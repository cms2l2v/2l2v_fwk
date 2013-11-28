#include <iostream>
#include <boost/shared_ptr.hpp>

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"



#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
//#include "UserCode/llvv_fwk/interface/DataEventSummaryHandler.h"
#include "UserCode/llvv_fwk/interface/llvvObjects.h"
#include "UserCode/llvv_fwk/interface/TMVAUtils.h"
#include "UserCode/llvv_fwk/interface/LeptonEfficiencySF.h"
#include "UserCode/llvv_fwk/interface/PDFInfo.h"
#include "UserCode/llvv_fwk/interface/MuScleFitCorrector.h"
#include "UserCode/llvv_fwk/interface/GammaWeightsHandler.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TEventList.h"
#include "TROOT.h"
#include "TNtuple.h"
#include <Math/VectorUtil.h>

using namespace std;



typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >::BetaVector BetaVector;


//
float getAngle(LorentzVector &a, LorentzVector &b)
{
  TVector3 mom1(a.px(),a.py(),a.pz());
  TVector3 mom2(b.px(),b.py(),b.pz());
  double cosine = mom1.Dot(mom2)/(mom1.Mag()*mom2.Mag());
  return acos(cosine);
}

//
std::vector<TString> getDijetCategories(double mjj,double etajj,std::vector<TString> &curTags, TString &mjjCat)
{
  if(mjj<250)               mjjCat="mjjq016";
  if(mjj>=250 && mjj<350)   mjjCat="mjjq033";
  if(mjj>=350 && mjj<450)   mjjCat="mjjq049";
  if(mjj>=450 && mjj<550)   mjjCat="mjjq066";
  if(mjj>=550 && mjj<750)   mjjCat="mjjq083";
  if(mjj>=750 && mjj<1000)  mjjCat="mjjq092";
  if(mjj>=1000)             mjjCat="mjjq100";
  
  //include new tags
  std::vector<TString> selTags;
  for(size_t i=0; i<curTags.size(); i++)
    {
      TString itag=curTags[i];
      selTags.push_back(itag);
      selTags.push_back(itag+mjjCat);
    }
  return selTags;
}


int main(int argc, char* argv[])
{
  //##############################################
  //########    GLOBAL INITIALIZATION     ########
  //##############################################

  // check arguments
  if(argc<2){ std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl; exit(0); }
  
  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();
  
  // configure the process
  const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");

  bool isMC = runProcess.getParameter<bool>("isMC");
  double xsec = runProcess.getParameter<double>("xsec");
  int mctruthmode=runProcess.getParameter<int>("mctruthmode");
  bool runPhotonSelection(mctruthmode==22 || mctruthmode==111);
  bool runLoosePhotonSelection(false);

  GammaWeightsHandler *gammaWgtHandler=0;
  if(runPhotonSelection) gammaWgtHandler=new GammaWeightsHandler(runProcess);
  
  float minJetPtToApply(30);

  std::vector<int> jacknifeCfg=runProcess.getParameter<std::vector<int> >("jacknife");
  int jacknife(jacknifeCfg[0]), jacks(jacknifeCfg[1]);
  if(jacknife>0 && jacks>0) cout << "Jacknife will be applied to every " << jacknife << " out of " << jacks << " events" << endl;
  
  std::vector<std::string> urls=runProcess.getParameter<std::vector<std::string> >("input");
  TString url = TString(urls[0]);
  TString baseDir    = runProcess.getParameter<std::string>("dirName");
  TString outFileUrl(gSystem->BaseName(url));
  outFileUrl.ReplaceAll(".root","");
  if(mctruthmode!=0) { outFileUrl += "_filt"; outFileUrl += mctruthmode; }
  TString outdir=runProcess.getParameter<std::string>("outdir");
  TString outUrl( outdir );
  gSystem->Exec("mkdir -p " + outUrl);
  bool filterOnlyEE(false), filterOnlyMUMU(false);
  if(!isMC)
    {
      if(url.Contains("DoubleEle")) filterOnlyEE=true;
      if(url.Contains("DoubleMu"))  filterOnlyMUMU=true;
    }
  bool isSingleMuPD(!isMC && url.Contains("SingleMu"));  
  bool isV0JetsMC(isMC && (url.Contains("DYJetsToLL_50toInf") || url.Contains("WJets")));
  bool isSignal(isMC && (url.Contains("VBFNLO") || url.Contains("lljj")) );

  TString outTxtUrl= outUrl + "/" + outFileUrl + ".txt";
  FILE* outTxtFile = NULL;
  if(!isMC)outTxtFile = fopen(outTxtUrl.Data(), "w");
  printf("TextFile URL = %s\n",outTxtUrl.Data());

  //residual corrections from Z+1 jet recoil
  std::vector<TGraph *> recoilResidualsGr;
  TFile *fIn=TFile::Open("/afs/cern.ch/user/p/psilva/work/CMSSW_5_3_3_patch2/src/CMGTools/HtoZZ2l2nu/data/recoilBalance.root");
  if(fIn){
    TGraph *gr=(TGraph *)fIn->Get("mumurecoilbalancevseta50toInfdydata2mc");
    if(gr) recoilResidualsGr.push_back( gr );
    gr=(TGraph *)fIn->Get("mumurecoilbalancevseta50toInfgdata2mc");
    if(gr) recoilResidualsGr.push_back( gr );
    fIn->Close();
    cout << "Read " << recoilResidualsGr.size() << " residual recoil systematics" << endl;
  }
  
  //Q^2 variations
  fIn=TFile::Open("/afs/cern.ch/user/p/psilva/work/CMSSW_5_3_3_patch2/src/CMGTools/HtoZZ2l2nu/data/vbfnloQ2weights.root");
  std::vector<TF1 *> Q2weightsGr;
  if(fIn && isSignal){
    TF1 *gr = (TF1 *)fIn->Get("q2upWgt_func");
    if(gr) Q2weightsGr.push_back(gr);
    gr = (TF1 *)fIn->Get("q2downWgt_func");
    if(gr) Q2weightsGr.push_back(gr);
    fIn->Close();
    cout << "Read " << Q2weightsGr.size() << " Q^2-scale re-weighting functions" << endl; 
  }
  
  //summary ntuple
  TString summaryTupleVarNames("ch:weight:nInitEvent:mjj:detajj:spt:setajj:dphijj:ystar:hardpt:fisher:llr:mva:ystar3:maxcjpt:ncjv:htcjv:ncjv15:htcjv15");
  TNtuple *summaryTuple = new TNtuple("ewkzp2j","ewkzp2j",summaryTupleVarNames);
  Float_t summaryTupleVars[summaryTupleVarNames.Tokenize(":")->GetEntriesFast()];
  summaryTuple->SetDirectory(0);

  //MVA
  bool useMVA = runProcess.getParameter<bool>("useMVA");
  TMVA::Reader *tmvaReader = 0;
  std::vector<string> tmvaMethods;
  std::vector<Float_t> tmvaDiscrVals(3,0);
  std::vector<std::string> tmvaVarNames;
  std::vector<Float_t> tmvaVars;
  if(useMVA)
    {
      edm::ParameterSet tmvaInput = runProcess.getParameter<edm::ParameterSet>("tmvaInput");
      std::string weightsDir      = tmvaInput.getParameter<std::string>("weightsDir");
      tmvaMethods                 = tmvaInput.getParameter<std::vector<std::string> >("methodList");
      tmvaDiscrVals.resize(tmvaMethods.size(),0);
      tmvaVarNames                = tmvaInput.getParameter<std::vector<std::string> >("varsList");
      tmvaVars.resize( tmvaVarNames.size(), 0 );

      //start the reader for the variables and methods
      tmvaReader = new TMVA::Reader( "!Color:!Silent" );
      for(size_t ivar=0; ivar<tmvaVarNames.size(); ivar++)   tmvaReader->AddVariable( tmvaVarNames[ivar], &tmvaVars[ivar] );

      //open the file with the method description
      for(size_t im=0; im<tmvaMethods.size(); im++)
	{
	  TString weightsFile = weightsDir + "/TMVAClassification_" + tmvaMethods[im] + TString(".weights.xml");
	  gSystem->ExpandPathName(weightsFile);
	  tmvaReader->BookMVA( tmvaMethods[im], weightsFile);
	}
    }

  //lepton efficiencies
  LeptonEfficiencySF lepEff;

  //systematics
  bool runSystematics                        = runProcess.getParameter<bool>("runSystematics");
  std::vector<TString> varNames(1,"");
  size_t nvarsToInclude(1);
  if(runSystematics && isMC && !runPhotonSelection)
    {
      varNames.push_back("_jerup"); varNames.push_back("_jerdown");
      varNames.push_back("_jesup"); varNames.push_back("_jesdown");
      varNames.push_back("_puup"); varNames.push_back("_pudown");
      if(isSignal)
	{
	  varNames.push_back("_q2up"); varNames.push_back("_q2down");
	  varNames.push_back("_pdfup"); varNames.push_back("_pdfdown");
	  varNames.push_back("_balanceup"); varNames.push_back("_balancedown");
	}
      nvarsToInclude=varNames.size();
      cout << nvarsToInclude << " systematics will be computed for this analysis" << endl;
    }
  if(runSystematics && runPhotonSelection) runLoosePhotonSelection=true;


  //##############################################
  //########    INITIATING HISTOGRAMS     ########
  //##############################################
  SmartSelectionMonitor mon;

  TH1 *h=mon.addHistogram( new TH1F ("eventflow", ";;Events", 5,0,5) );
  h->GetXaxis()->SetBinLabel(1,"#geq 2 leptons");
  h->GetXaxis()->SetBinLabel(2,"|M-M_{Z}|<15");
  h->GetXaxis()->SetBinLabel(3,"p_{T}^{ll}>50");
  h->GetXaxis()->SetBinLabel(4, "#eta^{ll}<1.44");
  h->GetXaxis()->SetBinLabel(5,"#geq 2 jets"); 

  mon.addHistogram( new TH1F("pthat",";#hat{p}_{T} [GeV];Events",50,0,1000) );
  mon.addHistogram( new TH1F("nup",";NUP;Events",10,0,10) );
  mon.addHistogram( new TH1F("nupfilt",";NUP;Events",10,0,10) );

  //pileup control
  mon.addHistogram( new TH1F( "nvtx",";Vertices;Events",50,0,50) ); 
  mon.addHistogram( new TH1F( "nvtxraw",";Vertices;Events",50,0,50) ); 
  mon.addHistogram( new TH1F( "rho",";#rho;Events",50,0,25) ); 
  mon.addHistogram( new TH1F( "rho25",";#rho(#eta<2.5);Events",50,0,25) ); 

  //lepton control
  mon.addHistogram( new TH1F( "leadpt", ";p_{T}^{l};Events", 50,0,500) );
  mon.addHistogram( new TH1F( "leadeta", ";#eta^{l};Events", 50,-2.6,2.6) );
  mon.addHistogram( new TH1F( "trailerpt", ";p_{T}^{l};Events", 50,0,500) );
  mon.addHistogram( new TH1F( "trailereta", ";#eta^{l};Events", 50,-2.6,2.6) );

  //balance histograms
  mon.addHistogram( new TH2F("ptllvsdphi",     ";p_{T}(ll) [GeV];#Delta #phi(ll,jet)",50,0,200,25,0,3.2) );
  mon.addHistogram( new TH2F("ptllvsdphipu",   ";p_{T}(ll) [GeV];#Delta #phi(ll,jet)",50,0,200,25,0,3.2) );
  mon.addHistogram( new TH2F("ptllvsdphitrue", ";p_{T}(ll) [GeV];#Delta #phi(ll,jet)",50,0,200,25,0,3.2) );
  for(size_t ireg=0; ireg<2; ireg++){
    TString regStr("lt15collinear");
    if(ireg==1) regStr="gt50back2back";
    mon.addHistogram( new TH1F("recoilbalancepumva"+regStr, "; PU discriminator; Jets", 100,-1,1) );
    mon.addHistogram( new TH1F("recoilbalancedrmean"+regStr, "; <#Delta R>; Jets", 100,0,0.5) );
    mon.addHistogram( new TH1F("recoilbalancebeta"+regStr, "; #beta; Jets", 100,0,1) );
    mon.addHistogram( new TH1F("recoilbalanceptrms"+regStr, "; RMS p_{T} [GeV]; Jets", 100,0,0.1) );
    mon.addHistogram( new TH1F("recoilbalance"+regStr, "; p_{T}(jet)/p_{T}; Jets", 100,0,5) );
    mon.addHistogram( new TH2F("recoilbalancevseta"+regStr, "; #eta(jet); <p_{T}(jet)/p_{T}>", 50,0,5,100,0,5) );
    TH2 *idH=(TH2 *)mon.addHistogram( new TH2F("recoilbalanceid"+regStr, "; Pseudo-rapidity; ID", 50,0,5,6,0,6) );
    idH->GetYaxis()->SetBinLabel(1,"no id");  
    idH->GetYaxis()->SetBinLabel(2,"PF loose");
    idH->GetYaxis()->SetBinLabel(3,"PF+PU");
    idH->GetYaxis()->SetBinLabel(4,"PF+sPU");
    idH->GetYaxis()->SetBinLabel(5,"PU");
    idH->GetYaxis()->SetBinLabel(6,"sPU");
    if(ireg==0)
      {
	mon.addHistogram( (TH2 *)idH->Clone("truejetsid") );
	mon.addHistogram( (TH2 *)idH->Clone("pujetsid") );
      }
  }
  
  //boson control
  mon.addHistogram( new TH1F( "qt",      ";p_{T}^{#gamma} [GeV];Events",500,0,1500));
  mon.addHistogram( new TH1F( "zpt",     ";p_{T}^{ll};Events", 50,0,500) );
  mon.addHistogram( new TH1F( "zptNM1",  ";p_{T}^{ll};Events", 50,0,500) );
  mon.addHistogram( new TH1F( "zeta",    ";#eta^{ll};Events", 50,-10,10) );
  mon.addHistogram( new TH1F( "zetaNM1", ";#eta^{ll};Events", 50,-10,10) );
  mon.addHistogram( new TH1F( "zy",      ";y^{ll};Events", 50,-6,6) );
  mon.addHistogram( new TH1F( "rapidity",";y^{ll};Events", 50,0,2) );
  mon.addHistogram( new TH1F( "zyNM1",   ";y^{ll};Events", 50,-6,6) );
  mon.addHistogram( new TH1F( "zmass",   ";M^{ll};Events", 100,40,250) );

  //jet control
  mon.addHistogram( new TH2F("jetptvseta", ";p_{T} [GeV];|#eta|;Events",50,0,400,25,0,5) );
  mon.addHistogram( new TH1F("jetgt3pt",        ";p_{T} [GeV];Events",50,0,400) );
  mon.addHistogram( new TH1F("jetgt3nhf",       ";Neutral hadron fraction;Events",50,0,1) );
  mon.addHistogram( new TH1F("jetgt3nemf",      ";Neutral e.m. fraction;Events",50,0,1) );
  mon.addHistogram( new TH1F("jetgt3chf",       ";Charged hadron fraction;Events",50,0,1) );
  mon.addHistogram( new TH1F("jetgt3ptrms",     ";p_{T} RMS [GeV];Events",50,0,0.2) );
  for(size_t ijesvar=0; ijesvar<3; ijesvar++)
    {
      TString pf("");
      if(ijesvar==1) pf+="_jesup";
      if(ijesvar==2) pf+="_jesdown";
      mon.addHistogram( new TH1F("jetpt"+pf       , ";p_{T} [GeV];Events",50,0,400) );
      mon.addHistogram( new TH1F("jeteta"+pf       , ";|#eta|;Events",25,0,5) );
      h=mon.addHistogram( new TH1F ("njets"+pf, ";Jet multiplicity;Events", 5,0,5) );
      h->GetXaxis()->SetBinLabel(1,"=0 jets");
      h->GetXaxis()->SetBinLabel(2,"=1 jets");
      h->GetXaxis()->SetBinLabel(3,"=2 jets");
      h->GetXaxis()->SetBinLabel(4,"=3 jets");
      h->GetXaxis()->SetBinLabel(5,"#geq 4 jets"); 
    }

  //vbf control
  mon.addHistogram( new TH2F("njetsvsavginstlumi",  ";;Jet multiplicity;Events",5,0,5,10,0,5000) );
  mon.addHistogram( new TH1F("vbfcandjeteta"     , ";#eta;Jets",                                 50,0,5) );
  mon.addHistogram( new TH1F("vbfcandjetpt"      , ";p_{T} [GeV];Jets",                        50,0,500) );
  mon.addHistogram( new TH1F("vbfcandjet1eta"    , ";#eta;Jets",                                 50,0,5) );
  mon.addHistogram( new TH1F("vbfcandjet1pt"     , ";p_{T} [GeV];Jets",                        50,0,500) );
  mon.addHistogram( new TH1F("vbfcandjet2eta"    , ";#eta;Jets",                                 50,0,5) );
  mon.addHistogram( new TH1F("vbfcandjet2pt"     , ";p_{T} [GeV];Jets",                        50,0,500) );
  mon.addHistogram( new TH1F("vbfcandjetdeta"    , ";|#Delta #eta|;Jets",                50,0,10) );  
  mon.addHistogram( new TH1F("vbfcandjetseta"    , ";#Sigma |#eta|;Jets",                50,0,15) );  
  mon.addHistogram( new TH1F("vbfcandjetetaprod" , ";#eta_{1} . #eta_{2};Jets",       100,-25,25) );
  mon.addHistogram( new TH1F("vbfhardpt"         , ";Hard p_{T} [GeV];Events",                   25,0,250) );
  mon.addHistogram( new TH1F("vbfspt"         , ";S_{p_{T}};Events",                   50,0,1) );
  h=mon.addHistogram( new TH1F("vbfcjv"          , ";Central jet count;Events",                     5,0,5) );
  h->GetXaxis()->SetBinLabel(1,"=0 jets");
  h->GetXaxis()->SetBinLabel(2,"=1 jets");
  h->GetXaxis()->SetBinLabel(3,"=2 jets");
  h->GetXaxis()->SetBinLabel(4,"=3 jets");
  h->GetXaxis()->SetBinLabel(5,"#geq 4 jets");
  mon.addHistogram( (TH1F *) h->Clone("vbfcjv15") );
  mon.addHistogram( (TH1F *) h->Clone("vbfcjv20") );
  mon.addHistogram( new TH1F("vbfhtcjv15"          , ";H_{T}(p_{T}>15) [GeV];Events",50,0,250) );
  mon.addHistogram( new TH1F("vbfhtcjv20"          , ";H_{T}(p_{T}>20) [GeV];Events",50,0,250) );
  mon.addHistogram( new TH1F("vbfhtcjv"            , ";H_{T}(p_{T}>30) [GeV];Events",50,0,250) );
  mon.addHistogram( new TH1F("vbfmaxcjvjpt"        , ";Central jet gap [GeV];Events",50,0,100) );
  mon.addHistogram( new TH1F("vbfystar3"           , ";y_{j3}-(y_{j1}+y_{j2})/2;Events",100,0,5) );

  Double_t mjjaxis[32];
  mjjaxis[0]=0.01;
  for(size_t i=1; i<20; i++)  mjjaxis[i]   =50*i;        //0-1000
  for(size_t i=0; i<5; i++)   mjjaxis[20+i]=1000+100*i; //1000-1500
  for(size_t i=0; i<=5; i++)   mjjaxis[25+i]=1500+300*i; //1500-5000  
  mjjaxis[31]=5000;
  mon.addHistogram( new TH1F("vbfmjj"            , ";M_{jj} [GeV];Events",               31,mjjaxis) );
  mon.addHistogram( new TH1F("vbfdphijj"         , ";#Delta#phi_{jj};Events",            20,0,3.5) );
  mon.addHistogram( new TH1F("vbfystar"           , ";#eta_{ll}-(#eta_{j1}+#eta_{j2})/2;Events",       50,0,5) );
  mon.addHistogram( new TH1F("vbfpt"             , ";p_{T}(visible) [GeV];Dijets",       50,0,500) );
  mon.addHistogram( new TH1F("met"               , ";E_{T}^{miss} [GeV]; Events",        50,0,500) );
  mon.addHistogram( new TH1F("metL"              , ";Axial E_{T}^{miss} [GeV]; Events",  50,-50,200) );
  
  std::vector<TH1 *> tmvaH;
  for(size_t im=0; im<tmvaMethods.size(); im++)
    tmvaH.push_back( mon.addHistogram( tmva::getHistogramForDiscriminator( tmvaMethods[im] ) ) );

  //statistical analysis
  std::vector<double> optim_Cuts2_jet_pt1; 
  std::vector<double> optim_Cuts2_jet_pt2; 
  for(double jet_pt1=minJetPtToApply;jet_pt1<=130;jet_pt1+=5)
    {
      for(double jet_pt2=minJetPtToApply;jet_pt2<=jet_pt1;jet_pt2+=5)
	{
	  optim_Cuts2_jet_pt1.push_back(jet_pt1);
	  optim_Cuts2_jet_pt2.push_back(jet_pt2);
	} 
    }
  TH2F* Hoptim_cuts2  =(TH2F*)mon.addHistogram(new TProfile2D("optim_cut2",      ";cut index;variable",       optim_Cuts2_jet_pt1.size(),0,optim_Cuts2_jet_pt1.size(), 2, 0, 2)) ;
  Hoptim_cuts2->GetYaxis()->SetBinLabel(1, "jpt1>");
  Hoptim_cuts2->GetYaxis()->SetBinLabel(2, "jpt2>");
  for(unsigned int index=0;index<optim_Cuts2_jet_pt1.size();index++){
    Hoptim_cuts2->Fill(index,0.0,optim_Cuts2_jet_pt1[index]); 
    Hoptim_cuts2->Fill(index,1.0,optim_Cuts2_jet_pt2[index]); 
  }

  TH1F* Hoptim_systs     =  (TH1F*) mon.addHistogram( new TH1F ("optim_systs"    , ";syst;", nvarsToInclude,0,nvarsToInclude) ) ;
  for(size_t ivar=0; ivar<nvarsToInclude; ivar++)
  {
    Hoptim_systs->GetXaxis()->SetBinLabel(ivar+1, varNames[ivar]);
    mon.addHistogram( new TH2F (TString("dijet_deta_shapes")+varNames[ivar],";cut index;|#Delta #eta|;Events",optim_Cuts2_jet_pt1.size(),0,optim_Cuts2_jet_pt1.size(),25,0,10) );
    if(tmvaH.size()) 
      for(size_t im=0; im<tmvaH.size(); im++)
	{
	  TString hname(tmvaMethods[im].c_str());
	  mon.addHistogram( new TH2F(hname+"_shapes" +varNames[ivar],";cut index;"+TString(tmvaH[im]->GetXaxis()->GetTitle())+";Events",
				     optim_Cuts2_jet_pt1.size(),0,optim_Cuts2_jet_pt1.size(),
				     tmvaH[im]->GetXaxis()->GetNbins(),tmvaH[im]->GetXaxis()->GetXmin(),tmvaH[im]->GetXaxis()->GetXmax()) );
	}
  }
  
  //##############################################
  //######## GET READY FOR THE EVENT LOOP ########
  //##############################################

  //open the file and get events tree
//  DataEventSummaryHandler evSummary;
//  TFile *file = TFile::Open(url);
//  printf("Looping on %s\n",url.Data());
//  if(file==0) return -1;
//  if(file->IsZombie()) return -1;
//  if( !evSummary.attach( (TTree *) file->Get(baseDir+"/data") ) ) { file->Close();  return -1; }
//  const Int_t totalEntries= evSummary.getEntries();

  fwlite::ChainEvent ev(urls);
  const Int_t totalEntries= ev.size();
 
  //MC normalization (to 1/pb)
  float nInitEvent=1.0;
  if(isMC){
     nInitEvent = (float)utils::getMergeableCounterValue(urls, "startCounter");
  }
  double xsecWeight = xsec/nInitEvent;
  if(!isMC) xsecWeight=1.0;

  //jet energy scale and uncertainties 
  TString jecDir = runProcess.getParameter<std::string>("jecDir");
  gSystem->ExpandPathName(jecDir);
  FactorizedJetCorrector *jesCor        = utils::cmssw::getJetCorrector(jecDir,isMC);
  JetCorrectionUncertainty *totalJESUnc = new JetCorrectionUncertainty((jecDir+"/MC_Uncertainty_AK5PFchs.txt").Data());

  //muon energy scale and uncertainties
  MuScleFitCorrector *muCor=getMuonCorrector(jecDir,url);
    
  //pdf info
  PDFInfo *mPDFInfo=0;
  if(isMC)
    {
      TString pdfUrl(url);
      pdfUrl.ReplaceAll(".root","_pdf.root");
      pdfUrl.ReplaceAll("/MC","/pdf/MC");
      mPDFInfo=new PDFInfo(pdfUrl,"cteq66.LHgrid");
      for(int i=0; i<mPDFInfo->numberPDFs(); i++)
	{
	  TString var("_"); var+=i;
	  mon.addHistogram( new TH1F("vbfcandjetdeta"+var    , ";|#Delta #eta|;Jets",                        50,0,10) );
	  mon.addHistogram( new TH1F("vbfcandjet1eta"+var    , ";#eta;Jets",                                 50,0,5) );
	  mon.addHistogram( new TH1F("vbfcandjet1pt"+var     , ";p_{T} [GeV];Jets",                        50,0,500) );
	  mon.addHistogram( new TH1F("vbfcandjet2eta"+var    , ";#eta;Jets",                                 50,0,5) );
	  mon.addHistogram( new TH1F("vbfcandjet2pt"+var     , ";p_{T} [GeV];Jets",                        50,0,500) );
	}
    }

  //pileup weighting: based on vtx for now...
  edm::LumiReWeighting* LumiWeights = NULL;
  utils::cmssw::PuShifter_t PuShifters;
  double PUNorm[] = {1,1,1};
  if(isMC){
     std::vector<double> dataPileupDistributionDouble = runProcess.getParameter< std::vector<double> >("datapileup");
     std::vector<float> dataPileupDistribution; for(unsigned int i=0;i<dataPileupDistributionDouble.size();i++){dataPileupDistribution.push_back(dataPileupDistributionDouble[i]);}
     std::vector<float> mcPileupDistribution;
     utils::getMCPileupDistribution(ev,dataPileupDistribution.size(), mcPileupDistribution);
     while(mcPileupDistribution.size()<dataPileupDistribution.size())  mcPileupDistribution.push_back(0.0);
     while(mcPileupDistribution.size()>dataPileupDistribution.size())dataPileupDistribution.push_back(0.0);

     LumiWeights= new edm::LumiReWeighting(mcPileupDistribution,dataPileupDistribution);
     PuShifters=utils::cmssw::getPUshifters(dataPileupDistribution,0.05);
     utils::getPileupNormalization(mcPileupDistribution, PUNorm, LumiWeights, PuShifters);
  }

  gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE

  //##############################################
  //########           EVENT LOOP         ########
  //##############################################
  //loop on all the events
  printf("Progressing Bar     :0%%       20%%       40%%       60%%       80%%       100%%\n");
  printf("Scanning the ntuple :");
  DuplicatesChecker duplicatesChecker;
  int nDuplicates(0);
  int step(totalEntries/50); 
  for( int iev=0; iev<totalEntries; iev++ ) 
    {
      if(iev%step==0){printf(".");fflush(stdout);}
      if(!isMC && jacknife>0 && jacks>0 && iev%jacks==jacknife) continue;
      
      //##############################################   EVENT LOOP STARTS   ##############################################
      //load the event content from tree
//      evSummary.getEntry(iev);
      ev.to(iev);
//      DataEventSummary &ev = evSummary.getEvent();
      if(!isMC && duplicatesChecker.isDuplicate( ev.eventAuxiliary().run(),ev.eventAuxiliary().luminosityBlock(),ev.eventAuxiliary().event()) ) { nDuplicates++; continue; }

      //get the collection of generated Particles
      fwlite::Handle< llvvGenEvent > genEventHandle;
      genEventHandle.getByLabel(ev, "llvvObjectProducersUsed");
      if(!genEventHandle.isValid()){printf("llvvGenEvent Object NotFound\n");continue;}
      llvvGenEvent genEv = *genEventHandle;


      if(isV0JetsMC){
	mon.fillHisto("nup","",genEv.nup,1);
	if(genEv.nup>5) continue;
	mon.fillHisto("nupfilt","",genEv.nup,1);
      }

      fwlite::Handle< llvvGenParticleCollection > genPartCollHandle;
      genPartCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
      if(!genPartCollHandle.isValid()){printf("llvvGenParticleCollection Object NotFound\n");continue;}
      llvvGenParticleCollection gen = *genPartCollHandle;

      fwlite::Handle< llvvLeptonCollection > leptonCollHandle;
      leptonCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
      if(!leptonCollHandle.isValid()){printf("llvvLeptonCollection Object NotFound\n");continue;}
      llvvLeptonCollection leptons = *leptonCollHandle;

      fwlite::Handle< llvvElectronInfoCollection > electronInfoCollHandle;
      electronInfoCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
      if(!electronInfoCollHandle.isValid()){printf("llvvElectronInfoCollection Object NotFound\n");continue;}

      fwlite::Handle< llvvMuonInfoCollection > muonInfoCollHandle;
      muonInfoCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
      if(!muonInfoCollHandle.isValid()){printf("llvvMuonInfoCollection Object NotFound\n");continue;}

      fwlite::Handle< llvvPhotonCollection > photonCollHandle;
      photonCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
      if(!photonCollHandle.isValid()){printf("llvvPhotonCollection Object NotFound\n");continue;}
      llvvPhotonCollection photons = *photonCollHandle;

      fwlite::Handle< llvvJetCollection > jetCollHandle;
      jetCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
      if(!jetCollHandle.isValid()){printf("llvvJetCollection Object NotFound\n");continue;}
      llvvJetExtCollection jets;
      for(unsigned int i=0;i<jetCollHandle->size();i++){jets.push_back(llvvJetExt((*jetCollHandle)[i]));}

      fwlite::Handle< llvvMet > metHandle;
      metHandle.getByLabel(ev, "llvvObjectProducersUsed", "pfMETPFlow"); 
      if(!metHandle.isValid()){printf("llvvMet Object NotFound\n");continue;}
      llvvMet met = *metHandle;

      fwlite::Handle< std::vector<bool> > triggerBitsHandle;
      triggerBitsHandle.getByLabel(ev, "llvvObjectProducersUsed", "triggerBits");
      if(!triggerBitsHandle.isValid()){printf("triggerBits Object NotFound\n");continue;}
      std::vector<bool> triggerBits = *triggerBitsHandle;

      fwlite::Handle< std::vector<int> > triggerPrescalesHandle;
      triggerPrescalesHandle.getByLabel(ev, "llvvObjectProducersUsed", "triggerPrescales");
      if(!triggerPrescalesHandle.isValid()){printf("triggerPrescales Object NotFound\n");continue;}
      std::vector<int> triggerPrescales = *triggerPrescalesHandle;

      fwlite::Handle< double > rhoHandle;
      rhoHandle.getByLabel(ev, "kt6PFJets", "rho");
      if(!rhoHandle.isValid()){printf("rho Object NotFound\n");continue;}
      double rho = *rhoHandle;

      fwlite::Handle< double > rho25Handle;
      rho25Handle.getByLabel(ev, "kt6PFJetsCentral", "rho");
      if(!rho25Handle.isValid()){printf("rho25 Object NotFound\n");continue;}
      double rho25 = *rho25Handle;

      fwlite::Handle< int > nvtxHandle;
      nvtxHandle.getByLabel(ev, "llvvObjectProducersUsed", "nvtx");
      if(!nvtxHandle.isValid()){printf("nvtx Object NotFound\n");continue;}
      int nvtx = *nvtxHandle;

      //require compatibilitiy of the event with the PD
      bool eeTrigger          = triggerBits[0];
      bool muTrigger          = triggerBits[6];
      bool mumuTrigger        = triggerBits[2] || triggerBits[3] || muTrigger;
      if(filterOnlyEE)   { mumuTrigger=false; }
      if(filterOnlyMUMU) { eeTrigger=false;   }
      if(isSingleMuPD)   { eeTrigger=false; if( mumuTrigger || !muTrigger ) mumuTrigger= false;  }
      
      bool hasPhotonTrigger(false);
      float triggerPrescale(1.0),triggerThreshold(0);
      if(runPhotonSelection){
	  eeTrigger=false; mumuTrigger=false;
	  for(size_t itrig=10; itrig>=7; itrig--)
	    {
	      if(!triggerBits[itrig]) continue;
	      hasPhotonTrigger=true;
	      triggerPrescale=triggerPrescales[itrig];
	      if(itrig==10) triggerThreshold=90;
	      if(itrig==9)  triggerThreshold=75;
	      if(itrig==8)  triggerThreshold=50;
	      if(itrig==7)  triggerThreshold=36;
	      break;
	    }
      }

      //pileup weight
      float weight = 1.0;
      double TotalWeight_plus = 1.0;
      double TotalWeight_minus = 1.0;
      float puWeight(1.0);
      if(isMC){
        puWeight          = LumiWeights->weight(genEv.ngenITpu) * PUNorm[0];
	weight            = xsecWeight*puWeight;
        TotalWeight_plus  = PuShifters[utils::cmssw::PUUP  ]->Eval(genEv.ngenITpu) * (PUNorm[2]/PUNorm[0]);
        TotalWeight_minus = PuShifters[utils::cmssw::PUDOWN]->Eval(genEv.ngenITpu) * (PUNorm[1]/PUNorm[0]);
      }

      //
      //
      // BELOW FOLLOWS THE ANALYSIS OF THE MAIN SELECTION WITH N-1 PLOTS
      //
      //

      //
      // photon selection
      //
      llvvPhotonCollection selPhotons;
      if(runPhotonSelection)
	{
	  //filter out number of photons generated in the hard process (no double counting)
	  int ngenpho(0);
	  for(size_t igen=0; igen<gen.size(); igen++)
	    {
	      if(gen[igen].id!=22) continue;
	      if(gen[igen].status!=3) continue;
	      ngenpho++;
	    }
	  //if(mctruthmode==111 && ngenpho>0) continue;
	  //if(mctruthmode==22 && ngenpho==0) continue;

	  //select the photons
	  for(size_t ipho=0; ipho<photons.size(); ipho++)
	    {
	      double pt=photons[ipho].pt();
	      double eta=photons[ipho].sceta;

	      //if systematics are active loosen the selection to the medium working point
	      Int_t idbits( photons[ipho].idbits );
	      bool hasPhotonId( (idbits >> 2 ) & 0x1 );
	      //double gIso    = photons[ipho].getVal("gIso03");
	      // double gArea   = utils::cmssw::getEffectiveArea(22,eta,3,"gIso");	      
	      double chIso   = photons[ipho].chIso03;
	      // 	      double chArea  = utils::cmssw::getEffectiveArea(22,eta,3,"chIso");
	      // 	      double nhIso   = photons[ipho].getVal("nhIso03");
	      // 	      double nhArea  = utils::cmssw::getEffectiveArea(22,eta,3,"nhIso");
	      
	      //select the photon
	      if(pt<triggerThreshold || fabs(eta)>1.4442 ) continue;
	      bool passId(true);
	      if(!hasPhotonId) passId=false;
	      if( photons[ipho].r9<0.9 ) passId=false;
	      if(!passId) continue;
	      bool passIso(true);
	      if(runLoosePhotonSelection){
		passIso &= ((chIso/pt)<0.20);
		// passIso &= (TMath::Max(chIso-chArea*rho,0.0) < 1.5); 
		// passIso &= (TMath::Max(nhIso-nhArea*rho,0.0) < 1.0+0.04*pt); 
		// passIso &= (TMath::Max(gIso-gArea*rho,  0.0) < 0.7+0.005*pt); 
	      }
	      else{
		passIso &= ((chIso/pt)<0.10);
	      }
	      if(!passIso) continue; 
	      selPhotons.push_back(photons[ipho]);
	    }
	}


      //
      // LEPTON ANALYSIS
      //
      llvvLeptonCollection selLeptons;
      for(size_t ilep=0; ilep<leptons.size(); ilep++)
	{
	  bool passKin(true),passId(true),passIso(true);
	  int lid=leptons[ilep].id;

	  //apply muon corrections
	  if(lid==13 && muCor){
            TLorentzVector p4(leptons[ilep].px(), leptons[ilep].py(), leptons[ilep].pz(), leptons[ilep].energy());
	    muCor->applyPtCorrection(p4 , lid<0 ? -1 :1 );
	    if(isMC) muCor->applyPtSmearing(p4, lid<0 ? -1 : 1, false);
            leptons[ilep].SetPxPyPzE(p4.Px(), p4.Py(), p4.Pz(), p4.Energy());
	  }

	  //no need for charge info any longer
	  lid=abs(lid);
	  TString lepStr( lid==13 ? "mu" : "e");

	  //veto nearby photon (loose electrons are many times photons...)
	  double minDRlg(9999.);
	  for(size_t ipho=0; ipho<selPhotons.size(); ipho++)
	    minDRlg=TMath::Min(minDRlg,deltaR(leptons[ilep],selPhotons[ipho]));
	  if(minDRlg<0.1) continue;

	  
	  //kinematics
//	  float leta = lid==11 ? leptons[ilep].getVal("sceta") : leptons[ilep].eta();
        float leta = lid==11 ? leptons[ilep].electronInfoRef->sceta : leptons[ilep].eta();
	  if(leptons[ilep].pt()<20)                   passKin=false;
	  if(leta> (lid==11 ? 2.5 : 2.4) )            passKin=false;
	  if(lid==11 && (leta>1.4442 && leta<1.5660)) passKin=false;

	  //id
	  Int_t idbits = leptons[ilep].idbits;
	  if(lid==11){
	    if(leptons[ilep].electronInfoRef->isConv)              passId=false;
	    bool isLoose = ((idbits >> 4) & 0x1);
	    if(!isLoose)                                   passId=false;
 	  }
	  else{
	    bool isLoose    = ((idbits >> 8) & 0x1);
	    if(!isLoose)                                   passId=false;
	  }

	  //isolation
	  Float_t gIso    = lid==11 ? leptons[ilep].gIso03    : leptons[ilep].gIso04;
	  Float_t chIso   = lid==11 ? leptons[ilep].chIso03   : leptons[ilep].chIso04;
	  Float_t puchIso = lid==11 ? leptons[ilep].puchIso03 : leptons[ilep].puchIso04;  
	  Float_t nhIso   = lid==11 ? leptons[ilep].nhIso03   : leptons[ilep].nhIso04;
	  float relIso= lid==11 ?
	    (TMath::Max(nhIso+gIso-rho*utils::cmssw::getEffectiveArea(11,leptons[ilep].electronInfoRef->sceta),double(0.))+chIso)/leptons[ilep].pt() :
	    (TMath::Max(nhIso+gIso-0.5*puchIso,0.)+chIso)/leptons[ilep].pt()
	    ;
	  if(lid==11){
	    if(relIso>0.15)                                passIso=false;
	  }
	  else{
	    if(relIso>0.20)                                passIso=false;
	  }
	  
	  if(!passId || !passIso || !passKin) continue;
	  selLeptons.push_back(leptons[ilep]);
	}
      std::sort(selLeptons.begin(), selLeptons.end(), sort_llvvObjectByPt);

      //at this point check if it's worth continuig
      if(runPhotonSelection && (selLeptons.size()!=0 || selPhotons.size()!=1)) continue;
      if(!runPhotonSelection && selLeptons.size()<2) continue;

      //apply data/mc correction factors
      //prepare the tag's vectors for histo filling
      std::vector<TString> chTags;
      int dilId(1);
      if(!runPhotonSelection)
	{
 	  for(size_t ilep=0; ilep<2; ilep++)
	    {
	      dilId *= selLeptons[ilep].id;
	      int id(abs(selLeptons[ilep].id));
	      weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[ilep].pt(), selLeptons[ilep].eta(), id,  id ==11 ? "loose" : "loose" ).first : 1.0;
	    }
     
	  //check the channel
	  if( abs(dilId)==121 && eeTrigger)   chTags.push_back("ee");
	  if( abs(dilId)==169 && mumuTrigger) chTags.push_back("mumu"); 
	}
      else{
	if(hasPhotonTrigger) {
	  dilId=22;
	  chTags.push_back("ee");
	  chTags.push_back("mumu");
	  weight *= triggerPrescale;
	}
      }
      if(chTags.size()==0) continue;
    
      //
      // DILEPTON ANALYSIS
      //
      LorentzVector leadingLep, trailerLep, zll;
      if(runPhotonSelection){leadingLep = selPhotons[0]; trailerLep=selPhotons[0]; zll=selPhotons[0];
      }else if(selLeptons[0].pt()>selLeptons[1].pt()){ leadingLep=selLeptons[0]; trailerLep=selLeptons[1]; zll=leadingLep+trailerLep;
      }else{                                           leadingLep=selLeptons[1]; trailerLep=selLeptons[0]; zll=leadingLep+trailerLep;
      } 
      float zy(zll.Rapidity());
      bool passZmass(runPhotonSelection || (fabs(zll.mass()-91)<15));
      bool passZpt(zll.pt()>50);
      bool passZeta(fabs(zy)<1.4442);

      //generator level
      LorentzVector genll;
      for(size_t ig=0; ig<gen.size(); ig++)
	{
	  int pid=abs(gen[ig].id);
	  if(pid!=11 && pid!=13) continue;
	  genll += gen[ig];
	}
      
      //
      //JET/MET ANALYSIS
      //
      llvvJetExtCollection selJets, selJetsNoId;
      int njets(0);
      for(size_t ijet=0; ijet<jets.size(); ijet++) 
	{
	  //correct jet
	  float toRawSF=jets[ijet].torawsf;
	  LorentzVector rawJet(jets[ijet]*toRawSF);
	  jesCor->setJetEta(rawJet.eta());
	  jesCor->setJetPt(rawJet.pt());
	  jesCor->setJetA(jets[ijet].area);
	  jesCor->setRho(rho);
	  float newJECSF=jesCor->getCorrection();
	  jets[ijet].SetPxPyPzE(rawJet.px(),rawJet.py(),rawJet.pz(),rawJet.energy());
	  jets[ijet] *= newJECSF;
	  jets[ijet].torawsf = 1./newJECSF;
	  if(jets[ijet].pt()<15 || fabs(jets[ijet].eta())>4.7 ) continue;
	  
	  //cross-clean with selected leptons and photons
	  double minDRlj(9999.),minDRlg(9999.);
          for(size_t ilep=0; ilep<selLeptons.size(); ilep++)
            minDRlj = TMath::Min( minDRlj, deltaR(jets[ijet],selLeptons[ilep]) );
	  for(size_t ipho=0; ipho<selPhotons.size(); ipho++)
	    minDRlg = TMath::Min( minDRlg, deltaR(jets[ijet],selPhotons[ipho]) );
	  if(minDRlj<0.4 || minDRlg<0.4) continue;
	  
	  //jet id
	  // float pumva=jets[ijet].puMVA;
	  Int_t idbits=jets[ijet].idbits;
	  bool passPFloose( ((idbits>>0) & 0x1));
	  int puId( ( idbits >>3 ) & 0xf );
	  bool passLoosePuId( ( puId >> 2) & 0x1);
	  int simplePuId( ( idbits >>7 ) & 0xf );
	  bool passLooseSimplePuId(  ( simplePuId >> 2) & 0x1);
	  TString jetType( jets[ijet].genj.pt()>0 ? "truejetsid" : "pujetsid" );
	  if(jets[ijet].pt()>30)
	    {
	      mon.fillHisto(jetType,chTags,fabs(jets[ijet].eta()),0);
	      if(passPFloose)                        mon.fillHisto(jetType,chTags,fabs(jets[ijet].eta()),1);
	      if(passPFloose && passLoosePuId)       mon.fillHisto(jetType,chTags,fabs(jets[ijet].eta()),2);
	      if(passPFloose && passLooseSimplePuId) mon.fillHisto(jetType,chTags,fabs(jets[ijet].eta()),3);
	      if(passLoosePuId)                      mon.fillHisto(jetType,chTags,fabs(jets[ijet].eta()),4);
	      if(passLooseSimplePuId)                mon.fillHisto(jetType,chTags,fabs(jets[ijet].eta()),5);
	    }
	  	
	  //add scale/resolution uncertainties
	  std::vector<float> smearPt=utils::cmssw::smearJER(jets[ijet].pt(),jets[ijet].eta(),jets[ijet].genj.pt());
	  jets[ijet].jer     = isMC ? smearPt[0] : jets[ijet].pt();
	  jets[ijet].jerup   = isMC ? smearPt[1] : jets[ijet].pt();
	  jets[ijet].jerdown = isMC ? smearPt[2] : jets[ijet].pt();
	  smearPt=utils::cmssw::smearJES(jets[ijet].pt(),jets[ijet].eta(), totalJESUnc);
	  jets[ijet].jesup   = isMC ? smearPt[0] : jets[ijet].pt();
	  jets[ijet].jesdown = isMC ? smearPt[1] : jets[ijet].pt();

	  selJetsNoId.push_back(jets[ijet]);
	  if(passPFloose && passLooseSimplePuId){
	    selJets.push_back(jets[ijet]);
	    if(jets[ijet].pt()>minJetPtToApply) njets++;
	  }
	}
      std::sort(selJets.begin(), selJets.end(), sort_llvvObjectByPt);
      
      //analyze dijets and activity in the dijet rapidity distance
      float maxPt(0), minPt(0), maxEta(0), minEta(0), maxAbsEta(0), minAbsEta(0);
      float dphijj(0), detajj(0), setajj(0), etaprod(0), ystar(0), mjj(0),ptjj(0),spt(0);
      float hardpt(0);
      int ncjv(0), ncjv15(0),ncjv20(0), htcjv(0), htcjv15(0),htcjv20(0);
      float pt3(0), ystar3(0);
      float ptmiss(0),metL(0);
      if(njets>=2)
	{
	  llvvJet jet1=selJets[0];
	  llvvJet jet2=selJets[1];
	  maxPt=max(jet1.pt(),jet2.pt());
	  minPt=min(jet1.pt(),jet2.pt());
	  maxAbsEta=max(fabs(jet1.eta()),fabs(jet2.eta()));
	  minAbsEta=min(fabs(jet1.eta()),fabs(jet2.eta()));
	  maxEta=max(jet1.eta(),jet2.eta());
	  minEta=min(jet1.eta(),jet2.eta());
	  detajj=fabs(maxEta-minEta);
	  setajj=fabs(maxEta)+fabs(minEta);
	  dphijj=deltaPhi(jet1.phi(),jet2.phi());
	  etaprod=jet1.eta()*jet2.eta(); 
	  ystar=zll.Rapidity()-0.5*(jet1.Rapidity()+jet2.Rapidity());
	  
	  LorentzVectorF vbfSyst=jet1+jet2;
	  mjj=vbfSyst.mass();
	  ptjj=vbfSyst.pt();
	  spt=vbfSyst.pt()/(jet1.pt()+jet2.pt());

	  LorentzVectorF hardSyst=vbfSyst+zll; 
	  hardpt=hardSyst.pt();
	  
	  TVector2 boson2(zll.px(),zll.py());
	  TVector2 met2(met.px(),met.py());
	  ptmiss=met.pt();
	  metL=boson2*met2; metL /= -zll.pt();
	  
	  //visible system rest frame
	  //BetaVector vbfBoost         = hardSyst.BoostToCM();
	  //LorentzVector jet1_cm    = ROOT::Math::VectorUtil::boost( jet1, vbfBoost );
	  //LorentzVector jet2_cm    = ROOT::Math::VectorUtil::boost( jet2, vbfBoost );
	  //LorentzVector jj_cm      = ROOT::Math::VectorUtil::boost( vbfSyst, vbfBoost );
	  //LorentzVector z_cm       = ROOT::Math::VectorUtil::boost( zll,  vbfBoost );
	  
	  for(size_t iotherjet=2; iotherjet<selJets.size(); iotherjet++){
	    bool isInRapGap(selJets[iotherjet].eta()>minEta && selJets[iotherjet].eta()<maxEta);
	    if(!isInRapGap) continue;
	    float ipt=selJets[iotherjet].pt();
	    if(ipt>pt3)
	      {
		pt3=ipt;
		ystar3=selJets[iotherjet].Rapidity()-0.5*(jet1.Rapidity()+jet2.Rapidity());
	      }
	    
	    ncjv15++;
	    htcjv15 += ipt;
	    if(ipt>20)
	      {
		ncjv20++;
		htcjv20 += ipt;
		if(ipt>30)
		  {
		    ncjv++;
		    htcjv += ipt;
		  }
	      }
	  }
	}

      //set the variables to be used in the MVA evaluation (independently of its use)
      for(size_t ivar=0; ivar<tmvaVarNames.size(); ivar++) 
	{
	  std::string variable     = tmvaVarNames[ivar];
	  if(variable=="mjj")        tmvaVars[ivar]=mjj;
	  if(variable=="detajj")     tmvaVars[ivar]=detajj;
          if(variable=="spt")        tmvaVars[ivar]=spt;
	}
      if(tmvaReader) {
	for(size_t im=0; im<tmvaMethods.size(); im++)
	  tmvaDiscrVals[im]=tmvaReader->EvaluateMVA( tmvaMethods[im] );
      }
            
      //
      // NOW FOR THE CONTROL PLOTS
      //
      //start analysis
      for(size_t ich=0; ich<chTags.size(); ich++)
	{
	  std::vector<TString> tags(1,chTags[ich]);

	  mon.fillHisto("eventflow",tags,0,weight);
	  if(passZmass)                                      mon.fillHisto("eventflow",tags,1,weight);
	  if(passZmass && passZpt)                           mon.fillHisto("eventflow",tags,2,weight);
	  if(passZmass && passZpt && passZeta)               mon.fillHisto("eventflow",tags,3,weight);
	  if(passZmass && passZpt && passZeta && njets>1)    mon.fillHisto("eventflow",tags,4,weight);

	  mon.fillHisto("zmass",    tags, zll.mass(), weight);  
	  if(passZmass){
	
	    //pu control
            mon.fillHisto("nvtx"     ,   tags, nvtx,      weight);
	    mon.fillHisto("nvtxraw"  ,   tags, nvtx,      weight/puWeight);
	    mon.fillHisto("rho"      ,   tags, rho,       weight);
            mon.fillHisto("rho25"    ,   tags, rho25,     weight);
	
	    //Z kinematics control
	    mon.fillHisto("leadpt"      ,   tags, leadingLep.pt(), weight);      
	    mon.fillHisto("trailerpt"   ,   tags, trailerLep.pt(), weight);      
	    mon.fillHisto("leadeta"     ,   tags, TMath::Max(fabs(leadingLep.eta()),fabs(trailerLep.eta())), weight);      
	    mon.fillHisto("trailereta"  ,   tags, TMath::Min(fabs(leadingLep.eta()),fabs(trailerLep.eta())), weight);      

	    mon.fillHisto("zpt"      , tags, zll.pt(),  weight);      
	    mon.fillHisto("zeta"     , tags, zll.eta(), weight);
	    mon.fillHisto("zy"       , tags, zy,        weight);
	    
	    //balance control
	    if(njets==1)
	      {
		//set as pu if no matched gen jet
		bool isPUjet( selJetsNoId[0].genj.pt()==0 ); 
		
		//kinematics
		float balance = selJetsNoId[0].pt()/zll.pt();
		float dphi    = fabs( deltaPhi(selJetsNoId[0].phi(),zll.phi()) );
		mon.fillHisto("ptllvsdphi",                                    tags,zll.pt(),dphi,weight);
		mon.fillHisto(TString("ptllvsdphi")+(isPUjet ? "pu" : "true"), tags,zll.pt(),dphi,weight);
	
		TString regStr("");
		if(dphi<1   && zll.pt()<15) regStr="lt15collinear";
		if(dphi>2.7 && zll.pt()>50) regStr="gt50back2back";
		if(regStr!="")
		  {
		    float drmean( selJetsNoId[0].dRMean) ;
		    float beta( selJetsNoId[0].beta );
		    float ptrms( selJetsNoId[0].ptRMS );
		    float pumva( selJetsNoId[0].puMVA );

		    //ids
		    Int_t idbits=selJetsNoId[0].idbits;
		    bool passPFloose ( ((idbits>>0) & 0x1) );
		    int puId((idbits>>3) & 0xf);
		    bool passLoosePuId( ( puId >> 2) & 0x1);
		    int simplePuId( ( idbits >>7 ) & 0xf );
		    bool passLooseSimplePuId(  ( simplePuId >> 2) & 0x1);

		    mon.fillHisto("recoilbalancepumva"+regStr, tags, pumva, weight);
		    mon.fillHisto("recoilbalancedrmean"+regStr, tags, drmean, weight);
		    mon.fillHisto("recoilbalancebeta"+regStr,   tags, beta,   weight);
		    mon.fillHisto("recoilbalanceptrms"+regStr,  tags, ptrms,  weight);

		    mon.fillHisto("recoilbalance"+regStr,      tags,balance, weight);
		    mon.fillHisto("recoilbalancevseta"+regStr, tags,fabs(selJetsNoId[0].eta()), balance, weight);
		    mon.fillHisto("recoilbalanceid"+regStr,    tags,fabs(selJetsNoId[0].eta()),0, weight);
		    if(passPFloose){
		      mon.fillHisto("recoilbalanceid"+regStr,tags,fabs(selJetsNoId[0].eta()),1, weight);
		      if(passLoosePuId)     mon.fillHisto("recoilbalanceid"+regStr,tags,fabs(selJetsNoId[0].eta()),2, weight);
		      if(passLooseSimplePuId)  mon.fillHisto("recoilbalanceid"+regStr,tags,fabs(selJetsNoId[0].eta()),3, weight);
		    }
		    if(passLoosePuId)     mon.fillHisto("recoilbalanceid"+regStr,tags,fabs(selJetsNoId[0].eta()),4, weight);
		    if(passLooseSimplePuId)  mon.fillHisto("recoilbalanceid"+regStr,tags,fabs(selJetsNoId[0].eta()),5, weight);
		  }
	      }
	    //end balance control
	    
	    if(passZpt && passZeta){
	  
	      //analyze dilepton kinematics
	      mon.fillHisto("leadeta"   ,  tags, leadingLep.eta(), weight);
	      mon.fillHisto("leadpt"    ,  tags, leadingLep.pt(),  weight);
	      mon.fillHisto("trailereta",  tags, trailerLep.eta(), weight);
	      mon.fillHisto("trailerpt",   tags, trailerLep.pt(),  weight);
	    
	      //analyze jet kinematics
	      for(size_t ijesvar=0; ijesvar<3; ijesvar++)
		{
		  TString pf("");
		  int njets_ivar(0);
		  for(size_t ijet=0; ijet<selJets.size(); ijet++)
		    {
		      float pt( selJets[ijet].pt() );
		      if(ijesvar==1) { pf="_jesup";   pt=selJets[ijet].jesup; }
		      if(ijesvar==2) { pf="_jesdown"; pt=selJets[ijet].jesdown; }
		      if(pt>minJetPtToApply){
			njets_ivar++;
			if(ijesvar==0){
			  mon.fillHisto("jetptvseta",  tags,pt,fabs(selJets[ijet].eta()),weight);
			  mon.fillHisto("jetgt3pt",    tags,pt,weight);
			  mon.fillHisto("jetgt3nhf",   tags,selJets[ijet].neutHadFrac,weight);
			  mon.fillHisto("jetgt3nemf",  tags,selJets[ijet].neutEmFrac, weight);
			  mon.fillHisto("jetgt3chf",   tags,selJets[ijet].chHadFrac,  weight);
			  mon.fillHisto("jetgt3ptrms", tags,selJets[ijet].ptRMS,      weight);
			}
			mon.fillHisto("jetpt"+pf,  tags, pt, weight);
			mon.fillHisto("jeteta"+pf, tags, fabs(selJets[ijet].eta()), weight);
		      }
		    }
		  mon.fillHisto("njets"+pf,tags, njets_ivar, weight);
		}
		
	      //signal region
	      float photonWeight(1.0);
	      if(njets>=2)
		{
		  TString mjjCat("");
		  std::vector<TString> selTags;
		  selTags = getDijetCategories(mjj,detajj,tags,mjjCat);

		  //re-weight for photons if needed
		  if(gammaWgtHandler!=0) {
		    mjjCat.ReplaceAll("mjjq100","mjjq092");
		    std::vector<Float_t> photonVars;
		    photonVars.push_back(selPhotons[0].pt());
		    photonWeight = gammaWgtHandler->getWeightFor(photonVars,chTags[ich]+mjjCat);
		  }
		  float catWeight=weight*photonWeight;

		  //veto events with very large weights in simulation
		  if(isMC && catWeight>5) catWeight=0;

		  //save for further analysis
		  if(mjj>200) {
		    summaryTupleVars[0]=dilId;  
		    summaryTupleVars[1]=catWeight*xsecWeight;    
		    summaryTupleVars[2]=nInitEvent;
		    summaryTupleVars[3]=mjj;     
		    summaryTupleVars[4]=detajj;               
		    summaryTupleVars[5]=spt;
		    summaryTupleVars[6]=setajj;  
		    summaryTupleVars[7]=dphijj;  
		    summaryTupleVars[8]=ystar;
		    summaryTupleVars[9]=hardpt; 
		    for(size_t im=0; im<tmvaMethods.size(); im++) summaryTupleVars[10+im]=tmvaDiscrVals[im];
		    summaryTupleVars[13]=ystar3; 
		    summaryTupleVars[14]=pt3; 
		    summaryTupleVars[15]=ncjv;
		    summaryTupleVars[16]=htcjv;  
		    summaryTupleVars[17]=ncjv15;  
		    summaryTupleVars[18]=htcjv15;
		    summaryTuple->Fill(summaryTupleVars);
		    for(size_t im=0; im<tmvaMethods.size(); im++) mon.fillHisto(tmvaMethods[im], selTags, tmvaDiscrVals[im], catWeight);
		  } 
	      
		  mon.fillHisto("qt",                 selTags, zll.pt(), catWeight,true);      
		  mon.fillHisto("rapidity"     ,      selTags, fabs(zy),    catWeight);
//FIXME		  mon.fillHisto("njetsvsavginstlumi", selTags, njets,ev.eventAuxiliary().luminosityBlock().avgInsDelLumi(),catWeight);
		  mon.fillHisto("vbfcandjetpt",       selTags, maxPt,catWeight);
		  mon.fillHisto("vbfcandjetpt",       selTags, minPt,catWeight);
		  mon.fillHisto("vbfcandjet1pt",      selTags, maxPt,catWeight);
		  mon.fillHisto("vbfcandjet2pt",      selTags, minPt,catWeight);
		  mon.fillHisto("vbfcandjet1eta",     selTags, maxAbsEta, catWeight);
		  mon.fillHisto("vbfcandjet2eta",     selTags, minAbsEta, catWeight);
		  mon.fillHisto("vbfcandjeteta",      selTags, maxAbsEta, catWeight);
		  mon.fillHisto("vbfcandjeteta",      selTags, minAbsEta, catWeight);
		  mon.fillHisto("vbfcandjetdeta",     selTags, detajj,catWeight);
		  mon.fillHisto("vbfcandjetseta",     selTags, setajj,catWeight);
		  mon.fillHisto("vbfcandjetetaprod",  selTags, etaprod,catWeight);
		  mon.fillHisto("vbfmjj",             selTags, mjj,catWeight,true);
		  mon.fillHisto("vbfhardpt",          selTags, hardpt,catWeight);
		  mon.fillHisto("vbfspt",             selTags, spt,catWeight);
		  mon.fillHisto("vbfdphijj",          selTags, fabs(dphijj),catWeight);
		  mon.fillHisto("vbfystar",           selTags, ystar,catWeight);
		  mon.fillHisto("vbfpt",              selTags, ptjj,catWeight);
		  mon.fillHisto("met",                selTags, ptmiss,catWeight);
		  mon.fillHisto("metL",               selTags, metL,catWeight);
		  if(ncjv15){
		    mon.fillHisto("vbfmaxcjvjpt",        selTags, pt3,catWeight);
		    mon.fillHisto("vbfystar3",           selTags, ystar3,catWeight);
		  }
		  mon.fillHisto("vbfcjv",              selTags, ncjv,catWeight);
		  mon.fillHisto("vbfcjv15",            selTags, ncjv15,catWeight);
		  mon.fillHisto("vbfcjv20",            selTags, ncjv20,catWeight);
		  mon.fillHisto("vbfhtcjv",            selTags, htcjv,catWeight);
		  mon.fillHisto("vbfhtcjv15",          selTags, htcjv15,catWeight);
		  mon.fillHisto("vbfhtcjv20",          selTags, htcjv20,catWeight);
	      
		  if(isMC && mPDFInfo)
		    {
		      std::vector<float> wgts=mPDFInfo->getWeights(iev);
		      for(size_t ipw=0; ipw<wgts.size(); ipw++) 
			{
			  TString var("_"); var+=ipw;
			  mon.fillHisto("vbfcandjetdeta"+var,     selTags, fabs(detajj),catWeight*wgts[ipw]);
			  mon.fillHisto("vbfcandjet1pt"+var,      tags, maxPt,catWeight*wgts[ipw]);
			  mon.fillHisto("vbfcandjet2pt"+var,      tags, minPt,catWeight*wgts[ipw]);
			  mon.fillHisto("vbfcandjet1eta"+var,     tags, maxAbsEta,catWeight*wgts[ipw]);
			  mon.fillHisto("vbfcandjet2eta"+var,     tags, minAbsEta,catWeight*wgts[ipw]);
			}
		    }
	      	      
		}//end nAJetsLoose
  
	      //STATISTICAL ANALYSIS
	      float Q2Weight_plus(1.0), Q2Weight_down(1.0);
	      float PDFWeight_plus(1.0), PDFWeight_down(1.0);
	      if(isSignal)
		{
		  if(mPDFInfo)
		    {
		      std::vector<float> wgts=mPDFInfo->getWeights(iev);
		      for(size_t ipw=0; ipw<wgts.size(); ipw++)
			{
			  PDFWeight_plus = TMath::Max(PDFWeight_plus,wgts[ipw]);
			  PDFWeight_down = TMath::Min(PDFWeight_down,wgts[ipw]);
			}
		    }
		  if(Q2weightsGr.size()==2)
		    {
		      Q2Weight_plus = Q2weightsGr[0]->Eval(genll.pt());
		      Q2Weight_down = Q2weightsGr[1]->Eval(genll.pt());
		    }
		}

	      for(size_t ivar=0; ivar<nvarsToInclude; ivar++){
		float iweight = weight*photonWeight;                                   //nominal
		if(ivar==5)                        iweight *= TotalWeight_plus;        //pu up
		if(ivar==6)                        iweight *= TotalWeight_minus;       //pu down
		if(ivar==7)                        iweight *= Q2Weight_plus;
		if(ivar==8)                        iweight *= Q2Weight_down;
		if(ivar==9)                        iweight *= PDFWeight_plus;
		if(ivar==10)                       iweight *= PDFWeight_down;
	    
		llvvJetExtCollection localSelJets;
		for(size_t ijet=0; ijet<jets.size(); ijet++){
	      
		  float rawpt=jets[ijet].pt();
		  float pt=rawpt;
		  if(ivar==1) pt=jets[ijet].jesup;
		  if(ivar==2) pt=jets[ijet].jesdown;
		  if(ivar==3) pt=jets[ijet].jerup;
		  if(ivar==4) pt=jets[ijet].jerdown;
		  if(pt<minJetPtToApply || fabs(jets[ijet].eta())>4.7) continue;
	      
		  Int_t idbits=jets[ijet].idbits;
		  bool passPFloose ( ((idbits>>0) & 0x1) );
		  //int puId((idbits>>3) & 0xf);
		  //bool passLoosePuId( ( puId >> 2) & 0x1);
		  int simplePuId( ( idbits >>7 ) & 0xf );
		  bool passLooseSimplePuId(  ( simplePuId >> 2) & 0x1);
		  if(!passPFloose || !passLooseSimplePuId) continue;

		  llvvJetExt iSelJet(jets[ijet]);
		  iSelJet *= pt/rawpt;
		  localSelJets.push_back( iSelJet );
		}
		if(localSelJets.size()<2)  continue;
		std::sort(localSelJets.begin(), localSelJets.end(),  sort_llvvObjectByPt);

		//recoil residuals uncertainty
		if( (ivar==11 || ivar==12) && recoilResidualsGr.size()==2)
		  {
		    for(size_t ijet=0; ijet<2; ijet++)
		      {
			float abseta=TMath::Abs(localSelJets[ijet].eta());
			float sfEnvelope=TMath::Max( TMath::Abs(1-recoilResidualsGr[0]->Eval(abseta)),
						     TMath::Abs(1-recoilResidualsGr[1]->Eval(abseta)) );
			localSelJets[ijet] *= (ivar==11 ? 1+sfEnvelope : 1-sfEnvelope);
		      }
		  }

		//re-assign the event category;
		std::vector<TString> locTags(1,chTags[ich]);
		for(unsigned int index=0; index<optim_Cuts2_jet_pt1.size();index++)
		  {
		    float minJetPt1=optim_Cuts2_jet_pt1[index];
		    float minJetPt2=optim_Cuts2_jet_pt2[index];

		    bool passLocalJet1Pt(localSelJets[0].pt()>minJetPt1);
		    bool passLocalJet2Pt(localSelJets[1].pt()>minJetPt2);
		    if(!passLocalJet1Pt || !passLocalJet2Pt) continue; 
		
		    LorentzVectorF vbfSyst=localSelJets[0]+localSelJets[1];
		    float mjj=vbfSyst.M();
		    float detajj=fabs(localSelJets[0].eta()-localSelJets[1].eta());
		    float spt=vbfSyst.pt()/(localSelJets[0].pt()+localSelJets[1].pt());
		    
		    TString mjjCat("");
		    std::vector<TString> localSelTags=getDijetCategories(mjj,detajj,locTags,mjjCat);
		    if(gammaWgtHandler!=0) {
		      mjjCat.ReplaceAll("mjjq100","mjjq092");
		      std::vector<Float_t> photonVars(1,selPhotons[0].pt());
		      iweight *= gammaWgtHandler->getWeightFor(photonVars,chTags[ich]+mjjCat);
		    }
		    mon.fillHisto(TString("dijet_deta_shapes")+varNames[ivar],localSelTags,index,detajj,iweight);
		    
		    //set the variables to be used in the MVA evaluation (independently of its use)
		    for(size_t mvar=0; mvar<tmvaVarNames.size(); mvar++) 
		      {
			std::string variable        = tmvaVarNames[mvar];
			if(variable=="mjj")        tmvaVars[mvar]=mjj;
			if(variable=="detajj")     tmvaVars[mvar]=detajj;
			if(variable=="spt")        tmvaVars[mvar]=spt;
		      }
		    if(tmvaReader)
		      for(size_t im=0; im<tmvaMethods.size(); im++)
			{
			  float iTmvaDiscrVal=tmvaReader->EvaluateMVA( tmvaMethods[im] );
			  mon.fillHisto(TString(tmvaMethods[im]+"_shapes")+varNames[ivar],localSelTags,index,iTmvaDiscrVal,iweight);
			}
		  }
	      }
	    }//end passZpt && passZeta
	
	    //
	    //N-1 CONTROL
	    //
	    if(           passZeta && (njets>=2))   { mon.fillHisto("zptNM1"      ,   tags, zll.pt(),     weight); }
	    if(passZpt             && (njets>=2))   { mon.fillHisto("zetaNM1"     ,   tags, zll.eta(),    weight);  mon.fillHisto("zyNM1"     , tags, zy,    weight);}
	    if(passZpt && passZeta && (njets>=2))
	      {
		mon.fillHisto  ("mjjvsdetajj",       tags, detajj, mjj, weight);
		mon.fillHisto  ("detajjvsmjj",       tags, mjj,    detajj, weight);
	      }

	  }//end passZmass
	}

  }
  
  printf("\n"); 
//  file->Close();
  
  //##############################################
  //########     SAVING HISTO TO FILE     ########
  //##############################################
  //save control plots to file
  outUrl += "/";
  outUrl += outFileUrl + ".root";
  printf("Results save in %s\n", outUrl.Data());

  //save all to the file
  TFile *ofile=TFile::Open(outUrl, "recreate");
  mon.Write();
  ofile->Close();

  //save summary tuple
  outUrl.ReplaceAll(".root","_summary.root");
  ofile=TFile::Open(outUrl,"recreate");
  summaryTuple->SetDirectory(ofile);
  summaryTuple->Write();
  ofile->Close();
  if(outTxtFile)fclose(outTxtFile);
}  





