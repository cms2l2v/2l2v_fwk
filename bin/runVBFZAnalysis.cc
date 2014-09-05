#include <iostream>
#include <boost/shared_ptr.hpp>

#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/llvv_fwk/interface/DataEventSummaryHandler.h"
#include "UserCode/llvv_fwk/interface/TMVAUtils.h"
#include "UserCode/llvv_fwk/interface/LeptonEfficiencySF.h"
#include "UserCode/llvv_fwk/interface/PDFInfo.h"
#include "UserCode/llvv_fwk/interface/MuScleFitCorrector.h"
#include "UserCode/llvv_fwk/interface/GammaWeightsHandler.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

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
bool passPhaseSpaceAcceptance(data::PhysicsObjectCollection_t &leptons, data::PhysicsObjectCollection_t &jets,bool isReco=false){
  if(leptons.size()<2 || jets.size()<2) return false;

  bool passLeptons(true);
  LorentzVector sumL(0,0,0,0);
  for(size_t i=0; i<2; i++){
    float pt=leptons[i].pt();
    float eta=leptons[i].eta();
    sumL+=leptons[i]; 
    if(isReco) passLeptons&=(fabs(eta)<2.5 && fabs(pt)>20); 
    else       passLeptons&=(fabs(eta)<2.5 && fabs(pt)>15); 
  }
  bool passDilepton( fabs(sumL.mass()-91)<30 );
  if(isReco) passDilepton=( fabs(sumL.mass()-91)<15);
  
  bool passJets(true);
  LorentzVector sumJ(0,0,0,0);
  for(size_t i=0; i<2; i++){
    float pt=jets[i].pt();
    float eta=jets[i].eta();
    sumJ+=jets[i]; 
    if(isReco) passJets&=(fabs(eta)<4.7 && fabs(pt)>30); 
    else       passJets&=(fabs(eta)<5.0 && fabs(pt)>25); 
  }
  bool passDijet( sumJ.mass()>120 );
  
  //build the final result
  //bool toReturn(passJets && passDijet);
  //bool toReturn(passJets && passDijet && passLeptons);
  bool toReturn(passJets && passDijet && passLeptons && passDilepton);
  return toReturn;
}


//
bool passGenAcceptance(data::PhysicsObjectCollection_t &gen,bool tauFilt=false){
  
  data::PhysicsObjectCollection_t jets,leptons;
  for(size_t i=0; i<gen.size(); i++){
    int status=gen[i].get("status");
    if(status!=3) continue;
    int pid=gen[i].get("id");
    if(fabs(pid)<6)                    jets.push_back(gen[i]);
    if(fabs(pid)==11 || fabs(pid)==13) leptons.push_back(gen[i]);
  }

  //select only ee or mumu events
  if(tauFilt) return (leptons.size()==2);

  //only last two quarks are outgoing
  std::reverse(jets.begin(),jets.end());
  return passPhaseSpaceAcceptance(leptons,jets);
}




//
float getAngle(LorentzVector &a, LorentzVector &b)
{
  TVector3 mom1(a.px(),a.py(),a.pz());
  TVector3 mom2(b.px(),b.py(),b.pz());
  double cosine = mom1.Dot(mom2)/(mom1.Mag()*mom2.Mag());
  return acos(cosine);
}

//
std::vector<TString> getDijetCategories(double mjj,double hardpt,std::vector<TString> &curTags, TString &mjjCat)
{
  std::vector<TString> mjjCats;
  if(mjj<250)               mjjCats.push_back("mjjq016");
  if(mjj>=250 && mjj<350)   mjjCats.push_back("mjjq033");
  if(mjj>=350 && mjj<450)   mjjCats.push_back("mjjq049");
  if(mjj>=450 && mjj<550)   mjjCats.push_back("mjjq066");
  if(mjj>=550 && mjj<750)   mjjCats.push_back("mjjq083");
  if(mjj>=750 && mjj<1000)  mjjCats.push_back("mjjq092");
  if(mjj>=1000)             mjjCats.push_back("mjjq100");
  if(mjj>=1250)             mjjCats.push_back("highmjj");
  if(mjj>=750)              mjjCats.push_back("mjjgt092");  
  if(mjj>350)
    {
      if(hardpt<50)            mjjCats.push_back("lowhardpt");
      else                     mjjCats.push_back("highhardpt");
    }
  mjjCat=mjjCats[0];
  
  //include new tags
  std::vector<TString> selTags;
  for(size_t i=0; i<curTags.size(); i++)
    {
      TString itag=curTags[i];
      selTags.push_back(itag);
      for(size_t j=0; j<mjjCats.size(); j++)
	selTags.push_back(itag+mjjCats[j]);
    }

  //all done here
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
  TGraph *recoilResidualsGr=0;
  TFile *fIn=TFile::Open("${CMSSW_BASE}/src/UserCode/llvv_fwk/data/jec/jecResponse.root");
  if(fIn){
    recoilResidualsGr=(TGraph *)fIn->Get("mumurecoilbalancevsetagt50back2backdydata2mc");
    fIn->Close();
    if(recoilResidualsGr) cout << "Read residual recoil systematics" << endl;
  }
  
  //Q^2 variations
  fIn=TFile::Open("${CMSSW_BASE}/src/UserCode/llvv_fwk/data/weights/vbfnloQ2weights.root");
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
  TString summaryTupleVarNames("ch:weight:cnorm:mjj:detajj:setajj:pt1:pt2:eta1:eta2:qg1:qg2:spt:ystar:hardpt:ncjv15:htcjv15:pt3:ystar3:mva");
  TNtuple *summaryTuple = new TNtuple("ewkzp2j","ewkzp2j",summaryTupleVarNames);
  Float_t summaryTupleVars[summaryTupleVarNames.Tokenize(":")->GetEntriesFast()];
  summaryTuple->SetDirectory(0);

  //MVA
  bool useMVA = runProcess.getParameter<bool>("useMVA");
  TMVA::Reader *tmvaReader = 0;
  std::vector<string> tmvaMethods;
  std::vector<Float_t> tmvaDiscrVals(3,0.0);
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

  TH1F* Hcutflow  = (TH1F*) mon.addHistogram(  new TH1F ("cutflow"    , "cutflow"    ,6,0,6) ) ;
  TH1 *h=mon.addHistogram( new TH1F ("eventflow", ";;Events", 5,0,5) );
  h->GetXaxis()->SetBinLabel(1,"#geq 2 leptons");
  h->GetXaxis()->SetBinLabel(2,"|M-M_{Z}|<15");
  h->GetXaxis()->SetBinLabel(3,"p_{T}^{ll}>50");
  h->GetXaxis()->SetBinLabel(4, "#eta^{ll}<1.44");
  h->GetXaxis()->SetBinLabel(5,"#geq 2 jets"); 

  h=mon.addHistogram( new TH1F ("triggerProj", ";;Events", 7,0,7) );
  h->GetXaxis()->SetBinLabel(1,"8 TeV baseline");
  h->GetXaxis()->SetBinLabel(2,"1/2 rate");
  h->GetXaxis()->SetBinLabel(3,"p^{boson}_{T}>30");
  h->GetXaxis()->SetBinLabel(4,"p^{boson}_{T}>50");
  h->GetXaxis()->SetBinLabel(5,"p_{T}^{jet}>50");
  h->GetXaxis()->SetBinLabel(6,"M_{jj}>120");
  h->GetXaxis()->SetBinLabel(7,"#splitline{p^{boson}_{T}>30}{M_{jj}>120}");

  h=mon.addHistogram( new TH1F ("psacceptance", ";;Events", 5,0,5) );
  h->GetXaxis()->SetBinLabel(1,"Gen");
  h->GetXaxis()->SetBinLabel(2,"Gen Acc");
  h->GetXaxis()->SetBinLabel(3,"Rec Acc");
  h->GetXaxis()->SetBinLabel(4,"Gen+Rec Acc");
  h->GetXaxis()->SetBinLabel(5,"~Gen Acc+Rec Acc");

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
  for(size_t ireg=0; ireg<3; ireg++){
    TString regStr("lt15collinear");
    if(ireg==1) regStr="gt50back2back";
    if(ireg==2) regStr="lt50back2back";
    mon.addHistogram( new TH1F("qgmva"+regStr, "; Quark/gluon discriminator; Jets", 100,0,1) );
    mon.addHistogram( new TH1F("ptd"+regStr, "; p_{T}D; Jets", 100,0,1) );
    mon.addHistogram( new TH1F("betastar"+regStr, "; #beta^{*}; Jets", 100,0,1) );
    mon.addHistogram( new TH1F("dr2mean"+regStr, "; <#DeltaR^{2}>; Jets", 100,0,1) );
    mon.addHistogram( new TH1F("recoilbalance"+regStr, "; p_{T}(jet)/p_{T}; Jets", 100,0,5) );
    mon.addHistogram( new TH2F("recoilbalancevseta"+regStr, "; #eta(jet); <p_{T}(jet)/p_{T}>", 50,0,5,100,0,5) );
    TH2 *idH=(TH2 *)mon.addHistogram( new TH2F("recoilbalanceid"+regStr, "; Pseudo-rapidity; ID", 50,0,5,4,0,4) );
    idH->GetYaxis()->SetBinLabel(1,"no id");  
    idH->GetYaxis()->SetBinLabel(2,"PF");
    idH->GetYaxis()->SetBinLabel(3,"PU");
    idH->GetYaxis()->SetBinLabel(4,"PF+PU");
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
  mon.addHistogram( new TH2F("njetsvsavginstlumi", ";Inst. luminosity;Jet multiplicity;Events",5,0,5,10,0,5000) );
  mon.addHistogram( new TH1F("vbfcandjeteta"     , ";Tag  jet #eta;Events x 2",                                 50,0,5) );
  mon.addHistogram( new TH1F("vbfcandjetpt"      , ";Tag jet p_{T} [GeV];Events x 2",                        50,0,500) );

  //  Double_t jetptaxis[]={30.,40.,50.,60.,70.,80.,90.,100.,150.,200.,250.,300.,400.,500.,750.,1000.};
  //Int_t nptBins=sizeof(jetptaxis)/sizeof(Double_t)-1;
  //mon.addHistogram( new TH1F("vbfcandjet1pt"     , ";Leading jet p_{T} [GeV];Events",                        nptBins,jetptaxis) );
  //mon.addHistogram( new TH1F("vbfcandjet2pt"     , ";Trailer jet p_{T} [GeV];Events",                        nptBins,jetptaxis) );

  mon.addHistogram( new TH1F("vbfcandjet1pt"     , ";Leading jet p_{T} [GeV];Events",     50,0,1000) );
  mon.addHistogram( new TH1F("vbfcandjet2pt"     , ";Trailer jet p_{T} [GeV];Events",     50,0,500) );
  mon.addHistogram( new TH1F("vbfqgmva1",             "; Quark/gluon discriminator; Events", 100,0,1) );
  mon.addHistogram( new TH1F("vbfqgmva2",             "; Quark/gluon discriminator; Events", 100,0,1) );
  mon.addHistogram( new TH1F("vbfcandjet1eta"    , ";Forward jet #eta;Events",                                 25,0,5) );
  mon.addHistogram( new TH1F("vbfcandjet2eta"    , ";Central jet #eta;Events",                                 25,0,5) );
  mon.addHistogram( new TH1F("vbfcandjetdeta"    , ";Dijet pseudo-rapidity distance (#Delta#eta);Events",   50,0,10) );  
  mon.addHistogram( new TH1F("vbfcandjetseta"    , ";#Sigma |#eta(j)|;Jets",                50,0,15) );  
  mon.addHistogram( new TH1F("vbfcandjetetaprod" , ";#eta_{j1} x #eta_{j2};Events",       50,-25,25) );
  mon.addHistogram( new TH1F("vbfhardpt"         , ";Hard process p_{T} [GeV];Events",                   25,0,250) );
  mon.addHistogram( new TH1F("vbfspt"            , ";Relative balance (#Delta^{rel}_{p_{T}});Events",                   50,0,1) );
  h=mon.addHistogram( new TH1F("vbfcjv"          , ";Central jet count;Events",                     5,0,5) );
  h->GetXaxis()->SetBinLabel(1,"=0 jets");
  h->GetXaxis()->SetBinLabel(2,"=1 jets");
  h->GetXaxis()->SetBinLabel(3,"=2 jets");
  h->GetXaxis()->SetBinLabel(4,"=3 jets");
  h->GetXaxis()->SetBinLabel(5,"#geq 4 jets");
  mon.addHistogram( (TH1F *) h->Clone("vbfcjv15") );
  mon.addHistogram( (TH1F *) h->Clone("vbfcjv20") );
  mon.addHistogram( new TH1F("vbfhtcjv15"          , ";H_{T}(p_{T}>15) [GeV];Events",10,0,250) );
  mon.addHistogram( new TH1F("vbfhtcjv20"          , ";H_{T}(p_{T}>20) [GeV];Events",10,0,250) );
  mon.addHistogram( new TH1F("vbfhtcjv"            , ";H_{T}(p_{T}>30) [GeV];Events",10,0,250) );
  mon.addHistogram( new TH1F("vbfmaxcjvjpt"        , ";Third jet p_{T} [GeV];Events", 15,0,300) );
  mon.addHistogram( new TH1F("vbfystar3"           , ";y_{j3}-(y_{j1}+y_{j2})/2;Events",10,0,5) );

  Double_t mjjaxis[32];
  mjjaxis[0]=0.01;
  for(size_t i=1; i<20; i++)  mjjaxis[i]   =50*i;        //0-1000
  for(size_t i=0; i<5; i++)   mjjaxis[20+i]=1000+100*i; //1000-1500
  for(size_t i=0; i<=5; i++)   mjjaxis[25+i]=1500+300*i; //1500-5000  
  mjjaxis[31]=5000;
  mon.addHistogram( new TH1F("vbfmjj"            , ";Dijet invariant mass [GeV];Events",               31,mjjaxis) );
  mon.addHistogram( new TH1F("vbfdphijj"         , ";Dijet azimuthal opening (#Delta#phi_{jj});Events", 20,0,3.5) );
  mon.addHistogram( new TH1F("vbfystar"          , ";y^{*}_{ll}=#eta_{ll}-(#eta_{j1}+#eta_{j2})/2;Events",       10,0,5) );
  mon.addHistogram( new TH1F("vbfpt"             , ";Dijet p_{T} [GeV];Dijets",       50,0,500) );
  mon.addHistogram( new TH1F("met"               , ";E_{T}^{miss} [GeV]; Events",        50,0,500) );
  mon.addHistogram( new TH1F("metL"              , ";Axial E_{T}^{miss} [GeV]; Events",  50,-50,200) );
  mon.addHistogram( new TH2F("vbfystarvsmjj",";Dijet invariant mass [GeV];y^{*}_{ll}=#eta_{ll}-(#eta_{j1}+#eta_{j2})/2;Events", 30,0,3000,25,0,5));

  //profiles for soft hadronic activity
  for(size_t i=0; i<2; i++)
    {
      TString softHadType(i==0 ? "soft" : "softin" );
      mon.addHistogram( new TH2F(softHadType+"jetsvsnvtx",";Number of vertices;Soft jet multiplicity;Events", 10,0,50,10,0,10));
      mon.addHistogram( new TH2F(softHadType+"htvsnvtx",  ";Number of vertices;Soft H_{T} [GeV];Events",      10,0,50,25,0,250));
      
      mon.addHistogram( new TH2F(softHadType+"jetsvsmjj", ";Dijet invariant mass [GeV];Soft jet multiplicity;Events", 30,0,3000,10,0,10));
      mon.addHistogram( new TH2F(softHadType+"htvsmjj",   ";Dijet invariant mass [GeV];Soft H_{T} [GeV];Events",      30,0,3000,25,0,250));
      
      mon.addHistogram( new TH2F(softHadType+"jetsvsdetajj",";Dijet pseudo-rapidity distance (#Delta#eta);Soft jet multiplicity;Events", 25,0,10,10,0,10));
      mon.addHistogram( new TH2F(softHadType+"htvsdetajj",  ";Dijet pseudo-rapidity distance (#Delta#eta);Soft H_{T} [GeV];Events",      25,0,10,25,0,250));
      
      mon.addHistogram( new TH2F(softHadType+"jetsvsdphijj",";Dijet azimuthal opening (#Delta#phi_{jj});Soft jet multiplicity;Events", 10,0,3.5,10,0,10));
      mon.addHistogram( new TH2F(softHadType+"htvsdphijj",  ";Dijet azimuthal opening (#Delta#phi_{jj});Soft H_{T} [GeV];Events",      10,0,3.5,25,0,250));
    }

  
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
  DataEventSummaryHandler evSummary;
  TFile *file = TFile::Open(url);
  printf("Looping on %s\n",url.Data());
  if(file==0) return -1;
  if(file->IsZombie()) return -1;
  if( !evSummary.attach( (TTree *) file->Get(baseDir+"/data") , false) ) { file->Close();  return -1; }
  const Int_t totalEntries= evSummary.getEntries();
 
  //MC normalization (to 1/pb)
  float cnorm=1.0;
  if(isMC){
    TH1F* cutflowH = (TH1F *) file->Get(baseDir+"/cutflow");
    if(cutflowH) cnorm=cutflowH->GetBinContent(1);
    printf("cnorm = %f\n",cnorm);
  }
  Hcutflow->SetBinContent(1,cnorm);
  double xsecWeight = xsec;
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
	  mon.addHistogram( new TH1F("vbfcandjetdeta"+var    , ";|#Delta #eta|;Jets",                             50,0,10) );
	  mon.addHistogram( new TH1F("vbfcandjet1eta"+var    , ";#eta;Jets",                                      50,0,5) );
	  mon.addHistogram( new TH1F("vbfcandjet1pt"+var     , ";p_{T} [GeV];Jets",                               50,0,1000) ); //nptBins,jetptaxis) );
	  mon.addHistogram( new TH1F("vbfcandjet2eta"+var    , ";#eta;Jets",                                      50,0,5) );
	  mon.addHistogram( new TH1F("vbfcandjet2pt"+var     , ";p_{T} [GeV];Jets",                               50,0,500) ); //nptBins,jetptaxis) );
	  mon.addHistogram( new TH1F("vbfspt"+var            , ";Relative balance (#Delta^{rel}_{p_{T}});Events", 50,0,1) );
	}
      cout << "Readout " << mPDFInfo->numberPDFs() << " pdf variations" << endl;
    }

  //pileup weighting
  std::vector<double> dataPileupDistributionDouble = runProcess.getParameter< std::vector<double> >("datapileup");
  std::vector<float> dataPileupDistribution; for(unsigned int i=0;i<dataPileupDistributionDouble.size();i++){dataPileupDistribution.push_back(dataPileupDistributionDouble[i]);}
  std::vector<float> mcPileupDistribution;
  if(isMC){
    TString puDist(baseDir+"/pileup");
    TH1F* histo = (TH1F *) file->Get(puDist);
    if(!histo) std::cout<<"pileup histogram is null!!!\n";
    for(int i=1;i<=histo->GetNbinsX();i++){mcPileupDistribution.push_back(histo->GetBinContent(i));}
    delete histo;
  }
  while(mcPileupDistribution.size()<dataPileupDistribution.size())  mcPileupDistribution.push_back(0.0);
  while(mcPileupDistribution.size()>dataPileupDistribution.size())dataPileupDistribution.push_back(0.0);

  gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE
  edm::LumiReWeighting *LumiWeights= isMC ? new edm::LumiReWeighting(mcPileupDistribution,dataPileupDistribution): 0;
  utils::cmssw::PuShifter_t PuShifters;
  if(isMC) { PuShifters=utils::cmssw::getPUshifters(dataPileupDistribution,0.05); }

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
      evSummary.getEntry(iev);
      DataEventSummary &ev = evSummary.getEvent();
      if(!isMC && duplicatesChecker.isDuplicate( ev.run, ev.lumi, ev.event) ) { nDuplicates++; continue; }

      if(isV0JetsMC){
	mon.fillHisto("nup","",ev.nup,1);
	if(ev.nup>5) continue;
	mon.fillHisto("nupfilt","",ev.nup,1);
      }


      data::PhysicsObjectCollection_t photons=evSummary.getPhysicsObject(DataEventSummaryHandler::PHOTONS);
      data::PhysicsObjectCollection_t leptons=evSummary.getPhysicsObject(DataEventSummaryHandler::LEPTONS);
      data::PhysicsObjectCollection_t jets=evSummary.getPhysicsObject(DataEventSummaryHandler::JETS);
      data::PhysicsObjectCollection_t recoMet=evSummary.getPhysicsObject(DataEventSummaryHandler::MET);
      data::PhysicsObjectCollection_t gen=evSummary.getPhysicsObject(DataEventSummaryHandler::GENPARTICLES);      

      //require compatibilitiy of the event with the PD
      bool eeTrigger          = ev.t_bits[0];
      bool muTrigger          = ev.t_bits[6];
      bool mumuTrigger        = ev.t_bits[2] || ev.t_bits[3] || muTrigger;
      if(filterOnlyEE)   { mumuTrigger=false; }
      if(filterOnlyMUMU) { eeTrigger=false;   }
      if(isSingleMuPD)   { eeTrigger=false; if( mumuTrigger || !muTrigger ) mumuTrigger= false;  }
      
      bool hasPhotonTrigger(false);
      float triggerPrescale(1.0),triggerThreshold(0);
      TString phoTrigCat("");
      if(runPhotonSelection)
	{
	  eeTrigger=false; mumuTrigger=false;
	  for(size_t itrig=10; itrig>=7; itrig--)
	    {
	      if(!ev.t_bits[itrig]) continue;
	      hasPhotonTrigger=true;
	      triggerPrescale=ev.t_prescale[itrig];
	      if(itrig==10) {triggerThreshold=90; phoTrigCat="Photon90";}
	      if(itrig==9)  {triggerThreshold=75; phoTrigCat="Photon75";}
	      if(itrig==8)  {triggerThreshold=50; phoTrigCat="Photon50";}
	      if(itrig==7)  {triggerThreshold=36; phoTrigCat="Photon36";}
	      break;
	    }
	}

      //pileup weight
      float weight = 1.0;
      double TotalWeight_plus = 1.0;
      double TotalWeight_minus = 1.0;
      float puWeight(1.0);
      if(isMC){
        puWeight          = LumiWeights->weight(ev.ngenITpu);
	weight            = puWeight;
        TotalWeight_plus  = PuShifters[utils::cmssw::PUUP]->Eval(ev.ngenITpu);
        TotalWeight_minus = PuShifters[utils::cmssw::PUDOWN]->Eval(ev.ngenITpu);
      }
      Hcutflow->Fill(1,1);
      Hcutflow->Fill(2,weight);
      Hcutflow->Fill(3,weight*TotalWeight_minus);
      Hcutflow->Fill(4,weight*TotalWeight_plus);
      Hcutflow->Fill(5,weight);

      //
      //
      // BELOW FOLLOWS THE ANALYSIS OF THE MAIN SELECTION WITH N-1 PLOTS
      //
      //

      //
      // photon selection
      //
      data::PhysicsObjectCollection_t selPhotons;
      if(runPhotonSelection)
	{
	  //filter out number of prompt photons to avoid double counting
	  int ngenpho(0);
	  for(size_t igen=0; igen<gen.size(); igen++)
	    {
	      if(gen[igen].get("id")!=22 || gen[igen].get("status")!=1) continue;
	      float lxy=gen[igen].getVal("lxy");
	      if(lxy>0) continue;
	      ngenpho++;
	    }
	  if(mctruthmode==111 && ngenpho>0) continue;
	  //if(mctruthmode==22 && ngenpho==0) continue;

	  //select the photons
	  for(size_t ipho=0; ipho<photons.size(); ipho++)
	    {
	      double pt=photons[ipho].pt();
	      double eta=photons[ipho].getVal("sceta");

	      //if systematics are active loosen the selection to the medium working point
	      Int_t idbits( photons[ipho].get("id") );
	      bool hasLoosePhotonId( (idbits >> 0 ) & 0x1 );
	      bool hasTightPhotonId( (idbits >> 2 ) & 0x1 );
	      double gIso    = photons[ipho].getVal("gIso03");
	      double gArea   = utils::cmssw::getEffectiveArea(22,eta,3,"gIso");	      
	      double chIso   = photons[ipho].getVal("chIso03");
	      double chArea  = utils::cmssw::getEffectiveArea(22,eta,3,"chIso");
	      double nhIso   = photons[ipho].getVal("nhIso03");
	      double nhArea  = utils::cmssw::getEffectiveArea(22,eta,3,"nhIso");
	      
	      //select the photon
	      if(pt<triggerThreshold || fabs(eta)>1.4442 ) continue;
	      bool passId(true);
	      if( photons[ipho].getVal("r9")<0.9 ) passId=false;
	      if(!passId) continue;
	      bool passIso(true);
	      if(runLoosePhotonSelection){
		if(!hasLoosePhotonId) passId=false;
		passIso &= (TMath::Max(chIso-chArea*ev.rho,0.0) < 2.6); 
		passIso &= (TMath::Max(nhIso-nhArea*ev.rho,0.0) < 3.5+0.04*pt); 
		passIso &= (TMath::Max(gIso-gArea*ev.rho,  0.0) < 1.3+0.005*pt); 
	      }
	      else{
		if(!hasTightPhotonId) passId=false;
		passIso &= (TMath::Max(chIso-chArea*ev.rho,0.0) < 0.7); 
		passIso &= (TMath::Max(nhIso-nhArea*ev.rho,0.0) < 0.4+0.04*pt); 
		passIso &= (TMath::Max(gIso-gArea*ev.rho,  0.0) < 0.5+0.005*pt); 
	      }

	      if(!passIso) continue; 
	      selPhotons.push_back(photons[ipho]);
	    }
	}


      //
      // LEPTON ANALYSIS
      //
      data::PhysicsObjectCollection_t selLeptons;
      for(size_t ilep=0; ilep<leptons.size(); ilep++)
	{
	  bool passKin(true),passId(true),passIso(true);
	  int lid=leptons[ilep].get("id");

	  //apply muon corrections
	  if(abs(lid)==13 && muCor){
	    TLorentzVector p4(leptons[ilep].px(),leptons[ilep].py(),leptons[ilep].pz(),leptons[ilep].energy());
	    muCor->applyPtCorrection(p4 , lid<0 ? -1 :1 );
	    if(isMC) muCor->applyPtSmearing(p4, lid<0 ? -1 : 1, false);
	    leptons[ilep].SetPxPyPzE(p4.Px(),p4.Py(),p4.Pz(),p4.E());
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
	  float leta = lid==11 ? leptons[ilep].getVal("sceta") : leptons[ilep].eta();
	  if(leptons[ilep].pt()<20)                   passKin=false;
	  if(leta> (lid==11 ? 2.5 : 2.4) )            passKin=false;
	  if(lid==11 && (leta>1.4442 && leta<1.5660)) passKin=false;

	  //id
	  Int_t idbits = leptons[ilep].get("idbits");
	  if(lid==11){
	    if(leptons[ilep].getFlag("isconv"))              passId=false;
	    bool isLoose = ((idbits >> 4) & 0x1);
	    if(!isLoose)                                   passId=false;
 	  }
	  else{
	    bool isLoose    = ((idbits >> 8) & 0x1);
	    if(!isLoose)                                   passId=false;
	  }

	  //isolation
	  Float_t gIso    = leptons[ilep].getVal(lid==11 ? "gIso03"    : "gIso04");
	  Float_t chIso   = leptons[ilep].getVal(lid==11 ? "chIso03"   : "chIso04");
	  Float_t puchIso = leptons[ilep].getVal(lid==11 ? "puchIso03" : "puchIso04");  
	  Float_t nhIso   = leptons[ilep].getVal(lid==11 ? "nhIso03"   : "nhIso04");
	  float relIso= lid==11 ?
	    (TMath::Max(nhIso+gIso-ev.rho*utils::cmssw::getEffectiveArea(11,leptons[ilep].getVal("sceta")),Float_t(0.))+chIso)/leptons[ilep].pt() :
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

      std::sort(selLeptons.begin(), selLeptons.end(), data::PhysicsObject_t::sortByPt);


      
      //
      //JET/MET ANALYSIS
      //
      //add scale/resolution uncertainties and propagate to the MET
      utils::cmssw::updateJEC(jets,jesCor,totalJESUnc,ev.rho,ev.nvtx,isMC);
      std::vector<LorentzVector> met=utils::cmssw::getMETvariations(recoMet[0],jets,selLeptons,isMC);

      //select the jets
      data::PhysicsObjectCollection_t selJets, selJetsNoId;
      int njets(0);
      for(size_t ijet=0; ijet<jets.size(); ijet++) 
	{
	  if(jets[ijet].pt()<15 || fabs(jets[ijet].eta())>4.7 ) continue;

	  //mc truth for this jet
	  const data::PhysicsObject_t &genJet=jets[ijet].getObject("genJet");
	  TString jetType( genJet.pt()>0 ? "truejetsid" : "pujetsid" );
	  
	  //cross-clean with selected leptons and photons
	  double minDRlj(9999.),minDRlg(9999.);
          for(size_t ilep=0; ilep<selLeptons.size(); ilep++)
            minDRlj = TMath::Min( minDRlj, deltaR(jets[ijet],selLeptons[ilep]) );
	  for(size_t ipho=0; ipho<selPhotons.size(); ipho++)
	    minDRlg = TMath::Min( minDRlg, deltaR(jets[ijet],selPhotons[ipho]) );
	  if(minDRlj<0.4 || minDRlg<0.4) continue;
	  
	  //jet id
	  // float pumva=jets[ijet].getVal("puMVA");
	  Int_t idbits=jets[ijet].get("idbits");
	  bool passPFloose( ((idbits>>0) & 0x1));
	  //int puId( ( idbits >>3 ) & 0xf );
	  //bool passLoosePuId( ( puId >> 2) & 0x1);
	  int simplePuId( ( idbits >>7 ) & 0xf );
	  bool passLooseSimplePuId(  ( simplePuId >> 2) & 0x1);
	  if(jets[ijet].pt()>30)
	    {
	      mon.fillHisto(jetType,"",fabs(jets[ijet].eta()),0);
	      if(passPFloose)                        mon.fillHisto(jetType,"",fabs(jets[ijet].eta()),1);
	      if(passLooseSimplePuId)                mon.fillHisto(jetType,"",fabs(jets[ijet].eta()),2);
	      if(passPFloose && passLooseSimplePuId) mon.fillHisto(jetType,"",fabs(jets[ijet].eta()),3);
	    }
	  
	  selJetsNoId.push_back(jets[ijet]);
	  if(passPFloose && passLooseSimplePuId){
	    selJets.push_back(jets[ijet]);
	    if(jets[ijet].pt()>minJetPtToApply) njets++;
	  }
	}
      std::sort(selJets.begin(), selJets.end(), data::PhysicsObject_t::sortByPt);


      //check acceptance
      bool isGenEEorMuMu=passGenAcceptance(gen,true);
      bool isInGenAcceptance=passGenAcceptance(gen);
      bool isInRecoAcceptance=(passPhaseSpaceAcceptance(selLeptons,selJets,true) && (eeTrigger||mumuTrigger));
      if(isGenEEorMuMu){
	float weight=1.0;
	mon.fillHisto("psacceptance","",0,weight);
	if(isInGenAcceptance) mon.fillHisto("psacceptance","",1,weight);
	if(isInRecoAcceptance) mon.fillHisto("psacceptance","",2,weight);
	if(isInGenAcceptance && isInRecoAcceptance) mon.fillHisto("psacceptance","",3,weight);
	if(!isInGenAcceptance && isInRecoAcceptance) mon.fillHisto("psacceptance","",4,weight);
      }

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
	      dilId *= selLeptons[ilep].get("id");
	      int id(abs(selLeptons[ilep].get("id")));
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
      LorentzVector leadingLep(runPhotonSelection ? selPhotons[0] : selLeptons[0].pt()>selLeptons[1].pt() ? selLeptons[0]: selLeptons[1]);
      LorentzVector trailerLep(runPhotonSelection ? selPhotons[0] : selLeptons[0].pt()>selLeptons[1].pt() ? selLeptons[1]: selLeptons[0]);
      LorentzVector zll(runPhotonSelection ? selPhotons[0] : leadingLep+trailerLep);
      float zy(zll.Rapidity());
      bool passZmass(runPhotonSelection || (fabs(zll.mass()-91)<15));
      bool passZpt(zll.pt()>50);
      bool passZeta(fabs(zy)<1.4442);

      //generator level
      LorentzVector genll;
      for(size_t ig=0; ig<gen.size(); ig++)
	{
	  int pid=abs(gen[ig].get("id"));
	  if(pid!=11 && pid!=13) continue;
	  int status=abs(gen[ig].get("status"));
	  if(status!=3) continue;
	  genll += gen[ig];
	}
      
      //analyze dijets and activity in the dijet rapidity distance
      float maxPt(0), minPt(0), maxEta(0), minEta(0), maxAbsEta(0), minAbsEta(0);
      float dphijj(0), detajj(0), setajj(0), etaprod(0), ystar(0), mjj(0),ptjj(0),spt(0);
      float hardpt(0);
      int ncjv(0), ncjv15(0),ncjv20(0), htcjv(0), htcjv15(0),htcjv20(0);
      int nj15(0),htj15(0);
      float pt3(0), ystar3(0);
      float ptmiss(0),metL(0);
      if(njets>=2)
	{
	  LorentzVector jet1=selJets[0];
	  LorentzVector jet2=selJets[1];
	  maxPt=jet1.pt();
	  minPt=jet2.pt();
	  maxAbsEta=max(fabs(jet1.eta()),fabs(jet2.eta()));
	  minAbsEta=min(fabs(jet1.eta()),fabs(jet2.eta()));
	  maxEta=max(jet1.eta(),jet2.eta());
	  minEta=min(jet1.eta(),jet2.eta());
	  detajj=fabs(maxEta-minEta);
	  setajj=fabs(maxEta)+fabs(minEta);
	  dphijj=deltaPhi(jet1.phi(),jet2.phi());
	  etaprod=jet1.eta()*jet2.eta(); 
	  ystar=zll.Rapidity()-0.5*(jet1.Rapidity()+jet2.Rapidity());
	  if(detajj<2 && mjj>1000)
	    cout << jet1.pt() << "," << jet1.eta() << " | " << jet2.pt() << "," << jet2.eta() << " | " << zll.pt() << "," << zll.eta() << endl;
	    
	  LorentzVector vbfSyst=jet1+jet2;
	  mjj=vbfSyst.mass();
	  ptjj=vbfSyst.pt();
	  spt=vbfSyst.pt()/(jet1.pt()+jet2.pt());

	  LorentzVector hardSyst=vbfSyst+zll; 
	  hardpt=hardSyst.pt();
	  
	  TVector2 boson2(zll.px(),zll.py());
	  TVector2 met2(met[0].px(),met[0].py());
	  ptmiss=met[0].pt();
	  metL=boson2*met2; metL /= -zll.pt();
	  
	  //visible system rest frame
	  //BetaVector vbfBoost         = hardSyst.BoostToCM();
	  //LorentzVector jet1_cm    = ROOT::Math::VectorUtil::boost( jet1, vbfBoost );
	  //LorentzVector jet2_cm    = ROOT::Math::VectorUtil::boost( jet2, vbfBoost );
	  //LorentzVector jj_cm      = ROOT::Math::VectorUtil::boost( vbfSyst, vbfBoost );
	  //LorentzVector z_cm       = ROOT::Math::VectorUtil::boost( zll,  vbfBoost );
	  
	  for(size_t iotherjet=2; iotherjet<selJets.size(); iotherjet++){
	    float ipt=selJets[iotherjet].pt();
	    nj15++;
	    htj15+=ipt; 
	    bool isInRapGap(selJets[iotherjet].eta()>minEta && selJets[iotherjet].eta()<maxEta);
	    if(isInRapGap)
	      {
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
	}
      
      { 
	//trigger studies
	bool passHalfRate(false);
	bool passZpt30( zll.pt()>30 );
	bool passZpt50( zll.pt()>50 );
	bool passJet50( maxPt>50 );
	bool passMjj( mjj>120 );
	bool passBase( abs(dilId)==22 ? passZpt50 : true );
	
	if(abs(dilId)==11*11) passHalfRate=(leadingLep.pt()>30 && trailerLep.pt()>27 );
	if(abs(dilId)==13*13) passHalfRate=(leadingLep.pt()>23 && trailerLep.pt()>10 );
	if(abs(dilId)==22)    passHalfRate=(zll.pt()>155);
	
	if(passBase)
	  {
	    mon.fillHisto("triggerProj",chTags,0,1);
	    if(passHalfRate)              mon.fillHisto("triggerProj",chTags,1,weight);
	    if(passZpt30)                 mon.fillHisto("triggerProj",chTags,2,weight);
	    if(passZpt50)                 mon.fillHisto("triggerProj",chTags,3,weight);
	    if(passJet50)                 mon.fillHisto("triggerProj",chTags,4,weight);
	    if(passMjj)                   mon.fillHisto("triggerProj",chTags,5,weight);
	    if(passZpt30 && passMjj)      mon.fillHisto("triggerProj",chTags,6,weight);
	  }
      }


      
      //set the variables to be used in the MVA evaluation (independently of its use)
      if(njets>1)
	{
	  for(size_t ivar=0; ivar<tmvaVarNames.size(); ivar++) 
	    {
	      std::string variable     = tmvaVarNames[ivar];
	      if(variable=="mjj")        tmvaVars[ivar]=mjj;
	      if(variable.find("detajj")!=string::npos)     tmvaVars[ivar]=detajj;
	      if(variable=="spt")        tmvaVars[ivar]=spt;
	      if(variable=="setajj")     tmvaVars[ivar]=setajj;
	      if(variable=="pt1")        tmvaVars[ivar]=selJets[0].pt();
	      if(variable=="pt2")        tmvaVars[ivar]=selJets[1].pt();
	      if(variable=="eta1")       tmvaVars[ivar]=selJets[0].eta();
	      if(variable=="eta2")       tmvaVars[ivar]=selJets[1].eta();
	      if(variable=="qg1")        tmvaVars[ivar]=selJets[0].getVal("qgMVA");
	      if(variable=="qg2")        tmvaVars[ivar]=selJets[1].getVal("qgMVA");
	      if(variable=="hardpt")     tmvaVars[ivar]=hardpt;
	    }
	  if(tmvaReader) {
	    for(size_t im=0; im<tmvaMethods.size(); im++)
	      tmvaDiscrVals[im]=tmvaReader->EvaluateMVA( tmvaMethods[im] );
	  }
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
	    mon.fillHisto("nvtx"     ,   tags, ev.nvtx,      weight);
	    mon.fillHisto("nvtxraw"  ,   tags, ev.nvtx,      weight/puWeight);
	    mon.fillHisto("rho"      ,   tags, ev.rho,       weight);
	
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
		bool isPUjet( selJetsNoId[0].getObject("genJet").pt()==0 ); 
		
		//kinematics
		float balance = selJetsNoId[0].pt()/zll.pt();
		float dphi    = fabs( deltaPhi(selJetsNoId[0].phi(),zll.phi()) );
		mon.fillHisto("ptllvsdphi",                                    tags,zll.pt(),dphi,weight);
		mon.fillHisto(TString("ptllvsdphi")+(isPUjet ? "pu" : "true"), tags,zll.pt(),dphi,weight);
	
		TString regStr("");
		if(dphi<1   && zll.pt()<15) regStr="lt15collinear";
		if(dphi>2.7 && zll.pt()>50) regStr="gt50back2back";
		if(dphi>2.7 && zll.pt()<=50) regStr="lt50back2back";
		if(regStr!="")
		  {
		    //ids
		    Int_t idbits=selJetsNoId[0].get("idbits");
		    bool passPFloose ( ((idbits>>0) & 0x1) );
		    //int puId((idbits>>3) & 0xf);
		    //bool passLoosePuId( ( puId >> 2) & 0x1);
		    int simplePuId( ( idbits >>7 ) & 0xf );
		    bool passLooseSimplePuId(  ( simplePuId >> 2) & 0x1);

		    mon.fillHisto("recoilbalance"+regStr,      tags,balance, weight);
		    mon.fillHisto("recoilbalancevseta"+regStr, tags,fabs(selJetsNoId[0].eta()), balance, weight);
		    mon.fillHisto("recoilbalanceid"+regStr,    tags,fabs(selJetsNoId[0].eta()),0, weight);
		    if(passPFloose)                         mon.fillHisto("recoilbalanceid"+regStr,tags,fabs(selJetsNoId[0].eta()),1, weight);
		    if(passLooseSimplePuId)                 mon.fillHisto("recoilbalanceid"+regStr,tags,fabs(selJetsNoId[0].eta()),2, weight);
		    if(passPFloose && passLooseSimplePuId) {
		      mon.fillHisto("recoilbalanceid"+regStr,tags,fabs(selJetsNoId[0].eta()),      3,    weight);
		      mon.fillHisto("qgmva"+regStr,         tags,selJetsNoId[0].getVal("qgMVA"),         weight);
		      mon.fillHisto("ptd"+regStr,           tags,selJetsNoId[0].getVal("ptD"),           weight);
		      mon.fillHisto("betastar"+regStr,      tags,selJetsNoId[0].getVal("betaStar"),      weight);
		      mon.fillHisto("dr2mean"+regStr,       tags,selJetsNoId[0].getVal("dR2Mean"),       weight);
		    }
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
		      if(ijesvar==1) { pf="_jesup";   pt=selJets[ijet].getVal("jesup"); }
		      if(ijesvar==2) { pf="_jesdown"; pt=selJets[ijet].getVal("jesdown"); }
		      if(pt>minJetPtToApply){
			njets_ivar++;
			if(ijesvar==0){
			  mon.fillHisto("jetptvseta",  tags,pt,fabs(selJets[ijet].eta()),weight);
			  mon.fillHisto("jetgt3pt",    tags,pt,weight);
			  mon.fillHisto("jetgt3nhf",   tags,selJets[ijet].getVal("neutHadFrac"),weight);
			  mon.fillHisto("jetgt3nemf",  tags,selJets[ijet].getVal("neutEmFrac"), weight);
			  mon.fillHisto("jetgt3chf",   tags,selJets[ijet].getVal("chHadFrac"),  weight);
			  mon.fillHisto("jetgt3ptrms", tags,selJets[ijet].getVal("ptRMS"),      weight);
			}
			mon.fillHisto("jetpt"+pf,  tags, pt, weight);
			mon.fillHisto("jeteta"+pf, tags, fabs(selJets[ijet].eta()), weight);
		      }
		    }
		  mon.fillHisto("njets"+pf,tags, njets_ivar, weight);
		}
		
	      //signal region
	      float photonWeight(1.0);
	      if(njets>1)
		{
		  TString mjjCat("");
		  std::vector<TString> selTags;
		  selTags = getDijetCategories(mjj,hardpt,tags,mjjCat);
		  if(phoTrigCat!="") selTags.push_back(phoTrigCat);
	  
		  //re-weight for photons if needed
		  if(gammaWgtHandler!=0) {
		    std::vector<Float_t> gammaVars(1,selPhotons[0].pt());
		    photonWeight = gammaWgtHandler->getWeightFor(gammaVars,chTags[ich]+mjjCat);
		  }
		  float catWeight=weight*photonWeight;

		  //veto events with very large weights in simulation
		  if(isMC && catWeight>5) catWeight=0;

		  //save for further analysis
		  Int_t summaryDilId=dilId;
		  if(gammaWgtHandler && chTags[ich].Contains("ee")   ) summaryDilId=-11*11;
		  if(gammaWgtHandler && chTags[ich].Contains("mumu") ) summaryDilId=-13*13;
		  if(mjj>200 && (abs(summaryDilId)==11*11 || abs(summaryDilId)==13*13) ) {
		    //  if(mjj>2000 && ev.nvtx<18 && njets==2){
		    // 	  fprintf(outTxtFile,"--------- CANDIDATE EVENT ---------\n");
		    // 	  fprintf(outTxtFile,"%d:%d:%d    mjj=%f  ystar=%f spt=%f\n",ev.run,ev.lumi,ev.event,mjj,ystar,spt);
		    // 	  fprintf(outTxtFile,"j1 (%f,%f,%f)\n",selJets[0].pt(),selJets[0].eta(),selJets[0].phi());
		    // 	  fprintf(outTxtFile,"j2 (%f,%f,%f)\n",selJets[1].pt(),selJets[1].eta(),selJets[1].phi());
		    // 	  fprintf(outTxtFile,"z (%f,%f,%f) m=%f\n",zll.pt(),zll.eta(),zll.phi(),zll.mass());
		    //	}
		    summaryTupleVars[0]=summaryDilId;
		    summaryTupleVars[1]=catWeight*xsecWeight;    
		    summaryTupleVars[2]=cnorm;
		    summaryTupleVars[3]=mjj;     
		    summaryTupleVars[4]=fabs(detajj);               
		    summaryTupleVars[5]=fabs(setajj);
		    summaryTupleVars[6]=selJets[0].pt();
		    summaryTupleVars[7]=selJets[1].pt();
		    summaryTupleVars[8]=fabs(selJets[0].eta());
		    summaryTupleVars[9]=fabs(selJets[1].eta());
		    summaryTupleVars[10]=selJets[0].getVal("qgMVA");
		    summaryTupleVars[11]=selJets[1].getVal("qgMVA");
		    summaryTupleVars[12]=spt;
		    summaryTupleVars[13]=ystar;
		    summaryTupleVars[14]=hardpt;
		    summaryTupleVars[15]=ncjv15;
		    summaryTupleVars[16]=htcjv15;  
		    summaryTupleVars[17]=pt3; 
		    summaryTupleVars[18]=ystar3;
		    if(tmvaMethods.size() ) summaryTupleVars[19]=tmvaDiscrVals[ tmvaMethods.size()-1 ];

		    summaryTuple->Fill(summaryTupleVars);
		    for(size_t im=0; im<tmvaMethods.size(); im++) mon.fillHisto(tmvaMethods[im], selTags, tmvaDiscrVals[im], catWeight);
		  } 
	      
		  mon.fillHisto("qt",                 selTags, zll.pt(), catWeight,true);      
		  mon.fillHisto("rapidity"     ,      selTags, fabs(zy),    catWeight);
		  mon.fillHisto("njetsvsavginstlumi", selTags, njets,ev.instLumi,catWeight);
		  mon.fillHisto("vbfcandjetpt",       selTags, maxPt,catWeight);
		  mon.fillHisto("vbfcandjetpt",       selTags, minPt,catWeight);
		  mon.fillHisto("vbfcandjet1pt",      selTags, maxPt,catWeight);//,true);
		  mon.fillHisto("vbfcandjet2pt",      selTags, minPt,catWeight);//,true);
		  mon.fillHisto("vbfcandjet1eta",     selTags, maxAbsEta, catWeight);
		  mon.fillHisto("vbfcandjet2eta",     selTags, minAbsEta, catWeight);
		  mon.fillHisto("vbfqgmva1",          selTags, selJets[0].getVal("qgMVA"), catWeight);
		  mon.fillHisto("vbfqgmva2",          selTags, selJets[1].getVal("qgMVA"), catWeight);
		  mon.fillHisto("vbfcandjeteta",      selTags, maxAbsEta, catWeight);
		  mon.fillHisto("vbfcandjeteta",      selTags, minAbsEta, catWeight);
		  mon.fillHisto("vbfcandjetdeta",     selTags, detajj,catWeight);
		  mon.fillHisto("vbfcandjetseta",     selTags, setajj,catWeight);
		  mon.fillHisto("vbfcandjetetaprod",  selTags, etaprod,catWeight);
		  mon.fillHisto("vbfmjj",             selTags, mjj,catWeight,true);
		  mon.fillHisto("vbfhardpt",          selTags, hardpt,catWeight);
		  mon.fillHisto("vbfspt",             selTags, spt,catWeight);
		  mon.fillHisto("vbfdphijj",          selTags, fabs(dphijj),catWeight);
		  mon.fillHisto("vbfystar",           selTags, fabs(ystar),catWeight);
		  mon.fillHisto("vbfystarvsmjj",      selTags, mjj, fabs(ystar),catWeight);
		  mon.fillHisto("vbfpt",              selTags, ptjj,catWeight);
		  mon.fillHisto("met",                selTags, ptmiss,catWeight);
		  mon.fillHisto("metL",               selTags, metL,catWeight);
		  if(ncjv15){
		    mon.fillHisto("vbfmaxcjvjpt",        selTags, pt3,catWeight);
		    mon.fillHisto("vbfystar3",           selTags, fabs(ystar3),catWeight);
		  }
		  mon.fillHisto("vbfcjv",                selTags, ncjv,catWeight);
		  mon.fillHisto("vbfcjv15",              selTags, ncjv15,catWeight);
		  mon.fillHisto("vbfcjv20",              selTags, ncjv20,catWeight);
		  if(ncjv)   mon.fillHisto("vbfhtcjv",   selTags, htcjv,catWeight);
		  if(ncjv15) mon.fillHisto("vbfhtcjv15", selTags, htcjv15,catWeight);
		  if(ncjv20) mon.fillHisto("vbfhtcjv20", selTags, htcjv20,catWeight);
	      
		  
		  //profile soft hadronic activity
		  for(int isoft=0; isoft<2; isoft++)
		    {
		      int nsoftjets(isoft==0 ? nj15 : ncjv15);
		      int softht(isoft==0 ? htj15 : htcjv15);
		      TString softHadType(isoft==0 ? "soft" : "softin" );
		      
		      mon.fillHisto(softHadType+"jetsvsnvtx",    selTags, ev.nvtx,      nsoftjets, catWeight);
		      mon.fillHisto(softHadType+"jetsvsmjj",     selTags, mjj,          nsoftjets, catWeight);
		      mon.fillHisto(softHadType+"jetsvsdetajj",  selTags, fabs(detajj), nsoftjets, catWeight);
		      mon.fillHisto(softHadType+"jetsvsdphijj",  selTags, fabs(dphijj), nsoftjets, catWeight);
		      if(nsoftjets>0)
			{
			  mon.fillHisto(softHadType+"htvsnvtx",    selTags, ev.nvtx,       softht, catWeight);
			  mon.fillHisto(softHadType+"htvsmjj",     selTags, mjj,           softht, catWeight);
			  mon.fillHisto(softHadType+"htvsdetajj",  selTags, detajj,        softht, catWeight);
			  mon.fillHisto(softHadType+"htvsdphijj",  selTags, fabs(dphijj),  softht, catWeight);
			}
		    }
		  
		  
		  //stability of pdf variations 
		  if(isMC && mPDFInfo)
		    {
		      std::vector<float> wgts=mPDFInfo->getWeights(iev);
		      for(size_t ipw=0; ipw<wgts.size(); ipw++) 
			{
			  TString var("_"); var+=ipw;
			  mon.fillHisto("vbfcandjetdeta"+var,     selTags, detajj,    catWeight*wgts[ipw]);
			  mon.fillHisto("vbfcandjet1pt"+var,      selTags, maxPt,     catWeight*wgts[ipw]);//,true);
			  mon.fillHisto("vbfcandjet2pt"+var,      selTags, minPt,     catWeight*wgts[ipw]);//,true);
			  mon.fillHisto("vbfcandjet1eta"+var,     selTags, maxAbsEta, catWeight*wgts[ipw]);
			  mon.fillHisto("vbfcandjet2eta"+var,     selTags, minAbsEta, catWeight*wgts[ipw]);
			  mon.fillHisto("vbfspt"+var,            selTags, spt,       catWeight*wgts[ipw]);
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
		float iweight = weight;                                                //nominal
		if(ivar==5)                        iweight *= TotalWeight_plus;        //pu up
		if(ivar==6)                        iweight *= TotalWeight_minus;       //pu down
		if(ivar==7)                        iweight *= Q2Weight_plus;
		if(ivar==8)                        iweight *= Q2Weight_down;
		if(ivar==9)                        iweight *= PDFWeight_plus;
		if(ivar==10)                       iweight *= PDFWeight_down;
	    
		data::PhysicsObjectCollection_t localSelJets;
		for(size_t ijet=0; ijet<jets.size(); ijet++){
	      
		  float rawpt=jets[ijet].pt();
		  float pt=rawpt;
		  if(ivar==1) pt=jets[ijet].getVal("jesup");
		  if(ivar==2) pt=jets[ijet].getVal("jesdown");
		  if(ivar==3) pt=jets[ijet].getVal("jerup");
		  if(ivar==4) pt=jets[ijet].getVal("jerdown");
		  if(pt<minJetPtToApply || fabs(jets[ijet].eta())>4.7) continue;
	      
		  Int_t idbits=jets[ijet].get("idbits");
		  bool passPFloose ( ((idbits>>0) & 0x1) );
		  //int puId((idbits>>3) & 0xf);
		  //bool passLoosePuId( ( puId >> 2) & 0x1);
		  int simplePuId( ( idbits >>7 ) & 0xf );
		  bool passLooseSimplePuId(  ( simplePuId >> 2) & 0x1);
		  if(!passPFloose || !passLooseSimplePuId) continue;

		  data::PhysicsObject_t iSelJet(jets[ijet]);
		  iSelJet *= pt/rawpt;
		  localSelJets.push_back( iSelJet );
		}
		if(localSelJets.size()<2)  continue;
		std::sort(localSelJets.begin(), localSelJets.end(),  data::PhysicsObject_t::sortByPt);
		
		//recoil residuals uncertainty
		if( (ivar==11 || ivar==12) && recoilResidualsGr)
		  {
		    for(size_t ijet=0; ijet<2; ijet++)
		      {
			float abseta=TMath::Abs(localSelJets[ijet].eta());
			float sfEnvelope=TMath::Abs(1-recoilResidualsGr->Eval(abseta));
			localSelJets[ijet] *= (ivar==11 ? 1+sfEnvelope : 1-sfEnvelope);
		      }
		  }
		
		//recompute the discriminator variables
		LorentzVector vbfSyst=localSelJets[0]+localSelJets[1];
		float mjj=vbfSyst.M();
		float detajj=fabs(localSelJets[0].eta()-localSelJets[1].eta());
		float setajj=fabs(localSelJets[0].eta())+fabs(localSelJets[1].eta());
		float spt=vbfSyst.pt()/(localSelJets[0].pt()+localSelJets[1].pt());
		std::vector<Float_t> localTmvaDiscrVals(tmvaDiscrVals.size(),0.0);
		for(size_t iMvaVar=0; iMvaVar<tmvaVarNames.size(); iMvaVar++) 
		  {
		    std::string variable     = tmvaVarNames[iMvaVar];
		    if(variable=="mjj")        tmvaVars[iMvaVar]=mjj;
		    if(variable.find("detajj")!=string::npos)     tmvaVars[iMvaVar]=detajj;
		    if(variable=="spt")        tmvaVars[iMvaVar]=spt;
		    if(variable=="setajj")     tmvaVars[iMvaVar]=setajj;
		    if(variable=="pt1")        tmvaVars[iMvaVar]=localSelJets[0].pt();
		    if(variable=="pt2")        tmvaVars[iMvaVar]=localSelJets[1].pt();
		    if(variable=="eta1")       tmvaVars[iMvaVar]=localSelJets[0].eta();
		    if(variable=="eta2")       tmvaVars[iMvaVar]=localSelJets[1].eta();
		    if(variable=="qg1")        tmvaVars[iMvaVar]=localSelJets[0].getVal("qgMVA");
		    if(variable=="qg2")        tmvaVars[iMvaVar]=localSelJets[1].getVal("qgMVA");
		  }
		if(tmvaReader) {
		  for(size_t im=0; im<tmvaMethods.size(); im++)
		    localTmvaDiscrVals[im]=tmvaReader->EvaluateMVA( tmvaMethods[im] );
		}
		
		//re-assign the event category;
		std::vector<TString> locTags(1,chTags[ich]);
		TString mjjCat("");
		std::vector<TString> localSelTags=getDijetCategories(mjj,hardpt,locTags,mjjCat);
		float finalWeight(iweight);
		if(gammaWgtHandler!=0) {
		  std::vector<Float_t> gammaWgtVars(1,selPhotons[0].pt());
		  finalWeight *= gammaWgtHandler->getWeightFor(gammaWgtVars,chTags[ich]+mjjCat);
		}

		//re-select the event and fill the shapes
		for(unsigned int index=0; index<optim_Cuts2_jet_pt1.size();index++)
		  {
		    float minJetPt1=optim_Cuts2_jet_pt1[index];
		    float minJetPt2=optim_Cuts2_jet_pt2[index];
		    bool passLocalJet1Pt(localSelJets[0].pt()>minJetPt1);
		    bool passLocalJet2Pt(localSelJets[1].pt()>minJetPt2);
		    if(!passLocalJet1Pt || !passLocalJet2Pt) continue; 
		    
		    mon.fillHisto(TString("dijet_deta_shapes")+varNames[ivar],localSelTags,index,detajj,finalWeight);
		    if(tmvaReader)
		      for(size_t im=0; im<tmvaMethods.size(); im++)
			mon.fillHisto(TString(tmvaMethods[im]+"_shapes")+varNames[ivar],localSelTags,index,localTmvaDiscrVals[im],finalWeight);
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
  file->Close();
  
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





