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
  bool runLoosePhotonSelection(false);

  bool Cut_tautau_MVA_iso = true;//currently hardcoded

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
  if(runSystematics && isMC)
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

  //##############################################
  //########    INITIATING HISTOGRAMS     ########
  //##############################################
  SmartSelectionMonitor mon;

  TH1 *h=mon.addHistogram( new TH1F ("eventflow", ";;Events", 13,0,13) );
  h->GetXaxis()->SetBinLabel(1,"#geq 2 leptons");
  h->GetXaxis()->SetBinLabel(2,"|M-M_{Z}|<15");
  h->GetXaxis()->SetBinLabel(3,"Btag Veto");
  h->GetXaxis()->SetBinLabel(4, "lepton Veto");
  h->GetXaxis()->SetBinLabel(5,"HiggsCandidate"); 
  h->GetXaxis()->SetBinLabel(6,"ee_em");
  h->GetXaxis()->SetBinLabel(7,"ee_et");
  h->GetXaxis()->SetBinLabel(8,"ee_mt");
  h->GetXaxis()->SetBinLabel(9,"ee_tt");
  h->GetXaxis()->SetBinLabel(10,"mm_em");
  h->GetXaxis()->SetBinLabel(11,"mm_et");
  h->GetXaxis()->SetBinLabel(12,"mm_mt");
  h->GetXaxis()->SetBinLabel(13,"mm_tt");

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

  //tau control
  mon.addHistogram( new TH1F( "ntaus",      ";ntaus;Events", 6,0,6) );
  mon.addHistogram( new TH1F( "tauleadpt",  ";p_{T}^{#tau};Events", 50,0,500) );
  mon.addHistogram( new TH1F( "tauleadeta", ";#eta^{#tau};Events", 50,-2.6,2.6) );

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

  //higgs control
  mon.addHistogram( new TH1F( "higgspt",      ";p_{T}^{#higgs} [GeV];Events",500,0,1500));
  mon.addHistogram( new TH1F( "higgsmass",    ";M^{#higgs} [GeV];Events",500,0,1500));


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
//  PDFInfo *mPDFInfo=0;
//  if(isMC)
//    {
//      TString pdfUrl(url);
//      pdfUrl.ReplaceAll(".root","_pdf.root");
//      pdfUrl.ReplaceAll("/MC","/pdf/MC");
//      mPDFInfo=new PDFInfo(pdfUrl,"cteq66.LHgrid");
//      for(int i=0; i<mPDFInfo->numberPDFs(); i++)
//	{
//	  TString var("_"); var+=i;
//	  mon.addHistogram( new TH1F("vbfcandjetdeta"+var    , ";|#Delta #eta|;Jets",                        50,0,10) );
//	  mon.addHistogram( new TH1F("vbfcandjet1eta"+var    , ";#eta;Jets",                                 50,0,5) );
//	  mon.addHistogram( new TH1F("vbfcandjet1pt"+var     , ";p_{T} [GeV];Jets",                        50,0,500) );
//	  mon.addHistogram( new TH1F("vbfcandjet2eta"+var    , ";#eta;Jets",                                 50,0,5) );
//	  mon.addHistogram( new TH1F("vbfcandjet2pt"+var     , ";p_{T} [GeV];Jets",                        50,0,500) );
//	}
//    }

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

      fwlite::Handle< llvvTauCollection > tauCollHandle;
      tauCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
      if(!tauCollHandle.isValid()){printf("llvvLeptonCollection Object NotFound\n");continue;}
      llvvTauCollection taus = *tauCollHandle;

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

	  
	  //kinematics
//	  float leta = lid==11 ? leptons[ilep].getVal("sceta") : leptons[ilep].eta();
          float leta = lid==11 ? leptons[ilep].electronInfoRef->sceta : leptons[ilep].eta();
//	  if(leptons[ilep].pt()<20)                   passKin=false;
          if(leptons[ilep].pt()<10)                   passKin=false;
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
          float relIso = utils::cmssw::relIso(leptons[ilep], rho);
//	  if( (lid==11 && relIso>0.15) || (lid!=11 && relIso>0.20) ) passIso=false;
          if( (lid==11 && relIso>0.40) || (lid!=11 && relIso>0.40) ) passIso=false;
	  
	  if(!passId || !passIso || !passKin) continue;
	  selLeptons.push_back(leptons[ilep]);
	}
      std::sort(selLeptons.begin(), selLeptons.end(), sort_llvvObjectByPt);

      //at this point check if it's worth continuig
      if(selLeptons.size()<2) continue;
    
      //
      // DILEPTON ANALYSIS
      //
      LorentzVector leadingLep, trailerLep, zll;
      int dilLep1=-1, dilLep2=-1;
      double BestMass=0;
      //identify the best lepton pair
      for(unsigned int l1=0   ;l1<selLeptons.size();l1++){
         float relIso1 = utils::cmssw::relIso(selLeptons[l1], rho);
         if( (abs(selLeptons[l1].id)==11 && relIso1>0.15) || (abs(selLeptons[l1].id)!=11 && relIso1>0.20) ) continue;
         for(unsigned int l2=l1+1;l2<selLeptons.size();l2++){
            if(selLeptons[l1].pt()<20 || selLeptons[l2].pt()<20)continue;
            float relIso2 = utils::cmssw::relIso(selLeptons[l2], rho);
            if( (abs(selLeptons[l2].id)==11 && relIso2>0.15) || (abs(selLeptons[l2].id)!=11 && relIso2>0.20) ) continue;
            if(fabs(selLeptons[l1].id)!=fabs(selLeptons[l2].id))continue; //only consider same flavor lepton pairs
            if(fabs(BestMass-91.2)>((selLeptons[l1]+selLeptons[l2]).mass() - 91.2)){
               dilLep1 = l1; 
               dilLep2 = l2;
               zll=selLeptons[l1]+selLeptons[l2];
               BestMass=zll.mass();
            }
         }
      }
      if(selLeptons[dilLep1].pt()>selLeptons[dilLep2].pt()){ leadingLep=selLeptons[dilLep1]; trailerLep=selLeptons[dilLep2]; 
      }else{                                                 leadingLep=selLeptons[dilLep2]; trailerLep=selLeptons[dilLep1]; 
      } 
      float zy(zll.Rapidity());
      bool passZmass(fabs(zll.mass()-91.2)<30);//15);
      bool passZpt(true);//zll.pt()>50);
      bool passZeta(true);//fabs(zll.eta())<1.4442);

      //apply data/mc correction factors
      int dilId = selLeptons[dilLep1].id * selLeptons[dilLep2].id;
      weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[dilLep1].pt(), selLeptons[dilLep1].eta(), abs(selLeptons[dilLep1].id),  abs(selLeptons[dilLep1].id) ==11 ? "loose" : "loose" ).first : 1.0;
      weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[dilLep2].pt(), selLeptons[dilLep2].eta(), abs(selLeptons[dilLep2].id),  abs(selLeptons[dilLep2].id) ==11 ? "loose" : "loose" ).first : 1.0;

      //check the channel
      //prepare the tag's vectors for histo filling
      std::vector<TString> chTags;
      chTags.push_back("all");
      if( abs(dilId)==121 && eeTrigger  ) chTags.push_back("ee");
      if( abs(dilId)==169 && mumuTrigger) chTags.push_back("mumu"); 
      if(chTags.size()==0) continue;


      //generator level
      LorentzVector genll;
      for(size_t ig=0; ig<gen.size(); ig++){
	  int pid=abs(gen[ig].id);
	  if(pid!=11 && pid!=13) continue;
	  genll += gen[ig];
      }

      //
      //JET/MET ANALYSIS
      //
      llvvJetExtCollection selJets, selJetsNoId;
      int njets(0), nbjets(0);
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

          //bjets
          if(jets[ijet].pt()>20 && fabs(jets[ijet].eta())<2.4 && jets[ijet].origcsv>0.244)nbjets++;

	  selJetsNoId.push_back(jets[ijet]);
	  if(passPFloose && passLooseSimplePuId){
	    selJets.push_back(jets[ijet]);
	    if(jets[ijet].pt()>minJetPtToApply) njets++;
	  }
	}
      std::sort(selJets.begin(), selJets.end(), sort_llvvObjectByPt);




      //
      // TAU ANALYSIS
      //
      llvvTauCollection selTaus;
      for(size_t itau=0; itau<taus.size(); itau++){
         llvvTau& tau = taus[itau];
         if(tau.pt()>15.0 && fabs(tau.eta()) < 2.3 && 
           (tau.idbits&(1<<llvvTAUID::againstElectronLoose))!=0 && 
           (tau.idbits&(1<<llvvTAUID::againstMuonLoose2))!=0 && 
           (tau.idbits&(1<<llvvTAUID::decayModeFinding))!=0){

//            if(tau.id==0)tau.id = 15;;

            bool overlap=false;
//            for(size_t ilep=0; ilep<selLeptons.size(); ilep++){
//               if(deltaR(tau,selLeptons[ilep])<0.1){overlap=true; break;}
//            }
            if(!overlap)selTaus.push_back(tau);
         }
      }

      //
      // HIGGS ANALYSIS
      //

      LorentzVector higgsCand;
      int higgsCandId=0,  higgsCandMu=-1, higgsCandEl=-1, higgsCandT1=-1, higgsCandT2=-1;

      printf("NLeptons=%i  NTaus = %i CheckEvent=%i\n", (int)selLeptons.size()-2, (int)selTaus   .size(),  ((int)selLeptons.size()-2)+(int)selTaus   .size());

      //Check if the event is compatible with a Mu-El candidate
      for(int l1=0   ;!higgsCandId && l1<(int)selLeptons.size();l1++){
      for(int l2=l1+1;!higgsCandId && l2<(int)selLeptons.size();l2++){
         if(l1==dilLep1 || l1==dilLep2 || l2==dilLep1 || l2==dilLep2)continue; //lepton already used in the dilepton pair
         if(selLeptons[l1].id*selLeptons[l2].id!=-143)continue;//Only consider opposite sign, opposite flavor pairs
 
         int muId, elId;
         if(abs(selLeptons[l1].id)==13){muId=l1; elId=l2;}else{muId=l2; elId=l1;}

         printf("Check if it's a EMU event: Id1 %i Id2 %i Iso1 %f Iso2 %f DR %f SumPt %f\n", selLeptons[l1].id, selLeptons[l2].id, utils::cmssw::relIso(selLeptons[muId], rho), utils::cmssw::relIso(selLeptons[elId], rho), deltaR(selLeptons[muId], selLeptons[elId]), (selLeptons[muId].pt()+selLeptons[elId].pt()) );

         if(utils::cmssw::relIso(selLeptons[muId], rho)<0.25 && utils::cmssw::relIso(selLeptons[elId], rho)<0.25 && deltaR(selLeptons[muId], selLeptons[elId])>0.1 && (selLeptons[muId].pt()+selLeptons[elId].pt())>25 ){        
           printf("Is a candidate\n");
           higgsCand=selLeptons[muId]+selLeptons[elId]; higgsCandId=selLeptons[muId].id*selLeptons[elId].id;  higgsCandMu=muId; higgsCandEl=elId;
           break;//we found a candidate, stop the loop
         }
      }}

      //Check if the event is compatible with a Lep-Tau candidate
      for(int l1=0;!higgsCandId && l1<(int)selLeptons.size();l1++){
      for(int t1=0;!higgsCandId && t1<(int)selTaus   .size();t1++){
         if(l1==dilLep1 || l1==dilLep2)continue; //lepton already used in the dilepton pair
         if(selLeptons[l1].id*selTaus[t1].id>0)continue;//Only consider opposite sign pairs
         if(selTaus[t1].pt()<15 || fabs(selTaus[t1].eta())>2.3)continue;

         float lepIso = utils::cmssw::relIso(selLeptons[l1], rho);
         bool pasLepIso = (abs(selLeptons[l1].id)==13 && lepIso<0.25) || (abs(selLeptons[l1].id)==11 && lepIso<0.15);

         bool passTauIso = true;
         if(Cut_tautau_MVA_iso){passTauIso&=(selTaus[t1].idbits&(1<<llvvTAUID::byLooseIsolationMVA))!=0;}else{passTauIso&=(selTaus[t1].idbits&(1<<llvvTAUID::byLooseCombinedIsolationDeltaBetaCorr3Hits))!=0;}
         if(abs(selLeptons[l1].id)==11){passTauIso&=(selTaus[t1].idbits&(1<<againstElectronTightMVA3))!=0;}else{passTauIso&=(selTaus[t1].idbits&(1<<llvvTAUID::againstMuonTight2))!=0;}

         printf("Check if it's a LTau event: Id1 %i Id2 %i Iso1 %i Iso2 %i DR %f SumPt %f\n", selLeptons[l1].id, selTaus[t1].id, (int)pasLepIso, (int)passTauIso, deltaR(selLeptons[l1], selTaus[t1]), (selLeptons[l1].pt()+selTaus[t1].pt()));

         if(pasLepIso && passTauIso && deltaR(selLeptons[l1], selTaus[t1])>0.3 && (selLeptons[l1].pt()+selTaus[t1].pt())>45 ){
           higgsCand=selLeptons[l1]+selTaus[t1]; higgsCandId=selLeptons[l1].id*selTaus[t1].id;  if(abs(selLeptons[l1].id)==11){higgsCandEl=l1;}else{higgsCandEl=l1;} higgsCandT1=t1;
           printf("Is a candidate\n");
           break;//we found a candidate, stop the loop
         }
      }}

      //Check if the event is compatible with a Tau-Tau candidate
      for(int t1=0   ;!higgsCandId && t1<(int)selTaus   .size();t1++){
      for(int t2=t1+1;!higgsCandId && t2<(int)selTaus   .size();t2++){
         if(selTaus[t1].id*selTaus[t2].id>0)continue;//Only consider opposite sign pairs
         if(selTaus[t1].pt()<15 || fabs(selTaus[t1].eta())>2.3)continue;
         if(selTaus[t2].pt()<15 || fabs(selTaus[t2].eta())>2.3)continue;

         bool passTauIso = true;
         if(Cut_tautau_MVA_iso){passTauIso&=(selTaus[t1].idbits&(1<<llvvTAUID::byLooseIsolationMVA))!=0;}else{passTauIso&=(selTaus[t1].idbits&(1<<llvvTAUID::byLooseCombinedIsolationDeltaBetaCorr3Hits))!=0;}
         if(Cut_tautau_MVA_iso){passTauIso&=(selTaus[t2].idbits&(1<<llvvTAUID::byLooseIsolationMVA))!=0;}else{passTauIso&=(selTaus[t2].idbits&(1<<llvvTAUID::byLooseCombinedIsolationDeltaBetaCorr3Hits))!=0;}

         printf("Check if it's a TauTau event: Id1 %i Id2 %i Iso %i DR %f SumPt %f\n", selTaus[t1].id, selTaus[t2].id, (int)passTauIso, deltaR(selTaus[t1], selTaus[t2]), (selTaus[t1].pt()+selTaus[t2].pt()));

         if(passTauIso && deltaR(selTaus[t1], selTaus[t2])>0.3 && (selTaus[t1].pt()+selTaus[t2].pt())>75 ){
           higgsCand=selTaus[t1]+selTaus[t2]; higgsCandId=selTaus[t1].id*selTaus[t2].id;  higgsCandT1=t1; higgsCandT2=t2;
           printf("Is a candidate\n");
           break;//we found a candidate, stop the loop
         }
      }}

      //apply data/mc correction factors
      if(higgsCandMu!=-1)weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[higgsCandMu].pt(), selLeptons[higgsCandMu].eta(), abs(selLeptons[higgsCandMu].id),  abs(selLeptons[higgsCandMu].id) ==11 ? "loose" : "loose" ).first : 1.0;
      if(higgsCandEl!=-1)weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[higgsCandEl].pt(), selLeptons[higgsCandEl].eta(), abs(selLeptons[higgsCandEl].id),  abs(selLeptons[higgsCandEl].id) ==11 ? "loose" : "loose" ).first : 1.0;

      //check the channel
      //prepare the tag's vectors for histo filling
      bool passHiggs = abs(higgsCandId)>0;  int HiggsShortId = abs(selLeptons[dilLep1].id)==11?0:4;
           if( abs(higgsCandId)==143 ){ chTags.push_back(chTags[chTags.size()-1] + string("_elmu")); HiggsShortId+=0;}
      else if( abs(higgsCandId)==165 ){ chTags.push_back(chTags[chTags.size()-1] + string("_elha")); HiggsShortId+=1;}
      else if( abs(higgsCandId)==195 ){ chTags.push_back(chTags[chTags.size()-1] + string("_muha")); HiggsShortId+=2;}
      else if( abs(higgsCandId)==225 ){ chTags.push_back(chTags[chTags.size()-1] + string("_haha")); HiggsShortId+=3;}
      else                             chTags.push_back(chTags[chTags.size()-1] + string("_none"));
      printf("event is %+4i %s\n", higgsCandId, chTags[chTags.size()-1].Data());

      bool passBJetVeto = (nbjets==0);
      bool passLepVeto  = true;
      for(int l1=0;l1<(int)selLeptons.size();l1++){
         if(l1==dilLep1 || l1==dilLep2 || l1==higgsCandMu || l1==higgsCandEl)continue; //lepton already used in the dilepton pair or higgs candidate
         passLepVeto = false; break;
      }
      for(int t1=0;passLepVeto && t1<(int)selTaus   .size();t1++){
         if(t1==higgsCandT1 || t1==higgsCandT2)continue; //lepton already used in the dilepton pair or higgs candidate
         if(selTaus[t1].pt()>20 && ((selTaus[t1].idbits&(1<<llvvTAUID::byMediumIsolationMVA))!=0) && ((selTaus[t1].idbits&(1<<againstElectronLooseMVA3))!=0) && ((selTaus[t1].idbits&(1<<llvvTAUID::againstMuonLoose2))!=0) ){
            passLepVeto = false; break; 
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
	  if(passZmass && passBJetVeto)                           mon.fillHisto("eventflow",tags,2,weight);
	  if(passZmass && passBJetVeto && passLepVeto)               mon.fillHisto("eventflow",tags,3,weight);
	  if(passZmass && passBJetVeto && passLepVeto && passHiggs)  mon.fillHisto("eventflow",tags,4,weight);
          if(passZmass && passBJetVeto && passLepVeto && passHiggs)  mon.fillHisto("eventflow",tags,5+HiggsShortId,weight);

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
	    
//	    if(passZpt && passZeta){
          if(passBJetVeto && passLepVeto){
	  
	      //analyze dilepton kinematics
	      mon.fillHisto("leadeta"   ,  tags, leadingLep.eta(), weight);
	      mon.fillHisto("leadpt"    ,  tags, leadingLep.pt(),  weight);
	      mon.fillHisto("trailereta",  tags, trailerLep.eta(), weight);
	      mon.fillHisto("trailerpt",   tags, trailerLep.pt(),  weight);
 

              mon.fillHisto("ntaus"        ,  tags, selTaus.size(), weight);
              mon.fillHisto("tauleadpt"    ,  tags, selTaus.size()>0?selTaus[0].pt():-1,  weight);
              mon.fillHisto("tauleadeta"   ,  tags, selTaus.size()>0?selTaus[0].eta():-10, weight);

              if(passHiggs){
                 mon.fillHisto("higgspt"      , tags, higgsCand.pt(),    weight);
                 mon.fillHisto("higgsmass"    , tags, higgsCand.mass(),  weight);
              }




 
	      //STATISTICAL ANALYSIS
	      float Q2Weight_plus(1.0), Q2Weight_down(1.0);
	      float PDFWeight_plus(1.0), PDFWeight_down(1.0);
	      if(isSignal){
//		  if(mPDFInfo){
//		      std::vector<float> wgts=mPDFInfo->getWeights(iev);
//		      for(size_t ipw=0; ipw<wgts.size(); ipw++){
//			  PDFWeight_plus = TMath::Max(PDFWeight_plus,wgts[ipw]);
//			  PDFWeight_down = TMath::Min(PDFWeight_down,wgts[ipw]);
//			}
//		  }
	      }

	      for(size_t ivar=0; ivar<nvarsToInclude; ivar++){
		float iweight = weight;                                                //nominal
		if(ivar==5)                        iweight *= TotalWeight_plus;        //pu up
		if(ivar==6)                        iweight *= TotalWeight_minus;       //pu down
		if(ivar==7)                        iweight *= Q2Weight_plus;
		if(ivar==8)                        iweight *= Q2Weight_down;
		if(ivar==9)                        iweight *= PDFWeight_plus;
		if(ivar==10)                       iweight *= PDFWeight_down;

/*	    
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
		    mon.fillHisto(TString("dijet_deta_shapes")+varNames[ivar],localSelTags,index,detajj,iweight);		    
		  }
*/
	      }
	    }//end passZpt && passZeta

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





