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
  bool filterOnlyEE(false), filterOnlyMUMU(false), filterOnlyMuEG(false);
  if(!isMC)
    {
      if(url.Contains("DoubleEle")) filterOnlyEE=true;
      if(url.Contains("DoubleMu"))  filterOnlyMUMU=true;
      if(url.Contains("MuEG"))      filterOnlyMuEG=true;
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

  //lepton efficiencies
  LeptonEfficiencySF lepEff;

  //systematics
  bool runSystematics                        = runProcess.getParameter<bool>("runSystematics");
  std::vector<TString> varNames(1,"");
  size_t nvarsToInclude(1);
  if(runSystematics && isMC)
    {
//      varNames.push_back("_jerup"); varNames.push_back("_jerdown");
//      varNames.push_back("_jesup"); varNames.push_back("_jesdown");
//      varNames.push_back("_puup"); varNames.push_back("_pudown");
      if(isSignal)
	{
//	  varNames.push_back("_q2up"); varNames.push_back("_q2down");
//	  varNames.push_back("_pdfup"); varNames.push_back("_pdfdown");
//	  varNames.push_back("_balanceup"); varNames.push_back("_balancedown");
	}
      nvarsToInclude=varNames.size();
      cout << nvarsToInclude << " systematics will be computed for this analysis" << endl;
    }

  //##############################################
  //########    INITIATING HISTOGRAMS     ########
  //##############################################
  SmartSelectionMonitor mon;

  TH1 *h=mon.addHistogram( new TH1F ("eventflow", ";;Events", 4,0,4) );
  h->GetXaxis()->SetBinLabel(1 ,"#eq 2 leptons");
  h->GetXaxis()->SetBinLabel(2 ,"MET>=0");
  h->GetXaxis()->SetBinLabel(3 ,"lepton Veto");


  h=mon.addHistogram( new TH1F ("eventflowBNV", ";;Events", 10,0,10) );
  h->GetXaxis()->SetBinLabel(1 ,"ee OS");
  h->GetXaxis()->SetBinLabel(2 ,"ee SS-");
  h->GetXaxis()->SetBinLabel(3 ,"ee SS+");   
  h->GetXaxis()->SetBinLabel(4 ,"mm OS");
  h->GetXaxis()->SetBinLabel(5 ,"mm SS-");
  h->GetXaxis()->SetBinLabel(6, "mm SS+");
  h->GetXaxis()->SetBinLabel(7, "em OS");
  h->GetXaxis()->SetBinLabel(8, "em SS-");
  h->GetXaxis()->SetBinLabel(9, "em SS+");

  h=mon.addHistogram( new TH1F ("eventflowBNV2", ";;Events", 7,0,7) );
  h->GetXaxis()->SetBinLabel(1 ,"ee SS-");
  h->GetXaxis()->SetBinLabel(2 ,"ee SS+");
  h->GetXaxis()->SetBinLabel(3 ,"mm SS-");
  h->GetXaxis()->SetBinLabel(4, "mm SS+");
  h->GetXaxis()->SetBinLabel(5, "em SS-");
  h->GetXaxis()->SetBinLabel(6, "em SS+");

  mon.addHistogram( new TH1F( "llmass",    ";M^{ll} [GeV];Events",100,0,500));
  mon.addHistogram( new TH1F( "llpt1",     ";l1^{pT} [GeV];Events",100,0,500));
  mon.addHistogram( new TH1F( "llpt2",     ";l2^{pT} [GeV];Events",100,0,500));
  mon.addHistogram( new TH1F( "lleta1",     ";l1^{#eta};Events",100,-2.6,2.6));
  mon.addHistogram( new TH1F( "lleta2",     ";l2^{#eta};Events",100,-2.6,2.6));
  mon.addHistogram( new TH1F( "lldR",     ";ll^{dR};Events",100,-6,6));
  mon.addHistogram( new TH2F( "lleta1vseta2",";l1^{#eta};l2^{#eta};Events",25,-2.6,2.6, 25,-2.6,2.6));
  mon.addHistogram( new TH2F( "lliso1vsiso2",";l1^{iso};l2^{iso};Events",25,0,0.15, 25,0,0.15));
  mon.addHistogram( new TH2F( "llpt1vspt2",";l1^{pT};l2^{pT};Events",25,0,200, 25,0,200));





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
  mon.addHistogram( new TH1F( "leppt", ";p_{T}^{l};Events", 50,0,500) );

  //tau control
  mon.addHistogram( new TH1F( "ntaus",      ";ntaus;Events", 6,0,6) );
  mon.addHistogram( new TH1F( "tauleadpt",  ";p_{T}^{#tau};Events", 50,0,500) );
  mon.addHistogram( new TH1F( "tauleadeta", ";#eta^{#tau};Events", 50,-2.6,2.6) );
  mon.addHistogram( new TH1F( "taupt",  ";p_{T}^{#tau};Events", 50,0,500) );
  mon.addHistogram( new TH1F( "taucharge",  ";p_{T}^{#tau};Events", 5,-2,2) );
  mon.addHistogram( new TH1F( "taudz",      ";dz^{#tau};Events", 50,0,10) );
  mon.addHistogram( new TH1F( "tauvz",      ";vz^{#tau};Events", 50,0,10) );


  //bjets control
  mon.addHistogram( new TH1F( "nbjets",      ";ntaus;Events", 6,0,6) );
  mon.addHistogram( new TH1F( "bjetpt",  ";p_{T}^{bjet};Events", 50,0,500) );
  mon.addHistogram( new TH1F( "bjetcsv", ";#eta^{#tau};Events", 50,0, 1) );

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
  mon.addHistogram( new TH1F( "higgsmass",    ";M^{#higgs} [GeV];Events",100,0,500));
  mon.addHistogram( new TH1F( "higgsmasssvfit",    ";M^{#higgs} [GeV];Events",100,0,500));
  mon.addHistogram( new TH1F( "higgsmet",    ";MET [GeV];Events",100,0,500));
  mon.addHistogram( new TH1F( "higgsnjets",   ";NJets;Events",10,0,10));

  //statistical analysis
  std::vector<double> optim_Cuts_sumpt; 
  std::vector<double> optim_Cuts_lepIso; 
  for(double sumpt=10;sumpt<=100;sumpt+=5)
    {
      for(double lepIso=0.4;lepIso>=0.05;lepIso-=0.05)
	{
	  optim_Cuts_sumpt.push_back(sumpt);
	  optim_Cuts_lepIso.push_back(lepIso);
	} 
    }
  TH2F* Hoptim_cuts  =(TH2F*)mon.addHistogram(new TProfile2D("optim_cut",      ";cut index;variable",       optim_Cuts_sumpt.size(),0,optim_Cuts_sumpt.size(), 2, 0, 2)) ;
  Hoptim_cuts->GetYaxis()->SetBinLabel(1, "jpt1>");
  Hoptim_cuts->GetYaxis()->SetBinLabel(2, "jpt2>");
  for(unsigned int index=0;index<optim_Cuts_sumpt.size();index++){
    Hoptim_cuts->Fill(index,0.0,optim_Cuts_sumpt[index]); 
    Hoptim_cuts->Fill(index,1.0,optim_Cuts_lepIso[index]); 
  }

  TH1F* Hoptim_systs     =  (TH1F*) mon.addHistogram( new TH1F ("optim_systs"    , ";syst;", nvarsToInclude,0,nvarsToInclude) ) ;
  for(size_t ivar=0; ivar<nvarsToInclude; ivar++)
  {
    Hoptim_systs->GetXaxis()->SetBinLabel(ivar+1, varNames[ivar]);
    mon.addHistogram( new TH2F (TString("svfit_shapes")+varNames[ivar],";cut index;|#Delta #eta|;Events",optim_Cuts_sumpt.size(),0,optim_Cuts_sumpt.size(),25,0,250) );
  }


  //tau fakeRate
  mon.addHistogram( new TH2F( "taufakerate",    ";CutIndex;p_{T}^{#tau} [GeV]",optim_Cuts_sumpt.size(),0,optim_Cuts_sumpt.size(),20,0,100));

  
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
                                //some work as to be done on the triggers
      bool eeTrigger          = true;//triggerBits[0];
      bool muTrigger          = true;//triggerBits[6];
      bool mumuTrigger        = true;//triggerBits[2] || triggerBits[3] || muTrigger;
      bool emTrigger          = true;//!(eeTrigger||mumuTrigger);
      if(filterOnlyEE)   { mumuTrigger=false; emTrigger=false; }
      if(filterOnlyMUMU) { eeTrigger=false;   emTrigger=false; }
      if(filterOnlyMuEG) { mumuTrigger=false; eeTrigger=false; }
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
	  if(leptons[ilep].pt()<25)                   passKin=false;
//          if(leptons[ilep].pt()<10)                   passKin=false;
	  if(leta> (lid==11 ? 2.5 : 2.4) )            passKin=false;
	  if(lid==11 && (leta>1.4442 && leta<1.5660)) passKin=false;

	  //id
	  Int_t idbits = leptons[ilep].idbits;
	  if(lid==11){
	    if(leptons[ilep].electronInfoRef->isConv)              passId=false;
	    bool isLoose = ((idbits >> 4) & 0x1);
            bool isTight = ((idbits >> 6) & 0x1);
	    if(!isTight)                                   passId=false;
 	  }
	  else{
	    bool isLoose    = ((idbits >> 8) & 0x1);
            bool isTight    = ((idbits >> 10) & 0x1);
	    if(!isTight)                                   passId=false;
	  }

	  //isolation
          float relIso = utils::cmssw::relIso(leptons[ilep], rho);
//	  if( (lid==11 && relIso>0.15) || (lid!=11 && relIso>0.20) ) passIso=false;
        if( (lid==11 && relIso>0.1) || (lid!=11 && relIso>0.12) ) passIso=false;
//          if( (lid==11 && relIso>0.40) || (lid!=11 && relIso>0.40) ) passIso=false;
	  
	  if(!passId || !passIso || !passKin) continue;
	  selLeptons.push_back(leptons[ilep]);
	}
      std::sort(selLeptons.begin(), selLeptons.end(), sort_llvvObjectByPt);

      //at this point check if it's worth continuig
      if(selLeptons.size()!=2) continue;
    
/*      //
      // DILEPTON ANALYSIS
      //
      LorentzVector leadingLep, trailerLep, zll;
      int dilLep1=-1, dilLep2=-1;
      double BestMass=0;
      //identify the best lepton pair
      for(unsigned int l1=0   ;l1<selLeptons.size();l1++){
         float relIso1 = utils::cmssw::relIso(selLeptons[l1], rho);
         if( relIso1>0.30 ) continue;
         for(unsigned int l2=l1+1;l2<selLeptons.size();l2++){
            if(fabs(selLeptons[l1].id)!=fabs(selLeptons[l2].id))continue; //only consider same flavor lepton pairs
            if(selLeptons[l1].id*selLeptons[l2].id>=0)continue; //only consider opposite charge lepton pairs
            if( !(selLeptons[l1].pt()>=20 && selLeptons[l2].pt()>=10) || !(selLeptons[l1].pt()>=10 && selLeptons[l2].pt()>=20))continue;
            float relIso2 = utils::cmssw::relIso(selLeptons[l2], rho);
            if( relIso2>0.30 ) continue;

            if(deltaR(selLeptons[l1], selLeptons[l2])<0.1)continue;
            if(fabs(BestMass-91.2)>fabs((selLeptons[l1]+selLeptons[l2]).mass() - 91.2)){
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
*/
      //apply data/mc correction factors
//      int dilId = selLeptons[dilLep1].id * selLeptons[dilLep2].id;
      weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[0].pt(), selLeptons[0].eta(), abs(selLeptons[0].id),  abs(selLeptons[0].id) ==11 ? "loose" : "loose" ).first : 1.0;
      weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[1].pt(), selLeptons[1].eta(), abs(selLeptons[1].id),  abs(selLeptons[1].id) ==11 ? "loose" : "loose" ).first : 1.0;

      //check the channel
      //prepare the tag's vectors for histo filling
      std::vector<TString> chTags;
      chTags.push_back("all");

/*
      //generator level
      LorentzVector genll;
      for(size_t ig=0; ig<gen.size(); ig++){
	  int pid=abs(gen[ig].id);
	  if(pid!=11 && pid!=13) continue;
	  genll += gen[ig];
      }
*/
      //
      //JET/MET ANALYSIS
      //
      llvvJetExtCollection selJets, selJetsNoId, selBJets;
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
	  if(jets[ijet].pt()<30 || fabs(jets[ijet].eta())>4.7 ) continue;


          //bjets
          mon.fillHisto("bjetpt"    ,  chTags, jets[ijet].pt(),  weight);
          mon.fillHisto("bjetcsv"   ,  chTags, jets[ijet].origcsv,  weight);
//          if(jets[ijet].pt()>20 && fabs(jets[ijet].eta())<2.4 && jets[ijet].origcsv>0.244){
          if(jets[ijet].pt()>20 && fabs(jets[ijet].eta())<2.4 && jets[ijet].origcsv>0.898){
             selBJets.push_back(jets[ijet]);  
             nbjets++;
          }


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

	  selJetsNoId.push_back(jets[ijet]);
	  if(passPFloose && passLooseSimplePuId){
	    selJets.push_back(jets[ijet]);
	    if(jets[ijet].pt()>minJetPtToApply) njets++;
	  }
	}
      std::sort(selJets.begin(), selJets.end(), sort_llvvObjectByPt);
      std::sort(selBJets.begin(), selBJets.end(), sort_llvvObjectByPt);
      mon.fillHisto("nbjets"    ,  chTags, nbjets,  weight);




/*
      //
      // TAU ANALYSIS
      //
      llvvTauCollection selTaus;
      for(size_t itau=0; itau<taus.size(); itau++){
         llvvTau& tau = taus[itau];
         if(tau.pt()<15.0 || fabs(tau.eta()) > 2.3)continue; 


         bool overalWithLepton=false;
         for(int l1=0   ;l1<(int)selLeptons.size();l1++){
            if(deltaR(tau, selLeptons[l1])<0.1){overalWithLepton=true; break;}
         }
         if(overalWithLepton)continue;

//         printf("TauId: "); for(unsigned int i=0;i<64;i++){printf("%i ", (int) ((tau.idbits>>i)&1));}printf("\n");

         if(!tau.passId(llvvTAUID::againstElectronLoose))continue;
         if(!tau.passId(llvvTAUID::againstMuonLoose2))continue; 
         if(!tau.passId(llvvTAUID::decayModeFinding))continue;
         if(!tau.passId(llvvTAUID::byLooseCombinedIsolationDeltaBetaCorr3Hits))continue;

         selTaus.push_back(tau);         
      }

      //
      // HIGGS ANALYSIS
      //

      LorentzVector higgsCand;
      int higgsCandId=0,  higgsCandMu=-1, higgsCandEl=-1, higgsCandT1=-1, higgsCandT2=-1;

      bool passNLep = (selLeptons.size()+selTaus.size())>=4;

      //Check if the event is compatible with a Mu-El candidate
      for(int l1=0   ;!higgsCandId && l1<(int)selLeptons.size();l1++){
      for(int l2=l1+1;!higgsCandId && l2<(int)selLeptons.size();l2++){
         if(l1==dilLep1 || l1==dilLep2 || l2==dilLep1 || l2==dilLep2)continue; //lepton already used in the dilepton pair
         if(selLeptons[l1].id*selLeptons[l2].id!=-143)continue;//Only consider opposite sign, opposite flavor pairs

         //check all deltaRs
         if(deltaR(selLeptons[l1], selLeptons[dilLep1])<0.1)continue;
         if(deltaR(selLeptons[l1], selLeptons[dilLep2])<0.1)continue;
         if(deltaR(selLeptons[l2], selLeptons[dilLep1])<0.1)continue;
         if(deltaR(selLeptons[l2], selLeptons[dilLep2])<0.1)continue;
         if(deltaR(selLeptons[l1], selLeptons[l2     ])<0.1)continue;

//         if(utils::cmssw::relIso(selLeptons[l1], rho)>0.25)continue;
//         if(utils::cmssw::relIso(selLeptons[l2], rho)>0.25)continue;

//         if((selLeptons[l1].pt() + selLeptons[l2].pt())<25)continue;
 
         int muId, elId;
         if(abs(selLeptons[l1].id)==13){muId=l1; elId=l2;}else{muId=l2; elId=l1;}

         higgsCand=selLeptons[muId]+selLeptons[elId]; higgsCandId=selLeptons[muId].id*selLeptons[elId].id;  higgsCandMu=muId; higgsCandEl=elId;
         break;//we found a candidate, stop the loop
      }}

      //Check if the event is compatible with a Lep-Tau candidate
      for(int l1=0;!higgsCandId && l1<(int)selLeptons.size();l1++){
      for(int t1=0;!higgsCandId && t1<(int)selTaus   .size();t1++){
         if(l1==dilLep1 || l1==dilLep2)continue; //lepton already used in the dilepton pair
         if(selLeptons[l1].id*selTaus[t1].id>=0)continue;//Only consider opposite sign pairs
         if(selTaus[t1].pt()<15 || fabs(selTaus[t1].eta())>2.3)continue;

         //check all deltaRs
         if(deltaR(selLeptons[l1], selLeptons[dilLep1])<0.1)continue;
         if(deltaR(selLeptons[l1], selLeptons[dilLep2])<0.1)continue;
         if(deltaR(selTaus[t1]   , selLeptons[dilLep1])<0.1)continue;
         if(deltaR(selTaus[t1]   , selLeptons[dilLep2])<0.1)continue;
         if(deltaR(selLeptons[l1], selTaus[t1        ])<0.1)continue;

//         if(utils::cmssw::relIso(selLeptons[l1], rho)>0.25)continue;

         if(abs(selLeptons[l1].id)==11 && (!selTaus[t1].passId(llvvTAUID::decayModeFinding) || !selTaus[t1].passId(llvvTAUID::againstElectronTightMVA3) || !selTaus[t1].passId(llvvTAUID::againstMuonLoose2) || !selTaus[t1].passId(llvvTAUID::byLooseCombinedIsolationDeltaBetaCorr3Hits) ))continue;
         if(abs(selLeptons[l1].id)!=11 && (!selTaus[t1].passId(llvvTAUID::decayModeFinding) || !selTaus[t1].passId(llvvTAUID::againstElectronLoose    ) || !selTaus[t1].passId(llvvTAUID::againstMuonTight2) || !selTaus[t1].passId(llvvTAUID::byLooseCombinedIsolationDeltaBetaCorr3Hits) ))continue;

//         if((selLeptons[l1].pt()+selTaus[t1].pt())<45 )continue;

         higgsCand=selLeptons[l1]+selTaus[t1]; higgsCandId=selLeptons[l1].id*selTaus[t1].id;  if(abs(selLeptons[l1].id)==11){higgsCandEl=l1;}else{higgsCandMu=l1;} higgsCandT1=t1;
         break;//we found a candidate, stop the loop         
      }}

      //Check if the event is compatible with a Tau-Tau candidate
      for(int t1=0   ;!higgsCandId && t1<(int)selTaus   .size();t1++){
      for(int t2=t1+1;!higgsCandId && t2<(int)selTaus   .size();t2++){

         //printf("Tau %6.2f %6.2f %+3i %i %i %i %i\n", selTaus[t1].pt(), fabs(selTaus[t1].eta()), selTaus[t1].id, (int)selTaus[t1].passId(llvvTAUID::decayModeFinding), (int)selTaus[t1].passId(llvvTAUID::againstElectronLoose), (int) selTaus[t1].passId(llvvTAUID::againstMuonLoose2), (int) selTaus[t1].passId(llvvTAUID::byLooseCombinedIsolationDeltaBetaCorr3Hits) );
         //printf("    %6.2f %6.2f %+3i %i %i %i %i\n", selTaus[t2].pt(), fabs(selTaus[t2].eta()), selTaus[t2].id, (int)selTaus[t2].passId(llvvTAUID::decayModeFinding), (int)selTaus[t2].passId(llvvTAUID::againstElectronLoose), (int) selTaus[t2].passId(llvvTAUID::againstMuonLoose2), (int) selTaus[t2].passId(llvvTAUID::byLooseCombinedIsolationDeltaBetaCorr3Hits) );

         if(selTaus[t1].pt()<15 || fabs(selTaus[t1].eta())>2.3)continue;
         if(selTaus[t2].pt()<15 || fabs(selTaus[t2].eta())>2.3)continue;

         //check all deltaRs
         if(deltaR(selTaus[t1]   , selLeptons[dilLep1])<0.1)continue;
         if(deltaR(selTaus[t1]   , selLeptons[dilLep2])<0.1)continue;
         if(deltaR(selTaus[t2]   , selLeptons[dilLep1])<0.1)continue;
         if(deltaR(selTaus[t2]   , selLeptons[dilLep2])<0.1)continue;
         if(deltaR(selTaus[t1]   , selTaus[t2        ])<0.1)continue;

         if(!selTaus[t1].passId(llvvTAUID::decayModeFinding) || !selTaus[t1].passId(llvvTAUID::againstElectronLoose) || !selTaus[t1].passId(llvvTAUID::againstMuonLoose2) || !selTaus[t1].passId(llvvTAUID::byLooseCombinedIsolationDeltaBetaCorr3Hits) )continue;
         if(!selTaus[t2].passId(llvvTAUID::decayModeFinding) || !selTaus[t2].passId(llvvTAUID::againstElectronLoose) || !selTaus[t2].passId(llvvTAUID::againstMuonLoose2) || !selTaus[t2].passId(llvvTAUID::byLooseCombinedIsolationDeltaBetaCorr3Hits) )continue;

//         if((selTaus[t1].pt()+selTaus[t2].pt())<75)continue;

         if(selTaus[t1].id*selTaus[t2].id<0){        
            higgsCand=selTaus[t1]+selTaus[t2]; higgsCandId=selTaus[t1].id*selTaus[t2].id;  higgsCandT1=t1; higgsCandT2=t2;
         }else{
            higgsCand=selTaus[t1]+selTaus[t2]; higgsCandId=selTaus[t1].id;  higgsCandT1=t1; higgsCandT2=t2;            
         }
         break;//we found a candidate, stop the loop
      }}


      //apply data/mc correction factors
      if(higgsCandMu!=-1)weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[higgsCandMu].pt(), selLeptons[higgsCandMu].eta(), abs(selLeptons[higgsCandMu].id),  abs(selLeptons[higgsCandMu].id) ==11 ? "loose" : "loose" ).first : 1.0;
      if(higgsCandEl!=-1)weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[higgsCandEl].pt(), selLeptons[higgsCandEl].eta(), abs(selLeptons[higgsCandEl].id),  abs(selLeptons[higgsCandEl].id) ==11 ? "loose" : "loose" ).first : 1.0;

      //check the channel
      //prepare the tag's vectors for histo filling
      bool passHiggs = abs(higgsCandId)>100;  int HiggsShortId = abs(selLeptons[dilLep1].id)==13?0:4;
           if( abs(higgsCandId)==143 ){ chTags.push_back(chTags[chTags.size()-1] + string("_elmu")); HiggsShortId+=0;}
      else if( abs(higgsCandId)==165 ){ chTags.push_back(chTags[chTags.size()-1] + string("_elha")); HiggsShortId+=1;}
      else if( abs(higgsCandId)==195 ){ chTags.push_back(chTags[chTags.size()-1] + string("_muha")); HiggsShortId+=2;}
      else if( abs(higgsCandId)==225 ){ chTags.push_back(chTags[chTags.size()-1] + string("_haha")); HiggsShortId+=3;}
      else if( abs(higgsCandId)== 15 ){ chTags.push_back(chTags[chTags.size()-1] + string("_haCtrl"));}
      else                              chTags.push_back(chTags[chTags.size()-1] + string("_none"));

      bool passLepVeto  = true;
      for(int l1=0;l1<(int)selLeptons.size();l1++){
         mon.fillHisto("leppt"    ,  chTags, selLeptons[l1].pt(),  weight);
         if(l1==dilLep1 || l1==dilLep2 || l1==higgsCandMu || l1==higgsCandEl)continue; //lepton already used in the dilepton pair or higgs candidate
         passLepVeto = false; break;
      }
      for(int t1=0;passLepVeto && t1<(int)selTaus   .size();t1++){         
         mon.fillHisto("taupt"    ,  chTags, selTaus[t1].pt(),  weight);
         mon.fillHisto("taucharge",  chTags, selTaus[t1].id/15.0, weight);
         mon.fillHisto("taudz"    ,  chTags, selTaus[t1].vz,  weight);
         mon.fillHisto("tauvz"    ,  chTags, selTaus[t1].z_expo,  weight);

         if(t1==higgsCandT1 || t1==higgsCandT2)continue; //lepton already used in the dilepton pair or higgs candidate
         if(selTaus[t1].pt()<20)continue;
         if(!selTaus[t1].passId(llvvTAUID::againstElectronLoose) || !selTaus[t1].passId(llvvTAUID::againstMuonLoose2) || !selTaus[t1].passId(llvvTAUID::byMediumCombinedIsolationDeltaBetaCorr3Hits)  )continue;
         passLepVeto = false; break;          
      } 

      bool passBJetVeto = true;;
      for(int j1=0;j1<(int)selBJets.size();j1++){
         if(                   deltaR(selBJets[j1]   , selLeptons[dilLep1    ])>0.4){passBJetVeto=false; break;}
         if(                   deltaR(selBJets[j1]   , selLeptons[dilLep2    ])>0.4){passBJetVeto=false; break;}
         if(higgsCandMu!=-1 && deltaR(selBJets[j1]   , selLeptons[higgsCandMu])>0.4){passBJetVeto=false; break;}
         if(higgsCandEl!=-1 && deltaR(selBJets[j1]   , selLeptons[higgsCandEl])>0.4){passBJetVeto=false; break;}
         if(higgsCandT1!=-1 && deltaR(selBJets[j1]   , selTaus   [higgsCandT1])>0.4){passBJetVeto=false; break;}
         if(higgsCandT2!=-1 && deltaR(selBJets[j1]   , selTaus   [higgsCandT2])>0.4){passBJetVeto=false; break;}
      }

      int NCleanedJet=0;
      for(int j1=0;j1<(int)selJets.size();j1++){
         if(                   deltaR(selJets[j1]   , selLeptons[dilLep1    ])<0.4)continue;
         if(                   deltaR(selJets[j1]   , selLeptons[dilLep2    ])<0.4)continue;
         if(higgsCandMu!=-1 && deltaR(selJets[j1]   , selLeptons[higgsCandMu])<0.4)continue;
         if(higgsCandEl!=-1 && deltaR(selJets[j1]   , selLeptons[higgsCandEl])<0.4)continue;
         if(higgsCandT1!=-1 && deltaR(selJets[j1]   , selTaus   [higgsCandT1])<0.4)continue;
         if(higgsCandT2!=-1 && deltaR(selJets[j1]   , selTaus   [higgsCandT2])<0.4)continue;
         NCleanedJet++;
      }



      //SVFIT MASS
      double diTauMass = -1;
      if(passHiggs && passLepVeto && passBJetVeto){
            //taken from https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorkingSummer2013
            TMatrixD covMET(2, 2); // PFMET significance matrix
            covMET[0][0] = met.sigx2;
            covMET[0][1] = met.sigxy;
            covMET[1][0] = met.sigxy;
            covMET[1][1] = met.sigy2;
            std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
            if(higgsCandMu!=-1)measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, NSVfitStandalone::LorentzVector(selLeptons[higgsCandMu].px(), selLeptons[higgsCandMu].py(), selLeptons[higgsCandMu].pz(), selLeptons[higgsCandMu].E()) ));
            if(higgsCandEl!=-1)measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, NSVfitStandalone::LorentzVector(selLeptons[higgsCandEl].px(), selLeptons[higgsCandEl].py(), selLeptons[higgsCandEl].pz(), selLeptons[higgsCandEl].E()) ));
            if(higgsCandT1!=-1)measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kHadDecay, NSVfitStandalone::LorentzVector(selTaus   [higgsCandT1].px(), selTaus   [higgsCandT1].py(), selTaus   [higgsCandT1].pz(), selTaus   [higgsCandT1].E()) ));
            if(higgsCandT2!=-1)measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kHadDecay, NSVfitStandalone::LorentzVector(selTaus   [higgsCandT2].px(), selTaus   [higgsCandT2].py(), selTaus   [higgsCandT2].pz(), selTaus   [higgsCandT2].E()) ));
            NSVfitStandaloneAlgorithm algo(measuredTauLeptons, NSVfitStandalone::Vector(met.px(), met.py(), 0) , covMET, 0);
            algo.addLogM(false);
            algo.integrateMarkovChain();
            //algo.integrateVEGAS(); ////Use this instead for VEGAS integration
            if(algo.isValidSolution()){
               diTauMass     = algo.getMass(); 
               double diTauMassErr = algo.massUncert();
            }
      }
*/




      //check the channel
      //prepare the tag's vectors for histo filling
      int ChannelId = abs(selLeptons[0].id * selLeptons[1].id) + int(selLeptons[0].id>0?1:0) + int(selLeptons[1].id>0?1:0);
      int ChannelShortId = -1, ChannelShortId2 = -1;
           if( (ChannelId == 11*11 + 1) && eeTrigger){ chTags.push_back(chTags[chTags.size()-1] + string("_eeOS")) ; ChannelShortId=0;}
      else if( (ChannelId == 11*11 + 0) && eeTrigger){ chTags.push_back(chTags[chTags.size()-1] + string("_eeSSp")); ChannelShortId=2; ChannelShortId2=1;}
      else if( (ChannelId == 11*11 + 2) && eeTrigger){ chTags.push_back(chTags[chTags.size()-1] + string("_eeSSm")); ChannelShortId=1; ChannelShortId2=0;}
      else if( (ChannelId == 13*13 + 1) && mumuTrigger){ chTags.push_back(chTags[chTags.size()-1] + string("_mmOS")) ; ChannelShortId=3;}
      else if( (ChannelId == 13*13 + 0) && mumuTrigger){ chTags.push_back(chTags[chTags.size()-1] + string("_mmSSp")); ChannelShortId=5; ChannelShortId2=3;}
      else if( (ChannelId == 13*13 + 2) && mumuTrigger){ chTags.push_back(chTags[chTags.size()-1] + string("_mmSSm")); ChannelShortId=4; ChannelShortId2=2;}
      else if( (ChannelId == 11*13 + 1 && emTrigger) ){ chTags.push_back(chTags[chTags.size()-1] + string("_emOS")) ; ChannelShortId=6;} 
      else if( (ChannelId == 11*13 + 0 && emTrigger) ){ chTags.push_back(chTags[chTags.size()-1] + string("_emSSp")); ChannelShortId=8; ChannelShortId2=5;}
      else if( (ChannelId == 11*13 + 2 && emTrigger) ){ chTags.push_back(chTags[chTags.size()-1] + string("_emSSm")); ChannelShortId=7; ChannelShortId2=4;}
      else                                 chTags.push_back(chTags[chTags.size()-1] + string("_none"));

//      printf("%+2i - %+2i --> %3i %3i --> %3i  --> %s\n", selLeptons[0].id,  selLeptons[1].id, abs(selLeptons[0].id * selLeptons[1].id), int(selLeptons[0].id>0?1:0) + int(selLeptons[1].id>0?1:0), ChannelId, chTags[chTags.size()-1].Data());

      if(selJets.size()==0)chTags.push_back(chTags[chTags.size()-1] + "eq0jet");
      if(selJets.size()==1)chTags.push_back(chTags[chTags.size()-1] + "eq1jet");
      if(selJets.size()==2)chTags.push_back(chTags[chTags.size()-1] + "eq2jet");
      if(selJets.size()>=3)chTags.push_back(chTags[chTags.size()-1] + "geq3jet");
     
      bool passLepVeto  = deltaR(selLeptons[0]   , selLeptons[1])>0.3;
      bool passMET = (met.pt()>30);
             
      //
      // NOW FOR THE CONTROL PLOTS
      //
      //start analysis
//      for(size_t ich=0; ich<chTags.size(); ich++)
//	{
//	  std::vector<TString> tags(1,chTags[ich]);

            	                                                   mon.fillHisto("eventflow",chTags,0,weight);
	  if(passMET)                                              mon.fillHisto("eventflow",chTags,1,weight);
	  if(passMET && passLepVeto)                               mon.fillHisto("eventflow",chTags,2,weight);
//          if(false)                                              mon.fillHisto("eventflow",chTags,3,weight);
          if(passMET && passLepVeto && ChannelShortId>=0)          mon.fillHisto("eventflowBNV",chTags,ChannelShortId,weight);
          if(passMET && passLepVeto && ChannelShortId2>=0)          mon.fillHisto("eventflowBNV2",chTags,ChannelShortId2,weight);

          if(passMET && passLepVeto && ChannelShortId>=0){
             mon.fillHisto("llmass",chTags,(selLeptons[0]+selLeptons[1]).pt(),weight);
             mon.fillHisto("llpt1",chTags,selLeptons[0].pt(),weight);
             mon.fillHisto("llpt2",chTags,selLeptons[1].pt(),weight);
             mon.fillHisto("lleta1",chTags,selLeptons[0].eta(),weight);
             mon.fillHisto("lleta2",chTags,selLeptons[1].eta(),weight);
             mon.fillHisto("lldR",chTags,deltaR(selLeptons[0],selLeptons[1]),weight);
             mon.fillHisto("lleta1vseta2",chTags,selLeptons[0].eta(), selLeptons[1].eta(),weight);
             mon.fillHisto("llpt1vspt2",chTags,selLeptons[0].pt(), selLeptons[1].pt(),weight);
             mon.fillHisto("lliso1vsiso2",chTags,utils::cmssw::relIso(selLeptons[0], rho), utils::cmssw::relIso(selLeptons[1], rho),weight);




          }

/*
	  mon.fillHisto("zmass",    chTags, zll.mass(), weight);  

        if(passMET && passBJetVeto && passLepVeto && !passHiggs && abs(higgsCandId)== 15){
            mon.fillHisto("higgsmass"    , chTags, higgsCand.mass(),  weight);
           
            mon.fillHisto("taufakerate"     ,   chTags, 1,  selTaus   [higgsCandT1].jet.pt(),    weight);
            mon.fillHisto("taufakerate"     ,   chTags, 1,  selTaus   [higgsCandT2].jet.pt(),    weight);

            for(size_t itau=0; itau<taus.size(); itau++){
               llvvTau& tau = taus[itau];
               if(tau.pt()<15.0 || fabs(tau.eta()) > 2.3)continue;
               mon.fillHisto("taufakerate"     ,   chTags, 0,  tau.jet.pt(),    weight);
            }

        }


//	  if(passMET){
        if(passMET && passBJetVeto && passLepVeto && passHiggs){
	
	    //pu control
            mon.fillHisto("nvtx"     ,   chTags, nvtx,      weight);
	    mon.fillHisto("nvtxraw"  ,   chTags, nvtx,      weight/puWeight);
	    mon.fillHisto("rho"      ,   chTags, rho,       weight);
            mon.fillHisto("rho25"    ,   chTags, rho25,     weight);
	
	    //Z kinematics control
	    mon.fillHisto("leadpt"      ,   chTags, leadingLep.pt(), weight);      
	    mon.fillHisto("trailerpt"   ,   chTags, trailerLep.pt(), weight);      
	    mon.fillHisto("leadeta"     ,   chTags, TMath::Max(fabs(leadingLep.eta()),fabs(trailerLep.eta())), weight);      
	    mon.fillHisto("trailereta"  ,   chTags, TMath::Min(fabs(leadingLep.eta()),fabs(trailerLep.eta())), weight);      

	    mon.fillHisto("zpt"      , chTags, zll.pt(),  weight);      
	    mon.fillHisto("zeta"     , chTags, zll.eta(), weight);
	    mon.fillHisto("zy"       , chTags, zy,        weight);
	    
//	    if(passZpt && passZeta){
          if(passBJetVeto && passLepVeto){
	  
	      //analyze dilepton kinematics
	      mon.fillHisto("leadeta"   ,  chTags, leadingLep.eta(), weight);
	      mon.fillHisto("leadpt"    ,  chTags, leadingLep.pt(),  weight);
	      mon.fillHisto("trailereta",  chTags, trailerLep.eta(), weight);
	      mon.fillHisto("trailerpt",   chTags, trailerLep.pt(),  weight);
 

              mon.fillHisto("ntaus"        ,  chTags, selTaus.size(), weight);
              mon.fillHisto("tauleadpt"    ,  chTags, selTaus.size()>0?selTaus[0].pt():-1,  weight);
              mon.fillHisto("tauleadeta"   ,  chTags, selTaus.size()>0?selTaus[0].eta():-10, weight);

              if(passHiggs){
                 mon.fillHisto("higgspt"      , chTags, higgsCand.pt(),    weight);
                 mon.fillHisto("higgsmass"    , chTags, higgsCand.mass(),  weight);
                 mon.fillHisto("higgsmasssvfit", chTags, diTauMass,  weight);
                 mon.fillHisto("higgsnjets"   , chTags, NCleanedJet      , weight); 
                 mon.fillHisto("higgsmet"     , chTags, met.pt()         , weight);
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

		//re-assign the event category;
		std::vector<TString> locTags = chTags;

                double sumPt = 0;                  
                if(higgsCandMu!=-1) sumPt+= selLeptons[higgsCandMu].pt();
                if(higgsCandEl!=-1) sumPt+= selLeptons[higgsCandEl].pt();
                if(higgsCandT1!=-1) sumPt+= selTaus   [higgsCandT1].pt();
                if(higgsCandT2!=-1) sumPt+= selTaus   [higgsCandT2].pt();

                double Lep1Iso=-1, Lep2Iso=-1;
                if(higgsCandMu!=-1)Lep1Iso = utils::cmssw::relIso(selLeptons[higgsCandMu], rho);
                if(higgsCandEl!=-1)Lep2Iso = utils::cmssw::relIso(selLeptons[higgsCandEl], rho);

		for(unsigned int index=0; index<optim_Cuts_sumpt.size();index++)
		  {
                    bool passSumPt = sumPt>=optim_Cuts_sumpt[index];
                    bool passIso   = Lep1Iso<=optim_Cuts_lepIso[index] && Lep2Iso<=optim_Cuts_lepIso[index];

//		    TString mjjCat("");
//		    std::vector<TString> localSelTags=getDijetCategories(mjj,detajj,locTags,mjjCat);
                    std::vector<TString> localSelTags=locTags;
                    if(passSumPt && passIso)mon.fillHisto(TString("svfit_shapes")+varNames[ivar],localSelTags,index,diTauMass,iweight);		    
		  }

	      }
	    }//end passZpt && passZeta

	  }//end passMET
*/
	}

//  }
  
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





