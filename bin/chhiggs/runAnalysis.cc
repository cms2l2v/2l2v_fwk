#include <iostream>
#include <boost/shared_ptr.hpp>

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

//Load here all the dataformat that we will need
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
//#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h" //for svfit

#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/HiggsUtils.h"
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/llvv_fwk/interface/TMVAUtils.h"
#include "UserCode/llvv_fwk/interface/LeptonEfficiencySF.h"
#include "UserCode/llvv_fwk/interface/PDFInfo.h"
#include "UserCode/llvv_fwk/interface/MuScleFitCorrector.h"
#include "UserCode/llvv_fwk/interface/BtagUncertaintyComputer.h"
#include "UserCode/llvv_fwk/interface/GammaWeightsHandler.h"

#include "UserCode/llvv_fwk/interface/PatUtils.h"


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

  TString suffix=runProcess.getParameter<std::string>("suffix");
  std::vector<std::string> urls=runProcess.getParameter<std::vector<std::string> >("input");
  TString baseDir    = runProcess.getParameter<std::string>("dirName");
  TString url = TString(argv[1]);
  TString outFileUrl(gSystem->BaseName(url));
  outFileUrl.ReplaceAll("_cfg.py","");
  if(mctruthmode!=0) { outFileUrl += "_filt"; outFileUrl += mctruthmode; }
  TString outdir=runProcess.getParameter<std::string>("outdir");
  TString outUrl( outdir );
  gSystem->Exec("mkdir -p " + outUrl);


  bool filterOnlyEE(false), filterOnlyMUMU(false), filterOnlyEMU(false);
  if(!isMC)
    {
      if(url.Contains("DoubleEle")) filterOnlyEE=true;
      if(url.Contains("DoubleMu"))  filterOnlyMUMU=true;
      if(url.Contains("MuEG"))      filterOnlyEMU=true;
    }
  bool isSingleMuPD(!isMC && url.Contains("SingleMu"));  
  bool isV0JetsMC(isMC && (url.Contains("DYJetsToLL_50toInf") || url.Contains("WJets")));
  bool isWGmc(isMC && url.Contains("WG"));
  bool isZGmc(isMC && url.Contains("ZG"));
  bool isMC_GG  = isMC && ( string(url.Data()).find("GG" )  != string::npos);
  bool isMC_VBF = isMC && ( string(url.Data()).find("VBF")  != string::npos);
  bool isMC_125OnShell = isMC && (mctruthmode==521);
  if(isMC_125OnShell) mctruthmode=125;
  bool isMC_ZZ  = isMC && ( string(url.Data()).find("TeV_ZZ")  != string::npos);
  bool isMC_WZ  = isMC && ( string(url.Data()).find("TeV_WZ")  != string::npos);

  TString outTxtUrl= outUrl + "/" + outFileUrl + ".txt";
  FILE* outTxtFile = NULL;
  if(!isMC)outTxtFile = fopen(outTxtUrl.Data(), "w");
  printf("TextFile URL = %s\n",outTxtUrl.Data());

  //tree info
  TString dirname = runProcess.getParameter<std::string>("dirName");

  //systematics
  bool runSystematics                        = runProcess.getParameter<bool>("runSystematics");
  std::vector<TString> varNames(1,"");
  if(runSystematics){
    varNames.push_back("_jerup");    varNames.push_back("_jerdown");
    varNames.push_back("_jesup");    varNames.push_back("_jesdown");  
    varNames.push_back("_umetup");   varNames.push_back("_umetdown");  
    varNames.push_back("_lesup");    varNames.push_back("_lesdown");  
    varNames.push_back("_puup");     varNames.push_back("_pudown");  
    varNames.push_back("_btagup");   varNames.push_back("_btagdown");
    if(isMC_ZZ)             { varNames.push_back("_zzptup");   varNames.push_back("_zzptdown");     }
    if(isMC_WZ)             { varNames.push_back("_wzptup");   varNames.push_back("_wzptdown");     }
    if(isMC_GG || isMC_VBF) { varNames.push_back("_lshapeup"); varNames.push_back("_lshapedown"); }
  }
  size_t nvarsToInclude=varNames.size();
  
  std::vector<std::string> allWeightsURL=runProcess.getParameter<std::vector<std::string> >("weightsFile");
  std::string weightsDir( allWeightsURL.size() ? allWeightsURL[0] : "");

  //  //shape uncertainties for dibosons
  //  std::vector<TGraph *> vvShapeUnc;
  //  if(isMC_ZZ || isMC_WZ)
  //    {
  //      TString weightsFile=weightsDir+"/zzQ2unc.root";
  //      TString dist("zzpt");
  //      if(isMC_WZ) { weightsFile.ReplaceAll("zzQ2","wzQ2"); dist.ReplaceAll("zzpt","wzpt"); }
  //      gSystem->ExpandPathName(weightsFile);
  //      TFile *q2UncF=TFile::Open(weightsFile);
  //      if(q2UncF){
  //	vvShapeUnc.push_back( new TGraph( (TH1 *)q2UncF->Get(dist+"_up") ) );
  //	vvShapeUnc.push_back( new TGraph( (TH1 *)q2UncF->Get(dist+"_down") ) );
  //	q2UncF->Close();
  //      }
  //    }
  
  
  //##############################################
  //########    INITIATING HISTOGRAMS     ########
  //##############################################
  SmartSelectionMonitor mon;

  //generator level control : add an underflow entry to make sure the histo is kept
  //((TH1F*)mon.addHistogram( new TH1F( "higgsMass_raw",     ";Higgs Mass [GeV];Events", 500,0,1500) ))->Fill(-1.0,0.0001);

  //event selection
  TH1F *h=(TH1F*) mon.addHistogram( new TH1F ("eventflow", ";;Events", 9,0,9) );
  h->GetXaxis()->SetBinLabel(1,"raw");
  h->GetXaxis()->SetBinLabel(2,"#geq 2 iso leptons");
  h->GetXaxis()->SetBinLabel(3,"|M-91|<15");
  h->GetXaxis()->SetBinLabel(4,"p_{T}>55");
  h->GetXaxis()->SetBinLabel(5,"3^{rd}-lepton veto");
  h->GetXaxis()->SetBinLabel(6,"b-veto"); 
  h->GetXaxis()->SetBinLabel(7,"#Delta #phi(jet,E_{T}^{miss})>0.5");
  h->GetXaxis()->SetBinLabel(8,"E_{T}^{miss}>80");

  //pu control
  mon.addHistogram( new TH1F( "nvtx",";Vertices;Events",50,0,50) ); 
  mon.addHistogram( new TH1F( "nvtxraw",";Vertices;Events",50,0,50) ); 
  mon.addHistogram( new TH1F( "rho",";#rho;Events",50,0,25) ); 

  //tau control
  mon.addHistogram( new TH1F( "leadtaupt",     ";Transverse momentum [GeV];Events", 50,0,500) );
  TH1 *htaus=mon.addHistogram( new TH1F("ntaus",  ";Tau multiplicity;Events",5,0,5) );
  for(int ibin=1; ibin<=htaus->GetXaxis()->GetNbins(); ibin++)
    {
      TString label("");
      if(ibin==h->GetXaxis()->GetNbins()) label +="#geq";
      else                                label +="=";
      label += (ibin-1);
      htaus->GetXaxis()->SetBinLabel(ibin,label);
    } 

  //lepton control
  mon.addHistogram( new TH1F( "leadpt",     ";Transverse momentum [GeV];Events", 50,0,500) );
  mon.addHistogram( new TH1F( "leadeta",    ";Pseudo-rapidity;Events", 50,0,2.6) );
  mon.addHistogram( new TH1F( "trailerpt",  ";Transverse momentum [GeV];Events", 50,0,500) );
  mon.addHistogram( new TH1F( "trailereta", ";Pseudo-rapidity;Events", 50,0,2.6) );
  mon.addHistogram( new TH1F( "zy",         ";Rapidity;Events", 50,0,3) );
  mon.addHistogram( new TH1F( "zmass",      ";Mass [GeV];Events", 100,40,250) );
  mon.addHistogram( new TH1F( "zpt",        ";Transverse momentum [GeV];Events",100,0,1500));
  mon.addHistogram( new TH1F( "qmass",      ";Mass [GeV];Events / (1 GeV)",100,76,106));
  mon.addHistogram( new TH1F( "qt",         ";Transverse momentum [GeV];Events / (1 GeV)",1500,0,1500));
  mon.addHistogram( new TH1F( "qtraw",      ";Transverse momentum [GeV];Events / (1 GeV)",1500,0,1500));

  //extra leptons in the event
  mon.addHistogram( new TH1F( "nextraleptons", ";Extra leptons;Events",4,0,4) );
  mon.addHistogram( new TH1F( "thirdleptonpt", ";Transverse momentum;Events", 50,0,500) );
  mon.addHistogram( new TH1F( "thirdleptoneta", ";Pseudo-rapidity;Events", 50,0,2.6) );
  mon.addHistogram( new TH1F( "thirdleptonmt", ";Transverse mass(3^{rd} lepton,E_{T}^{miss}) [GeV];Events", 50,0,500) );


  mon.addHistogram( new TH1F("csv",      ";Combined Secondary Vertex;Jets",50,0.,1.) );
  mon.addHistogram( new TH1F("csvb",     ";Combined Secondary Vertex;Jets",50,0.,1.) );
  mon.addHistogram( new TH1F("csvc",     ";Combined Secondary Vertex;Jets",50,0.,1.) );
  mon.addHistogram( new TH1F("csvothers",";Combined Secondary Vertex;Jets",50,0.,1.) );
  TH1 *hbtags=mon.addHistogram( new TH1F("nbtags",   ";b-tag multiplicity;Events",5,0,5) );
  TH1 *hbtagsJP=mon.addHistogram( new TH1F("nbtagsJP",   ";b-tag multiplicity;Events",5,0,5) );
  mon.addHistogram( new TH1F("leadjetpt",    ";Transverse momentum [GeV];Events",50,0,1000) );
  mon.addHistogram( new TH1F("trailerjetpt", ";Transverse momentum [GeV];Events",50,0,1000) );
  mon.addHistogram( new TH1F("fwdjeteta",    ";Pseudo-rapidity;Events",25,0,5) );
  mon.addHistogram( new TH1F("cenjeteta",       ";Pseudo-rapidity;Events",25,0,5) );
  Double_t mjjaxis[32];
  mjjaxis[0]=0.01;
  for(size_t i=1; i<20; i++)  mjjaxis[i]   =50*i;        //0-1000
  for(size_t i=0; i<5; i++)   mjjaxis[20+i]=1000+100*i; //1000-1500
  for(size_t i=0; i<=5; i++)   mjjaxis[25+i]=1500+300*i; //1500-5000  
  mjjaxis[31]=5000;
  mon.addHistogram( new TH1F("vbfmjj"       , ";Dijet invariant mass [GeV];Events",31,mjjaxis) );
  mon.addHistogram( new TH1F("vbfdphijj"    , ";Azimuthal angle difference;Events",20,0,3.5) );
  mon.addHistogram( new TH1F("vbfdetajj"    , ";Pseudo-rapidity span;Events",20,0,10) );
  TH1 *hjets=mon.addHistogram( new TH1F("njets",  ";Jet multiplicity;Events",5,0,5) );
  for(int ibin=1; ibin<=hjets->GetXaxis()->GetNbins(); ibin++)
    {
      TString label("");
      if(ibin==h->GetXaxis()->GetNbins()) label +="#geq";
      else                                label +="=";
      label += (ibin-1);
      hjets->GetXaxis()->SetBinLabel(ibin,label);
      hbtags->GetXaxis()->SetBinLabel(ibin,label);
      hbtagsJP->GetXaxis()->SetBinLabel(ibin,label);
    } 

  mon.addHistogram( new TH1F( "mindphijmet",  ";min #Delta#phi(jet,E_{T}^{miss});Events",40,0,4) );
  mon.addHistogram( new TH1F( "mindphijmetNM1",  ";min #Delta#phi(jet,E_{T}^{miss});Events",40,0,4) );
  mon.addHistogram( new TH1D( "balance",      ";E_{T}^{miss}/q_{T};Events", 25,0,2.5) );
  mon.addHistogram( new TH1D( "balanceNM1",   ";E_{T}^{miss}/q_{T};Events", 25,0,2.5) );
  mon.addHistogram( new TH1F( "axialmet",     ";Axial missing transvere energy [GeV];Events", 50,-100,400) );
  mon.addHistogram( new TH1F( "axialmetNM1",   ";Axial missing transvere energy [GeV];Events", 50,-100,400) );
  Double_t metaxis[]={0,5,10,15,20,25,30,35,40,45,50,55,60,70,80,90,100,125,150,175,200,250,300,400,500};
  Int_t nmetAxis=sizeof(metaxis)/sizeof(Double_t);
  mon.addHistogram( new TH1F( "met",          ";Missing transverse energy [GeV];Events",nmetAxis-1,metaxis) ); //50,0,1000) );
  mon.addHistogram( new TH1F( "metNM1",        ";Missing transverse energy [GeV];Events",nmetAxis-1,metaxis) ); //50,0,1000) );
  Double_t mtaxis[]={100,120,140,160,180,200,220,240,260,280,300,325,350,375,400,450,500,600,700,800,900,1000,2000};
  Int_t nmtAxis=sizeof(mtaxis)/sizeof(Double_t);
  mon.addHistogram( new TH1F( "mt"  ,         ";Transverse mass;Events",nmtAxis-1,mtaxis) );
  mon.addHistogram( new TH1F( "mtNM1"  ,       ";Transverse mass;Events",nmtAxis-1,mtaxis) );
  mon.addHistogram( new TH1F( "mtresponse",   ";Transverse mass response;Events", 100,0,2) );
  mon.addHistogram( new TH1F( "mtcheckpoint"  ,         ";Transverse mass [GeV];Events",160,150,1750) );
  mon.addHistogram( new TH1F( "metcheckpoint" ,         ";Missing transverse energy [GeV];Events",100,0,500) );


  //
  // STATISTICAL ANALYSIS
  //
  TH1F* Hoptim_systs     =  (TH1F*) mon.addHistogram( new TH1F ("optim_systs"    , ";syst;", nvarsToInclude,0,nvarsToInclude) ) ;
  for(size_t ivar=0; ivar<nvarsToInclude; ivar++)
    Hoptim_systs->GetXaxis()->SetBinLabel(ivar+1, varNames[ivar]);
    

  
  //##############################################
  //######## GET READY FOR THE EVENT LOOP ########
  //##############################################

  fwlite::ChainEvent ev(urls);
  const size_t totalEntries= ev.size();

  //MC normalization (to 1/pb)
  double xsecWeight = xsec/totalEntries;
  if(!isMC) xsecWeight=1.0;

  //jet energy scale and uncertainties 
  TString jecDir = runProcess.getParameter<std::string>("jecDir");
  gSystem->ExpandPathName(jecDir);
  FactorizedJetCorrector *jesCor        = utils::cmssw::getJetCorrector(jecDir,isMC);
  JetCorrectionUncertainty *totalJESUnc = new JetCorrectionUncertainty((jecDir+"/MC_Uncertainty_AK5PFchs.txt").Data());
  
  //muon energy scale and uncertainties
  MuScleFitCorrector *muCor=getMuonCorrector(jecDir,url);

  //lepton efficiencies
  LeptonEfficiencySF lepEff;

  //b-tagging: beff and leff must be derived from the MC sample using the discriminator vs flavor
  //the scale factors are taken as average numbers from the pT dependent curves see:
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2012_Data_and_MC_EPS13_prescript
  BTagSFUtil btsfutil;
  float beff(0.68), sfb(0.99), sfbunc(0.015);
  float leff(0.13), sfl(1.05), sflunc(0.12);

  //pileup weighting
  edm::LumiReWeighting* LumiWeights = NULL;
  utils::cmssw::PuShifter_t PuShifters;
  double PUNorm[] = {1,1,1};
  if(isMC){
          std::vector<double> dataPileupDistributionDouble = runProcess.getParameter< std::vector<double> >("datapileup");
          std::vector<float> dataPileupDistribution; for(unsigned int i=0;i<dataPileupDistributionDouble.size();i++){dataPileupDistribution.push_back(dataPileupDistributionDouble[i]);}
          std::vector<float> mcPileupDistribution;
          utils::getMCPileupDistributionFromMiniAOD(ev,dataPileupDistribution.size(), mcPileupDistribution);
          while(mcPileupDistribution.size()<dataPileupDistribution.size())  mcPileupDistribution.push_back(0.0);
          while(mcPileupDistribution.size()>dataPileupDistribution.size())dataPileupDistribution.push_back(0.0);
          gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE
          LumiWeights = new edm::LumiReWeighting(mcPileupDistribution,dataPileupDistribution);
          PuShifters=utils::cmssw::getPUshifters(dataPileupDistribution,0.05);
          utils::getPileupNormalization(mcPileupDistribution, PUNorm, LumiWeights, PuShifters);
  }
 
  gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE


  //higgs::utils::EventCategory eventCategoryInst(higgs::utils::EventCategory::EXCLUSIVE2JETSVBF); //jet(0,>=1)+vbf binning


  //##############################################
  //########           EVENT LOOP         ########
  //##############################################
  //loop on all the events
  printf("Progressing Bar     :0%%       20%%       40%%       60%%       80%%       100%%\n");
  printf("Scanning the ntuple :");
  int treeStep(totalEntries/50);
  //DuplicatesChecker duplicatesChecker;
  //int nDuplicates(0);
  for( size_t iev=0; iev<totalEntries; iev++){
      if(iev%treeStep==0){printf(".");fflush(stdout);}

       //##############################################   EVENT LOOP STARTS   ##############################################
       ev.to(iev); //load the event content from the EDM file
       //if(!isMC && duplicatesChecker.isDuplicate( ev.run, ev.lumi, ev.event) ) { nDuplicates++; continue; }

       //apply trigger and require compatibilitiy of the event with the PD
       edm::TriggerResultsByName tr = ev.triggerResultsByName("HLT");
       if(!tr.isValid())return false;

      bool eeTrigger          = utils::passTriggerPatterns(tr, "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
      bool muTrigger          = utils::passTriggerPatterns(tr, "HLT_IsoMu24_eta2p1_v*");
      bool mumuTrigger        = utils::passTriggerPatterns(tr, "HLT_Mu17_Mu8_v*", "HLT_Mu17_TkMu8_v*"); 
      bool emuTrigger         = utils::passTriggerPatterns(tr, "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*", "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
      if(filterOnlyEE)   { mumuTrigger=false; emuTrigger=false;  }
      if(filterOnlyMUMU) { eeTrigger=false;   emuTrigger=false;  }
      if(isSingleMuPD)   { eeTrigger=false;   emuTrigger=false;  if( muTrigger && !mumuTrigger) mumuTrigger=true; else mumuTrigger=false; }
      if(filterOnlyEMU)  { eeTrigger=false;   mumuTrigger=false; }

      if(!(eeTrigger || muTrigger || mumuTrigger || emuTrigger ))continue;  //ONLY RUN ON THE EVENTS THAT PASS OUR TRIGGERS

       //##############################################   EVENT PASSED THE TRIGGER   #######################################


       //load all the objects we will need to access
       reco::VertexCollection vtx;
       fwlite::Handle< reco::VertexCollection > vtxHandle; 
       vtxHandle.getByLabel(ev, "offlineSlimmedPrimaryVertices");
       if(vtxHandle.isValid()){ vtx = *vtxHandle;}


       double rho = 0;
       fwlite::Handle< double > rhoHandle;
       rhoHandle.getByLabel(ev, "fixedGridRhoFastjetAll");
       if(rhoHandle.isValid()){ rho = *rhoHandle;}

       reco::GenParticleCollection gen;
       fwlite::Handle< reco::GenParticleCollection > genHandle;
       genHandle.getByLabel(ev, "prunedGenParticles");
       if(genHandle.isValid()){ gen = *genHandle;}

       pat::MuonCollection muons;
       fwlite::Handle< pat::MuonCollection > muonsHandle;
       muonsHandle.getByLabel(ev, "slimmedMuons");
       if(muonsHandle.isValid()){ muons = *muonsHandle;}

       pat::ElectronCollection electrons;
       fwlite::Handle< pat::ElectronCollection > electronsHandle;
       electronsHandle.getByLabel(ev, "slimmedElectrons");
       if(electronsHandle.isValid()){ electrons = *electronsHandle;}

       pat::JetCollection jets;
       fwlite::Handle< pat::JetCollection > jetsHandle;
       jetsHandle.getByLabel(ev, "slimmedJets");
       if(jetsHandle.isValid()){ jets = *jetsHandle;}

       pat::PhotonCollection photons;
       fwlite::Handle< pat::PhotonCollection > photonsHandle;
       photonsHandle.getByLabel(ev, "slimmedPhotons");
       if(photonsHandle.isValid()){ photons = *photonsHandle;}
       
       pat::METCollection mets;
       fwlite::Handle< pat::METCollection > metsHandle;
       metsHandle.getByLabel(ev, "slimmedMETs");
       if(metsHandle.isValid()){ mets = *metsHandle;}
       LorentzVector met = mets[0].p4(); 

       pat::TauCollection taus;
       fwlite::Handle< pat::TauCollection > tausHandle;
       tausHandle.getByLabel(ev, "slimmedTaus");
       if(tausHandle.isValid()){ taus = *tausHandle;}

      if(isV0JetsMC){
         fwlite::Handle< LHEEventProduct > lheEPHandle;
         lheEPHandle.getByLabel(ev, "externalLHEProducer");

 	 mon.fillHisto("nup","",lheEPHandle->hepeup().NUP,1);
 	 if(lheEPHandle->hepeup().NUP>5) continue;
	 mon.fillHisto("nupfilt","",lheEPHandle->hepeup().NUP,1);
      }



      //
      // DERIVE WEIGHTS TO APPLY TO SAMPLE
      //

       //pileup weight
       float weight = 1.0;
       double TotalWeight_plus = 1.0;
       double TotalWeight_minus = 1.0;
       float puWeight(1.0);

       if(isMC){          
          int ngenITpu = 0;

          fwlite::Handle< std::vector<PileupSummaryInfo> > puInfoH;
          puInfoH.getByLabel(ev, "addPileupInfo");
          for(std::vector<PileupSummaryInfo>::const_iterator it = puInfoH->begin(); it != puInfoH->end(); it++){
             if(it->getBunchCrossing()==0)      { ngenITpu += it->getPU_NumInteractions(); }
          }

          puWeight          = LumiWeights->weight(ngenITpu) * PUNorm[0];
          weight            = xsecWeight*puWeight;
          TotalWeight_plus  = PuShifters[utils::cmssw::PUUP  ]->Eval(ngenITpu) * (PUNorm[2]/PUNorm[0]);
          TotalWeight_minus = PuShifters[utils::cmssw::PUDOWN]->Eval(ngenITpu) * (PUNorm[1]/PUNorm[0]);
      }


      //
      //
      // BELOW FOLLOWS THE ANALYSIS OF THE MAIN SELECTION WITH N-1 PLOTS
      //
      //

      //
      // LEPTON ANALYSIS
      //
      
      //start by merging electrons and muons
      std::vector<patUtils::GenericLepton> leptons;
      for(size_t l=0;l<electrons.size();l++){leptons.push_back(patUtils::GenericLepton(electrons[l]));}      
      for(size_t l=0;l<muons    .size();l++){leptons.push_back(patUtils::GenericLepton(muons    [l]));}      
      std::sort(leptons.begin(),   leptons.end(), utils::sort_CandidatesByPt);

      LorentzVector muDiff(0,0,0,0);
      std::vector<patUtils::GenericLepton> selLeptons, extraLeptons;
      for(size_t ilep=0; ilep<leptons.size(); ilep++)
	{
	  bool passKin(true),passId(true),passIso(true);
	  bool passLooseLepton(true), passSoftMuon(true), passSoftElectron(true), passVetoElectron(true);

	  int lid=leptons[ilep].pdgId();

	  //apply muon corrections
	  if(abs(lid)==13)
	    {
	      passSoftMuon=false;
	      if(muCor){
		TLorentzVector p4(leptons[ilep].px(),leptons[ilep].py(),leptons[ilep].pz(),leptons[ilep].energy());
		muCor->applyPtCorrection(p4 , lid<0 ? -1 :1 );
		if(isMC) muCor->applyPtSmearing(p4, lid<0 ? -1 : 1, false);
		muDiff -= leptons[ilep].p4();
                leptons[ilep].setP4(LorentzVector(p4.Px(),p4.Py(),p4.Pz(),p4.E() ) );
		muDiff += leptons[ilep].p4();
	      }
	    }

	  //no need for charge info any longer
	  lid=abs(lid);
	  TString lepStr( lid==13 ? "mu" : "e");

	  // don't want to mess with photon ID // //veto nearby photon (loose electrons are many times photons...)
	  // don't want to mess with photon ID // double minDRlg(9999.);
	  // don't want to mess with photon ID // for(size_t ipho=0; ipho<selPhotons.size(); ipho++)
	  // don't want to mess with photon ID //   minDRlg=TMath::Min(minDRlg,deltaR(leptons[ilep].p4(),selPhotons[ipho].p4()));
	  // don't want to mess with photon ID // if(minDRlg<0.1) continue;

	  //kinematics
	  float leta = fabs(lid==11 ?  leptons[ilep].el.superCluster()->eta() : leptons[ilep].eta());
	  if(leta> (lid==11 ? 2.5 : 2.4) )            passKin=false;
	  if(lid==11 && (leta>1.4442 && leta<1.5660)) passKin=false;
	  passLooseLepton &= passKin;
	  passSoftMuon    &= passKin;
	  if(lid==13){
	    if(leptons[ilep].pt()<10) passLooseLepton=false;
	    if(leptons[ilep].pt()<3)  passSoftMuon=false;
	  }
	  else if(lid==11){
	    if(leptons[ilep].pt()<10) passLooseLepton=false;
	  }
	  if(leptons[ilep].pt()<20) passKin=false;

          //Cut based identification 
          passId = lid==11?patUtils::passId(leptons[ilep].el, vtx[0], patUtils::llvvElecId::Tight) : patUtils::passId(leptons[ilep].mu, vtx[0], patUtils::llvvMuonId::Tight);

	  //isolation
	  passIso = lid==11?patUtils::passIso(leptons[ilep].el,  patUtils::llvvElecIso::Tight) : patUtils::passIso(leptons[ilep].mu,  patUtils::llvvMuonIso::Tight);
	  
	  if(passId && passIso && passKin)          selLeptons.push_back(leptons[ilep]);
	  else if(passLooseLepton || passSoftMuon)  extraLeptons.push_back(leptons[ilep]);

	}
        std::sort(selLeptons.begin(),   selLeptons.end(), utils::sort_CandidatesByPt);
        std::sort(extraLeptons.begin(), extraLeptons.end(), utils::sort_CandidatesByPt);
        LorentzVector recoMET = met - muDiff;

      //
      //JET/MET ANALYSIS
      //
      //add scale/resolution uncertainties and propagate to the MET      
      //utils::cmssw::updateJEC(jets,jesCor,totalJESUnc,rho,vtx.size(),isMC);  //FIXME if still needed
      //std::vector<LorentzVector> met=utils::cmssw::getMETvariations(recoMet,jets,selLeptons,isMC); //FIXME if still needed

      //select the jets
      pat::JetCollection selJets;
      int njets(0),nbtags(0),nbtagsJP(0);
      float mindphijmet(9999.);
      for(size_t ijet=0; ijet<jets.size(); ijet++){
	  if(jets[ijet].pt()<15 || fabs(jets[ijet].eta())>4.7 ) continue;

	  //mc truth for this jet
	  const reco::GenJet* genJet=jets[ijet].genJet();
	  TString jetType( genJet && genJet->pt()>0 ? "truejetsid" : "pujetsid" );
	  
	  //cross-clean with selected leptons and photons
	  float minDRlj(9999.),minDRlg(9999.);
          for(size_t ilep=0; ilep<selLeptons.size(); ilep++)
            minDRlj = TMath::Min( minDRlj, deltaR(jets[ijet],selLeptons[ilep]) );
	  // don't want to mess with photon ID // for(size_t ipho=0; ipho<selPhotons.size(); ipho++)
	  // don't want to mess with photon ID //   minDRlg = TMath::Min( minDRlg, deltaR(jets[ijet],selPhotons[ipho]) );
	  if(minDRlj<0.4 /*|| minDRlg<0.4*/) continue;
	  
	  //jet id
	  bool passPFloose = true; //FIXME --> Need to be updated according to te latest recipe;
	  float PUDiscriminant = jets[ijet].userFloat("pileupJetId:fullDiscriminant");
	  bool passLooseSimplePuId = true; //FIXME --> Need to be updated according to the latest recipe
	  if(jets[ijet].pt()>30){
	      mon.fillHisto(jetType,"",fabs(jets[ijet].eta()),0);
	      if(passPFloose)                        mon.fillHisto(jetType,"",fabs(jets[ijet].eta()),1);
	      if(passLooseSimplePuId)                mon.fillHisto(jetType,"",fabs(jets[ijet].eta()),2);
	      if(passPFloose && passLooseSimplePuId) mon.fillHisto(jetType,"",fabs(jets[ijet].eta()),3);
	  }
	  if(!passPFloose || !passLooseSimplePuId) continue;
	  selJets.push_back(jets[ijet]);
	  if(jets[ijet].pt()>30) {
	    njets++;
	    float dphijmet=fabs(deltaPhi(met.phi(), jets[ijet].phi()));
	    if(dphijmet<mindphijmet) mindphijmet=dphijmet;
	    if(fabs(jets[ijet].eta())<2.5){
	      bool hasCSVtag(jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags")>0.405);
	      //update according to the SF measured by BTV
	      if(isMC){
		  int flavId=jets[ijet].partonFlavour();
		  if(abs(flavId)==5)        btsfutil.modifyBTagsWithSF(hasCSVtag,sfb,beff);
		  else if(abs(flavId)==4)   btsfutil.modifyBTagsWithSF(hasCSVtag,sfb/5,beff);
		  else		            btsfutil.modifyBTagsWithSF(hasCSVtag,sfl,leff);
              }
	      nbtags   += hasCSVtag;
	    }
	  }
	}
      std::sort(selJets.begin(), selJets.end(), utils::sort_CandidatesByPt);

      //select the taus
      pat::TauCollection selTaus;
      int ntaus(0);
      for(size_t itau=0; itau<taus.size(); ++itau){
	pat::Tau& tau = taus[itau];
	if(tau.pt()<20. || fabs(tau.eta()) >2.3) continue;
	
	//	bool overlapWithLepton(false);
	//	for(int l1=0; l1<(int)selLeptons.size();++l1){
	//	  if(deltaR(tau, selLeptons[l1])<0.1){overlapWithLepton=true; break;}
	//	}
	//	if(overlapWithLepton) continue;
	
	//	if(!tau.isPFTau()) continue; // Only PFTaus
	//	if(tau.emFraction() >=2.) continue;
	
	if(!tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"))continue;
	if(!tau.tauID("againstMuonTight3"))continue; 
	if(!tau.tauID("againstElectronMediumMVA5"))continue;
	
	selTaus.push_back(tau);
	ntaus++;
      }
      std::sort(selTaus.begin(), selTaus.end(), utils::sort_CandidatesByPt);
      

      //
      // ASSIGN CHANNEL
      //
      std::vector<TString> chTags;
      int dilId(1);
      LorentzVector boson(0,0,0,0);
      if(selLeptons.size()==2)
	{
 	  for(size_t ilep=0; ilep<2; ilep++)
	    {
	      dilId *= selLeptons[ilep].pdgId();
	      int id(abs(selLeptons[ilep].pdgId()));
	      weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[ilep].pt(), selLeptons[ilep].eta(), id,  id ==11 ? "loose" : "loose" ).first : 1.0;
	      boson += selLeptons[ilep].p4();
	    }
     
	  //check the channel
	  if( abs(dilId)==121 && eeTrigger)   chTags.push_back("ee");
	  if( abs(dilId)==169 && mumuTrigger) chTags.push_back("mumu"); 
	  if( abs(dilId)==143 && emuTrigger) chTags.push_back("emu"); 
	}

      // keep in mind the eventCategory thingy for more refined categorization // TString evCat=eventCategoryInst.GetCategory(selJets,boson);
      std::vector<TString> tags(1,"all");
      for(size_t ich=0; ich<chTags.size(); ich++){
	tags.push_back( chTags[ich] );
	//tags.push_back( chTags[ich]+evCat );
      }

      mon.fillHisto("eventflow",  tags,0,weight);
      if(chTags.size()==0) continue;

      //
      // BASELINE SELECTION
      //
      bool passMass(fabs(boson.mass()-91)<15);
      bool passQt(boson.pt()>55);
      bool passThirdLeptonVeto( selLeptons.size()==2 && extraLeptons.size()==0 );
      bool passBtags(nbtags==0);
      bool passMinDphijmet( njets==0 || mindphijmet>0.5);

      mon.fillHisto("eventflow",  tags,1,weight);
      mon.fillHisto("nvtxraw",  tags,vtx.size(),weight/puWeight);
      mon.fillHisto("nvtx",  tags,vtx.size(),weight);
      mon.fillHisto("rho",  tags,rho,weight);
      mon.fillHisto("leadpt",      tags,selLeptons[0].pt(),weight); 
      mon.fillHisto("trailerpt",   tags,selLeptons[1].pt(),weight); 
      mon.fillHisto("leadeta",     tags,fabs(selLeptons[0].eta()),weight); 
      mon.fillHisto("trailereta",  tags,fabs(selLeptons[1].eta()),weight); 
      
      mon.fillHisto("ntaus", tags, ntaus,weight);
      if(ntaus>0) mon.fillHisto("leadtaupt", tags, selTaus[0].pt(),weight);
      
      mon.fillHisto("zmass", tags,boson.mass(),weight); 
      mon.fillHisto("zy",    tags,fabs(boson.Rapidity()),weight); 

      if(passMass){

	mon.fillHisto("eventflow",tags, 2,weight);
	mon.fillHisto("zpt",      tags, boson.pt(),weight);

	//these two are used to reweight photon -> Z, the 3rd is a control
	mon.fillHisto("qt",       tags, boson.pt(),weight,true); 
	///	mon.fillHisto("qtraw",    tags, boson.pt(),weight/triggerPrescale,true); 

	if(passQt){
	  mon.fillHisto("eventflow",tags,3,weight);
	  int nExtraLeptons((selLeptons.size()-2)+extraLeptons.size());
	  mon.fillHisto("nextraleptons",tags,nExtraLeptons,weight);
	  if(nExtraLeptons>0){
	    LorentzVector thirdLepton(selLeptons.size()>2 ?  selLeptons[1].p4() : extraLeptons[0].p4());
	    double dphi=fabs(deltaPhi(thirdLepton.phi(),met.phi()));
	    double mt=TMath::Sqrt(2*thirdLepton.pt()*met.pt()*(1-TMath::Cos(dphi)));
	    mon.fillHisto("thirdleptonpt",tags,thirdLepton.pt(),weight);
	    mon.fillHisto("thirdleptoneta",tags,fabs(thirdLepton.eta()),weight);
	    mon.fillHisto("thirdleptonmt",tags,mt,weight);
	  }
	  if(passThirdLeptonVeto){
	    
	    mon.fillHisto("eventflow",tags,4,weight);
	    for(size_t ijet=0; ijet<selJets.size(); ijet++){
	      if(selJets[ijet].pt()<30 || fabs(selJets[ijet].eta())>2.5) continue;

	      float csv(selJets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"));
	      mon.fillHisto( "csv",tags,csv,weight);
	      if(!isMC) continue;
	      int flavId=selJets[ijet].partonFlavour();
	      TString jetFlav("others");
	      if(abs(flavId)==5)      jetFlav="b";
	      else if(abs(flavId)==4) jetFlav="c";
	      mon.fillHisto( "csv"+jetFlav,tags,csv,weight);
	    }
	    mon.fillHisto( "nbtags",tags,nbtags,weight);
	    mon.fillHisto( "nbtagsJP",tags,nbtagsJP,weight);
	    
	    if(passBtags){
	      mon.fillHisto("eventflow",tags,5,weight);

	      //include photon prediction from this point forward
	      //would require looping tag by tag as weights are category-specific

	    }
	  }
	}        
      }


      //
      // HISTOS FOR STATISTICAL ANALYSIS (include systematic variations)
      //
      //Fill histogram for posterior optimization, or for control regions
      for(size_t ivar=0; ivar<nvarsToInclude; ivar++){
	float iweight = weight;                            //nominal
	
	//energy scale/resolution
	bool varyJesUp( varNames[ivar]=="_jesup" );
	bool varyJesDown( varNames[ivar]=="_jesdown" );
	bool varyJerUp( varNames[ivar]=="_jerup" );
	bool varyJerDown( varNames[ivar]=="_jerdown" );
	bool varyUmetUp( varNames[ivar]=="_umetup" );
	bool varyUmetDown( varNames[ivar]=="_umetdown" );
	bool varyLesUp( varNames[ivar]=="_lesup" );
	bool varyLesDown( varNames[ivar]=="_lesdown" );
		
	//pileup variations
	if(varNames[ivar]=="_puup") iweight *=TotalWeight_plus;
	if(varNames[ivar]=="_pudown") iweight *=TotalWeight_minus;
	
	//btag
	bool varyBtagUp( varNames[ivar]=="_btagup" );
	bool varyBtagDown( varNames[ivar]=="_btagdown" );
	
	//Here were the Q^2 variations on VV pT spectum


	//recompute MET/MT if JES/JER was varied
	LorentzVector    zvv = mets[0].p4();
	//FIXME
	//	if(varyJesUp)    zvv = mets[0].shiftedP4(pat::MET::METUncertainty::JetEnUp);
	//	if(varyJesDown)  zvv = mets[0].shiftedP4(pat::MET::METUncertainty::JetEnDown);
	//	if(varyJerUp)    zvv = mets[0].shiftedP4(pat::MET::METUncertainty::JetResUp);
	//	if(varyJerDown)  zvv = mets[0].shiftedP4(pat::MET::METUncertainty::JetResDown);
	//	if(varyUmetUp)   zvv = mets[0].shiftedP4(pat::MET::METUncertainty::UnclusteredEnUp);
	//	if(varyUmetDown) zvv = mets[0].shiftedP4(pat::MET::METUncertainty::UnclusteredEnDown);
	//	if(varyLesUp)    zvv = met[utils::cmssw::LESUP]; //FIXME  must vary all leptons separately: MuonEnUp/MuonEnDown/ElectronEnUp/ElectronEnDown/TauEnUp/TauEnDown
	//	if(varyLesDown)  zvv = met[utils::cmssw::LESDOWN];
	
        pat::JetCollection tightVarJets;
	bool passLocalBveto(passBtags);
 	for(size_t ijet=0; ijet<jets.size(); ijet++){
	  
	  float eta=jets[ijet].eta();
	  if( fabs(eta)>4.7 ) continue;
	  float pt=jets[ijet].pt();
          //FIXME
//	  if(varyJesUp)    pt=jets[ijet].getVal("jesup");
//	  if(varyJesDown)  pt=jets[ijet].getVal("jesdown");
//	  if(varyJerUp)    pt=jets[ijet].getVal("jerup");
//	  if(varyJerDown)  pt=jets[ijet].getVal("jerdown");
	  if(pt<30) continue;

	  //cross-clean with selected leptons and photons
	  double minDRlj(9999.),minDRlg(9999.);
          for(size_t ilep=0; ilep<selLeptons.size(); ilep++)
            minDRlj = TMath::Min( minDRlj, deltaR(jets[ijet].p4(),selLeptons[ilep].p4()) );
	  // don't want to mess with photon ID // for(size_t ipho=0; ipho<selPhotons.size(); ipho++)
	  // don't want to mess with photon ID //   minDRlg = TMath::Min( minDRlg, deltaR(jets[ijet].p4(),selPhotons[ipho].p4()) );
	  if(minDRlj<0.4 /*|| minDRlg<0.4*/) continue;
	  
	  //jet id
	  bool passPFloose = true;//FIXME
	  int simplePuId = true;//FIXME
	  bool passLooseSimplePuId = true;//FIXME
	  if(!passPFloose || !passLooseSimplePuId) continue;
	 
	  //jet is selected
	  tightVarJets.push_back(jets[ijet]);

	  //check b-tag
	  if(pt<30 || fabs(eta)>2.5) continue;
	  if(!isMC) continue;
	  if(!varyBtagUp && !varyBtagDown) continue;
	  int flavId=jets[ijet].partonFlavour();
	  bool hasCSVtag (jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags")>0.405);
 	  if(varyBtagUp) {
	    if(abs(flavId)==5)        btsfutil.modifyBTagsWithSF(hasCSVtag,sfb+sfbunc,beff);
	    else if(abs(flavId)==4)   btsfutil.modifyBTagsWithSF(hasCSVtag,sfb/5+2*sfbunc,beff);
	    else		      btsfutil.modifyBTagsWithSF(hasCSVtag,sfl+sflunc,leff);
	  }
 	  else if(varyBtagDown) {
	    if(abs(flavId)==5)        btsfutil.modifyBTagsWithSF(hasCSVtag,sfb-sfbunc,beff);
	    else if(abs(flavId)==4)   btsfutil.modifyBTagsWithSF(hasCSVtag,sfb/5-2*sfbunc,beff);
	    else		      btsfutil.modifyBTagsWithSF(hasCSVtag,sfl-sflunc,leff);
 	  }
	  passLocalBveto |= hasCSVtag;
 	}
	
	bool isZsideBand    ( (boson.mass()>40  && boson.mass()<70) || (boson.mass()>110 && boson.mass()<200) );
	bool isZsideBandPlus( (boson.mass()>110 && boson.mass()<200) );
 	bool passPreselection                 (passMass && passQt && passThirdLeptonVeto && passMinDphijmet && passLocalBveto);
 	bool passPreselectionMbvetoMzmass     (            passQt && passThirdLeptonVeto && passMinDphijmet                  );          
	
 	//re-assign the event category to take migrations into account
	// 	TString evCat  = eventCategoryInst.GetCategory(tightVarJets,boson);
	for(size_t ich=0; ich<chTags.size(); ich++){
	  
	  TString tags_full=chTags[ich];//+evCat;
	  float chWeight(iweight);

	  //update weight and mass for photons
	  LorentzVector iboson(boson);
	  
	  //updet the transverse mass
	  float mt =higgs::utils::transverseMass(iboson,zvv,true);

	}
      }
  }
  printf("\n"); 
  
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

  if(outTxtFile)fclose(outTxtFile);
}  







