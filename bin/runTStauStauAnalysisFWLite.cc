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

#include "TauAnalysis/CandidateTools/interface/NSVfitStandaloneAlgorithm.h" //for svfit


#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
//#include "UserCode/llvv_fwk/interface/DataEventSummaryHandler.h"
#include "UserCode/llvv_fwk/interface/llvvObjects.h"
#include "UserCode/llvv_fwk/interface/TMVAUtils.h"
#include "UserCode/llvv_fwk/interface/LeptonEfficiencySF.h"
#include "UserCode/llvv_fwk/interface/PDFInfo.h"
#include "UserCode/llvv_fwk/interface/MuScleFitCorrector.h"
#include "UserCode/llvv_fwk/interface/GammaWeightsHandler.h"
#include "UserCode/llvv_fwk/interface/BtagUncertaintyComputer.h"
#include "UserCode/llvv_fwk/interface/GammaWeightsHandler.h"


#include "TROOT.h"
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
#include "TGraphErrors.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/shared_ptr.hpp>
#include <Math/VectorUtil.h>

/*****************************************************************************/
/* Return Codes:                                                             */
/*   0 - Everything OK                                                       */
/*   1 - Missing parameters_cfg.py configuration file                        */
/*****************************************************************************/

int main(int argc, char* argv[])
{
  /***************************************************************************/
  /*                          Global Initialization                          */
  /***************************************************************************/
  // Check arguments:
  if(argc<2)
    std::cout << "Usage: " << argv[0] << " parameters_cfg.py" << std::endl, exit(1);

  // Load FWLite:
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

  // Configure the process (aka "Load parameters from configuration file")
  const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");
  bool isMC = runProcess.getParameter<bool>("isMC");
  double xsec = runProcess.getParameter<double>("xsec");
  int mctruthmode = runProcess.getParameter<int>("mctruthmode"); // What is this one for?
  std::vector<std::string> urls = runProcess.getParameter<std::vector<std::string> >("input");
  std::string baseDir = runProcess.getParameter<std::string>("dirName");
  std::string outdir = runProcess.getParameter<std::string>("outdir");
  std::string jecDir = runProcess.getParameter<std::string>("jecDir");
  bool runSystematics = runProcess.getParameter<bool>("runSystematics");
  bool saveSummaryTree = runProcess.getParameter<bool>("saveSummaryTree");

  // Hardcoded configs
  bool runLoosePhotonSelection(false); // Specific to ZHTauTau?
  bool Cut_tautau_MVA_iso(true); // Specific to ZHTauTau?
  bool examineThisEvent(false); // ?
  double minElPt        = 35;
  double maxElEta       = 2.5;
  double ECALGap_MinEta = 1.4442;
  double ECALGap_MaxEta = 1.5660;
  double minMuPt        = 30;
  double maxMuEta       = 2.4;
  double minTauPt       = 20;
  double maxTauEta      =  2.3;
  double minJetPt       = 30;
  double maxJetEta      =  2.5;
  double minTauJetPt    = 20;
  double maxTauJetEta   =  2.5;

  // Lepton Efficiencies
  LeptonEfficiencySF lepEff;

  std::vector<std::string>  weightsFile = runProcess.getParameter<std::vector<std::string> >("weightsFile");
  TString wFile("");
  if(weightsFile.size())
    wFile = TString(gSystem->ExpandPathName(weightsFile[0].c_str()));

  // B-tag efficiencies
  std::map<std::pair<TString,TString>, std::pair<TGraphErrors*,TGraphErrors*> > btagEffCorr;
  if(weightsFile.size() && isMC)
  {
    TString btagEffCorrUrl(wFile); btagEffCorrUrl += "/btagEff.root";
    gSystem->ExpandPathName(btagEffCorrUrl);
    TFile *btagF=TFile::Open(btagEffCorrUrl);
    if(btagF!=0 && !btagF->IsZombie())
    {
      TList *dirs=btagF->GetListOfKeys();
      for(int itagger=0; itagger < dirs->GetEntries(); itagger++)
      {
        TString iDir(dirs->At(itagger)->GetName());
        btagEffCorr[std::pair<TString,TString>(iDir,"b")]    = std::pair<TGraphErrors*,TGraphErrors*>((TGraphErrors*) btagF->Get(iDir+"/beff"),(TGraphErrors*) btagF->Get(iDir+"/sfb"));
        btagEffCorr[std::pair<TString,TString>(iDir,"c")]    = std::pair<TGraphErrors*,TGraphErrors*>((TGraphErrors*) btagF->Get(iDir+"/ceff"),(TGraphErrors*) btagF->Get(iDir+"/sfc"));
        btagEffCorr[std::pair<TString,TString>(iDir,"udsg")] = std::pair<TGraphErrors*,TGraphErrors*>((TGraphErrors*) btagF->Get(iDir+"/udsgeff"),(TGraphErrors*) btagF->Get(iDir+"/sfudsg"));
      }
      std::cout << btagEffCorr.size() << " b-tag correction factors have been read" << std::endl;
    }
  }


  // Setting Up
  gSystem->Exec(("mkdir -p " + outdir).c_str());
  std::string url = urls[0];
  std::string outFileUrl(gSystem->BaseName(url.c_str()));
  while(outFileUrl.find(".root", 0) != std::string::npos)
    outFileUrl.replace(outFileUrl.find(".root", 0), 5, "");
  if(mctruthmode != 0)
  {
    outFileUrl += "_filt";
    outFileUrl += mctruthmode;
  }
  std::string outUrl = outdir;

  bool isSingleMuPD(!isMC && url.Contains("SingleMu"));



  /***************************************************************************/
  /*                         Initializing Histograms                         */
  /***************************************************************************/
  SmartSelectionMonitor mon;
  TH1F *cutFlow = (TH1F*)mon.addHistogram(new TH1F("cutFlow", ";;Events", 20, 0, 20));
  cutFlow->GetXaxis()->SetBinLabel(1, "All");
  cutFlow->GetXaxis()->SetBinLabel(2, "HLT");
  cutFlow->GetXaxis()->SetBinLabel(3, "> 1l");
  cutFlow->GetXaxis()->SetBinLabel(4, "B-veto");
  cutFlow->GetXaxis()->SetBinLabel(5, "> 1#tau");
  TH1F *fractions = (TH1F*)mon.addHistogram(new TH1F("fractions", ";;Events", 20, 0, 20));
  fractions->GetXaxis()->SetBinLabel(1, "All");
  fractions->GetXaxis()->SetBinLabel(2, "HLT");
  fractions->GetXaxis()->SetBinLabel(3, "> 1l");
  fractions->GetXaxis()->SetBinLabel(4, "B-veto");
  fractions->GetXaxis()->SetBinLabel(5, "> 1#tau");
  // ...
  // TH2D* hist = (TH2D*)mon.addHistogram(...);

  mon.addHistogram(new TH1F("nup",     ";NUP;Events", 10, 0, 10));
  //mon.addHistogram(new TH1F("nupfilt", ";NUP;Events", 10, 0, 10));

  // PU
  mon.addHistogram(new TH1F("nvtx",    ";Vertices;Events",       50, -0.5, 49.5));
  mon.addHistogram(new TH1F("nvtxraw", ";Vertices;Events",       50, -0.5, 49.5));
  mon.addHistogram(new TH1F("rho",     ";#rho;Events",           50,  0,   25));
  mon.addHistogram(new TH1F("rho25",   ";#rho(#eta<2.5);Events", 50,  0,   25));

  // Leptons
  mon.addHistogram(new TH1F("nlep",       ";nlep;Events",       10,  0,   10));
  mon.addHistogram(new TH1F("leadpt",     ";p_{T}^{l};Events", 100,  0,  100));
  mon.addHistogram(new TH1F("leadeta",    ";#eta^{l};Events",   50, -2.6,  2.6));
  mon.addHistogram(new TH1F("leadcharge", ";q^{l};Events",       5, -2,    2));
  TH1F *leptonCutFlow   = (TH1F*)mon.addHistogram(new TH1F("leptonCutFlow",   ";;Leptons", 4, 0, 4));
  TH1F *leptonFractions = (TH1F*)mon.addHistogram(new TH1F("leptonFractions", ";;Leptons", 4, 0, 4));
  leptonCutFlow->GetXaxis()->SetBinLabel(1, "All");
  leptonFractions->GetXaxis()->SetBinLabel(1, "All");
  leptonCutFlow->GetXaxis()->SetBinLabel(2, "ID");
  leptonFractions->GetXaxis()->SetBinLabel(2, "ID");
  leptonCutFlow->GetXaxis()->SetBinLabel(3, "Kin");
  leptonFractions->GetXaxis()->SetBinLabel(3, "Kin");
  leptonCutFlow->GetXaxis()->SetBinLabel(4, "Iso");
  leptonFractions->GetXaxis()->SetBinLabel(4, "Iso");

  // Lepton Isolation
  mon.addHistogram(new TH1F("isomu",  "RelIso(#mu);;Leptons", 100, -0.5, 9.5));
  mon.addHistogram(new TH1F("isoele", "RelIso(ele);;Leptons", 100, -0.5, 9.5));

  // Taus
  mon.addHistogram(new TH1F("ntaus",         ";ntaus;Events",         10, -0.5,  9.5));
  mon.addHistogram(new TH1F("tauleadpt",     ";p_{T}^{#tau};Events", 100,  0,  100));
  mon.addHistogram(new TH1F("tauleadeta",    ";#eta^{#tau};Events",   50, -2.6,  2.6));
  mon.addHistogram(new TH1F("tauleadcharge", ";q^{#tau};Events",       5, -2,    2));
  mon.addHistogram(new TH1F("taudz",         ";dz^{#tau};Events",     50,  0,    2));
  mon.addHistogram(new TH1F("tauvz",         ";vz^{#tau};Events",     50,  0,    2));
  mon.addHistogram(new TH1F("tauleademfrac", ";emf^{#tau};Events",    50,  0,    5));
  TH1F *tauCutFlow   = (TH1F*)mon.addHistogram(new TH1F("tauCutFlow",   ";;#tau", 6, 0, 6));
  TH1F *tauFractions = (TH1F*)mon.addHistogram(new TH1F("tauFractions", ";;#tau", 6, 0, 6));
  tauCutFlow->GetXaxis()->SetBinLabel(1, "All");
  tauFractions->GetXaxis()->SetBinLabel(1, "All");
  tauCutFlow->GetXaxis()->SetBinLabel(2, "PF");
  tauFractions->GetXaxis()->SetBinLabel(2, "PF");
  tauCutFlow->GetXaxis()->SetBinLabel(3, "ID");
  tauFractions->GetXaxis()->SetBinLabel(3, "ID");
  tauCutFlow->GetXaxis()->SetBinLabel(4, "Quality");
  tauFractions->GetXaxis()->SetBinLabel(4, "Quality");
  tauCutFlow->GetXaxis()->SetBinLabel(5, "Kin");
  tauFractions->GetXaxis()->SetBinLabel(5, "Kin");
  tauCutFlow->GetXaxis()->SetBinLabel(6, "Iso");
  tauFractions->GetXaxis()->SetBinLabel(6, "Iso");
  TH1F *tauID = (TH1F*)mon.addHistogram(new TH1F("tauID", ";;#tau", 6, 0, 6));
  tauID->GetXaxis()->SetBinLabel(1, "All");
  tauID->GetXaxis()->SetBinLabel(2, "Not e");
  tauID->GetXaxis()->SetBinLabel(3, "Not #mu");
  tauID->GetXaxis()->SetBinLabel(4, "Not decay mode");
  tauID->GetXaxis()->SetBinLabel(5, "Not medium comb iso");
  tauID->GetXaxis()->SetBinLabel(6, "Overall ID");

  // Jets
  mon.addHistogram(new TH1F("njets", ";njets;Events", 6, 0, 6));
  mon.addHistogram(new TH1F("nbjets", ";njets;Events", 6, 0, 6));
  mon.addHistogram(new TH1F("jetleadpt", ";p_{T}^{jet};Events", 50, 0, 500));
  mon.addHistogram(new TH1F("jetleadeta", ";#eta^{jet};Events", 50, -2.6, 2.6));
  mon.addHistogram(new TH1F("jetcsv", ";csv;jets", 50, 0, 1));
  TH1F *jetCutFlow   = (TH1F*)mon.addHistogram(new TH1F("jetCutFlow",   ";;jets", 6, 0, 6));
  TH1F *jetFractions = (TH1F*)mon.addHistogram(new TH1F("jetFractions", ";;jets", 6, 0, 6));
  jetCutFlow->GetXaxis()->SetBinLabel(1, "All");
  jetFractions->GetXaxis()->SetBinLabel(1, "All");
  jetCutFlow->GetXaxis()->SetBinLabel(2, "PF Loose");
  jetFractions->GetXaxis()->SetBinLabel(2, "PF Loose");
  jetCutFlow->GetXaxis()->SetBinLabel(3, "Pre-Selection");
  jetFractions->GetXaxis()->SetBinLabel(3, "Pre-Selection");
  jetCutFlow->GetXaxis()->SetBinLabel(4, "ID");
  jetFractions->GetXaxis()->SetBinLabel(4, "ID");
  jetCutFlow->GetXaxis()->SetBinLabel(5, "Iso");
  jetFractions->GetXaxis()->SetBinLabel(5, "Iso");
  jetCutFlow->GetXaxis()->SetBinLabel(6, "Kin");
  jetFractions->GetXaxis()->SetBinLabel(6, "Kin");



  /***************************************************************************/
  /*                          Prepare for Event Loop                         */
  /***************************************************************************/
  fwlite::ChainEvent ev(urls);
  const Int_t totalEntries = ev.size();

  // MC normalization to 1/pb
  double nInitEvent = 1.;
  if(isMC)
    nInitEvent = (double) utils::getMergeableCounterValue(urls, "startCounter");
  double xsecWeight = xsec/nInitEvent;
  if(!isMC)
    xsecWeight = 1.;

  // Jet Energy Scale and Uncertainties
  jecDir = gSystem->ExpandPathName(jecDir.c_str());
  FactorizedJetCorrector *jesCor = utils::cmssw::getJetCorrector(jecDir, isMC);
  JetCorrectionUncertainty *totalJESUnc = new JetCorrectionUncertainty((jecDir+"/MC_Uncertainty_AK5PFchs.txt"));

  // Muon Energy Scale and Uncertainties
  MuScleFitCorrector* muCor = getMuonCorrector(jecDir, url);

  // Pileup Weighting: Based on vtx (For now?)
  edm::LumiReWeighting* LumiWeights = NULL;
  utils::cmssw::PuShifter_t PuShifters;
  double PUNorm[] = {1, 1, 1};
  if(isMC)
  {
    std::vector<double> dataPileupDistributionDouble = runProcess.getParameter<std::vector<double> >("datapileup");
    std::vector<float> dataPileupDistribution;
      for(unsigned int i = 0; i < dataPileupDistributionDouble.size(); ++i)
        dataPileupDistribution.push_back(dataPileupDistributionDouble[i]);
    std::vector<float> mcPileupDistribution;

    utils::getMCPileupDistribution(ev, dataPileupDistribution.size(), mcPileupDistribution);
    while(mcPileupDistribution.size() < dataPileupDistribution.size())  mcPileupDistribution.push_back(0.0);
    while(mcPileupDistribution.size() > dataPileupDistribution.size())  dataPileupDistribution.push_back(0.0);

    LumiWeights= new edm::LumiReWeighting(mcPileupDistribution, dataPileupDistribution);
    PuShifters=utils::cmssw::getPUshifters(dataPileupDistribution, 0.05);
    utils::getPileupNormalization(mcPileupDistribution, PUNorm, LumiWeights, PuShifters);
  }

  gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE



  /***************************************************************************/
  /*                               Event Loop                                */
  /***************************************************************************/
  std::cout << "       Progress Bar:0%      20%       40%       60%       80%      100%" << std::endl;
  std::cout << "Scanning the ntuple:";

  DuplicatesChecker duplicatesChecker;
  int nDuplicates(0);
  int step(totalEntries/50);

  // Redirect stdout and stderr to a temporary buffer, then output buffer after event loop
  std::ostream myCout(std::cout.rdbuf());
  std::stringstream buffer;
  std::streambuf *coutbuf = std::cout.rdbuf();
  std::streambuf *cerrbuf = std::cerr.rdbuf();
  std::cout.rdbuf(buffer.rdbuf());
  std::cerr.rdbuf(buffer.rdbuf());

  // Loop on events
  for(int iev = 0; iev < totalEntries; ++iev)
  {
    if(iev%step == 0)
      myCout << "_" << std::flush;

    // Prepare tags to fill the histograms
    std::vector<TString> chTags;
    chTags.push_back("all");

    // Load the event content from tree
    ev.to(iev);


    /****     Get information/collections from the event     ****/
    // Number of vertexes
    int nvtx = 0;
    fwlite::Handle<int> nvtxHandle;
    nvtxHandle.getByLabel(ev, "llvvObjectProducersUsed", "nvtx");
    if(nvtxHandle.isValid()) nvtx = *nvtxHandle;

    // Collection of generated particles
    fwlite::Handle<llvvGenEvent> genEventHandle;
    genEventHandle.getByLabel(ev, "llvvObjectProducersUsed");
    if(!genEventHandle.isValid())
    {
      std::cout << "llvvGenEvent Object NotFound" << std::endl;
      continue;
    }
    llvvGenEvent genEv = *genEventHandle;

    fwlite::Handle<llvvGenParticleCollection> genPartCollHandle;
    genPartCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
    if(!genPartCollHandle.isValid())
    {
      std::cout << "llvvGenParticleCollection Object NotFound" << std::endl;
      continue;
    }
    llvvGenParticleCollection gen = *genPartCollHandle;

    // Collection of leptons
    fwlite::Handle<llvvLeptonCollection> leptonCollHandle;
    leptonCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
    if(!leptonCollHandle.isValid())
    {
      std::cout << "llvvLeptonCollection Object NotFound" << std::endl;
      continue;
    }
    llvvLeptonCollection leptons = *leptonCollHandle;

    // Electron Information Collection
    fwlite::Handle<llvvElectronInfoCollection> electronInfoCollHandle;
    electronInfoCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
    if(!electronInfoCollHandle.isValid())
    {
      std::cout << "llvvElectronInfoCollection Object NotFound" << std::endl;
      continue;
    }

    // Muon Information Collection
    fwlite::Handle<llvvMuonInfoCollection> muonInfoCollHandle;
    muonInfoCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
    if(!muonInfoCollHandle.isValid())
    {
      std::cout << "llvvMuonInfoCollection Object NotFound" << std::endl;
      continue;
    }

    // Tau Collection
    fwlite::Handle<llvvTauCollection> tauCollHandle;
    tauCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
    if(!tauCollHandle.isValid())
    {
      std::cout << "llvvLeptonCollection Object NotFound" << std::endl;
      continue;
    }
    llvvTauCollection taus = *tauCollHandle;

    // Jet Collection
    fwlite::Handle<llvvJetCollection> jetCollHandle;
    jetCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
    if(!jetCollHandle.isValid())
    {
      std::cout << "llvvJetCollection Object NotFound" << std::endl;
      continue;
    }
    llvvJetExtCollection jets;
    //for(unsigned int i = 0; i < jetCollHandle->size(); ++i)
      //jets.push_back(llvvJetExt((*jetCollHandle)[i]));
    for(auto i = jetCollHandle->begin(); i != jetCollHandle->end(); ++i)
      jets.push_back(llvvJetExt(*i));

    // MET Collection
    fwlite::Handle<llvvMet> metHandle;
    metHandle.getByLabel(ev, "llvvObjectProducersUsed", "pfMETPFlow");
    if(!metHandle.isValid())
    {
      std::cout << "llvvMet Object NotFound" << std::endl;
      continue;
    }
    llvvMet met = *metHandle;

    // Trigger Bits
    fwlite::Handle<std::vector<bool> > triggerBitsHandle;
    triggerBitsHandle.getByLabel(ev, "llvvObjectProducersUsed", "triggerBits");
    if(!triggerBitsHandle.isValid())
    {
      std::cout << "triggerBits Object NotFound" << std::endl;
      continue;
    }
    std::vector<bool> triggerBits = *triggerBitsHandle;

    // Trigger Prescales
    fwlite::Handle<std::vector<int> > triggerPrescalesHandle;
    triggerPrescalesHandle.getByLabel(ev, "llvvObjectProducersUsed", "triggerPrescales");
    if(!triggerPrescalesHandle.isValid())
    {
      std::cout << "triggerPrescales Object NotFound" << std::endl;
      continue;
    }
    std::vector<int> triggerPrescales = *triggerPrescalesHandle;

    // Rho
    fwlite::Handle<double> rhoHandle;
    rhoHandle.getByLabel(ev, "kt6PFJets", "rho");
    if(!rhoHandle.isValid())
    {
      std::cout << "rho Object NotFound" << std::endl;
      continue;
    }
    double rho = *rhoHandle;

    // Rho25
    fwlite::Handle<double> rho25Handle;
    rho25Handle.getByLabel(ev, "kt6PFJetsCentral", "rho");
    if(!rho25Handle.isValid())
    {
      std::cout << "rho25 Object NotFound" << std::endl;
      continue;
    }
    double rho25 = *rho25Handle;


    /****          Sort events acording to HLT Path          ****/
    bool singleETrigger  = triggerBits[13]; // HLT_Ele27_WP80_v*
    bool singleMuTrigger = triggerBits[15]; // HLT_IsoMu24_v*
    if(triggerBits.size() > 16)
    { // Add here my trigger bits for TauPlusX
    }
    bool triggeredOn = singleETrigger || singleMuTrigger;
    if(singleETrigger)
    {
      if(isSingleMuPD)
        continue;
      chTags.push_back("singleE");
    }
    if(singleMuTrigger)
      chTags.push_back("singleMu");


    // Pileup Weight
    double weight       = 1.;
    double weight_plus  = 1.;
    double weight_minus = 1.;
    double puWeight     = 1.;
    if(isMC)
    {
      puWeight     = LumiWeights->weight(genEv.ngenITpu) * PUNorm[0];
      weight       = xsecWeight*puWeight;
      weight_plus  = PuShifters[utils::cmssw::PUUP  ]->Eval(genEv.ngenITpu) * (PUNorm[2]/PUNorm[0]);
      weight_minus = PuShifters[utils::cmssw::PUDOWN]->Eval(genEv.ngenITpu) * (PUNorm[1]/PUNorm[0]);
    }


    // Get Leading Lepton
    llvvLeptonCollection selLeptons;
    for(size_t i = 0; i < leptons.size(); ++i)
    {
      int lepId = leptons[i].id;

      if(abs(lepId) == 13)
        mon.fillHisto("isomu", "all", utils::cmssw::relIso(leptons[i], rho), weight);
      else if(abs(lepId) == 11)
        mon.fillHisto("isoele", "all", utils::cmssw::relIso(leptons[i], rho), weight);

      if(lepId == 13 && muCor)
      {
        TLorentzVector p4(leptons[i].px(), leptons[i].py(), leptons[i].pz(), leptons[i].energy());
        muCor->applyPtCorrection(p4, (lepId>0)?1:-1);
        if(isMC)
          muCor->applyPtSmearing(p4, (lepId>0)?1:-1, false);
        leptons[i].SetPxPyPzE(p4.Px(), p4.Py(), p4.Pz(), p4.Energy());
      }

      lepId = abs(lepId);

      // Lepton Kinematics
      bool passKin = true;
      double eta = (lepId == 11)?(leptons[i].electronInfoRef->sceta):(leptons[i].eta());
      if(lepId == 11) // If Electron
      {
        if(leptons[i].pt() < minElPt)
          passKin = false;
        if(abs(eta) > maxElEta)
          passKin = false;
        if(abs(eta) > ECALGap_MinEta && abs(eta) < ECALGap_MaxEta)  // Remove electrons that fall in ECAL Gap
          passKin = false;
      }
      else            // If Muon
      {
        if(leptons[i].pt() < minMuPt)
          passKin = false;
        if(abs(eta) > maxMuEta)
          passKin = false;
      }

      // Lepton ID
      bool passID = true;
      Int_t idbits = leptons[i].idbits;
      if(lepId == 11)
      {
        if(leptons[i].electronInfoRef->isConv)
          passID = false;

        bool isLoose = ((idbits >> 4) & 0x1);
        if(!isLoose)
          passID = false;
      }
      else
      {
        bool isLoose = ((idbits >> 8) & 0x1);
        bool isTight = ((idbits >> 10) & 0x1);
        if(!isLoose)
          passID = false;
      }

      // Lepton Isolation
      bool passIso = true;
      double relIso = utils::cmssw::relIso(leptons[i], rho);
      if((lepId == 11 && relIso > 0.15) || (lepId == 13 && relIso > 0.12))
        passIso = false;

      // Keep desired leptons
      if(passKin && passID && passIso)
        selLeptons.push_back(leptons[i]);

      // Fill lepton control plots
      mon.fillHisto("leptonCutFlow", "all", 0, weight);
      mon.fillHisto("leptonFractions", "all", 0, weight);
      if(passID)
      {
        mon.fillHisto("leptonCutFlow", "all", 1, weight);
        mon.fillHisto("leptonFractions", "all", 1, weight);
        if(passKin)
        {
          mon.fillHisto("leptonCutFlow", "all", 2, weight);
          if(passIso)
            mon.fillHisto("leptonCutFlow", "all", 3, weight);
        }
      }
      if(passKin)
        mon.fillHisto("leptonFractions", "all", 2, weight);
      if(passIso)
        mon.fillHisto("leptonFractions", "all", 3, weight);
    }
    if(selLeptons.size() != 0)
      std::sort(selLeptons.begin(), selLeptons.end(), sort_llvvObjectByPt);

    // Get Jets
    llvvJetExtCollection selJets, selJetsNoId, selBJets;
    int nTauJets = 0;
    for(size_t i = 0; i < jets.size(); ++i)
    {
      // Apply jet corrections
      double toRawSF = jets[i].torawsf;
      LorentzVector rawJet(jets[i]*toRawSF);
      jesCor->setJetEta(rawJet.eta());
      jesCor->setJetPt(rawJet.pt());
      jesCor->setJetA(jets[i].area);
      jesCor->setRho(rho);

      double newJECSF(jesCor->getCorrection());
      jets[i].SetPxPyPzE(rawJet.px(),rawJet.py(),rawJet.pz(),rawJet.energy());
      jets[i] *= newJECSF;
      jets[i].torawsf = 1./newJECSF;

      // Jet Pre-Selection
      bool passPreSel = true;
      if(jets[i].pt() < 15)
        passPreSel = false;
      if(abs(jets[i].eta()) > 4.7)
        passPreSel = false;

      // Cross clean with selected leptons and taus
      bool passIso = true;
      double minDRlj = 9999.9;
      double minDRlg = 9999.9;
      double minDRtj = 9999.9;
      for(size_t j = 0; j < selLeptons.size(); ++j)
        minDRlj = TMath::Min(minDRlj, deltaR(jets[i], selLeptons[j]));
      for(size_t j = 0; j < taus.size(); ++j)
      {
        if(taus[j].pt() < minTauPt || abs(taus[j].eta()) > maxTauEta)
          continue;
        minDRtj = TMath::Min(minDRtj, deltaR(jets[i], taus[j]));
      }
      if(minDRlj < 0.4 || minDRlg < 0.4 || minDRtj < 0.4)
        passIso = false;

      // Jet ID
      Int_t idbits = jets[i].idbits;
      bool passPFLoose = (idbits & 0x01);
      int puID = (idbits >> 3) & 0x0f;
      bool passLoosePuID = ((puID >> 2) & 0x01);
      int simplePuID = (idbits >> 7) & 0x0f;
      bool passLooseSimplePuID = ((simplePuID >> 2) & 0x01);
      std::string jetType = ((jets[i].genj.pt() > 0)?("truejetsid"):("pujetsid"));
      bool passID = passLoosePuID;

      // Jet Kinematics
      bool passKin = true;
      if(jets[i].pt() < minJetPt)
        passKin = false;
      if(abs(jets[i].eta()) > maxJetEta)
        passKin = false;
      bool isTauJet = (jets[i].pt() > minTauJetPt) && (abs(jets[i].eta()) < maxTauJetEta);

      // B-jets
      bool hasCSVV1L = jets[i].csv > 0.405;
      bool hasBtagCorr = hasCSVV1L;
      if(isMC)
      {
        // Get a "unique" seed
        double bseed_sin_phi = sin(jets[i].phi()*1000000);
        double bseed = abs(static_cast<int>(bseed_sin_phi*100000));

        // Get Jet Flavour
        int bflavid = jets[i].genflav;

        // Init
        BTagSFUtil btsfutil(bseed);

        TString flavKey("udsg");
        if(abs(bflavid) == 4)
          flavKey = "c";
        if(abs(bflavid) == 5)
          flavKey = "b";
        std::pair<TString,TString> btagKey("csvL", flavKey);
        if(btagEffCorr.find(btagKey) != btagEffCorr.end())
        {
          TGraphErrors* mceffGr = btagEffCorr[btagKey].first;
          TGraphErrors* sfGr    = btagEffCorr[btagKey].second;
          if(mceffGr && sfGr)
          {
            double jetpt = jets[i].pt();
            double eff = mceffGr->Eval(jetpt);
            double sf  = sfGr->Eval(jetpt);
            // Systematics: btagup, btagdown, ...
            //if(var == "btagup" || var == "btagdown" || var == "unbtagup" || var == "unbtagdown")
            //{
            //}
            btsfutil.modifyBTagsWithSF(hasBtagCorr, sf, eff);
          }
        }
      }

      // Compute scale and resolution uncertainties
      if(isMC)
      {
        std::vector<float> smearPt = utils::cmssw::smearJER(jets[i].pt(),jets[i].eta(),jets[i].genj.pt());
        jets[i].jer     = smearPt[0];
        jets[i].jerup   = smearPt[1];
        jets[i].jerdown = smearPt[2];
        smearPt = utils::cmssw::smearJES(jets[i].pt(),jets[i].eta(), totalJESUnc);
        jets[i].jesup   = smearPt[0];
        jets[i].jesdown = smearPt[1];
      }
      else
      {
        jets[i].jer     = jets[i].pt();
        jets[i].jerup   = jets[i].pt();
        jets[i].jerdown = jets[i].pt();
        jets[i].jesup   = jets[i].pt();
        jets[i].jesdown = jets[i].pt();
      }

      // Save selected jets/counters
      if(passPreSel && passIso)
        selJetsNoId.push_back(jets[i]);
      if(passPreSel && passIso && passPFLoose && passLoosePuID && isTauJet)
        ++nTauJets;
      if(passPreSel && passIso && passPFLoose && passID && passKin)
      {
        selJets.push_back(jets[i]);
        mon.fillHisto("jetcsv", chTags, jets[i].origcsv, weight);
      }
      if(passPreSel && passIso && passPFLoose && passID && passKin && hasBtagCorr)
        selBJets.push_back(jets[i]);

      // Fill Jet control histograms
      mon.fillHisto("jetCutFlow", chTags, 0, weight);
      mon.fillHisto("jetFractions", chTags, 0, weight);
      if(passPFLoose)
      {
        mon.fillHisto("jetCutFlow", chTags, 1, weight);
        mon.fillHisto("jetFractions", chTags, 1, weight);
        if(passPreSel)
        {
          mon.fillHisto("jetCutFlow", chTags, 2, weight);
          if(passID)
          {
            mon.fillHisto("jetCutFlow", chTags, 3, weight);
            if(passIso)
            {
              mon.fillHisto("jetCutFlow", chTags, 4, weight);
              if(passKin)
                mon.fillHisto("jetCutFlow", chTags, 5, weight);
            }
          }
        }
      }
      if(passPreSel)
        mon.fillHisto("jetFractions", chTags, 2, weight);
      if(passID)
        mon.fillHisto("jetFractions", chTags, 3, weight);
      if(passIso)
        mon.fillHisto("jetFractions", chTags, 4, weight);
      if(passKin)
        mon.fillHisto("jetFractions", chTags, 5, weight);
    }
    if(selJets.size() != 0)
      std::sort(selJets.begin(), selJets.end(), sort_llvvObjectByPt);
    if(selBJets.size() != 0)
      std::sort(selBJets.begin(), selBJets.end(), sort_llvvObjectByPt);

    // Get taus
    llvvTauCollection selTaus;
    for(size_t i = 0; i < taus.size(); ++i)
    {
      llvvTau& tau = taus[i];

      // Tau Kinematics
      bool passKin = true;
      if(tau.pt() < minTauPt)
        passKin = false;
      if(abs(tau.eta()) > maxTauEta)
        passKin = false;

      // Tau overlap with leptons
      bool passIso = true;
      for(size_t lep = 0; lep < selLeptons.size(); ++lep)
      {
        if(deltaR(tau, selLeptons[lep]) < 0.1)
        {
          passIso = false;
          break;
        }
      }

      bool passQual = true;
      if(abs(tau.dZ) > 0.5)
        passQual = false;
      if(tau.emfraction >= 2.0)
        passQual = false;

      // Tau ID
      bool passID = true;
      mon.fillHisto("tauID", chTags, 0, weight);
      //if(!tau.passId(llvvTAUID::againstElectronMediumMVA3))                   passID = false;
      if(!tau.passId(llvvTAUID::againstElectronMediumMVA5))                   passID = false;
      else mon.fillHisto("tauID", chTags, 1, weight);
      if(!tau.passId(llvvTAUID::againstMuonTight3))                           passID = false;
      else mon.fillHisto("tauID", chTags, 2, weight);
      if(!tau.passId(llvvTAUID::decayModeFinding))                            passID = false;
      else mon.fillHisto("tauID", chTags, 3, weight);
      if(!tau.passId(llvvTAUID::byMediumCombinedIsolationDeltaBetaCorr3Hits)) passID = false;
      else mon.fillHisto("tauID", chTags, 4, weight);
      if(passID) mon.fillHisto("tauID", chTags, 5, weight);

      // Keep selected taus
      if(passID && passKin && passIso && passQual && tau.isPF)
        selTaus.push_back(tau);

      // Fill control histograms
      mon.fillHisto("tauCutFlow", chTags, 0, weight);
      mon.fillHisto("tauFractions", chTags, 0, weight);
      if(tau.isPF)
      {
        mon.fillHisto("tauCutFlow", chTags, 1, weight);
        mon.fillHisto("tauFractions", chTags, 1, weight);
        if(passID)
        {
          mon.fillHisto("tauCutFlow", chTags, 2, weight);
          if(passQual)
          {
            mon.fillHisto("tauCutFlow", chTags, 3, weight);
            if(passKin)
            {
              mon.fillHisto("tauCutFlow", chTags, 4, weight);
              if(passIso)
                mon.fillHisto("tauCutFlow", chTags, 5, weight);
            }
          }
        }
      }
      if(passIso)
        mon.fillHisto("tauFractions", chTags, 5, weight);
      if(passKin)
        mon.fillHisto("tauFractions", chTags, 4, weight);
      if(passQual)
        mon.fillHisto("tauFractions", chTags, 3, weight);
      if(passID)
        mon.fillHisto("tauFractions", chTags, 2, weight);
    }
    if(selTaus.size() != 0)
      std::sort(selTaus.begin(), selTaus.end(), sort_llvvObjectByPt);


    mon.fillHisto("fractions", chTags, 0, weight);
    if(triggeredOn)
      mon.fillHisto("fractions", chTags, 1, weight);
    if(selLeptons.size() > 0)
      mon.fillHisto("fractions", chTags, 2, weight);
    if(selBJets.size() == 0)
      mon.fillHisto("fractions", chTags, 3, weight);
    if(selTaus.size() > 0)
      mon.fillHisto("fractions", chTags, 4, weight);

    mon.fillHisto("cutFlow", chTags, 0, weight);
    if(triggeredOn)
    {
      mon.fillHisto("cutFlow", chTags, 1, weight);
      if(selLeptons.size() > 0)
      {
        mon.fillHisto("cutFlow", chTags, 2, weight);
        mon.fillHisto("nbjets", chTags, selBJets.size(), weight);
        if(selBJets.size() == 0)
        {
          mon.fillHisto("cutFlow", chTags, 3, weight);
          if(selTaus.size() > 0)
          {
            mon.fillHisto("cutFlow", chTags, 4, weight);

            mon.fillHisto("nvtx", chTags, nvtx, weight);
            mon.fillHisto("nvtxraw", chTags, nvtx, weight/puWeight);
            mon.fillHisto("nup", "", genEv.nup, 1);

            mon.fillHisto("rho", chTags, rho, weight);
            mon.fillHisto("rho25", chTags, rho25, weight);

            mon.fillHisto("nlep", chTags, selLeptons.size(), weight);
            if(selLeptons.size() != 0)
            {
              mon.fillHisto("leadeta", chTags, selLeptons[0].eta(), weight);
              mon.fillHisto("leadpt", chTags, selLeptons[0].pt(), weight);
              mon.fillHisto("leadcharge", chTags, ((selLeptons[0].id > 0)?(-1):(1)), weight);
            }

            mon.fillHisto("njets", chTags, selJets.size(), weight);
            if(selJets.size() != 0)
            {
              mon.fillHisto("jetleadpt", chTags, selJets[0].pt(), weight);
              mon.fillHisto("jetleadeta", chTags, selJets[0].eta(), weight);
            }

            mon.fillHisto("ntaus", chTags, selTaus.size(), weight);
            if(selTaus.size() != 0)
            {
              mon.fillHisto("tauleadpt", chTags, selTaus[0].pt(), weight);
              mon.fillHisto("tauleadeta", chTags, selTaus[0].eta(), weight);
              mon.fillHisto("tauleadcharge", chTags, ((selTaus[0].id > 0)?(-1):(1)), weight);
              mon.fillHisto("tauleademfrac", chTags, selTaus[0].emfraction, weight);
              mon.fillHisto("taudz", chTags, selTaus[0].dZ, weight);
            }
          }
        }
      }
    }

//    break;
  }

  // Output temporary buffer and restore cout and cerr behaviour
  std::cout.rdbuf(coutbuf);
  std::cerr.rdbuf(cerrbuf);
  std::cout << std::endl;
  std::cout << buffer.str();

  std::cout << "totalEntries: " << totalEntries << "; vs nInitEvent: " << nInitEvent << ";" << std::endl;


  /***************************************************************************/
  /*                        Saving Histograms to File                        */
  /***************************************************************************/
  outUrl += "/";
  outUrl += outFileUrl + ".root";
  std::cout << "Saving results in " << outUrl << std::endl;
  TFile* outfile = new TFile(outUrl.c_str(), "RECREATE");
  mon.Write();
  outfile->Close();
  delete outfile;

  if(saveSummaryTree)
  {
    outUrl.replace(outUrl.find(".root", 0), 5, "_summary.root");
    std::cout << "Saving summary results in " << outUrl << std::endl;
    outfile = new TFile(outUrl.c_str(), "RECREATE");
    // Write Tuple/Tree
    //summaryTuple->SetDirectory(ofile);
    //summaryTuple->Write();
    outfile->Close();
    delete outfile;
  }

  return 0;
}
