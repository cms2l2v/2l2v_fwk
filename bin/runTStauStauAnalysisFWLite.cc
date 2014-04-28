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
  float minJetPtToApply(30); // Min Jet Pt to accept jet?
  bool examineThisEvent(false); // ?
  double minTauPt = 20;
  double maxTauEta = 2.3;

  // Lepton Efficiencies
  LeptonEfficiencySF lepEff;


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



  /***************************************************************************/
  /*                         Initializing Histograms                         */
  /***************************************************************************/
  SmartSelectionMonitor mon;
  TH1F *h = (TH1F*)mon.addHistogram(new TH1F("cutFlow", ";;Events", 20, 0, 20));
  h->GetXaxis()->SetBinLabel(1, "All");
  h->GetXaxis()->SetBinLabel(2, "HLT");
  h->GetXaxis()->SetBinLabel(3, "> 1l");
  h->GetXaxis()->SetBinLabel(4, "> 1#tau");
  // ...
  // TH2D* hist = (TH2D*)mon.addHistogram(...);

  mon.addHistogram(new TH1F("nup",     ";NUP;Events", 10, 0, 10));
  mon.addHistogram(new TH1F("nupfilt", ";NUP;Events", 10, 0, 10));

  // PU
  mon.addHistogram(new TH1F("nvtx",    ";Vertices;Events",       50, -0.5, 49.5));
  mon.addHistogram(new TH1F("nvtxraw", ";Vertices;Events",       50, -0.5, 49.5));
  mon.addHistogram(new TH1F("rho",     ";#rho;Events",           50,  0,   25));
  mon.addHistogram(new TH1F("rho25",   ";#rho(#eta<2.5);Events", 50,  0,   25));

  // Leptons
  mon.addHistogram(new TH1F("nlep",       ";nlep;Events",       10,  0,   10));
  mon.addHistogram(new TH1F("leadpt",     ";p_{T}^{l};Events", 100,  0,  100));
  mon.addHistogram(new TH1F("leadeta",    ";#eta^{l};Events",   50, -2.6,  2.6));
  mon.addHistogram(new TH1F("leadcharge", ";q^{l};Events",   5, -2,    2));

  // Lepton Isolation
  mon.addHistogram(new TH1F("isomu",  "RelIso(#mu)", 100, -0.5, 9.5));
  mon.addHistogram(new TH1F("isoele", "RelIso(ele)", 100, -0.5, 9.5));

  // Taus
  mon.addHistogram(new TH1F("ntaus",         ";ntaus;Events",         10, -0.5,  9.5));
  mon.addHistogram(new TH1F("tauleadpt",     ";p_{T}^{#tau};Events", 100,  0,  100));
  mon.addHistogram(new TH1F("tauleadeta",    ";#eta^{#tau};Events",   50, -2.6,  2.6));
  mon.addHistogram(new TH1F("tauleadcharge", ";q^{#tau};Events",       5, -2,    2));
  mon.addHistogram(new TH1F("taudz",         ";dz^{#tau};Events",     50,  0,   10));
  mon.addHistogram(new TH1F("tauvz",         ";vz^{#tau};Events",     50,  0,   10));
  mon.addHistogram(new TH1F("tauleademfrac", ";emf^{#tau};Events",    50,  0,    5));

  // Jets
  mon.addHistogram(new TH1F("njets", ";njets;Events", 6, 0, 6));



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


    /****         Filter events acording to HLT Path         ****/
    bool singleETrigger  = triggerBits[13]; // HLT_Ele27_WP80_v*
    bool singleMuTrigger = triggerBits[15]; // HLT__IsoMu24_v*
    if(singleETrigger)
      chTags.push_back("singleE");
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

    // Fill the all events bin in the cutflow:
    mon.fillHisto("cutFlow", chTags, 0, weight);

    if(!(singleETrigger || singleMuTrigger))
      continue;
    mon.fillHisto("cutFlow", chTags, 1, weight);

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
      if(leptons[i].pt() < ((lepId == 11)?(35):(30)))  // Remove low Pt leptons
        passKin = false;
      if(abs(eta) > ((lepId == 11)?(2.5):(2.4))) // Only keep leptons inside detector acceptance (different for el and mu)
        passKin = false;
      if(lepId == 11 && (eta > 1.4442 && eta < 1.5660)) // Remove electrons that fall in ECAL "hole"
        passKin = false;

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
    }
    if(selLeptons.size() == 0)
      continue;
    std::sort(selLeptons.begin(), selLeptons.end(), sort_llvvObjectByPt);
    mon.fillHisto("cutFlow", chTags, 2, weight);

    // Get Jets
    llvvJetExtCollection selJets, selJetsNoId, selBJets;
    int nJets = 0;
    int nbJets = 0;
    int nTrueJets = 0;
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

      // Jet Kinematics
      if(jets[i].pt() < 15)
        continue;
      if(abs(jets[i].eta()) > 4.7)
        continue;

      // Cross clean with selected leptons and taus
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
        continue;
    }

    // Get taus
    llvvTauCollection selTaus;
    for(size_t i = 0; i < taus.size(); ++i)
    {
      llvvTau& tau = taus[i];

      // Tau Kinematics
      if(tau.pt() < minTauPt)
        continue;
      if(abs(tau.eta()) > maxTauEta)
        continue;

      // Tau overlap with leptons
      bool overlapWithLepton = false;
      for(size_t lep = 0; lep < selLeptons.size(); ++lep)
      {
        if(deltaR(tau, selLeptons[lep]) < 0.1)
        {
          overlapWithLepton = true;
          break;
        }
      }
      if(overlapWithLepton)
        continue;

      if(!tau.isPF)  // Only keep PF taus
        continue;
      if(abs(tau.dZ) > 0.5)
        continue;
      if(tau.emfraction >= 2.0)
        continue;

      // Tau ID
      if(!tau.passId(llvvTAUID::againstElectronMediumMVA5))                   continue;
      if(!tau.passId(llvvTAUID::againstMuonTight2))                           continue;
      if(!tau.passId(llvvTAUID::decayModeFinding))                            continue;
      if(!tau.passId(llvvTAUID::byMediumCombinedIsolationDeltaBetaCorr3Hits)) continue;

      selTaus.push_back(tau);
    }
    if(selTaus.size() == 0)
      continue;
    std::sort(selTaus.begin(), selTaus.end(), sort_llvvObjectByPt);
    mon.fillHisto("cutFlow", chTags, 3, weight);


    mon.fillHisto("nvtx", chTags, nvtx, weight);
    mon.fillHisto("nvtxraw", chTags, nvtx, weight/puWeight);
    mon.fillHisto("nup", "", genEv.nup, 1);

    mon.fillHisto("rho", chTags, rho, weight);
    mon.fillHisto("rho25", chTags, rho25, weight);

    mon.fillHisto("nlep", chTags, selLeptons.size(), weight);
    mon.fillHisto("leadeta", chTags, selLeptons[0].eta(), weight);
    mon.fillHisto("leadpt", chTags, selLeptons[0].pt(), weight);
    mon.fillHisto("leadcharge", chTags, ((selLeptons[0].id > 0)?(-1):(1)), weight);

    mon.fillHisto("ntaus", chTags, selTaus.size(), weight);
    mon.fillHisto("tauleadpt", chTags, selTaus[0].pt(), weight);
    mon.fillHisto("tauleadeta", chTags, selTaus[0].eta(), weight);
    mon.fillHisto("tauleadcharge", chTags, ((selTaus[0].id > 0)?(-1):(1)), weight);

//    break;
  }

  // Output temporary buffer and restore cout and cerr behaviour
  std::cout.rdbuf(coutbuf);
  std::cerr.rdbuf(cerrbuf);
  std::cout << std::endl;
  std::cout << buffer.str();


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
