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

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h" //for SVfit


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
#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TEventList.h"
#include "TNtuple.h"
#include "TGraphErrors.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRotation.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/shared_ptr.hpp>
#include <Math/VectorUtil.h>
#include <bitset>
#include <cctype>
#include <cmath>

// Include MT2 library:
// http://particle.physics.ucdavis.edu/hefti/projects/doku.php?id=wimpmass    ** Code from here
// http://www.hep.phy.cam.ac.uk/~lester/mt2/    ** Other libraries
#include "UserCode/llvv_fwk/interface/mt2_bisect.h"


#ifndef DEBUG_EVENT
//#define DEBUG_EVENT true
#endif

#define NAN_WARN(X) if(std::isnan(X)) std::cout << "  Warning: " << #X << " is nan" << std::endl;



enum ID_Type {LooseID, MediumID, TightID};
enum TAU_E_ID {antiELoose, antiEMedium, antiETight, antiEMva, antiEMva3Loose, antiEMva3Medium, antiEMva3Tight, antiEMva3VTight, antiEMva5Medium};

double stauCrossSec(double stauM, double neutM);
bool electronMVAID(double mva, llvvLepton& lepton, ID_Type id);
int old_main(int argc, char* argv[]);
double tauSF(llvvTau& tau, llvvGenParticleCollection& genPartColl, TAU_E_ID eId);
double leptonIdAndIsoScaleFactor(llvvLepton& lepton);
double leptonTauTriggerScaleFactor(llvvLepton& lepton, llvvTau& tau);
double efficiency(double m, double m0, double sigma, double alpha, double n, double norm);

/*****************************************************************************/
/* Return Codes:                                                             */
/*   0 - Everything OK                                                       */
/*   1 - Missing parameters_cfg.py configuration file                        */
/*****************************************************************************/
int main(int argc, char* argv[])
{
  if(argc < 2)
    std::cout << "Usage: " << argv[0] << " parameters_cfg.py" << std::endl, exit(1);

  size_t limit = 0;

  #if defined(DEBUG_EVENT)
  bool debugEvent = false;
  int skipEvents = 0;
  #endif

  int fileIndex = 1;
  if(argc > 2)
  {
    std::stringstream parser;

    for(int i = 1; i < argc; ++i)
    {
      if(argv[i][0] != '-')
      {
        fileIndex = i;
        break;
      }

      std::string arg = argv[i];
      if(arg.find("--limit") != std::string::npos)
      {
        char first = argv[i+1][0];
        if(!isdigit(first))
          continue;

        parser << argv[i+1];
        parser >> limit;

        ++i;
        continue;
      }

      #if defined(DEBUG_EVENT)
      if(arg.find("--debugEvent") != std::string::npos)
      {
        debugEvent = true;
        continue;
      }
      if(arg.find("--skipEvents") != std::string::npos)
      {
        char first = argv[i+1][0];
        if(!isdigit(first))
          continue;

        parser << argv[i+1];
        parser >> skipEvents;

        ++i;
        continue;
      }
      #endif
    }
  }

  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

  // Read parameters from the configuration file
  const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[fileIndex])->getParameter<edm::ParameterSet>("runProcess");
  // If we are supposed to run the old analysis, run the old analysis
  bool doOldAnalysis = false;
  if(runProcess.exists("doOldAnalysis"))
    doOldAnalysis = runProcess.getParameter<bool>("doOldAnalysis");
  if(doOldAnalysis)
  {
    return old_main(argc, argv); // Old analysis has been deprecated, only kept here for historical reasons
  }

  bool isMC = runProcess.getParameter<bool>("isMC");
  double xsec = runProcess.getParameter<double>("xsec");
  std::vector<std::string> urls = runProcess.getParameter<std::vector<std::string> >("input");
  std::string baseDir = runProcess.getParameter<std::string>("dirName");
  std::string outdir = runProcess.getParameter<std::string>("outdir");
  std::string jecDir = runProcess.getParameter<std::string>("jecDir");
  bool runSystematics = runProcess.getParameter<bool>("runSystematics");
  bool saveSummaryTree = runProcess.getParameter<bool>("saveSummaryTree");
  bool exclusiveRun = runProcess.getParameter<bool>("exclusiveRun");
  double stauMtoPlot = 120;
  double neutralinoMtoPlot = 20;  // TStauStau mass point to plot by default
  if(runProcess.exists("stauMtoPlot"))
    stauMtoPlot = runProcess.getParameter<double>("stauMtoPlot");
  if(runProcess.exists("neutralinoMtoPlot"))
    neutralinoMtoPlot = runProcess.getParameter<double>("neutralinoMtoPlot");
  bool doSVfit = false;
  if(runProcess.exists("doSVfit"))
    doSVfit = runProcess.getParameter<bool>("doSVfit");
  bool applyScaleFactors = false;
  if(runProcess.exists("applyScaleFactors"))
    applyScaleFactors = runProcess.getParameter<bool>("applyScaleFactors");
  bool debug = false;
  if(runProcess.exists("debug"))
    debug = runProcess.getParameter<bool>("debug");
  bool doTightTauID = false;
  if(runProcess.exists("doTightTauID"))
    debug = runProcess.getParameter<bool>("doTightTauID");

  if(debug)
    std::cout << "Finished loading config file" << std::endl;

  // Hardcoded Values
  double sqrtS          =  8;      // Center of mass energy
  double minElPt        = 24;      // Selected electron pT and eta
  double maxElEta       =  2.1;
  double ECALGap_MinEta =  1.4442; // ECAL gap parameters
  double ECALGap_MaxEta =  1.5660;
  double minMuPt        = 20;      // Selected muon pT and eta
  double maxMuEta       =  2.1;
  double minTauPt       = 20;      // Selected tau pT and eta (I was using 25)
  double maxTauEta      =  2.3;
  double maxJetEta      =  4.7;    // Selected jet eta

  // Setting up -------------------------------------------------------------------------
  if(debug)
    std::cout << "Setting up" << std::endl;
  gSystem->Exec(("mkdir -p " + outdir).c_str());
  std::string url = urls[0];
  std::string outFileUrl(gSystem->BaseName(url.c_str()));
  while(outFileUrl.find(".root", 0) != std::string::npos)
    outFileUrl.replace(outFileUrl.find(".root", 0), 5, "");
  std::string outUrl = outdir;
  outUrl += "/";
  outUrl += outFileUrl + ".root";

  TString turl(url);
  bool isV0JetsMC(isMC && (turl.Contains("DYJetsToLL_50toInf") || turl.Contains("WJets")));
  bool isStauStau(isMC && turl.Contains("TStauStau"));

  TTree* summaryTree = NULL;
  TFile* summaryOutFile = NULL;
  if(saveSummaryTree)
  {
    TDirectory* cwd = gDirectory;

    std::string summaryOutUrl = outUrl;
    summaryOutUrl.replace(summaryOutUrl.find(".root", 0), 5, "_summary.root");
    std::cout << "Saving summary results in " << summaryOutUrl << std::endl;
    summaryOutFile = new TFile(summaryOutUrl.c_str(), "RECREATE");

    summaryTree = new TTree("Events", "Events");

    summaryTree->SetDirectory(summaryOutFile);

    cwd->cd();
  }



  /***************************************************************************/
  /*                         Initializing Histograms                         */
  /***************************************************************************/
  if(debug)
    std::cout << "Initializing histograms" << std::endl;
  SmartSelectionMonitor mon;
  TH1D *eventflow = (TH1D*)mon.addHistogram(new TH1D("eventflow", ";;Events", 8, 0, 8));
  eventflow->GetXaxis()->SetBinLabel(1, "HLT");
  eventflow->GetXaxis()->SetBinLabel(2, "MET > 30");
  eventflow->GetXaxis()->SetBinLabel(3, "> 1l");
  eventflow->GetXaxis()->SetBinLabel(4, "> 1#tau");
  eventflow->GetXaxis()->SetBinLabel(5, "B-veto");
  eventflow->GetXaxis()->SetBinLabel(6, "OS");
  eventflow->GetXaxis()->SetBinLabel(7, "lep veto");
  eventflow->GetXaxis()->SetBinLabel(8, "SVfit");

  mon.addHistogram(new TH1D("nup", ";NUP;Events", 10, 0, 10));

  // Pile Up
//  mon.addHistogram(new TH1D("nvtxAll", ";Vertices;Events", 50, -0.5, 49.5));
  mon.addHistogram(new TH1D("nvtx", ";Vertices;Events", 50, -0.5, 49.5));
  mon.addHistogram(new TH1D("nvtxraw", ";Vertices;Events", 50, -0.5, 49.5));
  mon.addHistogram(new TH1D("rho", ";#rho;Events", 25, 0, 25));
  mon.addHistogram(new TH1D("rho25", ";#rho(#eta<2.5);Events", 25, 0, 25));


  // Leptons
  mon.addHistogram(new TH1D("nlep", ";nlep;Events", 10, 0, 10));
  mon.addHistogram(new TH1D("ptSelectedLep", ";p_{T}^{l};Events", 50, 0, 100));
  mon.addHistogram(new TH1D("etaSelectedLep", ";#eta^{l};Events", 25, -2.6, 2.6));
  mon.addHistogram(new TH1D("chargeSelectedLep", ";q^{l};Events", 5, -2, 2));
  TH1D *leptonCutFlow = (TH1D*)mon.addHistogram(new TH1D("leptonCutFlow", ";;Leptons", 4, 0, 4));
  leptonCutFlow->GetXaxis()->SetBinLabel(1, "All");
  leptonCutFlow->GetXaxis()->SetBinLabel(2, "ID");
  leptonCutFlow->GetXaxis()->SetBinLabel(3, "Kin");
  leptonCutFlow->GetXaxis()->SetBinLabel(4, "Iso");

  // Lepton Isolation
  mon.addHistogram(new TH1D("isomu", "RelIso(#mu);;Leptons", 100, -0.5, 9.5));
  mon.addHistogram(new TH1D("isoele", "RelIso(ele);;Leptons", 100, -0.5, 9.5));

  // Taus
  mon.addHistogram(new TH1D("ntaus", ";ntaus;Events", 10, 0, 10));
  mon.addHistogram(new TH1D("ptSelectedTau", ";p_{T}^{#tau};Events", 50, 0, 100));
  mon.addHistogram(new TH1D("ptSelectedTauExtended", ";p_{T}^{#tau};Events", 50, 0, 250));
  mon.addHistogram(new TH1D("etaSelectedTau", ";#eta^{#tau};Events", 25, -2.6, 2.6));
  mon.addHistogram(new TH1D("chargeSelectedTau", ";q^{#tau};Events", 5, -2, 2));
  mon.addHistogram(new TH1D("dzSelectedTau", ";dz^{#tau};Events", 25, 0, 2));
  mon.addHistogram(new TH1D("emfracSelectedTau", ";emf^{#tau};Events", 25, 0, 5));
  TH1D *tauCutFlow = (TH1D*)mon.addHistogram(new TH1D("tauCutFlow", ";;#tau", 6, 0, 6));
  tauCutFlow->GetXaxis()->SetBinLabel(1, "All");
  tauCutFlow->GetXaxis()->SetBinLabel(2, "PF");
  tauCutFlow->GetXaxis()->SetBinLabel(3, "ID");
  tauCutFlow->GetXaxis()->SetBinLabel(4, "Quality");
  tauCutFlow->GetXaxis()->SetBinLabel(5, "Kin");
  tauCutFlow->GetXaxis()->SetBinLabel(6, "Iso");
  TH1D *tauID = (TH1D*)mon.addHistogram(new TH1D("tauID", ";;#tau", 5, 0, 5));
  tauID->GetXaxis()->SetBinLabel(1, "All");
  tauID->GetXaxis()->SetBinLabel(2, "Medium comb iso");
  tauID->GetXaxis()->SetBinLabel(3, "Decay mode");
  tauID->GetXaxis()->SetBinLabel(4, "Not e");
  tauID->GetXaxis()->SetBinLabel(5, "Not #mu");

  // Jets
  mon.addHistogram(new TH1D("njets", ";njets;Events", 6, 0, 6));
  mon.addHistogram(new TH1D("nbjets", ";njets;Events", 6, 0, 6));
  mon.addHistogram(new TH1D("jetleadpt", ";p_{T}^{jet};Events", 25, 0, 500));
  mon.addHistogram(new TH1D("jetleadeta", ";#eta^{jet};Events", 50, -5, 5));
  mon.addHistogram(new TH1D("jetcsv", ";csv;jets", 25, 0, 1));
  TH1D *jetCutFlow = (TH1D*)mon.addHistogram(new TH1D("jetCutFlow", ";;jets", 4, 0, 4));
  jetCutFlow->GetXaxis()->SetBinLabel(1, "All");
  jetCutFlow->GetXaxis()->SetBinLabel(2, "PF Loose");
  jetCutFlow->GetXaxis()->SetBinLabel(3, "ID");
  jetCutFlow->GetXaxis()->SetBinLabel(4, "Kin");

  // MET
  mon.addHistogram(new TH1D("MET", ";MET [GeV];Events", 25, 0, 200));

  // MT
  mon.addHistogram(new TH1D("MT", ";MT [GeV];Events", 25, 0, 200));

  // Deconstructed MT: https://indico.cern.ch/event/344807/
  mon.addHistogram(new TH1D("Q80", ";Q_{80};Events", 30, -2, 1));
  mon.addHistogram(new TH1D("Q100", ";Q_{100};Events", 30, -2, 1));
  mon.addHistogram(new TH1D("cosPhi", ";cos#Phi;Events", 30, -1, 1));

  // MT2
  mon.addHistogram(new TH1D("MT2", ";M_{T2} [GeV];Events", 25, 0, 500));

  // SVFit Mass
  if(doSVfit)
    mon.addHistogram(new TH1D("SVFitMass", ";M_{SVFit};Events", 50, 0, 500));

  // Invariant Mass
  mon.addHistogram(new TH1D("InvMass", ";M_{l-#tau};Events", 50, 0, 500));

  // Angles
  mon.addHistogram(new TH1D("deltaAlphaLepTau", ";#Delta#alpha_{l-#tau}(Lab);Events", 30, 0, TMath::Pi()));
  mon.addHistogram(new TH1D("deltaPhiLepTauMET", ";#Delta#phi_{l#tau-MET}(Lab);Events", 30, 0, TMath::Pi()));
  mon.addHistogram(new TH1D("deltaPhiLepTau", ";#Delta#phi_{l-#tau}(Lab);Events", 30, 0, TMath::Pi()));
  mon.addHistogram(new TH1D("cosThetaTau", ";cos#theta_{#tau}(Lab);Events", 30, -1, 1));
  mon.addHistogram(new TH1D("cosThetaLep", ";cos#theta_{l}(Lab);Events", 30, -1, 1));
  mon.addHistogram(new TH1D("deltaPhiLepMETCS", ";#Delta#phi_{l-MET}(CS);Events", 30, 0, TMath::Pi()));
  mon.addHistogram(new TH1D("cosThetaCS", ";cos#theta(CS);Events", 30, -1, 1));
  mon.addHistogram(new TH1D("minDeltaPhiMETJetPt40", ";min(#Delta#phi_{MET-Jet40});Events", 20, 0, TMath::Pi()));

  // 2D variables
  mon.addHistogram(new TH2D("metVSPTl", ";p_{T}(l);MET", 50, 0, 100, 25, 0, 200));
  mon.addHistogram(new TH2D("metVSPTtau", ";p_{T}(#tau);MET", 50, 0, 100, 25, 0, 200));
  //  Deconstructed MT 2D Plots:
  mon.addHistogram(new TH2D("Q80VsCosPhi", ";cos#Phi;Q_{80}", 20, -1, 1, 20, -2, 1));
  mon.addHistogram(new TH2D("Q100VsCosPhi", ";cos#Phi;Q_{100}", 20, -1, 1, 20, -2, 1));



  /***************************************************************************/
  /*                          Prepare for Event Loop                         */
  /***************************************************************************/
  if(debug)
    std::cout << "Preparing for event loop" << std::endl;
  fwlite::ChainEvent ev(urls);
  const size_t totalEntries = ev.size();

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
    while(mcPileupDistribution.size() < dataPileupDistribution.size()) mcPileupDistribution.push_back(0.0);
    while(mcPileupDistribution.size() > dataPileupDistribution.size()) dataPileupDistribution.push_back(0.0);

    LumiWeights= new edm::LumiReWeighting(mcPileupDistribution, dataPileupDistribution);
    PuShifters=utils::cmssw::getPUshifters(dataPileupDistribution, 0.05);
    utils::getPileupNormalization(mcPileupDistribution, PUNorm, LumiWeights, PuShifters);
  }

  gROOT->cd(); //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE

  DuplicatesChecker duplicatesChecker;
//  int nDuplicates(0);
  int step = int(totalEntries/50);

  // Redirect stdout and stderr to a temporary buffer, then output buffer after event loop
  std::ostream myCout(std::cout.rdbuf());
  std::stringstream buffer;
  std::streambuf *coutbuf = std::cout.rdbuf();
  std::streambuf *cerrbuf = std::cerr.rdbuf();
  std::cout.rdbuf(buffer.rdbuf());
  std::cerr.rdbuf(buffer.rdbuf());

  // Variables used in loop
  if(debug)
    myCout << "  Declaring all variables used in loop" << std::endl;
  int nvtx = 0;
  std::vector<TString> chTags;
  std::vector<bool> triggerBits;
  bool triggeredOn = false;
  llvvMet met;
  double rho = 0;
  double rho25 = 0;
  double crossSection = xsec;
  double weight = 1.;
  double weight_plus = 1.;
  double weight_minus = 1.;
  double puWeight = 1.;
  double triggerSF = 1.;
  double leptonIdIsoSF = 1.;
  double tauSF = 1.;
  llvvLeptonCollection selLeptons;
  llvvJetExtCollection selJets;
  llvvJetCollection selJetsOut;
  llvvJetCollection selBJets;
  llvvTauCollection selTaus;
  int nJets = 0;
  int nBJets = 0;
  bool selected = false;
  int tauIndex = -1, leptonIndex = -1;
  bool isOS = false;
  bool isMultilepton = false;
  bool isSVfit = true;
  double mass = -1;
  double invMass = -1;
  double mt = -1;
  double Q80 = 2;
  double Q100 = 2;
  double cosPhi = -10;
  double mt2 = -1;
  double stauMass = 0;
  double neutralinoMass = 0;
  double deltaAlphaLepTau = 0;
  double deltaPhiLepTauMET = 0;
  double deltaPhiLepTau = 0;
  double cosThetaTau = 0;
  double cosThetaLep = 0;
  double deltaPhiLepMETCS = 0;
  double cosThetaCS = 0;
  double minDeltaPhiMETJetPt40 = 0;
  double tauLeadPt = 0;
  double lepLeadPt = 0;
  double maxPtSum = 0;
//  int nTauJets = 0;

  // Prepare summary tree
  if(saveSummaryTree)
  {
    if(debug)
      std::cout << "  Defining all branches in output root file" << std::endl;

    TDirectory* cwd = gDirectory;
    summaryOutFile->cd();

    // Dataset specific variables
    summaryTree->Branch("isMC", &isMC);
    summaryTree->Branch("xSecWeight", &xsecWeight);

    // Event specific variables
    summaryTree->Branch("selected", &selected);
    summaryTree->Branch("chTags", &chTags);
    summaryTree->Branch("nvtx", &nvtx);
    summaryTree->Branch("triggerBits", &triggerBits);
    summaryTree->Branch("triggeredOn", &triggeredOn);
    summaryTree->Branch("rho", &rho);
    summaryTree->Branch("rho25", &rho25);
    summaryTree->Branch("met", &met);
    summaryTree->Branch("crossSection", &crossSection);
    summaryTree->Branch("weight", &weight);
    summaryTree->Branch("weight_plus", &weight_plus);
    summaryTree->Branch("weight_minus", &weight_minus);
    summaryTree->Branch("puWeight", &puWeight);
    summaryTree->Branch("triggerSF", &triggerSF);
    summaryTree->Branch("leptonIdIsoSF", &leptonIdIsoSF);
    summaryTree->Branch("tauSF", &tauSF);
    summaryTree->Branch("selLeptons", &selLeptons);
    summaryTree->Branch("selTaus", &selTaus);
//    summaryTree->Branch("selJets", &selJetsOut);
//    summaryTree->Branch("selBJets", &selBJets);
    summaryTree->Branch("nJets", &nJets);
    summaryTree->Branch("nBJets", &nBJets);
    summaryTree->Branch("isOS", &isOS);
    summaryTree->Branch("isMultilepton", &isMultilepton);
    summaryTree->Branch("isSVfit", &isSVfit);
    summaryTree->Branch("tauIndex", &tauIndex);
    summaryTree->Branch("leptonIndex", &leptonIndex);
    if(doSVfit)
    {
      summaryTree->Branch("SVFitMass", &mass);
    }
    summaryTree->Branch("InvariantMass", &invMass);
    summaryTree->Branch("MT", &mt);
    summaryTree->Branch("Q80", &Q80);
    summaryTree->Branch("Q100", &Q100);
    summaryTree->Branch("cosPhi", &cosPhi);
    summaryTree->Branch("MT2", &mt2);
    summaryTree->Branch("stauMass", &stauMass);
    summaryTree->Branch("neutralinoMass", &neutralinoMass);
    summaryTree->Branch("deltaAlphaLepTau", &deltaAlphaLepTau);
    summaryTree->Branch("deltaPhiLepTauMET", &deltaPhiLepTauMET);
    summaryTree->Branch("deltaPhiLepTau", &deltaPhiLepTau);
    summaryTree->Branch("cosThetaTau", &cosThetaTau);
    summaryTree->Branch("cosThetaLep", &cosThetaLep);
    summaryTree->Branch("deltaPhiLepMETCS", &deltaPhiLepMETCS);
    summaryTree->Branch("cosThetaCS", &cosThetaCS);
    summaryTree->Branch("minDeltaPhiMETJetPt40", &minDeltaPhiMETJetPt40);
    summaryTree->Branch("tauLeadPt", &tauLeadPt);
    summaryTree->Branch("lepLeadPt", &lepLeadPt);

    cwd->cd();
  }



  /***************************************************************************/
  /*                                Event Loop                               */
  /***************************************************************************/
  myCout << "       Progress Bar:0%      20%       40%       60%       80%      100%" << std::endl;
  myCout << "Scanning the ntuple:";

  // Loop on events
  for(size_t iev = 0; iev < totalEntries; ++iev)
  {
    #if defined(DEBUG_EVENT)
    if(iev < skipEvents)
      continue;
    if(debugEvent)
      myCout << "## Event " << iev << std::endl;
    #endif
    if(iev%step == 0)
      myCout << "_" << std::flush;

    // Init variables
    deltaAlphaLepTau = 0;
    deltaPhiLepTauMET = 0;
    deltaPhiLepTau = 0;
    cosThetaTau = 0;
    cosThetaLep = 0;
    deltaPhiLepMETCS = 0;
    cosThetaCS = 0;
    minDeltaPhiMETJetPt40 = 0;
    nvtx = 0;
    selected = false;
    weight = 1.;
    weight_plus = 1.;
    weight_minus = 1.;
    puWeight = 1.;
    triggerSF = 1.;
    leptonIdIsoSF = 1.;
    tauSF = 1.;
    chTags.clear();
    selLeptons.clear();
    nJets = 0;
    nBJets = 0;
    selJets.clear();
    selJetsOut.clear();
    selBJets.clear();
    selTaus.clear();
    tauIndex = -1, leptonIndex = -1;
    isOS = false;
    isMultilepton = false;
    isSVfit = false;
    mass = -1;
    invMass = -1;
    mt = -1;
    Q80 = 2;
    Q100 = 2;
    cosPhi = -10;
    mt2 = -1;
    stauMass = -1;
    neutralinoMass = -1;
    tauLeadPt = 0;
    lepLeadPt = 0;
    maxPtSum = 0;
//    nTauJets = 0;

    // Prepare tags to fill the histograms
    chTags.push_back("all");

    // Load the event content from tree
    ev.to(int(iev));


    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Loading collections" << std::endl;
    #endif

    /**** Get information/collections from the event ****/
    // Number of vertexes
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
    /**** Remove double counting if running on exclusive samples ****/
    if(exclusiveRun && isV0JetsMC)
    {
      if(genEv.nup > 5) // Drop V+{1,2,3,4}Jets from VJets samples to avoid double counting (but keep V+0Jets) [V = W,Z]
        continue;
    }

    /**** Get LHE comments with mass info ****/
    if(isStauStau)
    {
      fwlite::Handle<LHEEventProduct> LHEHandle;
      LHEHandle.getByLabel(ev, "source");
      if(!LHEHandle.isValid())
      {
        std::cout << "LHEEventProduct Object not Found" << std::endl;
        continue;
      }
      if(LHEHandle->comments_size() == 0)
        continue;

      for(auto comment = LHEHandle->comments_begin(); comment != LHEHandle->comments_end(); ++comment)
      {
        auto modelPos = comment->find("# model TStauStau_");
        if(modelPos != std::string::npos)
        {
          std::stringstream tmp;
          auto numPos = comment->find_first_of("1234567890", modelPos);

          tmp << comment->substr(numPos, comment->find("_", numPos)-numPos);
          tmp >> stauMass;
          tmp.clear();

          numPos = comment->find("_", numPos);
          numPos = comment->find_first_of("1234567890", numPos);
          tmp << comment->substr(numPos, comment->find("\n", numPos)-numPos);
          tmp >> neutralinoMass;

          break;
        }
      }
    }

    // Trigger Bits
    fwlite::Handle<std::vector<bool> > triggerBitsHandle;
    triggerBitsHandle.getByLabel(ev, "llvvObjectProducersUsed", "triggerBits");
    if(!triggerBitsHandle.isValid())
    {
      std::cout << "triggerBits Object NotFound" << std::endl;
      continue;
    }
    triggerBits = *triggerBitsHandle;

    //bool singleETrigger = triggerBits[13]; // HLT_Ele27_WP80_v*
    //bool singleMuTrigger = triggerBits[15]; // HLT_IsoMu24_v*
    bool TauPlusE2012A = triggerBits[18]; // HLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v*
    bool TauPlusMu2012A = triggerBits[22]; // HLT_IsoMu18_eta2p1_LooseIsoPFTau20_v*
    //bool TauPlusE2012A = triggerBits[19]; // HLT_Ele22_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v
    //bool TauPlusMu2012A = triggerBits[23]; // HLT_IsoMu20_eta2p1_LooseIsoPFTau20_v*
    bool TauPlusE2012B = triggerBits[17]; // HLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_v*
    bool TauPlusMu2012B = triggerBits[21]; // HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v*
    //bool TauPlusE2012D = triggerBits[16]; // HLT_Ele13_eta2p1_WP90Rho_LooseIsoPFTau20_v*
    //bool TauPlusMu2012D = triggerBits[20]; // HLT_IsoMu8_eta2p1_LooseIsoPFTau20_v*
    bool TauPlusETrigger = TauPlusE2012A || TauPlusE2012B;
    bool TauPlusMuTrigger = TauPlusMu2012A || TauPlusMu2012B;

    triggeredOn = TauPlusETrigger || TauPlusMuTrigger;
    if(TauPlusETrigger)
      chTags.push_back("TauPlusE");
    if(TauPlusMuTrigger)
      chTags.push_back("TauPlusMu");

    // Rest of Gen Particles
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
      std::cout << "llvvTauCollection Object NotFound" << std::endl;
      continue;
    }
    llvvTauCollection taus = *tauCollHandle;

    // Boosted tau Collection
    fwlite::Handle<llvvTauCollection> boostedTauCollHandle;
    boostedTauCollHandle.getByLabel(ev, "llvvObjectProducersUsed", "boosted");
    if(!boostedTauCollHandle.isValid())
    {
      std::cout << "llvvTauCollection Boosted Object NotFound" << std::endl;
    }
    else
    {
      llvvTauCollection boostedTaus = *boostedTauCollHandle;
      //for(size_t i = 0; i < boostedTaus.size(); ++i)
      //  taus.push_back(boostedTaus[i]);
    }

    // Jet Collection
    fwlite::Handle<llvvJetCollection> jetCollHandle;
    jetCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
    if(!jetCollHandle.isValid())
    {
      std::cout << "llvvJetCollection Object NotFound" << std::endl;
      continue;
    }
    llvvJetCollection jets_ = *jetCollHandle;
    llvvJetExtCollection jets;
    for(auto i = jetCollHandle->begin(); i != jetCollHandle->end(); ++i)
      jets.push_back(llvvJetExt(*i));

    // MET Collection
    fwlite::Handle<llvvMet> metHandle;
    metHandle.getByLabel(ev, "llvvObjectProducersUsed", "pfMETPFlow");
//    metHandle.getByLabel(ev, "llvvObjectProducersUsed", "pfType1CorrectedMet");
    if(!metHandle.isValid())
    {
      std::cout << "llvvMet Object NotFound" << std::endl;
      continue;
    }
    met = *metHandle;

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
    rho = *rhoHandle;

    // Rho25
    fwlite::Handle<double> rho25Handle;
    rho25Handle.getByLabel(ev, "kt6PFJetsCentral", "rho");
    if(!rho25Handle.isValid())
    {
      std::cout << "rho25 Object NotFound" << std::endl;
      continue;
    }
    rho25 = *rho25Handle;



    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Finished loading collections, now computing PU weight and trigger scale factor" << std::endl;
    #endif

    // Pileup Weight
    if(isMC)
    {
      if(isStauStau)
      {
        int nEvents = 10000;
        double xsec = stauCrossSec(stauMass, neutralinoMass);
        crossSection = xsec;
        xsecWeight  = xsec/nEvents;
      }
      puWeight     = LumiWeights->weight(genEv.ngenITpu) * PUNorm[0];
      weight       = xsecWeight*puWeight;
      weight_plus  = PuShifters[utils::cmssw::PUUP ]->Eval(genEv.ngenITpu) * (PUNorm[2]/PUNorm[0]);
      weight_minus = PuShifters[utils::cmssw::PUDOWN]->Eval(genEv.ngenITpu) * (PUNorm[1]/PUNorm[0]);
    }

    // Get trigger Scale Factor
    triggerSF = 1;
    if(isMC)
    {
      #if defined(DEBUG_EVENT)
      if(debugEvent)
      {
        myCout << " Event";
        if(TauPlusETrigger)
          myCout << ", it is a TauPlusE event";
        if(TauPlusMuTrigger)
          myCout << ", it is a TauPlusMu event";
        myCout << std::endl;

        if(triggeredOn)
        {
          myCout << "  Looping on leptons:" << std::endl;
          for(auto lep = leptons.begin(); lep != leptons.end(); ++lep)
            myCout << "    Lepton (" << lep->id << ", pT=" << lep->pt() << ") trigger bits: " << bitset<8*sizeof(int)>(lep->Tbits) << std::endl;

          myCout << "  Looping on taus:" << std::endl;
          for(auto tau = taus.begin(); tau != taus.end(); ++tau)
            myCout << "    Tau (pT=" << tau->pt() << ") trigger bits: " << bitset<8*sizeof(int)>(tau->Tbits) << std::endl;
        }
      }
      #endif

      if(TauPlusETrigger)
      {
        llvvTau* trigTau = NULL, *leadTau = NULL;
        llvvLepton* trigE = NULL, *leadE = NULL;

        // This is working, but sometimes it can't find the triggered lepton
        //so, when it can't be found, we select the leading pt lepton of the right type
        for(auto lep = leptons.begin(); lep != leptons.end(); ++lep)
        {
          if(lep->Tbits & (3 << 17))
          {
            if(trigE == NULL)
              trigE = &(*lep);
            else
              if(lep->pt() > trigE->pt())
                trigE = &(*lep);
          }

          if(abs(lep->id) == 11)
          {
            if(leadE == NULL)
              leadE = &(*lep);
            else
              if(lep->pt() > leadE->pt())
                leadE = &(*lep);
          }
        }
        if(trigE == NULL)
          trigE = leadE;

        // Tau trigger matching has not yet been enabled in the nTuple production (Tbits is filled with random data)
        //so we use the leading pt pt
        for(auto tau = taus.begin(); tau != taus.end(); ++tau)
        {
          //if(tau->Tbits & (3 << 17))
          //{
          //  trigTau = &(*tau);
          //  break;
          //}

          if(leadTau == NULL)
            leadTau = &(*tau);
          else
            if(tau->pt() > leadTau->pt())
              leadTau = &(*tau);
        }
        if(trigTau == NULL)
          trigTau = leadTau;

        if(trigTau != NULL && trigE != NULL)
        {
          triggerSF *= leptonTauTriggerScaleFactor(*trigE, *trigTau);
        }
        else
        {
          #if defined(DEBUG_EVENT)
          if(debugEvent)
          {
            if(trigE == NULL)
              myCout << " TauPlusE trigSF: Unable to find triggered electron" << std::endl;
            if(trigTau == NULL)
              myCout << " TauPlusE trigSF: Unable to find triggered tau" << std::endl;
          }
          #endif
        }
      }

      if(TauPlusMuTrigger)
      {
        llvvTau* trigTau = NULL, *leadTau = NULL;
        llvvLepton* trigMu = NULL, *leadMu = NULL;

        // This is working, but sometimes it can't find the triggered lepton
        //so, when it can't be found, we select the leading pt lepton of the right type
        for(auto lep = leptons.begin(); lep != leptons.end(); ++lep)
        {
          if(lep->Tbits & (3 << 21))
          {
            if(trigMu == NULL)
              trigMu = &(*lep);
            else
              if(lep->pt() > trigMu->pt())
                trigMu = &(*lep);
          }

          if(abs(lep->id) == 13)
          {
            if(leadMu == NULL)
              leadMu = &(*lep);
            else
              if(lep->pt() > leadMu->pt())
                leadMu = &(*lep);
          }
        }
        if(trigMu == NULL)
          trigMu = leadMu;

        // Tau trigger matching has not yet been enabled in the nTuple production (Tbits is filled with random data)
        //so we use the leading pt pt
        for(auto tau = taus.begin(); tau != taus.end(); ++tau)
        {
          //if(tau->Tbits & (3 << 21))
          //{
          //  trigTau = &(*tau);
          //  break;
          //}

          if(leadTau == NULL)
            leadTau = &(*tau);
          else
            if(tau->pt() > leadTau->pt())
              leadTau = &(*tau);
        }
        if(trigTau == NULL)
          trigTau = leadTau;

        if(trigTau != NULL && trigMu != NULL)
        {
          triggerSF *= leptonTauTriggerScaleFactor(*trigMu, *trigTau);
        }
        else
        {
          #if defined(DEBUG_EVENT)
          if(debugEvent)
          {
            if(trigMu == NULL)
              myCout << " TauPlusMu trigSF: Unable to find triggered muon" << std::endl;
            if(trigTau == NULL)
              myCout << " TauPlusMu trigSF: Unable to find triggered tau" << std::endl;
          }
          #endif
        }
      }

      #if defined(DEBUG_EVENT)
      if(debugEvent)
        myCout << "  Computed trigger SF: " << triggerSF << std::endl;
      #endif
    }
    if(applyScaleFactors && isMC)
      weight *= triggerSF;



    #if defined(DEBUG_EVENT)
    if(debugEvent)
    {
      myCout << " Finished computing PU weight and trigger scale factors" << std::endl;
      myCout << " Getting leading lepton" << std::endl;
    }
    #endif
    // Get Leading Lepton
    for(size_t i = 0; i < leptons.size(); ++i)
    {
      int lepId = leptons[i].id;

      if(abs(lepId) == 13 && muCor)
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
      bool keepKin = true; // We want to keep leptons with a looser selection than the one for the ID so we can then remove events with multiple leptons
      double eta = (lepId == 11)?(leptons[i].electronInfoRef->sceta):(leptons[i].eta());
      if(lepId == 11) // If Electron
      {
        if(leptons[i].pt() < minElPt)
          passKin = false;
        if(abs(eta) > maxElEta)
          passKin = false;

        if(leptons[i].pt() < 10)
          keepKin = false;
        if(abs(eta) > 2.4)
          keepKin = false;

        if(abs(eta) > ECALGap_MinEta && abs(eta) < ECALGap_MaxEta) // Remove electrons that fall in ECAL Gap
        {
          passKin = false;
          keepKin = false;
        }
      }
      else // If Muon
      {
        if(leptons[i].pt() < minMuPt)
          passKin = false;
        if(abs(eta) > maxMuEta)
          passKin = false;

        if(leptons[i].pt() < 10)
          keepKin = false;
        if(abs(eta) > 2.4)
          keepKin = false;
      }

      // Lepton ID
      bool passID = true;
      Int_t idbits = leptons[i].idbits;
      if(lepId == 11)
      {
        // bool isTight = electronMVAID(leptons[i].electronInfoRef->mvanontrigv0, leptons[i], TightID);
        // bool isLoose = electronMVAID(leptons[i].electronInfoRef->mvanontrigv0, leptons[i], LooseID);
        // bool isLoose = ((idbits >> 4) & 0x1);
        // bool isTight = ((idbits >> 6) & 0x1);
        passID = electronMVAID(leptons[i].electronInfoRef->mvanontrigv0, leptons[i], LooseID);
        if(leptons[i].d0 > 0.045)
          passID = false;
        if(leptons[i].dZ > 0.2)
          passID = false;
        if(leptons[i].electronInfoRef->isConv)
          passID = false;
        if(leptons[i].trkLostInnerHits > 0)
          passID = false;
      }
      else
      {
        // bool isLoose = ((idbits >> 8) & 0x1);
        // bool isTight = ((idbits >> 10) & 0x1);
        passID = ((idbits >> 8) & 0x1);
        if(leptons[i].d0 > 0.045)
          passID = false;
        if(leptons[i].dZ > 0.2)
          passID = false;
      }

      // Lepton Isolation
      bool passIso = true;
      double relIso = utils::cmssw::relIso(leptons[i], rho);
      if(lepId == 11)
      {
        if(relIso > 0.3)
          passIso = false;
      }
      else
      {
        if(relIso > 0.3)
          passIso = false;
      }

      // Keep desired leptons
      if(keepKin && passID && passIso)
        selLeptons.push_back(leptons[i]);
      if(!triggeredOn)
        continue;

      mon.fillHisto("leptonCutFlow", chTags, 0, weight);
      if(passID)
      {
        mon.fillHisto("leptonCutFlow", chTags, 1, weight);
        if(passKin)
        {
          mon.fillHisto("leptonCutFlow", chTags, 2, weight);
          if(passIso)
            mon.fillHisto("leptonCutFlow", chTags, 3, weight);
        }
      }
    }

    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Getting taus" << std::endl;
    #endif
    // Get taus
    for(size_t i = 0; i < taus.size(); ++i)
    {
      llvvTau& tau = taus[i];

      // Tau Kinematics
      bool passKin = true;
      if(tau.pt() < minTauPt)
        passKin = false;
      if(abs(tau.eta()) > maxTauEta)
        passKin = false;

      // Tau overlap with selected leptons
      bool passIso = true;
      for(size_t lep = 0; lep < selLeptons.size(); ++lep)
      {
        int lepId = abs(selLeptons[lep].id);
        if(lepId == 11)
        {
          if(selLeptons[lep].pt() < minElPt)
            continue;
          if(abs(selLeptons[lep].dZ) > 0.1)
            continue;
          double eta = selLeptons[lep].electronInfoRef->sceta;
          if(abs(eta) > maxElEta)
            continue;
          double relIso = utils::cmssw::relIso(selLeptons[lep], rho);
          if(relIso > 0.1)
            continue;
        }
        else
        {
          if(selLeptons[lep].pt() < minMuPt)
            continue;
          if(abs(selLeptons[lep].eta()) > maxMuEta)
            continue;
          double relIso = utils::cmssw::relIso(selLeptons[lep], rho);
          if(relIso > 0.1)
            continue;
          Int_t idbits = selLeptons[lep].idbits;
          bool isTight = ((idbits >> 10) & 0x1);
          if(!isTight)
            continue;
        }
        if(deltaR(tau, selLeptons[lep]) < 0.1)
        {
          passIso = false;
          break;
        }
      }

      bool passQual = true;
      if(abs(tau.dZ) > 0.5)
        passQual = false;

      // Tau ID
      bool passID = true;
      if(!tau.passId(llvvTAUID::decayModeFinding)) passID = false;
      if(!doTightTauID)
      {
        if(!tau.passId(llvvTAUID::byMediumCombinedIsolationDeltaBetaCorr3Hits)) passID = false;
      }
      else
      {
        if(!tau.passId(llvvTAUID::byTightCombinedIsolationDeltaBetaCorr3Hits)) passID = false;
      }
      if(!tau.passId(llvvTAUID::againstMuonTight3)) passID = false;
      if(!tau.passId(llvvTAUID::againstElectronMediumMVA5)) passID = false;

      if(passID && passKin && passIso && passQual && tau.isPF)
        selTaus.push_back(tau);
      if(!triggeredOn)
        continue;

      // Fill control histograms
      mon.fillHisto("tauCutFlow", chTags, 0, weight);
      if(tau.isPF)
      {
        mon.fillHisto("tauCutFlow", chTags, 1, weight);
        mon.fillHisto("tauID", chTags, 0, weight);
        if((doTightTauID && tau.passId(llvvTAUID::byTightCombinedIsolationDeltaBetaCorr3Hits)) || (!doTightTauID && tau.passId(llvvTAUID::byMediumCombinedIsolationDeltaBetaCorr3Hits)))
        {
          mon.fillHisto("tauID", chTags, 1, weight);
          if(tau.passId(llvvTAUID::decayModeFinding))
          {
            mon.fillHisto("tauID", chTags, 2, weight);
            if(tau.passId(llvvTAUID::againstElectronMediumMVA5))
            {
              mon.fillHisto("tauID", chTags, 3, weight);
              if(tau.passId(llvvTAUID::againstMuonTight3))
                mon.fillHisto("tauID", chTags, 4, weight);
            }
          }
        }
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
    }

    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Getting jets" << std::endl;
    #endif
    // Get Jets
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

      // Jet ID
      bool passID = true;
      Int_t idbits = jets[i].idbits;
      bool passPFLoose = (idbits & 0x01);
      int fullPuId = (idbits >> 3) & 0x0f;
      bool passLooseFullPuId = ((fullPuId >> 2) & 0x01);
      passID = passLooseFullPuId;

      // Jet Kinematics
      bool passKin = true;
      if(abs(jets[i].eta()) > maxJetEta)
        passKin = false;

      // B-jets
      bool isBJet = false;
//      bool hasBtagCorr = false;
      if(jets[i].csv > 0.679)
      {
        isBJet = true;
//        hasBtagCorr = true;
      }

      if(isMC)
      {
      }

      // Compute scale and resolution uncertainties
      if(isMC)
      {
        std::vector<float> smearPt = utils::cmssw::smearJER(jets[i].pt(),jets[i].eta(),jets[i].genj.pt());
        jets[i].jer = smearPt[0];
        jets[i].jerup = smearPt[1];
        jets[i].jerdown = smearPt[2];
        smearPt = utils::cmssw::smearJES(jets[i].pt(),jets[i].eta(), totalJESUnc);
        jets[i].jesup = smearPt[0];
        jets[i].jesdown = smearPt[1];
      }
      else
      {
        jets[i].jer = jets[i].pt();
        jets[i].jerup = jets[i].pt();
        jets[i].jerdown = jets[i].pt();
        jets[i].jesup = jets[i].pt();
        jets[i].jesdown = jets[i].pt();
      }

      if(passPFLoose && passID && passKin)
      {
        selJets.push_back(jets[i]);
        selJetsOut.push_back(jets_[i]);
      }
      if(passPFLoose && passID && passKin && isBJet)
        selBJets.push_back(jets_[i]);
      if(!triggeredOn)
        continue;

      // Fill Jet control histograms
      mon.fillHisto("jetCutFlow", chTags, 0, weight);
      if(passPFLoose)
      {
        mon.fillHisto("jetCutFlow", chTags, 1, weight);
        if(passID)
        {
          mon.fillHisto("jetCutFlow", chTags, 2, weight);
          if(passKin)
            mon.fillHisto("jetCutFlow", chTags, 5, weight);
        }
      }
    }

    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Sorting leptons, taus and jets" << std::endl;
    #endif
    if(selLeptons.size() != 0)
    {
      std::sort(selLeptons.begin(), selLeptons.end(), sort_llvvObjectByPt);

      if(abs(selLeptons[0].id) == 11)
        chTags.push_back("leadingE");
      else
        chTags.push_back("leadingMu");
    }

    if(selTaus.size() != 0)
      std::sort(selTaus.begin(), selTaus.end(), sort_llvvObjectByPt);

    nBJets = selBJets.size();
    nJets = selJets.size();
    if(nJets != 0)
      std::sort(selJets.begin(), selJets.end(), sort_llvvObjectByPt);



    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Requiring an opposite sign pair" << std::endl;
    #endif
    // Opposite Sign requirements
    maxPtSum = 0;
    tauIndex = -1;
    leptonIndex = -1;
    for(size_t i = 0; i < selLeptons.size(); ++i)
    {
      if(abs(selLeptons[i].id) == 11) // Electron
      {
        if(selLeptons[i].pt() < minElPt)
          continue;
        if(abs(selLeptons[i].dZ) > 0.1)
          continue;
        double eta = selLeptons[i].electronInfoRef->sceta;
        if(abs(eta) > maxElEta)
          continue;
        double relIso = utils::cmssw::relIso(selLeptons[i], rho);
        if(relIso > 0.1)
          continue;
      }
      else
      {
        if(selLeptons[i].pt() < minMuPt)
          continue;
        if(abs(selLeptons[i].eta()) > maxMuEta)
          continue;
        double relIso = utils::cmssw::relIso(selLeptons[i], rho);
        if(relIso > 0.1)
          continue;
        Int_t idbits = selLeptons[i].idbits;
        bool isTight = ((idbits >> 10) & 0x1);
        if(!isTight)
          continue;
      }

      for(size_t j = 0; j < selTaus.size(); ++j)
      {
        double PtSum = selLeptons[i].pt() + selTaus[j].pt();
        if(selLeptons[i].id * selTaus[j].id < 0)
        {
          if(PtSum > maxPtSum)
          {
            maxPtSum = PtSum;
            tauIndex = j;
            leptonIndex = i;
            isOS = true;
            tauLeadPt = selTaus[tauIndex].pt();
            lepLeadPt = selLeptons[leptonIndex].pt();
//            triggerSF = leptonTauTriggerScaleFactor(selLeptons[leptonIndex], selTaus[tauIndex]);
            leptonIdIsoSF = leptonIdAndIsoScaleFactor(selLeptons[leptonIndex]);
            tauSF = ::tauSF(selTaus[tauIndex], gen, antiEMva5Medium);
          }
        }
        if(PtSum < maxPtSum) // Skip a few iterations if it is not expected that we will find a better candidate
          break;
      }
    }
    if(!isMC)
    {
//      triggerSF = 1;
      leptonIdIsoSF = 1;
      tauSF = 1;
    }
    if(isOS && applyScaleFactors)
    {
//      weight *= triggerSF;
      weight *= leptonIdIsoSF;
      weight *= tauSF;
    }

    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Rejecting event if multilepton" << std::endl;
    #endif
    // Reject events with more leptons
    isMultilepton = false;
    if(isOS)
    {
      if(abs(selLeptons[leptonIndex].id) == 11) //Electron tau channel
      {
        for(size_t i  = 0; i < selLeptons.size(); ++i)
        {
          if(i == (size_t)leptonIndex)
            continue;
          if(abs(selLeptons[i].id) != 11)
            continue;
          if(selLeptons[i].pt() < 20)
            continue;
          if(abs(selLeptons[i].electronInfoRef->sceta) > 2.1)
            continue;
          isMultilepton = true;
          break;
        }
      }
      else // Muon tau channel
      {
        for(size_t i = 0; i < selLeptons.size(); ++i)
        {
          if(i == (size_t)leptonIndex)
            continue;
          isMultilepton = true;
          break;
        }
      }

      // Set up channels
      if(!isMultilepton)
      {
        if(abs(selLeptons[leptonIndex].id) == 11)
          chTags.push_back("etau");
        else
          chTags.push_back("mutau");

        #if defined(DEBUG_EVENT)
        if(debugEvent || true)
        {
          bool isMuTau = false;
          bool isAll = false;
          for(auto tag = chTags.begin(); tag != chTags.end(); ++tag)
          {
            if(*tag == "all")
              isAll = true;
            if(*tag == "mutau")
              isMuTau = true;
          }

          if(!isAll)
          {
            myCout << " Found a ";
            if(isMuTau)
              myCout << "mu-tau";
            else
              myCout << "e-tau";
            myCout << " event";
            if(isAll)
              myCout << " and it is tagged as all" << std::endl;
            else
              myCout << " but it is not tagged as all" << std::endl;
          }
        }
        #endif
      }
    }

    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Computing angular variables" << std::endl;
    #endif
    // Get angular variables
    if(isOS && !isMultilepton)
    {
      /**    LAB FRAME    **/
      TLorentzVector lep(selLeptons[leptonIndex].Px(), selLeptons[leptonIndex].Py(), selLeptons[leptonIndex].Pz(), selLeptons[leptonIndex].E());
      TLorentzVector tau(selTaus[tauIndex].Px(), selTaus[tauIndex].Py(), selTaus[tauIndex].Pz(), selTaus[tauIndex].E());
      TLorentzVector Tmet(met.Px(), met.Py(), met.Pz(), met.E());

      deltaAlphaLepTau = lep.Angle(tau.Vect());
      deltaPhiLepTauMET = Tmet.DeltaPhi(lep + tau);
      deltaPhiLepTau = lep.DeltaPhi(tau);

      minDeltaPhiMETJetPt40 = 5;
      for(auto &jet : selJets)
      {
        TLorentzVector tJet(jet.Px(), jet.Py(), jet.Pz(), jet.E());

        if(tJet.Pt() < 40)
          break;

        double temp = Tmet.DeltaPhi(tJet);
        if(temp < minDeltaPhiMETJetPt40)
          minDeltaPhiMETJetPt40 = temp;
      }

      double posSign = Tmet.CosTheta();
      cosThetaTau = tau.CosTheta();
      cosThetaLep = lep.CosTheta();
      if(posSign < 0)
      {
        cosThetaTau = -cosThetaTau;
        cosThetaLep = -cosThetaLep;
      }

      /**   CS FRAME   **/
      TLorentzVector tauSystem = tau + lep;
      TLorentzVector tauCS = tau, lepCS = lep, metCS = Tmet;
      TLorentzVector beam1(0, 0, 0, 0);
      TLorentzVector beam2(0, 0, 0, 0);
      double energy = sqrtS * 500.; // sqrtS / 2 * 1000
      double mom = sqrt(energy*energy + 0.938*0.938);
      if(posSign > 0)
      {
        beam1.SetPxPyPzE(0, 0,  mom, energy);
        beam2.SetPxPyPzE(0, 0, -mom, energy);
      }
      else
      {
        beam1.SetPxPyPzE(0, 0, -mom, energy);
        beam2.SetPxPyPzE(0, 0,  mom, energy);
      }

      TVector3 boost = -tauSystem.BoostVector();
      //tauSystem.Boost(boost); // By construction, this will be 0
      tauCS.Boost(boost);
      lepCS.Boost(boost);
      metCS.Boost(boost);
      beam1.Boost(boost);
      beam2.Boost(boost);

      TLorentzVector SQA = beam1 - beam2; // Spin quantization axis
      TRotation rotation;

      TVector3 newZAxis = SQA.Vect().Unit();
      TVector3 targetZaxis(0, 0, 1);
      TVector3 rotAxis = targetZaxis.Cross(newZAxis);
      double rotAngle = targetZaxis.Angle(newZAxis);
      rotation.Rotate(-rotAngle, rotAxis);

      SQA.Transform(rotation);
      tauCS.Transform(rotation);
      lepCS.Transform(rotation);
      metCS.Transform(rotation);
      beam1.Transform(rotation);
      beam2.Transform(rotation);

      cosThetaCS = lepCS.CosTheta();
      deltaPhiLepMETCS = metCS.DeltaPhi(lepCS);
    }

    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Computing pair mass" << std::endl;
    #endif
    // Tau-Lepton pair mass calculation
    isSVfit = doSVfit;
    if(isOS && !isMultilepton)
    {
      auto selLepton = selLeptons[leptonIndex];
      auto selTau    = selTaus[tauIndex];
      invMass = (selLepton+selTau).M(); // Invariant mass

      if(doSVfit)
      {
        TMatrixD covMET(2, 2);
        std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
        svFitStandalone::Vector measuredMET(met.px(), met.py(), 0);

        covMET[0][0] = met.sigx2;
        covMET[0][1] = met.sigxy;
        covMET[1][0] = met.sigxy;
        covMET[1][1] = met.sigy2;

        if(abs(selLepton.id) == 11)
          measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToElecDecay, svFitStandalone::LorentzVector(selLepton.px(), selLepton.py(), selLepton.pz(), selLepton.E())));
        else
          measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToMuDecay, svFitStandalone::LorentzVector(selLepton.px(), selLepton.py(), selLepton.pz(), selLepton.E())));
        measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToHadDecay, svFitStandalone::LorentzVector(selTau.px(), selTau.py(), selTau.pz(), selTau.E())));

        SVfitStandaloneAlgorithm SVfit_algo(measuredTauLeptons, measuredMET, covMET, 0);
        //SVfit_algo.maxObjFunctionCalls(10000) // To change the max number of iterations before minimization is terminated, default 5000
        SVfit_algo.addLogM(false); // To not use the LogM penalty, it is used by default
        //SVfit_algo.metPower(0.5); // Additional power to enhance MET likelihood, default is 1.
        //SVfit_algo.fit();
        //SVfit_algo.integrate();
        SVfit_algo.integrateVEGAS();
        //SVfit_algo.integrateMarkovChain();
        if(SVfit_algo.isValidSolution())
          mass = SVfit_algo.mass();
        else
          isSVfit = false;
      }
    }

    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Computing MT and deconstructed MT" << std::endl;
    #endif
    // MT and deconstructed MT calculation
    if(isOS && !isMultilepton && (!doSVfit || isSVfit))
    {
      auto& selLepton = selLeptons[leptonIndex];
      double cosDeltaPhi = cos(deltaPhi(selLepton.phi(), met.phi()));
      double fac = 2 * met.pt() * selLepton.pt();

      mt = sqrt(fac * (1 - cosDeltaPhi));
      Q80 = 1 - (80.0*80.0) / fac;
      Q100 = 1 - (100.0*100.0) / fac;
      cosPhi = cosDeltaPhi;
    }

    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Computing MT2" << std::endl;
    #endif
    // MT2 calculation
    if(isOS && !isMultilepton && (!doSVfit || isSVfit))
    {
      auto selLepton = selLeptons[leptonIndex];
      auto selTau    = selTaus[tauIndex];
      double pa[3];
      double pb[3];
      double pmiss[3];
      double mn;

      pa[0] = selLepton.M();
      pa[1] = selLepton.px();
      pa[2] = selLepton.py();
      pb[0] = selTau.M();
      pb[1] = selTau.px();
      pb[2] = selTau.py();
      pmiss[0] = 0;
      pmiss[1] = met.px();
      pmiss[2] = met.py();
      mn = 0;

      mt2_bisect::mt2 mt2_event;
      mt2_event.set_momenta(pa,pb,pmiss);
      mt2_event.set_mn(mn);
      mt2 = mt2_event.get_mt2();

/*      mn = 50;
      mt2_event.set_momenta(pa,pb,pmiss);
      mt2_event.set_mn(mn);
      mt2_50 = mt2_event.get_mt2();

      mn = 150;
      mt2_event.set_momenta(pa,pb,pmiss);
      mt2_event.set_mn(mn);
      mt2_150 = mt2_event.get_mt2(); // */
    }

    bool stauPlot = false;
    if(stauMass == stauMtoPlot && neutralinoMass == neutralinoMtoPlot)
      stauPlot = true;

    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Filling histograms" << std::endl;
    #endif
    bool plotThisEvent = !isStauStau || stauPlot;
    if(plotThisEvent)
    {
      //mon.fillHisto("nvtxAll", chTags, nvtx, weight);
      if(triggeredOn)
      {
        mon.fillHisto("eventflow", chTags, 0, weight);
        if(met.pt() > 30)
        {
          mon.fillHisto("eventflow", chTags, 1, weight);
          if(selLeptons.size() > 0)
          {
            mon.fillHisto("eventflow", chTags, 2, weight);
            mon.fillHisto("nbjets", chTags, selBJets.size(), weight);
            if(selTaus.size() > 0)
            {
              mon.fillHisto("eventflow", chTags, 3, weight);
              if(selBJets.size() == 0)
              {
                mon.fillHisto("eventflow", chTags, 4, weight);
                if(isOS)
                {
                  mon.fillHisto("eventflow", chTags, 5, weight);
                  if(!isMultilepton)
                  {
                    mon.fillHisto("eventflow", chTags, 6, weight);
                    if(!doSVfit || isSVfit)
                    {
                      mon.fillHisto("eventflow", chTags, 7, weight);

                      mon.fillHisto("nvtx", chTags, nvtx, weight);
                      mon.fillHisto("nvtxraw", chTags, nvtx, weight/puWeight);
                      mon.fillHisto("nup", "", genEv.nup, 1);

                      mon.fillHisto("rho", chTags, rho, weight);
                      mon.fillHisto("rho25", chTags, rho25, weight);

                      mon.fillHisto("MET", chTags, met.pt(), weight);

                      mon.fillHisto("MT", chTags, mt, weight);
                      mon.fillHisto("Q80", chTags, Q80, weight);
                      mon.fillHisto("Q100", chTags, Q100, weight);
                      mon.fillHisto("cosPhi", chTags, cosPhi, weight);
                      mon.fillHisto("Q80VsCosPhi", chTags, cosPhi, Q80, weight);
                      mon.fillHisto("Q100VsCosPhi", chTags, cosPhi, Q100, weight);

                      mon.fillHisto("MT2", chTags, mt2, weight);
                      if(doSVfit)
                        mon.fillHisto("SVFitMass", chTags, mass, weight);
                      mon.fillHisto("InvMass", chTags, invMass, weight);

                      mon.fillHisto("deltaAlphaLepTau", chTags, deltaAlphaLepTau, weight);
                      mon.fillHisto("deltaPhiLepTauMET", chTags, deltaPhiLepTauMET, weight);
                      mon.fillHisto("deltaPhiLepTau", chTags, deltaPhiLepTau, weight);
                      mon.fillHisto("cosThetaTau", chTags, cosThetaTau, weight);
                      mon.fillHisto("cosThetaLep", chTags, cosThetaLep, weight);
                      mon.fillHisto("cosThetaCS", chTags, cosThetaCS, weight);
                      mon.fillHisto("deltaPhiLepMETCS", chTags, deltaPhiLepMETCS, weight);
                      mon.fillHisto("minDeltaPhiMETJetPt40", chTags, minDeltaPhiMETJetPt40, weight);

                      mon.fillHisto("metVSPTl", chTags, selLeptons[leptonIndex].pt(), met.pt(), weight);
                      mon.fillHisto("metVSPTtau", chTags, selTaus[tauIndex].pt(), met.pt(), weight);

                      mon.fillHisto("nlep", chTags, selLeptons.size(), weight);
                      double eta = selLeptons[leptonIndex].eta();
                      if(abs(selLeptons[leptonIndex].id) == 11) eta = selLeptons[leptonIndex].electronInfoRef->sceta;
                      mon.fillHisto("etaSelectedLep", chTags, eta, weight);
                      mon.fillHisto("ptSelectedLep", chTags, selLeptons[leptonIndex].pt(), weight);
                      mon.fillHisto("chargeSelectedLep", chTags, (selLeptons[leptonIndex].id > 0)?(-1):(1), weight);

                      mon.fillHisto("ntaus", chTags, selTaus.size(), weight);
                      mon.fillHisto("ptSelectedTau", chTags, selTaus[tauIndex].pt(), weight);
                      mon.fillHisto("ptSelectedTauExtended", chTags, selTaus[tauIndex].pt(), weight);
                      mon.fillHisto("etaSelectedTau", chTags, selTaus[tauIndex].eta(), weight);
                      mon.fillHisto("chargeSelectedTau", chTags, (selTaus[tauIndex].id > 0)?(-1):(1), weight);
                      mon.fillHisto("emfracSelectedTau", chTags, selTaus[tauIndex].emfraction, weight);
                      mon.fillHisto("dzSelectedTau", chTags, selTaus[tauIndex].dZ, weight);

                      mon.fillHisto("njets", chTags, selJets.size(), weight);
                      if(selJets.size() != 0)
                      {
                        mon.fillHisto("jetleadpt", chTags, selJets[0].pt(), weight);
                        mon.fillHisto("jetleadeta", chTags, selJets[0].eta(), weight);
                      }

                      for(auto i = selJets.begin(); i != selJets.end(); ++i)
                      {
                        mon.fillHisto("jetcsv", chTags, i->origcsv, weight);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    if(triggeredOn && met.pt() > 30 && selLeptons.size() > 0 && selBJets.size() == 0 && selTaus.size() > 0 && isOS && !isMultilepton && (!doSVfit || isSVfit))
      selected = true;

    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Filling TTree" << std::endl;
    NAN_WARN(rho)
    NAN_WARN(rho25)
    NAN_WARN(weight)
    NAN_WARN(weight_plus)
    NAN_WARN(weight_minus)
    NAN_WARN(puWeight)
    NAN_WARN(triggerSF)
    NAN_WARN(leptonIdIsoSF)
    NAN_WARN(tauSF)
    NAN_WARN(mass)
    NAN_WARN(invMass)
    NAN_WARN(mt)
    NAN_WARN(Q80)
    NAN_WARN(Q100)
    NAN_WARN(cosPhi)
    NAN_WARN(mt2)
    NAN_WARN(stauMass)
    NAN_WARN(neutralinoMass)
    NAN_WARN(deltaAlphaLepTau)
    NAN_WARN(deltaPhiLepTauMET)
    NAN_WARN(deltaPhiLepTau)
    NAN_WARN(cosThetaTau)
    NAN_WARN(cosThetaLep)
    NAN_WARN(deltaPhiLepMETCS)
    NAN_WARN(cosThetaCS)
    NAN_WARN(minDeltaPhiMetJet40)
    NAN_WARN(tauLeadPt)
    NAN_WARN(lepLeadPt)
    NAN_WARN(maxPtSum)
    #endif
    if(saveSummaryTree)
    {
      TDirectory* cwd = gDirectory;
      summaryOutFile->cd();
      summaryTree->Fill();
      cwd->cd();
    }

    #if defined(DEBUG_EVENT)
    if(debugEvent)
      myCout << " Finished event " << iev << std::endl;
    #endif
//    break;
    if(limit != 0)
      if(iev >= limit - 1)
        break;
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
  std::cout << "Saving results in " << outUrl << std::endl;
  TFile* outfile = new TFile(outUrl.c_str(), "RECREATE");
  mon.Write();
  outfile->Close();
  delete outfile;

  if(saveSummaryTree)
  {
//    std::string summaryOutUrl = outUrl;
//    summaryOutUrl.replace(summaryOutUrl.find(".root", 0), 5, "_summary.root");
//    std::cout << "Saving summary results in " << summaryOutUrl << std::endl;
//    summaryOutFile = new TFile(summaryOutUrl.c_str(), "RECREATE");
    // Write Tuple/Tree
//    summaryTree->SetDirectory(summaryOutFile);
    TDirectory* cwd = gDirectory;
    summaryOutFile->cd();
    summaryTree->Write();
    summaryOutFile->Close();
    delete summaryOutFile;
    cwd->cd();
  }

  return 0;
}

int old_main(int argc, char* argv[])
{
  /***************************************************************************/
  /*                          Global Initialization                          */
  /***************************************************************************/
  // Check arguments:
  if(argc<2)
    std::cout << "Usage: " << argv[0] << " parameters_cfg.py" << std::endl, exit(1);

  mt2_bisect::mt2 mt2_event;
  // Format: M, px, py
  double pa[3] = { 0.106, 39.0, 12.0 };
  double pb[3] = { 0.106, 119.0, -33.0 };
  double pmiss[3] = { 0, -29.9, 35.9 };
  double mn    = 50.; mn = 0; // Using same neutralino mass as IPM
  //mt2_event.set_momenta(pa,pb,pmiss);
  //mt2_event.set_mn(mn);
  //std::cout << std::endl << " mt2 = " << mt2_event.get_mt2() << std::endl;

  // Load FWLite:
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

  // Configure the process (aka "Load parameters from configuration file")
  const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");

  bool isMC = runProcess.getParameter<bool>("isMC");
  double xsec = runProcess.getParameter<double>("xsec");
  std::vector<std::string> urls = runProcess.getParameter<std::vector<std::string> >("input");
  std::string baseDir = runProcess.getParameter<std::string>("dirName");
  std::string outdir = runProcess.getParameter<std::string>("outdir");
  std::string jecDir = runProcess.getParameter<std::string>("jecDir");
//  bool runSystematics = runProcess.getParameter<bool>("runSystematics");
  bool saveSummaryTree = runProcess.getParameter<bool>("saveSummaryTree");
  bool exclusiveRun = runProcess.getParameter<bool>("exclusiveRun");
  std::string periodHLT = "2012B";  // This variable should hold the period of the HLT triggers wished to be used,
                                    //possible values are: 2012A, 2012B, 2012C, 2012D. 2012B and 2012C are equal.
  if(runProcess.exists("periodHLT"))
  {
    periodHLT = runProcess.getParameter<std::string>("periodHLT");
    if(periodHLT == "2012C")
      periodHLT = "2012B";
    if(periodHLT != "2012A" && periodHLT != "2012B" && periodHLT != "2012D")
      periodHLT = "2012B";

    std::cout << "Using the HLT from " << periodHLT << std::endl;
  }
  double stauMtoPlot = 120;
  double neutralinoMtoPlot = 20;  // TStauStau mass point to plot by default
  if(runProcess.exists("stauMtoPlot"))
    stauMtoPlot = runProcess.getParameter<double>("stauMtoPlot");
  if(runProcess.exists("neutralinoMtoPlot"))
    neutralinoMtoPlot = runProcess.getParameter<double>("neutralinoMtoPlot");
  std::string selection = "LIP";
  if(runProcess.exists("selection"))
    selection = runProcess.getParameter<std::string>("selection");
  if(selection != "LIP")
  {
    // All recognised selection schemes should be added here so thay are valid options later in the code
    if(selection != "IPM")
      selection = "LIP";
  }
  bool doSVfit = false;
  if(runProcess.exists("doSVfit"))
    doSVfit = runProcess.getParameter<bool>("doSVfit");
  bool doLIPTauID = false;
  if(runProcess.exists("doLIPTauID"))
    doLIPTauID = runProcess.getParameter<bool>("doLIPTauID");
  bool doLIPBVeto = false;
  if(runProcess.exists("doLIPBVeto"))
    doLIPBVeto = runProcess.getParameter<bool>("doLIPBVeto");

  bool doCombined = false;
  if(runProcess.exists("doCombined"))
    doCombined = runProcess.getParameter<bool>("doCombined");
  if(doCombined)
    doLIPTauID = true;
  bool doTauAgainstElectronLoose = false;
  if(runProcess.exists("doTauAgainstElectronLoose"))
    doTauAgainstElectronLoose = runProcess.getParameter<bool>("doTauAgainstElectronLoose");
  bool noTauOverlap = false;
  if(runProcess.exists("noTauOverlap"))
    noTauOverlap = runProcess.getParameter<bool>("noTauOverlap");
  //bool doQuickFix = false;
  //if(runProcess.exists("doQuickFix"))
  //  doQuickFix = runProcess.getParameter<bool>("doQuickFix");

  // Hardcoded configs
  double sqrtS          =  8;
  double minElPt        = 25;
  double maxElEta       =  2.5;
  double ECALGap_MinEta =  1.4442;
  double ECALGap_MaxEta =  1.5660;
  double minMuPt        = 20;
  double maxMuEta       =  2.4;
  double minTauPt       = 25;
  double maxTauEta      =  2.3;
  double minJetPt       = 30;
  double maxJetEta      =  2.5;
  double minTauJetPt    = 20;
  double maxTauJetEta   =  2.5;
  if(periodHLT == "2012D")
  {
    minMuPt = 12;
    minElPt = 17;
  }
  if(selection == "IPM")
  {
    minElPt      = 20;
    maxElEta     =  2.1;
    minMuPt      = 20;
    maxMuEta     =  2.1;
    minTauPt     = 20;
    maxTauEta    =  2.3;
    //minJetPt     = ;
    maxJetEta    =  4.7;
    //minTauJetPt  = ;
    //maxTauJetEta = ;
  }


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
  std::string outUrl = outdir;

  TString turl(url);
//  bool isSingleMuPD(!isMC && turl.Contains("SingleMu"));
  bool isV0JetsMC(isMC && (turl.Contains("DYJetsToLL_50toInf") || turl.Contains("WJets")));
  bool isStauStau(isMC && turl.Contains("TStauStau"));

  TTree* summaryTree = NULL;
  if(saveSummaryTree)
    summaryTree = new TTree("Events", "Events");



  /***************************************************************************/
  /*                         Initializing Histograms                         */
  /***************************************************************************/
  SmartSelectionMonitor mon;
  TH1D *eventflow = (TH1D*)mon.addHistogram(new TH1D("eventflow", ";;Events", 5, 0, 5));
  eventflow->GetXaxis()->SetBinLabel(1, "HLT");
  eventflow->GetXaxis()->SetBinLabel(2, "> 1l");
  eventflow->GetXaxis()->SetBinLabel(3, "B-veto");
  eventflow->GetXaxis()->SetBinLabel(4, "> 1#tau");
  eventflow->GetXaxis()->SetBinLabel(5, "OS");
  // ...
  // TH2D* hist = (TH2D*)mon.addHistogram(...);

  mon.addHistogram(new TH1D("nup",     ";NUP;Events", 10, 0, 10));
  //mon.addHistogram(new TH1D("nupfilt", ";NUP;Events", 10, 0, 10));

  // PU
  mon.addHistogram(new TH1D("nvtxAll",    ";Vertices;Events",       50, -0.5, 49.5));
  mon.addHistogram(new TH1D("nvtxHLT",    ";Vertices;Events",       50, -0.5, 49.5));
  mon.addHistogram(new TH1D("nvtx1l",    ";Vertices;Events",       50, -0.5, 49.5));
  mon.addHistogram(new TH1D("nvtxBveto",    ";Vertices;Events",       50, -0.5, 49.5));
  mon.addHistogram(new TH1D("nvtx1tau",    ";Vertices;Events",       50, -0.5, 49.5));
  mon.addHistogram(new TH1D("nvtxOS",    ";Vertices;Events",       50, -0.5, 49.5));

  mon.addHistogram(new TH1D("nvtx",    ";Vertices;Events",       50, -0.5, 49.5));
  mon.addHistogram(new TH1D("nvtxraw", ";Vertices;Events",       50, -0.5, 49.5));
  mon.addHistogram(new TH1D("rho",     ";#rho;Events",           25,  0,   25));
  mon.addHistogram(new TH1D("rho25",   ";#rho(#eta<2.5);Events", 25,  0,   25));


  // Leptons
  mon.addHistogram(new TH1D("nlep",       ";nlep;Events",       10,  0,   10));
  mon.addHistogram(new TH1D("leadpt",     ";p_{T}^{l};Events",  50,  0,  100));
  mon.addHistogram(new TH1D("leadeta",    ";#eta^{l};Events",   25, -2.6,  2.6));
  mon.addHistogram(new TH1D("leadcharge", ";q^{l};Events",       5, -2,    2));
  TH1D *leptonCutFlow   = (TH1D*)mon.addHistogram(new TH1D("leptonCutFlow",   ";;Leptons", 4, 0, 4));
  leptonCutFlow->GetXaxis()->SetBinLabel(1, "All");
  leptonCutFlow->GetXaxis()->SetBinLabel(2, "ID");
  leptonCutFlow->GetXaxis()->SetBinLabel(3, "Kin");
  leptonCutFlow->GetXaxis()->SetBinLabel(4, "Iso");

  // Lepton Isolation
  mon.addHistogram(new TH1D("isomu",  "RelIso(#mu);;Leptons", 100, -0.5, 9.5));
  mon.addHistogram(new TH1D("isoele", "RelIso(ele);;Leptons", 100, -0.5, 9.5));

  // Taus
  mon.addHistogram(new TH1D("ntaus",         ";ntaus;Events",         10,  0, 10));
  mon.addHistogram(new TH1D("tauleadpt",     ";p_{T}^{#tau};Events",  50,  0,  100));
  mon.addHistogram(new TH1D("tauleadeta",    ";#eta^{#tau};Events",   25, -2.6,  2.6));
  mon.addHistogram(new TH1D("tauleadcharge", ";q^{#tau};Events",       5, -2,    2));
  mon.addHistogram(new TH1D("tauleaddz",     ";dz^{#tau};Events",     25,  0,    2));
  mon.addHistogram(new TH1D("tauleademfrac", ";emf^{#tau};Events",    25,  0,    5));
  TH1D *tauCutFlow   = (TH1D*)mon.addHistogram(new TH1D("tauCutFlow",   ";;#tau", 6, 0, 6));
  tauCutFlow->GetXaxis()->SetBinLabel(1, "All");
  tauCutFlow->GetXaxis()->SetBinLabel(2, "PF");
  tauCutFlow->GetXaxis()->SetBinLabel(3, "ID");
  tauCutFlow->GetXaxis()->SetBinLabel(4, "Quality");
  tauCutFlow->GetXaxis()->SetBinLabel(5, "Kin");
  tauCutFlow->GetXaxis()->SetBinLabel(6, "Iso");
  TH1D *tauID = (TH1D*)mon.addHistogram(new TH1D("tauID", ";;#tau", 5, 0, 5));
  tauID->GetXaxis()->SetBinLabel(1, "All");
  tauID->GetXaxis()->SetBinLabel(2, "Not e");
  tauID->GetXaxis()->SetBinLabel(3, "Not #mu");
  tauID->GetXaxis()->SetBinLabel(4, "Decay mode");
  tauID->GetXaxis()->SetBinLabel(5, "Medium comb iso");

  // Jets
  mon.addHistogram(new TH1D("njets", ";njets;Events", 6, 0, 6));
  mon.addHistogram(new TH1D("nbjets", ";njets;Events", 6, 0, 6));
  mon.addHistogram(new TH1D("jetleadpt", ";p_{T}^{jet};Events", 25, 0, 500));
  mon.addHistogram(new TH1D("jetleadeta", ";#eta^{jet};Events", 25, -2.6, 2.6));
  mon.addHistogram(new TH1D("jetcsv", ";csv;jets", 25, 0, 1));
  TH1D *jetCutFlow   = (TH1D*)mon.addHistogram(new TH1D("jetCutFlow",   ";;jets", 6, 0, 6));
  jetCutFlow->GetXaxis()->SetBinLabel(1, "All");
  jetCutFlow->GetXaxis()->SetBinLabel(2, "PF Loose");
  jetCutFlow->GetXaxis()->SetBinLabel(3, "Pre-Selection");
  jetCutFlow->GetXaxis()->SetBinLabel(4, "ID");
  jetCutFlow->GetXaxis()->SetBinLabel(5, "Iso");
  jetCutFlow->GetXaxis()->SetBinLabel(6, "Kin");

  // MET
  mon.addHistogram(new TH1D("MET", ";MET [GeV];Events", 25, 0, 200));

  // MT2
  mon.addHistogram(new TH1D("MT2",     ";M_{T2} [GeV];Events", 25, 0, 500));
  mon.addHistogram(new TH1D("MT2_50",  ";M_{T2(50)} [GeV];Events", 25, 0, 500));
  mon.addHistogram(new TH1D("MT2_150", ";M_{T2(150)} [GeV];Events", 25, 0, 500));

  // SVFit Mass
  mon.addHistogram(new TH1D("SVFitMass", ";M_{SVFit};Events", 50, 0, 500));
  // Invariant Mass
  mon.addHistogram(new TH1D("InvMass", ";M_{l-#tau};Events", 50, 0, 500));

  // Angles
  mon.addHistogram(new TH1D("deltaAlphaTauTau", ";#Delta#alpha_{#tau-#tau}(Lab);Events", 30, 0, TMath::Pi()));
  mon.addHistogram(new TH1D("deltaPhiTauTauMET", ";#Delta#phi_{#tau#tau-MET}(Lab);Events", 30, 0, TMath::Pi()));
  mon.addHistogram(new TH1D("deltaPhiTauTau", ";#Delta#phi_{#tau-#tau}(Lab);Events", 30, 0, TMath::Pi()));
  mon.addHistogram(new TH1D("cosThetaTauH", ";cos#theta_{#tau}(Lab);Events", 30, -1, 1));
  mon.addHistogram(new TH1D("cosThetaTauL", ";cos#theta_{l}(Lab);Events", 30, -1, 1));
  mon.addHistogram(new TH1D("deltaPhiLepMETCS", ";#Delta#phi_{l-MET}(CS);Events", 30, 0, TMath::Pi()));
  mon.addHistogram(new TH1D("cosThetaCS", ";cos#theta(CS);Events", 30, -1, 1));
  mon.addHistogram(new TH2D("metVSPTl", ";p_{T}(l);MET", 50, 0, 100, 25, 0, 200));
  mon.addHistogram(new TH2D("metVSPTtau", ";p_{T}(#tau);MET", 50, 0, 100, 25, 0, 200));



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
    // Double lepton PU distribution
    //std::vector<double> dataPileupDistributionDouble = runProcess.getParameter<std::vector<double> >("datapileup");
    // Single lepton PU distribution
    std::vector<double> dataPileupDistributionDouble = runProcess.getParameter<std::vector<double> >("datapileupSingleLep");
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
  DuplicatesChecker duplicatesChecker;
//  int nDuplicates(0);
  int step(totalEntries/50);

  // Redirect stdout and stderr to a temporary buffer, then output buffer after event loop
  std::ostream myCout(std::cout.rdbuf());
  std::stringstream buffer;
  std::streambuf *coutbuf = std::cout.rdbuf();
  std::streambuf *cerrbuf = std::cerr.rdbuf();
  std::cout.rdbuf(buffer.rdbuf());
  std::cerr.rdbuf(buffer.rdbuf());

  // Variables used in loop
  int nvtx = 0;
  std::vector<TString> chTags;
  std::vector<bool> triggerBits;
  bool triggeredOn = false;
  llvvMet met;
  double rho = 0;
  double rho25 = 0;
  double weight       = 1.;
  double weight_plus  = 1.;
  double weight_minus = 1.;
  double puWeight     = 1.;
  llvvLeptonCollection selLeptons;
  llvvJetExtCollection selJets, selJetsNoId, selBJets;
  llvvTauCollection selTaus;
  bool selected = false;
  int tauIndex = -1, leptonIndex = -1;
  bool isOS = false;
  double mass = -1;
  double invMass = -1;
  double mt2 = -1;
  double mt2_50 = -1;
  double mt2_150 = -1;
  double stauMass = 0;
  double neutralinoMass = 0;
  double deltaAlphaTauTau = 0;
  double deltaPhiTauTauMET = 0;
  double deltaPhiTauTau = 0;
  double cosThetaTauH = 0;
  double cosThetaTauL = 0;
  double deltaPhiLepMETCS = 0;
  double cosThetaCS = 0;
  double tauLeadPt = 0;
  double lepLeadPt = 0;
  double maxPtSum = 0;
  int nTauJets = 0;

  // Prepare summary tree
  if(saveSummaryTree)
  {
    // Dataset specific variables
    summaryTree->Branch("isMC", &isMC);
    summaryTree->Branch("xSecWeight", &xsecWeight);

    // Event specific variables
    summaryTree->Branch("selected", &selected);
    summaryTree->Branch("chTags", &chTags);
    summaryTree->Branch("nvtx", &nvtx);
    summaryTree->Branch("triggerBits", &triggerBits);
    summaryTree->Branch("triggeredOn", &triggeredOn);
    summaryTree->Branch("rho", &rho);
    summaryTree->Branch("rho25", &rho25);
    summaryTree->Branch("met", &met);
    summaryTree->Branch("weight", &weight);
    summaryTree->Branch("weight_plus", &weight_plus);
    summaryTree->Branch("weight_minus", &weight_minus);
    summaryTree->Branch("puWeight", &puWeight);
    summaryTree->Branch("selLeptons", &selLeptons);
    summaryTree->Branch("selJets", &selJets);
//    summaryTree->Branch("selJetsNoId", &selJetsNoId);
    summaryTree->Branch("selBJets", &selBJets);
    summaryTree->Branch("selTaus", &selTaus);
    summaryTree->Branch("isOS", &isOS);
    summaryTree->Branch("SVFitMass", &mass);
    summaryTree->Branch("InvariantMass", &invMass);
    summaryTree->Branch("MT2", &mt2);
    summaryTree->Branch("MT2_50", &mt2_50);
    summaryTree->Branch("MT2_150", &mt2_150);
    summaryTree->Branch("tauIndex", &tauIndex);
    summaryTree->Branch("leptonIndex", &leptonIndex);
    summaryTree->Branch("stauMass", &stauMass);
    summaryTree->Branch("neutralinoMass", &neutralinoMass);
    summaryTree->Branch("deltaAlphaTauTau", &deltaAlphaTauTau);
    summaryTree->Branch("deltaPhiTauTauMET", &deltaPhiTauTauMET);
    summaryTree->Branch("deltaPhiTauTau", &deltaPhiTauTau);
    summaryTree->Branch("cosThetaTauH", &cosThetaTauH);
    summaryTree->Branch("cosThetaTauL", &cosThetaTauL);
    summaryTree->Branch("deltaPhiLepMETCS", &deltaPhiLepMETCS);
    summaryTree->Branch("cosThetaCS", &cosThetaCS);
    summaryTree->Branch("tauLeadPt", &tauLeadPt);
    summaryTree->Branch("lepLeadPt", &lepLeadPt);
  }

  myCout << "       Progress Bar:0%      20%       40%       60%       80%      100%" << std::endl;
  myCout << "Scanning the ntuple:";


  // Loop on events
  for(int iev = 0; iev < totalEntries; ++iev)
  {
    if(iev%step == 0)
      myCout << "_" << std::flush;

    // Init variables
    deltaAlphaTauTau = 0;
    deltaPhiTauTauMET = 0;
    deltaPhiTauTau = 0;
    cosThetaTauH = 0;
    cosThetaTauL = 0;
    deltaPhiLepMETCS = 0;
    cosThetaCS = 0;
    nvtx = 0;
    selected = false;
    weight       = 1.;
    weight_plus  = 1.;
    weight_minus = 1.;
    puWeight     = 1.;
    chTags.clear();
    llvvLeptonCollection selLowLeptons;
    selLeptons.clear();
    selJets.clear(), selJetsNoId.clear(), selBJets.clear();
    selTaus.clear();
    tauIndex = -1, leptonIndex = -1;
    isOS = false;
    mass = -1;
    invMass = -1;
    mt2 = -1;
    stauMass = -1;
    neutralinoMass = -1;
    tauLeadPt = 0;
    lepLeadPt = 0;
    maxPtSum = 0;
    nTauJets = 0;

    // Prepare tags to fill the histograms
    chTags.push_back("all");

    // Load the event content from tree
    ev.to(iev);


    /****     Get information/collections from the event     ****/
    // Number of vertexes
    fwlite::Handle<int> nvtxHandle;
    nvtxHandle.getByLabel(ev, "llvvObjectProducersUsed", "nvtx");
    if(nvtxHandle.isValid()) nvtx = *nvtxHandle;
    //if(nvtx > 19 && doQuickFix)
    //  continue;

    // Collection of generated particles
    fwlite::Handle<llvvGenEvent> genEventHandle;
    genEventHandle.getByLabel(ev, "llvvObjectProducersUsed");
    if(!genEventHandle.isValid())
    {
      std::cout << "llvvGenEvent Object NotFound" << std::endl;
      continue;
    }
    llvvGenEvent genEv = *genEventHandle;
    /**** Remove double counting if running on exclusive samples ****/
    if(exclusiveRun && isV0JetsMC)
    {
      if(genEv.nup > 5) // Drop V+{1,2,3,4}Jets from VJets samples to avoid double counting (but keep V+0Jets) [V = W,Z]
        continue;
    }

    /****           Get LHE comments with mass info          ****/
    if(isStauStau)
    {
      fwlite::Handle<LHEEventProduct> LHEHandle;
      LHEHandle.getByLabel(ev, "source");
      if(!LHEHandle.isValid())
      {
        std::cout << "LHEEventProduct Object not Found" << std::endl;
        continue;
      }
      if(LHEHandle->comments_size() == 0)
        continue;

      for(auto comment = LHEHandle->comments_begin(); comment != LHEHandle->comments_end(); ++comment)
      {
        auto modelPos = comment->find("# model TStauStau_");
        if(modelPos != std::string::npos)
        {
          std::stringstream tmp;
          auto numPos = comment->find_first_of("1234567890", modelPos);

          tmp << comment->substr(numPos, comment->find("_", numPos)-numPos);
          tmp >> stauMass;
          tmp.clear();

          numPos = comment->find("_", numPos);
          numPos = comment->find_first_of("1234567890", numPos);
          tmp << comment->substr(numPos, comment->find("\n", numPos)-numPos);
          tmp >> neutralinoMass;

          break;
        }
      }
    }

    // Trigger Bits
    fwlite::Handle<std::vector<bool> > triggerBitsHandle;
    triggerBitsHandle.getByLabel(ev, "llvvObjectProducersUsed", "triggerBits");
    if(!triggerBitsHandle.isValid())
    {
      std::cout << "triggerBits Object NotFound" << std::endl;
      continue;
    }
    triggerBits = *triggerBitsHandle;
    /****          Sort events acording to HLT Path          ****/
//    bool singleETrigger   = triggerBits[13]; // HLT_Ele27_WP80_v*
//    bool singleMuTrigger  = triggerBits[15]; // HLT_IsoMu24_v*
    bool TauPlusE2012A    = triggerBits[18]; // HLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v*
    bool TauPlusMu2012A   = triggerBits[22]; // HLT_IsoMu18_eta2p1_LooseIsoPFTau20_v*
    //bool TauPlusE2012A    = triggerBits[19]; // HLT_Ele22_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v
    //bool TauPlusMu2012A   = triggerBits[23]; // HLT_IsoMu20_eta2p1_LooseIsoPFTau20_v*
    bool TauPlusE2012B    = triggerBits[17]; // HLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_v*
    bool TauPlusMu2012B   = triggerBits[21]; // HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v*
    bool TauPlusE2012D    = triggerBits[16]; // HLT_Ele13_eta2p1_WP90Rho_LooseIsoPFTau20_v*
    bool TauPlusMu2012D   = triggerBits[20]; // HLT_IsoMu8_eta2p1_LooseIsoPFTau20_v*
    bool TauPlusETrigger  = false;
    bool TauPlusMuTrigger = false;
    if(periodHLT == "2012A")
    {
      TauPlusETrigger  = TauPlusE2012A;
      TauPlusMuTrigger = TauPlusMu2012A;
    }
    if(periodHLT == "2012B")
    {
      TauPlusETrigger  = TauPlusE2012B;
      TauPlusMuTrigger = TauPlusMu2012B;
    }
    if(periodHLT == "2012D")
    {
      TauPlusETrigger  = TauPlusE2012D;
      TauPlusMuTrigger = TauPlusMu2012D;
    }
    triggeredOn = TauPlusETrigger || TauPlusMuTrigger;
    if(TauPlusETrigger)
      chTags.push_back("TauPlusE");
    if(TauPlusMuTrigger)
      chTags.push_back("TauPlusMu");
    /*triggeredOn = singleETrigger || singleMuTrigger;
    if(singleETrigger)
    {
      if(isSingleMuPD)  // Remove double counting between SingleE and SingleMu
        continue;
      chTags.push_back("singleE");
    }
    if(singleMuTrigger)
      chTags.push_back("singleMu"); // */

    // Rest of Gen Particles
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
      std::cout << "llvvTauCollection Object NotFound" << std::endl;
      continue;
    }
    llvvTauCollection taus = *tauCollHandle;

    // Boosted tau Collection
//    fwlite::Handle<llvvTauCollection> boostedTauCollHandle;
//    boostedTauCollHandle.getByLabel(ev, "llvvObjectProducersUsed", "boosted");
//    if(!boostedTauCollHandle.isValid())
//    {
//      std::cout << "llvvTauCollection Boosted Object NotFound" << std::endl;
//      continue;
//    }
//    llvvTauCollection boostedTaus = *boostedTauCollHandle;
    //if(boostedTaus.size() > 0)
    //  continue;
    //for(size_t i = 0; i < boostedTaus.size(); ++i)
    //  taus.push_back(boostedTaus[i]);

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
    if(selection == "IPM")
      metHandle.getByLabel(ev, "llvvObjectProducersUsed", "pfMet");
    else
      metHandle.getByLabel(ev, "llvvObjectProducersUsed", "pfMETPFlow");
    if(!metHandle.isValid())
    {
      std::cout << "llvvMet Object NotFound" << std::endl;
      continue;
    }
    met = *metHandle;

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
    rho = *rhoHandle;

    // Rho25
    fwlite::Handle<double> rho25Handle;
    rho25Handle.getByLabel(ev, "kt6PFJetsCentral", "rho");
    if(!rho25Handle.isValid())
    {
      std::cout << "rho25 Object NotFound" << std::endl;
      continue;
    }
    rho25 = *rho25Handle;


    // Pileup Weight
    if(isMC)
    {
      if(isStauStau)
      {
        int nEvents = 10000;
        double xsec = stauCrossSec(stauMass, neutralinoMass);
        xsecWeight = xsec/nEvents;
      }
      puWeight     = LumiWeights->weight(genEv.ngenITpu) * PUNorm[0];
      weight       = xsecWeight*puWeight;
      weight_plus  = PuShifters[utils::cmssw::PUUP  ]->Eval(genEv.ngenITpu) * (PUNorm[2]/PUNorm[0]);
      weight_minus = PuShifters[utils::cmssw::PUDOWN]->Eval(genEv.ngenITpu) * (PUNorm[1]/PUNorm[0]);
    }

    //goto BoostedTausQuickFix;

    // Get Leading Lepton
    for(size_t i = 0; i < leptons.size(); ++i)
    {
      int lepId = leptons[i].id;

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
        if(selection == "IPM")
        {
          if(leptons[i].pt() < 10)
            passKin = false;
          if(abs(eta) > 2.4)
            passKin = false;
        }
        else
        {
          if(leptons[i].pt() < minElPt)
            passKin = false;
          if(abs(eta) > maxElEta)
            passKin = false;
        }
        if(abs(eta) > ECALGap_MinEta && abs(eta) < ECALGap_MaxEta)  // Remove electrons that fall in ECAL Gap
          passKin = false;
      }
      else            // If Muon
      {
        if(selection == "IPM")
        {
          if(leptons[i].pt() < 15)
            passKin = false;
          if(abs(eta) > 2.4)
            passKin = false;
        }
        else
        {
          if(leptons[i].pt() < minMuPt)
            passKin = false;
          if(abs(eta) > maxMuEta)
            passKin = false;
        }
      }

      // Lepton ID
      bool passID = true;
      bool isLoose = false;
      bool isTight = false;
      Int_t idbits = leptons[i].idbits;
      if(lepId == 11)
      {
        if(selection == "IPM")
        {
          if(leptons[i].pt() < 20)
          {
            if(abs(eta) < 0.8)
            {
              if(leptons[i].electronInfoRef->mvanontrigv0 > 0.925)
                isLoose = true;
            }
            else
            {
              if(abs(eta) < 1.479)
              {
                if(leptons[i].electronInfoRef->mvanontrigv0 > 0.915)
                  isLoose = true;
              }
              else
              {
                if(leptons[i].electronInfoRef->mvanontrigv0 > 0.965)
                  isLoose = true;
              }
            }
          }
          else
          {
            if(abs(eta) < 0.8)
            {
              if(leptons[i].electronInfoRef->mvanontrigv0 > 0.905)
                isLoose = true;
              if(leptons[i].electronInfoRef->mvanontrigv0 > 0.925)
                isTight = true;
            }
            else
            {
              if(abs(eta) < 1.479)
              {
                if(leptons[i].electronInfoRef->mvanontrigv0 > 0.955)
                  isLoose = true;
                if(leptons[i].electronInfoRef->mvanontrigv0 > 0.975)
                  isTight = true;
              }
              else
              {
                if(leptons[i].electronInfoRef->mvanontrigv0 > 0.975)
                  isLoose = true;
                if(leptons[i].electronInfoRef->mvanontrigv0 > 0.985)
                  isTight = true;
              }
            }
          }

          passID = isTight;
          passID = isLoose;
          if(leptons[i].d0 > 0.045)
            passID = false;
          if(leptons[i].dZ > 0.2)
            passID = false;
          if(leptons[i].electronInfoRef->isConv)
            passID = false;
          if(leptons[i].trkLostInnerHits > 0)
            passID = false;
        }
        else
        {
          if(leptons[i].electronInfoRef->isConv)
            passID = false;

          isLoose = ((idbits >> 4) & 0x1);
          isTight = ((idbits >> 6) & 0x1);
          if(!isLoose)
            passID = false;
        }
      }
      else  // Muon
      {
        isLoose = ((idbits >> 8) & 0x1);
        isTight = ((idbits >> 10) & 0x1);
        if(selection == "IPM")
        {
          if(!isLoose)
            passID = false;
          if(leptons[i].d0 > 0.045)
            passID = false;
          if(leptons[i].dZ > 0.2)
            passID = false;
        }
        else
        {
          if(!isLoose)
            passID = false;
        }
      }

      // Lepton Isolation
      bool passIso = true;
      double relIso = utils::cmssw::relIso(leptons[i], rho);
      if(lepId == 11)
      {
        if(selection == "IPM")
        {
          if(relIso > 0.3)
            passIso = false;
        }
        else
        {
          if(relIso > 0.15)
            passIso = false;
        }
      }
      else
      {
        if(selection == "IPM")
        {
          if(relIso > 0.3)
            passIso = false;
        }
        else
        {
          if(relIso > 0.12)
            passIso = false;
        }
      }

      // Keep desired leptons
      if(passKin && passID && passIso)
        selLeptons.push_back(leptons[i]);
      if(!triggeredOn)
        continue;

      // Fill lepton control plots
      mon.fillHisto("leptonCutFlow", chTags, 0, weight);
      if(passID)
      {
        mon.fillHisto("leptonCutFlow", chTags, 1, weight);
        if(passKin)
        {
          mon.fillHisto("leptonCutFlow", chTags, 2, weight);
          if(passIso)
            mon.fillHisto("leptonCutFlow", chTags, 3, weight);
        }
      }
    }
    if(selLeptons.size() != 0)
    {
      std::sort(selLeptons.begin(), selLeptons.end(), sort_llvvObjectByPt);

      if(abs(selLeptons[0].id) == 11)
        chTags.push_back("leadingE");
      else
        chTags.push_back("leadingMu");
    }

    // Get taus
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
      if(!noTauOverlap)
      {
        if(selection == "LIP" || doLIPTauID)
        {
          for(size_t lep = 0; lep < selLeptons.size(); ++lep)
          {
            if(deltaR(tau, selLeptons[lep]) < 0.1)
            {
              passIso = false;
              break;
            }
          }
        }
      }

      bool passQual = true;
      if(selection == "LIP" || doLIPTauID)
      {
        if(abs(tau.dZ) > 0.5)
          passQual = false;
//        if(tau.emfraction >= 2.0)
//          passQual = false;
      }

      // Tau ID
      bool passID = true;
      if(!tau.passId(llvvTAUID::decayModeFinding))  passID = false;
      if(!tau.passId(llvvTAUID::againstMuonTight3)) passID = false;
      if(selection == "LIP" || doLIPTauID)
      {
        if(!tau.passId(llvvTAUID::byMediumCombinedIsolationDeltaBetaCorr3Hits)) passID = false;
        //if(!tau.passId(llvvTAUID::againstElectronMediumMVA3))                   passID = false;
        if(doTauAgainstElectronLoose)
        {
          if(!tau.passId(llvvTAUID::againstElectronLoose))                      passID = false;
        }
        else
        {
          if(!tau.passId(llvvTAUID::againstElectronMediumMVA5))                 passID = false;
        }
      }
      else
      {
        if(!tau.passId(llvvTAUID::byLooseCombinedIsolationDeltaBetaCorr3Hits))
          passID = false;
        if(!tau.passId(llvvTAUID::againstElectronLoose))
          passID = false;
      }

      // Keep selected taus
      if(passID && passKin && passIso && passQual && tau.isPF)
        selTaus.push_back(tau);
      if(!triggeredOn)
        continue;


      // Fill control histograms
      mon.fillHisto("tauCutFlow", chTags, 0, weight);
      if(tau.isPF)
      {
        mon.fillHisto("tauCutFlow", chTags, 1, weight);
        mon.fillHisto("tauID", chTags, 0, weight);
        if(tau.passId(llvvTAUID::againstElectronMediumMVA5))
        {
          mon.fillHisto("tauID", chTags, 1, weight);
          if(tau.passId(llvvTAUID::againstMuonTight3))
          {
            mon.fillHisto("tauID", chTags, 2, weight);
            if(tau.passId(llvvTAUID::decayModeFinding))
            {
              mon.fillHisto("tauID", chTags, 3, weight);
              if(tau.passId(llvvTAUID::byMediumCombinedIsolationDeltaBetaCorr3Hits))
              {
                mon.fillHisto("tauID", chTags, 4, weight);
              }
            }
          }
        }
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
    }
    if(selTaus.size() != 0)
      std::sort(selTaus.begin(), selTaus.end(), sort_llvvObjectByPt);

    // Get Jets
    nTauJets = 0;
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
      if(selection == "IPM")
      {
      }
      else
      {
        if(jets[i].pt() < 15)
          passPreSel = false;
        if(abs(jets[i].eta()) > 4.7)
          passPreSel = false;
      }

      // Cross clean with selected leptons and taus
      bool passIso = true;
      if(selection == "IPM" && !doCombined)
      {
      }
      else
      {
        double minDRlj = 9999.9;
        double minDRlg = 9999.9;
        double minDRtj = 9999.9;
        for(size_t j = 0; j < selLeptons.size(); ++j)
          minDRlj = TMath::Min(minDRlj, deltaR(jets[i], selLeptons[j]));
        for(size_t j = 0; j < selTaus.size(); ++j)
          minDRtj = TMath::Min(minDRtj, deltaR(jets[i], selTaus[j]));
        if(minDRlj < 0.4 || minDRlg < 0.4 || minDRtj < 0.4)
          passIso = false;
      }

      // Jet ID
      bool passID = true;
      Int_t idbits = jets[i].idbits;
      bool passPFLoose = (idbits & 0x01);
      if(selection == "IPM")
      {
        int fullPuId = (idbits >> 3) & 0x0f;
        bool passLooseFullPuId = ((fullPuId >> 2) & 0x01);
        passID = passLooseFullPuId;
      }
      else
      {
        int puID = (idbits >> 3) & 0x0f;
        bool passLoosePuID = ((puID >> 2) & 0x01);
//        int simplePuID = (idbits >> 7) & 0x0f;
//        bool passLooseSimplePuID = ((simplePuID >> 2) & 0x01);
        std::string jetType = ((jets[i].genj.pt() > 0)?("truejetsid"):("pujetsid"));
        passID = passLoosePuID;
      }

      // Jet Kinematics
      bool passKin = true;
      bool isTauJet = false;
      if(selection == "IPM")
      {
        if(abs(jets[i].eta()) > maxJetEta)
          passKin = false;
      }
      else
      {
        if(jets[i].pt() < minJetPt)
          passKin = false;
        if(abs(jets[i].eta()) > maxJetEta)
          passKin = false;
        isTauJet = (jets[i].pt() > minTauJetPt) && (abs(jets[i].eta()) < maxTauJetEta);
      }

      // B-jets
      bool isBJet = false;
      bool hasBtagCorr = false;
      if(selection == "IPM")
      {
        if((jets[i].csv > 0.679 && !doLIPBVeto) || (doLIPBVeto && jets[i].csv > 0.405))
        {
          isBJet = true;
          hasBtagCorr = true;
        }
      }
      else
      {
        if(jets[i].csv > 0.405)
        {
          isBJet = true;
          hasBtagCorr = true;
        }
      }

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
      if(passPreSel && passIso && passPFLoose && passID && isTauJet)
        ++nTauJets;
      if(passPreSel && passIso && passPFLoose && passID && passKin)
      {
        selJets.push_back(jets[i]);
      }
      if(passPreSel && passIso && passPFLoose && passID && passKin && isBJet)
        selBJets.push_back(jets[i]);
      if(!triggeredOn)
        continue;

      // Fill Jet control histograms
      mon.fillHisto("jetCutFlow", chTags, 0, weight);
      if(passPFLoose)
      {
        mon.fillHisto("jetCutFlow", chTags, 1, weight);
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
    }
    if(selJets.size() != 0)
      std::sort(selJets.begin(), selJets.end(), sort_llvvObjectByPt);
    if(selBJets.size() != 0)
      std::sort(selBJets.begin(), selBJets.end(), sort_llvvObjectByPt);

    // Opposite Sign requirements
    maxPtSum = 0;
    tauIndex = -1;
    leptonIndex = -1;
    for(size_t i = 0; i < selLeptons.size(); ++i)
    {
      if(selection == "IPM")
      {
      if(abs(selLeptons[i].id) == 11) // Electron
      {
        if(abs(selLeptons[i].dZ) > 0.1)
          continue;
        if((selLeptons[i].pt() < 20 && !doCombined) || (selLeptons[i].pt() < 25 && doCombined))
          continue;
        double eta = selLeptons[i].electronInfoRef->sceta;
        if(abs(eta) > 2.1)
          continue;
        double relIso = utils::cmssw::relIso(selLeptons[i], rho);
        if(relIso > 0.1)
          continue;
      }
      else // Muon
      {
        if(selLeptons[i].pt() < 20)
          continue;
        if(selLeptons[i].eta() > 2.1)
          continue;
        double relIso = utils::cmssw::relIso(selLeptons[i], rho);
        if(relIso > 0.1)
          continue;
        Int_t idbits = selLeptons[i].idbits;
        bool isTight = ((idbits >> 10) & 0x1);
        if(!isTight)
          continue;
      }
      }
      for(size_t j = 0; j < selTaus.size(); ++j)
      {
        if(selLeptons[i].id * selTaus[j].id < 0)
        {
          double PtSum = selLeptons[i].pt() + selTaus[j].pt();
          if(PtSum > maxPtSum)
          {
            maxPtSum = PtSum;
            tauIndex = j;
            leptonIndex = i;
            isOS = true;
            tauLeadPt = selTaus[tauIndex].pt();
            lepLeadPt = selLeptons[leptonIndex].pt();
          }
          else  // This allows us to jump a few loops if there are many taus/leptons
            break;
        }
      }
    }
    if(isOS && selection == "IPM") // Reject events with more leptons
    {
      if(abs(selLeptons[leptonIndex].id) == 11) //Electron tau channel
      {
        for(size_t i  = 0; i < selLeptons.size(); ++i)
        {
          if(i == (size_t)leptonIndex)
            continue;
          if(abs(selLeptons[i].id) != 11)
            continue;
          if(selLeptons[i].pt() < 20)
            continue;
          if(selLeptons[i].electronInfoRef->sceta > 2.1)
            continue;
          isOS = false;
          break;
        }
      }
      else // Muon tau channel
      {
        for(size_t i = 0; i < selLeptons.size(); ++i)
        {
          if(i == (size_t)leptonIndex)
            continue;
          isOS = false;
          break;
        }
      }
    }
    if(isOS)
    {
      if(abs(selLeptons[leptonIndex].id) == 11)
        chTags.push_back("etau");
      else
        chTags.push_back("mutau");
    }

    // Get angular variables
    if(isOS)
    {
      /**    LAB FRAME    **/
      TLorentzVector lep(selLeptons[leptonIndex].Px(), selLeptons[leptonIndex].Py(), selLeptons[leptonIndex].Pz(), selLeptons[leptonIndex].E());
      TLorentzVector tau(selTaus[tauIndex].Px(), selTaus[tauIndex].Py(), selTaus[tauIndex].Pz(), selTaus[tauIndex].E());
      TLorentzVector Tmet(met.Px(), met.Py(), met.Pz(), met.E());

      deltaAlphaTauTau = lep.Angle(tau.Vect());
      deltaPhiTauTauMET = Tmet.DeltaPhi(lep + tau);
      deltaPhiTauTau = lep.DeltaPhi(tau);

      double posSign = Tmet.CosTheta();
      cosThetaTauH = tau.CosTheta();
      cosThetaTauL = lep.CosTheta();
      if(posSign < 0)
      {
        cosThetaTauH = -cosThetaTauH;
        cosThetaTauL = -cosThetaTauL;
      }

      /**   CS FRAME   **/
      TLorentzVector tauSystem = tau + lep;
      TLorentzVector tauCS = tau, lepCS = lep, metCS = Tmet;
      TLorentzVector beam1(0, 0, 0, 0);
      TLorentzVector beam2(0, 0, 0, 0);
      double energy = sqrtS * 500.; // sqrtS / 2 * 1000
      double mom = sqrt(energy*energy + 0.938*0.938);
      if(posSign > 0)
      {
        beam1.SetPxPyPzE(0, 0,  mom, energy);
        beam2.SetPxPyPzE(0, 0, -mom, energy);
      }
      else
      {
        beam1.SetPxPyPzE(0, 0, -mom, energy);
        beam2.SetPxPyPzE(0, 0,  mom, energy);
      }

      TVector3 boost = -tauSystem.BoostVector();
      //tauSystem.Boost(boost); // By construction, this will be 0
      tauCS.Boost(boost);
      lepCS.Boost(boost);
      metCS.Boost(boost);
      beam1.Boost(boost);
      beam2.Boost(boost);

      TLorentzVector SQA = beam1 - beam2; // Spin quantization axis
      TRotation rotation;

      TVector3 newZAxis = SQA.Vect().Unit();
      TVector3 targetZaxis(0, 0, 1);
      TVector3 rotAxis = targetZaxis.Cross(newZAxis);
      double rotAngle = targetZaxis.Angle(newZAxis);
      rotation.Rotate(-rotAngle, rotAxis);

      SQA.Transform(rotation);
      tauCS.Transform(rotation);
      lepCS.Transform(rotation);
      metCS.Transform(rotation);
      beam1.Transform(rotation);
      beam2.Transform(rotation);

      cosThetaCS = lepCS.CosTheta();
      deltaPhiLepMETCS = metCS.DeltaPhi(lepCS);
    }

    // Tau-Lepton pair mass calculation
    if(isOS)
    {
      auto selLepton = selLeptons[leptonIndex];
      auto selTau    = selTaus[tauIndex];

      invMass = (selLepton+selTau).M(); // Invariant mass
      if(doSVfit)
      {
        TMatrixD covMET(2, 2);
        std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
        svFitStandalone::Vector measuredMET(met.px(), met.py(), 0);

        covMET[0][0] = met.sigx2;
        covMET[0][1] = met.sigxy;
        covMET[1][0] = met.sigxy;
        covMET[1][1] = met.sigy2;

        if(abs(selLepton.id) == 11)
          measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToElecDecay, svFitStandalone::LorentzVector(selLepton.px(), selLepton.py(), selLepton.pz(), selLepton.E())));
        else
          measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToMuDecay, svFitStandalone::LorentzVector(selLepton.px(), selLepton.py(), selLepton.pz(), selLepton.E())));
        measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToHadDecay, svFitStandalone::LorentzVector(selTau.px(), selTau.py(), selTau.pz(), selTau.E())));

        SVfitStandaloneAlgorithm SVfit_algo(measuredTauLeptons, measuredMET, covMET, 0);
        //SVfit_algo.maxObjFunctionCalls(10000) // To change the max number of iterations before minimization is terminated, default 5000
        SVfit_algo.addLogM(false); // To not use the LogM penalty, it is used by default
        //SVfit_algo.metPower(0.5); // Additional power to enhance MET likelihood, default is 1.
        //SVfit_algo.fit();
        //SVfit_algo.integrate();
        SVfit_algo.integrateVEGAS();
        //SVfit_algo.integrateMarkovChain();
        if(SVfit_algo.isValidSolution())
          mass = SVfit_algo.mass();
      }
    }

    // MT2 calculation
    if(mass >= 0 || (!doSVfit && isOS))
    {
      auto selLepton = selLeptons[leptonIndex];
      auto selTau    = selTaus[tauIndex];
      pa[0] = selLepton.M();
      pa[1] = selLepton.px();
      pa[2] = selLepton.py();
      pb[0] = selTau.M();
      pb[1] = selTau.px();
      pb[2] = selTau.py();
      pmiss[1] = met.px();
      pmiss[2] = met.py();

      mn = 0;
      mt2_event.set_momenta(pa,pb,pmiss);
      mt2_event.set_mn(mn);
      mt2 = mt2_event.get_mt2();

      mn = 50;
      mt2_event.set_momenta(pa,pb,pmiss);
      mt2_event.set_mn(mn);
      mt2_50 = mt2_event.get_mt2();

      mn = 150;
      mt2_event.set_momenta(pa,pb,pmiss);
      mt2_event.set_mn(mn);
      mt2_150 = mt2_event.get_mt2();
    }

//    BoostedTausQuickFix:
    bool stauPlot = false;
    if(stauMass == stauMtoPlot && neutralinoMass == neutralinoMtoPlot)
      stauPlot = true;

    if(!isStauStau || stauPlot)
    {
      mon.fillHisto("nvtxAll", chTags, nvtx, weight);
    }
    if(triggeredOn)
    {
      if(!isStauStau || stauPlot)
      {
        mon.fillHisto("eventflow", chTags, 0, weight);
        mon.fillHisto("nvtxHLT", chTags, nvtx, weight);
      }
      if(selLeptons.size() > 0)
      {
        if(!isStauStau || stauPlot)
        {
          mon.fillHisto("eventflow", chTags, 1, weight);
          mon.fillHisto("nvtx1l", chTags, nvtx, weight);
          mon.fillHisto("nbjets", chTags, selBJets.size(), weight);
        }
        if(selBJets.size() == 0)
        {
          if(!isStauStau || stauPlot)
          {
            mon.fillHisto("eventflow", chTags, 2, weight);
            mon.fillHisto("nvtxBveto", chTags, nvtx, weight);
          }
          if(selTaus.size() > 0)
          {
            if(!isStauStau || stauPlot)
            {
              mon.fillHisto("eventflow", chTags, 3, weight);
              mon.fillHisto("nvtx1tau", chTags, nvtx, weight);
            }
            if(isOS)
            {
              if(!isStauStau || stauPlot)
              {
                mon.fillHisto("eventflow", chTags, 4, weight);
                mon.fillHisto("nvtxOS", chTags, nvtx, weight);
              }
              if(mass >= 0 || !doSVfit)
              {
                //if(!isStauStau || stauPlot) mon.fillHisto("eventflow", chTags, 5, weight);

            selected = true;

            if(!isStauStau || stauPlot)
            {
            mon.fillHisto("nvtx", chTags, nvtx, weight);
            mon.fillHisto("nvtxraw", chTags, nvtx, weight/puWeight);
            mon.fillHisto("nup", "", genEv.nup, 1);

            mon.fillHisto("rho", chTags, rho, weight);
            mon.fillHisto("rho25", chTags, rho25, weight);

            mon.fillHisto("MET", chTags, met.pt(), weight);

            mon.fillHisto("MT2", chTags, mt2, weight);
            mon.fillHisto("MT2_50", chTags, mt2_50, weight);
            mon.fillHisto("MT2_150", chTags, mt2_150, weight);
            mon.fillHisto("SVFitMass", chTags, mass, weight);
            mon.fillHisto("InvMass", chTags, invMass, weight);

            mon.fillHisto("deltaAlphaTauTau", chTags, deltaAlphaTauTau, weight);
            mon.fillHisto("deltaPhiTauTauMET", chTags, deltaPhiTauTauMET, weight);
            mon.fillHisto("deltaPhiTauTau", chTags, deltaPhiTauTau, weight);
            mon.fillHisto("cosThetaTauH", chTags, cosThetaTauH, weight);
            mon.fillHisto("cosThetaTauL", chTags, cosThetaTauL, weight);
            mon.fillHisto("cosThetaCS", chTags, cosThetaCS, weight);
            mon.fillHisto("deltaPhiLepMETCS", chTags, deltaPhiLepMETCS, weight);

            mon.fillHisto("metVSPTl", chTags, selLeptons[leptonIndex].pt(), met.pt(), weight);
            mon.fillHisto("metVSPTtau", chTags, selTaus[tauIndex].pt(), met.pt(), weight);

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
              mon.fillHisto("tauleaddz", chTags, selTaus[0].dZ, weight);
            }

            for(auto i = selLeptons.begin(); i != selLeptons.end(); ++i)
            {
              int lepId = abs(i->id);
              if(lepId == 13)
                mon.fillHisto("isomu", chTags, utils::cmssw::relIso(*i, rho), weight);
              else if(lepId == 11)
                mon.fillHisto("isoele", chTags, utils::cmssw::relIso(*i, rho), weight);
            }

            for(auto i = selJets.begin(); i != selJets.end(); ++i)
            {
              mon.fillHisto("jetcsv", chTags, i->origcsv, weight);
            }
            }
              }
            }
          }
        }
      }
    }

    if(saveSummaryTree)
      summaryTree->Fill();

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
    summaryTree->SetDirectory(outfile);
    summaryTree->Write();
    outfile->Close();
    delete outfile;
  }

  return 0;
}

double stauCrossSec(double stauM, double neutM)
{
  //Points taken from  http://arxiv.org/abs/1204.2379
  //Mstau == 100 => 0.1
  //Mstau == 125 => 0.05
  //Mstau == 145 => 0.03
  //Mstau == 195 => 0.01
  //Mstau == 240 => 0.005
  //Mstau == 275 => 0.003
  //Mstau == 300 => 0.002
  //Mstau == 360 => 0.001
  //Mstau == 425 => 0.0005
  double a = 0.2979;
  double b = 17.626;
  double c = 67.632;
  double d = 3.463;
  return a / (1 + std::pow((stauM - b) / c, d));
}

bool electronMVAID(double mva, llvvLepton& lepton, ID_Type id)
{
  if(id == MediumID)
    id = TightID;

  double eta = lepton.electronInfoRef->sceta;
  bool pass = false;

  switch(id)
  {
  case LooseID:
    if(lepton.pt() < 20)
    {
      if(abs(eta) < 0.8)
      {
        if(mva > 0.925)
          pass = true;
      }
      else
      {
        if(abs(eta) < 1.479)
        {
          if(mva > 0.915)
            pass = true;
        }
        else
        {
          if(mva > 0.965)
            pass = true;
        }
      }
    }
    else
    {
      if(abs(eta) < 0.8)
      {
        if(mva > 0.905)
          pass = true;
      }
      else
      {
        if(abs(eta) < 1.479)
        {
          if(mva > 0.955)
            pass = true;
        }
        else
        {
          if(mva > 0.975)
            pass = true;
        }
      }
    }
    break;
  case TightID:
  default:
    if(lepton.pt() >= 20)
    {
      if(abs(eta) < 0.8)
      {
        if(mva > 0.925)
          pass = true;
      }
      else
      {
        if(abs(eta) < 1.479)
        {
          if(mva > 0.975)
            pass = true;
        }
        else
        {
          if(mva > 0.985)
            pass = true;
        }
      }
    }
    break;
  }

  return pass;
}

double tauSF(llvvTau& tau, llvvGenParticleCollection& genPartColl, TAU_E_ID eId)
{
  double scaleFactor = 1;

  // No correction necessary for tau ID
  // No correction necessary for normalization of Jet->Tau fake (if doing shape analysis, should investigate pt dependence of this)
  // No correction necessary for mu->Tau fake if using tight muon discriminator
  // Hadronic tau energy scale, no correction necessary
  // Tau charge misidentification rate, no correction necessary
  // This leaves only e->Tau fake, which must be corrected according to the anti-e discriminator used

  bool isElectronFakingTau = false;

  for(auto genPart = genPartColl.begin(); genPart != genPartColl.end(); ++genPart)
  {
    if(abs(genPart->id) == 11 && genPart->status == 3) // If the gen particle is a stable electron
    {
      if(deltaR(tau, *genPart) < 0.3)
      {
        isElectronFakingTau = true;
        break;
      }
    }
  }

  if(isElectronFakingTau)
  {
    double barrelSF = 1;
    double endcapSF = 1;

    switch(eId)
    {
    case antiELoose: // Both are 1, so no change
      break;
    case antiEMedium:
      barrelSF = 0.95;
      endcapSF = 0.75;
      break;
    case antiETight:
      barrelSF = 0.90;
      endcapSF = 0.70;
      break;
    case antiEMva:
      barrelSF = 0.85;
      endcapSF = 0.65;
      break;
    case antiEMva3Loose:
      barrelSF = 1.4; // +- 0.3
      endcapSF = 0.8; // +- 0.3
      break;
    case antiEMva3Medium:
      barrelSF = 1.6; // +- 0.3
      endcapSF = 0.8; // +- 0.3
      break;
    case antiEMva3Tight:
      barrelSF = 2.0; // +- 0.4
      endcapSF = 1.2; // +- 0.4
      break;
    case antiEMva3VTight:
      barrelSF = 2.4; // +- 0.5
      endcapSF = 1.2; // +- 0.5
      break;
    case antiEMva5Medium: // 1.6 +/- 0.3 for the barrel (abs(tauEta) < 1.5) and 1.1 +/- 0.3 for the endcap.
    default:
      barrelSF = 1.6;
      endcapSF = 1.1;
      break;
    }

    if(tau.eta() < 1.5)
    {
      scaleFactor = barrelSF;
    }
    else
    {
      scaleFactor = endcapSF;
    }
  }

  return scaleFactor;
}

double leptonIdAndIsoScaleFactor(llvvLepton& lepton)
{
  double scaleFactor = 1;

  if(abs(lepton.id) == 11) // If an electron
  {
    double isoSF = 0;
    double idSF  = 0;

    double pt = lepton.pt();
    double eta = lepton.electronInfoRef->sceta;
    if(abs(eta) < 1.479)  // Electron in barrel
    {
      if(pt < 30)
      {
        idSF  = 0.8999; // +- 0.0018
        isoSF = 0.9417; // +- 0.0019
      }
      else
      {
        idSF  = 0.9486; // +- 0.0003
        isoSF = 0.9804; // +- 0.0003
      }
    }
    else // Electron in endcap
    {
      if(pt < 30)
      {
        idSF  = 0.7945; // +- 0.0055
        isoSF = 0.9471; // +- 0.0037
      }
      else
      {
        idSF  = 0.8866; // +- 0.0001
        isoSF = 0.9900; // +- 0.0002
      }
    }

    scaleFactor = isoSF * idSF;
  }
  else // If a muon
  {
    double isoSF = 0;
    double idSF  = 0;

    double eta = lepton.eta();
    double pt  = lepton.pt();
    if(abs(eta) < 0.8) // Barrel muons
    {
      if(pt < 30)
      {
        idSF  = 0.9818; // +- 0.0005
        isoSF = 0.9494; // +- 0.0015
      }
      else
      {
        idSF  = 0.9852; // +- 0.0001
        isoSF = 0.9883; // +- 0.0003
      }
    }
    else
    {
      if(abs(eta) < 1.2) // Transition muons
      {
        if(pt < 30)
        {
          idSF  = 0.9829; // +- 0.0009
          isoSF = 0.9835; // +- 0.0020
        }
        else
        {
          idSF  = 0.9852; // +- 0.0002
          isoSF = 0.9937; // +- 0.0004
        }
      }
      else // Endcap muons
      {
        if(pt < 30)
        {
          idSF  = 0.9869; // +- 0.0007
          isoSF = 0.9923; // +- 0.0013
        }
        else
        {
          idSF  = 0.9884; // +- 0.0001
          isoSF = 0.9996; // +- 0.0005
        }
      }
    }

    scaleFactor = isoSF * idSF;
  }

  return scaleFactor;
}

double leptonTauTriggerScaleFactor(llvvLepton& lepton, llvvTau& tau)
{
  double scaleFactor = 1;
  double m0[2], sigma[2], alpha[2], n[2], norm[2]; // Index 0 - Data; Index 1 - MC
  double pt, eta;

  if(abs(lepton.id) == 11) // eTau channel
  {
    // Electron leg
    eta = lepton.electronInfoRef->sceta;
    pt  = lepton.pt();
    if(abs(eta) < 1.479) // In barrel
    {
      m0[0]    = 22.9704;
      m0[1]    = 21.7243;
      sigma[0] = 1.0258;
      sigma[1] = 0.619015;
      alpha[0] = 1.26889;
      alpha[1] = 0.739301;
      n[0]     = 1.31024;
      n[1]     = 1.34903;
      norm[0]  = 1.06409;
      norm[1]  = 1.02594;
    }
    else // In endcap
    {
      m0[0] = 21.9816;
      m0[1] = 22.1217;
      sigma[0] = 1.40993;
      sigma[1] = 1.34054;
      alpha[0] = 0.978597;
      alpha[1] = 1.8885;
      n[0] = 2.33144;
      n[1] = 1.01855;
      norm[0] = 0.937552;
      norm[1] = 4.7241;
    }

    double electronSF = 1;
    if(pt >= 20) // Do not apply for electrons with pt below threshold
    {
      double electronDataEff = efficiency(pt, m0[0], sigma[0], alpha[0], n[0], norm[0]);
      double electronMCEff   = efficiency(pt, m0[1], sigma[1], alpha[1], n[1], norm[1]);
      electronSF = electronDataEff/electronMCEff;
    }

    // Tau leg
    eta = tau.eta();
    pt  = tau.pt();
    if(abs(eta) < 1.5) // In barrel
    {
      m0[0]    = 18.538229;
      m0[1]    = 18.605055;
      sigma[0] = 0.651562;
      sigma[1] = 0.264062;
      alpha[0] = 0.324869;
      alpha[1] = 0.139561;
      n[0]     = 13.099048;
      n[1]     = 4.792849;
      norm[0]  = 0.902365;
      norm[1]  = 0.915035;
    }
    else // In endcap
    {
      m0[0]    = 18.756548;
      m0[1]    = 18.557810;
      sigma[0] = 0.230732;
      sigma[1] = 0.280908;
      alpha[0] = 0.142859;
      alpha[1] = 0.119282;
      n[0]     = 3.358497;
      n[1]     = 17.749043;
      norm[0]  = 0.851919;
      norm[1]  = 0.865756;
    }

    double tauSF = 1;
    if(pt >= 20)
    {
      double tauDataEff = efficiency(pt, m0[0], sigma[0], alpha[0], n[0], norm[0]);
      double tauMCEff   = efficiency(pt, m0[1], sigma[1], alpha[1], n[1], norm[1]);
      tauSF = tauDataEff/tauMCEff;
    }

    scaleFactor = electronSF * tauSF;
  }
  else // muTau channel
  {
    // Muon leg
    eta = lepton.eta();
    pt  = lepton.pt();
    if(eta < -1.2)
    {
      m0[0]    = 15.9977;
      m0[1]    = 16.0051;
      sigma[0] = 7.64004e-05;
      sigma[1] = 2.45144e-05;
      alpha[0] = 6.4951e-08;
      alpha[1] = 4.3335e-09;
      n[0]     = 1.57403;
      n[1]     = 1.66134;
      norm[0]  = 0.865325;
      norm[1]  = 0.87045;
    }
    else if(eta < -0.8)
    {
      m0[0]    = 17.3974;
      m0[1]    = 17.3135;
      sigma[0] = 0.804001;
      sigma[1] = 0.747636;
      alpha[0] = 1.47145;
      alpha[1] = 1.21803;
      n[0]     = 1.24295;
      n[1]     = 1.40611;
      norm[0]  = 0.928198;
      norm[1]  = 0.934983;
    }
    else if(eta < 0)
    {
      m0[0]    = 16.4307;
      m0[1]    = 15.9556;
      sigma[0] = 0.226312;
      sigma[1] = 0.0236127;
      alpha[0] = 0.265553;
      alpha[1] = 0.00589832;
      n[0]     = 1.55756;
      n[1]     = 1.75409;
      norm[0]  = 0.974462;
      norm[1]  = 0.981338;
    }
    else if(eta < 0.8)
    {
      m0[0]    = 17.313;
      m0[1]    = 15.9289;
      sigma[0] = 0.662731;
      sigma[1] = 0.0271317;
      alpha[0] = 1.3412;
      alpha[1] = 0.00448573;
      n[0]     = 1.05778;
      n[1]     = 1.92101;
      norm[0]  = 1.26624;
      norm[1]  = 0.978625;
    }
    else if(eta < 1.2)
    {
      m0[0]    = 16.9966;
      m0[1]    = 16.5678;
      sigma[0] = 0.550532;
      sigma[1] = 0.328333;
      alpha[0] = 0.807863;
      alpha[1] = 0.354533;
      n[0]     = 1.55402;
      n[1]     = 1.67085;
      norm[0]  = 0.885134;
      norm[1]  = 0.916992;
    }
    else
    {
      m0[0]    = 15.9962;
      m0[1]    = 15.997;
      sigma[0] = 0.000106195;
      sigma[1] = 7.90069e-05;
      alpha[0] = 4.95058e-08;
      alpha[1] = 4.40036e-08;
      n[0]     = 1.9991;
      n[1]     = 1.66272;
      norm[0]  = 0.851294;
      norm[1]  = 0.884502;
    }

    double muonSF = 1;
    if(pt >= 17)
    {
      double muonDataEff = efficiency(pt, m0[0], sigma[0], alpha[0], n[0], norm[0]);
      double muonMCEff   = efficiency(pt, m0[1], sigma[1], alpha[1], n[1], norm[1]);
      muonSF = muonDataEff/muonMCEff;
    }

    // Tau leg
    eta = tau.eta();
    pt  = tau.pt();
    if(abs(eta) < 1.5) // In barrel
    {
      m0[0]    = 18.604910;
      m0[1]    = 18.532997;
      sigma[0] = 0.276042;
      sigma[1] = 1.027880;
      alpha[0] = 0.137039;
      alpha[1] = 2.262950;
      n[0]     = 2.698437;
      n[1]     = 1.003322;
      norm[0]  = 0.940721;
      norm[1]  = 5.297292;
    }
    else // In endcap
    {
      m0[0]    = 18.701715;
      m0[1]    = 18.212782;
      sigma[0] = 0.216523;
      sigma[1] = 0.338119;
      alpha[0] = 0.148111;
      alpha[1] = 0.122828;
      n[0]     = 2.245081;
      n[1]     = 12.577926;
      norm[0]  = 0.895320;
      norm[1]  = 0.893975;
    }

    double tauSF = 1;
    if(pt >= 20)
    {
      double tauDataEff = efficiency(pt, m0[0], sigma[0], alpha[0], n[0], norm[0]);
      double tauMCEff   = efficiency(pt, m0[1], sigma[1], alpha[1], n[1], norm[1]);
      tauSF = tauDataEff/tauMCEff;
    }

    scaleFactor = muonSF * tauSF;
  }

  return scaleFactor;
}

// Following function from: https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2012#ETau_MuTau_trigger_turn_on_Joshu
// it parametrizes a trigger efficiency turn on curve. m is the pT of the object
double efficiency(double m, double m0, double sigma, double alpha, double n, double norm)
{
  const double sqrtPiOver2 = 1.2533141373;
  const double sqrt2 = 1.4142135624;
  double sig = fabs((double) sigma);
  double t = (m - m0)/sig;
  if(alpha < 0)
    t = -t;
  double absAlpha = fabs(alpha/sig);
  double a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
  double b = absAlpha - n/absAlpha;
  double ApproxErf;
  double arg = absAlpha / sqrt2;
  if (arg > 5.) ApproxErf = 1;
  else if (arg < -5.) ApproxErf = -1;
  else ApproxErf = TMath::Erf(arg);
  double leftArea = (1 + ApproxErf) * sqrtPiOver2;
  double rightArea = ( a * 1/TMath::Power(absAlpha - b,n-1)) / (n - 1);
  double area = leftArea + rightArea;
  if( t <= absAlpha )
  {
    arg = t / sqrt2;
    if(arg > 5.) ApproxErf = 1;
    else if (arg < -5.) ApproxErf = -1;
    else ApproxErf = TMath::Erf(arg);
    return norm * (1 + ApproxErf) * sqrtPiOver2 / area;
  }
  else
  {
    return norm * (leftArea + a * (1/TMath::Power(t-b,n-1) - 1/TMath::Power(absAlpha - b,n-1)) / (1 - n)) / area;
  }
}
