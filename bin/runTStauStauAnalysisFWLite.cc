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
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
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

// Include MT2 library (do not forget to quote):
// http://particle.physics.ucdavis.edu/hefti/projects/doku.php?id=wimpmass    ** Code from here
// http://www.hep.phy.cam.ac.uk/~lester/mt2/    ** Other libraries
#include "UserCode/llvv_fwk/interface/mt2_bisect.h"


double stauCrossSec(double stauM, double neutM);

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
  //bool runSystematics = runProcess.getParameter<bool>("runSystematics");
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
  bool doQuickFix = false;
  if(runProcess.exists("doQuickFix"))
    doQuickFix = runProcess.getParameter<bool>("doQuickFix");

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
  //bool isSingleMuPD(!isMC && turl.Contains("SingleMu"));
  bool isV0JetsMC(isMC && (turl.Contains("DYJetsToLL_50toInf") || turl.Contains("WJets")));
  bool isStauStau(isMC && turl.Contains("TStauStau"));

  TTree* summaryTree = NULL;
  if(saveSummaryTree)
    summaryTree = new TTree("Events", "Events");



  /***************************************************************************/
  /*                         Initializing Histograms                         */
  /***************************************************************************/
  SmartSelectionMonitor mon;
  TH1F *eventflow = (TH1F*)mon.addHistogram(new TH1F("eventflow", ";;Events", 8, 0, 8));
  eventflow->GetXaxis()->SetBinLabel(1, "All");
  eventflow->GetXaxis()->SetBinLabel(2, "HLT");
  eventflow->GetXaxis()->SetBinLabel(3, "> 1l");
  eventflow->GetXaxis()->SetBinLabel(4, "B-veto");
  eventflow->GetXaxis()->SetBinLabel(5, "> 1#tau");
  eventflow->GetXaxis()->SetBinLabel(6, "OS");
  eventflow->GetXaxis()->SetBinLabel(7, "Mass");
  eventflow->GetXaxis()->SetBinLabel(8, "MT");
  // ...
  // TH2D* hist = (TH2D*)mon.addHistogram(...);

  mon.addHistogram(new TH1F("nup",     ";NUP;Events", 10, 0, 10));
  //mon.addHistogram(new TH1F("nupfilt", ";NUP;Events", 10, 0, 10));

  // PU
  mon.addHistogram(new TH1F("nvtxAll",    ";Vertices;Events",       50, -0.5, 49.5));
  mon.addHistogram(new TH1F("nvtxHLT",    ";Vertices;Events",       50, -0.5, 49.5));
  mon.addHistogram(new TH1F("nvtx1l",    ";Vertices;Events",       50, -0.5, 49.5));
  mon.addHistogram(new TH1F("nvtxBveto",    ";Vertices;Events",       50, -0.5, 49.5));
  mon.addHistogram(new TH1F("nvtx1tau",    ";Vertices;Events",       50, -0.5, 49.5));
  mon.addHistogram(new TH1F("nvtxOS",    ";Vertices;Events",       50, -0.5, 49.5));

  mon.addHistogram(new TH1F("nvtx",    ";Vertices;Events",       50, -0.5, 49.5));
  mon.addHistogram(new TH1F("nvtxraw", ";Vertices;Events",       50, -0.5, 49.5));
  mon.addHistogram(new TH1F("rho",     ";#rho;Events",           25,  0,   25));
  mon.addHistogram(new TH1F("rho25",   ";#rho(#eta<2.5);Events", 25,  0,   25));


  // Leptons
  mon.addHistogram(new TH1F("nlep",       ";nlep;Events",       10,  0,   10));
  mon.addHistogram(new TH1F("leadpt",     ";p_{T}^{l};Events",  50,  0,  100));
  mon.addHistogram(new TH1F("leadeta",    ";#eta^{l};Events",   25, -2.6,  2.6));
  mon.addHistogram(new TH1F("leadcharge", ";q^{l};Events",       5, -2,    2));
  TH1F *leptonCutFlow   = (TH1F*)mon.addHistogram(new TH1F("leptonCutFlow",   ";;Leptons", 4, 0, 4));
  leptonCutFlow->GetXaxis()->SetBinLabel(1, "All");
  leptonCutFlow->GetXaxis()->SetBinLabel(2, "ID");
  leptonCutFlow->GetXaxis()->SetBinLabel(3, "Kin");
  leptonCutFlow->GetXaxis()->SetBinLabel(4, "Iso");

  // Lepton Isolation
  mon.addHistogram(new TH1F("isomu",  "RelIso(#mu);;Leptons", 100, -0.5, 9.5));
  mon.addHistogram(new TH1F("isoele", "RelIso(ele);;Leptons", 100, -0.5, 9.5));

  // Taus
  mon.addHistogram(new TH1F("ntaus",         ";ntaus;Events",         10,  0, 10));
  mon.addHistogram(new TH1F("tauleadpt",     ";p_{T}^{#tau};Events",  50,  0,  100));
  mon.addHistogram(new TH1F("tauleadeta",    ";#eta^{#tau};Events",   25, -2.6,  2.6));
  mon.addHistogram(new TH1F("tauleadcharge", ";q^{#tau};Events",       5, -2,    2));
  mon.addHistogram(new TH1F("tauleaddz",     ";dz^{#tau};Events",     25,  0,    2));
  mon.addHistogram(new TH1F("tauleademfrac", ";emf^{#tau};Events",    25,  0,    5));
  TH1F *tauCutFlow   = (TH1F*)mon.addHistogram(new TH1F("tauCutFlow",   ";;#tau", 6, 0, 6));
  tauCutFlow->GetXaxis()->SetBinLabel(1, "All");
  tauCutFlow->GetXaxis()->SetBinLabel(2, "PF");
  tauCutFlow->GetXaxis()->SetBinLabel(3, "ID");
  tauCutFlow->GetXaxis()->SetBinLabel(4, "Quality");
  tauCutFlow->GetXaxis()->SetBinLabel(5, "Kin");
  tauCutFlow->GetXaxis()->SetBinLabel(6, "Iso");
  TH1F *tauID = (TH1F*)mon.addHistogram(new TH1F("tauID", ";;#tau", 5, 0, 5));
  tauID->GetXaxis()->SetBinLabel(1, "All");
  tauID->GetXaxis()->SetBinLabel(2, "Not e");
  tauID->GetXaxis()->SetBinLabel(3, "Not #mu");
  tauID->GetXaxis()->SetBinLabel(4, "Decay mode");
  tauID->GetXaxis()->SetBinLabel(5, "Medium comb iso");

  // Jets
  mon.addHistogram(new TH1F("njets", ";njets;Events", 6, 0, 6));
  mon.addHistogram(new TH1F("nbjets", ";njets;Events", 6, 0, 6));
  mon.addHistogram(new TH1F("jetleadpt", ";p_{T}^{jet};Events", 25, 0, 500));
  mon.addHistogram(new TH1F("jetleadeta", ";#eta^{jet};Events", 25, -2.6, 2.6));
  mon.addHistogram(new TH1F("jetcsv", ";csv;jets", 25, 0, 1));
  TH1F *jetCutFlow   = (TH1F*)mon.addHistogram(new TH1F("jetCutFlow",   ";;jets", 6, 0, 6));
  jetCutFlow->GetXaxis()->SetBinLabel(1, "All");
  jetCutFlow->GetXaxis()->SetBinLabel(2, "PF Loose");
  jetCutFlow->GetXaxis()->SetBinLabel(3, "Pre-Selection");
  jetCutFlow->GetXaxis()->SetBinLabel(4, "ID");
  jetCutFlow->GetXaxis()->SetBinLabel(5, "Iso");
  jetCutFlow->GetXaxis()->SetBinLabel(6, "Kin");

  // MET
  mon.addHistogram(new TH1F("MET", ";MET [GeV];Events", 25, 0, 200));

  // MT2
  mon.addHistogram(new TH1F("MT2",     ";M_{T2} [GeV];Events", 25, 0, 500));
  mon.addHistogram(new TH1F("MT2_50",  ";M_{T2(50)} [GeV];Events", 25, 0, 500));
  mon.addHistogram(new TH1F("MT2_150", ";M_{T2(150)} [GeV];Events", 25, 0, 500));

  // SVFit Mass
  mon.addHistogram(new TH1F("SVFitMass", ";M_{SVFit};Events", 50, 0, 500));

  // Angles
  mon.addHistogram(new TH1F("deltaAlphaTauTau", ";#Delta#alpha_{#tau-#tau}(Lab);Events", 30, 0, TMath::Pi()));
  mon.addHistogram(new TH1F("deltaPhiTauTauMET", ";#Delta#phi_{#tau#tau-MET}(Lab);Events", 30, 0, TMath::Pi()));
  mon.addHistogram(new TH1F("deltaPhiTauTau", ";#Delta#phi_{#tau-#tau}(Lab);Events", 30, 0, TMath::Pi()));
  mon.addHistogram(new TH1F("cosThetaTauH", ";cos#theta_{#tau}(Lab);Events", 30, -1, 1));
  mon.addHistogram(new TH1F("cosThetaTauL", ";cos#theta_{l}(Lab);Events", 30, -1, 1));
  mon.addHistogram(new TH1F("deltaPhiLepMETCS", ";#Delta#phi_{l-MET}(CS);Events", 30, 0, TMath::Pi()));
  mon.addHistogram(new TH1F("cosThetaCS", ";cos#theta(CS);Events", 30, -1, 1));
  mon.addHistogram(new TH2F("metVSPTl", ";p_{T}(l);MET", 50, 0, 100, 25, 0, 200));
  mon.addHistogram(new TH2F("metVSPTtau", ";p_{T}(#tau);MET", 50, 0, 100, 25, 0, 200));



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
  //int nDuplicates(0);
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
    if(nvtx > 20 && doQuickFix)
      continue;

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
    //bool singleETrigger   = triggerBits[13]; // HLT_Ele27_WP80_v*
    //bool singleMuTrigger  = triggerBits[15]; // HLT_IsoMu24_v*
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
    fwlite::Handle<llvvTauCollection> boostedTauCollHandle;
    boostedTauCollHandle.getByLabel(ev, "llvvObjectProducersUsed", "boosted");
    if(!boostedTauCollHandle.isValid())
    {
      std::cout << "llvvTauCollection Boosted Object NotFound" << std::endl;
      continue;
    }
    llvvTauCollection boostedTaus = *boostedTauCollHandle;
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
      if(selection == "IPM")
      {
      }
      else
      {
        double minDRlj = 9999.9;
        double minDRlg = 9999.9;
        double minDRtj = 9999.9;
        for(size_t j = 0; j < selLeptons.size(); ++j)
          minDRlj = TMath::Min((double)minDRlj, (double)deltaR(jets[i], selLeptons[j]));
        for(size_t j = 0; j < taus.size(); ++j)
        {
          if(taus[j].pt() < minTauPt || abs(taus[j].eta()) > maxTauEta)
            continue;
          minDRtj = TMath::Min((double)minDRtj, (double)deltaR(jets[i], taus[j]));
        }
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
        //int simplePuID = (idbits >> 7) & 0x0f;
        //bool passLooseSimplePuID = ((simplePuID >> 2) & 0x01);
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
        if(jets[i].csv > 0.679)
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
      if(selection == "IPM")
      {
      }
      else
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

      bool passQual = true;
      if(selection == "IPM")
      {
      }
      else
      {
        if(abs(tau.dZ) > 0.5)
          passQual = false;
        if(tau.emfraction >= 2.0)
          passQual = false;
      }

      // Tau ID
      bool passID = true;
      if(selection == "IPM")
      {
        if(!tau.passId(llvvTAUID::decayModeFinding))
          passID = false;
        if(!tau.passId(llvvTAUID::byLooseCombinedIsolationDeltaBetaCorr3Hits))
          passID = false;
        if(!tau.passId(llvvTAUID::againstElectronLoose))
          passID = false;
        if(!tau.passId(llvvTAUID::againstMuonTight3))
          passID = false;
      }
      else
      {
        //if(!tau.passId(llvvTAUID::againstElectronMediumMVA3))                   passID = false;
        if(!tau.passId(llvvTAUID::againstElectronMediumMVA5))                   passID = false;
        if(!tau.passId(llvvTAUID::againstMuonTight3))                           passID = false;
        if(!tau.passId(llvvTAUID::decayModeFinding))                            passID = false;
        if(!tau.passId(llvvTAUID::byMediumCombinedIsolationDeltaBetaCorr3Hits)) passID = false;
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

    // Opposite Sign requirements
    maxPtSum = 0;
    tauIndex = -1;
    leptonIndex = -1;
    for(size_t i = 0; i < selLeptons.size(); ++i)
    {
      if(abs(selLeptons[i].id) == 11) // Electron
      {
        if(abs(selLeptons[i].dZ) > 0.1)
          continue;
        if(selLeptons[i].pt() < 20)
          continue;
        double eta = selLeptons[i].electronInfoRef->sceta;
        if(abs(eta) > 2.1)
          continue;
        double relIso = utils::cmssw::relIso(selLeptons[i], rho);
        if(relIso > 0.1)
          continue;
        bool isTight = false;
        if(selLeptons[i].pt() >= 20)
        {
          if(abs(eta) < 0.8)
          {
            if(selLeptons[i].electronInfoRef->mvanontrigv0 > 0.925)
              isTight = true;
          }
          else
          {
            if(abs(eta) < 1.479)
            {
              if(selLeptons[i].electronInfoRef->mvanontrigv0 > 0.975)
                isTight = true;
            }
            else
            {
              if(selLeptons[i].electronInfoRef->mvanontrigv0 > 0.985)
                isTight = true;
            }
          }
        }
        if(!isTight)
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
      TMatrixD covMET(2, 2);
      std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
      auto selLepton = selLeptons[leptonIndex];
      auto selTau    = selTaus[tauIndex];
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
      //mass = (selLepton+selTau).M(); // Temporary hold
    }

    // MT2 calculation
    if(mass >= 0)
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

    //BoostedTausQuickFix:
    bool stauPlot = false;
    if(stauMass == stauMtoPlot && neutralinoMass == neutralinoMtoPlot)
      stauPlot = true;

    if(!isStauStau || stauPlot)
    {
      mon.fillHisto("eventflow", chTags, 0, weight);
      mon.fillHisto("nvtxAll", chTags, nvtx, weight);
    }
    if(triggeredOn)
    {
      if(!isStauStau || stauPlot)
      {
        mon.fillHisto("eventflow", chTags, 1, weight);
        mon.fillHisto("nvtxHLT", chTags, nvtx, weight);
      }
      if(selLeptons.size() > 0)
      {
        if(!isStauStau || stauPlot)
        {
          mon.fillHisto("eventflow", chTags, 2, weight);
          mon.fillHisto("nvtx1l", chTags, nvtx, weight);
          mon.fillHisto("nbjets", chTags, selBJets.size(), weight);
        }
        if(selBJets.size() == 0)
        {
          if(!isStauStau || stauPlot)
          {
            mon.fillHisto("eventflow", chTags, 3, weight);
            mon.fillHisto("nvtxBveto", chTags, nvtx, weight);
          }
          if(selTaus.size() > 0)
          {
            if(!isStauStau || stauPlot)
            {
              mon.fillHisto("eventflow", chTags, 4, weight);
              mon.fillHisto("nvtx1tau", chTags, nvtx, weight);
            }
            if(isOS)
            {
              if(!isStauStau || stauPlot)
              {
                mon.fillHisto("eventflow", chTags, 5, weight);
                mon.fillHisto("nvtxOS", chTags, nvtx, weight);
              }
              if(mass >= 0)
              {
                if(!isStauStau || stauPlot) mon.fillHisto("eventflow", chTags, 6, weight);

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
  // TODO: Get cross section as a function of mass, for now use placeholder value of 1 pb
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
