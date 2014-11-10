#include <iostream>
#include <math.h>
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
#include "UserCode/llvv_fwk/interface/GammaWeightsHandler.h"
#include "UserCode/llvv_fwk/interface/BtagUncertaintyComputer.h"

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
#include <TRandom3.h>
#include <TMath.h>

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
 
  std::vector<std::string> allWeightsURL=runProcess.getParameter<std::vector<std::string> >("weightsFile");
  std::string weightsDir( allWeightsURL.size() ? allWeightsURL[0] : "");
    
  //##############################################
  //########    INITIATING HISTOGRAMS     ########
  //##############################################
  SmartSelectionMonitor mon;
  
  int evcounter = 0; 
  //Muon Efficency Vs Eta
  //std::cout << "=> Definition of the Histo for the Eff Vs Eta" << std::endl;
  mon.addHistogram( new TH1F("Probe_eta_Mu", "Probe Vs #eta - Muon", 480, -2.4, 2.4));
  mon.addHistogram( new TH1F("PassProbe_eta_Th_Mu", "PassProbe Vs #eta - Tight Selection", 480, -2.4, 2.4));
  mon.addHistogram( new TH1F("PassProbe_eta_Sf_Mu", "PassProbe Vs #eta - Soft Selection",  480, -2.4, 2.4));
  mon.addHistogram( new TH1F("PassProbe_eta_Lo_Mu", "PassProbe Vs #eta - Loose Selection", 480, -2.4, 2.4));
  //Muon Efficency Vs Pt
  //std::cout << "=> Definition of the Histo for the Eff Vs Pt" << std::endl;
  mon.addHistogram( new TH1F("Probe_pT_Mu", "Probe Vs pT - Muon", 300, 0, 300));
  mon.addHistogram( new TH1F("PassProbe_pT_Th_Mu", "PassProbe Vs pT - Tight Selection", 300, 0, 300));
  mon.addHistogram( new TH1F("PassProbe_pT_Sf_Mu", "PassProbe Vs pT - Soft Selection",  300, 0, 300));
  mon.addHistogram( new TH1F("PassProbe_pT_Lo_Mu", "PassProbe Vs pT - Loose Selection", 300, 0, 300));  

  //Muon Kinematic
  mon.addHistogram(new TH1F("pT_Led_Mu","pT - Leading Mu", 1000, 0, 1000));
  mon.addHistogram(new TH1F("pT_SubLed_Mu", "pT - SubLeading Mu", 1000, 0, 1000));

  //Electron Efficency Vs Eta
  mon.addHistogram( new TH1F("Probe_eta_Ele", "Probe Vs #eta - Electron", 480, -2.4, 2.4));
  mon.addHistogram( new TH1F("PassProbe_eta_Th_Ele", "PassProbe Vs #eta - Tight Selection",  480, -2.4, 2.4));
  mon.addHistogram( new TH1F("PassProbe_eta_Md_Ele", "PassProbe Vs #eta - Medium Selection", 480, -2.4, 2.4));
  mon.addHistogram( new TH1F("PassProbe_eta_Lo_Ele", "PassProbe Vs #eta - Loose Selection",  480, -2.4, 2.4));
  //Electron Efficency Vs Pt
  mon.addHistogram( new TH1F("Probe_pT_Ele", "Probe Vs pT - Electron", 300, 0, 300));
  mon.addHistogram( new TH1F("PassProbe_pT_Th_Ele", "PassProbe Vs pT - Tight Selection",   300, 0, 300));
  mon.addHistogram( new TH1F("PassProbe_pT_Md_Ele", "PassProbe Vs pT - Medium  Selection", 300, 0, 300));
  mon.addHistogram( new TH1F("PassProbe_pT_Lo_Ele", "PassProbe Vs pT - Loose Selection",   300, 0, 300));

  //Electron Kinematic
  mon.addHistogram(new TH1F("dEtaIn", "dEtaIn", 500, 0, 1));
  mon.addHistogram(new TH1F("dPhiIn", "dPhiIn", 500, 0, 1));
  mon.addHistogram(new TH1F("SiEta", "SiEta", 500, 0, 1));
  mon.addHistogram(new TH1F("HE","HE", 500, 0, 1));
  mon.addHistogram(new TH1F("dx", "dx", 500, 0, 1));
  mon.addHistogram(new TH1F("dz", "dz", 500, 0, 1));
  mon.addHistogram(new TH1F("InvEInvP", "InvEInvP", 500, 0, 1));
  mon.addHistogram(new TH1F("Iso", "Iso", 500, 0, 1));
  mon.addHistogram(new TH1F("pT_Led_Ele","pT - Leading Ele", 1000, 0, 1000));
  mon.addHistogram(new TH1F("pT_SubLed_Ele", "pT - SubLeading Ele", 1000, 0, 1000));
 
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


  //##############################################
  //########           EVENT LOOP         ########
  //##############################################
  //loop on all the events
  printf("Progressing Bar     :0%%       20%%       40%%       60%%       80%%       100%%\n");
  printf("Scanning the ntuple :");
  int treeStep(totalEntries/50);
  //DuplicatesChecker duplicatesChecker;
  //int nDuplicates(0);
  //totalEntries->10
  for( size_t iev=0; iev<totalEntries; iev++){
      if(iev%treeStep==0){printf(".");fflush(stdout);}

       //##############################################   EVENT LOOP STARTS   ##############################################
       ev.to(iev); //load the event content from the EDM file
       //if(!isMC && duplicatesChecker.isDuplicate( ev.run, ev.lumi, ev.event) ) { nDuplicates++; continue; }

       //apply trigger and require compatibilitiy of the event with the PD
       edm::TriggerResultsByName tr = ev.triggerResultsByName("HLT");
       if(!tr.isValid())return false;
      
       //std::cout << "Definition of the Booleian variable" << std::endl;
 
       bool passOnlyMuon(false);
       bool passOnlyEle(false); 

       bool eeTrigger          = utils::passTriggerPatterns(tr, "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
       bool muTrigger          = utils::passTriggerPatterns(tr, "HLT_IsoMu24_eta2p1_v*");
       bool mumuTrigger        = utils::passTriggerPatterns(tr, "HLT_Mu17_Mu8_v*", "HLT_Mu17_TkMu8_v*"); 
       bool emuTrigger         = utils::passTriggerPatterns(tr, "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*", "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
       if(filterOnlyEE)   { mumuTrigger=false; emuTrigger=false;  }
       if(filterOnlyMUMU) { eeTrigger=false;   emuTrigger=false;  }
       if(isSingleMuPD)   { eeTrigger=false;   emuTrigger=false;  if( muTrigger && !mumuTrigger) mumuTrigger=true; else mumuTrigger=false; }
       if(filterOnlyEMU)  { eeTrigger=false;   mumuTrigger=false; }

       bool hasPhotonTrigger(false);
       float triggerPrescale(1.0),triggerThreshold(0);
       bool runPhotonSelection(mctruthmode==22 || mctruthmode==111);
       if(runPhotonSelection){
	  eeTrigger=false; mumuTrigger=false;

          std::string successfulPath="";
          if(     utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_v*")){ hasPhotonTrigger=true; triggerThreshold=92; }
          else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_v*")){ hasPhotonTrigger=true; triggerThreshold=77; }
          else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_v*")){ hasPhotonTrigger=true; triggerThreshold=50; }
          else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_v*")){ hasPhotonTrigger=true; triggerThreshold=36; }
          
          if(successfulPath!=""){ //get the prescale associated to it
             fwlite::Handle< pat::PackedTriggerPrescales > prescalesHandle;
             prescalesHandle.getByLabel(ev, "patTrigger");
             pat::PackedTriggerPrescales prescales = *prescalesHandle;
             const edm::TriggerResults& trResults =  prescales.triggerResults();
             prescales.setTriggerNames( ev.triggerNames(trResults) );
             triggerPrescale = prescales.getPrescaleForName(successfulPath);
          }
      }
       if(!(eeTrigger || muTrigger || mumuTrigger || emuTrigger || hasPhotonTrigger))continue;  //ONLY RUN ON THE EVENTS THAT PASS OUR TRIGGERS

       //std::cout << "The trigger selection is passed!" << std::endl;
       //##############################################   EVENT PASSED THE TRIGGER   #######################################

       evcounter++;
       TRandom3 *rndm  = new TRandom3(evcounter);
       //std::cout << "The random variable is defined" << std::endl;
       
       //load all the objects we will need to access
       reco::VertexCollection vtx;
       fwlite::Handle< reco::VertexCollection > vtxHandle; 
       vtxHandle.getByLabel(ev, "offlineSlimmedPrimaryVertices");
       if(vtxHandle.isValid()){ vtx = *vtxHandle;}

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

       ///////////////////////
       ///                 ///
       /// LEPTON ANALYSIS ///
       ///                 ///
       ///////////////////////


       //start by merging electrons and muons
       //std::cout << "=>Merging Leptons" << std::endl;
       std::vector<patUtils::GenericLepton> leptons;
       for(size_t l=0;l<electrons.size();l++){leptons.push_back(patUtils::GenericLepton(electrons[l]));}      
       for(size_t l=0;l<muons    .size();l++){leptons.push_back(patUtils::GenericLepton(muons    [l]));}      
       std::sort(leptons.begin(),   leptons.end(), utils::sort_CandidatesByPt);

       LorentzVector muDiff(0,0,0,0);
       std::vector<patUtils::GenericLepton> selLeptons, extraLeptons;
       //Request exactly two leptons
       //std::cout << "=> Cut on the number of leptons in the event" << std::endl;
       if(leptons.size() != 2) continue;
       //std::cout << "=> The dimension of the vector containing the leptons is: " << leptons.size() << std::endl;
       //std::cout << "=> PdgIds are => " << leptons[0].pdgId() << ", the second is " << leptons[1].pdgId() << std::endl;
       if(abs(leptons[0].pdgId()) == 13 && abs(leptons[1].pdgId()) == 13) passOnlyMuon = true;
       if(abs(leptons[0].pdgId()) == 11 && abs(leptons[1].pdgId()) == 11) passOnlyEle = true;
     
       //std::cout << "=> Passed the selection of the events containing only muon or ele" << std::endl; 
       bool passKinMu(false),passIdMu(false),passIsoMu(false);
       bool passKinEle(false),passIdEle(false),passIsoEle(false);

       /// MUON TAG ///
       //Logical tag for Tag and Probe Muon
       bool Tag(false), Probe(false);
       //Logical tag for Id Muons
       bool passIdTh(false), passIdSf(false), passIdLo(false);
       //Logical tag for Pass Probe Muons  
       bool PassProbeTh(false), PassProbeSf(false), PassProbeLo(false);

       /// ELE TAG ///
       //Logical tag for Tag and Probe Ele
       bool TagEle(false), ProbeEle(false);
       //Logical tag for Id Ele
       bool passIdEleTh(false), passIdEleMd(false), passIdEleLo(false);
       //Logical tag for Pass Probe Ele
       bool PassProbeEleTh(false), PassProbeEleMd(false), PassProbeEleLo(false);
       bool ZPick(false);

       //////////////////////
       ///                ///
       /// MUON EFFICENCY ///
       ///                ///
       //////////////////////

       if(passOnlyMuon){
         //std::cout << "=> Analyzing a couple of Muons" << std::endl; 
         double first  = rndm->Rndm();
         double second = rndm->Rndm();
         int select_first;
         int select_second;
         mon.fillHisto("pT_Led_Mu","",leptons[0].pt(),weight);
         mon.fillHisto("pT_SubLed_Mu","",leptons[1].pt(),weight);
         //std::cout << "=> Filling Random variables " << std::endl;
         //std::cout << " First Random Number: " << first << ", Second Number: " << second << std::endl;
         if(first >= second){
            select_first  = 0;
            select_second = 1;
         } else {
            select_first  = 1;
            select_second = 0;
         } 

         //std::cout << "First random number is greater than the second one" << std::endl;
         double etaf = leptons[select_first].eta();
         double ptf  = leptons[select_first].pt();
         double etas = leptons[select_second].eta();
         double pts  = leptons[select_second].pt();

         //std::cout << "The Eta of the First Lep is: " << etaf << ", of the Second Lep is: " << etas << std::endl;
         //std::cout << "The Pt of the First Lep is: " << ptf << ", of the Second Lep is: " << pts << std::endl;

         //Select the Tag muon
         passKinMu = (abs(etaf) < 2.4);
         passKinMu = (ptf > 20);
         passIdMu  = patUtils::passId(leptons[select_first].mu, vtx[0], patUtils::llvvMuonId::Tight);      
         passIsoMu = patUtils::passIso(leptons[select_first].mu,  patUtils::llvvMuonIso::Loose); 
         if(passIdMu && passKinMu) Tag = true;
         if(passIdMu){
           std::cout << "I Have passed the tight selection of the muon" << std::endl;
         }
         //std::cout << "=> Selected the Tag Muon" << std::endl;
         //Select the Probe Muon
         Probe = (pts > 10);
         TLorentzVector lep1(leptons[select_first].px(),leptons[select_first].py(),leptons[select_first].pz(),leptons[select_first].energy());
         TLorentzVector lep2(leptons[select_second].px(),leptons[select_second].py(),leptons[select_second].pz(),leptons[select_second].energy());
         double mass = (lep1+lep2).M();
         if((mass > 70 && mass < 110) && Tag && Probe) {
            mon.fillHisto("Probe_eta_Mu","",etas,weight);
            mon.fillHisto("Probe_pT_Mu","",pts,weight);
            ZPick = true;
         }
         //Select the PassProbe
         if(ZPick){
           //std::cout << "=> Selection of the Id of the Particle" << std::endl;
           passIdTh = patUtils::passId(leptons[select_second].mu, vtx[0], patUtils::llvvMuonId::Tight);
           passIdLo = patUtils::passId(leptons[select_second].mu, vtx[0], patUtils::llvvMuonId::Loose);
           passIdSf = patUtils::passId(leptons[select_second].mu, vtx[0], patUtils::llvvMuonId::Soft);
           if(passIdTh) PassProbeTh = true;
           if(passIdLo) PassProbeLo = true;
           if(passIdSf) PassProbeSf = true;
           if(passIdTh){
             mon.fillHisto("PassProbe_eta_Th_Mu","",etas,weight);
             mon.fillHisto("PassProbe_pT_Th_Mu","",pts,weight);
           }
           if(passIdLo){
             mon.fillHisto("PassProbe_eta_Lo_Mu","",etas,weight);
             mon.fillHisto("PassProbe_pT_Lo_Mu","",pts,weight);
           }         
           if(passIdSf){
             mon.fillHisto("PassProbe_eta_Sf_Mu","",etas,weight);
             mon.fillHisto("PassProbe_pT_Sf_Mu","",pts,weight);
           }
         }

      //////////////////////////
      ///                    ///
      /// ELECTRON EFFICENCY ///
      ///                    ///
      ////////////////////////// 

       } else if(passOnlyEle){
         double first  = rndm->Rndm();
         double second = rndm->Rndm();
         int select_first;
         int select_second;
         for(int i=0; i<leptons.size(); i++){
            float dEtaln         = fabs(leptons[i].el.deltaEtaSuperClusterTrackAtVtx());
            float dPhiln         = fabs(leptons[i].el.deltaPhiSuperClusterTrackAtVtx());
            float sigmaletaleta  = leptons[i].el.sigmaIetaIeta();
            float hem            = leptons[i].el.hadronicOverEm();
            double resol         = fabs((1/leptons[i].el.ecalEnergy())-(leptons[i].el.eSuperClusterOverP()/leptons[i].el.ecalEnergy()));
            double dxy           = fabs(leptons[i].el.gsfTrack()->dxy(vtx[0].position()));
            double dz            = fabs(leptons[i].el.gsfTrack()->dz(vtx[0].position()));
            mon.fillHisto("dEtaIn","",dEtaln,weight);
            mon.fillHisto("dPhiIn","",dPhiln,weight);
            mon.fillHisto("SiEta","",sigmaletaleta,weight);
            mon.fillHisto("HE","",hem,weight);
            mon.fillHisto("dx","",dxy,weight);
            mon.fillHisto("dz","",dz,weight);
            mon.fillHisto("InvEInvP","",resol,weight);
         }
         mon.fillHisto("pT_Led_Ele","",leptons[0].pt(),weight);
         mon.fillHisto("pT_SubLed_Ele","",leptons[1].pt(),weight);   
         if(first >= second){
           select_first  = 0;
           select_second = 1;
         } else {
           select_first  = 1;
           select_second = 0;
         }
         double etaf = leptons[select_first].el.superCluster()->eta();
         double ptf  = leptons[select_first].pt();
         double etas = leptons[select_second].el.superCluster()->eta();
         double pts  = leptons[select_second].pt();

         std::cout << "Electron Eta: " << etaf << ", " << etas << std::endl;
         std::cout << "Electron Pt: " << ptf << ", " << pts << std::endl;
         //Selection of the Tag
         passKinEle = (ptf > 20);
         passIdEle  = patUtils::passId(leptons[select_first].el, vtx[0], patUtils::llvvElecId::Loose); 
         passIsoEle = patUtils::passIso(leptons[select_first].el,  patUtils::llvvElecIso::Loose);
         if(passKinEle && passIdEle) TagEle = true; 
         if(TagEle){
            std::cout << "Pass Tag Selection, First > Second" << std::endl;
         }
         //Selection of the Probe
         ProbeEle = (pts > 10);
         if(ProbeEle){
           std::cout << "Pass Probe Selection" << std::endl;
         }
         TLorentzVector lep1(leptons[select_first].px(),leptons[select_first].py(),leptons[select_first].pz(),leptons[select_first].energy());
         TLorentzVector lep2(leptons[select_second].px(),leptons[select_second].py(),leptons[select_second].pz(),leptons[select_second].energy());
         double mass = (lep1+lep2).M();
         //std::cout << "L'invariant Mass of the two Ele is: " << mass << std::endl;
         if((mass > 70 && mass < 110) && TagEle && ProbeEle){
           std::cout << "The Invariant Mass is: " << mass << std::endl;
           mon.fillHisto("Probe_eta_Ele","",etas,weight);
           mon.fillHisto("Probe_pT_Ele","",pts,weight);
           ZPick = true;
           std::cout << "Counting the Probe" << std::endl;
       
         }

         //Counting the Passing Prob
         if(ZPick){
           std::cout << "Passed the ZPick" << std::endl;
           passIdEleTh = patUtils::passId(leptons[select_second].el, vtx[0], patUtils::llvvElecId::Tight);
           passIdEleMd = patUtils::passId(leptons[select_second].el, vtx[0], patUtils::llvvElecId::Medium);
           passIdEleLo = patUtils::passId(leptons[select_second].el, vtx[0], patUtils::llvvElecId::Loose);

           if(passIdEleTh) PassProbeEleTh = true;
           if(passIdEleLo) PassProbeEleLo = true;
           if(passIdEleMd) PassProbeEleMd = true;
           std::cout << "Counting the PassProbe" << std::endl;
           if(passIdEleTh){
             mon.fillHisto("PassProbe_eta_Th_Ele","",etas,weight);
             mon.fillHisto("PassProbe_pT_Th_Ele","",pts,weight);
           }
           if(passIdEleLo){
             mon.fillHisto("PassProbe_eta_Lo_Ele","",etas,weight);
             mon.fillHisto("PassProbe_pT_Lo_Ele","",pts,weight);
           }
           if(passIdEleMd){
             mon.fillHisto("PassProbe_eta_Md_Ele","",etas,weight);
             mon.fillHisto("PassProbe_pT_Md_Ele","",pts,weight);
           }

         }

       }
         
  }

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
