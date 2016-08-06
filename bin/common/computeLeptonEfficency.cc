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
  TString dtag=runProcess.getParameter<std::string>("dtag");

  TString suffix=runProcess.getParameter<std::string>("suffix");
  std::vector<std::string> urls=runProcess.getUntrackedParameter<std::vector<std::string> >("input");
  TString outUrl     = runProcess.getParameter<std::string>("outfile");

  //Get Good Lumi From Data
  lumiUtils::GoodLumiFilter goodLumiFilter(runProcess.getUntrackedParameter<std::vector<edm::LuminosityBlockRange> >("lumisToProcess", std::vector<edm::LuminosityBlockRange>()));

  bool filterOnlyEE(false), filterOnlyMUMU(false), filterOnlyEMU(false);
  if(!isMC)
    {
      if(dtag.Contains("DoubleEle")) filterOnlyEE=true;
      if(dtag.Contains("DoubleMu"))  filterOnlyMUMU=true;
      if(dtag.Contains("MuEG"))      filterOnlyEMU=true;
    }
  bool isSingleMuPD(!isMC && dtag.Contains("SingleMu"));  

  TString outTxtUrl= outUrl + ".txt";
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
 
  //Efficency Vs Eta
  double NewEtaMuBins[9] = {-2.4,-2.1,-1.2,-0.9,0.,0.9,1.2,2.1,2.4};
  mon.addHistogram( new TH1F("Probe_eta_Mu", "Efficency Vs #eta", 8, NewEtaMuBins));
  double NewEtaEleBins[9] = {-2.5,-2.0,-1.4,-0.8,0.,0.8,1.4,2.0,2.5};
  mon.addHistogram( new TH1F("Probe_eta_Ele", "Efficency Vs #eta", 8, NewEtaEleBins));

  // Efficency Vs Pt
  double NewPtMuBins[12] = {0.,10.,20.,25.,30.,35.,40.,50.,60.,90.,140.,300.};
  mon.addHistogram( new TH1F("Probe_pT_Mu", "Efficency Vs pT", 11, NewPtMuBins));
  double NewPtEleBins[7] = {0.,10.,20.,30.,40.,50.,100.};
  mon.addHistogram( new TH1F("Probe_pT_Ele", "Efficency Vs pT", 6, NewPtEleBins));


  //Efficency in sub-Eta region
  mon.addHistogram( new TH2F("Probe_pT_Eta_Mu",  "Probe Vs pT - Eta SubRegion", 8,  NewEtaMuBins, 11,  NewPtMuBins));
  mon.addHistogram( new TH2F("Probe_pT_Eta_Ele", "Probe Vs pT - Eta SubRegion", 8, NewEtaEleBins,  6, NewPtEleBins));

  //Lepton Kinematic
  mon.addHistogram(new TH1F("Z_mass", "Z Pick", 1000,  0,  1000));
  mon.addHistogram(new TH1F("pT",   "pT - Lepton", 1000,    0, 1000));
  mon.addHistogram(new TH1F("Eta", "Eta - Lepton",  480, -2.4,  2.4));
 
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
  //MuScleFitCorrector *muCor=getMuonCorrector(jecDir,dtag);

  //pileup weighting
  edm::LumiReWeighting* LumiWeights = NULL;
  utils::cmssw::PuShifter_t PuShifters;
  double PUNorm[] = {1,1,1};
  if(isMC){
          std::vector<double> dataPileupDistributionDouble = runProcess.getParameter< std::vector<double> >("datapileup");
          std::vector<float> dataPileupDistribution; for(unsigned int i=0;i<dataPileupDistributionDouble.size();i++){dataPileupDistribution.push_back(dataPileupDistributionDouble[i]);}
          std::vector<float> mcPileupDistribution;
          //utils::getMCPileupDistributionFromMiniAOD(urls,dataPileupDistribution.size(), mcPileupDistribution);
          utils::getMCPileupDistributionFromMiniAOD(urls,dataPileupDistribution.size(), mcPileupDistribution);
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
       if(!goodLumiFilter.isGoodLumi(ev.eventAuxiliary().run(),ev.eventAuxiliary().luminosityBlock()))continue;

       //apply trigger and require compatibilitiy of the event with the PD
       edm::TriggerResultsByName tr = ev.triggerResultsByName("HLT");
       if(!tr.isValid())return false;
       
       bool passOnlyMuon(false);
       bool passOnlyEle(false); 

       bool eeTrigger          = utils::passTriggerPatterns(tr, "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
       bool muTrigger          = utils::passTriggerPatterns(tr, "HLT_Mu34_TrkIsoVVL_v*");
       bool mumuTrigger        = utils::passTriggerPatterns(tr, "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");
       bool emuTrigger         = utils::passTriggerPatterns(tr, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*", "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*");

       if(filterOnlyEE)   { mumuTrigger=false; emuTrigger=false;  }
       if(filterOnlyMUMU) { eeTrigger=false;   emuTrigger=false;  }
       if(isSingleMuPD)   { eeTrigger=false;   emuTrigger=false;  if( muTrigger && !mumuTrigger) mumuTrigger=true; else mumuTrigger=false; }
       if(filterOnlyEMU)  { eeTrigger=false;   mumuTrigger=false; }

       float triggerPrescale(1.0),triggerThreshold(0);

       if(!(eeTrigger || muTrigger || mumuTrigger || emuTrigger ))continue;  //ONLY RUN ON THE EVENTS THAT PASS OUR TRIGGERS

       //std::cout << "The trigger selection is passed!" << std::endl;
       //##############################################   EVENT PASSED THE TRIGGER   #######################################
              
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

       //Computing of the different weights
       float weight(1.0);
       float puWeight(1.0);
       float shapeWeight(1.0);

       if(isMC){ 

          //Pileup
          int ngenITpu = 0;
          //fwlite::Handle< std::vector<PileupSummaryInfo> > puInfoH;
          //puInfoH.getByLabel(ev, "addPileupInfo");
          //for(std::vector<PileupSummaryInfo>::const_iterator it = puInfoH->begin(); it != puInfoH->end(); it++){
          //   if(it->getBunchCrossing()==0)      { ngenITpu += it->getPU_NumInteractions(); }
          //}
          ngenITpu  = vtx.size();
          puWeight  = LumiWeights->weight(ngenITpu) * PUNorm[0];

          //Shape
	  fwlite::Handle< GenEventInfoProduct > genEventInfoHandle;
          genEventInfoHandle.getByLabel(ev, "generator");
          if(genEventInfoHandle.isValid()){ if(genEventInfoHandle->weight()<0){shapeWeight*=-1;}  }

          //Global Weight
          weight = xsecWeight * puWeight * shapeWeight;

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
       
       //PRE-SELECTION based to pt value
       for(unsigned int j=0; j<leptons.size(); j++){
         if(leptons[j].pt() > 20) selLeptons.push_back(leptons[j]);
       }

       std::sort(selLeptons.begin(),   selLeptons.end(), utils::sort_CandidatesByPt);
       if(selLeptons.size() != 2) continue;
       if(abs(selLeptons[0].pdgId()) == 13 && abs(selLeptons[1].pdgId()) == 13) passOnlyMuon = true;
       if(abs(selLeptons[0].pdgId()) == 11 && abs(selLeptons[1].pdgId()) == 11) passOnlyEle = true;
     
       //std::cout << "=> Passed the selection of the events containing only muon or ele" << std::endl; 
       bool passKinMu(false),passIdMu(false),passIsoMu(false);
       bool passKinEle(false),passIdEle(false),passIsoEle(false);

       /// MUON TAG ///
       //Logical tag for Tag and Probe Muon
       bool Tag(false), ProbeMu(false), ProbeMuIso(false), ProbeMuKin(false);
       //Logical tag for Id Muons
       bool passIdTh(false), passIdMd(false), passIdLo(false), passIsoProbeMu(false);
       //Logical tag for Pass Probe Muons  
       bool PassProbeTh(false), PassProbeMd(false), PassProbeLo(false);

       /// ELE TAG ///
       //Logical tag for Tag and Probe Ele
       bool TagEle(false), ProbeEle(false), ProbeEleKin(false), ProbeEleEta(false), ProbeEleIso(false);
       //Logical tag for Id Ele
       bool passIdEleTh(false), passIdEleMd(false), passIdEleLo(false), passIsoProbeEle(false);
       //Logical tag for Pass Probe Ele
       bool PassProbeEleTh(false), PassProbeEleMd(false), PassProbeEleLo(false);
       bool ZPick(false);

       //////////////////////
       ///                ///
       /// MUON EFFICENCY ///
       ///                ///
       //////////////////////

       if(passOnlyMuon){

         int first  = rand()%2;
         int second = (first+1)%2;

         mon.fillHisto("pT",     "Mu_Led", selLeptons[0].pt(),  weight);
         mon.fillHisto("pT",  "Mu_SubLed", selLeptons[1].pt(),  weight);
         mon.fillHisto("Eta",    "Mu_Led", selLeptons[0].eta(), weight);
         mon.fillHisto("Eta", "Mu_SubLed", selLeptons[1].eta(), weight);

         double etaf = selLeptons[first].eta();
         double ptf  = selLeptons[first].pt();
         double etas = selLeptons[second].eta();
         double pts  = selLeptons[second].pt();

         //Select the Tag muon 
         passKinMu = (ptf > 20);
         passIdMu  = patUtils::passId(selLeptons[first].mu, vtx[0], patUtils::llvvMuonId::Tight, patUtils::CutVersion::ICHEP16Cut);      
         passIsoMu = patUtils::passIso(selLeptons[first].mu,  patUtils::llvvMuonIso::Tight, patUtils::CutVersion::ICHEP16Cut); 
         if(passIdMu && passKinMu && passIsoMu) Tag = true;
         if(Tag){
           mon.fillHisto("pT",  "Mu_Tag",  ptf, weight);
           mon.fillHisto("Eta", "Mu_Tag", etaf, weight);
         }         

         //Select the Probe Muon
         ProbeMuKin = (pts > 20);
         ProbeMuIso = patUtils::passIso(selLeptons[second].mu,  patUtils::llvvMuonIso::Tight, patUtils::CutVersion::ICHEP16Cut);
         if( ProbeMuKin ) ProbeMu = true;
         if(ProbeMu){
           mon.fillHisto("pT",  "Mu_NoCutProbe",  pts, weight);
           mon.fillHisto("Eta", "Mu_NoCutProbe", etas, weight);
         }
         TLorentzVector lep1(selLeptons[first].px(),selLeptons[first].py(),selLeptons[first].pz(),selLeptons[first].energy());
         TLorentzVector lep2(selLeptons[second].px(),selLeptons[second].py(),selLeptons[second].pz(),selLeptons[second].energy());
         double mass = (lep1+lep2).M();
         mon.fillHisto("Z_mass","Mu",mass,weight);
         if((mass > 70 && mass < 110) && Tag && ProbeMu) {
            mon.fillHisto("Probe_eta_Mu",     "",    etas, weight);
            mon.fillHisto("Probe_pT_Mu",      "",     pts, weight);
            mon.fillHisto("Probe_pT_Eta_Mu",  "",   etas, pts, weight);
            ZPick = true;
         }
         //Select the PassProbe
         if(ZPick){
           passIdTh = patUtils::passId(selLeptons[second].mu, vtx[0], patUtils::llvvMuonId::Tight, patUtils::CutVersion::ICHEP16Cut);
           passIdLo = patUtils::passId(selLeptons[second].mu, vtx[0], patUtils::llvvMuonId::Loose, patUtils::CutVersion::ICHEP16Cut);
           passIdMd = patUtils::passId(selLeptons[second].mu, vtx[0], patUtils::llvvMuonId::Soft, patUtils::CutVersion::ICHEP16Cut);
           passIsoProbeMu = patUtils::passIso(selLeptons[second].mu,  patUtils::llvvMuonIso::Tight, patUtils::CutVersion::ICHEP16Cut);
           if(passIdTh && passIsoProbeMu) PassProbeTh = true;
           if(passIdLo && passIsoProbeMu) PassProbeLo = true;
           if(passIdMd && passIsoProbeMu) PassProbeMd = true;
           if(PassProbeTh){
             mon.fillHisto("Probe_eta_Mu",     "Th_Pass",  etas, weight);
             mon.fillHisto("Probe_pT_Mu",      "Th_Pass",   pts, weight);
             mon.fillHisto("Probe_pT_Eta_Mu",  "Th_Pass",  etas, pts, weight);
           }
           if(PassProbeLo){
             mon.fillHisto("Probe_eta_Mu",     "Lo_Pass",  etas, weight);
             mon.fillHisto("Probe_pT_Mu",      "Lo_Pass",   pts, weight);
             mon.fillHisto("Probe_pT_Eta_Mu",  "Lo_Pass",  etas, pts, weight);
           }         
           if(PassProbeMd){
             mon.fillHisto("Probe_eta_Mu",     "Md_Pass",  etas, weight);
             mon.fillHisto("Probe_pT_Mu",      "Md_Pass",   pts, weight);
             mon.fillHisto("Probe_pT_Eta_Mu",  "Md_Pass",  etas, pts, weight);
           }
         }

      //////////////////////////
      ///                    ///
      /// ELECTRON EFFICENCY ///
      ///                    ///
      ////////////////////////// 

       } else if(passOnlyEle){
        
         int first  = rand()%2;
         int second = (first+1)%2;
         
         mon.fillHisto("pT",     "Ele_Led",   selLeptons[0].pt(), weight);
         mon.fillHisto("pT",  "Ele_SubLed",   selLeptons[1].pt(), weight); 
         mon.fillHisto("Eta",    "Ele_Led",  selLeptons[0].eta(), weight);  
         mon.fillHisto("Eta", "Ele_SubLed",  selLeptons[1].eta(), weight);         
         mon.fillHisto("Eta",    "Ele_Led_Cluster",  selLeptons[0].el.superCluster()->eta(), weight);
         mon.fillHisto("Eta", "Ele_SubLed_Cluster",  selLeptons[1].el.superCluster()->eta(), weight);

         double etaf = selLeptons[first].eta();
         double ptf  = selLeptons[first].pt();
         double etas = selLeptons[second].eta();         
         double pts  = selLeptons[second].pt();
       
         //Selection of the Tag
         passKinEle = (ptf > 20);
         passKinEle = ((abs(etas) >= 0 && abs(etas) <= 1.4442) || (abs(etas) >= 1.5660 && abs(etas) <= 2.5));
         passIdEle  = patUtils::passId(selLeptons[first].el, vtx[0], patUtils::llvvElecId::Tight, patUtils::CutVersion::ICHEP16Cut);  
         passIsoEle = patUtils::passIso(selLeptons[first].el,  patUtils::llvvElecIso::Tight, patUtils::CutVersion::ICHEP16Cut, 0.);
         if(passKinEle && passIdEle && passIsoEle) TagEle = true; 
         if(TagEle) {
           mon.fillHisto("passTag", "Ele",    1,     1.); 
           mon.fillHisto("pT",  "Ele_Tag",  ptf, weight);
           mon.fillHisto("Eta", "Ele_Tag", etaf, weight);
         }

         //Selection of the Probe
         ProbeEleKin = (pts > 20);
         ProbeEleEta = ((abs(etas) >= 0 && abs(etas) <= 1.4442) || (abs(etas) >= 1.5660 && abs(etas) <= 2.5));
         ProbeEleIso = patUtils::passIso(selLeptons[first].el,  patUtils::llvvElecIso::Tight, patUtils::CutVersion::ICHEP16Cut, 0.);
         if( ProbeEleKin && ProbeEleEta ) ProbeEle = true; 
         if(ProbeEle){
           mon.fillHisto("pT",  "Ele_NoCutProbe",  pts, weight);
           mon.fillHisto("Eta", "Ele_NoCutProbe", etas, weight);
         }
         TLorentzVector lep1(selLeptons[first].px(),selLeptons[first].py(),selLeptons[first].pz(),selLeptons[first].energy());
         TLorentzVector lep2(selLeptons[second].px(),selLeptons[second].py(),selLeptons[second].pz(),selLeptons[second].energy());
         double mass = (lep1+lep2).M();
         mon.fillHisto("Z_mass","Ele",mass,weight);
         if((mass > 70 && mass < 110) && TagEle && ProbeEle){ 
           mon.fillHisto("Probe_eta_Ele",    "",   etas, weight);
           mon.fillHisto("Probe_pT_Ele",     "",    pts, weight);
           mon.fillHisto("Probe_pT_Eta_Ele", "",   etas, pts, weight);
           ZPick = true;
         }

         //Counting the Passing Prob
         if(ZPick){
           passIdEleTh = patUtils::passId(selLeptons[second].el, vtx[0], patUtils::llvvElecId::Tight, patUtils::CutVersion::ICHEP16Cut);
           passIdEleMd = patUtils::passId(selLeptons[second].el, vtx[0], patUtils::llvvElecId::Medium, patUtils::CutVersion::ICHEP16Cut);
           passIdEleLo = patUtils::passId(selLeptons[second].el, vtx[0], patUtils::llvvElecId::Loose, patUtils::CutVersion::ICHEP16Cut);
           passIsoProbeEle = patUtils::passIso(selLeptons[first].el,  patUtils::llvvElecIso::Tight, patUtils::CutVersion::ICHEP16Cut, 0.);
           if(passIdEleTh && passIsoProbeEle) PassProbeEleTh = true;
           if(passIdEleLo && passIsoProbeEle) PassProbeEleLo = true;
           if(passIdEleMd && passIsoProbeEle) PassProbeEleMd = true;
           if(PassProbeEleTh){
             mon.fillHisto("Probe_eta_Ele",    "Th_Pass",  etas, weight);
             mon.fillHisto("Probe_pT_Ele",     "Th_Pass",   pts, weight);
             mon.fillHisto("Probe_pT_Eta_Ele", "Th_Pass",  etas, pts, weight);
           }
           if(PassProbeEleLo){ 
             mon.fillHisto("Probe_eta_Ele",    "Lo_Pass",  etas, weight);
             mon.fillHisto("Probe_pT_Ele",     "Lo_Pass",   pts, weight);
             mon.fillHisto("Probe_pT_Eta_Ele", "Lo_Pass",  etas, pts, weight);
           }
           if(PassProbeEleMd){ 
             mon.fillHisto("Probe_eta_Ele",    "Md_Pass",  etas, weight);
             mon.fillHisto("Probe_pT_Ele",     "Md_Pass",   pts, weight);
             mon.fillHisto("Probe_pT_Eta_Ele", "Md_Pass",  etas, pts, weight);
           }

         }

       }
         
  }

  //##############################################
  //########     SAVING HISTO TO FILE     ########
  //##############################################
  //save control plots to file
  printf("Results save in %s\n", outUrl.Data());
  
  //save all to the file
  TFile *ofile=TFile::Open(outUrl, "recreate");
  mon.Write();
  ofile->Close();

  if(outTxtFile)fclose(outTxtFile);

  //Now that everything is done, dump the list of lumiBlock that we processed in this job
  if(!isMC){
     goodLumiFilter.FindLumiInFiles(urls);
     goodLumiFilter.DumpToJson(((outUrl.ReplaceAll(".root",""))+".json").Data());
  }

}
