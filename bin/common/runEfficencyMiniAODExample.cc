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
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h" ////////////////////////
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "DataFormats/Math/interface/deltaR.h"

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
#include "UserCode/llvv_fwk/interface/LumiUtils.h"

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



double DeltaRMC(const pat::Electron& el, std::vector<reco::GenParticle> object)
   {
    std::vector<double> coll;
    for (unsigned int i = 0 ; i < object.size() ; i++){
        float dR = deltaR(el.superCluster()->eta(),el.phi(),object[i].eta(),object[i].phi()) ;

        coll.push_back(dR);
   }   
        std::sort(coll.begin(),coll.end());
        if(coll.size() <=  0 ) return 999.0 ;
        else return coll.at(0);
   }   

double DeltaRMCmu(const pat::Muon& mu, std::vector<reco::GenParticle> object)
   {
    std::vector<double> coll;
    for (unsigned int i = 0 ; i < object.size() ; i++){
        float dR = deltaR(mu.eta(),mu.phi(),object[i].eta(),object[i].phi()) ;

        coll.push_back(dR);
   }   
        std::sort(coll.begin(),coll.end());
        if(coll.size() <=  0 ) return 999.0 ;
        else return coll.at(0);
}

bool passFilter(pat::TriggerObjectStandAlone& Triggerobj, std::string& FilterName)
   {   
     for (unsigned h = 0; h < Triggerobj.filterLabels().size(); ++h){
       if(FilterName.compare(Triggerobj.filterLabels()[h])==0) return true;
   }   
       return false;

   }   

double DeltaRtrig(const pat::Electron& el, std::vector<pat::TriggerObjectStandAlone> object)
   {   
    std::vector<double> coll;
    for (unsigned int i = 0 ; i < object.size() ; i++){
        float dR = deltaR(el.superCluster()->eta(),el.phi(),object[i].eta(),object[i].phi()) ;
        coll.push_back(dR);
   }
        std::sort(coll.begin(),coll.end());
        if(coll.size() <=  0 ) return 999.0 ;
        else return coll.at(0);
   }

double DeltaRtrigMu(const pat::Muon& mu, std::vector<pat::TriggerObjectStandAlone> object)
   {   
    std::vector<double> coll;
    for (unsigned int i = 0 ; i < object.size() ; i++){
        float dR = deltaR(mu.eta(),mu.phi(),object[i].eta(),object[i].phi()) ;
        coll.push_back(dR);
   }
        std::sort(coll.begin(),coll.end());
        if(coll.size() <=  0 ) return 999.0 ;
        else return coll.at(0);
   }

bool hasZMother(const reco::GenParticle  p)
{   
    bool foundZ(false);                 
    const reco::Candidate  *part = (p.mother());
    // loop on the mother particles to check if is has a W has mother
    while ((part->numberOfMothers()>0)) {
        const reco::Candidate  *MomPart =part->mother();
        if ((fabs(MomPart->pdgId())==23)){
             foundZ = true;
            break;   
        }            
        part = MomPart;
    }                
    return foundZ;   
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
  TString dtag=runProcess.getParameter<std::string>("dtag");

  TString suffix=runProcess.getParameter<std::string>("suffix");
  std::vector<std::string> urls=runProcess.getUntrackedParameter<std::vector<std::string> >("input");

  TString outUrl = runProcess.getParameter<std::string>("outfile");

// good lumi mask
  lumiUtils::GoodLumiFilter goodLumiFilter(runProcess.getUntrackedParameter<std::vector<edm::LuminosityBlockRange> >("lumisToProcess", std::vector<edm::LuminosityBlockRange>()));

  bool filterOnlyEE(false), filterOnlyMUMU(false), filterOnlyEMU(false);
/*  if(!isMC)
    {
      if(dtag.Contains("DoubleEle")) filterOnlyEE=true;
      if(dtag.Contains("DoubleMu"))  filterOnlyMUMU=true;
      if(dtag.Contains("MuEG"))      filterOnlyEMU=true;
    }
*/
  bool isSingleMuPD(!isMC && dtag.Contains("SingleMu"));  
  bool isSingleElePD(!isMC && dtag.Contains("SingleEle"));  
  bool isV0JetsMC(isMC && (dtag.Contains("DYJetsToLL"))); // #FIXME should be reactivated as soon as we have exclusive jet samples
  bool isWGmc(isMC && dtag.Contains("WG"));
  bool isZGmc(isMC && dtag.Contains("ZG"));
  bool isMC_GG  = isMC && ( string(dtag.Data()).find("GG" )  != string::npos);
  bool isMC_VBF = isMC && ( string(dtag.Data()).find("VBF")  != string::npos);
  bool isMC_125OnShell = isMC && (mctruthmode==521);
  if(isMC_125OnShell) mctruthmode=125;
  bool isMC_ZZ  = isMC && ( string(dtag.Data()).find("TeV_ZZ")  != string::npos);
  bool isMC_WZ  = isMC && ( string(dtag.Data()).find("TeV_WZ")  != string::npos);

//  bool isData_DoubleEle = string(dtag.Data()).find("Data13TeV_DoubleElectron");


  //##############################################
  //########    INITIATING HISTOGRAMS     ########
  //##############################################
  SmartSelectionMonitor mon;
  
  int evcounter = 0;
  int counter1 = 0;
  int counter2 = 0; 
  //Efficency Vs Eta
  double NewEtaBins[5] = {-2.4,-1.4,0.,1.4,2.4};
  mon.addHistogram( new TH1F("Probe_eta", "Efficency Vs #eta", 480,-2.4,2.4));

  // Efficency Vs Pt
  double NewPtBins[7] = {0.,10.,20.,30.,40.,50.,100.};
  mon.addHistogram( new TH1F("Probe_pT", "Efficency Vs pT", 1000,0,1000));

  //Efficency in sub-Eta region
  mon.addHistogram( new TH2F("Probe_pT_Eta", "Probe Vs pT - Eta SubRegion", 6, NewPtBins, 4, NewEtaBins));
  
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
//  TString jecDir = runProcess.getParameter<std::string>("jecDir");
//  gSystem->ExpandPathName(jecDir);
//  FactorizedJetCorrector *jesCor        = utils::cmssw::getJetCorrector(jecDir,isMC);
//  JetCorrectionUncertainty *totalJESUnc = new JetCorrectionUncertainty((jecDir+"/MC_Uncertainty_AK5PFchs.txt").Data());
  
  //muon energy scale and uncertainties
//  MuScleFitCorrector *muCor=getMuonCorrector(jecDir,dtag);

  //lepton efficiencies
//  LeptonEfficiencySF lepEff;

  //pileup weighting
  edm::LumiReWeighting* LumiWeights = NULL;
  utils::cmssw::PuShifter_t PuShifters;
  double PUNorm[] = {1,1,1};
  if(isMC){
          std::vector<double> dataPileupDistributionDouble = runProcess.getParameter< std::vector<double> >("datapileup");
          std::vector<float> dataPileupDistribution; for(unsigned int i=0;i<dataPileupDistributionDouble.size();i++){dataPileupDistribution.push_back(dataPileupDistributionDouble[i]);}
          std::vector<float> mcPileupDistribution;
          // Temporary hack for nvtx-based pileup           
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
//    for(unsigned int f=0;f<urls.size();f++){
//     TFile* file = TFile::Open(urls[f].c_str() );
//     fwlite::Event ev(file);
//     printf("Scanning the ntuple %2i/%2i :", (int)f+1, (int)urls.size());
     int iev=0;
     int treeStep(ev.size()/50);

     for( size_t iev=0; iev<totalEntries; iev++){
      if(iev%treeStep==0){printf(".");fflush(stdout);}

       //##############################################   EVENT LOOP STARTS   ##############################################
       ev.to(iev); //load the event content from the EDM file
       //if(!isMC && duplicatesChecker.isDuplicate( ev.run, ev.lumi, ev.event) ) { nDuplicates++; continue; }

	++evcounter;
       //apply trigger and require compatibilitiy of the event with the PD
       edm::TriggerResultsByName tr = ev.triggerResultsByName("HLT");
       if(!tr.isValid())return false;
      
       //std::cout << "Definition of the Booleian variable" << std::endl;
 
       bool passOnlyMuon(false);
       bool passOnlyEle(false); 

              
       //load all the objects we will need to access
       reco::VertexCollection vtx;
       fwlite::Handle< reco::VertexCollection > vtxHandle; 
       vtxHandle.getByLabel(ev, "offlineSlimmedPrimaryVertices");
       if(vtxHandle.isValid()){ vtx = *vtxHandle;}

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

       fwlite::Handle<edm::TriggerResults> triggerBits;
       triggerBits.getByLabel(ev,"TriggerResults","","HLT");
       const edm::TriggerNames &names = ev.triggerNames(*triggerBits);

       fwlite::Handle< pat::TriggerObjectStandAloneCollection > triggerObjects;
       triggerObjects.getByLabel(ev, "selectedPatTrigger");       

       std::vector<pat::TriggerObjectStandAlone> passTagObj;
       std::vector<pat::TriggerObjectStandAlone> passTagObjMuon;

       


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
       
       //PRE-SELECTION based to pt value
       for(unsigned int j=0; j<leptons.size(); j++){
         if(leptons[j].pt() > 8) selLeptons.push_back(leptons[j]);
       }

       std::sort(selLeptons.begin(),   selLeptons.end(), utils::sort_CandidatesByPt);
       if(selLeptons.size() != 2) continue;

	 // opposite charge condition
	 

       if(abs(selLeptons[0].pdgId()) == 13 && abs(selLeptons[1].pdgId()) == 13) passOnlyMuon = true;
       if(abs(selLeptons[0].pdgId()) == 11 && abs(selLeptons[1].pdgId()) == 11) passOnlyEle = true;
     
       //std::cout << "=> Passed the selection of the events containing only muon or ele" << std::endl; 
       bool passKinMu(false),passIdMu(false),passIsoMu(false);
       bool passKinEle(false),passIdEle(false),passIsoEle(false);

       /// MUON TAG ///
       //Logical tag for Tag and Probe Muon
       bool Tag(false), Probe(false),ProbeMuKin(false),ProbeId(false); 
       //Logical tag for Id Muons
       bool passIdTh(false), passIdSf(false), passIdLo(false), passIdMd(false),passIsoProbeMu(false);
       //Logical tag for Pass Probe Muons  
       bool PassProbeTh(false), PassProbeSf(false), PassProbeLo(false),PassProbeMd(false);

       /// ELE TAG ///
       //Logical tag for Tag and Probe Ele
       bool TagEle(false), ProbeEle(false), ProbeEleKin(false), passLegforTag(false);
       bool TagMuon(false), passLegforTagMuon(false);
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
//       if(false)

         int first  = rand()%2;
         int second = (first+1)%2;

       if (! (selLeptons[first].charge()*selLeptons[second].charge() < 0) ) continue;
         mon.fillHisto("pT",     "Mu_Led", selLeptons[0].pt(),  weight);
         mon.fillHisto("pT",  "Mu_SubLed", selLeptons[1].pt(),  weight);
         mon.fillHisto("Eta",    "Mu_Led", selLeptons[0].eta(), weight);
         mon.fillHisto("Eta", "Mu_SubLed", selLeptons[1].eta(), weight);

         double etaf = selLeptons[first].eta();
         double ptf  = selLeptons[first].pt();
         double etas = selLeptons[second].eta();
         double pts  = selLeptons[second].pt();

    unsigned int eleIndex = 0;
        
        // if in MC then do the matching with MC particles
        bool  isFromZMu1=false;
        bool  isFromZMu2=false;
        if (isMC){


       fwlite::Handle<reco::GenParticleCollection> genParticles;
       genParticles.getByLabel(ev, "prunedGenParticles");
       int  theNbOfGenParticles = genParticles->size();       
            int theGenPartRef1 = -1;
            float minDr1 = 1000;
            int iteMinDr1=-1;
            int theGenPartRef2 = -1;
            float minDr2 = 1000;
            int iteMinDr2=-1;
            for (int iteGen = 0 ; iteGen < theNbOfGenParticles ; iteGen++){
                const reco::GenParticle & genMuon = (*genParticles)[iteGen];
                if (fabs(genMuon.pdgId())!=13) continue;
                float deltaR1 = sqrt(pow((etaf-genMuon.eta()),2)+ pow(acos(cos(selLeptons[first].phi()-genMuon.phi())),2)) ;
                float deltaR2 = sqrt(pow((etas-genMuon.eta()),2)+ pow(acos(cos(selLeptons[second].phi()-genMuon.phi())),2)) ;
                if (deltaR1<minDr1){
                    minDr1 = deltaR1;
                    iteMinDr1 = iteGen;
                }
                if (deltaR2<minDr2){
                    minDr2 = deltaR2;
                    iteMinDr2 = iteGen;
                }
            }
            if (minDr1<0.1) theGenPartRef1 = iteMinDr1;
            if (minDr2<0.1) theGenPartRef2 = iteMinDr2;
            
 
            if (theGenPartRef1>=0){
                const reco::GenParticle & genMuon1 = (*genParticles)[theGenPartRef1];
		isFromZMu1 = hasZMother(genMuon1);
				}

            if (theGenPartRef2>=0){
                const reco::GenParticle & genMuon2 = (*genParticles)[theGenPartRef2];
		isFromZMu2 = hasZMother(genMuon2);
				}

		}
	   
//            if ((isFromZMu1 && isFromZMu2) != true) continue;

/*         for(reco::GenParticle Genobj : *pruned){
	 	genMuons.push_back(Genobj);
                }
	 bool passMCmatching(false);
   //      std::cout<<"I found MC sample, So I'm In!"<<endl;
         double dRMC1 = DeltaRMCmu(selLeptons[first].mu,genMuons);
         double dRMC2 = DeltaRMCmu(selLeptons[second].mu,genMuons);
         if(dRMC1 < 0.1 && dRMC2 < 0.1) passMCmatching = true;
         if (!passMCmatching){std::cout<<"didn't pass MC matching. :("<<endl; continue;}
*/

             for ( pat::TriggerObjectStandAlone obj: *triggerObjects ) { 
             obj.unpackPathNames(names);
            // std::string passTagFilterMuon("hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09");
             std::string passTagFilterMuon("hltL3fL1sMu16Eta2p1L1f0Tkf20QL3trkIsoFiltered0p09");
             passLegforTagMuon   = passFilter(obj,passTagFilterMuon);
            if  (passLegforTagMuon){
             passTagObjMuon.push_back(obj);
                }
                }

         bool   passTagdRMuonCut = false;
         double dRTagMuon =DeltaRtrigMu(selLeptons[first].mu,passTagObjMuon);
         if (dRTagMuon < 0.1) passTagdRMuonCut = true;
 
	//Select the Tag muon 
         passKinMu = (ptf > 20 && abs(etaf) < 2.1);
         passIdMu  = patUtils::passId(selLeptons[first].mu, vtx[0], patUtils::llvvMuonId::StdTight, patUtils::CutVersion::ICHEP16Cut);      
         passIsoMu = patUtils::passIso(selLeptons[first].mu,  patUtils::llvvMuonIso::Tight, patUtils::CutVersion::ICHEP16Cut); 
         if(passIdMu && passKinMu && passIsoMu && passTagdRMuonCut) Tag = true;
         if(!Tag) continue;
           mon.fillHisto("pT",  "Mu_Tag",  ptf, weight);
           mon.fillHisto("Eta", "Mu_Tag", etaf, weight);
                  

         //Select all Probe Muon
         ProbeMuKin = (pts > 20);
//         ProbeId  = patUtils::passId(selLeptons[second].mu, vtx[0], patUtils::llvvMuonId::StdLoose);
         if( (/*ProbeId &&*/ ProbeMuKin) ){

         TLorentzVector lep1(selLeptons[first].px(),selLeptons[first].py(),selLeptons[first].pz(),selLeptons[first].energy());
         TLorentzVector lep2(selLeptons[second].px(),selLeptons[second].py(),selLeptons[second].pz(),selLeptons[second].energy());
         double mass = (lep1+lep2).M();
         ZPick = mass > 70 && mass < 110 ;
         mon.fillHisto("Z_mass","Mu",mass,weight);
	 if (!ZPick) continue;

	   mon.fillHisto("pT",  "Mu_Probe",  pts, weight);
           mon.fillHisto("Eta", "Mu_Probe", etas, weight);
           mon.fillHisto("Probe_pT_Eta",  "Mu_Probe",  pts, etas, weight);
         
         //Select the PassProbe
        
           passIdTh = patUtils::passId(selLeptons[second].mu, vtx[0], patUtils::llvvMuonId::StdTight, patUtils::CutVersion::ICHEP16Cut);
           passIdMd = patUtils::passId(selLeptons[second].mu, vtx[0], patUtils::llvvMuonId::StdMedium, patUtils::CutVersion::ICHEP16Cut);
           passIdLo = patUtils::passId(selLeptons[second].mu, vtx[0], patUtils::llvvMuonId::StdLoose, patUtils::CutVersion::ICHEP16Cut);
           passIsoProbeMu = patUtils::passIso(selLeptons[second].mu,  patUtils::llvvMuonIso::Tight, patUtils::CutVersion::ICHEP16Cut);
           if(passIdTh && passIsoProbeMu) PassProbeTh = true;
           if(passIdLo && passIsoProbeMu) PassProbeLo = true;
           if(passIdMd && passIsoProbeMu) PassProbeMd = true;
           if(PassProbeTh){

             mon.fillHisto("Probe_eta",     "Mu_Th_Pass", etas, weight);
             mon.fillHisto("Probe_pT",      "Mu_Th_Pass",  pts, weight);
             mon.fillHisto("Probe_pT_Eta",  "Mu_Th_Pass",  pts, etas, weight);
           }
           if(PassProbeMd){
             mon.fillHisto("Probe_eta",     "Mu_Md_Pass", etas, weight);
             mon.fillHisto("Probe_pT",      "Mu_Md_Pass",  pts, weight);
             mon.fillHisto("Probe_pT_Eta",  "Mu_Md_Pass",  pts, etas, weight);
            }

           if(PassProbeLo){
             mon.fillHisto("Probe_eta",     "Mu_Lo_Pass", etas, weight);
             mon.fillHisto("Probe_pT",      "Mu_Lo_Pass",  pts, weight);
             mon.fillHisto("Probe_pT_Eta",  "Mu_Lo_Pass",  pts, etas, weight);
           }         
	  }
	}
      //////////////////////////
      ///                    ///
      /// ELECTRON EFFICENCY ///
      ///                    ///
      ////////////////////////// 

       else if(passOnlyEle){
//        else if(false){
        
         int first  = rand()%2;
         int second = (first+1)%2;
         
       if (! (selLeptons[first].charge()*selLeptons[second].charge() < 0) ) continue;
         mon.fillHisto("pT",     "Ele_Led",   selLeptons[0].pt(), weight);
         mon.fillHisto("pT",  "Ele_SubLed",   selLeptons[1].pt(), weight); 
         mon.fillHisto("Eta",    "Ele_Led",  selLeptons[0].eta(), weight);  
         mon.fillHisto("Eta", "Ele_SubLed",  selLeptons[1].eta(), weight);         
         mon.fillHisto("Eta",    "Ele_LedSC",  selLeptons[0].el.superCluster()->eta(), weight);
         mon.fillHisto("Eta", "Ele_SubLedSC",  selLeptons[1].el.superCluster()->eta(), weight);

         double etaf = selLeptons[first].eta();
         double ptf  = selLeptons[first].pt();
         double etas = selLeptons[second].eta();         
         double pts  = selLeptons[second].pt();
         double etaSCf = selLeptons[first].superCluster()->eta();
         double etaSCs = selLeptons[second].superCluster()->eta();
      
        bool  isFromZEle1=false;
        bool  isFromZEle2=false;
        if (isMC){
       fwlite::Handle<reco::GenParticleCollection> genParticles;
       genParticles.getByLabel(ev, "prunedGenParticles");
       int  theNbOfGenParticles = genParticles->size();       
            int theGenPartRef1 = -1; 
            float minDr1 = 1000;
            int iteMinDr1=-1;
            int theGenPartRef2 = -1; 
            float minDr2 = 1000;
            int iteMinDr2=-1;
            for (int iteGen = 0 ; iteGen < theNbOfGenParticles ; iteGen++){
                const reco::GenParticle & genElectron = (*genParticles)[iteGen];
                if (fabs(genElectron.pdgId())!=11) continue;
                float deltaR1 = sqrt(pow((etaSCf-genElectron.eta()),2)+ pow(acos(cos(selLeptons[first].phi()-genElectron.phi())),2));
                float deltaR2 = sqrt(pow((etaSCs-genElectron.eta()),2)+ pow(acos(cos(selLeptons[second].phi()-genElectron.phi())),2));
                if (deltaR1<minDr1){
                    minDr1 = deltaR1;
                    iteMinDr1 = iteGen;
                }
                if (deltaR2<minDr2){
                    minDr2 = deltaR2;
                    iteMinDr2 = iteGen;
                }
            }
            if (minDr1<0.1) theGenPartRef1 = iteMinDr1;
            if (minDr2<0.1) theGenPartRef2 = iteMinDr2;
    
 
            if (theGenPartRef1>=0){
                const reco::GenParticle & genElectron1 = (*genParticles)[theGenPartRef1];
                isFromZEle1 = hasZMother(genElectron1);
                                }

            if (theGenPartRef2>=0){
                const reco::GenParticle & genElectron2 = (*genParticles)[theGenPartRef2];
                isFromZEle2 = hasZMother(genElectron2);
                                }

                }
    
//            if ((isFromZEle1 && isFromZEle2) != true) continue;



 
/*         for(reco::GenParticle Genobj : *pruned){
	 	genElectrons.push_back(Genobj);
                }
	 bool passMCmatching(false);
         std::cout<<"I found MC sample, So I'm In!"<<endl;
         double dRMC1 = DeltaRMC(selLeptons[first].el,genElectrons);
         double dRMC2 = DeltaRMC(selLeptons[second].el,genElectrons);
//cout<<"---------------------------------------1-------------------------------"<<endl;
         if(dRMC1 < 0.1 && dRMC2 < 0.1) passMCmatching = true;
         if (!passMCmatching){std::cout<<"didn't pass MC matching. :("<<endl; continue;}
*/
//cout<<"---------------------------------------2-------------------------------"<<endl;
             for ( pat::TriggerObjectStandAlone obj: *triggerObjects ) { 
             obj.unpackPathNames(names);
         //    std::string passTagFilter("hltSingleEle22WP75GsfTrackIsoFilter");
         //    std::string passTagFilter("hltSingleEle22WPTightGsfTrackIsoFilter");
             std::string passTagFilter("hltEle27WPLooseGsfTrackIsoFilter");
           //  std::string passTagFilter("hltEle27WP75GsfTrackIsoFilter");

             passLegforTag   = passFilter(obj,passTagFilter);
            if  (passLegforTag){
             passTagObj.push_back(obj);
                }
		}
           
         bool   passTagdRCut = false;
         double dRTag =DeltaRtrig(selLeptons[first].el,passTagObj);
         if (dRTag < 0.1) passTagdRCut = true;

         passKinEle = (((abs(etaf) >= 0 && abs(etaf) <= 1.4442) || (abs(etaf) >= 1.5660 && abs(etaf) < 2.1)) && ptf > 30);
         passIdEle  = patUtils::passId(selLeptons[first].el, vtx[0], patUtils::llvvElecId::Tight, patUtils::CutVersion::ICHEP16Cut);  
         passIsoEle = patUtils::passIso(selLeptons[first].el,  patUtils::llvvElecIso::Tight, patUtils::CutVersion::ICHEP16Cut, 0.);
         if(passKinEle && passIdEle && passIsoEle && passTagdRCut) TagEle = true; 
         if(!TagEle)  continue;
           mon.fillHisto("passTag", "Ele",    1,     1.);
           mon.fillHisto("pT",  "Ele_Tag",  ptf, weight);
           mon.fillHisto("Eta", "Ele_Tag", etaf, weight);
         

         //Selection of the Probe
         ProbeEleKin = (((abs(etas) >= 0 && abs(etas) <= 1.4442) || (abs(etas) >= 1.5660 && abs(etas) <= 2.5)) && pts > 20);
//         ProbeEle = patUtils::passId(selLeptons[second].el,vtx[0],  patUtils::llvvElecId::Loose); 
//         if(! (ProbeEle && ProbeEleKin )) continue;
         if( (/*ProbeEle &&*/ ProbeEleKin )){
	  
         TLorentzVector lep1(selLeptons[first].px(),selLeptons[first].py(),selLeptons[first].pz(),selLeptons[first].energy());
         TLorentzVector lep2(selLeptons[second].px(),selLeptons[second].py(),selLeptons[second].pz(),selLeptons[second].energy());
         double mass = (lep1+lep2).M();

         mon.fillHisto("Z_mass","Ele",mass,weight);
         ZPick = mass > 70 && mass < 110;
			
 	 if(!ZPick)  continue;
           mon.fillHisto("pT",  "Ele_Probe",  pts, weight);
           mon.fillHisto("Eta", "Ele_Probe", etas, weight);
           mon.fillHisto("Probe_pT_Eta", "Ele_Probe",   pts, etas, weight);

	   ++counter2;
           passIdEleTh = patUtils::passId(selLeptons[second].el, vtx[0], patUtils::llvvElecId::Tight, patUtils::CutVersion::ICHEP16Cut);
           passIdEleMd = patUtils::passId(selLeptons[second].el, vtx[0], patUtils::llvvElecId::Medium, patUtils::CutVersion::ICHEP16Cut);
           passIdEleLo = patUtils::passId(selLeptons[second].el, vtx[0], patUtils::llvvElecId::Loose, patUtils::CutVersion::ICHEP16Cut);
           passIsoProbeEle = patUtils::passIso(selLeptons[first].el,  patUtils::llvvElecIso::Loose, patUtils::CutVersion::ICHEP16Cut, 0.);
           if(passIdEleTh && passIsoProbeEle) PassProbeEleTh = true;
           if(passIdEleLo && passIsoProbeEle) PassProbeEleLo = true;
           if(passIdEleMd && passIsoProbeEle) PassProbeEleMd = true;
           if(PassProbeEleTh){
             mon.fillHisto("Probe_eta",    "Ele_Th_Pass", etas, weight);
             mon.fillHisto("Probe_pT",     "Ele_Th_Pass",  pts, weight);
             mon.fillHisto("Probe_pT_Eta", "Ele_Th_Pass",  pts, etas, weight);
           }
           if(PassProbeEleLo){ 
             mon.fillHisto("Probe_eta",    "Ele_Lo_Pass", etas, weight);
             mon.fillHisto("Probe_pT",     "Ele_Lo_Pass",  pts, weight);
             mon.fillHisto("Probe_pT_Eta", "Ele_Lo_Pass",  pts, etas, weight);
           }
           if(PassProbeEleMd){ 
             mon.fillHisto("Probe_eta",    "Ele_Sf_Pass", etas, weight);
             mon.fillHisto("Probe_pT",     "Ele_Sf_Pass",  pts, weight);
             mon.fillHisto("Probe_pT_Eta", "Ele_Sf_Pass",  pts, etas, weight);
           }
	 } 

         }      

       }
         
 //   

  //##############################################
  //########     SAVING HISTO TO FILE     ########
  //##############################################
  //save control plots to file
  printf("Results save in %s\n", outUrl.Data());
  
  //save all to the file
  TFile *ofile=TFile::Open(outUrl, "recreate");
  mon.Write();
  ofile->Close();
std::cout<<"evcounter : "<<evcounter<<" counter1 : "<<counter1<<" counter2 : "<<counter2;
  //if(outTxtFile)fclose(outTxtFile);
}
