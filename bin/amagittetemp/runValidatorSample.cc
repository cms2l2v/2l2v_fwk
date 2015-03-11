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
  std::vector<std::string> allWeightsURL=runProcess.getParameter<std::vector<std::string> >("weightsFile");
  std::string weightsDir( allWeightsURL.size() ? allWeightsURL[0] : "");
    		  
  //##############################################
  //########    INITIATING HISTOGRAMS     ########
  //##############################################
  SmartSelectionMonitor mon;

     
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


  higgs::utils::EventCategory eventCategoryInst(higgs::utils::EventCategory::EXCLUSIVE2JETSVBF); //jet(0,>=1)+vbf binning


  //##############################################
  //########           EVENT LOOP         ########
  //##############################################
  //loop on all the events
  printf("Progressing Bar     :0%%       20%%       40%%       60%%       80%%       100%%\n");
  printf("Scanning the ntuple :");
  int treeStep(totalEntries/50);
  //DuplicatesChecker duplicatesChecker;
  //int nDuplicates(0);
  for( size_t iev=0; iev<totalEntries; iev++ ) { 
      if( iev%treeStep == 0 ) { printf("."); fflush(stdout); }

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

      /////////////////////////////////////////
      //                                     // 
      //  DERIVE WEIGHTS TO APPLY TO SAMPLE  //
      //                                     //
      /////////////////////////////////////////

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


      //Higgs specific weights
      float lShapeWeights[3]={ 1.0, 1.0, 1.0 };
      for( unsigned int nri=0; nri<NRparams.size(); nri++ ){ NRweights[nri] = 1.0; }
      if( isMC ){
	LorentzVector higgs(0,0,0,0);
	LorentzVector totLeptons(0,0,0,0);
	for(size_t igen=0; igen<gen.size(); igen++){
	  if(gen[igen].status()!=3) continue;
	  if(abs(gen[igen].pdgId())>=11 && abs(gen[igen].pdgId())<=16) totLeptons += gen[igen].p4();
	  if(gen[igen].pdgId()==25)                                      higgs=gen[igen].p4();
	}
	if(mctruthmode==125) {
	  higgs=totLeptons;
	  if(isMC_125OnShell && higgs.mass()>180) continue;
	  if(!isMC_125OnShell && higgs.mass()<=180) continue;
	}
	float shapeWeight(1.0);
        if((isMC_VBF || isMC_GG) && higgs.pt()>0){
	  {
	    //Line shape weights 
	    if(isMC_VBF || isMC_GG)
	      {
		std::vector<TGraph *> nominalShapeWgtGr=hLineShapeGrVec.begin()->second;
		for(size_t iwgt=0; iwgt<nominalShapeWgtGr.size(); iwgt++)
		  {
		    if(nominalShapeWgtGr[iwgt]==0) continue;
		    lShapeWeights[iwgt]=nominalShapeWgtGr[iwgt]->Eval(higgs.mass());
		  }
	      }
	    shapeWeight   = lShapeWeights[0];
	  }
	} 
	//final event weight
	weight = puWeight * shapeWeight;
      }

      ////////////////////////
      //                    //
      //  PHOTON SELECTION  //
      //                    //
      ////////////////////////

      pat::PhotonCollection selPhotons;
      if(runPhotonSelection)
        {
          for(size_t ipho=0; ipho<photons.size(); ipho++)
            {
              double pt=photons[ipho].pt();
              double eta=photons[ipho].superCluster()->eta();
              bool hasTightPhotonId = true;
              bool passId = (photons[ipho].r9()>=0.9 && hasTightPhotonId);
              bool passIso(true);
              float relIso = photons[ipho].particleIso() / photons[ipho].pt();
              if(relIso<0.4)passIso = true;
              if(pt<triggerThreshold || fabs(eta)>1.4442 ) continue;
              if(!passId) continue;
              if(!passIso) continue;
              selPhotons.push_back(photons[ipho]);
            }
        }

      ///////////////////////
      //                   //
      //  LEPTON ANALYSIS  //
      //                   //
      ///////////////////////
      
      //start by merging electrons and muons
      std::vector<patUtils::GenericLepton> leptons;
      for( size_t l=0; l<electrons.size(); l++ ){ leptons.push_back(patUtils::GenericLepton(electrons[l])); }      
      for( size_t l=0; l<muons    .size(); l++ ){ leptons.push_back(patUtils::GenericLepton(muons    [l])); }      
      std::sort(leptons.begin(),   leptons.end(), utils::sort_CandidatesByPt);

      LorentzVector muDiff(0,0,0,0);
      std::vector<patUtils::GenericLepton> selLeptons, extraLeptons;
      for(size_t ilep=0; ilep<leptons.size(); ilep++)
	{
	  bool passKin(true),passId(true),passIso(true);
	  bool passLooseLepton(true), passSoftMuon(true), passSoftElectron(true), passVetoElectron(true);

	  int lid=leptons[ilep].pdgId();

	  //apply muon corrections
	  if( abs(lid)==13 )
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

	  //veto nearby photon (loose electrons are many times photons...)
	  double minDRlg(9999.);
	  for(size_t ipho=0; ipho<selPhotons.size(); ipho++)
	    minDRlg=TMath::Min(minDRlg,deltaR(leptons[ilep].p4(),selPhotons[ipho].p4()));
	  if(minDRlg<0.1) continue;

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

	  //Isolation
	  passIso = lid==11?patUtils::passIso(leptons[ilep].el,  patUtils::llvvElecIso::Tight) : patUtils::passIso(leptons[ilep].mu,  patUtils::llvvMuonIso::Tight);
	  
	  if(passId && passIso && passKin)          selLeptons.push_back(leptons[ilep]);
	  else if(passLooseLepton || passSoftMuon)  extraLeptons.push_back(leptons[ilep]);

	}
        std::sort(selLeptons.begin(),   selLeptons.end(), utils::sort_CandidatesByPt);
        std::sort(extraLeptons.begin(), extraLeptons.end(), utils::sort_CandidatesByPt);
        LorentzVector recoMET = met - muDiff;
	
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

