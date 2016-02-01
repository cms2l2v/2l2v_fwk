#include "UserCode/llvv_fwk/interface/PatUtils.h"

#include "DataFormats/METReco/interface/HcalNoiseSummary.h"

namespace patUtils
{

  bool passId (VersionedPatElectronSelector id, edm::EventBase const & event, pat::Electron el){
    // This assumes an object to be created ( *before the event loop* ):
    // VersionedPatElectronSelector loose_id("some_VID_tag_including_the_WP");
    return id(el, event);
  }
  
  bool passId(pat::Electron& el,  reco::Vertex& vtx, int IdLevel){
    
    //for electron Id look here: //https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Spring15_selection_25ns 
    //for the meaning of the different cuts here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification
    float dEtaln         = fabs(el.deltaEtaSuperClusterTrackAtVtx());
    float dPhiln         = fabs(el.deltaPhiSuperClusterTrackAtVtx());
    float sigmaletaleta  = el.sigmaIetaIeta();
    float hem            = el.hadronicOverEm();
    double resol         = fabs((1/el.ecalEnergy())-(el.eSuperClusterOverP()/el.ecalEnergy()));
    double dxy           = fabs(el.gsfTrack()->dxy(vtx.position()));
    double dz            = fabs(el.gsfTrack()->dz(vtx.position())); 
    double mHits         = el.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
    bool conversionVeto = el.passConversionVeto(); 
    bool barrel = (fabs(el.superCluster()->eta()) <= 1.479);
    bool endcap = (!barrel && fabs(el.superCluster()->eta()) < 2.5);
    
    // Spring15 selection 
    switch(IdLevel){
    case llvvElecId::Veto :
         if(barrel                   &&
         dEtaln        < 0.0152   &&
         dPhiln        < 0.216    &&
         sigmaletaleta < 0.0114   &&
         hem           < 0.181    &&
         dxy           < 0.0564   &&
         dz            < 0.472    &&
         resol         < 0.207    &&
         mHits         <=2        &&
         conversionVeto            )
	return true;
         if(endcap                   &&
         dEtaln        < 0.0113   &&
         dPhiln        < 0.237    &&
         sigmaletaleta < 0.0352   &&
         hem           < 0.116    &&
         dxy           < 0.222    &&
         dz            < 0.921    &&
         resol         < 0.174    &&
         mHits <= 3               &&
         conversionVeto            )
        return true;
	break;
      
      case llvvElecId::Loose :
      if(barrel                   &&
         dEtaln        < 0.0105   &&
         dPhiln        < 0.115    &&
         sigmaletaleta < 0.0103   &&
         hem           < 0.104    &&
         dxy           < 0.0261   &&
         dz            < 0.41     &&
         resol         < 0.102    &&
         mHits         <= 2       &&
         conversionVeto            )
        return true;
      if(endcap                   &&
         dEtaln        < 0.00814  &&
         dPhiln        < 0.182    &&
         sigmaletaleta < 0.0301   &&
         hem           < 0.0897   &&
         dxy           < 0.118    &&
         dz            < 0.822    &&
         resol         < 0.126    &&
         mHits         <= 1       &&
         conversionVeto            )
        return true;
      break;

    case llvvElecId::Medium :
      if(barrel                     &&
         dEtaln          < 0.0103   &&
         dPhiln          < 0.0336   &&
         sigmaletaleta   < 0.0101   &&
         hem             < 0.0876   &&
         dxy             < 0.0118   &&
         dz              < 0.373    &&
         resol           < 0.0174   &&
         mHits           <= 2       &&
         conversionVeto             )
        return true;
      if(endcap                     &&
         dEtaln          < 0.00733  &&
         dPhiln          < 0.114    &&
         sigmaletaleta   < 0.0283   &&
         hem             < 0.0678   &&
         dxy             < 0.0739   &&
         dz              < 0.602    &&
         resol           < 0.0898   &&
         mHits            <= 1      &&
         conversionVeto             )
        return true;
      break;

    case llvvElecId::Tight :
      if(barrel                     &&
         dEtaln          < 0.00926  &&
         dPhiln          < 0.0336   &&
         sigmaletaleta   < 0.0101   &&
         hem             < 0.0597   &&
         dxy             < 0.0111   &&
         dz              < 0.0466   &&
         resol           < 0.012    &&
         mHits           <= 2       &&
         conversionVeto             )
        return true;
      if(endcap                     &&
         dEtaln          < 0.00724  &&
         dPhiln          < 0.0918   &&
         sigmaletaleta   < 0.0279   &&
         hem             < 0.0615   &&
         dxy             < 0.0351   &&
         dz              < 0.417    &&
         resol           < 0.00999  &&
         mHits           <= 1       &&
         conversionVeto             )
        return true;
      break;

    case llvvElecId::LooseMVA :
    case llvvElecId::MediumMVA :
    case llvvElecId::TightMVA :
      printf("FIXME: MVA ID not yet implemented for the electron\n");
      return false;
                  break;
                  
    default:
      printf("FIXME ElectronId llvvElecId::%i is unkown\n", IdLevel);
      return false;
      break;
    }
    return false;
  }
  
  bool passId(pat::Muon& mu,  reco::Vertex& vtx, int IdLevel){
    //for muon Id look here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#LooseMuon

    // Muon IDs for 74X are supposed to be already implemented in the standard pat::Muon methods (see https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015#Muons )
    // They are added here as "StdLoose", "StdTight" etc
    switch(IdLevel){
      
    case llvvMuonId::Loose :
      if(mu.isPFMuon() && (mu.isGlobalMuon() || mu.isTrackerMuon()))return true;
      break;
      
    case llvvMuonId::Soft :
      if(mu.isPFMuon() && mu.isTrackerMuon() && mu.muonID("TMOneStationTight") && mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 && mu.innerTrack()->hitPattern().pixelLayersWithMeasurement() > 1 &&
         fabs(mu.innerTrack()->dxy(vtx.position())) < 0.3 && fabs(mu.innerTrack()->dz(vtx.position())) < 20. && mu.innerTrack()->normalizedChi2() < 1.8) return true;
      break;
      
    case llvvMuonId::Tight :
      if( mu.isPFMuon() && mu.isGlobalMuon() && mu.globalTrack()->normalizedChi2() < 10. && mu.globalTrack()->hitPattern().numberOfValidMuonHits() > 0. && mu.numberOfMatchedStations() > 1 &&
          fabs(mu.muonBestTrack()->dxy(vtx.position())) < 0.2 && fabs(mu.muonBestTrack()->dz(vtx.position())) < 0.5 && mu.innerTrack()->hitPattern().numberOfValidPixelHits() > 0 &&
          mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5)return true;
      break;
      
    case llvvMuonId::StdLoose :
      if(mu.isLooseMuon()) return true;
      break;
      
    case llvvMuonId::StdSoft :
      if(mu.isSoftMuon(vtx)) return true;
      break;
      
    case llvvMuonId::StdTight :
      if(mu.isTightMuon(vtx)) return true;
      break;
      
    default:
      printf("FIXME MuonId llvvMuonId::%i is unkown\n", IdLevel);
      return false;
      break;
    }
    return false;
  }  
  
  bool passId(pat::Photon& photon, double rho, int IdLevel){
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonIdentificationRun2
    // CSA14 selection, conditions: 25ns, better detector alignment. 
    // Used Savvas Kyriacou's slides, mailed from Ilya. 
    
    //    bool elevto = photon.hasPixelSeed();
    
    // sigma ieta ieta
    // full5x5 is not ready in 720 yet 
    // float sigmaIetaIeta = photon.full5x5_sigmaIetaIeta();
    // taken from https://github.com/cms-sw/cmssw/blob/CMSSW_7_2_X/PhysicsTools/PatAlgos/plugins/PATPhotonSlimmer.cc#L119-L130
    
    // float sigmaIetaIeta = photon.sigmaIetaIeta(); 
    float sigmaIetaIeta = photon.userFloat("sigmaIetaIeta_NoZS"); 

    // H/E 
    float hoe = photon.hadTowOverEm();

    // isolation
    double pt=photon.pt();
    double eta=photon.superCluster()->eta();

    float chIso = photon.chargedHadronIso(); 
    //    float chArea = utils::cmssw::getEffectiveArea(22,eta,3,"chIso"); 

    float nhIso = photon.neutralHadronIso();
    float nhArea = utils::cmssw::getEffectiveArea(22,eta,3,"nhIso");

    float gIso = photon.photonIso();
    float gArea = utils::cmssw::getEffectiveArea(22,eta,3,"gIso");

    bool barrel = (fabs(eta) <= 1.479);
    bool endcap = (!barrel && fabs(eta) < 2.5);

    // PHYS14 selections PU20 bunch crossing 25 ns
    switch(IdLevel){
    case llvvPhotonId::Loose :
        
      if ( barrel
	   //	   && !elevto
	   && hoe < 0.05
	   && sigmaIetaIeta < 0.0102
	   && chIso < 3.32 
	   && TMath::Max(nhIso-nhArea*rho,0.0) < 1.92 + 0.014*pt + 0.000019*pt*pt
	   && TMath::Max(gIso-gArea*rho,  0.0) < 0.81 + 0.0053*pt  )
	return true; 
      if ( endcap
	   //	   && !elevto
	   && hoe < 0.05
	   && sigmaIetaIeta < 0.0274
	   && chIso < 1.97
	   && TMath::Max(nhIso-nhArea*rho,0.0) < 11.86 + 0.0139*pt+0.000025*pt*pt
	   && TMath::Max(gIso-gArea*rho,  0.0) < 0.83 + 0.0034*pt  )
	return true; 
            
      break;
      
    case llvvPhotonId::Medium :

      if ( barrel
	   //	   && !elevto
	   && hoe < 0.05
	   && sigmaIetaIeta < 0.0102
	   && chIso < 1.37 
	   && TMath::Max(nhIso-nhArea*rho,0.0) < 1.06 + 0.014*pt + 0.000019*pt*pt 
	   && TMath::Max(gIso-gArea*rho,  0.0) < 0.28 + 0.0053*pt  )
	return true; 
      if ( endcap
	   //	   && !elevto
	   && hoe < 0.05
	   && sigmaIetaIeta < 0.0268
	   && chIso < 1.10 
	   && TMath::Max(nhIso-nhArea*rho,0.0) < 2.69 + 0.0139*pt+0.000025*pt*pt
	   && TMath::Max(gIso-gArea*rho,  0.0) < 0.39 + 0.0034*pt  )
	return true; 
            
      break;
    case llvvPhotonId::Tight :
      
      if ( barrel    
	   && hoe < 0.05      
           && sigmaIetaIeta < 0.0100
	   && chIso < 0.76
	   && TMath::Max(nhIso-nhArea*rho,0.0) < 0.97 + 0.014*pt+0.000019*pt*pt
	   && TMath::Max(gIso-gArea*rho,  0.0) < 0.08 + 0.0053*pt )
	return true;

      if ( endcap
	   && hoe < 0.05 
	   && sigmaIetaIeta < 0.0268
	   && chIso < 0.56 
	   && TMath::Max(nhIso-nhArea*rho,0.0) < 2.09 + 0.0139*pt + 0.000025*pt*pt
	   && TMath::Max(gIso-gArea*rho,  0.0) < 0.16 + 0.0034*pt ) 
	return true; 

      break;          

    default:
      printf("FIXME PhotonId llvvPhotonId::%i is unkown\n", IdLevel);
      return false;
      break;
      
    }    
    
    return false; 
  }

  float relIso(patUtils::GenericLepton& lep, double rho){
    
    int lid=lep.pdgId();
    float relIso = 0.0; 
      
    if(lid==13){

      float  chIso   = lep.mu.pfIsolationR04().sumChargedHadronPt;
      float  nhIso   = lep.mu.pfIsolationR04().sumNeutralHadronEt;
      float  gIso    = lep.mu.pfIsolationR04().sumPhotonEt;
      float  puchIso = lep.mu.pfIsolationR04().sumPUPt;
      
      relIso  = (chIso + TMath::Max(0.,nhIso+gIso-0.5*puchIso)) / lep.mu.pt();
      
    } else if (lid==11){ 
      
      //https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Spring15_selection_25ns
      float  chIso   = lep.el.pfIsolationVariables().sumChargedHadronPt;
      float  nhIso   = lep.el.pfIsolationVariables().sumNeutralHadronEt;
      float  gIso    = lep.el.pfIsolationVariables().sumPhotonEt;
        
      if (rho == 0) {
	float  puchIso = lep.el.pfIsolationVariables().sumPUPt; 
	relIso  = (chIso + TMath::Max(0.,nhIso+gIso-0.5*puchIso)) / lep.el.pt();
      }
      else {
	float effArea = utils::cmssw::getEffectiveArea(11,lep.el.superCluster()->eta(),3);
	relIso  = (chIso + TMath::Max(0.,nhIso+gIso-rho*effArea)) / lep.el.pt();
      }
      
    }
    
    return relIso;

  }


  bool passIso (VersionedPatElectronSelector id, pat::Electron& el){
    // This assumes an object to be created ( *before the event loop* ):
    // VersionedPatElectronSelector loose_id("some_VID_tag_including_the_WP");
    return true; // Isolation is now embedded into the ID.
  }
  
  bool passIso(pat::Electron& el, int IsoLevel, double rho){
         //https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Spring15_selection_25ns
	  float  chIso   = el.pfIsolationVariables().sumChargedHadronPt;
          float  nhIso   = el.pfIsolationVariables().sumNeutralHadronEt;
          float  gIso    = el.pfIsolationVariables().sumPhotonEt;

	  float  relIso = 0.0; 
	  
	  if (rho == 0) {
	    float  puchIso = el.pfIsolationVariables().sumPUPt; 
	    relIso  = (chIso + TMath::Max(0.,nhIso+gIso-0.5*puchIso)) / el.pt();
	  }
	  else {
	    float effArea = utils::cmssw::getEffectiveArea(11,el.superCluster()->eta(),3);
	    relIso  = (chIso + TMath::Max(0.,nhIso+gIso-rho*effArea)) / el.pt();
	  }
	  
          bool barrel = (fabs(el.superCluster()->eta()) <= 1.479);
          bool endcap = (!barrel && fabs(el.superCluster()->eta()) < 2.5);

	  // Spring15 selection, conditions: PU20 bx25
          switch(IsoLevel){
               case llvvElecIso::Veto :
                  if( barrel && relIso < 0.126    ) return true;
                  if( endcap && relIso < 0.144    ) return true;
                  break;

               case llvvElecIso::Loose :
                  if( barrel && relIso < 0.0893   ) return true;
                  if( endcap && relIso < 0.121    ) return true;
                  break;

               case llvvElecIso::Medium :
                  if( barrel && relIso < 0.0766   ) return true;
                  if( endcap && relIso < 0.0678   ) return true;
                  break;

               case llvvElecIso::Tight :
                  if( barrel && relIso < 0.0354   ) return true;
                  if( endcap && relIso < 0.0646   ) return true;
                  break;

               default:
                  printf("FIXME MuonIso llvvMuonIso::%i is unkown\n", IsoLevel);
                  return false;
                  break;
          }
          return false;  
   }
  
  bool passIso(pat::Muon& mu, int IsoLevel){
    //https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Muon_Isolation
    float  chIso   = mu.pfIsolationR04().sumChargedHadronPt;
    float  nhIso   = mu.pfIsolationR04().sumNeutralHadronEt;
    float  gIso    = mu.pfIsolationR04().sumPhotonEt;
    float  puchIso = mu.pfIsolationR04().sumPUPt;
    float  relIso  = (chIso + TMath::Max(0.,nhIso+gIso-0.5*puchIso)) / mu.pt();
    
    switch(IsoLevel){
       case llvvMuonIso::Loose : 
          if( relIso < 0.20 ) return true;
	  break;
		 
       case llvvMuonIso::Tight :
         if( relIso < 0.15 ) return true;
         break;
      
       default:
         printf("FIXME MuonIso llvvMuonIso::%i is unkown\n", IsoLevel);
         return false;
         break;
    }
    return false;          
  }

  bool passPhotonTrigger(fwlite::Event &ev, float &triggerThreshold,
			 float &triggerPrescale ){
    edm::TriggerResultsByName tr = ev.triggerResultsByName("HLT");
    if( !tr.isValid() ) return false;

    bool hasPhotonTrigger(false);

    triggerPrescale = 1.0; 
    triggerThreshold = 0.0;
    

    std::string successfulPath="";
    if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon300_*")){
      hasPhotonTrigger=true;
      triggerThreshold=300;
    }
    else if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon250_*")){
      hasPhotonTrigger=true;
      triggerThreshold=250;
    }
    else if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon165_R9Id90_HE10_IsoM_*")){
      hasPhotonTrigger=true;
     triggerThreshold=165;
    }
    else if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon120_R9Id90_HE10_IsoM_*")){ //HE10_Iso40_EBOnly_*")){
      hasPhotonTrigger=true;
      triggerThreshold=120;
    }
    else if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon90_R9Id90_HE10_IsoM_*")){ // HE10_Iso40_EBOnly_*")){
      hasPhotonTrigger=true;
      triggerThreshold=90;
    }
    else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon75_R9Id90_HE10_IsoM_*")){ //HE10_Iso40_EBOnly_*")){
      hasPhotonTrigger=true;
      triggerThreshold=75;
    }
    else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon50_R9Id90_HE10_IsoM_*")){ //HE10_Iso40_EBOnly_*")){
      hasPhotonTrigger=true;
      triggerThreshold=50;
    }
    else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon36_R9Id90_HE10_IsoM_*")){ //HE10_Iso40_EBOnly_*")){
      hasPhotonTrigger=true;
      triggerThreshold=36;
    }
    else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon30_R9Id90_HE10_IsoM_*")){ //HE10_Iso40_EBOnly_*")){                       
      hasPhotonTrigger=true;                                                                                                                     
      triggerThreshold=30;                                                                                                                       
    } 
    else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon22_R9Id90_HE10_IsoM_*")){ //HE10_Iso40_EBOnly_*")){
      hasPhotonTrigger=true;
      triggerThreshold=22;
    }
      
    if(successfulPath!=""){ //get the prescale associated to it
      fwlite::Handle< pat::PackedTriggerPrescales > prescalesHandle;
      prescalesHandle.getByLabel(ev, "patTrigger");
      pat::PackedTriggerPrescales prescales = *prescalesHandle;
      const edm::TriggerResults& trResults =  prescales.triggerResults();
      prescales.setTriggerNames( ev.triggerNames(trResults) );
      triggerPrescale = prescales.getPrescaleForName(successfulPath);
    }

    return hasPhotonTrigger; 
  }

  /*
  bool passPhotonTrigger(fwlite::Event &ev, float &triggerThreshold,
			 float &triggerPrescale ){
    edm::TriggerResultsByName tr = ev.triggerResultsByName("HLT");
    if( !tr.isValid() ) return false;

    bool hasPhotonTrigger(false);
    // float triggerPrescale(1.0); 
    // float triggerThreshold(0);
    triggerPrescale = 1.0; 
    triggerThreshold = 0.0;

    std::string successfulPath="";
    if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon300_*")){
      hasPhotonTrigger=true;
      triggerThreshold=300;
    }
    else if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon250_*")){
      hasPhotonTrigger=true;
      triggerThreshold=250;
    }

    else if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_*")){
      hasPhotonTrigger=true;
      triggerThreshold=120;
    }
    else if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_*")){
      hasPhotonTrigger=true;
      triggerThreshold=90;
    }
    else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_*")){
      hasPhotonTrigger=true;
      triggerThreshold=75;
    }
    else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_*")){
      hasPhotonTrigger=true;
      triggerThreshold=50;
    }
    else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_*")){
      hasPhotonTrigger=true;
      triggerThreshold=36;
    }
    else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_*")){
      hasPhotonTrigger=true;
      triggerThreshold=22;
    }
      
    if(successfulPath!=""){ //get the prescale associated to it
      fwlite::Handle< pat::PackedTriggerPrescales > prescalesHandle;
      prescalesHandle.getByLabel(ev, "patTrigger");
      pat::PackedTriggerPrescales prescales = *prescalesHandle;
      const edm::TriggerResults& trResults =  prescales.triggerResults();
      prescales.setTriggerNames( ev.triggerNames(trResults) );
      triggerPrescale = prescales.getPrescaleForName(successfulPath);
    }

    return hasPhotonTrigger; 
  }
  */

  bool passPFJetID(std::string label,
                  pat::Jet jet){

    bool passID = false;

    float rawJetEn(jet.correctedJet("Uncorrected").energy() );
    // Note: All fractions are calculated with the raw/uncorrected energy of the jet (only then they add up to unity). So the PF JetID has to be applied before the jet energy corrections. 

    float nhf( (jet.neutralHadronEnergy() + jet.HFHadronEnergy())/rawJetEn );
    float nef( jet.neutralEmEnergy()/rawJetEn );
    float cef( jet.chargedEmEnergy()/rawJetEn );
    float chf( jet.chargedHadronEnergy()/rawJetEn );
    float nch    = jet.chargedMultiplicity();
    float nconst = jet.numberOfDaughters();
    //float muf( jet.muonEnergy()/rawJetEn);
    float NumNeutralParticles = jet.neutralMultiplicity();

    //From https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data
    if (label == "Loose"){
      if( fabs(jet.eta()) <= 3.0) passID = ( (nhf<0.99  && nef<0.99 && nconst>1) && ( fabs(jet.eta())>2.4 || (fabs(jet.eta()) <= 2.4 && chf>0 && nch>0 && cef<0.99) ) );
      if( fabs(jet.eta()) > 3.0) passID = ( nef<0.90 && NumNeutralParticles > 10);
    }
    if (label == "Tight"){ 
      if( fabs(jet.eta()) <= 3.0) passID = ( (nhf<0.90  && nef<0.90 && nconst>1) && ( fabs(jet.eta())>2.4||(fabs(jet.eta()) <= 2.4 && chf>0 && nch>0 && cef<0.99) ) );
      if( fabs(jet.eta()) > 3.0) passID = ( nef<0.90 && NumNeutralParticles > 10 );
    }
    return passID;
  }


  bool passPUJetID(pat::Jet j){

    double jpt = j.pt();
    double jeta = j.eta();

    //Recommendation of HZZ :https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4l2015#Jets
    //FIXME recommendation of HZZ are oudated (run 1)--->update them
    //FIXME Progress: Partially done (removed the cuts at eta>2.75) but have to wait for improvements from JetMET group
    //FIXME In miniAOD V2, the variable used is bugged ! Please avoir using it!
    float jpumva=0.;
    jpumva=j.userFloat("pileupJetId:fullDiscriminant");

    bool passPU = true;
    if(jpt>20){
      if(jeta>3.){
        //if(jpumva<=-0.45)passPU=false;
      }else if(jeta>2.75){
        //if(jpumva<=-0.55)passPU=false;
      }else if(jeta>2.5){
        if(jpumva<=-0.6)passPU=false;
      }else if(jpumva<=-0.63)passPU=false;
    }else{ //pt<20 : in the 2l2nu analysis, this means 15<pt<20
      if(jeta>3.){
        //if(jpumva<=-0.95)passPU=false;
      }else if(jeta>2.75){
        //if(jpumva<=-0.94)passPU=false;
      }else if(jeta>2.5){
        if(jpumva<=-0.96)passPU=false;
      }else if(jpumva<=-0.95)passPU=false;
    }
    return passPU;
  }

  bool exclusiveDataEventFilter(const double& run, const bool& isMC, const bool& isPromptReco)
  {
    bool passExclusiveDataEventFilter(false);
    
    if(isMC)
      passExclusiveDataEventFilter=true;
    else
      {
        bool isForPromptReco(run<251162 || run>251562);
        if(isPromptReco)  // Prompt reco keeps events outside of that range
          {
            if(isForPromptReco) passExclusiveDataEventFilter=true;
            else                passExclusiveDataEventFilter=false;
          }
        else // 17Jul15 ReReco keeps event inside that range
          {
            if(isForPromptReco) passExclusiveDataEventFilter=false;
            else                passExclusiveDataEventFilter=true;
          }
      }
    return passExclusiveDataEventFilter;
  }


void MetFilter::FillBadEvents(std::string path){
     unsigned int Run=0; unsigned int Lumi=1; unsigned int Event=2;
     //LOOP on the files and fill the map
     std::ifstream File;
     File.open(path.c_str());
     if(!File.good() ){
       std::cout<<"ERROR:: File "<<path<<" NOT FOUND!!"<<std::endl; 
       return;
     }

     std::string line;
     while ( File.good() ){
         getline(File, line);
         std::istringstream stringfile(line);
         stringfile >> Run >> Lumi >> Event >> std::ws;
         map[RuLuEv(Run, Lumi, Event)] += 1;
    }
    File.close();
}



int MetFilter::passMetFilterInt(const fwlite::Event& ev){  
     if(map.find( RuLuEv(ev.eventAuxiliary().run(), ev.eventAuxiliary().luminosityBlock(), ev.eventAuxiliary().event()))!=map.end())return 1;

    // Legacy: the different collection name was necessary with the early 2015B prompt reco        
//    edm::TriggerResultsByName metFilters = isPromptReco ? ev.triggerResultsByName("RECO") : ev.triggerResultsByName("PAT");
    edm::TriggerResultsByName metFilters = ev.triggerResultsByName("PAT");   //is present only if PAT (and miniAOD) is not run simultaniously with RECO
    if(!metFilters.isValid()){metFilters = ev.triggerResultsByName("RECO");} //if not present, then it's part of RECO
    if(!metFilters.isValid()){       
       printf("TriggerResultsByName for MET filters is not found in the process, as a consequence the MET filter is disabled for this event\n");    
    }


    // Documentation:    
    // -------- Full MET filters list (see bin/chhiggs/runAnalysis.cc for details on how to print it out ------------------
    // Flag_trackingFailureFilter
    // Flag_goodVertices        -------> Recommended by PAG
    // Flag_CSCTightHaloFilter  -------> Recommended by PAG
    // Flag_trkPOGFilters
    // Flag_trkPOG_logErrorTooManyClusters
    // Flag_EcalDeadCellTriggerPrimitiveFilter
    // Flag_ecalLaserCorrFilter
    // Flag_trkPOG_manystripclus53X
    // Flag_eeBadScFilter       -------> Recommended by PAG
    // Flag_METFilters
    // Flag_HBHENoiseFilter     -------> Recommended by PAG
    // Flag_trkPOG_toomanystripclus53X
    // Flag_hcalLaserEventFilter
    //
    // Notes (from https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015#ETmiss_filters ):
    // - For the RunIISpring15DR74 MC campaing, the process name in PAT.
    // - For Run2015B PromptReco Data, the process name is RECO.
    // - For Run2015B re-MiniAOD Data 17Jul2015, the process name is PAT.
    // - MET filters are available in PromptReco since run 251585; for the earlier run range (251162-251562) use the re-MiniAOD 17Jul2015
    // - Recommendations on how to use MET filers are given in https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2 . Note in particular that the HBHO noise filter must be re-run from MiniAOD instead of using the flag stored in the TriggerResults; this applies to all datasets (MC, PromptReco, 17Jul2015 re-MiniAOD)
    // -------------------------------------------------

    if(!utils::passTriggerPatterns(metFilters, "Flag_CSCTightHaloFilter")) return 2; 
    if(!utils::passTriggerPatterns(metFilters, "Flag_goodVertices"      )) return 3;
    if(!utils::passTriggerPatterns(metFilters, "Flag_eeBadScFilter"     )) return 4;

    // HBHE filter needs to be complemented with , because it is messed up in data (see documentation below)
    //if(!utils::passTriggerPatterns(metFilters, "Flag_HBHENoiseFilter"   )) return false; // Needs to be rerun for both data (prompt+reReco) and MC, for now.
    // C++ conversion of the python FWLITE example: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2

    HcalNoiseSummary summary;
    fwlite::Handle <HcalNoiseSummary> summaryHandle;
    summaryHandle.getByLabel(ev, "hcalnoise");
    if(summaryHandle.isValid()) summary=*summaryHandle;

    //HBHE NOISE
    if(summary.maxHPDHits() >= 17)return 5;
    if(summary.maxHPDNoOtherHits() >= 10)return 6;
    if( summary.maxZeros() >= 9e9) return 7;
    if(summary.HasBadRBXRechitR45Loose())return 8;  //for 25ns only!  CHeck the recipe for 50ns.
//    bool failCommon(summary.maxHPDHits() >= 17  || summary.maxHPDNoOtherHits() >= 10 || summary.maxZeros() >= 9e9 );
    // IgnoreTS4TS5ifJetInLowBVRegion is always false, skipping.
//    if((failCommon || summary.HasBadRBXRechitR45Loose() )) return 5;  //for 25ns only
   
     //HBHE ISO NOISE
     if(summary.numIsolatedNoiseChannels() >=10)return 9;
     if(summary.isolatedNoiseSumE() >=50) return 10;
     if(summary.isolatedNoiseSumEt() >=25) return 11;

     return 0;
   }


bool MetFilter::passMetFilter(const fwlite::Event& ev){
   return passMetFilterInt(ev)==0;
}

std::pair<double, double> scaleVariation(const fwlite::Event& ev){
        fwlite::Handle<LHEEventProduct> lheEPHandle;
        lheEPHandle.getByLabel( ev, "externalLHEProducer");
        double scaleUp = 1;
        double scaleDw = 1;
        std::vector<double> scaleVect;
        if( lheEPHandle.isValid() ){
                for (unsigned int i=0; i<lheEPHandle->weights().size(); i++) {
                        if( lheEPHandle->weights()[i].id != "1001" || lheEPHandle->weights()[i].id != "1006" || lheEPHandle->weights()[i].id != "1008" ){
                                double local_weight = 0;
                                local_weight = ( lheEPHandle->weights()[i].wgt / lheEPHandle->originalXWGTUP() );
                                scaleUp = std::max(scaleUp, local_weight);
                                scaleDw = std::min(scaleDw, local_weight);
                        }
                 }
         }
         return std::make_pair(scaleUp, scaleDw);
}

double pdfVariation(const fwlite::Event& ev){
        fwlite::Handle<LHEEventProduct> lheEPHandle;
        lheEPHandle.getByLabel( ev, "externalLHEProducer");
        std::vector<double> weight_vect;
        double pdfVar;
        double sum = 0;
        if( lheEPHandle.isValid() ){
                for (unsigned int i=0; i<lheEPHandle->weights().size(); i++) {
                        std::string::size_type sz;
                        double id = std::stod( lheEPHandle->weights()[i].id, &sz);
                        for( unsigned int i_id = 2002; i_id<2101; i_id++){
                                if( i_id == id ){
                                        sum += std::pow( (lheEPHandle->weights()[i].wgt / lheEPHandle->originalXWGTUP() - 1 ), 2);
                                        weight_vect.push_back(lheEPHandle->weights()[i].wgt);
                                }
                        }
                }
        }
        pdfVar = 1+sum/weight_vect.size(); //+1 variation
        return pdfVar;
}

double alphaVariation(const fwlite::Event& ev){
        fwlite::Handle<LHEEventProduct> lheEPHandle;
        lheEPHandle.getByLabel( ev, "externalLHEProducer");
        double alphaVar = 0;
        double local_alpha_one = 0;
        double local_alpha_two = 0;
        if( lheEPHandle.isValid() ){
                for (unsigned int i=0; i<lheEPHandle->weights().size(); i++) {
                        if( lheEPHandle->weights()[i].id != "2101" ){
                                local_alpha_one = ( lheEPHandle->weights()[i].wgt / lheEPHandle->originalXWGTUP() );
                        } else if( lheEPHandle->weights()[i].id != "2102" ){
                                local_alpha_two = ( lheEPHandle->weights()[i].wgt / lheEPHandle->originalXWGTUP() );
                        }
                }
        }
        alphaVar = 1+std::sqrt(0.75)*std::abs( local_alpha_one - local_alpha_two )*0.5; //+1 variation
        return alphaVar;
}

}
