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
    
    bool elevto = photon.hasPixelSeed();
    
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
    float chArea = utils::cmssw::getEffectiveArea(22,eta,3,"chIso"); 

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
	   && !elevto
	   && hoe < 0.028
	   && sigmaIetaIeta < 0.0107
	   && TMath::Max(chIso-chArea*rho,0.0) < 2.67 
	   && TMath::Max(nhIso-nhArea*rho,0.0) < 7.23 + exp(0.0028*pt + 0.5408) 
	   && TMath::Max(gIso-gArea*rho,  0.0) < 2.11 + 0.0014*pt )
	return true; 
      if ( endcap
	   && !elevto
	   && hoe < 0.093
	   && sigmaIetaIeta < 0.0272
	   && TMath::Max(chIso-chArea*rho,0.0) < 1.79
	   && TMath::Max(nhIso-nhArea*rho,0.0) < 8.89 + 0.01725*pt
	   && TMath::Max(gIso-gArea*rho,  0.0) < 3.09 + 0.0091*pt )
	return true; 
            
      break;
      
    case llvvPhotonId::Medium :

      if ( barrel
	   && !elevto
	   && hoe < 0.012
	   && sigmaIetaIeta < 0.0100
	   && TMath::Max(chIso-chArea*rho,0.0) < 1.79 
	   && TMath::Max(nhIso-nhArea*rho,0.0) < 0.16 + exp(0.0028*pt+0.5408) 
	   && TMath::Max(gIso-gArea*rho,  0.0) < 1.90 + 0.0014*pt )
	return true; 
      if ( endcap
	   && !elevto
	   && hoe < 0.023
	   && sigmaIetaIeta < 0.0267
	   && TMath::Max(chIso-chArea*rho,0.0) < 1.09 
	   && TMath::Max(nhIso-nhArea*rho,0.0) < 4.31 + 0.0172*pt
	   && TMath::Max(gIso-gArea*rho,  0.0) < 1.90 + 0.0091*pt )
	return true; 
            
      break;
    case llvvPhotonId::Tight :

      if ( barrel
	   && !elevto
	   && hoe < 0.010
	   && sigmaIetaIeta < 0.0100
	   && TMath::Max(chIso-chArea*rho,0.0) < 1.66 
	   && TMath::Max(nhIso-nhArea*rho,0.0) < 0.14 + exp(0.0028*pt+0.5408) 
	   && TMath::Max(gIso-gArea*rho,  0.0) < 1.40 + 0.0014*pt )
	return true; 
      if ( endcap
	   && !elevto
	   && hoe < 0.015
	   && sigmaIetaIeta < 0.0265
	   && TMath::Max(chIso-chArea*rho,0.0) < 1.04 
	   && TMath::Max(nhIso-nhArea*rho,0.0) < 3.89 + 0.0172*pt
	   && TMath::Max(gIso-gArea*rho,  0.0) < 1.40 + 0.0091*pt )
	return true; 
            
      break;
      
    default:
      printf("FIXME PhotonId llvvPhotonId::%i is unkown\n", IdLevel);
      return false;
      break;
      
    }    
    
    return false; 
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
                  if( relIso < 0.12 ) return true;
                  break;

               default:
                  printf("FIXME MuonIso llvvMuonIso::%i is unkown\n", IsoLevel);
                  return false;
                  break;
          }
          return false;          
   }

  bool passPhotonTrigger(fwlite::ChainEvent ev, float &triggerThreshold,
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
    // else if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon160_*")){
    //   hasPhotonTrigger=true;
    //   triggerThreshold=160;
    // }
    // else if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon150_*")){
    //   hasPhotonTrigger=true;
    //   triggerThreshold=150;
    // }
    // else if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon135_*")){
    //   hasPhotonTrigger=true;
    //   triggerThreshold=135;
    // }
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



  bool passMetFilters(const fwlite::Event& ev, const bool& isPromptReco){
    bool passMetFilter(false);

    edm::TriggerResultsByName metFilters = isPromptReco ? ev.triggerResultsByName("RECO") : ev.triggerResultsByName("PAT");
    
    bool CSC(     utils::passTriggerPatterns(metFilters, "Flag_CSCTightHaloFilter")); 
    bool GoodVtx( utils::passTriggerPatterns(metFilters, "Flag_goodVertices"      ));
    bool eeBadSc( utils::passTriggerPatterns(metFilters, "Flag_eeBadScFilter"     ));
    // HBHE filter needs to be complemented with , because it is messed up in data (see documentation below)
    // bool HBHE(    utils::passTriggerPatterns(metFilters, "Flag_HBHENoiseFilter"   )); // Needs to be rerun for both data (prompt+reReco) and MC, for now.
    // C++ conversion of the python FWLITE example: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
    bool HBHE(false);
    HcalNoiseSummary summary;
    fwlite::Handle <HcalNoiseSummary> summaryHandle;
    summaryHandle.getByLabel(ev, "hcalnoise");
    if(summaryHandle.isValid()) summary=*summaryHandle;
    bool failCommon(
                    summary.maxHPDHits() >= 17  ||
                    summary.maxHPDNoOtherHits() >= 10 ||
                    summary.maxZeros() >= 9e9
                    );
    // IgnoreTS4TS5ifJetInLowBVRegion is always false, skipping.
    HBHE = !(failCommon || summary.HasBadRBXTS4TS5());
    
    if(!HBHE || !CSC || !GoodVtx || !eeBadSc) passMetFilter=false;
    else                                      passMetFilter=true;
    
    return passMetFilter;
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

}
