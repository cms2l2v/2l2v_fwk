#include "UserCode/llvv_fwk/interface/PatUtils.h"

namespace patUtils
{
   bool passId(pat::Electron& el,  reco::Vertex& vtx, int IdLevel){

            //for electron Id look here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
            //for the meaning of the different cuts here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification
            float dEtaln         = fabs(el.deltaEtaSuperClusterTrackAtVtx());
            float dPhiln         = fabs(el.deltaPhiSuperClusterTrackAtVtx());
            float sigmaletaleta  = el.sigmaIetaIeta();
            float hem            = el.hadronicOverEm();
            double resol         = fabs((1/el.ecalEnergy())-(el.eSuperClusterOverP()/el.ecalEnergy()));
            double dxy           = fabs(el.gsfTrack()->dxy(vtx.position()));
            double dz            = fabs(el.gsfTrack()->dz(vtx.position())); 
            double mHits         = el.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);

            bool barrel = (fabs(el.superCluster()->eta()) <= 1.479);
            bool endcap = (!barrel && fabs(el.superCluster()->eta()) < 2.5);

            switch(IdLevel){
               case llvvElecId::Veto :
                  if(barrel && dEtaln < 0.0200 && dPhiln < 0.2579 && sigmaletaleta < 0.0125 && hem < 0.2564 && dxy < 0.0250 && dz < 0.5863 && resol < 0.1508 && mHits <=2)return true;
                  if(endcap && dEtaln < 0.0141 && dPhiln < 0.2591 && sigmaletaleta < 0.0371 && hem < 0.1335 && dxy < 0.2232 && dz < 0.9513 && resol < 0.1542 && mHits <= 2)return true;
                  break;

               case llvvElecId::Loose :
                  if(barrel && dEtaln < 0.0181 && dPhiln < 0.0936 && sigmaletaleta < 0.0123 && hem < 0.141 && dxy < 0.0166 && dz < 0.54342 && resol < 0.1043 && mHits <= 1)return true; 
                  if(endcap && dEtaln < 0.0124 && dPhiln < 0.0642 && sigmaletaleta < 0.035 && hem < 0.1115 && dxy < 0.0980 && dz < 0.91870 && resol < 0.1443 && mHits <= 1)return true; 
                  break;

                case llvvElecId::Medium :
                  if(barrel && dEtaln < 0.0106 && dPhiln < 0.0323 && sigmaletaleta < 0.0107 && hem < 0.067 && dxy < 0.0131 && dz < 0.22310 && resol < 0.1043 && mHits <= 1)return true; 
                  if(endcap && dEtaln < 0.0108 && dPhiln < 0.0455 && sigmaletaleta < 0.0318 && hem < 0.097 && dxy < 0.0845 && dz < 0.75230 && resol < 0.1201 && mHits <= 1)return true; 
                  break;
  
               case llvvElecId::Tight :
                  if(barrel && dEtaln < 0.0091 && dPhiln < 0.0310 && sigmaletaleta < 0.0106 && hem < 0.0532 && dxy < 0.0126 && dz < 0.0116 && resol < 0.0609 && mHits <= 1)return true; 
                  if(endcap && dEtaln < 0.0106 && dPhiln < 0.0359 && sigmaletaleta < 0.0305 && hem < 0.0835 && dxy < 0.0163 && dz < 0.5999 && resol < 0.1126 && mHits <= 1)return true; 
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
 
    switch(IdLevel){
    case llvvPhotonId::Loose :

      if ( barrel
	   && !elevto
	   && hoe < 0.032
	   && sigmaIetaIeta < 0.0100
	   && TMath::Max(chIso-chArea*rho,0.0) < 2.94 
	   && TMath::Max(nhIso-nhArea*rho,0.0) < 3.16 + 0.0023*pt
	   && TMath::Max(gIso-gArea*rho,  0.0) < 4.43 + 0.0004*pt )
	return true; 
      if ( endcap
	   && !elevto
	   && hoe < 0.023
	   && sigmaIetaIeta < 0.0270
	   && TMath::Max(chIso-chArea*rho,0.0) < 3.07 
	   && TMath::Max(nhIso-nhArea*rho,0.0) < 17.16 + 0.0116*pt
	   && TMath::Max(gIso-gArea*rho,  0.0) < 2.11 + 0.0037*pt )
	return true; 
            
      break;
      
    case llvvPhotonId::Medium :

      if ( barrel
	   && !elevto
	   && hoe < 0.020
	   && sigmaIetaIeta < 0.0099
	   && TMath::Max(chIso-chArea*rho,0.0) < 2.62 
	   && TMath::Max(nhIso-nhArea*rho,0.0) < 2.69 + 0.0023*pt
	   && TMath::Max(gIso-gArea*rho,  0.0) < 1.35 + 0.0004*pt )
	return true; 
      if ( endcap
	   && !elevto
	   && hoe < 0.011
	   && sigmaIetaIeta < 0.0269
	   && TMath::Max(chIso-chArea*rho,0.0) < 1.40 
	   && TMath::Max(nhIso-nhArea*rho,0.0) < 4.92 + 0.0116*pt
	   && TMath::Max(gIso-gArea*rho,  0.0) < 2.11 + 0.0037*pt )
	return true; 
            
      break;
    case llvvPhotonId::Tight :

      if ( barrel
	   && !elevto
	   && hoe < 0.012
	   && sigmaIetaIeta < 0.0098
	   && TMath::Max(chIso-chArea*rho,0.0) < 1.91 
	   && TMath::Max(nhIso-nhArea*rho,0.0) < 2.55 + 0.0023*pt
	   && TMath::Max(gIso-gArea*rho,  0.0) < 1.29 + 0.0004*pt )
	return true; 
      if ( endcap
	   && !elevto
	   && hoe < 0.011
	   && sigmaIetaIeta < 0.0264
	   && TMath::Max(chIso-chArea*rho,0.0) < 1.26 
	   && TMath::Max(nhIso-nhArea*rho,0.0) < 2.71 + 0.0116*pt
	   && TMath::Max(gIso-gArea*rho,  0.0) < 1.91 + 0.0037*pt )
	return true; 
            
      break;
      
    default:
      printf("FIXME PhotonId llvvPhotonId::%i is unkown\n", IdLevel);
      return false;
      break;
      
    }    
    
    return false; 
  }
  
   bool passIso(pat::Electron& el, int IsoLevel){
          //https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
          float  chIso   = el.pfIsolationVariables().sumChargedHadronPt;
          float  nhIso   = el.pfIsolationVariables().sumNeutralHadronEt;
          float  gIso    = el.pfIsolationVariables().sumPhotonEt;
          float  puchIso = el.pfIsolationVariables().sumPUPt; 
          float  relIso  = (chIso + TMath::Max(0.,nhIso+gIso-0.5*puchIso)) / el.pt();

          bool barrel = (fabs(el.superCluster()->eta()) <= 1.479);
          bool endcap = (!barrel && fabs(el.superCluster()->eta()) < 2.5);

          switch(IsoLevel){
               case llvvElecIso::Veto :
                  if( barrel && relIso < 0.3313 ) return true;
                  if( endcap && relIso < 0.3816 ) return true;
                  break;

               case llvvElecIso::Loose :
                  if( barrel && relIso < 0.2400 ) return true;
                  if( endcap && relIso < 0.3529 ) return true;
                  break;

               case llvvElecIso::Medium :
                  if( barrel && relIso < 0.2179 ) return true;
                  if( endcap && relIso < 0.2540 ) return true;
                  break;

               case llvvElecIso::Tight :
                  if( barrel && relIso < 0.1649 ) return true;
                  if( endcap && relIso < 0.2075 ) return true;
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


}
