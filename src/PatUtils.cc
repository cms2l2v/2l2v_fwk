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
                  if(barrel && relIso>0.3313) return true;
                  if(endcap && relIso>0.3816) return true;
                  break;

               case llvvElecIso::Loose :
                  if(barrel && relIso>0.2400) return true;
                  if(endcap && relIso>0.3529) return true;
                  break;

               case llvvElecIso::Medium :
                  if(barrel && relIso>0.2179) return true;
                  if(endcap && relIso>0.2540) return true;
                  break;

               case llvvElecIso::Tight :
                  if(barrel && relIso>0.1649) return true;
                  if(endcap && relIso>0.2075) return true;
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
                  if(relIso>0.20) return true;
                  break;

               case llvvMuonIso::Tight :
                  if(relIso>0.12) return true;
                  break;

               default:
                  printf("FIXME MuonIso llvvMuonIso::%i is unkown\n", IsoLevel);
                  return false;
                  break;
          }
          return false;          
   }
}
