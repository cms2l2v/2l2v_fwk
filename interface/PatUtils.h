#ifndef patutils_h
#define patutils_h


#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

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
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/PtrVector.h"

//need for the good lumi filter
#include "DataFormats/Provenance/interface/LuminosityBlockID.h"
#include "DataFormats/Provenance/interface/LuminosityBlockRange.h"
#include "FWCore/Utilities/interface/Algorithms.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/LumiUtils.h"

// Electron ID
#include "RecoEgamma/ElectronIdentification/interface/VersionedPatElectronSelector.h"

#include <vector>
#include "TVector3.h"
#include "TMath.h"
#include "TGraph.h"
#include <Math/VectorUtil.h>

namespace patUtils
{
   typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

   //define a generic container to hold information related to pat electrons, muons, taus
   //very weirdly, it seems it's the best way to do it...
   class GenericLepton  : public pat::GenericParticle {
      public:
      // constructor
      ~GenericLepton(){};
       GenericLepton(pat::Electron el_) : pat::GenericParticle(el_) {el = el_; };
       GenericLepton(pat::Muon     mu_) : pat::GenericParticle(mu_) {mu = mu_; };
       GenericLepton(pat::Tau     tau_) : pat::GenericParticle(tau_){tau = tau_; };
         pat::Electron el;
         pat::Muon     mu;
	 pat::Tau     tau;
   };

   namespace llvvElecId { enum ElecId  {Veto, Loose, Medium, Tight, LooseMVA, MediumMVA, TightMVA}; }
   namespace llvvMuonId { enum MuonId  {Loose, Soft, Tight, tkHighPT, TightAndTlkHighPt, StdLoose, StdSoft, StdMedium, StdTight}; }
   namespace llvvPhotonId { enum PhotonId  {Loose, Medium, Tight}; }
   namespace llvvElecIso{ enum ElecIso {Veto, Loose, Medium, Tight}; }
   namespace llvvMuonIso{ enum MuonIso {Loose,Tight, H4lWP, TightBoosted, TightAndTkRelatBoosted}; }
   namespace CutVersion { enum CutSet {Spring15Cut25ns, ICHEP16Cut, Moriond17Cut}; }

   unsigned int GainSeed (pat::Electron &electron, const EcalRecHitCollection* recHits);
   unsigned int GainSeed (pat::Photon &photon, const EcalRecHitCollection* recHits);
   bool passId (VersionedPatElectronSelector id, edm::EventBase const & event, pat::Electron el);
   bool passId (pat::Electron& el,  reco::Vertex& vtx, int IdLevel, int cutVersion, bool is2016=true); // Old PHYS14 ID
   bool passId (pat::Muon&     mu,  reco::Vertex& vtx, int IdLevel, int cutVersion);
   bool passId (pat::Photon& photon,  double rho, int IdLevel);
   float relIso(patUtils::GenericLepton& lep, double rho);
   bool passIso (VersionedPatElectronSelector id, pat::Electron& el);
   bool passIso(pat::Electron& el,  int IsoLevel, int cutVersion, double rho=0.0 ); // Old PHYS15 Iso
   bool passIso(pat::Muon&     mu,  int IsoLevel, int cutVersion);
   bool passIso(pat::Muon& mu, int IsoLevel, int cutVersion, std::vector<pat::PackedCandidate> thePats, pat::MuonCollection muons, reco::Vertex& vertex );
   bool passPhotonTrigger(fwlite::Event &ev, float &triggerThreshold, float &triggerPrescale, float& triggerThresholdHigh);
   bool passVBFPhotonTrigger(fwlite::Event &ev, float &triggerThreshold, float &triggerPrescale, float& triggerThresholdHigh);
   bool passPFJetID(std::string label, pat::Jet jet);
   bool passPUJetID(pat::Jet j);

   bool exclusiveDataEventFilter(const double&run, const bool& isMC, const bool& isPromptReco);

   std::pair<double, double> scaleVariation(const fwlite::Event& ev);
   double alphaVariation(const fwlite::Event& ev);
   double pdfVariation(const fwlite::Event& ev);

   double getHTScaleFactor(TString dtag, double lheHt);
  bool outInOnly(const reco::Muon &mu)  ;
  bool preselection(const reco::Muon &mu,bool selectClones_)  ;
  bool tighterId(const reco::Muon &mu)  ;
  bool tightGlobal(const reco::Muon &mu)  ;
  bool safeId(const reco::Muon &mu)  ;
  bool partnerId(const reco::Muon &mu)  ;
  float computePFphotonIsolation(std::vector<pat::PackedCandidate>, float , float , const pat::Muon &, pat::MuonCollection , reco::Vertex& vertex);
  float makeTkIsoCloseBySafe( float coneSize, float initialValue, const pat::Muon &theMuon, pat::MuonCollection muons);




   class MetFilter{
    private :
     struct RuLuEv {
        unsigned int Run;  unsigned int Lumi;  unsigned int Event;
        RuLuEv(unsigned int Run_, unsigned int Lumi_, unsigned int Event_){ Run = Run_; Lumi = Lumi_; Event = Event_;}
        bool operator==(const RuLuEv &other) const { return (Run == other.Run && Lumi == other.Lumi && Event == other.Event); }
     };
     struct RuLuEvHasher{
         std::size_t operator()(const RuLuEv& k) const{ using std::size_t; using std::hash;  using std::string;
            return ((hash<unsigned int>()(k.Run) ^ (hash<unsigned int>()(k.Lumi) << 1)) >> 1) ^ (hash<unsigned int>()(k.Event) << 1);
         }
     };

     typedef std::unordered_map<RuLuEv, int, RuLuEvHasher> MetFilterMap;
     MetFilterMap map;
    public :
     MetFilter(){}
     ~MetFilter(){}
     void Clear(){map.clear();}
     void FillBadEvents(std::string path);

     //     New Met Filters for 2016 Run II:
     bool passBadPFMuonFilter(const fwlite::Event& ev);
     bool passBadChargedCandidateFilter(const fwlite::Event& ev);

     int  passMetFilterInt(const fwlite::Event& ev);
     int  passMetFilterInt(const fwlite::Event& ev, bool is2016);
     bool passMetFilter(const fwlite::Event& ev);
     bool BadGlobalMuonTaggerFilter(const fwlite::Event& ev,std::unique_ptr<edm::PtrVector<reco::Muon>> &out, bool selectClones=false);
     bool BadGlobalMuonTaggerFilter(const fwlite::Event& ev,std::unique_ptr<std::vector<reco::Muon*>> &out, bool selectClones=false);
   };

}

#endif
