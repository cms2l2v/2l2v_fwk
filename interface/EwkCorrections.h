#ifndef EwkCorrections_h
#define EwkCorrections_h


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

namespace EwkCorrections
{
   //typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

   //define a generic container to hold information related to pat electrons, muons, taus
   //very weirdly, it seems it's the best way to do it...
   //class GenericLepton  : public pat::GenericParticle {
   //   public:
   //   // constructor
   //   ~GenericLepton(){};
   //   GenericLepton(pat::Electron el_) : pat::GenericParticle(el_){el = el_; };
   //    GenericLepton(pat::Muon     mu_) : pat::GenericParticle(mu_){mu = mu_; };
   //      pat::Electron el;
   //      pat::Muon     mu;
   //};

   //namespace llvvElecId { enum ElecId  {Veto, Loose, Medium, Tight, LooseMVA, MediumMVA, TightMVA}; }
   //namespace llvvMuonId { enum MuonId  {Loose, Soft, Tight, StdLoose, StdSoft, StdMedium, StdTight}; }
   //namespace llvvPhotonId { enum PhotonId  {Loose, Medium, Tight}; }
   //namespace llvvElecIso{ enum ElecIso {Veto, Loose, Medium, Tight}; }
   //namespace llvvMuonIso{ enum MuonIso {Loose,Tight}; }

  std::vector<std::vector<float>> readFile_and_loadEwkTable(TString url);
  std::vector<float> findCorrection(const std::vector<std::vector<float>> & Table_EWK, float sqrt_s_hat, float t_hat);
  double getEwkCorrections(TString url, reco::GenParticleCollection genParticles);


}

#endif
