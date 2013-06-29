#ifndef llvvObjects_H
#define llvvObjects_H

#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "Math/LorentzVector.h"
#include <vector>

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVectorF;

class llvvGenEvent 
{
   public:
   // constructor
   llvvGenEvent(){};
   ~llvvGenEvent(){};

   //member variables
   public:
  int ngenITpu, ngenOOTpu, ngenOOTpum1, ngenTruepu;
  float pthat, genWeight, qscale, x1,x2;
  int id1, id2, nup;
};
typedef  std::vector<llvvGenEvent> llvvGenEventCollection;
typedef  edm::Ref<llvvGenEventCollection> llvvGenEventRef;
typedef  edm::RefProd<llvvGenEventCollection> llvvGenEventRefProd;
typedef  edm::RefVector<llvvGenEventCollection> llvvGenEventRefVector;


class llvvGenParticle : public LorentzVectorF
{
   public:
   // constructor
   llvvGenParticle(){};
   ~llvvGenParticle(){};

   //member variables
   public:
   float lxy;
   int id, status;
};
typedef  std::vector<llvvGenParticle> llvvGenParticleCollection;
typedef  edm::Ref<llvvGenParticleCollection> llvvGenParticleRef;
typedef  edm::RefProd<llvvGenParticleCollection> llvvGenParticleRefProd;
typedef  edm::RefVector<llvvGenParticleCollection> llvvGenParticleRefVector;


class llvvMuonInfo {
   public:
   // constructor
   llvvMuonInfo(){};
   ~llvvMuonInfo(){};

   //member variables
   public:
   float nMatches, nMatchedStations, validMuonHits, innerTrackChi2, trkLayersWithMeasurement, pixelLayersWithMeasurement;

};
typedef  std::vector<llvvMuonInfo> llvvMuonInfoCollection;
typedef  edm::Ref<llvvMuonInfoCollection> llvvMuonInfoRef;
typedef  edm::RefProd<llvvMuonInfoCollection> llvvMuonInfoRefProd;
typedef  edm::RefVector<llvvMuonInfoCollection> llvvMuonInfoRefVector;


class llvvElectronInfo {
   public:
   // constructor
   llvvElectronInfo(){};
   ~llvvElectronInfo(){};

   //member variables
   public:
   bool  isConv;
   float hoe,h2te,dphiin,detain,sihih,sipip,sihip, eopin, eopout,r9,fbrem;
   float sce,sceta,scphi, ooemoop;
   float mvatrigv0, mvanontrigv0;
};
typedef  std::vector<llvvElectronInfo> llvvElectronInfoCollection;
typedef  edm::Ref<llvvElectronInfoCollection> llvvElectronInfoRef;
typedef  edm::RefProd<llvvElectronInfoCollection> llvvElectronInfoRefProd;
typedef  edm::RefVector<llvvElectronInfoCollection> llvvElectronInfoRefVector;


class llvvLepton : public LorentzVectorF
{
   public:
   // constructor
   llvvLepton(){};
   ~llvvLepton(){};

   //member variables
   public:
   int id,          idbits,    genid, Tbits;//, pid;
   int isPF;
   LorentzVectorF gen;
   float ecalIso03, hcalIso03, trkIso03;
   float gIso03,    chIso03,   puchIso03, nhIso03; 
   float ecalIso04, hcalIso04, trkIso04;
   float gIso04,    chIso04,   puchIso04, nhIso04; 
   float d0,        dZ,        ip3d,      ip3dsig;
   float trkchi2, trkValidPixelHits, trkValidTrackerHits, trkLostInnerHits, trkPtErr;
   LorentzVectorF trk;

   //Specific Lepton Information
   llvvMuonInfoRef     muonInfoRef;
   llvvElectronInfoRef     electronInfoRef;
};
typedef  std::vector<llvvLepton> llvvLeptonCollection;
typedef  edm::Ref<llvvLeptonCollection> llvvLeptonRef;
typedef  edm::RefProd<llvvLeptonCollection> llvvLeptonRefProd;
typedef  edm::RefVector<llvvLeptonCollection> llvvLeptonRefVector;

class llvvPhoton : public LorentzVectorF
{
   public:
   // constructor
   llvvPhoton(){};
   ~llvvPhoton(){};

   //member variables
   public:
   int idbits, pid;
   float ecalIso03, hcalIso03, trkIso03, gIso03,    chIso03,   puchIso03, nhIso03;
   float ecalIso04, hcalIso04, trkIso04, gIso04,    chIso04,   puchIso04, nhIso04;

   bool  isConv;
   float hoe,h2te,sihih,sipip,sihip,r9;
   float sce,sceta,scphi;
};
typedef  std::vector<llvvPhoton> llvvPhotonCollection;
typedef  edm::Ref<llvvPhotonCollection> llvvPhotonRef;
typedef  edm::RefProd<llvvPhotonCollection> llvvPhotonRefProd;
typedef  edm::RefVector<llvvPhotonCollection> llvvPhotonRefVector;

class llvvJet : public LorentzVectorF
{
   public:
   // constructor
   llvvJet(){};
   ~llvvJet(){};

   //member variables
   public:
   int idbits, pfstart, pfend;
   float torawsf;
   float neutHadFrac, neutEmFrac, chHadFrac, muFrac, area;
   float tchp, jp, origcsv, csv, jpcsv, slcsv, supercsv, ssvhe, ivf;
   float svxPx, svxPy, svxPz, svxM, svxNtrk, svxLxy, svxLxyErr;
   float ivfPx, ivfPy, ivfPz, ivfM, ivfNtrk, ivfLxy, ivfLxyErr;
   float puMVA, qgMVA;
   float beta,betaStar, dRMean, dR2Mean, ptRMS,ptD, etaW, phiW;
   int   genflav, genid;
   LorentzVectorF gen, genj;
};
typedef  std::vector<llvvJet> llvvJetCollection;
typedef  edm::Ref<llvvJetCollection> llvvJetRef;
typedef  edm::RefProd<llvvJetCollection> llvvJetRefProd;
typedef  edm::RefVector<llvvJetCollection> llvvJetRefVector;

class llvvMet : public LorentzVectorF
{
   public:
   // constructor
   llvvMet(){};
   ~llvvMet(){};

   //member variables
   public:
   float sig,sigx2,sigxy,sigy2;
};
typedef  std::vector<llvvMet> llvvMetCollection;
typedef  edm::Ref<llvvMetCollection> llvvMetRef;
typedef  edm::RefProd<llvvMetCollection> llvvMetRefProd;
typedef  edm::RefVector<llvvMetCollection> llvvMetRefVector;


class llvvPFParticle : public LorentzVectorF
{
   public:
   // constructor
   llvvPFParticle(){};
   ~llvvPFParticle(){};

   //member variables
   public:
   int id, charge;
};
typedef  std::vector<llvvPFParticle> llvvPFParticleCollection;
typedef  edm::Ref<llvvPFParticleCollection> llvvPFParticleRef;
typedef  edm::RefProd<llvvPFParticleCollection> llvvPFParticleRefProd;
typedef  edm::RefVector<llvvPFParticleCollection> llvvPFParticleRefVector;

#endif




