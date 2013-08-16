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
   llvvElectronInfoRef electronInfoRef;

   //functions

};
typedef  std::vector<llvvLepton> llvvLeptonCollection;
typedef  edm::Ref<llvvLeptonCollection> llvvLeptonRef;
typedef  edm::RefProd<llvvLeptonCollection> llvvLeptonRefProd;
typedef  edm::RefVector<llvvLeptonCollection> llvvLeptonRefVector;

class llvvTau  : public LorentzVectorF {
   public:
   // constructor
   llvvTau(){};
   ~llvvTau(){};

   //member variables
   public:
   int id,  genid, Tbits;
   bool isPF;
   uint64_t idbits;
   LorentzVectorF gen;
   float d0,        dZ,        ip3d,      ip3dsig;
   float trkchi2, trkValidPixelHits, trkValidTrackerHits, trkLostInnerHits, trkPtErr;
   std::vector<LorentzVectorF> tracks;
   std::vector<LorentzVectorF> pi0s;

   float vz, z_expo;
   float emfraction, hcalEnergy, ecalEnergy;
   LorentzVectorF jet;
   int   numChargedParticlesSigCone, numNeutralHadronsSigCone, numPhotonsSigCone, numPiZeroSigCone, numParticlesSigCone;
   int   numChargedParticlesIsoCone, numNeutralHadronsIsoCone, numPhotonsIsoCone,                   numParticlesIsoCone;
   float ptSumChargedParticlesIsoCone, ptSumPhotonsIsoCone;
   float mva_e_pi, mva_pi_mu, mva_e_mu;

   //functions
   bool passId(unsigned int IdBit){return ((idbits&(1<<IdBit))>0);}
};
typedef  std::vector<llvvTau> llvvTauCollection;
typedef  edm::Ref<llvvTauCollection> llvvTauRef;
typedef  edm::RefProd<llvvTauCollection> llvvTauRefProd;
typedef  edm::RefVector<llvvTauCollection> llvvTauRefVector;



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




///////////////////////////////////////////////////////////////////////////////////////////////////////////
// // //  NOT USED TO SAVE THE OBJECTS IN THE EDM FORMAT BUT CONVENIENT WHEN ANALYZING THE OBJECTS // // //
///////////////////////////////////////////////////////////////////////////////////////////////////////////

//Just an ugly way to save additional info to the object
struct llvvJetInfo{};
class llvvJetExt : public llvvJet{
   public:
   llvvJetExt(llvvJet jet_){
      SetPxPyPzE(jet_.px(), jet_.py(), jet_.pz(), jet_.energy());
      idbits=jet_.idbits; pfstart=jet_.pfstart; pfend=jet_.pfend;
      torawsf=jet_.torawsf;
      neutHadFrac=jet_.neutHadFrac; neutEmFrac=jet_.neutEmFrac; chHadFrac=jet_.chHadFrac; muFrac=jet_.muFrac; area=jet_.area;
      tchp=jet_.tchp; jp=jet_.jp; origcsv=jet_.origcsv; csv=jet_.csv; jpcsv=jet_.jpcsv; slcsv=jet_.slcsv; supercsv=jet_.supercsv; ssvhe=jet_.ssvhe; ivf=jet_.ivf;
      svxPx=jet_.svxPx; svxPy=jet_.svxPy; svxPz=jet_.svxPz; svxM=jet_.svxM; svxNtrk=jet_.svxNtrk; svxLxy=jet_.svxLxy; svxLxyErr=jet_.svxLxyErr;
      ivfPx=jet_.ivfPx; ivfPy=jet_.ivfPy; ivfPz=jet_.ivfPz; ivfM=jet_.ivfM; ivfNtrk=jet_.ivfNtrk; ivfLxy=jet_.ivfLxy; ivfLxyErr=jet_.ivfLxyErr;
      puMVA=jet_.puMVA; qgMVA=jet_.qgMVA;
      beta=jet_.beta; betaStar=jet_.betaStar; dRMean=jet_.dRMean; dR2Mean=jet_.dR2Mean; ptRMS=jet_.ptRMS; ptD=jet_.ptD; etaW=jet_.etaW; phiW=jet_.phiW;
      genflav=jet_.genflav; genid=jet_.genid;
      gen=jet_.gen; genj=jet_.genj;
   }
   ~llvvJetExt(){}
   public: 
   double jer; double jerup; double jerdown; double jesup; double jesdown;
};
typedef  std::vector<llvvJetExt> llvvJetExtCollection;

//ALL Sorting functions
inline bool sort_llvvObjectByPt(const LorentzVectorF &a, const LorentzVectorF &b)  { return a.pt()>b.pt(); }
inline bool sort_llvvJetByCSV(const llvvJet &a, const llvvJet &b) { return a.supercsv>b.supercsv; }

//ONLY ADD STUFF AT THE END... CAN HOST UP TO 64 VARIABLES
enum llvvTAUID {
	 decayModeFinding
	,byVLooseCombinedIsolationDeltaBetaCorr
	,byLooseCombinedIsolationDeltaBetaCorr
	,byMediumCombinedIsolationDeltaBetaCorr
	,byTightCombinedIsolationDeltaBetaCorr
	,byLooseCombinedIsolationDeltaBetaCorr3Hits
	,byMediumCombinedIsolationDeltaBetaCorr3Hits
	,byTightCombinedIsolationDeltaBetaCorr3Hits
	,byCombinedIsolationDeltaBetaCorrRaw3Hits
	,againstElectronLoose                    
	,againstElectronMedium
	,againstElectronTight
	,againstElectronMVA3category
        ,againstElectronMVA3raw
	,againstElectronLooseMVA3
	,againstElectronMediumMVA3
	,againstElectronTightMVA3
	,againstElectronVTightMVA3
	,againstMuonLoose
	,againstMuonMedium
	,againstMuonTight
	,againstMuonLoose2
	,againstMuonMedium2
	,againstMuonTight2
	,againstMuonLoose3
	,againstMuonTight3
	,byIsolationMVAraw
	,byLooseIsolationMVA
	,byMediumIsolationMVA
	,byTightIsolationMVA
	,byIsolationMVA2raw
	,byLooseIsolationMVA2
	,byMediumIsolationMVA2
	,byTightIsolationMVA2
};


#endif
