#ifndef HiggsUtils_h
#define HiggsUtils_h

#include<iostream>
#include<vector>

#include "TGraph.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TFile.h"
#include "TMath.h"

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"

#include "Math/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "ZZMatrixElement/MELA/interface/Mela.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;


namespace higgs{

  namespace utils{

    class EventCategory {
    public:

      enum CategoryType { INCLUSIVE, EXCLUSIVEVBF, EXCLUSIVE3JETS, EXCLUSIVE3JETSVBF, EXCLUSIVE2JETS, EXCLUSIVE2JETSVBF};      

      /// Constructor
      EventCategory(int mode);
      
      /// Destructor
      ~EventCategory();

      //return current configuration
      inline int GetMode() { return mode_; }

      //classify event

      TString GetCategory(pat::JetCollection& jets, LorentzVector &boson);
     
    private :
      int mode_;
    };

    //get tgraph to reweight
    TGraph* getWeightGraphFromShapes(TH1D* newLineShape, TH1D* originalLineShape,  double mH);       

    //automatically take cares of morphing and finding the right file
    TH1D* getHistoFromNRfile(std::string histoName, double mass, double Cprime, double BRnew, TFile *nrLineShapesFile);       

    //reweight the resonance
    TGraph* weightNarrowResonnance(bool isVBF, double mass, double Cprime, double BRnew, TFile *nrLineShapesFile, double& Norm, TString pf);
    TGraph* weightGGZZContinuum(TFile *nrLineShapesFile, double& Norm, TString pf);

    double weightNarrowResonnance_MELA( Mela& mela, bool isVBF, TString MelaMode, double Cprime, double resonance, fwlite::Event& eV);  
    float ComputeInterfWeight( Mela& mela, bool isVBF, TString MelaMode, double width, double mass, SimpleParticleCollection_t& daughters, SimpleParticleCollection_t& associated, SimpleParticleCollection_t& mothers);
    //float ComputeAllWeight( Mela& mela, bool isVBF, TString MelaMode, double kFactor, double width, double mass, SimpleParticleCollection_t& daughters, SimpleParticleCollection_t& associated, SimpleParticleCollection_t& mothers);
    TGraph* Get_NNLO_kFactors();

    double weightContinuum_MELA( bool isVBF, double CP, double heavyMass);
    TGraph* Get_CPS_weights(double mass);
    
    //reweight to H125 interference
    double weightToH125Interference(double mass,double width,TFile *intFile,TString var); 

    //transverse mass
    double transverseMass(const LorentzVector &visible, const LorentzVector &invisible, bool assumeSameMass);

    inline bool sort_CandidatesByPt_V2(const SimpleParticle_t &a, const SimpleParticle_t &b) { return a.second.Pt()>b.second.Pt(); }

  }  
}
#endif
