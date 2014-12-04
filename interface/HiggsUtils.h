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

#include "Math/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

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

    //reweight the resonance
    TGraph* weightNarrowResonnance(std::string SampleName, double m_gen, double mass, double Cprime, double BRnew, TGraph* hLineShapeNominal, TF1 *decayProbPdf, TFile *nrLineShapesFile=0,TString pf="");

    //reweight to H125 interference
    double weightToH125Interference(double mass,double width,TFile *intFile,TString var); 

    //transverse mass
    double transverseMass(LorentzVector &visible, LorentzVector &invisible, bool assumeSameMass);

  }  
}
#endif
