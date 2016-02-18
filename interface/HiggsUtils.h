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

    //get tgraph to reweight
    TGraph* getWeightGraphFromShapes(TH1D* newLineShape, TH1D* originalLineShape,  double mH);       

    //automatically take cares of morphing and finding the right file
    TH1D* getHistoFromNRfile(std::string histoName, double mass, double Cprime, double BRnew, TFile *nrLineShapesFile);       

    //reweight the resonance
    TGraph* weightNarrowResonnance(std::string SampleName, double mass, double Cprime, double BRnew, TFile *nrLineShapesFile, double& Norm, TString pf);
    TGraph* weightGGZZContinuum(std::string SampleName, TFile *nrLineShapesFile, double& Norm, TString pf);
  

    //reweight to H125 interference
    double weightToH125Interference(double mass,double width,TFile *intFile,TString var); 

    //transverse mass
    double transverseMass(const LorentzVector &visible, const LorentzVector &invisible, bool assumeSameMass);

  }  
}
#endif
