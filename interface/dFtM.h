#ifndef dftm_hh
#define dftm_hh

#include "UserCode/llvv_fwk/interface/tdrstyle.h"
#include "UserCode/llvv_fwk/interface/JSONWrapper.h"

#include "TFile.h"
#include "TRandom2.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TPaveText.h"
#include "THStack.h"

#include "RooGaussian.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooProdPdf.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooAddition.h"
#include "RooPlot.h"
#include "RooMinuit.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooLognormal.h"
#include "RooUniform.h"

#include "RooStats/ConfInterval.h"
#include "RooStats/PointSetInterval.h"
#include "RooStats/ConfidenceBelt.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/ModelConfig.h"

//
class dFtM
{
 public:

  enum FitTypes        { FIT_SFb, FIT_SFq, FIT_SFb_AND_SFq };
  
  struct FitResult_t
  {
    FitResult_t() : postFitNuisGr(0),plrGr(0) { };

    Bool_t status;
    Double_t poiFit,poiErr,poiFitLoLim,poiFitUpLim,poiFitStatLoLim,poiFitStatUpLim,poiStatErr;
    std::map<std::string,Double_t> postFitNuis, uncBreakup;
    TH1F *postFitNuisGr;
    TGraph *plrGr;
  };

  /**
     @short CTOR
   */
  dFtM(int fitType, TString flavCfg, TString effCfg, TString btagsUrl);

  /**
     @short init the PDFs for the fit
  */
  void initModel();

  /**
     @short dump fitter configuration 
  */
  void printConfiguration(std::ostream &os);

  /**
     @short performs a fit with the current configurations
  */
  dFtM::FitResult_t fit();
  
  /**
     @short getters
  */
  inline std::string getWP()         { return wp_; }
  inline std::string getSampleType() { return sampleType_; }
  inline std::map<std::string,FitResult_t> &getResults() { return curResults_; }
  inline std::string title()         { return fitTypeTitle_.Data(); }
  inline std::pair<Int_t,Int_t> getJetKinematicsForCat(std::string cat) { return catToJetKinematics_[cat]; }
    
  /**
     @short store the model
  */
  inline void save(TString url)
    {
      if(ws_==0)   return;
      if(mc_)   ws_->import(*mc_);
      if(data_) ws_->import(*data_);
      if( ws_->writeToFile(url,kTRUE) )  std::cout << "[Warn] Failed exporting dFtM workspace to " << url << std::endl;
    }

  /**
     @short gets the observed and expected b-tag multiplicity histograms per channel: ee, ee_exp, mumu, mumu_exp, emu, emu_exp
   */
  std::map<TString,TH1F *> getExtendedBtagMultiplicityHistograms();

  /**
     @short DTOR
  */
  ~dFtM() { }
  
 private:

  /**
     @short parses configuration file and adds/updates values to the model
  */
  void parseFitConfig(TString url);

  /**
     @short reads histograms from file and instantiates dataset
   */
  void readoutDataFrom(TString url);

  /**
     @short reads signal and background expectations from file and instantiates expected counts
   */
  void readoutExpectationsFrom(TString url);
    
  /**
     @short resets the current model values to the default ones
  */
  void reset();
  
  /**
     @short wrapper for the PLR analysis
  */
  FitResult_t plrFit(RooDataSet *data, RooStats::ModelConfig *mc,bool debug, bool systBreakup=false);


/*   /\** */
/*      @short the pdf has to be obtained category-by-category */
/*   *\/ */
/*   TH1F *getProjectedModel(TString tag,Double_t cts,TString name,TString title); */

/*   /\** */
/*      @short just a fancy jargon translator */
/*    *\/ */
/*   std::string getNuisanceTitle(std::string &name); */

/*   /\** */
/*      @short dumps canvas to file */
/*   *\/ */
/*   void saveGraphicalResult(TCanvas *c,std::string name); */

  int fitType_;
  TString fitTypeTitle_, fitTypeName_;  
  RooWorkspace *ws_;
  RooStats::ModelConfig *mc_;
  RooDataSet *data_;

  std::set<std::string> sampleCats_;
  std::map<std::string, std::pair<Int_t,Int_t> > catToJetKinematics_;
  std::string wp_,sampleType_;
  std::map<std::string,FitResult_t> curResults_;
};


#endif
