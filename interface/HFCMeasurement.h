#ifndef hfcmeasurement_hh
#define hfcmeasurement_hh

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
class HFCMeasurement
{
 public:

  enum FitTypes        { FIT_R, FIT_EB, FIT_R_AND_EB, FIT_VTB, FIT_GAMMAT, FIT_R_AND_MU, FIT_EB_AND_MU };
    
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
  HFCMeasurement(int fitType, TString fitConfig, TString wpConfig="",int maxJets=4);

  HFCMeasurement(RooWorkspace *ws,int maxJets,int fitType)
    {
      fitType_=fitType;
      switch(fitType_)
	{
	case FIT_EB:              fitTypeTitle_="#varepsilon_{b}";                    fitTypeName_="effb";                        break;
	case FIT_R_AND_EB:        fitTypeTitle_="R vs #varepsilon_{b}";               fitTypeName_="rvseffb";                     break;
	case FIT_VTB:             fitTypeTitle_="|V_{tb}|";                           fitTypeName_="vtb";                         break;
	default:                  fitTypeTitle_="R";                                  fitTypeName_="r";           fitType_=FIT_R; break;
	}

      maxJets_=maxJets;
      for(int ich=0; ich<=2; ich++)
	{
	  TString ch("ee");
	  if(ich==1) ch="mumu";
	  if(ich==2) ch="emu";
	  for(int ijet=2; ijet<=maxJets; ijet++)
	    {
	      TString cat(ch); cat+=ijet; cat+="jets";
	      sampleCats_.insert(cat.Data());
	    }
	}
      
      ws_=ws;
      mc_ = (RooStats::ModelConfig*) ws_->obj("mc");
      mc_->SetWorkspace(*ws_);
      data_= (RooDataSet *)ws_->data("data");
    }

  /**
     @short parses configuration file and adds/updates values to the model
  */
  void parseFitConfig(TString url);

  /**
     @short init the PDFs for the fit
  */
  void initHFCModel();

  /**
     @short dump fitter configuration
  */
  void printConfiguration(std::ostream &os);

  /**
     @short performs a fit with the current configurations (use with care)
   */
  HFCMeasurement::FitResult_t plrFit(TH1F *h);
  
  /**
     @short wrapper for the PLR analysis
  */
  FitResult_t plrFit(RooDataSet *data, RooStats::ModelConfig *mc,bool debug, bool systBreakup=false);

  
  /**
     @short steer the fit 
   */
  void fitHFCfrom(TH1 *, bool debug=false);
  
  inline std::string getWP()         { return wp_; }
  inline std::string getSampleType() { return sampleType_; }
  inline std::map<std::string,FitResult_t> &getResults() { return curResults_; }
  inline std::string title()         { return fitTypeTitle_.Data(); }

  /**
     @short generates the b-tag observed distribution from the current model parameters and data counts
   */
  TH1F *generateBtagObs();

  /**
     @short store the model
   */
  inline void saveWorkspace(TString url)
    {
      if(ws_==0)   return;
      if(mc_)   ws_->import(*mc_);
      if(data_) ws_->import(*data_);
      if( ws_->writeToFile(url,kTRUE) )  std::cout << "[Warn] Failed exporting HFC workspace to " << url << std::endl;
    }

  /**
     @short DTOR
  */
  ~HFCMeasurement() { }
 
 private:

  /**
     @short resets the current model values to the default ones
  */
  void resetModelValues();

  /**
     @short converts an histogram to a RooDataSet with categories
   */
  bool createDataset(TH1 *btagObs);

  /**
     @short fits the current dataset
   */
  bool fit(bool debug);

  /**
     @short the pdf has to be obtained category-by-category
  */
  TH1F *getProjectedModel(TString tag,Double_t cts,TString name,TString title);

  /**
     @short gets the probablity model as function of the POI
  */
  std::vector<TH1F *> getProbabilityModelFunction(TString tag,TString baseName);

  /**
     @short just a fancy jargon translator
   */
  std::string getNuisanceTitle(std::string &name);

  /**
     @short dumps canvas to file
  */
  void saveGraphicalResult(TCanvas *c,std::string name);

  /**
     @short adds single top information
   */
  void instantiateSingleTopContribution(RooWorkspace *);

  /**
     @short freeze/release nuisance parameters in the fit
   */
  enum FreezingMode_t { ALLNUISANCES, CORRELATEDNUISANCES, UNCORRELATEDNUISANCES };
  void freezeNuisances(RooStats::ModelConfig *mc,int mode,bool setConstant);

  int fitType_;
  TString fitTypeTitle_, fitTypeName_;  
  RooWorkspace *ws_;
  RooStats::ModelConfig *mc_;
  RooDataSet *data_;
  std::set<std::string> sampleCats_;
  std::string wp_,sampleType_;
  int maxJets_;
  std::map<std::string,FitResult_t> curResults_;
};


#endif
