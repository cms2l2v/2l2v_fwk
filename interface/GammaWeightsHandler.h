/* 
 * Wrapper for common operations on a gamma event
 * Get weights/mass shapes from file
 * Analyze event and assign trigger categories, weights and massive candidates
 * $Date: 2013/01/14 21:28:45 $
 * $Revision: 1.11 $
 * \author Pedro Silva
 */

#ifndef GammaWeightsHandler_H
#define GammaWeightsHandler_H

#include "TString.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1.h"
#include "TGraph.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "UserCode/llvv_fwk/interface/DataEventSummaryHandler.h"
#include "UserCode/llvv_fwk/interface/llvvObjects.h"

class GammaWeightsHandler
{
 public: 

  GammaWeightsHandler(const edm::ParameterSet &runProcess,TString ewkSupWgt="",bool forceAllToData=false);

  float getWeightFor(std::vector<Float_t> &vars, TString evCategoryLabel="");
  LorentzVector getMassiveP4(LorentzVector &gamma,TString evCategoryLabel="");
  LorentzVector getMassiveP4(LorentzVectorF &gamma,TString evCategoryLabel="");
  ~GammaWeightsHandler();

 private:

  std::vector<TString> dilCats_; 
  std::map<TString, std::vector<TGraph *> > wgtsH_;
  std::map<TString, TH1 *> zmassH_;  
};

#endif
