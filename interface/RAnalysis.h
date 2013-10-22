#ifndef _ranalysis_
#define _ranalysis_

#include "UserCode/llvv_fwk/interface/BtagUncertaintyComputer.h"
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/llvv_fwk/interface/DataEventSummaryHandler.h"

#include "TGraphErrors.h"
#include "TRandom2.h"
#include "TString.h"

#include <map>

class RAnalysis
{
  
 public:
  RAnalysis(SmartSelectionMonitor &mon,std::vector<TString> &vars);
  void setTaggingCorrectors(std::map<std::pair<TString,TString>, std::pair<TGraphErrors *,TGraphErrors *> > &btagEffCorr) { btagEffCorr_=&btagEffCorr; }
  void prepareAnalysis(data::PhysicsObjectCollection_t &leptons, data::PhysicsObjectCollection_t &jets);
  void analyze(data::PhysicsObjectCollection_t &leptons, data::PhysicsObjectCollection_t &jets, float weight, TString var, bool isTopMC);

  ~RAnalysis(){};
  
 private:
  
  SmartSelectionMonitor *mon_;
  data::PhysicsObjectCollection_t rotLeptons_;
  TRandom2 rndGen_;
  BTagSFUtil btsfutil_;
  std::map<std::pair<TString,TString>, std::pair<TGraphErrors *,TGraphErrors *> > *btagEffCorr_;
};

#endif
