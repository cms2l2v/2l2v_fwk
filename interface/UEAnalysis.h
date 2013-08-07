#ifndef _ueanalysis_
#define _ueanalysis_

#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/llvv_fwk/interface/DataEventSummaryHandler.h"

#include <vector>

#include "TString.h"
#include "TNtuple.h"

class UEAnalysis
{

public:
  UEAnalysis(SmartSelectionMonitor &mon);

  void analyze(data::PhysicsObjectCollection_t &leptons, 
	       data::PhysicsObjectCollection_t &jets,
	       LorentzVector &met, 
	       data::PhysicsObjectCollection_t &pf,
	       data::PhysicsObjectCollection_t &mctruth,
	       int nvtx,
	       float weight);

  void fillSummaryTuple(float evWeight) { summaryTupleVars_[2]=evWeight; summaryTuple_->Fill(summaryTupleVars_); }
  TNtuple *getSummaryTuple()            { return summaryTuple_; }

  ~UEAnalysis() { delete summaryTupleVars_; }
  
private:

  SmartSelectionMonitor *mon_;
  std::vector<TString> ueReg_;
  TNtuple *summaryTuple_;
  Float_t *summaryTupleVars_;
};

#endif
