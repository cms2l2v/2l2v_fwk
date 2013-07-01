#ifndef _ueanalysis_
#define _ueanalysis_

#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/llvv_fwk/interface/DataEventSummaryHandler.h"

#include <vector>

#include "TString.h"

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

private:

  SmartSelectionMonitor *mon_;
  std::vector<TString> ueReg_;

};

#endif
