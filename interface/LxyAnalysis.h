#ifndef _lxyanalysis_
#define _lxyanalysis_

#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/llvv_fwk/interface/DataEventSummaryHandler.h"

class LxyAnalysis
{

public:
  LxyAnalysis(SmartSelectionMonitor &mon,bool runSystematics);

  void analyze(data::PhysicsObjectCollection_t & leptons, 
	       data::PhysicsObjectCollection_t &jets,
	       data::PhysicsObject_t &met, 
	       data::PhysicsObjectCollection_t &mctruth,
	       float weight);

private:

  SmartSelectionMonitor *mon_;

};

#endif
