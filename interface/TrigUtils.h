#ifndef trigutils_h
#define trigutils_h


#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

//Load here all the dataformat that we will need

#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

//need for the good lumi filter
#include "DataFormats/Provenance/interface/LuminosityBlockID.h"
#include "DataFormats/Provenance/interface/LuminosityBlockRange.h"
#include "FWCore/Utilities/interface/Algorithms.h"

#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/LumiUtils.h"
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"

#include <vector>

namespace trigUtils
{
  
  bool applyPrescale(float prop);

  bool passPhotonTriggerMC(fwlite::ChainEvent& ev, float &triggerThreshold,
			 float &triggerPrescale, bool &prescale);
  void photonControlSample(fwlite::ChainEvent& iEvent, pat::Photon& photon,
			  SmartSelectionMonitor& mon, TString tag);
  void photonControlEff(fwlite::ChainEvent& iEvent, pat::Photon& photon, 
			SmartSelectionMonitor& mon, TString tag);
  
}

#endif
