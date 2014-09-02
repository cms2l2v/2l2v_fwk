#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"

//#include "DataFormats/FWLite/interface/Handle.h"
//#include "DataFormats/FWLite/interface/Event.h"
//#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"


#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
//#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
//#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

//#include "DataFormats/Candidate/interface/Candidate.h"
//#include "DataFormats/VertexReco/interface/VertexFwd.h"

//#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
//#include "TrackingTools/Records/interface/TransientTrackRecord.h"
//#include "TrackingTools/IPTools/interface/IPTools.h"

//#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
//#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

//#include "UserCode/llvv_fwk/interface/llvvObjects.h"


void nInitEvents()
{
  unsigned long Total = 0;
  std::vector<std::string> urls;

  urls.push_back("~/T3-area/14_05_12_2l2nu_EDMtuples_merged/MC8TeV_DYJetsToLL_50toInf_0.root");

  for(std::vector<std::string>::iterator fileN = urls.begin(); fileN != urls.end(); ++fileN)
  {
    TFile* file = TFile::Open(fileN->c_str());
    fwlite::LuminosityBlock ls(file);
    for(ls.toBegin(); !ls.alEnd(); ++ls)
    {
      fwlite::Handle<edm::MergeableCounter> nEventsTotalCounter;
      nEventsTotalCounter.getByLabel(ls, "startCounter");
      if(!nEventsTotalCounter.isValid())
        continue;
      Total += nEventsTotalCounter->value;
    }
  }

  std::cout << "Processed " << Total << " events." << std::endl;
}
