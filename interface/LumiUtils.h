#ifndef lumiutils_h
#define lumiutils_h


#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

#include "DataFormats/FWLite/interface/LuminosityBlock.h"
#include "DataFormats/FWLite/interface/Run.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/Provenance/interface/LuminosityBlockID.h"
#include "DataFormats/Provenance/interface/LuminosityBlockRange.h"
#include "FWCore/Utilities/interface/Algorithms.h"


#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TEventList.h"
#include "TROOT.h"

#include <iostream>
#include <vector>

using namespace std;

namespace lumiUtils
{

   class GoodLumiFilter {
      public:
      // constructor
      ~GoodLumiFilter(){};
       GoodLumiFilter(std::vector<edm::LuminosityBlockRange> lumisToProcess_){lumisToProcess = lumisToProcess_;  sortAndRemoveOverlaps(lumisToProcess); };
       bool isGoodLumi(edm::RunNumber_t run, edm::LuminosityBlockNumber_t lumi);
       void FindLumiInFiles(vector<string>& fileNames);
       void DumpToJson(string FileName);
       void DumpToTime(string FileName);
       void RemoveRunsAfter(unsigned int RunMax);
 
      private:
         struct stRun {unsigned int runId; std::vector<unsigned int> lumiId; };

         std::vector<stRun*> RunMap;
         std::map<unsigned int, unsigned int> timeMap;

         std::vector<edm::LuminosityBlockRange> lumisToProcess;
   };


}

#endif
