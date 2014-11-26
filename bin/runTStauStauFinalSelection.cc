// <author>Cristóvão B da Cruz e Silva</author>
// <email>c.beirao@cern.ch</email>
// <date>2014-11-24</date>
// <summary>Yes, this applies a final selection, and it is as analysis independent as possible</summary>

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"


#include "UserCode/llvv_fwk/interface/JSONWrapper.h"
#include "UserCode/llvv_fwk/interface/llvvObjects.h"
#include "UserCode/llvv_fwk/interface/tdrstyle.h"

#include "UserCode/llvv_fwk/interface/doubleWithUncertainty.h"


#include "TROOT.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TEventList.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TEventList.h"
#include "TNtuple.h"
#include "TGraphErrors.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRotation.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <cassert>
#include <map>
#include <algorithm>
#include <cmath>


#ifndef DEBUG_EVENT
//#define DEBUG_EVENT true
#endif

#define NAN_WARN(X) if(std::isnan(X)) std::cout << "  Warning: " << #X << " is nan" << std::endl;

class MyStyle
{
public:
  MyStyle()
  {
    lcolor_ = 1;//kBlack;
    mcolor_ = 1;//kBlack;
    fcolor_ = 0;//kWhite;
    lwidth_ = -1;
    lstyle_ = -1;
    marker_ = -1;
  };

  inline int& marker(){return marker_;};
  inline int& lcolor(){return lcolor_;};
  inline int& mcolor(){return mcolor_;};
  inline int& fcolor(){return fcolor_;};
  inline int& lwidth(){return lwidth_;};
  inline int& lstyle(){return lstyle_;};

private:
  int lcolor_;
  int mcolor_;
  int fcolor_;
  int lwidth_;
  int lstyle_;
  int marker_;

protected:
};

struct SampleFiles
{
  std::string name;
  TChain* chain;
  int nFiles;
};

struct ProcessFiles
{
  std::string name;
  std::vector<SampleFiles> samples;
  MyStyle style;
};

struct UserCutInfo
{
  std::string variable;
  std::string property;
  std::string direction;
  double value;
};

typedef std::map<std::string,std::map<std::string,doubleUnc>> CutInfo;

class SignalRegion
{
public:

private:
  std::string selection_;
  std::vector<UserCutInfo> cuts_;

protected:
};

class DatacardMaker
{
public:

private: // TODO: Add Channels
  double iLumi_;
  std::string inDir_;
  std::string outDir_;
  std::string jsonFile_;
  std::string customExtension_;
  std::string ttree_;
  std::string baseSelection_;
  std::string binningVariable_;
  int nBins_;
  double minVal_;
  double maxVal_;
  std::string signalPointVariable_;
  std::vector<SignalRegion> signalRegions_;

protected:
};


void printHelp(std::string binName);

/*****************************************************************************/
/* Return Codes:                                                             */
/*   0 - Everything OK                                                       */
/*   1 - Missing parameters_cfg.py configuration file                        */
/*****************************************************************************/
int main(int argc, char** argv)
{
//  AutoLibraryLoader::enable();
  setTDRStyle();

  std::string jsonFile;
  std::string outDir = "./OUT/";
  bool verbose = false;

  // Parse the command line options
  for(int i = 1; i < argc; ++i)
  {
    if(argv[i][0] != '-')
      continue;

    std::string arg(argv[i]);

    if(arg.find("--help") != std::string::npos)
    {
      printHelp(argv[0]);
      return 0;
    }
  }

  return 0;
}

void printHelp(std::string binName)
{
  std::cout << "Usage: " << binName << " [options]" << std::endl;

  std::cout << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "  --help  - print this help message" << std::endl;
  std::cout << "  --outDir  - set the output directory (default: ./OUT/)" << std::endl;
  std::cout << "  --json  - configuration file specifying the final selection in the several signal regions and other auxiliary data. Example file: $CMSSW_BASE/src/UserCode/llvv_fwk/test/TStauStau/finalSelection.json" << std::endl;
  std::cout << "  --verbose  - give more output text, can be useful for debugging" << std::endl;
  std::cout << "  --unblind  - include data yields in the output datacards (Not implemented yet)" << std::endl; // TODO: Implement the unblind option

  return;
}
