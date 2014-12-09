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

class SignalRegionCutInfo
{
public:
  SignalRegionCutInfo();
  SignalRegionCutInfo(JSONWrapper::Object& json);

  inline bool isValid() const {return isValid_;};

private:
  std::string variable_;
  std::string property_;
  std::string direction_;
  double value_;

  bool isValid_;

  bool loadJson(JSONWrapper::Object& json);

protected:
};

typedef std::map<std::string,std::map<std::string,doubleUnc>> CutInfo;

class SignalRegion
{
public:
  SignalRegion();
  SignalRegion(JSONWrapper::Object& json);

  inline bool isValid() const {return isValid_;};

private:
  bool isValid_;
  std::string selection_;
  std::vector<SignalRegionCutInfo> cuts_;

  bool loadJson(JSONWrapper::Object& json);

protected:
};

class DatacardMaker
{
public:
  DatacardMaker();
  DatacardMaker(const std::string& jsonFile);

  inline void setVerbose()   {verbose_ = true;};
  inline void clearVerbose() {verbose_ = false;};
  inline void setUnblind()   {unblind_ = true;};
  inline void clearUnblind() {unblind_ = false;};

  bool loadJson(const std::string& jsonFile);
  bool genDatacards();

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

  bool verbose_;
  bool unblind_;
  std::string jsonLoaded_;

  std::map<std::string,bool> FileExists_;
  std::vector<std::string> printOrder_;
  std::map<std::string,std::vector<ProcessFiles>> processes_;


  bool loadJson(std::vector<JSONWrapper::Object>& selection);
  void clearSamples();

protected:
};


void printHelp(std::string binName);

/*****************************************************************************/
/* Return Codes:                                                             */
/*   0 - Everything OK                                                       */
/*   1 - Invalid arguments                                                   */
/*   2 - Problem parsing the arguments                                       */
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

    if(arg.find("--json") != std::string::npos)
    {
      jsonFile = argv[i+1];
      // TODO: Check if the file exists
      ++i;
    }

    if(arg.find("--outDir") != std::string::npos)
    {
      outDir = argv[i+1];
      ++i;
    }

    if(arg.find("--verbose") != std::string::npos)
    {
      verbose = true;
    }
  }

  if(jsonFile == "")
  {
    std::cout << "You should define at least the following arguments: json" << std::endl;
    std::cout << "For more information, consult the help (\"" << argv[0] << " --help\")" << std::endl;
    return 2;
  }

  if(verbose)
    std::cout << "Creating an object of the DatacardMaker class..." << std::endl;

  DatacardMaker myDatacardMaker;
  if(verbose)  myDatacardMaker.setVerbose();
  myDatacardMaker.genDatacards();

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

DatacardMaker::DatacardMaker()
{
  verbose_ = false;
  unblind_ = false;
  jsonLoaded_ = "";
}

DatacardMaker::DatacardMaker(const std::string& jsonFile)
{
  verbose_ = false;
  unblind_ = false;
  jsonLoaded_ = "";
  loadJson(jsonFile);
}

bool DatacardMaker::loadJson(const std::string& jsonFile)
{
  if(jsonFile == jsonLoaded_)
    return true;

  if(verbose_)
    std::cout << "DatacardMaker::loadJson(): Loading JSON - " << jsonFile << std::endl;

  JSONWrapper::Object json(jsonFile, true);

  if(loadJson(json["finalSelection"].daughters()))
  {
    jsonLoaded_ = jsonFile;
    return true;
  }
  else
  {
    jsonLoaded_ = "";
    return false;
  }

  return false;
}

bool DatacardMaker::loadJson(std::vector<JSONWrapper::Object>& selection)
{
  if(selection.size() <= 0)
  {
    std::cout << "DatacardMaker::loadJson(): No selection is defined, please check that the json file is well formed. See --help for more information." << std::endl;
    return false;
  }
  if(selection.size() > 1 && verbose_)
  {
    std::cout << "DatacardMaker::loadJson(): More than one final selection was defined. Only using the first and ignoring any other ones." << std::endl;
  }

  auto& mySelection = selection[0];

  iLumi_ = mySelection.getDouble("iLumi", 0);
  if(iLumi_ <= 0)
  {
    std::cout << "DatacardMaker::loadJson(): Integrated luminosity should be positive and non-zero." << std::endl;
    return false;
  }
  inDir_ = mySelection.getString("inDir", "");
  if(inDir_ == "")
  {
    std::cout << "DatacardMaker::loadJson(): The input directory should be defined." << std::endl;
    return false;
  }
  jsonFile_ = mySelection.getString("jsonFile", "");
  if(jsonFile_ == "")
  {
    std::cout << "DatacardMaker::loadJson(): The json file describing the samples should be defined." << std::endl;
    return false;
  }
  customExtension_ = mySelection.getString("customExtension", "");
  ttree_           = mySelection.getString("ttree", "Events");
  baseSelection_   = mySelection.getString("baseSelection", "");
  binningVariable_ = mySelection.getString("binningVariable", "");
  if(binningVariable_ == "")
  {
    std::cout << "DatacardMaker::loadJson(): You must define a binning variable for the datacards." << std::endl;
    return false;
  }
  nBins_ = mySelection.getInt("numBins", 0);
  if(nBins_ <= 0)
  {
    std::cout << "DatacardMaker::loadJson(): You must define the number of bins for the binning variable. It must be a positive and non-zero number." << std::endl;
    return false;
  }
  minVal_ = mySelection.getDouble("minVal", 0);
  maxVal_ = mySelection.getDouble("maxVal", 0);
  if(maxVal_ - minVal_ <= 0)
  {
    std::cout << "DatacardMaker::loadJson(): The range [minVal,maxVal] must be a valid interval. ie: maxVal > minVal" << std::endl;
    return false;
  }
  signalPointVariable_ = mySelection.getString("signalPointVariable", "");
  if(signalPointVariable_ == "")
  {
    std::cout << "DatacardMaker::loadJson(): You must define the signalPointVariable so that the different signal points can be distinguished, even if you only have one signal point" << std::endl;
    return false;
  }

  if(verbose_)
    std::cout << "DatacardMaker::loadJson(): Finished loading basic info from JSON. Now loading the different signal regions." << std::endl << "    Do not forget about the systematics." << std::endl;

  auto signalRegions = mySelection["signalRegions"].daughters();
  signalRegions_.clear();
  for(auto signalRegion = signalRegions.begin(); signalRegion != signalRegions.end(); ++signalRegion)
  {
    if(signalRegion->getString("signalRegionSelection", "") == "")
    {
      std::cout << "DatacardMaker::loadJson(): The signal region must have a signalRegionSelection defined. In case there is only one signal region and you do not want to apply any extra selection, I suggest using one of the base selection variables so that this value is defined." << std::endl;
      continue;
    }
    signalRegions_.push_back(SignalRegion(*signalRegion));
  }
  if(signalRegions_.size() == 0)
  {
    std::cout << "DatacardMaker::loadJson(): At least one valid signal region must be defined." << std::endl;
    return false;
  }

  // TODO - Load the Json file with the samples
  clearSamples();

  return true;
}

void DatacardMaker::clearSamples()
{
  if(verbose_)
    std::cout << "DatacardMaker::clearSamples(): Clearing the datastructures to hold the samples." << std::endl;

  for(auto iterator = processes_.begin(); iterator != processes_.end(); ++iterator)
  {
    for(auto process = iterator->second.begin(); process != iterator->second.end(); ++process)
    {
      for(auto sample = process->samples.begin(); sample != process->samples.end(); ++sample)
      {
        delete sample->chain;
      }
    }
  }

  processes_.clear();
  printOrder_.clear();

  std::vector<ProcessFiles> temp;
  processes_["BG"] = temp;
  processes_["SIG"] = temp;
  processes_["Data"] = temp;

  return;
}

SignalRegion::SignalRegion():
  isValid_(false),
  selection_("")
{
}

SignalRegion::SignalRegion(JSONWrapper::Object& json):
  isValid_(false),
  selection_("")
{
  if(loadJson(json))
  {
    isValid_ = true;
  }
  else
  {
  }
}

bool SignalRegion::loadJson(JSONWrapper::Object& json)
{
  selection_ = json.getString("signalRegionSelection", "");
  if(selection_ == "")
    return false;

  auto cuts = json["cuts"].daughters();
  for(auto cut = cuts.begin(); cut != cuts.end(); ++cut)
  {
    if(cut->getString("variableExpression", "") == "")
    {
      std::cout << "SignalRegion::loadJson(): The cut must have a valid expression. Skipping this cut." << std::endl;
      continue;
    }

    cuts_.push_back(SignalRegionCutInfo(*cut));
  }

  return true;
}

SignalRegionCutInfo::SignalRegionCutInfo():
  variable_(""),
  property_(""),
  direction_(""),
  value_(0),
  isValid_(false)
{
}

SignalRegionCutInfo::SignalRegionCutInfo(JSONWrapper::Object& json):
  variable_(""),
  property_(""),
  direction_(""),
  value_(0),
  isValid_(false)
{
  if(loadJson(json))
    isValid_ = true;
}

bool SignalRegionCutInfo::loadJson(JSONWrapper::Object& json)
{
  variable_ = json.getString("variableExpression", "");
  if(variable_ == "")
    return false;

  property_ = json.getString("variableProperty", "");

  direction_ = json.getString("cutDirection", "");
  std::transform(direction_.begin(), direction_.end(), direction_.begin(), ::tolower);
  if(direction_ == "below" || direction_ == "<")
    direction_ = "<";
  if(direction_ == "above" || direction_ == ">" || direction_ == "")
    direction_ = ">";
  if(direction_ != "<" && direction_ != ">")
  {
    std::cout << "SignalRegionCutInfo::loadJson(): Please specify a valid direction for the cutDirection" << std::endl;
    return false;
  }

  value_ = json.getDouble("value", 0.0);

  return true;
}

bool DatacardMaker::genDatacards()
{
  return false;
}
