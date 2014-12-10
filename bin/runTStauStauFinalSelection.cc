// <author>Cristóvão B da Cruz e Silva</author>
// <email>c.beirao@cern.ch</email>
// <date>2014-11-24</date>
// <summary>Yes, this applies a final selection, and it is as analysis independent as possible</summary>

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
// TODO: Implement shape correctly, using the shape analysis from combine and providing histograms


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
  TTree* tree;
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
  inline std::string signalSelection() const {return (isValid_)?(selection_):("");};

private:
  bool isValid_;
  std::string selection_;
  std::vector<SignalRegionCutInfo> cuts_;

  bool loadJson(JSONWrapper::Object& json);

protected:
};

class ChannelInfo
{
public:
  ChannelInfo();
  ChannelInfo(JSONWrapper::Object& json);

  inline bool isValid() const {return isValid_;};
  inline std::string name() const {return name_;};
  inline std::string selection() const {return selection_;};

private:
  bool isValid_;
  std::string name_;
  std::string selection_;

  bool loadJson(JSONWrapper::Object& json);

protected:
};

class DatacardMaker
{
public:
  DatacardMaker();
  DatacardMaker(const std::string& jsonFile);
  ~DatacardMaker();

  inline void setVerbose()   {verbose_ = true;};
  inline void clearVerbose() {verbose_ = false;};
  inline void setUnblind()   {unblind_ = true;};
  inline void clearUnblind() {unblind_ = false;};

  bool loadJson(const std::string& jsonFile);
  inline void setOutDir(const std::string& outDir) {outDir_ = outDir;};
  bool genDatacards();

private: // TODO: Add Channels
  double iLumi_;
  std::string inDir_;
  std::string outDir_;
  std::string jsonFile_;
  std::string customExtension_;
  std::string ttree_;
  std::string baseSelection_;
  std::string shapeVariable_;
  int nBins_;
  double minVal_;
  double maxVal_;
  std::string signalPointVariable_;
  std::vector<SignalRegion> signalRegions_;
  std::vector<ChannelInfo> channels_;

  bool verbose_;
  bool unblind_;
  std::string jsonLoaded_;

  std::map<std::string,bool> FileExists_;
  std::vector<std::string> printOrder_;
  std::map<std::string,std::vector<ProcessFiles>> processes_;

  TFile* scratchArea_;


  bool loadJson(std::vector<JSONWrapper::Object>& selection);
  void clearSamples();
  std::vector<int> getSignalPoints(std::string currentSelection);

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
  AutoLibraryLoader::enable();
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
  myDatacardMaker.loadJson(jsonFile);
  myDatacardMaker.setOutDir(outDir);
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

DatacardMaker::~DatacardMaker()
{
  std::cout << "The list of ignored files, either missing or corrupt, can be found below:" << std::endl;
  for(auto key = FileExists_.begin(); key != FileExists_.end(); ++key)
  {
    if(!key->second)
      std::cout << "  " << key->first << std::endl;
  }
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
  shapeVariable_   = mySelection.getString("shapeVariable", "");
  if(shapeVariable_ == "")
  {
    std::cout << "DatacardMaker::loadJson(): No shape variable defined for the datacards." << std::endl;
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

  auto channels = mySelection["channels"].daughters();
  channels_.clear();
  for(auto channel = channels.begin(); channel != channels.end(); ++channel)
  {
    ChannelInfo temp(*channel);

    if(temp.isValid())
      channels_.push_back(temp);
    else
      std::cout << "DatacardMaker::loadJson(): The channel must have a name and a selection defined." << std::endl;
  }
  if(channels_.size() == 0)
    std::cout << "DatacardMaker::loadJson(): No channels were defined, assuming a default channel with name 'channel'." << std::endl;

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

  clearSamples();

  JSONWrapper::Object json(jsonFile_, true);
  std::vector<JSONWrapper::Object> processes = json["proc"].daughters();
  for(auto process = processes.begin(); process != processes.end(); ++process)
  {
    bool isData = (*process)["isdata"].toBool();
    bool isSig  = !isData && (*process).isTag("spimpose") && (*process)["spimpose"].toBool();
    bool isMC   = !isData && !isSig;

    std::string type = "Data";
    if(isMC)
      type = "BG";
    if(isSig)
      type = "SIG";

    ProcessFiles tempProc;
    tempProc.name = (*process).getString("tag", "Sample");
    printOrder_.push_back(tempProc.name);

    if(process->isTag("color"))
    {
      tempProc.style.lcolor() = process->getInt("color");
      tempProc.style.mcolor() = tempProc.style.lcolor();
      tempProc.style.fcolor() = tempProc.style.lcolor();
    }
    if(process->isTag("lcolor"))
      tempProc.style.lcolor() = process->getInt("lcolor");
    if(process->isTag("mcolor"))
      tempProc.style.mcolor() = process->getInt("mcolor");
    if(process->isTag("fcolor"))
      tempProc.style.fcolor() = process->getInt("fcolor");
    if(process->isTag("fill"))
      tempProc.style.fcolor() = process->getInt("fill");
    if(process->isTag("lwidth"))
      tempProc.style.lwidth() = process->getInt("lwidth");
    if(process->isTag("lstyle"))
      tempProc.style.lstyle() = process->getInt("lstyle");
    if(process->isTag("marker"))
      tempProc.style.marker() = process->getInt("marker");

    std::string filtExt;
    if((*process).isTag("mctruthmode"))
    {
      std::stringstream buf;
      buf << "_filt" << (*process)["mctruthmode"].toInt();
      buf >> filtExt;
    }
    std::vector<JSONWrapper::Object> samples = (*process)["data"].daughters();
    for(auto sample = samples.begin(); sample != samples.end(); ++sample)
    {
      SampleFiles tempSample;
      tempSample.nFiles = 0;
      tempSample.name = (*sample).getString("dtag", "");
      tempSample.chain = new TChain(ttree_.c_str(), ((*sample).getString("dtag", "") + (*sample).getString("suffix", "")).c_str());
      int nFiles = (*sample).getInt("split", 1);

      for(int filen = 0; filen < nFiles; ++filen)
      {
        std::string segmentExt;
        if(nFiles != 1)
        {
          std::stringstream buf;
          buf << "_" << filen;
          buf >> segmentExt;
        }

        std::string fileName = inDir_ + "/" + (*sample).getString("dtag", "") + (*sample).getString("suffix", "") + segmentExt + filtExt + customExtension_ + ".root";

        TFile* file = new TFile(fileName.c_str(), "READONLY");
        bool& fileExists = FileExists_[fileName];

        if(!file || file->IsZombie() || !file->IsOpen() || file->TestBit(TFile::kRecovered))
        {
          fileExists = false;
          file->Close();
          delete file;
          continue;
        }
        else
        {
          fileExists = true;
          file->Close();
          delete file;
        }

        tempSample.chain->Add(fileName.c_str());
        ++(tempSample.nFiles);
      }

      tempProc.samples.push_back(tempSample);
    }

    processes_[type].push_back(tempProc);
  }

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
//  std::cout << "SignalRegion::loadJson(): signal Region Selection: " << selection_  << std::endl;
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

//TDirectory* cwd = gDirectory;
//cwd->cd();
bool DatacardMaker::genDatacards()
{
  TDirectory* cwd = gDirectory;
  scratchArea_ = new TFile(".finalSelectionScratchArea.root", "RECREATE");

  for(auto signalRegion = signalRegions_.begin(); signalRegion != signalRegions_.end(); ++signalRegion)
  {
    std::vector<int> signalPoints = getSignalPoints(signalRegion->signalSelection());
    for(auto point = signalPoints.begin(); point != signalPoints.end(); ++point)
    {
      std::stringstream temp;
      std::string fileName;

      temp << outDir_ << "/SignalPoint_" << *point << ".txt";
      temp >> fileName;

      bool exists = std::ifstream(fileName).good();
      if(exists)
        std::cout << "DatacardMaker::genDatacards(): The defined signal regions overlap, program execution will continue but only the last signal region for a given signal point will be kept. Please check " << jsonFile_ << " for errors." << std::endl;
    }

    //Todo: Do whatever has to be done for each signal point in the signal region
    for(auto point = signalPoints.begin(); point != signalPoints.end(); ++point)
    {
      
    }
  }

  cwd->cd();

  return false;
}

std::vector<int> DatacardMaker::getSignalPoints(std::string currentSelection)
{
  std::vector<int> retVal;

  std::string selection = baseSelection_ + "&&" + currentSelection;

  int maxVal = 501*1000+1; //Change-me if you are having problems with multipart signal samples and it only finds some of them
  TH1D* signalPoints = new TH1D("signalPoints", "signalPoints", maxVal, 0, maxVal);

  for(auto process = processes_["SIG"].begin(); process != processes_["SIG"].end(); ++process)
  {
    for(auto sample = process->samples.begin(); sample != process->samples.end(); ++sample)
    {
      if(sample->nFiles == 0)
        continue;

      TH1D* temp_hist = new TH1D("temp_hist", "temp_hist", maxVal, 0, maxVal);
      sample->chain->Draw((signalPointVariable_+">>temp_hist").c_str(), selection.c_str(), "goff");
      signalPoints->Add(temp_hist);
      delete temp_hist;
    }
  }

  int nBins = signalPoints->GetNbinsX();
  for(int i = 1; i <= nBins; ++i)
    if(signalPoints->GetBinContent(i) != 0)
    {
//      std::cout << signalPoints->GetBinLowEdge(i) << std::endl;
      retVal.push_back(int(signalPoints->GetBinLowEdge(i)));
    }

  delete signalPoints;

  return retVal;
}

ChannelInfo::ChannelInfo():
  isValid_(false),
  name_(""),
  selection_("")
{
}

ChannelInfo::ChannelInfo(JSONWrapper::Object& json):
  isValid_(false),
  name_(""),
  selection_("")
{
  if(loadJson(json))
    isValid_ = true;
}

bool ChannelInfo::loadJson(JSONWrapper::Object& json)
{
  name_ = json.getString("name", "");
  selection_ = json.getString("selection", "");

  if(name_ == "" || selection_ == "")
    return false;

  return true;
}
