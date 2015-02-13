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
#include <cstdio>
#include <cctype>


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
  std::string label;
  std::vector<SampleFiles> samples;
  MyStyle style;
  bool reweight;
  bool isData;
  bool isSig;
  bool isMC;
};

class SignalRegionCutInfo
{
public:
  SignalRegionCutInfo();
  SignalRegionCutInfo(JSONWrapper::Object& json);

  inline bool isValid() const {return isValid_;};
  std::string cut() const;

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

  std::string cuts() const;

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
  ChannelInfo(std::string name, std::string selection);

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

  inline void setVerbose()        {verbose_ = true;};
  inline void clearVerbose()      {verbose_ = false;};
  inline void setUnblind()        {unblind_ = true;};
  inline void clearUnblind()      {unblind_ = false;};
  inline void setCrossSection()   {upperLimitCrossSection_ = true;};
  inline void clearCrossSection() {upperLimitCrossSection_ = false;};

  bool loadJson(const std::string& jsonFile);
  inline void setOutDir(const std::string& outDir) {outDir_ = outDir;};
  bool genDatacards();

private:
  std::string name_;
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
  bool upperLimitCrossSection_;
  std::string jsonLoaded_;

  std::map<std::string,bool> FileExists_;
  std::vector<std::string> printOrder_;
  std::map<std::string,std::vector<ProcessFiles>> processes_;

  TFile* scratchArea_;


  bool loadJson(std::vector<JSONWrapper::Object>& selection);
  void clearSamples();
  std::vector<int> getSignalPoints(std::string currentSelection);
  std::map<std::string,std::map<std::string,std::map<std::string,doubleUnc>>> applySelection(std::string type, std::vector<ProcessFiles> &processes, const SignalRegion &signalRegion, std::string additionalSelection = "", bool doSyst = false);

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
  bool unblind = false;
  bool upperLimitCrossSection = false;

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
      verbose = true;

    if(arg.find("--unblind") != std::string::npos)
      unblind = true;

    if(arg.find("--xsec") != std::string::npos)
      upperLimitCrossSection = true;
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
  if(verbose)                myDatacardMaker.setVerbose();
  if(unblind)                myDatacardMaker.setUnblind();
  if(upperLimitCrossSection) myDatacardMaker.setCrossSection();
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
  std::cout << "  --help     - print this help message" << std::endl;
  std::cout << "  --outDir   - set the output directory (default: ./OUT/)" << std::endl;
  std::cout << "  --json     - configuration file specifying the final selection in the several signal regions and other auxiliary data. Example file: $CMSSW_BASE/src/UserCode/llvv_fwk/test/TStauStau/finalSelection.json" << std::endl;
  std::cout << "  --verbose  - give more output text, can be useful for debugging" << std::endl;
  std::cout << "  --unblind  - include data yields in the output datacards (Not completely implemented yet)" << std::endl; // TODO: Finish implementing the unblind option
  std::cout << "  --xsec     - build datacards so that the result is a limit on the cross section and not a limit on the signal strength" << std::endl;

  return;
}

DatacardMaker::DatacardMaker()
{
  verbose_ = false;
  unblind_ = false;
  upperLimitCrossSection_ = false;
  jsonLoaded_ = "";
}

DatacardMaker::DatacardMaker(const std::string& jsonFile)
{
  verbose_ = false;
  unblind_ = false;
  upperLimitCrossSection_ = false;
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

  name_ = mySelection.getString("name", "");
  if(name_ == "")
  {
    std::cout << "DatacardMaker::loadJson(): You should define a name for the analysis." << std::endl;
    return false;
  }
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
  for(auto &channel : channels)
  {
    ChannelInfo temp(channel);

    if(temp.isValid())
      channels_.push_back(temp);
    else
      std::cout << "DatacardMaker::loadJson(): The channel must have a name and a selection defined." << std::endl;
  }
  if(channels_.size() == 0)
  {
    std::cout << "DatacardMaker::loadJson(): No channels were defined, assuming a default channel with name 'channel'." << std::endl;
    channels_.push_back(ChannelInfo("channel", baseSelection_));
  }

  auto signalRegions = mySelection["signalRegions"].daughters();
  signalRegions_.clear();
  for(auto &signalRegion : signalRegions)
  {
    SignalRegion temp(signalRegion);

    if(temp.isValid())
      signalRegions_.push_back(temp);
    else
      std::cout << "DatacardMaker::loadJson(): The signal region must have a signalRegionSelection defined. In case there is only one signal region and you do not want to apply any extra selection, I suggest using one of the base selection variables so that this value is defined." << std::endl;
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
    tempProc.isData = isData;
    tempProc.isSig = isSig;
    tempProc.isMC = isMC;
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

    tempProc.reweight = process->getBool("reweight", true);

    tempProc.label = process->getString("label", "");
    if(tempProc.label == "")
      tempProc.label = tempProc.name;

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

  jsonLoaded_ = "";

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
  for(auto &cut : cuts)
  {
    SignalRegionCutInfo temp(cut);

    if(temp.isValid())
      cuts_.push_back(temp);
    else
      std::cout << "SignalRegion::loadJson(): The cut must have a valid expression. Skipping this cut." << std::endl;
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

std::map<std::string,std::map<std::string,std::map<std::string,doubleUnc>>> DatacardMaker::applySelection(std::string type, std::vector<ProcessFiles> &processes, const SignalRegion &signalRegion, std::string additionalSelection, bool doSyst)
{
  std::map<std::string,std::map<std::string,std::map<std::string,doubleUnc>>> retVal;
  // retVal[channel][process][systematic]
  // without systematic is called "noSyst"

  // type == "S"  -> Signal sample
  // type == "B"  -> Background sample
  // type == "D"  -> Data sample
  if(type.size() != 0)
    type = toupper(type[0]);
  else
    type = "B"; // Assume background by default

  std::string SRSelection = signalRegion.cuts();
  std::string weight = "weight";
  if(upperLimitCrossSection_ && type == "S")
    weight = "(weight/crossSection)";
  for(auto &channel : channels_)
  {
    std::string channelSelection = channel.selection();


    std::string selection;
    selection  = "((" + baseSelection_ + ")";
    selection += " && (" + channelSelection + ")";
    if(SRSelection != "")
      selection += " && (" + SRSelection + ")";
    if(additionalSelection != "")
      selection += " && (" + additionalSelection + ")";
    selection += ")";

    std::map<std::string,std::map<std::string,doubleUnc>> channelYields;
    for(auto &process : processes)
    {
      std::map<std::string,doubleUnc> yields;

      TH1D processHist("processHist", "processHist", 1, 0, 20);
      processHist.Sumw2();
      for(auto &sample : process.samples)
      {
        TH1D tempHist("tempHist", "tempHist", 1, 0, 20);
        tempHist.Sumw2();

        sample.chain->Draw("weight>>tempHist", (selection+"*"+weight).c_str(), "goff");

        if(process.reweight && !process.isData)
          tempHist.Scale(1.0/sample.nFiles);

        processHist.Add(&tempHist);
      }
      if(!process.isData)
        processHist.Scale(iLumi_);

      doubleUnc yield(0,0);
      yield.setValue(processHist.GetBinContent(0) + processHist.GetBinContent(1) + processHist.GetBinContent(2));
      TArrayD* w2Vec = processHist.GetSumw2();
      yield.setUncertainty2(w2Vec->fArray[0] + w2Vec->fArray[1] + w2Vec->fArray[2]);

      yields["noSyst"] = yield;

      // Todo: Add here the other systematics retrieval
      if(doSyst)
      {
      }

      channelYields[process.name] = yields;
    }

    retVal[channel.name()] = channelYields;
  }

  return retVal;
}

//TDirectory* cwd = gDirectory;
//cwd->cd();
bool DatacardMaker::genDatacards()
{
  if(jsonLoaded_ == "")
    return false;

  std::string cmd = "if [[ -d \"";
  cmd += outDir_;
  cmd += "\" ]]; then exit 1; else exit 0; fi";
  if(system(cmd.c_str()) == 0)
  {
    cmd = "mkdir ";
    cmd += outDir_;
    system(cmd.c_str());
  }
  else
  {
    // TODO: The directory already exists, check if there are signal region datacards and if they exist, prompt the user to delete them or not.
  }

  TDirectory* cwd = gDirectory;
  scratchArea_ = new TFile(".finalSelectionScratchArea.root", "RECREATE");

  for(auto &signalRegion : signalRegions_)
  {
    // Todo: Load the background info with the cuts from this signal region. Do not forget about systematics (Done? systematics will be tricky)
    auto backgrounds = applySelection("B", processes_["BG"], signalRegion);
    std::map<std::string,doubleUnc> totalBackground, data;
    {
      for(auto &channel : channels_)
      {
        doubleUnc total(0,0);
        for(auto &process : backgrounds[channel.name()])
          total += process.second["noSyst"];

        totalBackground[channel.name()] = total;
      }
      // TODO: if unblind, process data and fill the data variable from above
    }

    if(verbose_)
    {
      for(auto &channel : channels_)
      {
        std::cout << "Background yields for channel " << channel.name() << std::endl;
        for(auto &process : backgrounds[channel.name()])
        {
          std::cout << "  " << process.first << ": " << process.second["noSyst"] << std::endl;
        }
      }
    }

    std::vector<int> signalPoints = getSignalPoints(signalRegion.signalSelection());
    for(auto &point : signalPoints)
    {
      std::stringstream temp;
      std::string fileName;

      temp << outDir_ << "/SignalPoint_" << point << ".txt";
      temp >> fileName;

      {
        bool exists = std::ifstream(fileName).good();
        if(exists)
          std::cout << "DatacardMaker::genDatacards(): The defined signal regions overlap, program execution will continue but only the last signal region for a given signal point will be kept. Please check " << jsonFile_ << " for errors." << std::endl;
      }


      temp.clear();
      std::string signalPointSelection;
      temp << "((" << signalPointVariable_ << ")==" << point << ")";
      temp >> signalPointSelection;
      auto signals = applySelection("S", processes_["SIG"], signalRegion, signalPointSelection);

      ofstream file(fileName, std::ios::out | std::ios::trunc | std::ios::binary);

      std::streambuf *coutbuf = std::cout.rdbuf();
      if(!file.is_open())
        std::cout << "DatacardMaker::genDatacards(): There was a problem opening the file '" << fileName << "'. The output of the file will be redirected to the terminal." << std::endl;
      else
        std::cout.rdbuf(file.rdbuf());

      std::string separator = "";
      for(int i = 0; i < 80; ++i)
        separator += "-";

      // Output datacard to file (or stdout if opening the file failed)
      std::cout << "# Datacard for " << name_ << " analysis for signal point " << point << std::endl;
      std::cout << "imax " << channels_.size() << " number of channels" << std::endl;
      std::cout << "jmax * number of backgrounds" << std::endl;
      std::cout << "kmax * number of nuisance parameters" << std::endl;
      std::cout << separator << std::endl;

      std::cout << "bin";
      for(auto &channel : channels_)
      {
        std::cout << " " << channel.name();
      }
      std::cout << std::endl;
      // Todo: Change the following output if we are unblinded
      std::cout << "# We are blinded right now so the following observation numbers are the sum of the MC processes" << std::endl;
      std::cout << "observation";
      for(auto &channel : channels_)
      {
        // TODO: Add here to unblind
        std::cout << " " << int(totalBackground[channel.name()].value());
      }
      std::cout << std::endl;
      std::cout << separator << std::endl;

      std::cout << "bin";
      for(auto &channel : channels_)
      {
        for(auto &process : processes_["SIG"])
          std::cout << " " << channel.name();
        for(auto &process : processes_["BG"])
          std::cout << " " << channel.name();
      }
      std::cout << std::endl;

      std::cout << "process";
      for(auto &channel : channels_)
      {
        for(auto &process : processes_["SIG"])
          std::cout << " " << process.label;
        for(auto &process : processes_["BG"])
          std::cout << " " << process.label;
      }
      std::cout << std::endl;

      std::cout << "process";
      for(auto &channel : channels_)
      {
        int index = 0;
        for(auto &process : processes_["SIG"])
        {
          std::cout << " " << index;
          --index;
        }
        index = 1;
        for(auto &process : processes_["BG"])
        {
          std::cout << " " << index;
          ++index;
        }
      }
      std::cout << std::endl;

      std::cout << "rate";
      for(auto &channel : channels_)
      {
        for(auto &process : processes_["SIG"])
          std::cout << " " << signals[channel.name()][process.name]["noSyst"].value();
        for(auto &process : processes_["BG"])
          std::cout << " " << backgrounds[channel.name()][process.name]["noSyst"].value();
      }
      std::cout << std::endl;
      std::cout << separator << std::endl;

      // TODO: Fixme, do not use a hardcoded value, get it from the json
      std::cout << "lumi_8TeV lnN";
      for(auto &channel : channels_)
      {
        for(auto &process : processes_["SIG"])
          std::cout << " 1.026";
        for(auto &process : processes_["BG"])
          std::cout << " 1.026";
      }
      std::cout << std::endl;

      // Statistical uncertainty
      for(auto &channel : channels_)
      {
        for(auto &process : processes_["SIG"])
        {
          if(signals[channel.name()][process.name]["noSyst"].value() == 0)
            continue;

          std::cout << "CMS_" << name_ << "_" << process.label << "_" << channel.name() << "_stat lnN";

          for(auto &channel2 : channels_)
          {
            for(auto &process2 : processes_["SIG"])
            {
              if(process.label != process2.label || channel.name() != channel2.name())
                std::cout << " -";
              else
                std::cout << " " << (signals[channel.name()][process.name]["noSyst"].value() + signals[channel.name()][process.name]["noSyst"].uncertainty())/(signals[channel.name()][process.name]["noSyst"].value());
            }
            for(auto &process2 : processes_["BG"])
              std::cout << " -";
          }

          std::cout << std::endl;
        }

        for(auto &process : processes_["BG"])
        {
          if(backgrounds[channel.name()][process.name]["noSyst"].value() == 0)
            continue;

          std::cout << "CMS_" << name_ << "_" << process.label << "_" << channel.name() << "_stat lnN";

          for(auto &channel2 : channels_)
          {
            for(auto &process2 : processes_["SIG"])
              std::cout << " -";
            for(auto &process2 : processes_["BG"])
            {
              if(process.label != process2.label || channel.name() != channel2.name())
                std::cout << " -";
              else
                std::cout << " " << (backgrounds[channel.name()][process.name]["noSyst"].value() + backgrounds[channel.name()][process.name]["noSyst"].uncertainty())/(backgrounds[channel.name()][process.name]["noSyst"].value());
            }
          }

          std::cout << std::endl;
        }
      }
      // TODO: Finish outputting uncertainties

      if(file.is_open())
      {
        std::cout.rdbuf(coutbuf);
        file.close();
      }
    }
  }

  cwd->cd();
  scratchArea_->Close();
  std::remove(".finalSelectionScratchArea.root");

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

ChannelInfo::ChannelInfo(std::string name, std::string selection):
  isValid_(true),
  name_(name),
  selection_(selection)
{
  if(name == "" || selection == "")
    isValid_ = false;
}

bool ChannelInfo::loadJson(JSONWrapper::Object& json)
{
  name_ = json.getString("name", "");
  selection_ = json.getString("selection", "");

  if(name_ == "" || selection_ == "")
    return false;

  return true;
}

std::string SignalRegion::cuts() const
{
  std::string retVal = "(";
  bool first = true;

  for(auto &cut : cuts_)
  {
    if(cut.isValid())
    {
      if(!first)
        retVal += "&&";
      else
        first = false;
      retVal += "(" + cut.cut() + ")";
    }
  }

  if(retVal == "(")
    retVal = "";
  else
    retVal += ")";

  return retVal;
}

std::string SignalRegionCutInfo::cut() const
{
  std::string retVal;

  std::stringstream temp;
  temp << variable_;
  if(property_ != "")
    temp << "." << property_;
  temp << direction_ << value_;
  temp >> retVal;

  return retVal;
}
