// <author>Cristóvão B da Cruz e Silva</author>
// <email>c.beirao@cern.ch</email>
// <date>2014-06-02</date>
// <summary>Script that runs an iterative cut optimization procedure. Cut values are scanned and the value with the maximum FOM is used.</summary>

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <cassert>
#include <map>
#include <algorithm>
#include <cmath>
#include <fstream>

#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TTree.h"
#include "TCut.h"
#include "TEventList.h"
#include "TInterpreter.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "THStack.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include "UserCode/llvv_fwk/interface/JSONWrapper.h"
#include "UserCode/llvv_fwk/interface/llvvObjects.h"

class OptimizationRoundInfo;

void printHelp();

class MyStyle
{
public:
  MyStyle()
  {
    marker_ = 1;
    lcolor_ = kBlack;
    mcolor_ = kBlack;
    fcolor_ = kWhite;
    lwidth_ = 1;
    lstyle_ = 1;
  };

  inline int& marker(){return marker_;};
  inline int& lcolor(){return lcolor_;};
  inline int& mcolor(){return mcolor_;};
  inline int& fcolor(){return fcolor_;};
  inline int& lwidth(){return lwidth_;};
  inline int& lstyle(){return lstyle_;};

private:
  int marker_;
  int lcolor_;
  int mcolor_;
  int fcolor_;
  int lwidth_;
  int lstyle_;

protected:
};

struct doubleUnc
{
  double value;
  double uncertainty;
};

std::ostream& operator << (std::ostream &o, doubleUnc& val)
{
  return o << val.value << " +- " << val.uncertainty;
}

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
  std::string direction;
  double value;
};

struct CutInfo
{
  std::string variable;
  std::string direction;
  double value;
  std::map<std::string,std::map<std::string,doubleUnc>> yield;
  double FOM;
  double FOMerr;
};

struct ReportInfo
{
  int round;
  std::vector<CutInfo> cuts;
};

class CutOptimizer
{
public:
  CutOptimizer(const std::string & jsonFile, const std::string & outDir, const std::vector<std::string> & plotExt);
  ~CutOptimizer();

  bool LoadJson();
  bool OptimizeAllRounds();
  bool OptimizeRound(size_t n);

  inline std::string getJson() const {return jsonFile_;};
  std::string& setJson(const std::string & jsonFile);
  inline std::string getOutDir() const {return outDir_;};
  std::string& setOutDir(const std::string & outDir);
  inline std::vector<std::string> getPlotExt() const {return plotExt_;};
  std::vector<std::string>& setPlotExt(const std::vector<std::string> & plotExt);
  inline size_t getNRounds() const {return roundInfo_.size();};
  inline void setVerbose() {verbose_ = true; return;};

private:
  std::string jsonFile_;
  std::string outDir_;
  std::vector<std::string> plotExt_;

  bool jsonLoaded;
  bool verbose_;

  bool isMultipointSignalSample_;
  int  nSignalPoints_;

  std::map<std::string,bool> FileExists_;
  std::vector<ReportInfo> report_;
//  JSONWrapper::Object* json;

  std::vector<OptimizationRoundInfo> roundInfo_;
  std::map<std::string,std::vector<ProcessFiles>> processes_;

  bool OptimizeRound_(size_t n);
  bool GetSamples(size_t n);
  int  GetSigPoints(size_t n);
  CutInfo GetBestCutAndMakePlots(size_t n, ReportInfo& report);
  std::vector<doubleUnc> GetFOM(std::map<std::string,std::map<std::string,doubleUnc>>& yield);
  std::map<std::string,std::map<std::string,doubleUnc>> GetYields(ReportInfo& report, TCut signalSelection, TCut currentSelection, double integratedLuminosity);
  void SaveGraph(std::string& name, std::vector<double>& xVals, std::vector<double>& xValsUnc, std::string& xTitle, std::vector<double>& yVals, std::vector<double>& yValsUnc, std::string& yTitle);
  std::map<std::string,std::map<std::string,TH1D*>> GetAndSaveHists(ReportInfo& report, TCut signalSelection, TCut currentSelection, double integratedLuminosity, std::string varName, std::string varExpression, std::string varLabel, std::map<std::string,double>& varParams);
  bool ClearSamples();

protected:
};

class OptimizationVariableInfo
{
public:
  OptimizationVariableInfo();
  ~OptimizationVariableInfo();

  // If using pointers the following lines should be uncommented and the respective functions defined
  // More info: http://www.cplusplus.com/articles/y8hv0pDG/
  //OptimizationVariableInfo(const OptimizationVariableInfo& other);
  //OptimizationVariableInfo& operator=(const OptimizationVariableInfo& rhs);

  inline std::string& name(){return _name;};
  inline std::string& expression(){return _expression;};
  inline double& minVal(){return _minVal;};
  inline double& maxVal(){return _maxVal;};
  inline double& step(){return _step;};
  inline double& bins(){return _bins;};
  inline std::string& label(){return _label;};

  friend bool CutOptimizer::LoadJson();

private:
  std::string _name;
  std::string _expression;
  double _minVal, _maxVal, _step, _bins;
  std::string _label;

protected:
};

class OptimizationRoundInfo
{
public:
  OptimizationRoundInfo();
  ~OptimizationRoundInfo();

  // If using pointers the following lines should be uncommented and the respective functions defined
  // More info: http://www.cplusplus.com/articles/y8hv0pDG/
  //OptimizationRoundInfo(const OptimizationRoundInfo& other);
  //OptimizationRoundInfo& operator=(const OptimizationRoundInfo& rhs);

  inline std::string& name(){return _name;};
  inline std::string& ttree(){return _ttree;};
  inline std::string& customExtension(){return _customExtension;};
  inline std::string& baseSelection(){return _baseSelection;};
  inline std::string& signalSelection(){return _signalSelection;};
  inline std::string& channel(){return _channel;};
  inline double& iLumi(){return _iLumi;};
  inline std::string& inDir(){return _inDir;};
  inline std::string& jsonFile(){return _jsonFile;};
  inline size_t nVars(){return _variables.size();};
  inline std::string pointVariable(){return _pointVariable;};
  inline bool isSelected(int pass){return (static_cast<size_t>(pass) < _UserCuts.size());};
  std::string getUserCut(int pass);
  inline std::string getUserCutVar(int pass){if(isSelected(pass)) return _UserCuts[pass].variable; else return "";};
  inline std::string getUserCutDir(int pass){if(isSelected(pass)) return _UserCuts[pass].direction; else return "";};
  inline double getUserCutValue(int pass){if(isSelected(pass)) return _UserCuts[pass].value; else return -9999;};
  inline UserCutInfo getUserCutInfo(int pass){if(isSelected(pass)) return _UserCuts[pass]; else return UserCutInfo{"","",-9999}; };
  inline std::string signalPoint(){return _signalPoint;};
  inline double sigCrossSection(){return _sigCrossSection;};
  inline int nInitEvents(){return _nInitEvents;};

  std::vector<std::string> getListOfVariables();
  std::map<std::string,std::string> getVariableExpressions();
  std::map<std::string,std::map<std::string,double>> getVariableParameterMap();
  std::map<std::string,std::string> getVariableLabels();

  friend bool CutOptimizer::LoadJson();

private:
  static size_t _counter;
  std::string _name;
  std::string _ttree;
  std::string _customExtension;
  std::string _baseSelection;
  std::string _signalSelection;
  std::string _signalPoint;
  std::string _channel;
  std::vector<OptimizationVariableInfo> _variables;
  double      _iLumi;
  std::string _inDir;
  std::string _jsonFile;
  std::string _pointVariable;
  std::vector<UserCutInfo> _UserCuts;
  double      _sigCrossSection;
  int         _nInitEvents;

protected:
};
size_t OptimizationRoundInfo::_counter = 0;


//
// Output codes:
// 0 - Everything has run fine
// 1 - Invalid arguments
// 2 - Problem parsing the arguments
//
int main(int argc, char** argv)
{
  AutoLibraryLoader::enable();

  std::string jsonFile;
  std::string outDir = "./OUT/";
  std::vector<std::string> plotExt;
  bool verbose = false;

//  gInterpreter->GenerateDictionary("llvvMet", "");
//  gROOT->ProcessLine("#include \"UserCode/llvv_fwk/interface/llvvObjects.h\"");

  // Parse the command line options
  for(int i = 1; i < argc; ++i)
  {
    std::string arg(argv[i]);

    // First see if help info is being requested, if so, print help and then exit
    if(arg.find("--help") != std::string::npos)
    {
      printHelp();
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

    if(arg.find("--plotExt") != std::string::npos)
    {
      plotExt.push_back(argv[i+1]);
      ++i;
    }

    if(arg.find("--verbose") != std::string::npos)
    {
      verbose = true;
    }
  }
  if(plotExt.size() == 0)
    plotExt.push_back(".png");

  if(jsonFile == "")
  {
    std::cout << "You should define at least the following arguments: json" << std::endl;
    std::cout << "For more information, consult the help (\"runCutOptimizer --help\")" << std::endl;
    return 2;
  }

  if(verbose)
    std::cout << "Creating an object of the CutOptimizer class..." << std::endl;
  CutOptimizer myOptimizer(jsonFile, outDir, plotExt);
  if(verbose)
    myOptimizer.setVerbose();
  myOptimizer.LoadJson();
  myOptimizer.OptimizeAllRounds();
//  myOptimizer.OptimizeRound(2);

  return 0;
}

OptimizationRoundInfo::OptimizationRoundInfo()
{
  std::stringstream buf;
  buf << "Round_" << _counter;
  buf >> _name;
  _ttree           = "Events";
  _customExtension = "";
  _baseSelection   = "selected";
  _signalSelection = "";
  _channel         = "";
  _iLumi           = 1;
  _inDir           = "";
  _jsonFile        = "";
  ++_counter;
}

OptimizationRoundInfo::~OptimizationRoundInfo()
{
  _variables.clear();
}

std::string OptimizationRoundInfo::getUserCut(int pass)
{
  if(!isSelected(pass))
    return "";

  std::string retVal;
  bool found = false;

  for(auto var = _variables.begin(); var != _variables.end(); ++var)
  {
    if(var->name() == _UserCuts[pass].variable)
    {
      found = true;
      if(_UserCuts[pass].direction != "above" && _UserCuts[pass].direction != "below")
        break;
      std::stringstream buf;
      buf << _UserCuts[pass].value;
      buf >> retVal;
      if(_UserCuts[pass].direction == "above")
        retVal = ">" + retVal;
      else
        retVal = "<" + retVal;
      retVal = var->expression() + retVal;
      break;
    }
  }

  if(!found)
    return "";
  return retVal;
}

std::vector<std::string> OptimizationRoundInfo::getListOfVariables()
{
  std::vector<std::string> retVal;

  for(auto variable = _variables.begin(); variable != _variables.end(); ++variable)
    retVal.push_back(variable->name());

  return retVal;
}

std::map<std::string,std::string> OptimizationRoundInfo::getVariableLabels()
{
  std::map<std::string,std::string> retVal;

  for(auto variable = _variables.begin(); variable != _variables.end(); ++variable)
  {
    retVal[variable->name()] = variable->label();
    if(retVal[variable->name()] == "")
      retVal[variable->name()] = variable->name();
  }

  return retVal;
}

std::map<std::string,std::map<std::string,double>> OptimizationRoundInfo::getVariableParameterMap()
{
  std::map<std::string,std::map<std::string,double>> retVal;

  for(auto variable = _variables.begin(); variable != _variables.end(); ++variable)
  {
    std::map<std::string,double> tempVal;

    tempVal["minVal"] = variable->minVal();
    tempVal["maxVal"] = variable->maxVal();
    tempVal["step"]   = variable->step();
    tempVal["bins"]   = variable->bins();

    retVal[variable->name()] = tempVal;
  }

  return retVal;
}

std::map<std::string,std::string> OptimizationRoundInfo::getVariableExpressions()
{
  std::map<std::string,std::string> retVal;

  for(auto variable = _variables.begin(); variable != _variables.end(); ++variable)
  {
    retVal[variable->name()] = variable->expression();
  }

  return retVal;
}

OptimizationVariableInfo::OptimizationVariableInfo()
{
  _name = "";
  _expression = "";
  _minVal = 0;
  _maxVal = 0;
  _step = 0;
  _bins = 1;
}

OptimizationVariableInfo::~OptimizationVariableInfo()
{
}

CutOptimizer::CutOptimizer(const std::string & jsonFile, const std::string & outDir, const std::vector<std::string> & plotExt):
  jsonFile_(jsonFile), outDir_(outDir), plotExt_(plotExt), jsonLoaded(false), verbose_(false)
{
//  json = new JSONWrapper::Object(jsonFile_, true);
}

CutOptimizer::~CutOptimizer()
{
  std::cout << "The list of ignored files, either missing or corrupt, can be found below:" << std::endl;
  for(auto key = FileExists_.begin(); key != FileExists_.end(); ++key)
  {
    if(!key->second)
      std::cout << "  " << key->first << std::endl;
  }
}

bool CutOptimizer::LoadJson()
{
  if(jsonLoaded)
    return true;

  if(verbose_)
    std::cout << "Loading JSON" << std::endl;

  roundInfo_.clear();
  JSONWrapper::Object json(jsonFile_, true);

  std::vector<JSONWrapper::Object> rounds = json["optim"].daughters();
  if(verbose_)
    std::cout << "  Found " << rounds.size() << " rounds of optimization to perform." << std::endl;
  for(auto round = rounds.begin(); round != rounds.end(); ++round)
  {
    OptimizationRoundInfo roundInfo;

    std::string name = round->getString("name", "");
    if(name != "")
      roundInfo._name = name;

    if(verbose_)
      std::cout << "  Round: " << roundInfo._name << std::endl;

    roundInfo._iLumi = round->getDouble("iLumi", 0);
    if(roundInfo._iLumi <= 0)
    {
      std::cout << roundInfo._name << ": Integrated luminosity should be positive and non-zero. Continuing..." << std::endl;
      continue;
    }
    roundInfo._sigCrossSection = round->getDouble("sigCrossSection", 0);
    if(roundInfo._sigCrossSection <= 0)
    {
      std::cout << roundInfo._name << ": Signal cross section should be positive and non-zero. Continuing..." << std::endl;
      continue;
    }
    roundInfo._inDir = round->getString("inDir", "");
    if(roundInfo._inDir == "")
    {
      std::cout << roundInfo._name << ": Input directory should be defined in the json for the optimization round. Continuing..." << std::endl;
      continue;
    }
    roundInfo._jsonFile = round->getString("jsonFile", "");
    if(roundInfo._jsonFile == "")
    {
      std::cout << roundInfo._name << ": JSON file must be specified for optimization round in cut optimization JSON. Continuing..." << std::endl;
      continue;
    }

    roundInfo._ttree = round->getString("ttree", roundInfo._ttree);
    roundInfo._customExtension = round->getString("customExtension", roundInfo._customExtension);
    roundInfo._baseSelection = round->getString("baseSelection", roundInfo._baseSelection);
    roundInfo._signalSelection = round->getString("signalSelection", roundInfo._signalSelection);
    roundInfo._signalPoint = round->getString("signalPoint", roundInfo._signalPoint);
    roundInfo._channel = round->getString("channel", roundInfo._channel);
    roundInfo._pointVariable = round->getString("pointVariable", roundInfo._pointVariable);
    roundInfo._nInitEvents = round->getInt("nInitEvents", roundInfo._nInitEvents);

    if(verbose_)
      std::cout << "    Loaded round info, now trying to load user-defined cuts, if any." << std::endl;

    auto userCuts = (*round)["userCuts"].daughters();
    for(auto userCut = userCuts.begin(); userCut != userCuts.end(); ++userCut)
    {
      UserCutInfo temp;
      temp.variable = userCut->getString("variable", "");
      temp.direction = userCut->getString("direction", "");
      temp.value = userCut->getDouble("value", -9999);
      if(temp.variable != "" && (temp.direction == "below" || temp.direction == "above"))
        roundInfo._UserCuts.push_back(temp);
    }

    if(verbose_)
    {
      std::cout << "    Loaded " << roundInfo._UserCuts.size() << " user-defined cuts." << std::endl;
      std::cout << "    Loading variable information." << std::endl;
    }

    auto variables = (*round)["variables"].daughters();
    for(auto variable = variables.begin(); variable != variables.end(); ++variable)
    {
      OptimizationVariableInfo variableInfo;

      variableInfo._name = variable->getString("name", "");
      if(variableInfo._name == "")
      {
        std::cout << roundInfo._name << ": All variables must have names. Continuing..." << std::endl;
        continue;
      }
      variableInfo._expression = variable->getString("expression", "");
      if(variableInfo._expression == "")
      {
        std::cout << roundInfo._name << "::" << variableInfo._name << ": This variable must have an expression, it must be a valid root expression. Continuing..." << std::endl;
        continue;
      }
      variableInfo._minVal = variable->getDouble("minVal", 0);
      variableInfo._maxVal = variable->getDouble("maxVal", 0);
      if(variableInfo._maxVal - variableInfo._minVal <= 0)
      {
        std::cout << roundInfo._name << "::" << variableInfo._name << ": maxVal and minVal must be specified and define a valid range of values. Continuing..." << std::endl;
        continue;
      }
      variableInfo._step = variable->getDouble("bins", 1);
      if(variableInfo._step <= 0)
      {
        std::cout << roundInfo._name << "::" << variableInfo._name << ": bins must be a valid (positive and non-zero) value. Continuing..." << std::endl;
        continue;
      }
      variableInfo._step = variable->getDouble("step", 0);
      if(variableInfo._step <= 0 || variableInfo._step > variableInfo._maxVal - variableInfo._minVal)
      {
        std::cout << roundInfo._name << "::" << variableInfo._name << ": step must be a resonable and valid value. Continuing..." << std::endl;
        continue;
      }
      variableInfo._label = variable->getString("label", "");

      roundInfo._variables.push_back(variableInfo);
    }

    if(verbose_)
    {
      std::cout << "    Loaded " << roundInfo._variables.size() << " variables." << std::endl;
      std::cout << "    Saving round informantion" << std::endl;
    }

    roundInfo_.push_back(roundInfo);
  }

  if(verbose_)
    std::cout << "Finished loading information from JSON file" << std::endl;

  if(roundInfo_.size() != 0)
    jsonLoaded = true;
  else
    std::cout << "Unable to load JSON, it might be empty/malformed or not present. Please check and try again." << std::endl;
  return jsonLoaded;
}

std::string& CutOptimizer::setJson(const std::string & jsonFile)
{
  if(verbose_)
    std::cout << "Changing json file to: " << jsonFile << std::endl;

  if(jsonFile_ != jsonFile)
  {
    jsonLoaded = false;
    ClearSamples();
  }
  jsonFile_ = jsonFile;
  return jsonFile_;
}

std::string& CutOptimizer::setOutDir(const std::string & outDir)
{
  if(verbose_)
    std::cout << "Changing output directory to: " << outDir << std::endl;

  outDir_ = outDir;
  return outDir_;
}

std::vector<std::string>& CutOptimizer::setPlotExt(const std::vector<std::string> & plotExt)
{
  plotExt_ = plotExt;
  return plotExt_;
}

bool CutOptimizer::OptimizeAllRounds()
{
  bool retVal = true;

  if(verbose_)
    std::cout << "Running optimizer on all rounds" << std::endl;

  if(!jsonLoaded)
  {
    std::cout << "The JSON file has not yet been loaded, attempting to load it" << std::endl;
    LoadJson();
    if(!jsonLoaded)
      return false;
  }

  for(size_t i = 0; i < roundInfo_.size(); ++i)
  {
    if(!OptimizeRound_(i))
      retVal = false;
  }

  return retVal;
}

bool CutOptimizer::OptimizeRound(size_t n)
{
  if(!jsonLoaded)
  {
    std::cout << "The JSON file has not yet been loaded, attempting to load it" << std::endl;
    LoadJson();
    if(!jsonLoaded)
      return false;
  }

  if(n > roundInfo_.size())
    return false;

  return OptimizeRound_(n);
}

bool CutOptimizer::OptimizeRound_(size_t n)
{
  bool retVal = false;

  if(verbose_)
    std::cout << "Optimizing round " << roundInfo_[n].name() << std::endl;

  if(roundInfo_[n].nVars() == 0)
  {
    if(verbose_)
      std::cout << "  There are no variables to optimize cuts on, skipping this round." << std::endl;
    return retVal;
  }

  system(("mkdir -p " + outDir_).c_str());
  std::string roundName = roundInfo_[n].name();
  if(roundName == "")
  {
    std::stringstream stream;
    stream << "Round_" << n;
    stream >> roundName;
  }

  std::streambuf *coutbuf = std::cout.rdbuf();
  std::streambuf *cerrbuf = std::cerr.rdbuf();
  if(verbose_)
    std::cout << "Opening log file: " << outDir_+"/"+roundName << std::endl;
  ofstream out(outDir_+"/"+roundName+".log", std::ios::out | std::ios::trunc);
  out << "This is the log file for round: " << roundName << std::endl;
  out.flush();
  if(verbose_)
    std::cout << "Redirecting standard output and error output to log file (if you still want to see the output, enable the tee option)" << std::endl;
  std::cout.rdbuf(out.rdbuf()); // Todo - enable tee
  std::cerr.rdbuf(out.rdbuf());

  std::cout << "Processing round: " << roundInfo_[n].name() << std::endl;
  std::cout << "\tReading from " << roundInfo_[n].jsonFile() << " and taking samples from " << roundInfo_[n].inDir() << " directory." << std::endl;
  std::cout << "\tUsing an integrated luminosity of " << roundInfo_[n].iLumi() << "." << std::endl;
  std::cout << "\tReading from ttree: " << roundInfo_[n].ttree();
  if(roundInfo_[n].baseSelection() != "")
    std::cout << ", with a base selection of \"" << roundInfo_[n].baseSelection() << "\"";
  if(roundInfo_[n].channel() != "")
    std::cout << " and performing cut optimization on the channel " << roundInfo_[n].channel();
  std::cout << "." << std::endl;
  std::cout << "\tThere are " << roundInfo_[n].nVars() << " variables to perform cut optimization on." << std::endl;

  GetSamples(n);

  std::cout << "\tFound " << processes_["BG"].size()  << " background processes:" << std::endl;
  for(auto process = processes_["BG"].begin(); process != processes_["BG"].end(); ++process)
  {
    std::cout << "\t  " << process->name << ":" << std::endl;
    for(auto sample = process->samples.begin(); sample != process->samples.end(); ++sample)
      std::cout << "\t    " << sample->chain->GetTitle() << " with " << sample->chain->GetEntries() << " entries in " << sample->nFiles << " files" << std::endl;
  }
  std::cout << "\tFound " << processes_["SIG"].size()  << " signal processes:" << std::endl;
  isMultipointSignalSample_ = false;
  for(auto process = processes_["SIG"].begin(); process != processes_["SIG"].end(); ++process)
  {
    std::cout << "\t  " << process->name << ":" << std::endl;
    if(process->name.find("TStauStau") != std::string::npos)
      isMultipointSignalSample_ = true;
    for(auto sample = process->samples.begin(); sample != process->samples.end(); ++sample)
      std::cout << "\t    " << sample->chain->GetTitle() << " with " << sample->chain->GetEntries() << " entries in " << sample->nFiles << " files" << std::endl;
  }

  if(processes_["BG"].size() != 0 && processes_["SIG"].size() != 0)
  {
    if(isMultipointSignalSample_)
    {
      nSignalPoints_ = GetSigPoints(n);
      std::cout << "Running on " << nSignalPoints_ << " signal points." << std::endl;
    }
    else
    {
      nSignalPoints_ = 1;
    }
    bool improve = true;

    ReportInfo myReport;
    myReport.round = n;

    CutInfo previousCut;
    previousCut.FOM = 0;
    previousCut.FOMerr = 0;

    while(improve && nSignalPoints_ != 0)
    {
      std::cout << roundInfo_[n].name() << ": Starting pass " << myReport.cuts.size() << std::endl;

      CutInfo thisCut = GetBestCutAndMakePlots(n, myReport);

      bool isSelected = roundInfo_[n].isSelected(myReport.cuts.size());

      if(thisCut.FOM == 0 && !isSelected)
      {
        improve = false;
      }
      else
      {
        if(thisCut.FOM - previousCut.FOM > previousCut.FOMerr || isSelected)
        {
          previousCut = thisCut;
          myReport.cuts.push_back(thisCut);
          std::cout << "Added cut on " << thisCut.variable << " with values " << thisCut.direction << " " << thisCut.value << ".";
          std::cout << "The cut has a Figure of Merit (FOM) of " << thisCut.FOM << " +- " << thisCut.FOMerr << "; with the following yields:" << std::endl;

          std::cout << "  Signal Processes:" << std::endl;
          doubleUnc sigTot{0,0};
          for(auto process = thisCut.yield["SIG"].begin(); process != thisCut.yield["SIG"].end(); ++process)
          {
            std::cout << "    " << process->first << ":" << std::endl;
            doubleUnc processTot = process->second;

            sigTot.value += processTot.value;
            sigTot.uncertainty += processTot.uncertainty*processTot.uncertainty;
            std::cout << "      Total: " << processTot << std::endl;
          }
          sigTot.uncertainty = std::sqrt(sigTot.uncertainty);
          std::cout << "    Total: " << sigTot << std::endl;

          std::cout << "  Background Processes:" << std::endl;
          doubleUnc bgTot{0,0};
          for(auto process = thisCut.yield["BG"].begin(); process != thisCut.yield["BG"].end(); ++process)
          {
            std::cout << "    " << process->first << ":" << std::endl;
            doubleUnc processTot = process->second;

            bgTot.value += processTot.value;
            bgTot.uncertainty += processTot.uncertainty*processTot.uncertainty;
            std::cout << "      Total: " << processTot << std::endl;
          }
          bgTot.uncertainty = std::sqrt(bgTot.uncertainty);
          std::cout << "    Total: " << bgTot << std::endl;
        }
        else
          improve = false;
      }
    }

    // TODO - Print the final selection info
    // TODO - Make the report (including the tex table with the yields)
  }
  else
  {
    if(processes_["BG"].size() == 0)
      std::cout << "\tIt doesn't make sense to perform cut optimization without any background samples. Skipping this round of optimization, please check the file " << roundInfo_[n].jsonFile() << "." << std::endl;
    if(processes_["SIG"].size() == 0)
      std::cout << "\tIt doesn't make sense to perform cut optimization without any signal samples. Skipping this round of optimization, please check the file " << roundInfo_[n].jsonFile() << "." << std::endl;
  }

  std::cout.rdbuf(coutbuf);
  std::cerr.rdbuf(cerrbuf);
  if(verbose_)
    std::cout << "Finished optimization round, restored standard output and standard error" << std::endl;
  out.close();

  return retVal;
}

CutInfo CutOptimizer::GetBestCutAndMakePlots(size_t n, ReportInfo& report)
{
  CutInfo retVal;
  std::vector<std::string> variables = roundInfo_[n].getListOfVariables(); // TODO - buffer this stuff
  std::map<std::string,std::string> variableExpressions = roundInfo_[n].getVariableExpressions();
  std::map<std::string,std::map<std::string,double>> variableParameterMap = roundInfo_[n].getVariableParameterMap();
  std::map<std::string,std::string> variableLabels = roundInfo_[n].getVariableLabels();

  TCut baseSelection = roundInfo_[n].baseSelection().c_str();
  TCut signalSelection = roundInfo_[n].signalSelection().c_str();
  TCut cumulativeSelection = "";

  // Clean looping vector of variables already cut on and get the cumulative selection
  for(auto pass = report.cuts.begin(); pass != report.cuts.end(); ++pass)
  {
    auto newEnd = std::remove(variables.begin(), variables.end(), pass->variable);
    variables.erase(newEnd, variables.end());

    std::string tmpStr;
    std::stringstream buf;
    buf << pass->value;
    buf >> tmpStr;
    if(pass->direction == "below")
    {
      tmpStr = "<" + tmpStr;
    }
    else
    {
      if(pass->direction == "above")
        tmpStr = ">" + tmpStr;
      else
        assert(pass->direction == "below" || pass->direction == "above");
    }
    tmpStr = variableExpressions[pass->variable] + tmpStr;
    cumulativeSelection = cumulativeSelection && tmpStr.c_str();
  }

  bool isSelected = roundInfo_[n].isSelected(report.cuts.size());

  // Loop on variables
  for(auto variableName = variables.begin(); variableName != variables.end(); ++variableName)
  {
//    double minVal = variableParameterMap[*variableName]["minVal"];
//    double maxVal = variableParameterMap[*variableName]["maxVal"];
//    double step   = variableParameterMap[*variableName]["step"]; // TODO: remove all mentions of step from the code since it is no longer used
    double bins   = variableParameterMap[*variableName]["bins"];
    std::cout << roundInfo_[n].name() << "::" << *variableName << " has started processing, with " << bins + 1 << " steps to be processed." << std::endl;

//TODO - enable and use this
    std::map<std::string,std::map<std::string,TH1D*>> hists = GetAndSaveHists(report, signalSelection, baseSelection && cumulativeSelection, roundInfo_[n].iLumi(), *variableName, variableExpressions[*variableName], variableLabels[*variableName], variableParameterMap[*variableName]);

    TH1D* exampleHist = NULL;
    std::vector<double> xValsAbove, yValsAbove, xValsUncAbove, yValsUncAbove;
    std::vector<double> signalYieldAbove, backgroundYieldAbove, signalYieldUncAbove, backgroundYieldUncAbove;
    std::vector<double> xValsBelow, yValsBelow, xValsUncBelow, yValsUncBelow;
    std::vector<double> signalYieldBelow, backgroundYieldBelow, signalYieldUncBelow, backgroundYieldUncBelow;

    // Set up yield data structure
    std::map<std::string,std::map<std::string,doubleUnc>> yieldsAbove, yieldsBelow;
    for(auto index = hists.begin(); index != hists.end(); ++index)
    {
      std::map<std::string,doubleUnc> tmp;
      for(auto process = index->second.begin(); process != index->second.end(); ++process)
      {
        doubleUnc val{0,0};
        tmp[process->first] = val;
        exampleHist = process->second;
      }
      yieldsAbove[index->first] = tmp;
      yieldsBelow[index->first] = tmp;
    }

    int nbins = exampleHist->GetNbinsX();//int(bins);//+2; //TODO: remember to fix this when adding under/overflow bins

    for(int cutBin = 1; cutBin <= nbins; ++cutBin)
    {
      double binLow = exampleHist->GetBinLowEdge(cutBin);
      double binHigh = exampleHist->GetBinLowEdge(cutBin+1);
//      double binWidth = exampleHist->GetBinWidth(cutBin);

      for(auto index = hists.begin(); index != hists.end(); ++index)
      {
        for(auto process = index->second.begin(); process != index->second.end(); ++process)
        {
          yieldsAbove[index->first][process->first].value       = 0;
          yieldsAbove[index->first][process->first].uncertainty = 0;
          yieldsBelow[index->first][process->first].value       = 0;
          yieldsBelow[index->first][process->first].uncertainty = 0;

          for(int bin = 1; bin <= nbins; ++bin)
          {
            if(bin < cutBin)
            {
              yieldsBelow[index->first][process->first].value += process->second->GetBinContent(bin);
              yieldsBelow[index->first][process->first].uncertainty += process->second->GetBinError(bin)*process->second->GetBinError(bin);
            }
            if(bin > cutBin)
            {
              yieldsAbove[index->first][process->first].value += process->second->GetBinContent(bin);
              yieldsAbove[index->first][process->first].uncertainty += process->second->GetBinError(bin)*process->second->GetBinError(bin);
            }
            if(bin == cutBin)
            {
              yieldsAbove[index->first][process->first].value += process->second->GetBinContent(bin);
              yieldsAbove[index->first][process->first].uncertainty += process->second->GetBinError(bin)*process->second->GetBinError(bin);
              yieldsBelow[index->first][process->first].value += process->second->GetBinContent(bin);
              yieldsBelow[index->first][process->first].uncertainty += process->second->GetBinError(bin)*process->second->GetBinError(bin);
            }
          }

          yieldsAbove[index->first][process->first].uncertainty = std::sqrt(yieldsAbove[index->first][process->first].uncertainty);
          yieldsBelow[index->first][process->first].uncertainty = std::sqrt(yieldsBelow[index->first][process->first].uncertainty);
        }
      }

      std::vector<doubleUnc> fomReportAbove = GetFOM(yieldsAbove);
      std::vector<doubleUnc> fomReportBelow = GetFOM(yieldsBelow);

      std::cout << "  The Cut: " << variableExpressions[*variableName] << " -  > " << binLow << ":  < " << binHigh << std::endl;
      std::cout << "    n(Sig): " << fomReportAbove[1] << ": " << fomReportBelow[1] << std::endl;
      std::cout << "    n(Bkg): " << fomReportAbove[2] << ": " << fomReportBelow[2] << std::endl;
      std::cout << "    FOM: " << fomReportAbove[0] << ": " << fomReportBelow[0] << std::endl;

      if(fomReportAbove[2].value != 0 && fomReportAbove[1].value != 0)
      {
        xValsAbove.push_back(binLow);
        xValsUncAbove.push_back(0);
        yValsAbove.push_back(fomReportAbove[0].value);
        yValsUncAbove.push_back(fomReportAbove[0].uncertainty);
        signalYieldAbove.push_back(fomReportAbove[1].value);
        signalYieldUncAbove.push_back(fomReportAbove[1].uncertainty);
        backgroundYieldAbove.push_back(fomReportAbove[2].value);
        backgroundYieldUncAbove.push_back(fomReportAbove[2].uncertainty);

        if(retVal.FOM < fomReportAbove[0].value)// todo: implement here a minimum limit for the number of signal events, implement an option in the json to specify this amount
        {
          retVal.FOM = fomReportAbove[0].value;
          retVal.FOMerr = fomReportAbove[0].uncertainty;
          retVal.variable = *variableName;
          retVal.direction = "above";
          retVal.value = binLow;
          retVal.yield = yieldsAbove;
        }
      }

      if(fomReportBelow[2].value != 0 && fomReportBelow[1].value != 0)
      {
        xValsBelow.push_back(binHigh);
        xValsUncBelow.push_back(0);
        yValsBelow.push_back(fomReportBelow[0].value);
        yValsUncBelow.push_back(fomReportBelow[0].uncertainty);
        signalYieldBelow.push_back(fomReportBelow[1].value);
        signalYieldUncBelow.push_back(fomReportBelow[1].uncertainty);
        backgroundYieldBelow.push_back(fomReportBelow[2].value);
        backgroundYieldUncBelow.push_back(fomReportBelow[2].uncertainty);

        if(retVal.FOM < fomReportBelow[0].value)// todo: implement here a minimum limit for the number of signal events, implement an option in the json to specify this amount
        {
          retVal.FOM = fomReportBelow[0].value;
          retVal.FOMerr = fomReportBelow[0].uncertainty;
          retVal.variable = *variableName;
          retVal.direction = "below";
          retVal.value = binHigh;
          retVal.yield = yieldsBelow;
        }
      }
    }

    if(xValsAbove.size() == 0)
    {
      std::cout << "  No points processed for \"above\"" << std::endl;
    }
    else
    {
      std::string xTitle = variableLabels[*variableName];
      std::string yTitle = "FOM";
      std::string name;
      std::stringstream buf;
      buf << outDir_ << "/";
      buf << roundInfo_[n].name() << "_Pass" << report.cuts.size() << "_" << *variableName;
      buf >> name;
      name = name + "_above";

      SaveGraph(name, xValsAbove, xValsUncAbove, xTitle, yValsAbove, yValsUncAbove, yTitle);
      yTitle = "SignalYield";
      std::string tmp = name+"_signalYield";
      SaveGraph(tmp, xValsAbove, xValsUncAbove, xTitle, signalYieldAbove, signalYieldUncAbove, yTitle);
      yTitle = "BackgroundYield";
      tmp = name+"_backgroundYield";
      SaveGraph(tmp, xValsAbove, xValsUncAbove, xTitle, backgroundYieldAbove, backgroundYieldUncAbove, yTitle);
    }

    if(xValsBelow.size() == 0)
    {
      std::cout << "  No points processed for \"below\"" << std::endl;
    }
    else
    {
      std::string xTitle = variableLabels[*variableName];
      std::string yTitle = "FOM";
      std::string name;
      std::stringstream buf;
      buf << outDir_ << "/";
      buf << roundInfo_[n].name() << "_Pass" << report.cuts.size() << "_" << *variableName;
      buf >> name;
      name = name + "_below";

      SaveGraph(name, xValsBelow, xValsUncBelow, xTitle, yValsBelow, yValsUncBelow, yTitle);
      yTitle = "SignalYield";
      std::string tmp = name+"_signalYield";
      SaveGraph(tmp, xValsBelow, xValsUncBelow, xTitle, signalYieldBelow, signalYieldUncBelow, yTitle);
      yTitle = "BackgroundYield";
      tmp = name+"_backgroundYield";
      SaveGraph(tmp, xValsBelow, xValsUncBelow, xTitle, backgroundYieldBelow, backgroundYieldUncBelow, yTitle);
    }

    //Clearing memory
    for(auto index = hists.begin(); index != hists.end(); ++index)
    {
      for(auto process = index->second.begin(); process != index->second.end(); ++process)
      {
        delete process->second;
        process->second = NULL;
      }
    }
  }

  if(isSelected)
  {
    std::string thisCutStr = roundInfo_[n].getUserCut(report.cuts.size());
    std::cout << "  Found a user-defined cut for this pass: " << thisCutStr << std::endl;

    TCut thisCut = thisCutStr.c_str();
    TCut currentSelection = baseSelection && cumulativeSelection && thisCut;

    std::map<std::string,std::map<std::string,doubleUnc>> yield = GetYields(report, signalSelection, currentSelection, roundInfo_[n].iLumi());
    std::vector<doubleUnc> fomReport = GetFOM(yield);
    doubleUnc fom  = fomReport[0];
    doubleUnc nSig = fomReport[1];
    doubleUnc nBg  = fomReport[2];

    std::cout << "    n(Sig): " << nSig << std::endl;
    std::cout << "    n(Bkg): " << nBg << std::endl;
    std::cout << "    FOM: " << fom << std::endl;

    retVal.FOM = fom.value;
    retVal.FOMerr = fom.uncertainty;
    retVal.variable = roundInfo_[n].getUserCutVar(report.cuts.size());;
    retVal.direction = roundInfo_[n].getUserCutDir(report.cuts.size());;
    retVal.value = roundInfo_[n].getUserCutValue(report.cuts.size());;
    retVal.yield = yield;
  }

  return retVal;
}

std::map<std::string,std::map<std::string,TH1D*>> CutOptimizer::GetAndSaveHists(ReportInfo& report, TCut signalSelection, TCut currentSelection, double integratedLuminosity, std::string varName, std::string varExpression, std::string varLabel, std::map<std::string,double>& varParams)
{
  gROOT->cd();

  double minVal = varParams["minVal"];
  double maxVal = varParams["maxVal"];
//  double step   = varParams["step"];
  double bins   = varParams["bins"];

  std::map<std::string,TH1D*> hists; // Hists for ratio
  std::map<std::string,std::map<std::string,TH1D*>> histsPlot; // Hists for plotting //TODO - Add an option in the variable definition in the JSON for "noCut" so that the variable is not cut on but just plotted
  std::string title = ";"+varLabel+";Events";

  for(auto index = processes_.begin(); index != processes_.end(); ++index)
  {
    hists[index->first] = new TH1D((index->first).c_str(), (index->first+title).c_str(), bins, minVal, maxVal);
    hists[index->first]->Sumw2();

    TCut cut = currentSelection;
    if(index->first == "SIG" && isMultipointSignalSample_)
      cut = cut && signalSelection;

    for(auto process = index->second.begin(); process != index->second.end(); ++process)
    {
      histsPlot[index->first][process->name] = new TH1D((process->name).c_str(), (process->name + title).c_str(), bins, minVal, maxVal);
      histsPlot[index->first][process->name]->Sumw2();

      for(auto sample = process->samples.begin(); sample != process->samples.end(); ++sample)
      {
        if(sample->nFiles == 0)
          continue;

        TH1D temp_hist("temp_hist", title.c_str(), bins, minVal, maxVal);
        temp_hist.Sumw2();

        if(nSignalPoints_ > 1 && index->first == "SIG")
          sample->chain->Draw((varExpression+">>temp_hist").c_str(), cut*"puWeight", "goff");
        else
          sample->chain->Draw((varExpression+">>temp_hist").c_str(), cut*"weight", "goff");

        if(index->first != "Data")
        {
          if(!isMultipointSignalSample_ || index->first != "SIG")
          {
            temp_hist.Scale(1.0/sample->nFiles);
          }
          else
          {
            if(nSignalPoints_ > 1)
              temp_hist.Scale(roundInfo_[report.round].sigCrossSection() / (roundInfo_[report.round].nInitEvents() * nSignalPoints_));
          }

          temp_hist.Scale(integratedLuminosity);
        }

        histsPlot[index->first][process->name]->Add(&temp_hist);
      }

      histsPlot[index->first][process->name]->SetLineColor  (process->style.lcolor());
      histsPlot[index->first][process->name]->SetMarkerColor(process->style.mcolor());
      histsPlot[index->first][process->name]->SetFillColor  (process->style.fcolor());
      histsPlot[index->first][process->name]->SetLineWidth  (process->style.lwidth());
      histsPlot[index->first][process->name]->SetLineStyle  (process->style.lstyle());
      histsPlot[index->first][process->name]->SetMarkerStyle(process->style.marker());

      hists[index->first]->Add(histsPlot[index->first][process->name]);
    }
  }

  // Todo: expand histos to show underflow/overflow if needed
  // Todo: what about plotting 2D histos, only with the no-cut option

  std::stringstream buf;
  std::string baseFileName;
  buf << roundInfo_[report.round].name() + "_Pass" << report.cuts.size();
  buf << "_" << varName;
  buf >> baseFileName;
  double maximum = hists["BG"]->GetMaximum();
  for(auto hist = histsPlot["SIG"].begin(); hist != histsPlot["SIG"].end(); ++hist)
    if(hist->second->GetMaximum() > maximum)
      maximum = hist->second->GetMaximum();

  // Todo: add option to unblind - Do not forget about implementing the ratio
  TCanvas c1("c1", "c1", 800, 600);

  THStack* bgStack = new THStack((varName+"_BG").c_str(), (varName+title).c_str());
  for(auto hist = histsPlot["BG"].begin(); hist != histsPlot["BG"].end(); ++hist)
    bgStack->Add(hist->second);

  bgStack->Draw();
  bgStack->SetMaximum(maximum);
  for(auto hist = histsPlot["SIG"].begin(); hist != histsPlot["SIG"].end(); ++hist)
    hist->second->Draw("same");

  c1.BuildLegend(0.88, 0.67, 1, 1);

  for(auto ext = plotExt_.begin(); ext != plotExt_.end(); ++ext)
    c1.SaveAs((outDir_ + "/" + baseFileName + *ext).c_str());

  delete bgStack;

  // Todo: normalized plots


  // Cleaning up:
  for(auto hist = hists.begin(); hist != hists.end(); ++hist)
  {
    delete hist->second;
    hist->second = NULL;
  }

  return histsPlot;
}

std::vector<doubleUnc> CutOptimizer::GetFOM(std::map<std::string,std::map<std::string,doubleUnc>>& yield)
{
  std::vector<doubleUnc> retVal;

  double systUnc = 0.15;  // Hard coded systematic uncertainty /////////////////////////////////////////////////////////////////////////////////////
  // TODO - Add option to run with full systematic uncertainty (ie: calculate full analysis at each point) - is this feasible and reasonable?
  doubleUnc nSig{0,0};
  doubleUnc nBg{0,0};
  doubleUnc fom{0,0};

  for(auto process = yield["SIG"].begin(); process != yield["SIG"].end(); ++process)
  {
    nSig.value += process->second.value;
    nSig.uncertainty += process->second.uncertainty*process->second.uncertainty;
  }
  nSig.uncertainty = std::sqrt(nSig.uncertainty);

  for(auto process = yield["BG"].begin(); process != yield["BG"].end(); ++process)
  {
    nBg.value += process->second.value;
    nBg.uncertainty += process->second.uncertainty*process->second.uncertainty;
  }
  nBg.uncertainty = std::sqrt(nBg.uncertainty);

  if(nBg.value == 0 || nSig.value == 0)
  {
    retVal.push_back(fom);
    retVal.push_back(nSig);
    retVal.push_back(nBg);
    return retVal;
  }

  double dividend = std::sqrt(nBg.value + (systUnc*nBg.value)*(systUnc*nBg.value));
  fom.value = nSig.value/dividend;
  {
    double SIGContrib = nSig.uncertainty/dividend;
    double BGContrib = nBg.uncertainty*nSig.value*(1+2*systUnc*systUnc*nBg.value)/(2*dividend*dividend*dividend);
    fom.uncertainty = std::sqrt(SIGContrib*SIGContrib + BGContrib*BGContrib);
  }

  retVal.push_back(fom);
  retVal.push_back(nSig);
  retVal.push_back(nBg);

  return retVal;
}

void CutOptimizer::SaveGraph(std::string& name, std::vector<double>& xVals, std::vector<double>& xValsUnc, std::string& xTitle, std::vector<double>& yVals, std::vector<double>& yValsUnc, std::string& yTitle)
{
  if(xVals.size() == 0 || yVals.size() == 0)
    return;

  gROOT->cd();
  TCanvas c1("c1", "c1", 800, 600);

  if(yValsUnc.size() == 0)
  {
    TGraph Graph(xVals.size(), xVals.data(), yVals.data());
    // http://root.cern.ch/root/html/TAttMarker.html
    //FOMGraph.SetLineColor(2);
    //FOMGraph.SetLineWidth(2);
    //FOMGraph.SetMarkerColor(4);
    //FOMGraph.SetMarkerStyle(21);
    //http://root.cern.ch/root/html/TAttFill.html
    //http://root.cern.ch/root/html/TColor.html
    Graph.SetFillColor(kBlue-9);
    Graph.SetFillStyle(3354);
    Graph.SetMarkerStyle(21);

    Graph.Draw("ALP");

    Graph.GetXaxis()->SetTitle(xTitle.c_str());
    Graph.GetYaxis()->SetTitle(yTitle.c_str());
    Graph.GetXaxis()->CenterTitle();

    for(auto ext = plotExt_.begin(); ext != plotExt_.end(); ++ext)
    {
      c1.SaveAs((name + *ext).c_str());
    }
  }
  else
  {
    TGraphErrors Graph(xVals.size(), xVals.data(), yVals.data(), xValsUnc.data(), yValsUnc.data());
    // http://root.cern.ch/root/html/TAttMarker.html
    //FOMGraph.SetLineColor(2);
    //FOMGraph.SetLineWidth(2);
    //FOMGraph.SetMarkerColor(4);
    //FOMGraph.SetMarkerStyle(21);
    //http://root.cern.ch/root/html/TAttFill.html
    //http://root.cern.ch/root/html/TColor.html
    Graph.SetFillColor(kBlue-9);
    Graph.SetFillStyle(3354);

    Graph.Draw("3A");
    Graph.Draw("same LP");

    Graph.GetXaxis()->SetTitle(xTitle.c_str());
    Graph.GetYaxis()->SetTitle(yTitle.c_str());
    Graph.GetXaxis()->CenterTitle();

    for(auto ext = plotExt_.begin(); ext != plotExt_.end(); ++ext)
    {
      c1.SaveAs((name + *ext).c_str());
    }
  }

  return;
}

std::map<std::string,std::map<std::string,doubleUnc>> CutOptimizer::GetYields(ReportInfo& report, TCut signalSelection, TCut currentSelection, double integratedLuminosity)
{
  gROOT->cd();
  std::map<std::string,std::map<std::string,doubleUnc>> retVal;

  for(auto index = processes_.begin(); index != processes_.end(); ++index)
  {
    std::map<std::string,doubleUnc> typeYieldMap;
    for(auto process = index->second.begin(); process != index->second.end(); ++process)
    {
      doubleUnc processYield;
      processYield.value = 0;
      processYield.uncertainty = 0;

      for(auto sample = process->samples.begin(); sample != process->samples.end(); ++sample)
      {
        if(sample->nFiles == 0)
          continue;
        doubleUnc sampleYield;

        TH1D temp_histo("temp_histo", "temp_histo", 1, 0, 20);
        temp_histo.Sumw2();

        TCut cut = currentSelection;
        if(index->first == "SIG" && isMultipointSignalSample_)
          cut = signalSelection && currentSelection;

        if(nSignalPoints_ > 1 && index->first == "SIG")
          sample->chain->Draw("puWeight>>temp_histo", cut*"puWeight", "goff");
        else
          sample->chain->Draw("weight>>temp_histo", cut*"weight", "goff");

        sampleYield.value = temp_histo.GetBinContent(0) + temp_histo.GetBinContent(1) + temp_histo.GetBinContent(2);
        //double a = temp_histo.GetBinError(0);
        //double b = temp_histo.GetBinError(1);
        //double c = temp_histo.GetBinError(2);
        //double uncertainty = std::sqrt(a*a + b*b + c*c);
        TArrayD* w2Vec = temp_histo.GetSumw2();
        sampleYield.uncertainty = std::sqrt(w2Vec->fArray[0] + w2Vec->fArray[1] + w2Vec->fArray[2]);

        if(index->first != "Data")
        {
          if(!isMultipointSignalSample_ || index->first != "SIG")
          {
            sampleYield.value       = sampleYield.value       /sample->nFiles;
            sampleYield.uncertainty = sampleYield.uncertainty /sample->nFiles;
          }
          else
          {
            if(nSignalPoints_ > 1)
            {
              sampleYield.value       = (sampleYield.value       * roundInfo_[report.round].sigCrossSection())/(roundInfo_[report.round].nInitEvents() * nSignalPoints_);
              sampleYield.uncertainty = (sampleYield.uncertainty * roundInfo_[report.round].sigCrossSection())/(roundInfo_[report.round].nInitEvents() * nSignalPoints_);
            }
          }

          sampleYield.value       = sampleYield.value * integratedLuminosity;
          sampleYield.uncertainty = sampleYield.uncertainty * integratedLuminosity;
        }

        processYield.value += sampleYield.value;
        processYield.uncertainty += sampleYield.uncertainty*sampleYield.uncertainty;
      }
      processYield.uncertainty = std::sqrt(processYield.uncertainty);

      typeYieldMap[process->name] = processYield;
    }
    retVal[index->first] = typeYieldMap;
  }

  return retVal;
}

int CutOptimizer::GetSigPoints(size_t n)
{
  gROOT->cd();
  int retVal = 0;
  int maxVal = 501*1000+1; //Change-me if you are having problems with multipart signal samples and it only finds some of them
  TH1D finalHist("finalHist", "finalHist", maxVal, 0, maxVal);

  std::string variable = roundInfo_[n].pointVariable();
  TCut baseSelection = roundInfo_[n].baseSelection().c_str();
  TCut signalSelection = roundInfo_[n].signalSelection().c_str();
  TCut cut = baseSelection && signalSelection;

  for(auto process = processes_["SIG"].begin(); process != processes_["SIG"].end(); ++process)
  {
    for(auto sample = process->samples.begin(); sample != process->samples.end(); ++sample)
    {
      if(sample->nFiles == 0)
        continue;

      TH1D* temp_hist = new TH1D("temp_hist", "temp_hist", maxVal, 0, maxVal);
      sample->chain->Draw((variable+">>temp_hist").c_str(), cut*"weight", "goff");
      finalHist.Add(temp_hist);
      delete temp_hist;
    }
  }

  int nBins = finalHist.GetNbinsX();
  for(int i = 1; i <= nBins; ++i)
    if(finalHist.GetBinContent(i) != 0)
      ++retVal;

  return retVal;
}

bool CutOptimizer::GetSamples(size_t n)
{
  ClearSamples();
  JSONWrapper::Object json(roundInfo_[n].jsonFile(), true);
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

    // Todo add possibility of plotting multiple signals, by extending the data json to include the cuts and a label for each signal sample
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
      tempSample.chain = new TChain(roundInfo_[n].ttree().c_str(), ((*sample).getString("dtag", "") + (*sample).getString("suffix", "")).c_str());
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

        std::string fileName = roundInfo_[n].inDir() + "/" + (*sample).getString("dtag", "") + (*sample).getString("suffix", "") + segmentExt + filtExt + roundInfo_[n].customExtension() + ".root";

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

bool CutOptimizer::ClearSamples()
{
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

  std::vector<ProcessFiles> temp;
  processes_["BG"] = temp;
  processes_["SIG"] = temp;
  processes_["Data"] = temp;

  return true;
}

void printHelp()
{
  std::cout << "runCutOptimizer help - There are the following options:" << std::endl << std::endl;

  std::cout << "--help    -->  Print this help message" << std::endl;
  std::cout << "--json    -->  Configuration file for cut optimization, should define which files to run on, where they are located, the integrated luminosity and the variables to optimize the cuts on" << std::endl;
  std::cout << "--outDir  -->  Path to the directory where to output plots and tables (will be created if it doesn't exist)" << std::endl;
  std::cout << "--plotExt -->  Extension format with which to save the plots, repeat this command if multiple formats are desired" << std::endl;
  std::cout << "--verbose -->  Set to verbose mode, the current step will be printed out to screen" << std::endl;

  std::cout << std::endl << "Example command:" << std::endl << "\trunCutOptimizer --json optimization_options.json --outDir ./OUT/" << std::endl;
  return;
}
