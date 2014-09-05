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
#include <unordered_map>
#include <algorithm>
#include <cmath>

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
#include "TGraph.h"
#include "TGraphErrors.h"

#include "UserCode/llvv_fwk/interface/JSONWrapper.h"
#include "UserCode/llvv_fwk/interface/llvvObjects.h"

// Black magic code to get root to be able to read our custom classes
//#pragma link C++ class llvvMet;

// REMEMBER TO IGNORE ANY DATA DEFINED IN THE JSON, YOU SHOULD NEVER OPTIMIZE CUTS ON DATA
// IT MIGHT BE OK TO USE THE DATA IF USING A METHOD WITH A DATA DRIVEN BACKGROUND ESTIMATION (CONFIRM THIS BEFORE USING IT HERE)

class OptimizationRoundInfo;

typedef std::vector<std::pair<std::string,std::vector<std::pair<int,TChain*>>>> fileChains;

void printHelp();
fileChains getChainsFromJSON(JSONWrapper::Object& json, std::string RootDir, std::string type="BG", std::string treename="Events", std::string customExtension="_selected");
std::vector<OptimizationRoundInfo> getRoundsFromJSON(JSONWrapper::Object& json);
int getNumberOfPoints(fileChains files, std::string variable, TCut cut);
double applyCut(fileChains files,  TCut cut, bool isSig = false, bool correctFiles = true, bool normalize = false);
double getSignalScale(fileChains files,  TCut cut, bool correctFiles = true);

bool gDoneScale = false;
double gScale = 0;
bool gDoneSignalScale = false;
double gSignalScale = 0;

class CutOptimizer
{
public:
private:
  std::string jsonFile_;
  std::vector<OptimizationRoundInfo> roundInfo_;

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
  inline std::string& cutDir(){return _cutDir;};
  inline std::string& label(){return _label;};

  friend std::vector<OptimizationRoundInfo> getRoundsFromJSON(JSONWrapper::Object& json);

private:
  std::string _name;
  std::string _expression;
  double _minVal, _maxVal, _step;
  std::string _cutDir;
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
  inline bool isSelected(int pass){return (static_cast<size_t>(pass) < _selected.size());};
  inline std::string getSelection(int pass){if(isSelected(pass)) return selection(pass); else return "";};
  inline std::string signalPoint(){return _signalPoint;};
  inline double sigCrossSection(){return _sigCrossSection;};

  std::vector<std::string> getListOfVariables();
  std::unordered_map<std::string,std::string> getVariableExpressions();
  std::unordered_map<std::string,std::unordered_map<std::string,double>> getVariableParameterMap();
  std::unordered_map<std::string,std::string> getVariableLabels();

  friend std::vector<OptimizationRoundInfo> getRoundsFromJSON(JSONWrapper::Object& json);

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
  std::vector<std::string> _selected;
  double      _sigCrossSection;

  inline std::string selection(int pass){return _selected[pass];};

protected:
};
size_t OptimizationRoundInfo::_counter = 0;

struct FOMInfo
{
  double FOM;
  double err;
  std::string var;
  double cutVal;
};


std::unordered_map<std::string,bool> FileExists;


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
  }
  if(plotExt.size() == 0)
    plotExt.push_back(".png");

  if(jsonFile == "")
  {
    std::cout << "You should define at least the following arguments: json" << std::endl;
    std::cout << "For more information, consult the help (\"runCutOptimizer --help\")" << std::endl;
    return 2;
  }

  system(("mkdir -p " + outDir).c_str());

  JSONWrapper::Object json(jsonFile, true);

  std::vector<OptimizationRoundInfo> rounds = getRoundsFromJSON(json);

  for(auto round = rounds.begin(); round != rounds.end(); ++round)
  {
    if(round->nVars() == 0)
    {
      std::cout << "There are no variables to optimize cuts on, continuing." << std::endl;
      continue;
    }
    std::cout << "Processing round: " << round->name() << std::endl;
    std::cout << "\tReading from " << round->jsonFile() << " and taking samples from " << round->inDir() << " directory." << std::endl;
    std::cout << "\tUsing an integrated luminosity of " << round->iLumi() << "." << std::endl;
    std::cout << "\tReading from ttree: " << round->ttree();
    if(round->baseSelection() != "")
      std::cout << ", with a base selection of \"" << round->baseSelection() << "\"";
    if(round->channel() != "")
      std::cout << " and performing cut optimization on the channel " << round->channel();
    std::cout << "." << std::endl;
    std::cout << "\tThere are " << round->nVars() << " variables to perform cut optimization on." << std::endl;


    JSONWrapper::Object json(round->jsonFile(), true);
    auto BG_samples  = getChainsFromJSON(json, round->inDir(),  "BG", round->ttree(), round->customExtension());
    auto SIG_samples = getChainsFromJSON(json, round->inDir(), "SIG", round->ttree(), round->customExtension());
    gDoneScale = false;
    gDoneSignalScale = false;

    std::cout << "\tFound " << BG_samples.size()  << " background processes:" << std::endl;
    for(auto process = BG_samples.begin(); process != BG_samples.end(); ++process)
    {
      std::cout << "\t  " << process->first << ":" << std::endl;
      for(auto sample = process->second.begin(); sample != process->second.end(); ++sample)
        std::cout << "\t    " << sample->second->GetTitle() << " with " << sample->second->GetEntries() << " entries in " << sample->first << " files" << std::endl;
    }
    std::cout << "\tFound " << SIG_samples.size()  << " signal processes:" << std::endl;
    bool isStauStau = false;
    for(auto process = SIG_samples.begin(); process != SIG_samples.end(); ++process)
    {
      std::cout << "\t  " << process->first << ":" << std::endl;
      if(process->first.find("TStauStau") != std::string::npos)
        isStauStau = true;
      for(auto sample = process->second.begin(); sample != process->second.end(); ++sample)
        std::cout << "\t    " << sample->second->GetTitle() << " with " << sample->second->GetEntries() << " entries in " << sample->first << " files" << std::endl;
    }


    // Cut Optimization performed here ----------------------------------------------------------------------------------------------------------------------
    if(BG_samples.size() == 0)
    {
      std::cout << "\tIt doesn't make sense to perform cut optimization without any background samples. Skipping this round of optimization, please check the file " << round->jsonFile() << "." << std::endl;
      continue;
    }
    if(SIG_samples.size() == 0)
    {
      std::cout << "\tIt doesn't make sense to perform cut optimization without any signal samples. Skipping this round of optimization, please check the file " << round->jsonFile() << "." << std::endl;
      continue;
    }

    std::vector<std::pair<std::string,double>> cutValues;
    FOMInfo highestFOM, previousFOM;
    std::pair<std::string,double> highestCut;

    highestFOM.FOM = 0;
    highestFOM.err = 0;
    highestFOM.var = "";
    highestFOM.cutVal = 0;
    previousFOM.FOM = 0;
    previousFOM.err = 0;
    previousFOM.var = "";
    previousFOM.cutVal = 0;

    std::vector<std::string> variables = round->getListOfVariables();
    std::unordered_map<std::string,std::string> variableExpressions = round->getVariableExpressions();
    std::unordered_map<std::string,std::unordered_map<std::string,double>> variableParameterMap = round->getVariableParameterMap();
    std::unordered_map<std::string,std::string> variableLabels = round->getVariableLabels();

    TCut baseSelection = round->baseSelection().c_str();
    TCut signalSelection = round->signalSelection().c_str();
    TCut cumulativeSelection = "";
    std::string pointVariable = round->pointVariable();

    int nSigPoints = getNumberOfPoints(SIG_samples, pointVariable, baseSelection&&signalSelection);
    std::cout << "Running on " << nSigPoints << " signal points." << std::endl;

    TCanvas c1("c1", "c1", 800, 600);

    bool improve = true;
    int nPass = -1;
    while(improve && nSigPoints!=0)
    {
      ++nPass;
      highestFOM.FOM = 0;
      std::cout << round->name() << ": Starting pass " << nPass << ", with " << variables.size() << " variables to optimize cuts on." << std::endl;
      bool isSelected = false;
      if(round->isSelected(nPass))
      {
        std::cout << "  Found a selection for this pass: " << round->getSelection(nPass) << std::endl;
        isSelected = true;
      }

      for(auto variableName = variables.begin(); variableName != variables.end(); ++variableName)
      {
        double minVal = variableParameterMap[*variableName]["minVal"];
        double maxVal = variableParameterMap[*variableName]["maxVal"];
        double step   = variableParameterMap[*variableName]["step"];
        double cutDir = variableParameterMap[*variableName]["cutDir"];
        std::cout << round->name() << "::" << *variableName << " has started processing, with " << (maxVal-minVal)/step + 1 << " steps to be processed." << std::endl;

        std::vector<double> xVals, yVals, xValsErr, yValsErr;
        std::vector<double> signalYield, backgroundYield;

        for(double cutVal = minVal; cutVal <= maxVal; cutVal+=step)
        {
          std::string thisCutStr;
          std::stringstream buf;
          //buf << variableExpressions[*variableName];
          if(cutDir > 0) // Positive values of cut dir means we want to remove events where the variable has a value above the cut
            buf << "<";  // since the cut expression is for the events we want to keep, the expression "seems" inverted.
          else
            buf << ">";
          buf << cutVal;
          buf >> thisCutStr;
          thisCutStr = variableExpressions[*variableName] + thisCutStr;
          std::cout << "  The Cut: " << thisCutStr << std::endl;
          TCut thisCut = thisCutStr.c_str();

          //double sigScale = getSignalScale(SIG_samples, round->signalPoint().c_str(), !isStauStau) * round->iLumi();
          double xSection = round->sigCrossSection();
          double nInitEvents = 10000;
          double nSIG = applyCut(SIG_samples, (baseSelection && signalSelection && cumulativeSelection && thisCut), true, !isStauStau); // Very compute intensive
          double nBG  = applyCut(BG_samples,  (baseSelection && cumulativeSelection && thisCut)) * round->iLumi(); // Very compute intensive
          nSIG = nSIG / nSigPoints;
          double unweightedYield = nSIG;
          nSIG = nSIG * round->iLumi() * xSection / nInitEvents;
          double systErr = 0.15;  // Hard coded systematic error /////////////////////////////////////////////////////////////////////////////////////
          std::cout << "    n(Sig): " << nSIG << std::endl;
          std::cout << "    n(Bkg): " << nBG  << std::endl;

          if(nBG == 0 || nSIG == 0)
            continue;
          double nBGErr = std::sqrt(nBG);
          double nSIGErr = std::sqrt(nSIG);
          if(1/std::sqrt(unweightedYield + nInitEvents) > 1/unweightedYield + 1/nInitEvents)
            nSIGErr = nSIG*std::sqrt(1/unweightedYield + 1/nInitEvents);
          else
            nSIGErr = nSIG*std::sqrt(1/unweightedYield + 1/nInitEvents - 1/std::sqrt(unweightedYield + nInitEvents));
          double dividend = std::sqrt(nBG + (systErr*nBG)*(systErr*nBG));

          double FOM    = nSIG/dividend;
          double FOMErr = 0;
          {
            double SIGContrib = nSIGErr/dividend;
            double BGContrib = nBGErr*nSIG*(1+2*systErr*systErr*nBG)/(2*dividend*dividend*dividend);
            FOMErr = std::sqrt(SIGContrib*SIGContrib + BGContrib*BGContrib);
          }
          std::cout << "    FOM: " << FOM << " +- " << FOMErr << std::endl;

          xVals.push_back(cutVal);
          xValsErr.push_back(0);
          yVals.push_back(FOM);
          yValsErr.push_back(FOMErr);
          signalYield.push_back(nSIG);
          backgroundYield.push_back(nBG);

          if(FOM > highestFOM.FOM)
          {
            highestFOM.FOM = FOM;
            highestFOM.err = FOMErr;
            highestFOM.var = *variableName;
            highestFOM.cutVal = cutVal;
          }
        }

        if(xVals.size() == 0)
        {
          std::cout << "  No points processed" << std::endl;
          continue;
        }
        TGraphErrors FOMGraph(xVals.size(), xVals.data(), yVals.data(), xValsErr.data(), yValsErr.data());
        // http://root.cern.ch/root/html/TAttMarker.html
        //FOMGraph.SetLineColor(2);
        //FOMGraph.SetLineWidth(2);
        //FOMGraph.SetMarkerColor(4);
        //FOMGraph.SetMarkerStyle(21);
        //http://root.cern.ch/root/html/TAttFill.html
        //http://root.cern.ch/root/html/TColor.html
        FOMGraph.SetFillColor(kBlue-9);
        FOMGraph.SetFillStyle(3354);
        FOMGraph.Draw("3A");
        FOMGraph.Draw("same LP");
        FOMGraph.GetXaxis()->SetTitle(variableLabels[*variableName].c_str());
        FOMGraph.GetXaxis()->CenterTitle();
        FOMGraph.GetYaxis()->SetTitle("FOM");

        std::stringstream buf;
        std::string plotName;
        buf << outDir << "/";
        buf << round->name() << "_Pass" << nPass << "_" << *variableName;
        buf >> plotName;

        for(auto ext = plotExt.begin(); ext != plotExt.end(); ++ext)
          c1.SaveAs((plotName + *ext).c_str());


        TGraph signalYieldGraph(xVals.size(), xVals.data(), signalYield.data());
        signalYieldGraph.SetFillColor(kBlue-9);
        signalYieldGraph.SetFillStyle(3354);
        signalYieldGraph.SetMarkerStyle(21);
        signalYieldGraph.Draw("ALP");
        signalYieldGraph.GetXaxis()->SetTitle(variableLabels[*variableName].c_str());
        signalYieldGraph.GetXaxis()->CenterTitle();

        buf.clear();
        plotName = "";
        buf << outDir << "/";
        buf << round->name() << "_Pass" << nPass << "_" << *variableName << "_signalYield";
        buf >> plotName;

        for(auto ext = plotExt.begin(); ext != plotExt.end(); ++ext)
          c1.SaveAs((plotName + *ext).c_str());


        TGraph backgroundYieldGraph(xVals.size(), xVals.data(), backgroundYield.data());
        backgroundYieldGraph.SetFillColor(kBlue-9);
        backgroundYieldGraph.SetFillStyle(3354);
        backgroundYieldGraph.SetMarkerStyle(21);
        backgroundYieldGraph.Draw("ALP");
        backgroundYieldGraph.GetXaxis()->SetTitle(variableLabels[*variableName].c_str());
        backgroundYieldGraph.GetXaxis()->CenterTitle();

        buf.clear();
        plotName = "";
        buf << outDir << "/";
        buf << round->name() << "_Pass" << nPass << "_" << *variableName << "_backgroundYield";
        buf >> plotName;

        for(auto ext = plotExt.begin(); ext != plotExt.end(); ++ext)
          c1.SaveAs((plotName + *ext).c_str());
      }

      if(isSelected)
      {
        improve = true;
        TCut tempCut = round->getSelection(nPass).c_str();
        cumulativeSelection = cumulativeSelection && tempCut;
      }
      else
      {
        if(highestFOM.FOM == 0)
          improve = false;
        else
        {
          if(highestFOM.FOM - previousFOM.FOM > previousFOM.err)
          {
            //Add cut to list of selected cuts
            previousFOM.FOM = highestFOM.FOM;
            previousFOM.err = highestFOM.err;

            std::stringstream buf;
            std::string tempStr;
            if(variableParameterMap[highestFOM.var]["cutDir"] > 0) // Positive values of cut dir means we want to remove events where the variable has a value above the cut
              buf << "<";  // since the cut expression is for the events we want to keep, the expression "seems" inverted.
            else
              buf << ">";
            buf << highestFOM.cutVal;
            buf >> tempStr;
            tempStr = variableExpressions[highestFOM.var] + tempStr;
            std::cout << "Adding cut: " << tempStr << std::endl;
            TCut tempCut = tempStr.c_str();
            cumulativeSelection = cumulativeSelection && tempCut;

            auto newEnd = std::remove(variables.begin(), variables.end(), highestFOM.var);
            variables.erase(newEnd, variables.end());
          }
          else
            improve = false;
        }
      }
    }
    std::cout << round->name() << ": Chose the cuts - " << cumulativeSelection << std::endl;
    // ------------------------------------------------------------------------------------------------------------------------------------------------------
  }

  std::cout << "The list of ignored files, either missing or corrupt, can be found below:" << std::endl;
  for(auto key = FileExists.begin(); key != FileExists.end(); ++key)
  {
    if(!key->second)
      std::cout << "  " << key->first << std::endl;
  }
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

std::vector<std::string> OptimizationRoundInfo::getListOfVariables()
{
  std::vector<std::string> retVal;

  for(auto variable = _variables.begin(); variable != _variables.end(); ++variable)
    retVal.push_back(variable->name());

  return retVal;
}

std::unordered_map<std::string,std::string> OptimizationRoundInfo::getVariableLabels()
{
  std::unordered_map<std::string,std::string> retVal;

  for(auto variable = _variables.begin(); variable != _variables.end(); ++variable)
  {
    retVal[variable->name()] = variable->label();
    if(retVal[variable->name()] == "")
      retVal[variable->name()] = variable->name();
  }

  return retVal;
}

std::unordered_map<std::string,std::unordered_map<std::string,double>> OptimizationRoundInfo::getVariableParameterMap()
{
  std::unordered_map<std::string,std::unordered_map<std::string,double>> retVal;

  for(auto variable = _variables.begin(); variable != _variables.end(); ++variable)
  {
    std::unordered_map<std::string,double> tempVal;

    tempVal["minVal"] = variable->minVal();
    tempVal["maxVal"] = variable->maxVal();
    tempVal["step"]   = variable->step();
    tempVal["cutDir"] = (variable->cutDir()=="below")?-1:1;

    retVal[variable->name()] = tempVal;
  }

  return retVal;
}

std::unordered_map<std::string,std::string> OptimizationRoundInfo::getVariableExpressions()
{
  std::unordered_map<std::string,std::string> retVal;

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
  _cutDir = "below";
}

OptimizationVariableInfo::~OptimizationVariableInfo()
{
}

int getNumberOfPoints(fileChains files, std::string variable, TCut cut)
{
  gROOT->cd();
  int retVal = 0;
  int maxVal = 501*1000+1;
  TH1D finalHist("finalHist", "finalHist", maxVal, 0, maxVal);

  for(auto process = files.begin(); process != files.end(); ++process)
  {
    for(auto sample = process->second.begin(); sample != process->second.end(); ++sample)
    {
      if(sample->first == 0)
        continue;

      TH1D* temp_hist = new TH1D("temp_hist", "temp_hist", maxVal, 0, maxVal);
      sample->second->Draw((variable+">>temp_hist").c_str(), cut*"weight", "goff");
      finalHist.Add(temp_hist);

      delete temp_hist;
    }
  }

  int nBins = finalHist.GetNbinsX();
  for(int i = 1; i <= nBins; ++i)
  {
    if(finalHist.GetBinContent(i) != 0)
      retVal++;
  }

  return retVal;
}

double getSignalScale(fileChains files,  TCut cut, bool correctFiles)
{
  double retVal = 0;

  if(!gDoneSignalScale)
  {
    for(auto process = files.begin(); process != files.end(); ++process)
    {
      for(auto sample = process->second.begin(); sample != process->second.end(); ++sample)
      {
        if(sample->first == 0)
          continue;
        TH1D temp_histo("temp_histo", "temp_histo", 1, 0, 20);
        sample->second->Draw("weight>>temp_histo", cut*"weight", "goff");

        double count = temp_histo.GetBinContent(0) + temp_histo.GetBinContent(1) + temp_histo.GetBinContent(2);

        if(correctFiles)
          count = count/sample->first;
        retVal += count;
      }
    }
  }
  else
    retVal = gSignalScale;

  gSignalScale = retVal;
  gDoneSignalScale = true;
  return retVal;
}

double applyCut(fileChains files,  TCut cut, bool isSig, bool correctFiles, bool normalize)
{
  gROOT->cd();
  double retVal = 0;

  if(!gDoneScale && isSig)
  {
    gScale = 0;
    for(auto process = files.begin(); process != files.end(); ++process)
    {
      for(auto sample = process->second.begin(); sample != process->second.end(); ++sample)
      {
        if(sample->first == 0)
          continue;
        TH1D temp_histo("temp_histo", "temp_histo", 1, 0, 20);
        sample->second->Draw("puWeight>>temp_histo", "puWeight", "goff");

        double count = temp_histo.GetBinContent(0) + temp_histo.GetBinContent(1) + temp_histo.GetBinContent(2);

        if(correctFiles)
          count = count/sample->first;
        gScale += count;
      }
    }

    gScale = 1/gScale;
    gDoneScale = true;
  }

  for(auto process = files.begin(); process != files.end(); ++process)
  {
    for(auto sample = process->second.begin(); sample != process->second.end(); ++sample)
    {
      if(sample->first == 0)
        continue;
//      TEventList selectedEvents("selectedList");
//      sample->second->Draw(">>selectedList", cut, "goff");
      TH1D temp_histo("temp_histo", "temp_histo", 1, 0, 20);
      if(isSig)
        sample->second->Draw("puWeight>>temp_histo", cut*"puWeight", "goff");
      else
        sample->second->Draw("weight>>temp_histo", cut*"weight", "goff");

      double count = temp_histo.GetBinContent(0) + temp_histo.GetBinContent(1) + temp_histo.GetBinContent(2);

      /*double weight = 0;
      double count = 0;
      sample->second->SetBranchAddress("weight", &weight);

      for(auto index = selectedEvents.GetN()-1; index >= 0; --index)
      {
        Long64_t entry = selectedEvents.GetEntry(index);
        sample->second->GetEntry(entry);
        count += weight;
      }
      sample->second->ResetBranchAddresses();// */

      if(correctFiles)
        count = count/sample->first;
      retVal += count;
    }
  }

  return retVal;
}

std::vector<OptimizationRoundInfo> getRoundsFromJSON(JSONWrapper::Object& json)
{
  std::vector<OptimizationRoundInfo> retVal;

  std::vector<JSONWrapper::Object> rounds = json["optim"].daughters();
  for(auto round = rounds.begin(); round != rounds.end(); ++round)
  {
    OptimizationRoundInfo roundInfo;
    //roundInfo._pointVariable = "stauMass*1000+neutralinoMass";

    std::string name = round->getString("name", "");
    if(name != "")
      roundInfo._name = name;
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

    auto selected = (*round)["selected"].daughters();
    for(auto selection = selected.begin(); selection != selected.end(); ++selection)
    {
      std::string temp = selection->getString("selection", "");
      if(temp != "")
        roundInfo._selected.push_back(temp);
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
      variableInfo._step = variable->getDouble("step", 0);
      if(variableInfo._step <= 0 || variableInfo._step > variableInfo._maxVal - variableInfo._minVal)
      {
        std::cout << roundInfo._name << "::" << variableInfo._name << ": step must be a resonable and valid value. Continuing..." << std::endl;
        continue;
      }
      variableInfo._cutDir = variable->getString("cutDir", "below");
      std::transform(variableInfo._cutDir.begin(), variableInfo._cutDir.end(), variableInfo._cutDir.begin(), ::tolower);
      if(variableInfo._cutDir != "below" && variableInfo._cutDir != "above")
      {
        std::cout << roundInfo._name << "::" << variableInfo._name << ": the cut direction (cutDir) must be either below or above. Continuing..." << std::endl;
        continue;
      }
      variableInfo._label = variable->getString("label", "");

      roundInfo._variables.push_back(variableInfo);
    }

    retVal.push_back(roundInfo);
  }

  return retVal;
}

fileChains getChainsFromJSON(JSONWrapper::Object& json, std::string RootDir, std::string type, std::string treename, std::string customExtension)
{
  std::vector<std::pair<std::string,std::vector<std::pair<int,TChain*>>>> retVal;
  std::pair<std::string,std::vector<std::pair<int,TChain*>>> tempProcess;
  std::vector<std::pair<int,TChain*>> tempSamples;
  std::pair<int,TChain*> tempSample;

  if(type != "BG" && type != "SIG")
  {
    std::cout << "Unknown sample type requested." << std::endl;
    assert(type == "BG" || type == "SIG");
  }

  std::vector<JSONWrapper::Object> processes = json["proc"].daughters();
  for(auto process = processes.begin(); process != processes.end(); ++process)
  {
    bool isData = (*process)["isdata"].toBool();
    bool isSig  = !isData && (*process).isTag("spimpose") && (*process)["spimpose"].toBool();
    bool isMC   = !isData && !isSig;

    if(isData) // Here we are enforcing for the data samples to not even be present, might not make sense for a data-driven background estimation
      continue;

    if(type == "SIG" && isMC) // Here we make sure we are only processing the requested processes
      continue;
    if(type == "BG" && isSig)
      continue;

    tempProcess.first = (*process).getString("tag", "Sample");

    std::string filtExt;
    if((*process).isTag("mctruthmode"))
    {
      std::stringstream buf;
      buf << "_filt" << (*process)["mctruthmode"].toInt();
      buf >> filtExt;
    }

    std::vector<JSONWrapper::Object> samples = (*process)["data"].daughters();
    tempSamples.clear();
    for(auto sample = samples.begin(); sample != samples.end(); ++sample)
    {
      int nFiles = (*sample).getInt("split", 1);
      tempSample.first = 0;
      tempSample.second = new TChain(treename.c_str(), ((*sample).getString("dtag", "") + (*sample).getString("suffix", "")).c_str());
      for(int f = 0; f < nFiles; ++f)
      {
        std::string segmentExt;
        if(nFiles != 1)
        {
          std::stringstream buf;
          buf << "_" << f;
          buf >> segmentExt;
        }

        std::string fileName = RootDir + "/" + (*sample).getString("dtag", "") + (*sample).getString("suffix", "") + segmentExt + filtExt + customExtension + ".root";
        TFile* file = new TFile(fileName.c_str());
        bool& fileExists = FileExists[fileName];

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

        // Chain the valid files together
        tempSample.second->Add(fileName.c_str());
        ++tempSample.first;
      }
      tempSamples.push_back(tempSample);
    }
    tempProcess.second = tempSamples;
    retVal.push_back(tempProcess);
  }

  return retVal;
}

void printHelp()
{
  std::cout << "runCutOptimizer help - There are the following options:" << std::endl << std::endl;

  std::cout << "--help    -->  Print this help message" << std::endl;
  std::cout << "--json    -->  Configuration file for cut optimization, should define which files to run on, where they are located, the integrated luminosity and the variables to optimize the cuts on" << std::endl;
  std::cout << "--outDir  -->  Path to the directory where to output plots and tables (will be created if it doesn't exist)" << std::endl;
  std::cout << "--plotExt -->  Extension format with which to save the plots, repeat this command if multiple formats are desired" << std::endl;

  std::cout << std::endl << "Example command:" << std::endl << "\trunCutOptimizer --json optimization_options.json --outDir ./OUT/" << std::endl;
  return;
}

