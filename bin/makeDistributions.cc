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
#include "TStyle.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TTree.h"
#include "TCut.h"
#include "TEventList.h"
#include "TInterpreter.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include "UserCode/llvv_fwk/interface/JSONWrapper.h"
#include "UserCode/llvv_fwk/interface/llvvObjects.h"


class MyVariable;
class MyStyle;
struct SampleFiles;
struct ProcessFiles;
struct ProcessHists;
class My2DPlot;

typedef std::vector<std::pair<std::string,std::vector<std::pair<int,TChain*>>>> fileChains;
typedef std::map<std::string,std::vector<ProcessFiles>> dataChains;
typedef std::map<std::string,std::vector<ProcessHists>> dataHists;

void printHelp();
std::vector<MyVariable> getVariables(JSONWrapper::Object& json);
fileChains old_getChainsFromJSON(JSONWrapper::Object& json, std::string RootDir, std::string type, std::string treename, std::string customExtension);
THStack* old_getStack(fileChains files,  TCut cut, MyVariable variable, bool correctFiles = true);
TH1D* getHist(fileChains files,  TCut cut, MyVariable variable, bool correctFiles = true);
dataChains getChainsFromJSON(std::string jsonFile, std::string ttreeName, std::string rootDir, std::string customExtension);
dataHists getHists(dataChains& chains, TCut cut, TCut sigCut, MyVariable& variable, double iLumi, double sigXSec = 0, int sigNInitEvents = 0, std::string pointVar = "", bool isMultipointSignalSample = false, bool normalize = false);

std::vector<My2DPlot> get2DPlots(JSONWrapper::Object& json);
void make2D(std::string outDir, std::vector<std::string> plotExt, fileChains files, TCut cut, My2DPlot TwoDPlot, bool correctFiles = true);
std::string cleanString(std::string str, std::string illegal, char replacement);


class MyStyle
{
public:
  MyStyle():marker_(1),lcolor_(kBlack),mcolor_(kBlack),fcolor_(kWhite),lwidth_(1),lstyle_(1){};

  inline int marker(){return marker_;};
  inline int lcolor(){return lcolor_;};
  inline int mcolor(){return mcolor_;};
  inline int fcolor(){return fcolor_;};
  inline int lwidth(){return lwidth_;};
  inline int lstyle(){return lstyle_;};

  friend fileChains old_getChainsFromJSON(JSONWrapper::Object& json, std::string RootDir, std::string type, std::string treename, std::string customExtension);
  friend dataChains getChainsFromJSON(std::string jsonFile, std::string ttreeName, std::string rootDir, std::string customExtension);

private:
  int marker_;
  int lcolor_;
  int mcolor_;
  int fcolor_;
  int lwidth_;
  int lstyle_;

protected:
};

class MyVariable
{
public:
//  OptimizationVariable();
//  ~OptimizationVariable();

  // If using pointers the following lines should be uncommented and the respective functions defined
  // More info: http://www.cplusplus.com/articles/y8hv0pDG/
  //OptimizationVariable(const OptimizationVariable& other);
  //OptimizationVariable& operator=(const OptimizationVariable& rhs);

  inline std::string& name(){return _name;};
  inline std::string& expression(){return _expression;};
  inline double& minVal(){return _minVal;};
  inline double& maxVal(){return _maxVal;};
  inline int& bins(){return _bins;};
  inline std::string& label(){return _label;};

  friend std::vector<MyVariable> getVariables(JSONWrapper::Object& json);
  friend std::vector<My2DPlot> get2DPlots(JSONWrapper::Object& json);

private:
  std::string _name;
  std::string _expression;
  double _minVal, _maxVal;
  int _bins;
  std::string _label;

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
  MyStyle style;
  std::vector<SampleFiles> samples;
};

struct ProcessHists
{
  std::string name;
  TH1* hist;
  doubleUnc yield;
};

class My2DPlot
{
public:

  inline MyVariable& xVar(){return xVar_;};
  inline MyVariable& yVar(){return yVar_;};

  friend std::vector<My2DPlot> get2DPlots(JSONWrapper::Object& json);

private:
  MyVariable xVar_;
  MyVariable yVar_;

protected:
};

std::unordered_map<std::string,bool> FileExists;
std::unordered_map<std::string,MyStyle> styles;

//
// Output codes:
// 0 - Everything has run fine
// 1 - Invalid arguments
// 2 - Problem parsing the arguments
//

int main(int argc, char** argv)
{
  AutoLibraryLoader::enable();

  std::string jsonFile, variablesFile;
  std::string outDir = "./OUT/";
  std::string inDir;
  std::vector<std::string> plotExt;
  std::string baseSelection = "selected";
  std::string signalSelection = "";
  std::string extraSignal = "";
  std::string ttree_name = "Events";
  std::string customExtension = "_summary";
  bool normalize = false;
  double iLumi = 0;
  double sigXSec = 0;
  int sigNInitEvents = 0;
  std::string pointVar = "";
  bool isMultipointSignalSample = false;
  bool unblind = false;
  bool printProcesses = false;

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

    if(arg.find("--iLumi") != std::string::npos)
    {
      std::stringstream buf;
      std::string temp = argv[i+1];
      buf << temp;
      buf >> iLumi;
      ++i;
    }

    if(arg.find("--printProcesses") != std::string::npos)
    {
      printProcesses = true;
    }

    if(arg.find("--sigXSec") != std::string::npos)
    {
      std::stringstream buf;
      std::string temp = argv[i+1];
      buf << temp;
      buf >> sigXSec;
      ++i;
    }

    if(arg.find("--sigNInitEvents") != std::string::npos)
    {
      std::stringstream buf;
      std::string temp = argv[i+1];
      buf << temp;
      buf >> sigNInitEvents;
      ++i;
    }

    if(arg.find("--pointVar") != std::string::npos)
    {
      pointVar = argv[i+1];
      ++i;
    }

    if(arg.find("--json") != std::string::npos)
    {
      jsonFile = argv[i+1];
      // TODO: Check if the file exists
      ++i;
    }

    if(arg.find("--variables") != std::string::npos)
    {
      variablesFile = argv[i+1];
      // TODO: Check if the file exists
      ++i;
    }

    if(arg.find("--normalize") != std::string::npos)
    {
      normalize = true;
    }

    if(arg.find("--unblind") != std::string::npos)
    {
      unblind = true;
    }

    if(arg.find("--outDir") != std::string::npos)
    {
      outDir = argv[i+1];
      ++i;
    }

    if(arg.find("--inDir") != std::string::npos)
    {
      inDir = argv[i+1];
      ++i;
    }

    if(arg.find("--plotExt") != std::string::npos)
    {
      plotExt.push_back(argv[i+1]);
      ++i;
    }

    if(arg.find("--baseSelection") != std::string::npos)
    {
      baseSelection = argv[i+1];
      ++i;
    }

    if(arg.find("--signalSelection") != std::string::npos)
    {
      signalSelection = argv[i+1];
      ++i;
    }

    if(arg.find("--extraSignal") != std::string::npos)
    {
      extraSignal = argv[i+1];
      ++i;
    }

    if(arg.find("--ttree") != std::string::npos)
    {
      ttree_name = argv[i+1];
      ++i;
    }

    if(arg.find("--customExtension") != std::string::npos)
    {
      customExtension = argv[i+1];
      ++i;
    }
  }
  if(plotExt.size() == 0)
    plotExt.push_back(".png");

  system(("mkdir -p " + outDir).c_str());

  if(iLumi <= 0)
  {
    if(!normalize)
    {
      std::cout << "You should define an integrated luminosity if not plotting the normalized histograms." << std::endl;
      return 1;
    }
    else
      iLumi = 1;
  }

  if(pointVar == "")
  {
    isMultipointSignalSample = false;
    sigXSec = 0;
    sigNInitEvents = 0;
  }
  else
  {
    isMultipointSignalSample = true;
    if(sigNInitEvents <= 0 || sigXSec <= 0)
    {
      std::cout << "If you are using a multipoint signal sample, by defining pointVar, you must define sigNInitEvents and sigXSec." << std::endl;
      return 1;
    }
  }

  if(jsonFile == "" || variablesFile == "" || inDir == "")
  {
    std::cout << "You should define at least the following arguments: json, variables, inDir" << std::endl;
    std::cout << "For more information, consult the help (\"makeDistributions --help\")" << std::endl;
    return 2;
  }

  JSONWrapper::Object json(jsonFile, true);
  JSONWrapper::Object variables_json(variablesFile, true);

  std::vector<MyVariable> variables = getVariables(variables_json);
  std::vector<My2DPlot> TwoDPlots = get2DPlots(variables_json);

  auto  BG_samples = old_getChainsFromJSON(json, inDir,  "BG", ttree_name, customExtension);
  auto SIG_samples = old_getChainsFromJSON(json, inDir, "SIG", ttree_name, customExtension);
  auto processes = getChainsFromJSON(jsonFile, ttree_name, inDir, customExtension);
  dataHists yields;

  TCut BGCut  = baseSelection.c_str();
  TCut SIGCut = BGCut && signalSelection.c_str();

  TCanvas c1("c1", "c1", 800, 600);
  c1.SetLogy();
  gStyle->SetOptStat(0);
  doubleUnc bgYield{0,0}, sigYield{0,0}, dataYield{0,0}, extraSigYield{0,0};
  for(auto variable = variables.begin(); variable != variables.end(); ++variable)
  {
    std::cout << "Processing " << variable->name() << std::endl;
    std::string label = variable->label();
    if(label == "")
      label = variable->name();

    dataHists hists = getHists(processes, BGCut, SIGCut, *variable, iLumi, sigXSec, sigNInitEvents, pointVar, isMultipointSignalSample, normalize);

    THStack* bgStack = new THStack("Background", ("Background;"+label+";% Events").c_str());
    TH1D* sigHist = NULL;
    TH1D* dataHist = NULL;
    TH1D* extraSigHist = NULL;
    bgYield.value = 0, sigYield.value = 0, dataYield.value = 0, extraSigYield.value = 0;
    bgYield.uncertainty = 0, sigYield.uncertainty = 0, dataYield.uncertainty = 0, extraSigYield.uncertainty = 0;

    for(auto process = hists["BG"].begin(); process != hists["BG"].end(); ++process)
    {
      bgStack->Add(process->hist);
      bgYield.value += process->yield.value;
      bgYield.uncertainty += process->yield.uncertainty*process->yield.uncertainty;
    }
    bgYield.uncertainty = std::sqrt(bgYield.uncertainty);
    for(auto process = hists["SIG"].begin(); process != hists["SIG"].end(); ++process)
    {
      if(sigHist == NULL)
        sigHist = static_cast<TH1D*>(process->hist->Clone());
      else
        sigHist->Add(process->hist);
      sigYield.value += process->yield.value;
      sigYield.uncertainty += process->yield.uncertainty*process->yield.uncertainty;
    }
    sigYield.uncertainty = std::sqrt(sigYield.uncertainty);
    for(auto process = hists["Data"].begin(); process != hists["Data"].end(); ++process)
    {
      if(dataHist == NULL)
        dataHist = static_cast<TH1D*>(process->hist->Clone());
      else
        dataHist->Add(process->hist);
      dataYield.value += process->yield.value;
      dataYield.uncertainty += process->yield.uncertainty*process->yield.uncertainty;
    }
    dataYield.uncertainty = std::sqrt(dataYield.uncertainty);

    bgStack->Draw("hist");
    double max = sigHist->GetMaximum();
    if(bgStack->GetMaximum() > max)
      max = bgStack->GetMaximum();
    if(dataHist != NULL && unblind)
      if(dataHist->GetMaximum() > max)
        max = dataHist->GetMaximum();

    if(printProcesses)
    {
      yields = hists;
    }

    if(extraSignal != "")
    {
      TCut SIG2Cut = extraSignal.c_str();
      dataChains extraSignalChains;
      extraSignalChains["SIG"] = processes["SIG"];
      dataHists extraHists = getHists(extraSignalChains, BGCut, SIG2Cut, *variable, iLumi, sigXSec, sigNInitEvents, pointVar, isMultipointSignalSample, normalize);

      for(auto process = extraHists["SIG"].begin(); process != extraHists["SIG"].end(); ++process)
      {
        if(extraSigHist == NULL)
          extraSigHist = static_cast<TH1D*>(process->hist->Clone());
        else
          extraSigHist->Add(process->hist);
        extraSigYield.value += process->yield.value;
        extraSigYield.uncertainty += process->yield.uncertainty*process->yield.uncertainty;
        delete process->hist;
      }
      extraSigYield.uncertainty = std::sqrt(extraSigYield.uncertainty);

      if(extraSigHist->GetMaximum() > max)
        max = extraSigHist->GetMaximum();

      bgStack->SetMaximum(max);
      sigHist->Draw("same");
      extraSigHist->SetLineColor(kBlue);
      extraSigHist->Draw("same");
      if(printProcesses)
      {
        yields["ExtraSig"] = extraHists["SIG"];
      }
    }
    else
    {
      bgStack->SetMaximum(max);
      sigHist->Draw("same");
    }

    if(unblind)
      dataHist->Draw("same");

    c1.BuildLegend(0.88, 0.67, 1, 1);

    for(auto ext = plotExt.begin(); ext != plotExt.end(); ++ext)
      c1.SaveAs((outDir + variable->name() + *ext).c_str());

    delete sigHist;
    delete bgStack;
    if(dataHist != NULL)
      delete dataHist;
    if(extraSigHist != NULL)
      delete extraSigHist;

    for(auto type = hists.begin(); type != hists.end(); ++type)
    {
      for(auto process = type->second.begin(); process != type->second.end(); ++process)
      {
        delete process->hist;
        process->hist = NULL;
      }
    }
  }

  if(printProcesses)
  {
    std::cout << "Process cutflow breakdown:" << std::endl;
    for(auto type = yields.begin(); type != yields.end(); ++type)
    {
      if(type->first == "Data" && !unblind)
        continue;
      std::cout << "  " << type->first << std::endl;
      for(auto process = type->second.begin(); process != type->second.end(); ++process)
      {
        std::cout << "    " << process->name << ": " << process->yield << std::endl;
      }
    }
  }

  std::cout << "Yields:" << std::endl;
  std::cout << "  MC: " << bgYield << std::endl;
  if(unblind)
    std::cout << "  Data: " << dataYield << std::endl;
  std::cout << "  Sig: " << sigYield << std::endl;
  if(extraSignal != "")
    std::cout << "  ExtraSig: " << extraSigYield << std::endl;

  for(auto TwoDPlot = TwoDPlots.begin(); TwoDPlot != TwoDPlots.end(); ++TwoDPlot)
  {
    make2D(outDir, plotExt, BG_samples, BGCut, *TwoDPlot);
    make2D(outDir, plotExt, SIG_samples, SIGCut, *TwoDPlot, false);
  }

  std::cout << "The list of ignored files, either missing or corrupt, can be found below:" << std::endl;
  for(auto key = FileExists.begin(); key != FileExists.end(); ++key)
  {
    if(!key->second)
      std::cout << "  " << key->first << std::endl;
  }

  return 0;
}

dataHists getHists(dataChains& chains, TCut cut, TCut sigCut, MyVariable& variable, double iLumi, double sigXSec, int sigNInitEvents, std::string pointVar, bool isMultipointSignalSample, bool normalize)
{
  dataHists retVal;

  std::string name = variable.name();
  std::string expression = variable.expression();
  std::string label = variable.label();
  if(label == "")
    label = name;


  int nSignalPoints = 1;
  if(isMultipointSignalSample)
  {
    nSignalPoints = 0;
    int maxVal = 501*1000+1; //Change-me if you are having problems with multipart signal samples and it only finds some of them
    TH1D tempPoints("points", "points", maxVal, 0, maxVal);
    for(auto process = chains["SIG"].begin(); process != chains["SIG"].end(); ++process)
    {
      for(auto sample = process->samples.begin(); sample != process->samples.end(); ++sample)
      {
        if(sample->nFiles == 0)
          continue;

        TH1D* temp_hist = new TH1D("temp_hist", "temp_hist", maxVal, 0, maxVal);
        sample->chain->Draw((pointVar+">>temp_hist").c_str(), (cut && sigCut)*"weight", "goff");
        tempPoints.Add(temp_hist);
        delete temp_hist;
      }
    }

    int nBins = tempPoints.GetNbinsX();
    for(int i = 1; i <= nBins; ++i)
      if(tempPoints.GetBinContent(i) != 0)
        ++nSignalPoints;
  }


  for(auto type = chains.begin(); type != chains.end(); ++type)
  {
    std::vector<ProcessHists> tempProcInfo;
    double scale = 0;

    for(auto process = type->second.begin(); process != type->second.end(); ++process)
    {
      ProcessHists procHist;

      procHist.name = process->name;
      procHist.hist = new TH1D((name+process->name).c_str(), (process->name+";"+label+";Events").c_str(), variable.bins(), variable.minVal(), variable.maxVal());
      procHist.hist->Sumw2();

      for(auto sample = process->samples.begin(); sample != process->samples.end(); ++sample)
      {
        if(sample->nFiles == 0)
          continue;

        TH1D tempHist("temp", (name+";"+label+";Events").c_str(), variable.bins(), variable.minVal(), variable.maxVal());
        tempHist.Sumw2();

        if(nSignalPoints > 1 && type->first == "SIG")
          sample->chain->Draw((expression+">>temp").c_str(), (cut && sigCut)*"puWeight", "goff");
        else
          sample->chain->Draw((expression+">>temp").c_str(), cut*"weight", "goff");

        if(type->first != "Data")
        {
          if(isMultipointSignalSample && type->first == "SIG")
          {
            if(nSignalPoints > 1)
              tempHist.Scale(iLumi*sigXSec/(nSignalPoints*sigNInitEvents));
          }
          else
          {
            tempHist.Scale(iLumi/sample->nFiles);
          }
        }

        tempHist.SetLineColor  (process->style.lcolor());
        tempHist.SetMarkerColor(process->style.mcolor());
        tempHist.SetFillColor  (process->style.fcolor());
        tempHist.SetLineWidth  (process->style.lwidth());
        tempHist.SetLineStyle  (process->style.lstyle());
        tempHist.SetMarkerStyle(process->style.marker());

        procHist.hist->Add(&tempHist);
      }

      procHist.hist->SetLineColor  (process->style.lcolor());
      procHist.hist->SetMarkerColor(process->style.mcolor());
      procHist.hist->SetFillColor  (process->style.fcolor());
      procHist.hist->SetLineWidth  (process->style.lwidth());
      procHist.hist->SetLineStyle  (process->style.lstyle());
      procHist.hist->SetMarkerStyle(process->style.marker());

      procHist.yield.value = procHist.hist->IntegralAndError(0, variable.bins()+1, procHist.yield.uncertainty);
      scale += procHist.yield.value;
      tempProcInfo.push_back(procHist);
    }

    if(normalize)
    {
      for(auto process = tempProcInfo.begin(); process != tempProcInfo.end(); ++process)
      {
        process->hist->Scale(1./scale);
      }
    }

    retVal[type->first] = tempProcInfo;
  }

  return retVal;
}

void make2D(std::string outDir, std::vector<std::string> plotExt, fileChains files, TCut cut, My2DPlot TwoDPlot, bool correctFiles)
{
  std::string xName = TwoDPlot.xVar().name();
  std::string xExpression = TwoDPlot.xVar().expression();
  std::string xLabel = TwoDPlot.xVar().label();
  if(xLabel == "")
    xLabel = xName;
  std::string yName = TwoDPlot.yVar().name();
  std::string yExpression = TwoDPlot.yVar().expression();
  std::string yLabel = TwoDPlot.yVar().label();
  if(yLabel == "")
    yLabel = yName;

  TCanvas c2("c2", "c2", 800, 800);
  c2.SetLogz();
  for(auto process = files.begin(); process != files.end(); ++process)
  {
    TH2D* processHist = new TH2D((xName+yName+process->first).c_str(), (process->first+";"+xLabel+";"+yLabel).c_str(), TwoDPlot.xVar().bins(), TwoDPlot.xVar().minVal(), TwoDPlot.xVar().maxVal(), TwoDPlot.yVar().bins(), TwoDPlot.yVar().minVal(), TwoDPlot.yVar().maxVal());
    for(auto sample = process->second.begin(); sample != process->second.end(); ++sample)
    {
      if(sample->first == 0)
        continue;

      TH2D tempHisto("temp", (";"+xLabel+";"+yLabel).c_str(), TwoDPlot.xVar().bins(), TwoDPlot.xVar().minVal(), TwoDPlot.xVar().maxVal(), TwoDPlot.yVar().bins(), TwoDPlot.yVar().minVal(), TwoDPlot.yVar().maxVal());
      sample->second->Draw((yExpression+":"+xExpression+">>temp").c_str(), cut*"weight", "goff");
      if(correctFiles)
        tempHisto.Scale(1./sample->first);
      processHist->Add(&tempHisto);
    }

    processHist->Draw("colz");
    processHist->SetMinimum(processHist->GetMinimum());

    for(auto ext = plotExt.begin(); ext != plotExt.end(); ++ext)
      c2.SaveAs((outDir + cleanString(process->first, "\\\"?|<>:/#+{} ", '_') + "_" + xName + "_" + yName + *ext).c_str());

    delete processHist;
  }

  return;
}

std::string cleanString(std::string str, std::string illegal, char replacement)
{
  for(auto it = str.begin(); it != str.end(); ++it)
  {
    bool found = illegal.find(*it) != std::string::npos;
    if(found)
      *it = replacement;
  }

  return str;
}

THStack* old_getStack(fileChains files,  TCut cut, MyVariable variable, bool correctFiles)
{
  std::string name = variable.name();
  std::string expression = variable.expression();
  std::string label = variable.label();
  if(label == "")
    label = name;

  THStack* retVal = new THStack(name.c_str(), (name+";"+label+";Events").c_str());
  TH1D* tempHist = getHist(files, cut, variable, correctFiles);
  double scale = tempHist->Integral();
  delete tempHist;

  for(auto process = files.begin(); process != files.end(); ++process)
  {
    TH1D* processHist = new TH1D((name+process->first).c_str(), (process->first+";"+label+";Events").c_str(), variable.bins(), variable.minVal(), variable.maxVal());
    for(auto sample = process->second.begin(); sample != process->second.end(); ++sample)
    {
      if(sample->first == 0)
        continue;

      TH1D tempHisto("temp", (name+";"+label+";Events").c_str(), variable.bins(), variable.minVal(), variable.maxVal());
      sample->second->Draw((expression+">>temp").c_str(), cut*"weight", "goff");
      if(correctFiles)
        tempHisto.Scale(1./sample->first);
      processHist->Add(&tempHisto);
    }

    processHist->Scale(1/scale);
    processHist->SetLineColor  (styles[process->first].lcolor());
    processHist->SetMarkerColor(styles[process->first].mcolor());
    processHist->SetFillColor  (styles[process->first].fcolor());
    processHist->SetLineWidth  (styles[process->first].lwidth());
    processHist->SetLineStyle  (styles[process->first].lstyle());
    processHist->SetMarkerStyle(styles[process->first].marker());

    retVal->Add(processHist);
  }

  return retVal;
}

TH1D* getHist(fileChains files,  TCut cut, MyVariable variable, bool correctFiles)
{
  std::string name = variable.name();
  std::string expression = variable.expression();
  std::string label = variable.label();
  if(label == "")
    label = name;

  TH1D* retVal = new TH1D(name.c_str(), (name+";"+label+";Events").c_str(), variable.bins(), variable.minVal(), variable.maxVal());

  for(auto process = files.begin(); process != files.end(); ++process)
  {
    MyStyle style = styles[process->first];
    for(auto sample = process->second.begin(); sample != process->second.end(); ++sample)
    {
      if(sample->first == 0)
        continue;

      TH1D tempHisto("temp", (name+";"+label+";Events").c_str(), variable.bins(), variable.minVal(), variable.maxVal());
      sample->second->Draw((expression+">>temp").c_str(), cut*"weight", "goff");
      if(correctFiles)
        tempHisto.Scale(1./sample->first);
      retVal->Add(&tempHisto);
    }

    retVal->SetLineColor  (style.lcolor());
    retVal->SetMarkerColor(style.mcolor());
    retVal->SetFillColor  (style.fcolor());
    retVal->SetLineWidth  (style.lwidth());
    retVal->SetLineStyle  (style.lstyle());
    retVal->SetMarkerStyle(style.marker());

  }


  return retVal;
}

fileChains old_getChainsFromJSON(JSONWrapper::Object& json, std::string RootDir, std::string type, std::string treename, std::string customExtension)
{
  std::vector<std::pair<std::string,std::vector<std::pair<int,TChain*>>>> retVal;
  std::pair<std::string,std::vector<std::pair<int,TChain*>>> tempProcess;
  std::vector<std::pair<int,TChain*>> tempSamples;
  std::pair<int,TChain*> tempSample;
  MyStyle style;

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

    style.lcolor_ = 1;
    style.mcolor_ = 1;
    style.fcolor_ = 0;
    style.lwidth_ = 1;
    style.lstyle_ = 1;
    style.marker_ = 1;
    if(process->isTag("color"))
    {
      style.lcolor_ = process->getInt("color");
      style.mcolor_ = style.lcolor_;
      style.fcolor_ = style.lcolor_;
    }
    if(process->isTag("lcolor"))
      style.lcolor_ = process->getInt("lcolor");
    if(process->isTag("mcolor"))
      style.mcolor_ = process->getInt("mcolor");
    if(process->isTag("fcolor"))
      style.fcolor_ = process->getInt("fcolor");
    if(process->isTag("fill"))
      style.fcolor_ = process->getInt("fill");
    if(process->isTag("lwidth"))
      style.lwidth_ = process->getInt("lwidth");
    if(process->isTag("lstyle"))
      style.lstyle_ = process->getInt("lstyle");
    if(process->isTag("marker"))
      style.marker_ = process->getInt("marker");
    styles[tempProcess.first] = style;

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
      for(int filen = 0; filen < nFiles; ++filen)
      {
        std::string segmentExt;
        if(nFiles != 1)
        {
          std::stringstream buf;
          buf << "_" << filen;
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

std::vector<My2DPlot> get2DPlots(JSONWrapper::Object& json)
{
  std::vector<My2DPlot> retVal;

  std::vector<JSONWrapper::Object> TwoDPlots = json["2DPlots"].daughters();
  for(auto TwoDPlot = TwoDPlots.begin(); TwoDPlot != TwoDPlots.end(); ++TwoDPlot)
  {
    My2DPlot TwoDPlotInfo;

    TwoDPlotInfo.xVar_._name = TwoDPlot->getString("xname", "");
    TwoDPlotInfo.yVar_._name = TwoDPlot->getString("yname", "");
    if(TwoDPlotInfo.xVar()._name == "" || TwoDPlotInfo.yVar()._name == "")
    {
      std::cout << "All variables must have names. Continuing..." << std::endl;
      continue;
    }
    TwoDPlotInfo.xVar_._expression = TwoDPlot->getString("xexpression", "");
    TwoDPlotInfo.yVar_._expression = TwoDPlot->getString("yexpression", "");
    if(TwoDPlotInfo.xVar_._expression == "")
    {
      std::cout << TwoDPlotInfo.yVar_._name << "_Vs_" << TwoDPlotInfo.xVar_._name << ": The variable " << TwoDPlotInfo.xVar_._name << " must have and expression, it must be a valid root exprtession. Continuing..." << std::endl;
      continue;
    }
    if(TwoDPlotInfo.yVar_._expression == "")
    {
      std::cout << TwoDPlotInfo.yVar_._name << "_Vs_" << TwoDPlotInfo.xVar_._name << ": The variable " << TwoDPlotInfo.yVar_._name << " must have and expression, it must be a valid root exprtession. Continuing..." << std::endl;
      continue;
    }
    TwoDPlotInfo.xVar_._minVal = TwoDPlot->getDouble("xminVal", 0);
    TwoDPlotInfo.xVar_._maxVal = TwoDPlot->getDouble("xmaxVal", 0);
    TwoDPlotInfo.yVar_._minVal = TwoDPlot->getDouble("yminVal", 0);
    TwoDPlotInfo.yVar_._maxVal = TwoDPlot->getDouble("ymaxVal", 0);
    if(TwoDPlotInfo.xVar_._maxVal - TwoDPlotInfo.xVar_._minVal <= 0)
    {
      std::cout << TwoDPlotInfo.yVar_._name << "_Vs_" << TwoDPlotInfo.xVar_._name << ": maxVal and minVal for " << TwoDPlotInfo.xVar_._name << " must be specified and define a valid range of values. Continuing..." << std::endl;
      continue;
    }
    if(TwoDPlotInfo.yVar_._maxVal - TwoDPlotInfo.yVar_._minVal <= 0)
    {
      std::cout << TwoDPlotInfo.yVar_._name << "_Vs_" << TwoDPlotInfo.xVar_._name << ": maxVal and minVal for " << TwoDPlotInfo.yVar_._name << " must be specified and define a valid range of values. Continuing..." << std::endl;
      continue;
    }
    TwoDPlotInfo.xVar_._bins = TwoDPlot->getInt("xbins", 0);
    TwoDPlotInfo.yVar_._bins = TwoDPlot->getInt("ybins", 0);
    if(TwoDPlotInfo.xVar_._bins <= 0)
    {
      std::cout << TwoDPlotInfo.yVar_._name << "_Vs_" << TwoDPlotInfo.xVar_._name << ": bins for " << TwoDPlotInfo.xVar_._name << " must be a resonable and valid value. Continuing..." << std::endl;
      continue;
    }
    if(TwoDPlotInfo.yVar_._bins <= 0)
    {
      std::cout << TwoDPlotInfo.yVar_._name << "_Vs_" << TwoDPlotInfo.xVar_._name << ": bins for " << TwoDPlotInfo.yVar_._name << " must be a resonable and valid value. Continuing..." << std::endl;
      continue;
    }
    TwoDPlotInfo.xVar_._label = TwoDPlot->getString("xlabel", "");
    TwoDPlotInfo.yVar_._label = TwoDPlot->getString("ylabel", "");

    retVal.push_back(TwoDPlotInfo);
  }

  return retVal;
}

std::vector<MyVariable> getVariables(JSONWrapper::Object& json)
{
  std::vector<MyVariable> retVal;

  std::vector<JSONWrapper::Object> variables = json["variables"].daughters();
  for(auto variable = variables.begin(); variable != variables.end(); ++variable)
  {
    MyVariable variableInfo;

    variableInfo._name = variable->getString("name", "");
    if(variableInfo._name == "")
    {
      std::cout << "All variables must have names. Continuing..." << std::endl;
      continue;
    }
    variableInfo._expression = variable->getString("expression", "");
    if(variableInfo._expression == "")
    {
      std::cout << variableInfo._name << ": This variable must have an expression, it must be a valid root expression. Continuing..." << std::endl;
      continue;
    }
    variableInfo._minVal = variable->getDouble("minVal", 0);
    variableInfo._maxVal = variable->getDouble("maxVal", 0);
    if(variableInfo._maxVal - variableInfo._minVal <= 0)
    {
      std::cout << variableInfo._name << ": maxVal and minVal must be specified and define a valid range of values. Continuing..." << std::endl;
      continue;
    }
    variableInfo._bins = variable->getInt("bins", 0);
    if(variableInfo._bins <= 0)
    {
      std::cout << variableInfo._name << ": bins must be a resonable and valid value. Continuing..." << std::endl;
      continue;
    }
    variableInfo._label = variable->getString("label", "");

    retVal.push_back(variableInfo);
  }

  return retVal;
}

dataChains getChainsFromJSON(std::string jsonFile, std::string ttreeName, std::string rootDir, std::string customExtension)
{
  dataChains retVal;
  {
    std::vector<ProcessFiles> temp;
    retVal["BG"] = temp;
    retVal["SIG"] = temp;
    retVal["Data"] = temp;
  }

  JSONWrapper::Object json(jsonFile, true);
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

    tempProc.style.lcolor_ = 1;
    tempProc.style.mcolor_ = 1;
    tempProc.style.fcolor_ = 0;
    tempProc.style.lwidth_ = 1;
    tempProc.style.lstyle_ = 1;
    tempProc.style.marker_ = 1;
    if(process->isTag("color"))
    {
      tempProc.style.lcolor_ = process->getInt("color");
      tempProc.style.mcolor_ = tempProc.style.lcolor_;
      tempProc.style.fcolor_ = tempProc.style.lcolor_;
    }
    if(process->isTag("lcolor"))
      tempProc.style.lcolor_ = process->getInt("lcolor");
    if(process->isTag("mcolor"))
      tempProc.style.mcolor_ = process->getInt("mcolor");
    if(process->isTag("fcolor"))
      tempProc.style.fcolor_ = process->getInt("fcolor");
    if(process->isTag("fill"))
      tempProc.style.fcolor_ = process->getInt("fill");
    if(process->isTag("lwidth"))
      tempProc.style.lwidth_ = process->getInt("lwidth");
    if(process->isTag("lstyle"))
      tempProc.style.lstyle_ = process->getInt("lstyle");
    if(process->isTag("marker"))
      tempProc.style.marker_ = process->getInt("marker");

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
      tempSample.chain = new TChain(ttreeName.c_str(), ((*sample).getString("dtag", "") + (*sample).getString("suffix", "")).c_str());
      int nFiles = (*sample).getInt("split", 1);

      for(int file = 0; file < nFiles; ++file)
      {
        std::string segmentExt;
        if(nFiles != 1)
        {
          std::stringstream buf;
          buf << "_" << file;
          buf >> segmentExt;
        }

        std::string fileName = rootDir + "/" + (*sample).getString("dtag", "") + (*sample).getString("suffix", "") + segmentExt + filtExt + customExtension + ".root";

        TFile* file = new TFile(fileName.c_str(), "READONLY");
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

        tempSample.chain->Add(fileName.c_str());
        ++(tempSample.nFiles);
      }

      tempProc.samples.push_back(tempSample);
    }

    retVal[type].push_back(tempProc);
  }

  return retVal;
}

void printHelp()
{
  std::cout << "makeDistributions help - There are the following options (default values, if any, in parenthesis):" << std::endl << std::endl;

  std::cout << "--help             -->  Print this help message" << std::endl;
  std::cout << "--json             -->  Configuration file with the samples" << std::endl;
  std::cout << "--inDir            -->  Path to the directory where to find the input files" << std::endl;
  std::cout << "--outDir           -->  Path to the directory where to output plots and tables (will be created if it doesn't exist, default: ./OUT/)" << std::endl;
  std::cout << "--baseSelection    -->  Base selection to apply to the tree to get the selected events (default: selected)" << std::endl;
  std::cout << "--ttree            -->  Name of the ttree containing the events (default: Events)" << std::endl;
  std::cout << "--customExtension  -->  Custom extension on the files with the ttrees (default: _summary)" << std::endl;
  std::cout << "--plotExt          -->  Extension format with which to save the plots, repeat this command if multiple formats are desired (default: png)" << std::endl;
  std::cout << "--variables        -->  JSON file with the variables which will be plotted" << std::endl;
  std::cout << "--normalize        -->  If the plots are to be normalized" << std::endl;
  std::cout << "--unblind          -->  If the data should be plotted and reported as well" << std::endl;
  std::cout << "--iLumi            -->  The integrated luminosity to be used in the plots" << std::endl;
  std::cout << "--pointVar         -->  To use a multipoint signal sample, define this value with the \"Draw\" option to use to differentiate signal points" << std::endl;
  std::cout << "--signalSelection  -->  When using a multipoint signal sample, extra selection to apply only to signal samples to get the selected events" << std::endl;
  std::cout << "--sigXSec          -->  When using a multipoint signal sample, you must define the signal cross section to use" << std::endl;
  std::cout << "--sigNInitEvents   -->  When using a multipoint signal sample, the number of initial events for each signal point" << std::endl;
  std::cout << "--extraSignal      -->  A second extra selection to apply only to signal" << std::endl;

  std::cout << std::endl << "Example command:" << std::endl << "\tmakeDistributions --json samples.json --outDir ./OUT/ --variables variables.json --inDir /directory" << std::endl;
  return;
}
