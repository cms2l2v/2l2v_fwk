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
class My2DPlot;

typedef std::vector<std::pair<std::string,std::vector<std::pair<int,TChain*>>>> fileChains;

void printHelp();
std::vector<MyVariable> getVariables(JSONWrapper::Object& json);
std::vector<My2DPlot> get2DPlots(JSONWrapper::Object& json);
fileChains getChainsFromJSON(JSONWrapper::Object& json, std::string RootDir, std::string type, std::string treename, std::string customExtension);
TH1D* getHist(fileChains files,  TCut cut, MyVariable variable, bool correctFiles = true);
THStack* getStack(fileChains files,  TCut cut, MyVariable variable, bool correctFiles = true);
void make2D(std::string outDir, std::vector<std::string> plotExt, fileChains files, TCut cut, My2DPlot TwoDPlot, bool correctFiles = true);
std::string cleanString(std::string str, std::string illegal, char replacement);


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

  inline int marker(){return marker_;};
  inline int lcolor(){return lcolor_;};
  inline int mcolor(){return mcolor_;};
  inline int fcolor(){return fcolor_;};
  inline int lwidth(){return lwidth_;};
  inline int lstyle(){return lstyle_;};

  friend fileChains getChainsFromJSON(JSONWrapper::Object& json, std::string RootDir, std::string type, std::string treename, std::string customExtension);

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

    if(arg.find("--variables") != std::string::npos)
    {
      variablesFile = argv[i+1];
      // TODO: Check if the file exists
      ++i;
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
  auto  BG_samples = getChainsFromJSON(json, inDir,  "BG", ttree_name, customExtension);
  auto SIG_samples = getChainsFromJSON(json, inDir, "SIG", ttree_name, customExtension);

  TCut BGCut  = baseSelection.c_str();
  TCut SIGCut = BGCut && signalSelection.c_str();

  TCanvas c1("c1", "c1", 800, 600);
  c1.SetLogy();
  gStyle->SetOptStat(0);
  for(auto variable = variables.begin(); variable != variables.end(); ++variable)
  {
    std::cout << "Processing " << variable->name() << std::endl;
    std::string label = variable->label();
    if(label == "")
      label = variable->name();
    THStack* final_bg  = getStack( BG_samples,   BGCut, *variable);
    final_bg->SetNameTitle("Background", ("Background;"+label+";% Events").c_str());
    TH1D* final_sig = getHist(SIG_samples,  SIGCut, *variable, false);
    final_sig->SetNameTitle("Signal", ("Signal;"+label+";% Events").c_str());

    //final_bg->Scale(1/final_bg->Integral());
    final_sig->Scale(1/final_sig->Integral());

    //final_bg->SetLineColor(kBlue);
    //final_sig->SetLineColor(kRed);
    final_bg->Draw();
    final_bg->SetMinimum(10e-6);
    double max = final_sig->GetMaximum();
    double min = final_sig->GetMinimum();
    if(final_bg->GetMaximum() > max)
      max = final_bg->GetMaximum();
    if(final_bg->GetMinimum() < min)
      min = final_bg->GetMinimum();

    if(extraSignal != "")
    {
      TCut SIG2Cut = BGCut && extraSignal.c_str();
      TH1D* sig_2 = getHist(SIG_samples,  SIG2Cut, *variable, false);
      sig_2->SetNameTitle("ExtraSignal", ("ExtraSignal;"+label+";% Events").c_str());
      sig_2->Scale(1/sig_2->Integral());
      sig_2->SetLineColor(kBlue);
      if(sig_2->GetMaximum() > max)
        max = sig_2->GetMaximum();
      if(sig_2->GetMinimum() < min)
        min = sig_2->GetMinimum();
      final_bg->SetMaximum(max);
//      final_bg->SetMinimum(min);
      final_sig->Draw("same");
      sig_2->Draw("same");
    }
    else
    {
      final_bg->SetMaximum(max);
//      final_bg->SetMinimum(min);
      final_sig->Draw("same");
    }

    c1.BuildLegend(0.88, 0.67, 1, 1);

    for(auto ext = plotExt.begin(); ext != plotExt.end(); ++ext)
      c1.SaveAs((outDir + variable->name() + *ext).c_str());

    delete final_bg;
    delete final_sig;
  }

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

THStack* getStack(fileChains files,  TCut cut, MyVariable variable, bool correctFiles)
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

fileChains getChainsFromJSON(JSONWrapper::Object& json, std::string RootDir, std::string type, std::string treename, std::string customExtension)
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

void printHelp()
{
  std::cout << "makeDistributions help - There are the following options (default values, if any, in parenthesis):" << std::endl << std::endl;

  std::cout << "--help             -->  Print this help message" << std::endl;
  std::cout << "--json             -->  Configuration file with the samples" << std::endl;
  std::cout << "--inDir            -->  Path to the directory where to find the input files" << std::endl;
  std::cout << "--outDir           -->  Path to the directory where to output plots and tables (will be created if it doesn't exist, default: ./OUT/)" << std::endl;
  std::cout << "--baseSelection    -->  Base selection to apply to the tree to get the selected events (default: selected)" << std::endl;
  std::cout << "--signalSelection  -->  Extra selection to apply only to signal samples to get the selected events" << std::endl;
  std::cout << "--ttree            -->  Name of the ttree containing the events (default: Events)" << std::endl;
  std::cout << "--customExtension  -->  Custom extension on the files with the ttrees (default: _summary)" << std::endl;
  std::cout << "--plotExt          -->  Extension format with which to save the plots, repeat this command if multiple formats are desired (default: png)" << std::endl;
  std::cout << "--variables        -->  JSON file with the variables which will be plotted" << std::endl;
  std::cout << "--extraSignal      -->  A second extra selection to apply only to signal" << std::endl;

  std::cout << std::endl << "Example command:" << std::endl << "\tmakeDistributions --json samples.json --outDir ./OUT/ --variables variables.json --inDir /directory" << std::endl;
  return;
}
