#ifndef FRANALYZER_H
#define FRANALYZER_H

#include "TCanvas.h"
#include "TFile.h"
#include "TKey.h"
#include "TBenchmark.h"
#include "TStyle.h"
#include "THStack.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TFrame.h"
#include "TF1.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TObject.h"
#include "TMath.h"
#include "TRandom.h"
#include "TGraphAsymmErrors.h"
#include "TDirectoryFile.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <sys/stat.h>

using namespace std;

class FRAnalyzer{

 public:
  //Constructors
  FRAnalyzer();
  FRAnalyzer(string File1, 
             string File2,
             int    Index,
             int    Ntoy,
             bool   Data,
             string Dir
             );

  //Destructor
  virtual ~FRAnalyzer();

  //Methods
  void CreateDir();
  void ReadFiles();
  void Initialize();
  TPaveText* GetPaveTextCMS();
  TPaveText* GetPaveTextSim(bool simulation);
  TPaveText* GetPaveTextLumi();
  vector<TLatex*> FRlatexText(int i, double x, double y);
  TLatex* ClosureTestLatexText(double x, double y, string channel);
  TH1F* Subtract(TH1F* h1, TH1F* h2);
  vector<TCanvas*> CreateCanvases(int nc, string name);
  vector<TH1F*> GetFRHistograms(bool numerator,string sample);
  vector<TH1F*> GetRatios(bool rebin);
  vector<TH1F*> GetRatiosData();
  vector<TH2F*> BuildShapeHistograms(string s, string cat);
  void Printing(int i,int j, double integral01, double integral10, double integral00, double uCR01, double uCR10, double uCR00);
  void PrintFinalNumbers(int i,int j, vector<TH2F*> shape,vector<TH2F*> shape1,vector<TH2F*> shap2);
  void Estimate(string dir,vector<TH1F*> ratios, 
                double elIso, double muIso, double tauIso, double sumPt, int cut_index,
                vector<TH2F*> h2DH_CR01,  vector<TH2F*> h2DH_CR10,  vector<TH2F*> h2DH_CR00,  vector<TH2F*> h2DH_CR11,
                vector<TH2F*> h2DA_CR01,  vector<TH2F*> h2DA_CR10,  vector<TH2F*> h2DA_CR00,  vector<TH2F*> h2DA_CR11);
  vector<string> GetStringVector(string s);
  vector<TFile*> Overall();
  void SetExtremes(vector<TH2F*> hCR);
  //void SetExtremes(vector<TH2F*> hCR01, vector<TH2F*> hCR10, vector<TH2F*> hCR00);
  vector<TH2F*> SumUpOSandSS(vector<TH2F*> hCR);
  void SetErrorForEmptyBins(vector<TH2F*> hCR01, vector<TH2F*> hCR10, vector<TH2F*> hCR00);
  vector<TH2F*> GetShapes(vector<TH2F*> hCR01, vector<TH2F*> hCR10, vector<TH2F*> hCR00);
  vector<TH2F*> SubtractShapes(vector<TH2F*> hData, vector<TH2F*> hMC);

  vector<TH2F*> GetShapeSignal(vector<TH2F*> hCR11, bool isMain);  
  void RescaleShape( vector<TH2F*> shape);
  void GetFCIntervals(vector<TH2F*> shape);
  vector<TFile*> FinalEstimate();
  void RetrieveResults();
  void RetrieveCTResults();
  void WriteFinalNumbers(vector<TH2F*> h1, vector<TH2F*> h2, vector<TH2F*> h3, vector<TH2F*> h4,vector<TH2F*> h5,vector<TH2F*> h6);
  vector<TH2F*> GetStraightMC();

  TFile* GetToyExperiments(vector<TH2F*>sCR01,vector<TH2F*>sCR10,vector<TH2F*>sCR00,vector<TH2F*>hCR01,vector<TH2F*>hCR10,vector<TH2F*>hCR00,
                vector<TH2F*>sDataFit,vector<TH2F*>sDataHis,vector<TH2F*>sMC,bool useBayesian);
  pair<double,double> LogNormalFit(TH1 *h1);
  vector<double> GetFRWeights();
  void DefineLikelihoodFunction();
  TF1* GetLikelihoodFunction(int count);
  vector<TH2F*> BuildSummaryHistograms(string s);
  double AutomatizeBinning(double content, double error, bool max);
  TF1* GetFitFunction(); 
  pair<double,double> GetClosureTestInfo(vector<TH2F*> shapeDD, vector<TH2F*> shapeMC, int ch1, int ch2);
  TGraphErrors* GetClosureTestGraph(pair<double,double> info);
  TFile* ClosureTest();


 protected:
  string m_File1;
  string m_File2;
  int    m_Index;
  int    m_Ntoy;
  bool   m_Data;
  bool   debug;
  string m_Dir;
  
  TFile* file1;
  TFile* file2;

  TF1* f1;
  TF1* f01;
  TF1* f02;
  TF1* f03;
  TF1* f04;
  TF1* f05;
  TF1* f06;
  TF1* f07;
  TF1* f08;
  TF1* f09;
  TF1* f10;
  TF1* f11;
  TF1* fitFun;
  TF1* tmpFitFun;

  stringstream cidx; string Idx;
  vector<string> shapeName;
  vector<string> finalStates;
  vector<string> forName;
  vector<string> channels;
  vector<string> cutIndexLabel;
  vector<string> hlabel;
  string sample;
  string dir_data;
  string dir_wz;
  string dir_ww;
  string dir_vg;
  string dir_ttZ;
  string dir_top;
  string dir_qcd;
  string dir_qcdmu;
  string dir_wjets;
  string dir_dy;
  string dir_zz2l2tau;
  string dir_zz2l2nu;
  string dir_zz4l;
  string dir_wwz;
  string dir_zzz;
  string dir_wzz;
  string dir_www;
  string dir_zh;
  
  vector<double> tauIso;
  vector<double> eleIso;
  vector<double> muoIso;

  TH1F* Hbins;
  TH2F* hBins;
  TH2F* h2ref;
  int   nbx;
  int   maxx;
  int   minx;
  int   nbins; 
  TH2F* h_eeemOS; TH2F* h_eeemSS;
  TH2F* h_eeetOS; TH2F* h_eeetSS;
  TH2F* h_eemtOS; TH2F* h_eemtSS;
  TH2F* h_eettOS; TH2F* h_eettSS;
  TH2F* h_mmemOS; TH2F* h_mmemSS;
  TH2F* h_mmetOS; TH2F* h_mmetSS;
  TH2F* h_mmmtOS; TH2F* h_mmmtSS;
  TH2F* h_mmttOS; TH2F* h_mmttSS;

};
#endif
