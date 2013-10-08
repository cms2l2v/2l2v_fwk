// Original Author:  Loic Quertenmont

#include <iostream>
#include <boost/shared_ptr.hpp>
#include "Math/GenVector/Boost.h"

#include "UserCode/llvv_fwk/src/tdrstyle.C"
#include "UserCode/llvv_fwk/src/JSONWrapper.cc"
#include "UserCode/llvv_fwk/interface/RootUtils.h"
#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "HiggsAnalysis/CombinedLimit/interface/th1fmorph.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TString.h"
#include "TList.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TObjArray.h"
#include "THStack.h"
#include "TGraphErrors.h"

#include<iostream>
#include<fstream>
#include<map>
#include<algorithm>
#include<vector>
#include<set>



using namespace std;
double NonResonnantSyst = 0.25;
double GammaJetSyst = 0.25;


TString signalSufix="";
TString histo(""), histoVBF("");
int rebinVal = 1;
double MCRescale = 1.0;
int mass;
bool shape = false;
TString postfix="";
TString systpostfix="";
bool runSystematics = false; 
std::vector<TString> Channels;
std::vector<TString> AnalysisBins;
double DDRescale = 1.0;
TString DYFile ="";

bool BackExtrapol  = false;
bool subNRB        = false;
bool MCclosureTest = false;

bool mergeWWandZZ = false;
bool skipWW = true;
bool fast = false;
bool skipGGH = false;
bool skipQQH = false;
bool subDY = false;
bool subWZ = false;
bool blindData = false;
bool blindWithSignal = false; 
TString inFileUrl(""),jsonFile("");
double shapeMin = 0;
double shapeMax = 9999;
double shapeMinVBF = 0;
double shapeMaxVBF = 9999;
bool doInterf = false;
double minSignalYield = 0;

bool dirtyFix1 = false;
bool dirtyFix2 = false;

int indexvbf = -1;
int indexcut   = -1, indexcutL=-1, indexcutR=-1;
int massL=-1, massR=-1;

void setTGraph(TString proc, TString suffix);
void initializeTGraph();
TGraph *ggH7TG_xsec=NULL, *ggH7TG_errp=NULL, *ggH7TG_errm=NULL, *ggH7TG_scap=NULL, *ggH7TG_scam=NULL, *ggH7TG_pdfp=NULL, *ggH7TG_pdfm=NULL;
TGraph *qqH7TG_xsec=NULL, *qqH7TG_errp=NULL, *qqH7TG_errm=NULL, *qqH7TG_scap=NULL, *qqH7TG_scam=NULL, *qqH7TG_pdfp=NULL, *qqH7TG_pdfm=NULL;
TGraph *ggH8TG_xsec=NULL, *ggH8TG_errp=NULL, *ggH8TG_errm=NULL, *ggH8TG_scap=NULL, *ggH8TG_scam=NULL, *ggH8TG_pdfp=NULL, *ggH8TG_pdfm=NULL;
TGraph *qqH8TG_xsec=NULL, *qqH8TG_errp=NULL, *qqH8TG_errm=NULL, *qqH8TG_scap=NULL, *qqH8TG_scam=NULL, *qqH8TG_pdfp=NULL, *qqH8TG_pdfm=NULL;
TGraph *    TG_xsec=NULL, *    TG_errp=NULL, *    TG_errm=NULL, *    TG_scap=NULL, *    TG_scam=NULL, *    TG_pdfp=NULL, *    TG_pdfm=NULL;

TGraph* TG_QCDScaleK0ggH0=NULL, *TG_QCDScaleK0ggH1=NULL, *TG_QCDScaleK1ggH1=NULL, *TG_QCDScaleK1ggH2=NULL, *TG_QCDScaleK2ggH2=NULL;
TGraph* TG_UEPSf0=NULL, *TG_UEPSf1=NULL, *TG_UEPSf2=NULL;



//wrapper for a projected shape for a given proc
class ShapeData_t
{
  public:
  std::map<string, double> uncScale;
  std::map<string, TH1*  > uncShape;
  TH1* fit;

  ShapeData_t(){
	  fit=NULL;
  }
  ~ShapeData_t(){}

  TH1* histo(){
     if(uncShape.find("")==uncShape.end())return NULL;
     return uncShape[""];
  }

  void makeStatUnc(string prefix="", string suffix=""){
     TH1* h = histo();
     if(!h || h->Integral()<=0)return;
     TH1* statU=(TH1 *)h->Clone(TString(h->GetName())+"StatU");
     TH1* statD=(TH1 *)h->Clone(TString(h->GetName())+"StatD");
     for(int ibin=1; ibin<=statU->GetXaxis()->GetNbins(); ibin++){
        statU->SetBinContent(ibin,std::min(2*h->GetBinContent(ibin), std::max(0.01*h->GetBinContent(ibin), statU->GetBinContent(ibin) + statU->GetBinError(ibin))));
        statD->SetBinContent(ibin,std::min(2*h->GetBinContent(ibin), std::max(0.01*h->GetBinContent(ibin), statD->GetBinContent(ibin) - statD->GetBinError(ibin))));
     }
     uncShape[prefix+"stat"+suffix+"Up"  ] = statU;
     uncShape[prefix+"stat"+suffix+"Down"] = statD;
  }

  double getScaleUncertainty(){
     double Total=0;
     for(std::map<string, double>::iterator unc=uncScale.begin();unc!=uncScale.end();unc++){
        if(unc->second<0)continue;
        Total+=pow(unc->second,2);
     }     
     return Total>0?sqrt(Total):-1;
  }

};

class ChannelInfo_t
{
  public:
        string bin;
        string channel;
	std::map<string, ShapeData_t> shapes;

	ChannelInfo_t(){}
	~ChannelInfo_t(){}
};

class ProcessInfo_t
{
  public:
  bool isData;
  bool isBckg;
  bool isSign;
  double xsec;
  double br;
  double mass;
  string shortName;
  std::map<string, ChannelInfo_t> channels;
  JSONWrapper::Object jsonObj;

  ProcessInfo_t(){xsec=0;}
  ~ProcessInfo_t(){}
};

class AllInfo_t
{
  public:
	std::map<string, ProcessInfo_t> procs;
        std::vector<string> sorted_procs;

	AllInfo_t(){};
	~AllInfo_t(){};

        // reorder the procs to get the backgrounds; total bckg, signal, data 
        void sortProc();

        // Sum up all background processes and add this as a total process
        void addProc(ProcessInfo_t& dest, ProcessInfo_t& src);

        // Sum up all background processes and add this as a total process
        void computeTotalBackground();

        // Replace the Data process by TotalBackground
        void blind();

        // Print the Yield table
         void getYieldsFromShape(FILE* pFile, std::vector<TString>& selCh, string histoName);

        // Dump efficiencies
        void getEffFromShape(FILE* pFile, std::vector<TString>& selCh, string histoName);

        // drop background process that have a negligible yield
        void dropSmallBckgProc(std::vector<TString>& selCh, string histoName, double threshold);

        // drop control channels
        void dropCtrlChannels(std::vector<TString>& selCh);

        // Make a summary plot
        void showShape(std::vector<TString>& selCh , TString histoName, TString SaveName);

        // Make a summary plot
        void saveHistoForLimit(string histoName, TFile* fout);
     
        // produce the datacards 
        void buildDataCards(string histoName, TString url);

        // Load histograms from root file and json to memory
        void getShapeFromFile(TFile* inF, std::vector<string> channelsAndShapes, int cutBin, JSONWrapper::Object &Root,  double minCut=0, double maxCut=9999, bool onlyData=false);

        // replace MC NonResonnant Backgrounds by DataDriven estimate
        void doBackgroundSubtraction(FILE* pFile, std::vector<TString>& selCh,TString ctrlCh,TString mainHisto, TString sideBandHisto);

        // replace MC Z+Jets Backgrounds by DataDriven Gamma+Jets estimate
        void doDYReplacement(FILE* pFile, std::vector<TString>& selCh,TString ctrlCh,TString mainHisto);

        // Rebin histograms to make sure that high mt/met region have no empty bins
        void rebinMainHisto(string histoName);

        // Interpollate the signal sample between two mass points 
        void SignalInterpolation(string histoName);

        // Rescale signal sample for the effect of the interference and propagate the uncertainty 
        void RescaleForInterference(string histoName);
};


void printHelp();
void printHelp()
{
  printf("Options\n");
  printf("--in        --> input file with from plotter\n");
  printf("--json      --> json file with the sample descriptor\n");
  printf("--histoVBF  --> name of histogram to be used for VBF\n");
  printf("--histo     --> name of histogram to be used\n");
  printf("--shapeMin  --> left cut to apply on the shape histogram\n");
  printf("--shapeMax  --> right cut to apply on the shape histogram\n");
  printf("--shapeMinVBF  --> left cut to apply on the shape histogram for Vbf bin\n");
  printf("--shapeMaxVBF  --> right cut to apply on the shape histogram for Vbf bin\n");
  printf("--indexvbf  --> index of selection to be used for the vbf bin (if unspecified same as --index)\n");
  printf("--index     --> index of selection to be used (Xbin in histogram to be used)\n");
  printf("--indexL    --> index of selection to be used (Xbin in histogram to be used) used for interpolation\n");
  printf("--indexR    --> index of selection to be used (Xbin in histogram to be used) used for interpolation\n");
  printf("--m         --> higgs mass to be considered\n");
  printf("--mL        --> higgs mass on the left  of the mass to be considered (used for interpollation\n");
  printf("--mR        --> higgs mass on the right of the mass to be considered (used for interpollation\n");
  printf("--syst      --> use this flag if you want to run systematics, default is no systematics\n");
  printf("--shape     --> use this flag if you want to run shapeBased analysis, default is cut&count\n");
  printf("--subNRB    --> use this flag if you want to subtract non-resonant-backgounds similarly to what was done in 2011 (will also remove H->WW)\n");
  printf("--subNRB12  --> use this flag if you want to subtract non-resonant-backgounds using a new technique that keep H->WW\n");
  printf("--subDY     --> histogram that contains the Z+Jets background estimated from Gamma+Jets)\n");
  printf("--subWZ     --> use this flag if you want to subtract WZ background by the 3rd lepton SB)\n");
  printf("--DDRescale --> factor to be used in order to multiply/rescale datadriven estimations\n");
  printf("--closure   --> use this flag if you want to perform a MC closure test (use only MC simulation)\n");
  printf("--bins      --> list of bins to be used (they must be comma separated without space)\n");
  printf("--HWW       --> use this flag to consider HWW signal)\n");
  printf("--skipGGH   --> use this flag to skip GGH signal)\n");
  printf("--skipQQH   --> use this flag to skip GGH signal)\n");
  printf("--blind     --> use this flag to replace observed data by total predicted background)\n");
  printf("--blindWithSignal --> use this flag to replace observed data by total predicted background+signal)\n");
  printf("--fast      --> use this flag to only do assymptotic prediction (very fast but inaccurate))\n");
  printf("--postfix    --> use this to specify a postfix that will be added to the process names)\n");
  printf("--systpostfix    --> use this to specify a syst postfix that will be added to the process names)\n");
  printf("--MCRescale    --> use this to specify a syst postfix that will be added to the process names)\n");
  printf("--interf     --> use this to rescale xsection according to WW interferences)\n");
  printf("--minSignalYield   --> use this to specify the minimum Signal yield you want in each channel)\n");
  printf("--signalSufix --> use this flag to specify a suffix string that should be added to the signal 'histo' histogram\n");
  printf("--rebin         --> rebin the histogram\n");
  printf("--BackExtrapol --> extrapollate background histograms to high mass/met\n");
}

//
int main(int argc, char* argv[])
{
  setTDRStyle();
  gStyle->SetPadTopMargin   (0.06);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadRightMargin (0.16);
  gStyle->SetPadLeftMargin  (0.14);
  gStyle->SetTitleSize(0.04, "XYZ");
  gStyle->SetTitleXOffset(1.1);
  gStyle->SetTitleYOffset(1.45);
  gStyle->SetPalette(1);
  gStyle->SetNdivisions(505);
  gStyle->SetOptStat(0);  
  gStyle->SetOptFit(0);

  //init the TGraphs
  initializeTGraph();

  //get input arguments
  for(int i=1;i<argc;i++){
    string arg(argv[i]);
    if(arg.find("--help")          !=string::npos) { printHelp(); return -1;} 
    else if(arg.find("--minSignalYield") !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%lf",&minSignalYield ); i++; printf("minSignalYield = %f\n", minSignalYield);}
    else if(arg.find("--subNRB")   !=string::npos) { subNRB=true; skipWW=true; printf("subNRB = True\n");}
    else if(arg.find("--subDY")    !=string::npos) { subDY=true; DYFile=argv[i+1];  i++; printf("Z+Jets will be replaced by %s\n",DYFile.Data());}
    else if(arg.find("--subWZ")    !=string::npos) { subWZ=true; printf("WZ will be estimated from 3rd lepton SB\n");}
    else if(arg.find("--DDRescale")!=string::npos && i+1<argc)  { sscanf(argv[i+1],"%lf",&DDRescale); i++;}
    else if(arg.find("--MCRescale")!=string::npos && i+1<argc)  { sscanf(argv[i+1],"%lf",&MCRescale); i++;}
    else if(arg.find("--HWW")      !=string::npos) { skipWW=false; printf("HWW = True\n");}
    else if(arg.find("--skipGGH")  !=string::npos) { skipGGH=true; printf("skipGGH = True\n");}
    else if(arg.find("--skipQQH")  !=string::npos) { skipQQH=true; printf("skipQQH = True\n");}
    else if(arg.find("--blindWithSignal")    !=string::npos) { blindData=true; blindWithSignal=true; printf("blindData = True; blindWithSignal = True\n");}
    else if(arg.find("--blind")    !=string::npos) { blindData=true; printf("blindData = True\n");}
    else if(arg.find("--closure")  !=string::npos) { MCclosureTest=true; printf("MCclosureTest = True\n");}
    else if(arg.find("--shapeMinVBF") !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%lf",&shapeMinVBF); i++; printf("Min cut on shape for VBF = %f\n", shapeMinVBF);}
    else if(arg.find("--shapeMaxVBF") !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%lf",&shapeMaxVBF); i++; printf("Max cut on shape for VBF = %f\n", shapeMaxVBF);}
    else if(arg.find("--shapeMin") !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%lf",&shapeMin); i++; printf("Min cut on shape = %f\n", shapeMin);}
    else if(arg.find("--shapeMax") !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%lf",&shapeMax); i++; printf("Max cut on shape = %f\n", shapeMax);}
    else if(arg.find("--interf")    !=string::npos) { doInterf=true; printf("doInterf = True\n");}
    else if(arg.find("--indexvbf") !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%i",&indexvbf); i++; printf("indexVBF = %i\n", indexvbf);}
    else if(arg.find("--indexL")    !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%i",&indexcutL); i++; printf("indexL = %i\n", indexcutL);}
    else if(arg.find("--indexR")    !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%i",&indexcutR); i++; printf("indexR = %i\n", indexcutR);}
    else if(arg.find("--index")    !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%i",&indexcut); i++; printf("index = %i\n", indexcut);}
    else if(arg.find("--in")       !=string::npos && i+1<argc)  { inFileUrl = argv[i+1];  i++;  printf("in = %s\n", inFileUrl.Data());  }
    else if(arg.find("--json")     !=string::npos && i+1<argc)  { jsonFile  = argv[i+1];  i++;  printf("json = %s\n", jsonFile.Data()); }
    else if(arg.find("--histoVBF") !=string::npos && i+1<argc)  { histoVBF  = argv[i+1];  i++;  printf("histoVBF = %s\n", histoVBF.Data()); }
    else if(arg.find("--histo")    !=string::npos && i+1<argc)  { histo     = argv[i+1];  i++;  printf("histo = %s\n", histo.Data()); }
    else if(arg.find("--mL")       !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%i",&massL ); i++; printf("massL = %i\n", massL);}
    else if(arg.find("--mR")       !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%i",&massR ); i++; printf("massR = %i\n", massR);}
    else if(arg.find("--m")        !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%i",&mass ); i++; printf("mass = %i\n", mass);}
    else if(arg.find("--bins")     !=string::npos && i+1<argc)  { char* pch = strtok(argv[i+1],",");printf("bins are : ");while (pch!=NULL){printf(" %s ",pch); AnalysisBins.push_back(pch);  pch = strtok(NULL,",");}printf("\n"); i++; }
    else if(arg.find("--channels") !=string::npos && i+1<argc)  { char* pch = strtok(argv[i+1],",");printf("channels are : ");while (pch!=NULL){printf(" %s ",pch); Channels.push_back(pch);  pch = strtok(NULL,",");}printf("\n"); i++; }
    else if(arg.find("--fast")     !=string::npos) { fast=true; printf("fast = True\n");}
    else if(arg.find("--postfix")   !=string::npos && i+1<argc)  { postfix = argv[i+1]; systpostfix = argv[i+1]; i++;  printf("postfix '%s' will be used\n", postfix.Data());  }
    else if(arg.find("--systpostfix")   !=string::npos && i+1<argc)  { systpostfix = argv[i+1];  i++;  printf("systpostfix '%s' will be used\n", systpostfix.Data());  }
    else if(arg.find("--syst")     !=string::npos) { runSystematics=true; printf("syst = True\n");}
    else if(arg.find("--shape")    !=string::npos) { shape=true; printf("shapeBased = True\n");}
    else if(arg.find("--dirtyFix2")    !=string::npos) { dirtyFix2=true; printf("dirtyFix2 = True\n");}
    else if(arg.find("--dirtyFix1")    !=string::npos) { dirtyFix1=true; printf("dirtyFix1 = True\n");}
    else if(arg.find("--signalSufix") !=string::npos) { signalSufix = argv[i+1]; i++; printf("signalSufix '%s' will be used\n", signalSufix.Data()); }
    else if(arg.find("--rebin")    !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%i",&rebinVal); i++; printf("rebin = %i\n", rebinVal);}
    else if(arg.find("--BackExtrapol")    !=string::npos) { BackExtrapol=true; printf("BackExtrapol = True\n");}
  }
  if(jsonFile.IsNull() || inFileUrl.IsNull() || histo.IsNull() || indexcut == -1 || mass==-1) { printHelp(); return -1; }
  if(AnalysisBins.size()==0)AnalysisBins.push_back("");
  if(Channels.size()==0){Channels.push_back("ee");Channels.push_back("mumu");}


///////////////////////////////////////////////


  //init the json wrapper
  JSONWrapper::Object Root(jsonFile.Data(), true);

  //init globalVariables
  TString massStr(""); if(mass>0)massStr += mass;
  std::vector<TString> allCh,allProcs;


  TString ch[]={"mumu","ee","emu"};
  const size_t nch=sizeof(ch)/sizeof(TString);
  std::vector<TString> sh;
  sh.push_back(histo);
  if(subNRB)sh.push_back(histo+"_NRBctrl");
  if(subWZ)sh.push_back(histo+"_3rdLepton");

  std::vector<string> channelsAndShapes;
  std::vector<string> vbfchannelsAndShapes;

  const size_t nsh=sh.size();
  for(size_t i=0; i<nch; i++){
     for(size_t b=0; b<AnalysisBins.size(); b++){
       int indexcut_ = indexcut; double cutMin=shapeMin; double cutMax=shapeMax;
        for(size_t j=0; j<nsh; j++){
           if(indexvbf>=0 && AnalysisBins[b].Contains("vbf")){
              vbfchannelsAndShapes.push_back((ch[i]+TString(";")+AnalysisBins[b]+TString(";")+sh[j]).Data());
           }else{
              channelsAndShapes.push_back((ch[i]+TString(";")+AnalysisBins[b]+TString(";")+sh[j]).Data());
           }
        }
     }
  }

  AllInfo_t allInfo;

  //open input file
  TFile* inF = TFile::Open(inFileUrl);
  if( !inF || inF->IsZombie() ){ printf("Invalid file name : %s\n", inFileUrl.Data());}
  gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE
  allInfo.getShapeFromFile(inF,    channelsAndShapes, indexcut, Root, shapeMin   , shapeMax   );
  allInfo.getShapeFromFile(inF, vbfchannelsAndShapes, indexvbf, Root, shapeMinVBF, shapeMaxVBF);
  inF->Close();
  printf("Loading all shapes... Done\n");

  allInfo.computeTotalBackground();
  if(MCclosureTest)allInfo.blind();

  FILE* pFile;

  //define vector for search
  std::vector<TString>& selCh = Channels;

  //remove the non-resonant background from data
  if(subNRB){
  pFile = fopen("NonResonnant.tex","w");
  allInfo.doBackgroundSubtraction(pFile, selCh,"emu",histo,histo+"_NRBctrl");
  fclose(pFile);
  }

  //replace Z+Jet background by Gamma+Jet estimates
  if(subDY){
  pFile = fopen("GammaJets.tex","w");
  allInfo.doDYReplacement(pFile, selCh,"gamma",histo);
  fclose(pFile);
  }

  //replace data by total MC background
  if(blindData)allInfo.blind();

  //interpollate signal sample if desired mass point is not available
  allInfo.SignalInterpolation(histo.Data());

  //rescale for interference
  if(doInterf || signalSufix!="")allInfo.RescaleForInterference(histo.Data());

  //extrapolate backgrounds toward higher mt/met region to make sure that there is no empty bins
  if(shape && BackExtrapol)allInfo.rebinMainHisto(histo.Data());

  //drop backgrounds with rate<1%
  allInfo.dropSmallBckgProc(selCh, histo.Data(), 0.01);

  //drop control channels
  allInfo.dropCtrlChannels(selCh);

  //print event yields from the mt shapes
  pFile = fopen("Yields.tex","w");
  allInfo.getYieldsFromShape(pFile, selCh, histo.Data());
  fclose(pFile);

  //print signal efficiency
  pFile = fopen("Efficiency.tex","w");
  allInfo.getEffFromShape(pFile, selCh, histo.Data());
  fclose(pFile);

  //produce a plot
  allInfo.showShape(selCh,histo,"plot");

  //prepare the output
  string limitFile=("hzz2l2v_"+massStr+systpostfix+".root").Data();
  TFile *fout=TFile::Open(limitFile.c_str(),"recreate");

  allInfo.saveHistoForLimit(histo.Data(), fout);

  allInfo.buildDataCards(histo.Data(), limitFile);

  //all done
  fout->Close();

}



void setTGraph(TString proc, TString suffix){
         if( suffix.Contains('8') &&  proc.Contains("qq")){ TG_xsec=qqH8TG_xsec;    TG_errp=qqH8TG_errp;    TG_errm=qqH8TG_errm;    TG_scap=qqH8TG_scap;    TG_scam=qqH8TG_scam;    TG_pdfp=qqH8TG_pdfp;    TG_pdfm=qqH8TG_pdfm;
   }else if( suffix.Contains('8') && !proc.Contains("qq")){ TG_xsec=ggH8TG_xsec;    TG_errp=ggH8TG_errp;    TG_errm=ggH8TG_errm;    TG_scap=ggH8TG_scap;    TG_scam=ggH8TG_scam;    TG_pdfp=ggH8TG_pdfp;    TG_pdfm=ggH8TG_pdfm;
   }else if(!suffix.Contains('8') &&  proc.Contains("qq")){ TG_xsec=qqH7TG_xsec;    TG_errp=qqH7TG_errp;    TG_errm=qqH7TG_errm;    TG_scap=qqH7TG_scap;    TG_scam=qqH7TG_scam;    TG_pdfp=qqH7TG_pdfp;    TG_pdfm=qqH7TG_pdfm;
   }else if(!suffix.Contains('8') && !proc.Contains("qq")){ TG_xsec=ggH7TG_xsec;    TG_errp=ggH7TG_errp;    TG_errm=ggH7TG_errm;    TG_scap=ggH7TG_scap;    TG_scam=ggH7TG_scam;    TG_pdfp=ggH7TG_pdfp;    TG_pdfm=ggH7TG_pdfm;
   }
}
void initializeTGraph(){
   double ggH7_mass [] = {90.0,95.0,100.0,105.0,110.0,110.5,111.0,111.5,112.0,112.5,113.0,113.5,114.0,114.5,115.0,115.5,116.0,116.5,117.0,117.5,118.0,118.5,119.0,119.5,120.0,120.5,121.0,121.5,122.0,122.5,123.0,123.5,124.0,124.5,125.0,125.5,126.0,126.5,127.0,127.5,128.0,128.5,129.0,129.5,130.0,130.5,131.0,131.5,132.0,132.5,133.0,133.5,134.0,134.5,135.0,135.5,136.0,136.5,137.0,137.5,138.0,138.5,139.0,139.5,140.0,141.0,142.0,143.0,144.0,145.0,146.0,147.0,148.0,149.0,150.0,151.0,152.0,153.0,154.0,155.0,156.0,157.0,158.0,159.0,160.0,162.0,164.0,166.0,168.0,170.0,172.0,174.0,176.0,178.0,180.0,182.0,184.0,186.0,188.0,190.0,192.0,194.0,196.0,198.0,200.0,202.0,204.0,206.0,208.0,210.0,212.0,214.0,216.0,218.0,220.0,222.0,224.0,226.0,228.0,230.0,232.0,234.0,236.0,238.0,240.0,242.0,244.0,246.0,248.0,250.0,252.0,254.0,256.0,258.0,260.0,262.0,264.0,266.0,268.0,270.0,272.0,274.0,276.0,278.0,280.0,282.0,284.0,286.0,288.0,290.0,295.0,300.0,305.0,310.0,315.0,320.0,325.0,330.0,335.0,340.0,345.0,350.0,360.0,370.0,380.0,390.0,400.0,420.0,440.0,460.0,480.0,500.0,520.0,540.0,560.0,580.0,600.0,620.0,640.0,660.0,680.0,700.0,720.0,740.0,760.0,780.0,800.0,820.0,840.0,860.0,880.0,900.0,920.0,940.0,960.0,980.0,1000.0};
   double ggH7_xsec [] = {29.51,26.51,24.00,21.77,19.84,19.66,19.48,19.31,19.13,18.96,18.79,18.63,18.46,18.30,18.14,17.98,17.83,17.67,17.52,17.37,17.22,17.08,16.93,16.79,16.65,16.51,16.37,16.23,16.10,15.97,15.84,15.71,15.58,15.45,15.32,15.20,15.08,14.96,14.85,14.73,14.62,14.50,14.38,14.27,14.16,14.05,13.94,13.83,13.72,13.62,13.51,13.41,13.31,13.21,13.11,13.01,12.91,12.81,12.72,12.62,12.53,12.44,12.35,12.26,12.18,12.00,11.82,11.65,11.49,11.33,11.18,11.02,10.87,10.72,10.58,10.43,10.29,10.16,10.02,9.886,9.754,9.624,9.487,9.349,9.202,8.830,8.519,8.246,8.009,7.786,7.578,7.389,7.212,7.041,6.869,6.696,6.522,6.349,6.179,6.017,5.865,5.725,5.598,5.483,5.377,5.277,5.188,5.106,5.009,4.922,4.833,4.758,4.695,4.608,4.528,4.449,4.381,4.321,4.245,4.177,4.114,4.056,3.990,3.924,3.854,3.789,3.726,3.667,3.611,3.555,3.501,3.449,3.398,3.349,3.301,3.255,3.211,3.167,3.125,3.083,3.044,3.006,2.970,2.934,2.900,2.866,2.833,2.803,2.773,2.744,2.677,2.616,2.563,2.516,2.478,2.443,2.418,2.403,2.398,2.407,2.431,2.428,2.408,2.362,2.283,2.175,2.049,1.776,1.507,1.263,1.050,0.8708,0.7211,0.5976,0.4960,0.4126,0.3444,0.2883,0.2422,0.2042,0.1728,0.1468,0.1252,0.1071,0.09190,0.07930,0.06850,0.05950,0.05180,0.04520,0.03970,0.03480,0.03070,0.02710,0.02410,0.02140,0.01900};
   double ggH7_errp [] = {+16.0,+15.8,+15.5,+15.4,+15.2,+15.2,+15.2,+15.2,+15.2,+15.1,+15.1,+15.1,+15.1,+15.1,+15.1,+15.1,+15.1,+15.0,+15.0,+15.0,+14.9,+14.9,+14.9,+14.8,+14.8,+14.8,+14.8,+14.8,+14.7,+14.7,+14.7,+14.7,+14.7,+14.7,+14.7,+14.7,+14.7,+14.7,+14.7,+14.7,+14.6,+14.6,+14.6,+14.6,+14.6,+14.6,+14.6,+14.6,+14.6,+14.6,+14.5,+14.5,+14.5,+14.5,+14.5,+14.5,+14.5,+14.5,+14.5,+14.5,+14.4,+14.4,+14.4,+14.4,+14.4,+14.4,+14.4,+14.3,+14.3,+14.3,+14.3,+14.3,+14.2,+14.2,+14.2,+14.2,+14.1,+14.1,+14.1,+14.0,+14.0,+14.0,+13.9,+13.9,+13.9,+13.9,+13.9,+13.9,+13.9,+13.8,+13.8,+13.7,+13.7,+13.7,+13.7,+13.7,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.5,+13.5,+13.6,+13.7,+13.8,+13.9,+14.1,+14.1,+14.0,+13.9,+13.8,+13.6,+13.6,+13.5,+13.5,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.6,+13.7,+13.7,+13.7,+13.7,+13.7,+13.7,+13.7,+13.7,+13.7,+13.7,+13.8,+13.8,+13.8,+13.8,+13.8,+13.8,+13.9,+13.9,+14.0,+14.0,+14.1,+14.1,+14.1,+14.1,+14.2,+14.2,+14.3,+14.3,+14.5,+14.7,+15.0,+15.1,+15.2,+15.3,+15.5,+15.6,+15.7,+15.8,+16.0,+16.2,+16.4,+16.5,+16.7,+16.8,+17.0,+17.1,+17.3,+17.4,+17.5,+17.7,+17.9,+18.2,+18.5,+18.8,+19.3,+19.7,+20.1,+20.5,+20.9,+21.2};
   double ggH7_errm [] = {-15.4,-15.3,-15.2,-15.2,-15.0,-15.0,-15.0,-15.0,-15.0,-15.0,-15.0,-15.0,-15.0,-15.0,-15.0,-15.0,-15.0,-15.0,-15.0,-15.0,-15.0,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.8,-14.8,-14.8,-14.8,-14.8,-14.8,-14.8,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.9,-14.8,-14.8,-14.8,-14.8,-14.8,-14.8,-14.8,-14.8,-14.9,-14.9,-14.9,-14.9,-14.9,-14.8,-14.8,-14.8,-14.8,-14.8,-14.7,-14.7,-14.7,-14.7,-14.7,-14.6,-14.6,-14.6,-14.6,-14.6,-14.6,-14.6,-14.6,-14.6,-14.6,-14.6,-14.6,-14.5,-14.5,-14.5,-14.5,-14.5,-14.5,-14.5,-14.5,-14.5,-14.5,-14.4,-14.4,-14.4,-14.4,-14.4,-14.4,-14.4,-14.4,-14.4,-14.4,-14.4,-14.4,-14.4,-14.4,-14.4,-14.3,-14.3,-14.3,-14.3,-14.3,-14.3,-14.3,-14.3,-14.3,-14.4,-14.4,-14.4,-14.4,-14.3,-14.3,-14.3,-14.3,-14.4,-14.4,-14.4,-14.3,-14.3,-14.3,-14.3,-14.4,-14.3,-14.2,-14.1,-14.0,-13.9,-13.9,-14.0,-14.0,-14.1,-14.1,-14.2,-14.3,-14.4,-14.6,-14.7,-14.9,-15.0,-15.1,-15.2,-15.3,-15.4,-15.6,-15.7,-15.8,-16.0,-16.3,-16.7,-17.0,-17.4,-17.8,-18.2,-18.6,-18.9,-19.2};
   double ggH7_scap [] = {+8.2,+8.0,+7.8,+7.7,+7.5,+7.5,+7.5,+7.5,+7.5,+7.4,+7.4,+7.4,+7.4,+7.4,+7.4,+7.4,+7.4,+7.4,+7.3,+7.3,+7.3,+7.3,+7.3,+7.2,+7.2,+7.2,+7.2,+7.2,+7.2,+7.1,+7.1,+7.1,+7.1,+7.1,+7.1,+7.1,+7.1,+7.1,+7.1,+7.1,+7.0,+7.0,+7.0,+7.0,+7.0,+7.0,+7.0,+7.0,+7.0,+7.0,+6.9,+6.9,+6.9,+6.9,+6.9,+6.9,+6.9,+6.9,+6.9,+6.9,+6.8,+6.8,+6.8,+6.8,+6.8,+6.8,+6.8,+6.7,+6.7,+6.7,+6.7,+6.7,+6.6,+6.6,+6.6,+6.6,+6.6,+6.5,+6.5,+6.5,+6.5,+6.5,+6.4,+6.4,+6.4,+6.4,+6.4,+6.4,+6.4,+6.3,+6.3,+6.2,+6.2,+6.2,+6.2,+6.2,+6.1,+6.1,+6.1,+6.1,+6.1,+6.1,+6.1,+6.1,+6.0,+6.0,+6.0,+6.0,+6.0,+6.0,+6.1,+6.1,+6.2,+6.4,+6.5,+6.5,+6.4,+6.3,+6.1,+5.9,+5.9,+5.8,+5.8,+5.9,+5.9,+5.9,+5.9,+5.9,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.8,+5.9,+5.9,+5.9,+5.9,+5.9,+5.9,+5.9,+6.0,+6.0,+6.0,+6.0,+6.0,+6.1,+6.1,+6.2,+6.2,+6.3,+6.3,+6.3,+6.4,+6.4,+6.5,+6.5,+6.5,+6.5,+6.5,+6.6,+6.7,+6.8,+6.8,+6.8,+6.9,+7.0};
   double ggH7_scam [] = {-8.7,-8.6,-8.4,-8.3,-8.1,-8.1,-8.1,-8.1,-8.1,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-7.9,-7.9,-7.9,-7.9,-7.9,-7.9,-7.9,-7.9,-7.9,-7.9,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.7,-7.7,-7.7,-7.7,-7.7,-7.7,-7.7,-7.7,-7.7,-7.7,-7.6,-7.6,-7.6,-7.6,-7.6,-7.6,-7.6,-7.6,-7.6,-7.6,-7.5,-7.5,-7.5,-7.5,-7.5,-7.5,-7.5,-7.5,-7.5,-7.5,-7.5,-7.5,-7.5,-7.4,-7.4,-7.4,-7.4,-7.3,-7.3,-7.3,-7.3,-7.3,-7.2,-7.2,-7.2,-7.2,-7.2,-7.2,-7.2,-7.1,-7.1,-7.0,-7.0,-7.0,-7.0,-7.0,-6.9,-6.9,-6.9,-6.9,-6.9,-6.8,-6.8,-6.8,-6.8,-6.8,-6.8,-6.8,-6.7,-6.7,-6.7,-6.7,-6.6,-6.6,-6.6,-6.6,-6.6,-6.5,-6.5,-6.5,-6.5,-6.5,-6.4,-6.4,-6.4,-6.4,-6.4,-6.3,-6.3,-6.3,-6.3,-6.3,-6.3,-6.3,-6.3,-6.3,-6.3,-6.3,-6.2,-6.2,-6.2,-6.2,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.1,-6.0,-6.0,-6.0,-6.0,-6.0,-6.0,-6.0,-5.9,-5.9,-5.9,-5.9,-5.9,-5.8,-5.6,-5.5,-5.4,-5.3,-5.3,-5.3,-5.2,-5.2,-5.2,-5.2,-5.2,-5.2,-5.2,-5.2,-5.2,-5.2,-5.3,-5.3,-5.3,-5.4,-5.4,-5.4,-5.4,-5.4,-5.5,-5.5,-5.6,-5.6,-5.6,-5.7,-5.7,-5.7,-5.7};
   double ggH7_pdfp [] = {+7.8,+7.8,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.5,+7.6,+7.6,+7.6,+7.6,+7.6,+7.6,+7.5,+7.5,+7.5,+7.5,+7.6,+7.6,+7.6,+7.6,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.7,+7.8,+7.8,+7.8,+7.8,+7.8,+7.8,+7.8,+7.8,+7.8,+7.8,+7.9,+7.9,+7.9,+7.9,+7.9,+7.9,+7.9,+7.9,+7.9,+7.9,+8.0,+8.0,+8.0,+8.0,+8.0,+8.0,+8.1,+8.1,+8.2,+8.2,+8.3,+8.3,+8.3,+8.3,+8.4,+8.4,+8.4,+8.4,+8.6,+8.8,+9.1,+9.2,+9.3,+9.4,+9.5,+9.6,+9.7,+9.8,+9.9,+10.1,+10.2,+10.4,+10.5,+10.6,+10.7,+10.8,+10.9,+11.0,+11.1,+11.2,+11.4,+11.7,+11.9,+12.3,+12.6,+13.0,+13.3,+13.7,+14.0,+14.2};
   double ggH7_pdfm [] = {-6.7,-6.7,-6.8,-6.9,-6.9,-6.9,-6.9,-6.9,-6.9,-6.9,-7.0,-7.0,-7.0,-7.0,-7.0,-7.0,-7.0,-7.0,-7.0,-7.0,-7.0,-7.0,-7.0,-7.0,-7.0,-7.0,-7.0,-7.0,-7.0,-7.0,-7.1,-7.1,-7.1,-7.1,-7.1,-7.1,-7.1,-7.1,-7.1,-7.2,-7.2,-7.2,-7.2,-7.2,-7.2,-7.2,-7.2,-7.2,-7.2,-7.3,-7.3,-7.3,-7.3,-7.3,-7.3,-7.3,-7.3,-7.3,-7.3,-7.3,-7.3,-7.3,-7.3,-7.3,-7.3,-7.3,-7.3,-7.4,-7.4,-7.4,-7.4,-7.4,-7.5,-7.5,-7.5,-7.5,-7.5,-7.5,-7.5,-7.5,-7.5,-7.5,-7.6,-7.6,-7.6,-7.6,-7.7,-7.7,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.8,-7.9,-7.9,-7.9,-7.9,-7.9,-7.9,-7.9,-7.9,-7.9,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.1,-8.1,-8.1,-8.1,-8.1,-8.1,-8.1,-8.1,-8.1,-8.1,-8.1,-8.1,-8.1,-8.1,-8.1,-8.2,-8.2,-8.2,-8.2,-8.2,-8.3,-8.3,-8.3,-8.3,-8.3,-8.3,-8.3,-8.4,-8.4,-8.4,-8.4,-8.4,-8.4,-8.4,-8.4,-8.5,-8.6,-8.6,-8.6,-8.6,-8.6,-8.7,-8.7,-8.8,-8.9,-9.0,-9.0,-9.1,-9.2,-9.4,-9.5,-9.7,-9.8,-9.8,-9.9,-10.0,-10.1,-10.2,-10.3,-10.4,-10.6,-10.9,-11.1,-11.5,-11.8,-12.2,-12.5,-12.9,-13.2,-13.5};


   double qqH7_mass [] = {90.0,95.0,100.0,105.0,110.0,110.5,111.0,111.5,112.0,112.5,113.0,113.5,114.0,114.5,115.0,115.5,116.0,116.5,117.0,117.5,118.0,118.5,119.0,119.5,120.0,120.5,121.0,121.5,122.0,122.5,123.0,123.5,124.0,124.5,125.0,125.5,126.0,126.5,127.0,127.5,128.0,128.5,129.0,129.5,130.0,130.5,131.0,131.5,132.0,132.5,133.0,133.5,134.0,134.5,135.0,135.5,136.0,136.5,137.0,137.5,138.0,138.5,139.0,139.5,140.0,141.0,142.0,143.0,144.0,145.0,146.0,147.0,148.0,149.0,150.0,151.0,152.0,153.0,154.0,155.0,156.0,157.0,158.0,159.0,160.0,162.0,164.0,165.0,166.0,168.0,170.0,172.0,174.0,175.0,176.0,178.0,180.0,182.0,184.0,185.0,186.0,188.0,190.0,192.0,194.0,195.0,196.0,198.0,200.0,202.0,204.0,206.0,208.0,210.0,212.0,214.0,216.0,218.0,220.0,222.0,224.0,226.0,228.0,230.0,232.0,234.0,236.0,238.0,240.0,242.0,244.0,246.0,248.0,250.0,252.0,254.0,256.0,258.0,260.0,262.0,264.0,266.0,268.0,270.0,272.0,274.0,276.0,278.0,280.0,282.0,284.0,286.0,288.0,290.0,295.0,300.0,305.0,310.0,315.0,320.0,325.0,330.0,335.0,340.0,345.0,350.0,360.0,370.0,380.0,390.0,400.0,420.0,440.0,450.0,460.0,480.0,500.0,520.0,540.0,550.0,560.0,580.0,600.0,620.0,640.0,650.0,660.0,680.0,700.0,720.0,740.0,750.0,760.0,780.0,800.0,820.0,840.0,850.0,860.0,880.0,900.0,920.0,940.0,950.0,960.0,980.0,1000.0};
   double qqH7_xsec [] = {1.723,1.639,1.557,1.478,1.410,1.404,1.396,1.391,1.382,1.375,1.369,1.363,1.356,1.349,1.344,1.335,1.330,1.324,1.317,1.310,1.304,1.297,1.292,1.286,1.279,1.275,1.269,1.263,1.257,1.251,1.246,1.241,1.234,1.227,1.222,1.219,1.211,1.206,1.199,1.194,1.187,1.184,1.178,1.173,1.168,1.161,1.157,1.152,1.147,1.142,1.136,1.133,1.127,1.121,1.117,1.112,1.107,1.103,1.097,1.092,1.087,1.082,1.078,1.074,1.069,1.059,1.050,1.040,1.032,1.023,1.015,1.005,0.9980,0.9880,0.9800,0.9720,0.9640,0.9561,0.9487,0.9415,0.9339,0.9270,0.9199,0.9126,0.9043,0.8906,0.8755,0.8694,0.8613,0.8473,0.8338,0.8201,0.8063,0.7998,0.7934,0.7809,0.7684,0.7561,0.7433,0.7375,0.7314,0.7195,0.7080,0.6960,0.6845,0.6790,0.6735,0.6629,0.6524,0.6429,0.6343,0.6262,0.6184,0.6108,0.6033,0.5955,0.5879,0.5802,0.5724,0.5646,0.5570,0.5493,0.5416,0.5341,0.5266,0.5190,0.5114,0.5038,0.4959,0.4882,0.4807,0.4733,0.4661,0.4588,0.4519,0.4452,0.4385,0.4320,0.4256,0.4193,0.4131,0.4069,0.4010,0.3951,0.3894,0.3837,0.3783,0.3729,0.3676,0.3623,0.3572,0.3521,0.3471,0.3422,0.3387,0.3350,0.3237,0.3130,0.3028,0.2929,0.2834,0.2745,0.2661,0.2585,0.2519,0.2380,0.2142,0.2039,0.1949,0.1859,0.1772,0.1601,0.1442,0.1368,0.1296,0.1164,0.1046,0.0940,0.0846,0.0803,0.0763,0.0689,0.0624,0.0567,0.0515,0.0492,0.0469,0.0429,0.0392,0.0360,0.0330,0.0317,0.0304,0.0281,0.0260,0.0241,0.0224,0.0216,0.0208,0.0194,0.0181,0.0170,0.0159,0.0154,0.0149,0.0140,0.0132};
   double qqH7_errp [] = {+2.7,+2.5,+2.6,+2.5,+2.8,+2.8,+2.7,+2.7,+2.7,+2.6,+2.6,+2.6,+2.6,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.7,+2.7,+2.7,+2.7,+2.7,+2.7,+2.7,+2.7,+2.7,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.9,+2.9,+2.9,+3.0,+3.0,+3.1,+3.1,+3.1,+3.1,+3.0,+3.0,+3.0,+2.9,+2.9,+2.9,+2.9,+2.8,+2.8,+2.8,+2.9,+3.0,+3.1,+3.1,+3.1,+3.0,+3.0,+2.9,+2.9,+2.9,+2.9,+3.1,+3.1,+3.1,+3.1,+3.0,+3.0,+2.9,+2.9,+2.9,+3.1,+3.1,+3.1,+3.2,+3.2,+3.2,+3.2,+3.2,+3.2,+3.2,+3.1,+3.2,+3.3,+3.4,+3.4,+3.4,+3.3,+3.3,+3.4,+3.4,+3.4,+3.4,+3.4,+3.4,+3.4,+3.5,+3.5,+3.5,+3.5,+3.5,+3.5,+3.5,+3.5,+3.5,+3.5,+3.7,+3.7,+3.7,+3.7,+3.7,+3.8,+3.8,+3.8,+3.8,+3.8,+3.9,+3.9,+3.9,+3.9,+4.0,+4.1,+4.2,+4.2,+4.2,+4.1,+4.2,+4.1,+4.1,+4.1,+4.2,+4.3,+4.4,+4.4,+4.4,+4.4,+4.5,+4.5,+4.5,+4.5,+4.6,+4.7,+4.7,+4.8,+4.9,+4.9,+5.0,+5.1,+5.1,+5.1,+5.2,+5.3,+5.5,+5.6,+5.7,+5.9,+6.2,+6.5,+6.6,+6.7,+7.1,+7.3,+7.5,+7.8,+7.9,+8.0,+8.3,+8.6,+8.8,+9.2,+9.3,+9.4,+9.7,+9.9,+10.2,+10.6,+10.7,+10.8,+11.1,+11.3,+11.5,+11.9,+12.0,+12.1,+12.4,+12.6,+12.9,+13.3,+13.5,+13.6,+13.9,+14.2};
   double qqH7_errm [] = {-2.3,-2.5,-2.4,-2.4,-2.3,-2.3,-2.3,-2.3,-2.3,-2.3,-2.3,-2.3,-2.3,-2.3,-2.3,-2.3,-2.3,-2.4,-2.4,-2.4,-2.4,-2.4,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.4,-2.4,-2.4,-2.4,-2.4,-2.4,-2.4,-2.4,-2.4,-2.4,-2.3,-2.3,-2.3,-2.3,-2.3,-2.3,-2.3,-2.3,-2.3,-2.3,-2.3,-2.2,-2.2,-2.2,-2.2,-2.2,-2.2,-2.2,-2.2,-2.2,-2.2,-2.3,-2.3,-2.3,-2.3,-2.3,-2.3,-2.2,-2.2,-2.1,-2.1,-2.1,-2.1,-2.2,-2.2,-2.2,-2.2,-2.2,-2.1,-2.1,-2.1,-2.1,-2.2,-2.2,-2.3,-2.3,-2.3,-2.2,-2.2,-2.2,-2.3,-2.3,-2.3,-2.2,-2.2,-2.2,-2.3,-2.4,-2.3,-2.1,-2.1,-2.1,-2.2,-2.2,-2.3,-2.4,-2.4,-2.4,-2.3,-2.2,-2.2,-2.2,-2.3,-2.3,-2.3,-2.3,-2.3,-2.4,-2.4,-2.4,-2.4,-2.4,-2.4,-2.4,-2.4,-2.4,-2.4,-2.5,-2.5,-2.5,-2.5,-2.5,-2.6,-2.6,-2.6,-2.6,-2.5,-2.5,-2.4,-2.4,-2.4,-2.5,-2.5,-2.6,-2.6,-2.6,-2.6,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.8,-2.8,-2.8,-2.7,-2.6,-2.6,-2.7,-2.7,-2.8,-2.8,-2.9,-2.9,-3.0,-3.0,-3.0,-3.0,-3.1,-3.1,-3.1,-3.1,-3.2,-3.3,-3.4,-3.4,-3.5,-3.5,-3.6,-3.6,-3.7,-3.8,-3.9,-3.9,-3.9,-3.9,-4.0,-4.1,-4.2,-4.2,-4.2,-4.3,-4.4,-4.5,-4.5,-4.5,-4.5,-4.6,-4.7,-4.7,-4.7,-4.7,-4.7,-4.8,-4.9};
   double qqH7_scap [] = {+0.6,+0.4,+0.4,+0.3,+0.5,+0.5,+0.4,+0.4,+0.4,+0.3,+0.3,+0.3,+0.3,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.4,+0.4,+0.4,+0.4,+0.4,+0.5,+0.5,+0.5,+0.5,+0.4,+0.4,+0.4,+0.3,+0.3,+0.3,+0.3,+0.2,+0.2,+0.2,+0.3,+0.3,+0.4,+0.4,+0.4,+0.3,+0.3,+0.2,+0.2,+0.2,+0.2,+0.3,+0.3,+0.3,+0.3,+0.2,+0.2,+0.1,+0.1,+0.1,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.1,+0.0,+0.1,+0.2,+0.3,+0.3,+0.2,+0.1,+0.1,+0.2,+0.2,+0.2,+0.1,+0.1,+0.1,+0.1,+0.1,+0.1,+0.1,+0.1,+0.1,+0.0,+0.0,+0.0,+0.0,+0.0,+0.1,+0.1,+0.1,+0.1,+0.1,+0.1,+0.1,+0.1,+0.1,+0.1,+0.1,+0.1,+0.1,+0.1,+0.2,+0.2,+0.3,+0.3,+0.3,+0.2,+0.2,+0.1,+0.1,+0.1,+0.1,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.4,+0.4,+0.4,+0.4,+0.5,+0.6,+0.6,+0.6,+0.7,+0.7,+0.7,+0.8,+0.8,+0.8,+0.9,+1.0,+1.0,+1.1,+1.1,+1.1,+1.2,+1.2,+1.3,+1.4,+1.4,+1.4,+1.5,+1.5,+1.5,+1.6,+1.6,+1.6,+1.7,+1.7,+1.8,+1.9,+2.0,+2.0,+2.1,+2.2};
   double qqH7_scam [] = {-0.2,-0.4,-0.3,-0.3,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.3,-0.3,-0.3,-0.3,-0.3,-0.4,-0.4,-0.4,-0.4,-0.4,-0.4,-0.4,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.1,-0.1,0.0,0.0,0.0,0.0,-0.1,-0.1,-0.1,-0.1,-0.1,0.0,0.0,0.0,0.0,-0.1,-0.1,-0.2,-0.2,-0.2,-0.1,-0.1,-0.1,-0.2,-0.2,-0.2,-0.1,-0.1,-0.1,-0.2,-0.3,-0.2,-0.1,-0.1,-0.1,-0.2,-0.2,-0.3,-0.4,-0.4,-0.4,-0.3,-0.2,-0.2,-0.2,-0.3,-0.3,-0.3,-0.3,-0.3,-0.4,-0.4,-0.4,-0.4,-0.4,-0.4,-0.4,-0.4,-0.4,-0.4,-0.5,-0.5,-0.5,-0.5,-0.5,-0.6,-0.6,-0.6,-0.6,-0.5,-0.5,-0.4,-0.4,-0.4,-0.5,-0.5,-0.6,-0.6,-0.6,-0.6,-0.7,-0.7,-0.7,-0.7,-0.7,-0.7,-0.7,-0.7,-0.8,-0.8,-0.8,-0.8,-0.7,-0.7,-0.8,-0.8,-0.9,-0.9,-1.0,-1.0,-1.1,-1.1,-1.1,-1.1,-1.2,-1.2,-1.3,-1.3,-1.4,-1.5,-1.6,-1.6,-1.7,-1.7,-1.8,-1.9,-2.0,-2.1,-2.2,-2.2,-2.2,-2.3,-2.4,-2.5,-2.6,-2.6,-2.6,-2.7,-2.8,-2.9,-3.0,-3.0,-3.0,-3.1,-3.2,-3.2,-3.3,-3.3,-3.3,-3.4,-3.5};
   double qqH7_pdfp [] = {+2.1,+2.1,+2.2,+2.2,+2.3,+2.3,+2.3,+2.3,+2.3,+2.3,+2.3,+2.3,+2.3,+2.3,+2.3,+2.3,+2.3,+2.3,+2.3,+2.3,+2.4,+2.4,+2.4,+2.4,+2.4,+2.4,+2.4,+2.4,+2.4,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.7,+2.7,+2.7,+2.7,+2.7,+2.7,+2.7,+2.7,+2.7,+2.7,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.9,+2.9,+2.9,+3.0,+3.0,+3.0,+3.0,+3.0,+3.0,+3.1,+3.1,+3.1,+3.1,+3.1,+3.1,+3.2,+3.2,+3.2,+3.2,+3.2,+3.2,+3.3,+3.3,+3.3,+3.3,+3.4,+3.4,+3.4,+3.4,+3.4,+3.5,+3.5,+3.5,+3.5,+3.5,+3.6,+3.6,+3.6,+3.6,+3.6,+3.7,+3.7,+3.7,+3.7,+3.7,+3.8,+3.8,+3.8,+3.8,+3.8,+3.9,+3.9,+3.9,+3.9,+3.9,+4.0,+4.0,+4.0,+4.0,+4.1,+4.1,+4.2,+4.2,+4.2,+4.2,+4.3,+4.3,+4.3,+4.3,+4.4,+4.5,+4.5,+4.5,+4.6,+4.6,+4.7,+4.8,+4.8,+4.8,+4.9,+5.0,+5.1,+5.2,+5.3,+5.5,+5.7,+5.9,+6.0,+6.1,+6.4,+6.6,+6.8,+7.0,+7.1,+7.2,+7.4,+7.6,+7.8,+8.1,+8.2,+8.3,+8.5,+8.7,+8.9,+9.2,+9.3,+9.4,+9.6,+9.8,+10.0,+10.3,+10.4,+10.5,+10.7,+10.9,+11.1,+11.4,+11.5,+11.6,+11.8,+12.0};
   double qqH7_pdfm [] = {-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-1.9,-1.9,-1.9,-1.9,-1.9,-1.9,-1.9,-1.9,-1.9,-1.9,-1.9,-1.9,-1.9,-1.9,-1.9,-1.8,-1.8,-1.8,-1.8,-1.8,-1.8,-1.8,-1.8,-1.8,-1.7,-1.7,-1.7,-1.7,-1.7,-1.7,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.5,-1.5,-1.5,-1.5,-1.5,-1.5,-1.4,-1.4,-1.4,-1.4,-1.4};

   double ggH8_mass [] = {80.0, 81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0, 90.0, 91.0, 92.0, 93.0, 94.0, 95.0, 96.0, 97.0, 98.0, 99.0, 100.0, 101.0, 102.0, 103.0, 104.0, 105.0, 106.0, 107.0, 108.0, 109.0, 110.0, 110.5, 111.0, 111.5, 112.0, 112.5, 113.0, 113.5, 114.0, 114.5, 115.0, 115.5, 116.0, 116.5, 117.0, 117.5, 118.0, 118.5, 119.0, 119.5, 120.0, 120.5, 121.0, 121.5, 122.0, 122.5, 123.0, 123.5, 124.0, 124.5, 125.0, 125.5, 126.0, 126.5, 127.0, 127.5, 128.0, 128.5, 129.0, 129.5, 130.0, 130.5, 131.0, 131.5, 132.0, 132.5, 133.0, 133.5, 134.0, 134.5, 135.0, 135.5, 136.0, 136.5, 137.0, 137.5, 138.0, 138.5, 139.0, 139.5, 140.0, 141.0, 142.0, 143.0, 144.0, 145.0, 146.0, 147.0, 148.0, 149.0, 150.0, 151.0, 152.0, 153.0, 154.0, 155.0, 156.0, 157.0, 158.0, 159.0, 160.0, 162.0, 164.0, 165.0, 166.0, 168.0, 170.0, 172.0, 174.0, 175.0, 176.0, 178.0, 180.0, 182.0, 184.0, 185.0, 186.0, 188.0, 190.0, 192.0, 194.0, 195.0, 196.0, 198.0, 200.0, 202.0, 204.0, 206.0, 208.0, 210.0, 212.0, 214.0, 216.0, 218.0, 220.0, 222.0, 224.0, 226.0, 228.0, 230.0, 232.0, 234.0, 236.0, 238.0, 240.0, 242.0, 244.0, 246.0, 248.0, 250.0, 252.0, 254.0, 256.0, 258.0, 260.0, 262.0, 264.0, 266.0, 268.0, 270.0, 272.0, 274.0, 276.0, 278.0, 280.0, 282.0, 284.0, 286.0, 288.0, 290.0, 295.0, 300.0, 305.0, 310.0, 315.0, 320.0, 325.0, 330.0, 335.0, 340.0, 345.0, 350.0, 360.0, 370.0, 380.0, 390.0, 400.0, 420.0, 440.0, 450.0, 460.0, 480.0, 500.0, 520.0, 540.0, 550.0, 560.0, 580.0, 600.0, 620.0, 640.0, 650.0, 660.0, 680.0, 700.0, 720.0, 740.0, 750.0, 760.0, 780.0, 800.0, 820.0, 840.0, 850.0, 860.0, 880.0, 900.0, 920.0, 940.0, 950.0, 960.0, 980.0, 1000.0};
   double ggH8_xsec [] = {46.12, 45.04, 43.99, 42.99, 42.01, 41.07, 40.17, 39.29, 38.44, 37.62, 36.80, 36.05, 35.30, 34.58, 33.87, 33.19, 32.53, 31.89, 31.27, 30.66, 30.12, 29.55, 28.99, 28.44, 27.92, 27.39, 26.89, 26.42, 25.95, 25.49, 25.04, 24.82, 24.60, 24.39, 24.18, 24.05, 23.76, 23.56, 23.36, 23.16, 22.96, 22.84, 22.58, 22.39, 22.20, 22.09, 21.90, 21.72, 21.55, 21.37, 21.20, 20.96, 20.86, 20.69, 20.53, 20.37, 20.21, 20.04, 19.89, 19.73, 19.57, 19.42, 19.27, 19.07, 18.97, 18.78, 18.67, 18.53, 18.39, 18.21, 18.11, 17.98, 17.84, 17.68, 17.55, 17.42, 17.29, 17.19, 17.07, 16.92, 16.82, 16.70, 16.60, 16.46, 16.36, 16.24, 16.13, 16.01, 15.89, 15.78, 15.67, 15.43, 15.20, 14.98, 14.79, 14.59, 14.40, 14.31, 14.11, 13.92, 13.65, 13.33, 13.12, 12.97, 12.86, 12.84, 12.92, 12.84, 12.52, 12.17, 11.98, 11.96, 11.27, 10.96, 10.69, 10.40, 10.15, 9.895, 9.637, 9.531, 9.429, 9.203, 8.923, 8.634, 8.410, 8.315, 8.219, 7.994, 7.984, 7.724, 7.443, 7.431, 7.430, 7.237, 7.127, 7.002, 6.885, 6.767, 6.644, 6.534, 6.421, 6.327, 6.249, 6.138, 6.038, 5.938, 5.851, 5.776, 5.680, 5.593, 5.513, 5.438, 5.352, 5.272, 5.183, 5.101, 5.020, 4.945, 4.874, 4.802, 4.734, 4.667, 4.602, 4.539, 4.479, 4.421, 4.366, 4.308, 4.253, 4.198, 4.149, 4.100, 4.053, 4.008, 3.964, 3.921, 3.880, 3.841, 3.803, 3.767, 3.683, 3.606, 3.539, 3.482, 3.434, 3.392, 3.364, 3.349, 3.349, 3.367, 3.405, 3.406, 3.390, 3.336, 3.235, 3.093, 2.924, 2.552, 2.180, 2.003, 1.837, 1.538, 1.283, 1.069, 0.8912, 0.8141, 0.7442, 0.6229, 0.5230, 0.4403, 0.3718, 0.3423, 0.3152, 0.2680, 0.2288, 0.1962, 0.1687, 0.1566, 0.1455, 0.1260, 0.1095, 0.09547, 0.08346, 0.07811, 0.07321, 0.06443, 0.05684, 0.05030, 0.04463, 0.04206, 0.03969, 0.03539, 0.03163};
   double ggH8_errp [] = {+16.7, +16.7, +16.6, +16.6, +16.5, +16.5, +16.4, +16.4, +16.3, +16.2, +16.1, +16.1, +16.0, +16.0, +15.9, +15.9, +15.9, +15.8, +15.8, +15.8, +15.7, +15.7, +15.7, +15.6, +15.6, +15.5, +15.4, +15.4, +15.4, +15.3, +15.3, +15.3, +15.2, +15.1, +15.1, +15.1, +15.1, +15.1, +15.1, +15.1, +15.0, +15.0, +14.9, +14.9, +14.9, +14.9, +14.9, +14.8, +14.8, +14.8, +14.8, +14.8, +14.8, +14.8, +14.8, +14.7, +14.7, +14.7, +14.7, +14.7, +14.7, +14.7, +14.7, +14.7, +14.6, +14.6, +14.6, +14.6, +14.6, +14.6, +14.6, +14.6, +14.6, +14.5, +14.5, +14.5, +14.4, +14.4, +14.4, +14.4, +14.4, +14.4, +14.3, +14.3, +14.3, +14.3, +14.3, +14.3, +14.3, +14.3, +14.3, +14.1, +14.1, +14.1, +14.1, +14.1, +14.0, +14.0, +14.0, +14.1, +14.1, +14.1, +14.0, +14.1, +14.1, +14.1, +14.1, +14.1, +14.0, +14.0, +14.0, +14.0, +13.9, +13.9, +13.9, +13.9, +13.9, +13.8, +13.7, +13.7, +13.7, +13.7, +13.6, +13.6, +13.6, +13.6, +13.5, +13.5, +13.5, +13.5, +13.5, +13.5, +13.4, +13.4, +13.4, +13.4, +13.4, +13.4, +13.4, +13.4, +13.4, +13.4, +13.3, +13.3, +13.2, +13.2, +13.3, +13.3, +13.3, +13.3, +13.3, +13.2, +13.2, +13.2, +13.2, +13.2, +13.2, +13.2, +13.2, +13.2, +13.2, +13.3, +13.3, +13.3, +13.4, +13.4, +13.4, +13.4, +13.4, +13.4, +13.4, +13.4, +13.4, +13.3, +13.3, +13.3, +13.3, +13.3, +13.3, +13.3, +13.3, +13.4, +13.4, +13.4, +13.4, +13.4, +13.4, +13.5, +13.6, +13.6, +13.6, +13.7, +13.8, +13.9, +13.9, +14.0, +14.0, +14.1, +14.3, +14.4, +14.5, +14.7, +14.9, +15.0, +15.2, +15.2, +15.3, +15.4, +15.4, +15.5, +15.6, +15.6, +15.7, +15.9, +16.1, +16.2, +16.5, +16.5, +16.6, +16.6, +16.7, +16.9, +17.0, +17.1, +17.1, +17.3, +17.4, +17.5, +17.8, +17.9, +18.0, +18.3, +18.6};
   double ggH8_errm [] = {-15.9, -15.8, -15.7, -15.6, -15.6, -15.6, -15.5, -15.5, -15.4, -15.4, -15.4, -15.3, -15.2, -15.2, -15.1, -15.1, -15.1, -15.0, -15.0, -15.0, -14.9, -14.9, -14.9, -14.9, -14.8, -14.8, -14.9, -14.9, -14.8, -14.9, -14.9, -14.9, -14.9, -14.8, -14.8, -14.8, -14.8, -14.9, -14.9, -14.9, -14.9, -14.8, -14.8, -14.8, -14.8, -14.8, -14.8, -14.8, -14.8, -14.8, -14.8, -14.8, -14.8, -14.8, -14.8, -14.8, -14.8, -14.8, -14.8, -14.8, -14.7, -14.7, -14.7, -14.7, -14.7, -14.7, -14.7, -14.7, -14.7, -14.7, -14.6, -14.6, -14.7, -14.7, -14.7, -14.7, -14.7, -14.7, -14.7, -14.7, -14.7, -14.6, -14.6, -14.6, -14.6, -14.6, -14.5, -14.5, -14.5, -14.5, -14.5, -14.5, -14.4, -14.4, -14.4, -14.4, -14.4, -14.4, -14.4, -14.4, -14.4, -14.4, -14.5, -14.5, -14.4, -14.4, -14.4, -14.4, -14.4, -14.4, -14.4, -14.4, -14.4, -14.5, -14.5, -14.5, -14.5, -14.5, -14.5, -14.5, -14.5, -14.5, -14.5, -14.5, -14.5, -14.5, -14.4, -14.4, -14.4, -14.5, -14.5, -14.5, -14.4, -14.5, -14.5, -14.5, -14.5, -14.4, -14.5, -14.5, -14.4, -14.4, -14.3, -14.2, -14.2, -14.2, -14.2, -14.2, -14.2, -14.2, -14.2, -14.2, -14.2, -14.2, -14.1, -14.1, -14.1, -14.1, -14.1, -14.1, -14.1, -13.9, -13.8, -13.7, -13.7, -13.8, -13.9, -14.0, -14.0, -14.1, -14.1, -14.1, -14.2, -14.2, -14.2, -14.2, -14.1, -14.1, -14.1, -14.1, -14.1, -14.0, -13.9, -14.0, -14.0, -14.0, -14.1, -14.1, -14.0, -14.0, -14.1, -14.1, -14.1, -14.0, -13.8, -13.7, -13.6, -13.6, -13.7, -13.6, -13.6, -13.7, -13.6, -13.7, -13.7, -13.8, -13.8, -13.8, -13.8, -13.9, -14.0, -14.0, -14.2, -14.3, -14.4, -14.6, -14.7, -14.8, -14.9, -15.0, -15.0, -15.0, -15.1, -15.2, -15.3, -15.4, -15.5, -15.7, -16.0, -16.1, -16.2, -16.4, -16.7};
   double ggH8_scap [] = {+8.8, +8.8, +8.7, +8.7, +8.6, +8.6, +8.5, +8.5, +8.4, +8.4, +8.3, +8.3, +8.2, +8.2, +8.1, +8.1, +8.1, +8.0, +8.0, +8.0, +7.9, +7.9, +7.9, +7.8, +7.8, +7.8, +7.7, +7.7, +7.7, +7.6, +7.6, +7.6, +7.6, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.3, +7.3, +7.3, +7.3, +7.3, +7.3, +7.3, +7.3, +7.2, +7.2, +7.2, +7.2, +7.2, +7.2, +7.2, +7.2, +7.2, +7.1, +7.1, +7.1, +7.1, +7.1, +7.1, +7.1, +7.1, +7.1, +7.0, +7.0, +7.0, +7.0, +7.0, +7.0, +7.0, +7.0, +7.0, +6.9, +6.9, +6.9, +6.9, +6.9, +6.9, +6.9, +6.9, +6.9, +6.8, +6.8, +6.8, +6.8, +6.8, +6.7, +6.7, +6.7, +6.7, +6.7, +6.7, +6.6, +6.6, +6.6, +6.6, +6.6, +6.6, +6.5, +6.5, +6.5, +6.5, +6.4, +6.4, +6.4, +6.4, +6.4, +6.3, +6.3, +6.3, +6.3, +6.3, +6.2, +6.2, +6.2, +6.2, +6.1, +6.1, +6.1, +6.1, +6.1, +6.1, +6.0, +6.0, +6.0, +6.0, +6.0, +6.0, +6.0, +6.0, +6.0, +6.0, +5.9, +5.9, +5.9, +5.9, +5.9, +5.9, +5.9, +5.9, +5.9, +5.9, +5.9, +5.9, +5.9, +5.9, +5.9, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.7, +5.7, +5.7, +5.7, +5.7, +5.7, +5.7, +5.7, +5.7, +5.7, +5.7, +5.7, +5.7, +5.7, +5.7, +5.7, +5.7, +5.7, +5.7, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.8, +5.9, +5.9, +5.9, +5.9, +5.9, +5.9, +5.9, +6.0, +6.0, +6.0, +6.1, +6.1, +6.1, +6.1, +6.1, +6.2, +6.2, +6.2, +6.2, +6.3, +6.3, +6.3, +6.4, +6.4, +6.4, +6.5, +6.5};
   double ggH8_scam [] = {-9.2, -9.1, -9.1, -9.0, -9.0, -9.0, -8.9, -8.9, -8.8, -8.8, -8.8, -8.7, -8.7, -8.7, -8.6, -8.6, -8.6, -8.5, -8.5, -8.5, -8.4, -8.4, -8.4, -8.4, -8.3, -8.3, -8.3, -8.3, -8.2, -8.2, -8.2, -8.2, -8.2, -8.1, -8.1, -8.1, -8.1, -8.1, -8.1, -8.1, -8.1, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -7.9, -7.9, -7.9, -7.9, -7.9, -7.9, -7.9, -7.9, -7.9, -7.9, -7.8, -7.8, -7.8, -7.8, -7.8, -7.8, -7.8, -7.8, -7.8, -7.8, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.6, -7.6, -7.6, -7.6, -7.6, -7.6, -7.6, -7.6, -7.6, -7.6, -7.6, -7.5, -7.5, -7.5, -7.5, -7.5, -7.5, -7.4, -7.4, -7.4, -7.4, -7.4, -7.4, -7.3, -7.3, -7.3, -7.3, -7.3, -7.3, -7.3, -7.2, -7.2, -7.2, -7.2, -7.2, -7.1, -7.1, -7.1, -7.1, -7.1, -7.0, -7.0, -7.0, -7.0, -7.0, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.8, -6.8, -6.8, -6.8, -6.8, -6.7, -6.7, -6.7, -6.7, -6.7, -6.6, -6.6, -6.6, -6.6, -6.6, -6.6, -6.5, -6.5, -6.5, -6.5, -6.5, -6.5, -6.4, -6.4, -6.4, -6.4, -6.4, -6.4, -6.4, -6.3, -6.3, -6.3, -6.3, -6.3, -6.3, -6.3, -6.2, -6.2, -6.2, -6.2, -6.2, -6.2, -6.2, -6.2, -6.1, -6.1, -6.1, -6.1, -6.1, -6.1, -6.0, -6.0, -6.0, -6.0, -6.0, -6.0, -5.9, -5.9, -5.9, -5.9, -5.9, -5.8, -5.6, -5.5, -5.4, -5.3, -5.3, -5.2, -5.2, -5.2, -5.1, -5.1, -5.1, -5.1, -5.1, -5.1, -5.0, -5.0, -5.0, -5.0, -5.1, -5.1, -5.1, -5.1, -5.1, -5.1, -5.2, -5.2, -5.2, -5.2, -5.2, -5.3, -5.3, -5.3, -5.3, -5.3, -5.4, -5.4, -5.4, -5.4, -5.4};
   double ggH8_pdfp [] = {+7.9, +7.9, +7.9, +7.9, +7.9, +7.9, +7.9, +7.9, +7.9, +7.8, +7.8, +7.8, +7.8, +7.8, +7.8, +7.8, +7.8, +7.8, +7.8, +7.8, +7.8, +7.8, +7.8, +7.8, +7.8, +7.7, +7.7, +7.7, +7.7, +7.7, +7.7, +7.7, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.3, +7.3, +7.3, +7.3, +7.3, +7.3, +7.3, +7.3, +7.4, +7.4, +7.4, +7.4, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.5, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.4, +7.3, +7.3, +7.4, +7.4, +7.4, +7.4, +7.4, +7.3, +7.3, +7.3, +7.3, +7.3, +7.3, +7.4, +7.4, +7.4, +7.4, +7.5, +7.5, +7.5, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.6, +7.7, +7.7, +7.7, +7.7, +7.7, +7.7, +7.8, +7.9, +7.9, +7.9, +8.0, +8.0, +8.1, +8.1, +8.2, +8.2, +8.3, +8.5, +8.6, +8.7, +8.9, +9.1, +9.2, +9.4, +9.4, +9.4, +9.5, +9.5, +9.6, +9.7, +9.7, +9.8, +9.9, +10.1, +10.2, +10.4, +10.4, +10.5, +10.5, +10.6, +10.7, +10.8, +10.9, +10.9, +11.0, +11.1, +11.2, +11.4, +11.5, +11.6, +11.8, +12.1};
   double ggH8_pdfm [] = {-6.7, -6.7, -6.6, -6.6, -6.6, -6.6, -6.6, -6.6, -6.6, -6.6, -6.6, -6.6, -6.5, -6.5, -6.5, -6.5, -6.5, -6.5, -6.5, -6.5, -6.5, -6.5, -6.5, -6.5, -6.5, -6.5, -6.6, -6.6, -6.6, -6.7, -6.7, -6.7, -6.7, -6.7, -6.7, -6.7, -6.7, -6.8, -6.8, -6.8, -6.8, -6.8, -6.8, -6.8, -6.8, -6.8, -6.8, -6.8, -6.8, -6.8, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -7.0, -7.0, -7.0, -7.0, -7.0, -7.0, -7.0, -7.0, -7.0, -7.0, -7.0, -7.0, -7.0, -7.0, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -6.9, -7.0, -7.0, -7.0, -7.0, -7.1, -7.1, -7.1, -7.1, -7.1, -7.1, -7.1, -7.1, -7.1, -7.2, -7.2, -7.3, -7.3, -7.3, -7.4, -7.4, -7.4, -7.4, -7.4, -7.5, -7.5, -7.5, -7.5, -7.5, -7.5, -7.5, -7.5, -7.6, -7.6, -7.6, -7.6, -7.7, -7.7, -7.7, -7.7, -7.7, -7.8, -7.8, -7.7, -7.7, -7.7, -7.6, -7.6, -7.6, -7.6, -7.6, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.7, -7.6, -7.5, -7.4, -7.4, -7.5, -7.6, -7.7, -7.8, -7.9, -7.9, -7.9, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -7.9, -7.9, -8.0, -8.0, -8.0, -8.1, -8.1, -8.1, -8.1, -8.2, -8.2, -8.2, -8.2, -8.2, -8.2, -8.2, -8.3, -8.4, -8.4, -8.4, -8.5, -8.5, -8.6, -8.6, -8.7, -8.7, -8.7, -8.8, -8.9, -9.0, -9.0, -9.1, -9.2, -9.3, -9.5, -9.6, -9.7, -9.7, -9.8, -9.8, -9.8, -9.9, -9.9, -10.0, -10.1, -10.2, -10.4, -10.6, -10.7, -10.8, -11.0, -11.3};

   double qqH8_mass [] = {80.0,81.0,82.0,83.0,84.0,85.0,86.0,87.0,88.0,89.0,90.0,91.0,92.0,93.0,94.0,95.0,96.0,97.0,98.0,99.0,100.0,101.0,102.0,103.0,104.0,105.0,106.0,107.0,108.0,109.0,110.0,110.5,111.0,111.5,112.0,112.5,113.0,113.5,114.0,114.5,115.0,115.5,116.0,116.5,117.0,117.5,118.0,118.5,119.0,119.5,120.0,120.5,121.0,121.5,122.0,122.5,123.0,123.5,124.0,124.5,125.0,125.5,126.0,126.5,127.0,127.5,128.0,128.5,129.0,129.5,130.0,130.5,131.0,131.5,132.0,132.5,133.0,133.5,134.0,134.5,135.0,135.5,136.0,136.5,137.0,137.5,138.0,138.5,139.0,139.5,140.0,141.0,142.0,143.0,144.0,145.0,146.0,147.0,148.0,149.0,150.0,151.0,152.0,153.0,154.0,155.0,156.0,157.0,158.0,159.0,160.0,162.0,164.0,165.0,166.0,168.0,170.0,172.0,174.0,175.0,176.0,178.0,180.0,182.0,184.0,185.0,186.0,188.0,190.0,192.0,194.0,195.0,196.0,198.0,200.0,202.0,204.0,206.0,208.0,210.0,212.0,214.0,216.0,218.0,220.0,222.0,224.0,226.0,228.0,230.0,232.0,234.0,236.0,238.0,240.0,242.0,244.0,246.0,248.0,250.0,252.0,254.0,256.0,258.0,260.0,262.0,264.0,266.0,268.0,270.0,272.0,274.0,276.0,278.0,280.0,282.0,284.0,286.0,288.0,290.0,295.0,300.0,305.0,310.0,315.0,320.0,325.0,330.0,335.0,340.0,345.0,350.0,360.0,370.0,380.0,390.0,400.0,420.0,440.0,450.0,460.0,480.0,500.0,520.0,540.0,550.0,560.0,580.0,600.0,620.0,640.0,650.0,660.0,680.0,700.0,720.0,740.0,750.0,760.0,780.0,800.0,820.0,840.0,850.0,860.0,880.0,900.0,920.0,940.0,950.0,960.0,980.0,1000.0};
   double qqH8_xsec [] = {2.424,2.399,2.364,2.346,2.326,2.300,2.283,2.258,2.240,2.209,2.191,2.170,2.153,2.129,2.108,2.084,2.068,2.046,2.027,2.004,1.988,1.967,1.945,1.933,1.914,1.897,1.877,1.862,1.841,1.826,1.809,1.799,1.791,1.784,1.780,1.771,1.764,1.753,1.743,1.735,1.729,1.719,1.714,1.704,1.699,1.688,1.683,1.675,1.666,1.659,1.649,1.643,1.636,1.631,1.623,1.615,1.608,1.598,1.595,1.587,1.578,1.573,1.568,1.558,1.552,1.543,1.540,1.531,1.525,1.513,1.511,1.504,1.497,1.492,1.485,1.479,1.473,1.466,1.462,1.455,1.448,1.444,1.436,1.429,1.423,1.417,1.412,1.407,1.400,1.396,1.389,1.377,1.365,1.354,1.344,1.333,1.321,1.311,1.302,1.291,1.280,1.268,1.259,1.250,1.240,1.231,1.222,1.212,1.204,1.194,1.185,1.171,1.152,1.141,1.132,1.114,1.098,1.080,1.062,1.055,1.047,1.031,1.015,0.9980,0.9830,0.9760,0.9690,0.9536,0.9387,0.9238,0.9090,0.9018,0.8953,0.8819,0.8685,0.8568,0.8456,0.8356,0.8259,0.8163,0.8067,0.7970,0.7873,0.7776,0.7677,0.7579,0.7481,0.7381,0.7287,0.7190,0.7095,0.6999,0.6903,0.6806,0.6703,0.6604,0.6506,0.6410,0.6319,0.6225,0.6136,0.6050,0.5964,0.5879,0.5797,0.5714,0.5636,0.5554,0.5477,0.5401,0.5328,0.5255,0.5184,0.5115,0.5045,0.4978,0.4911,0.4845,0.4780,0.4716,0.4562,0.4408,0.4266,0.4131,0.3999,0.3875,0.3753,0.3637,0.3526,0.3422,0.3303,0.3200,0.3028,0.2896,0.2776,0.2660,0.2543,0.2317,0.2103,0.2002,0.1905,0.1724,0.1561,0.1414,0.1283,0.1223,0.1166,0.1062,0.09688,0.08861,0.08121,0.07784,0.07459,0.06865,0.06330,0.05853,0.05420,0.05235,0.05032,0.04682,0.04365,0.04078,0.03815,0.03706,0.03579,0.03363,0.03164,0.02986,0.02820,0.02745,0.02669,0.02524,0.02399};
   double qqH8_errp [] = {+2.9,+3.0,+2.9,+3.0,+2.9,+2.9,+2.9,+2.9,+2.8,+2.9,+2.9,+2.9,+2.8,+2.9,+2.9,+2.8,+2.9,+2.9,+2.9,+2.8,+2.8,+2.9,+2.7,+2.8,+2.8,+2.8,+2.8,+2.7,+2.7,+2.8,+2.7,+2.7,+2.7,+2.7,+2.7,+2.8,+2.8,+2.8,+2.7,+2.7,+2.7,+2.8,+2.8,+2.8,+2.7,+2.8,+2.9,+2.8,+2.9,+2.8,+2.8,+2.8,+2.8,+2.8,+2.9,+2.8,+2.8,+2.8,+2.9,+2.8,+2.8,+2.8,+2.9,+2.8,+2.9,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.7,+2.7,+2.8,+2.8,+2.7,+2.7,+2.7,+2.7,+2.7,+2.7,+2.7,+2.8,+2.7,+2.7,+2.8,+2.7,+2.8,+2.7,+2.7,+2.8,+2.7,+2.7,+2.7,+2.7,+2.7,+2.8,+2.7,+2.7,+2.7,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.9,+2.8,+2.7,+2.7,+2.7,+2.8,+2.8,+2.7,+2.8,+2.7,+2.8,+2.8,+2.8,+2.7,+2.8,+2.8,+2.8,+2.7,+2.7,+2.7,+2.7,+2.6,+2.7,+2.9,+2.8,+2.8,+2.8,+2.8,+2.8,+2.7,+2.8,+2.7,+2.8,+2.8,+2.7,+2.8,+2.7,+2.8,+2.8,+2.8,+2.7,+2.8,+2.8,+2.8,+2.8,+2.8,+2.7,+2.8,+2.8,+2.8,+2.7,+2.6,+2.7,+2.7,+2.8,+2.7,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.7,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.8,+2.9,+2.9,+2.9,+3.1,+3.1,+3.1,+3.2,+3.2,+3.4,+3.4,+3.5,+3.5,+3.5,+3.7,+3.9,+3.9,+4.0,+4.0,+4.2,+4.3,+4.4,+4.4,+4.5,+4.7,+4.8,+4.9,+5.0,+5.0,+5.2,+5.4,+5.5,+5.5,+5.6,+5.7,+5.8,+6.0,+6.1,+6.2,+6.3};
   double qqH8_errm [] = {-3.3,-3.3,-3.1,-3.0,-3.0,-3.0,-3.0,-3.2,-3.2,-3.2,-2.9,-2.9,-3.2,-2.9,-2.9,-3.1,-3.1,-3.0,-3.1,-3.1,-3.1,-3.1,-3.0,-3.0,-3.1,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-2.8,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-2.9,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-2.9,-2.9,-2.9,-2.9,-2.9,-2.9,-2.9,-2.9,-2.9,-2.9,-2.9,-2.9,-2.8,-2.9,-2.9,-2.9,-2.9,-2.9,-2.9,-2.9,-2.9,-2.9,-2.9,-2.9,-2.9,-2.9,-2.9,-2.9,-2.9,-2.9,-2.9,-2.9,-2.9,-2.8,-2.8,-2.8,-2.8,-2.8,-2.9,-2.9,-2.8,-2.8,-2.8,-2.9,-2.7,-2.8,-2.8,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.6,-2.6,-2.6,-2.6,-2.8,-2.6,-2.8,-2.8,-2.8,-2.8,-2.8,-2.7,-2.8,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.7,-2.7,-2.8,-2.8,-2.8,-2.8,-2.8,-2.8,-2.8,-2.8,-2.9,-2.8,-2.8,-2.8,-2.9,-3.0,-3.0,-3.1,-3.2,-3.2,-3.2,-3.3,-3.6,-3.6,-3.6,-3.7,-3.9,-4.1,-4.1,-4.1,-4.3,-4.3,-4.4,-4.6,-4.8,-4.8,-4.9,-5.0,-5.3,-5.2,-5.2,-5.4,-5.6,-5.7,-6.0,-6.1,-6.0,-6.2,-6.3,-6.5,-6.8,-6.9,-6.8,-7.0,-7.1,-7.2};
   double qqH8_scap [] = {+0.2,+0.4,+0.3,+0.4,+0.3,+0.3,+0.3,+0.3,+0.2,+0.3,+0.3,+0.3,+0.2,+0.3,+0.3,+0.2,+0.3,+0.3,+0.3,+0.2,+0.2,+0.3,+0.2,+0.2,+0.2,+0.3,+0.3,+0.2,+0.2,+0.3,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.3,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.3,+0.2,+0.3,+0.2,+0.2,+0.2,+0.2,+0.2,+0.3,+0.2,+0.2,+0.2,+0.3,+0.2,+0.2,+0.2,+0.3,+0.2,+0.3,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.3,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.3,+0.2,+0.2,+0.3,+0.2,+0.3,+0.2,+0.2,+0.3,+0.2,+0.2,+0.2,+0.2,+0.2,+0.3,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.3,+0.3,+0.2,+0.2,+0.2,+0.3,+0.3,+0.2,+0.3,+0.2,+0.3,+0.3,+0.3,+0.2,+0.3,+0.3,+0.3,+0.2,+0.2,+0.2,+0.2,+0.2,+0.2,+0.3,+0.2,+0.2,+0.3,+0.3,+0.3,+0.2,+0.3,+0.2,+0.3,+0.3,+0.2,+0.3,+0.2,+0.3,+0.3,+0.3,+0.2,+0.3,+0.3,+0.3,+0.3,+0.3,+0.2,+0.3,+0.3,+0.3,+0.3,+0.2,+0.3,+0.3,+0.3,+0.2,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.3,+0.4,+0.4,+0.4,+0.4,+0.4,+0.4,+0.4,+0.4,+0.4,+0.4,+0.5,+0.5,+0.5,+0.5,+0.6,+0.6,+0.6,+0.6,+0.6,+0.6,+0.6,+0.6,+0.7,+0.7,+0.7};
   double qqH8_scam [] = {-0.3,-0.3,-0.3,-0.2,-0.2,-0.2,-0.2,-0.3,-0.3,-0.3,-0.2,-0.2,-0.3,-0.2,-0.2,-0.2,-0.2,-0.3,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.1,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.1,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.1,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.1,-0.1,-0.1,-0.1,-0.1,-0.2,-0.2,-0.1,-0.1,-0.1,-0.2,-0.1,-0.2,-0.2,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.4,-0.4,-0.4,-0.4,-0.5,-0.5,-0.5,-0.5,-0.6,-0.6,-0.6,-0.6,-0.7,-0.7,-0.7,-0.8,-0.8,-0.8,-0.8,-0.8,-0.9,-0.9,-0.9,-0.9,-1.0,-1.0,-1.1,-1.0,-1.1,-1.1,-1.1,-1.2,-1.2,-1.2,-1.3,-1.2,-1.3,-1.3};
   double qqH8_pdfp [] = {+2.7,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.5,+2.6,+2.6,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.6,+2.5,+2.6,+2.5,+2.5,+2.5,+2.6,+2.6,+2.6,+2.5,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.5,+2.6,+2.5,+2.5,+2.6,+2.6,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.6,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.4,+2.5,+2.6,+2.6,+2.6,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.4,+2.4,+2.4,+2.4,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.4,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.5,+2.6,+2.6,+2.6,+2.8,+2.8,+2.8,+2.9,+2.9,+3.1,+3.1,+3.2,+3.2,+3.2,+3.4,+3.5,+3.5,+3.6,+3.6,+3.8,+3.9,+4.0,+4.0,+4.1,+4.3,+4.3,+4.4,+4.5,+4.5,+4.6,+4.8,+4.9,+4.9,+5.0,+5.1,+5.2,+5.4,+5.4,+5.5,+5.6};
   double qqH8_pdfm [] = {-3.0,-3.0,-2.8,-2.8,-2.8,-2.8,-2.8,-2.9,-2.9,-2.9,-2.7,-2.7,-2.9,-2.7,-2.7,-2.9,-2.9,-2.7,-2.9,-2.9,-2.9,-2.9,-2.8,-2.8,-2.9,-2.8,-2.8,-2.8,-2.8,-2.8,-2.8,-2.6,-2.8,-2.8,-2.8,-2.8,-2.8,-2.8,-2.8,-2.8,-2.8,-2.8,-2.8,-2.8,-2.8,-2.8,-2.8,-2.8,-2.8,-2.8,-2.8,-2.8,-2.8,-2.8,-2.8,-2.8,-2.8,-2.8,-2.8,-2.8,-2.8,-2.8,-2.8,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.7,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.5,-2.5,-2.5,-2.5,-2.7,-2.5,-2.7,-2.7,-2.7,-2.7,-2.7,-2.6,-2.7,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.5,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.6,-2.5,-2.5,-2.5,-2.6,-2.7,-2.7,-2.8,-2.8,-2.8,-2.8,-2.9,-3.1,-3.1,-3.1,-3.2,-3.3,-3.5,-3.5,-3.5,-3.6,-3.6,-3.7,-3.8,-4.0,-4.0,-4.1,-4.2,-4.4,-4.3,-4.3,-4.5,-4.6,-4.7,-4.9,-5.1,-4.9,-5.1,-5.2,-5.3,-5.6,-5.7,-5.5,-5.8,-5.8,-5.9};

   ggH7TG_xsec = new TGraph(sizeof(ggH7_mass)/sizeof(double), ggH7_mass, ggH7_xsec);
   ggH7TG_errp = new TGraph(sizeof(ggH7_mass)/sizeof(double), ggH7_mass, ggH7_errp);
   ggH7TG_errm = new TGraph(sizeof(ggH7_mass)/sizeof(double), ggH7_mass, ggH7_errm);
   ggH7TG_scap = new TGraph(sizeof(ggH7_mass)/sizeof(double), ggH7_mass, ggH7_scap);
   ggH7TG_scam = new TGraph(sizeof(ggH7_mass)/sizeof(double), ggH7_mass, ggH7_scam);
   ggH7TG_pdfp = new TGraph(sizeof(ggH7_mass)/sizeof(double), ggH7_mass, ggH7_pdfp);
   ggH7TG_pdfm = new TGraph(sizeof(ggH7_mass)/sizeof(double), ggH7_mass, ggH7_pdfm);

   qqH7TG_xsec = new TGraph(sizeof(qqH7_mass)/sizeof(double), qqH7_mass, qqH7_xsec);
   qqH7TG_errp = new TGraph(sizeof(qqH7_mass)/sizeof(double), qqH7_mass, qqH7_errp);
   qqH7TG_errm = new TGraph(sizeof(qqH7_mass)/sizeof(double), qqH7_mass, qqH7_errm);
   qqH7TG_scap = new TGraph(sizeof(qqH7_mass)/sizeof(double), qqH7_mass, qqH7_scap);
   qqH7TG_scam = new TGraph(sizeof(qqH7_mass)/sizeof(double), qqH7_mass, qqH7_scam);
   qqH7TG_pdfp = new TGraph(sizeof(qqH7_mass)/sizeof(double), qqH7_mass, qqH7_pdfp);
   qqH7TG_pdfm = new TGraph(sizeof(qqH7_mass)/sizeof(double), qqH7_mass, qqH7_pdfm);

   ggH8TG_xsec = new TGraph(sizeof(ggH8_mass)/sizeof(double), ggH8_mass, ggH8_xsec);
   ggH8TG_errp = new TGraph(sizeof(ggH8_mass)/sizeof(double), ggH8_mass, ggH8_errp);
   ggH8TG_errm = new TGraph(sizeof(ggH8_mass)/sizeof(double), ggH8_mass, ggH8_errm);
   ggH8TG_scap = new TGraph(sizeof(ggH8_mass)/sizeof(double), ggH8_mass, ggH8_scap);
   ggH8TG_scam = new TGraph(sizeof(ggH8_mass)/sizeof(double), ggH8_mass, ggH8_scam);
   ggH8TG_pdfp = new TGraph(sizeof(ggH8_mass)/sizeof(double), ggH8_mass, ggH8_pdfp);
   ggH8TG_pdfm = new TGraph(sizeof(ggH8_mass)/sizeof(double), ggH8_mass, ggH8_pdfm);

   qqH8TG_xsec = new TGraph(sizeof(qqH8_mass)/sizeof(double), qqH8_mass, qqH8_xsec);
   qqH8TG_errp = new TGraph(sizeof(qqH8_mass)/sizeof(double), qqH8_mass, qqH8_errp);
   qqH8TG_errm = new TGraph(sizeof(qqH8_mass)/sizeof(double), qqH8_mass, qqH8_errm);
   qqH8TG_scap = new TGraph(sizeof(qqH8_mass)/sizeof(double), qqH8_mass, qqH8_scap);
   qqH8TG_scam = new TGraph(sizeof(qqH8_mass)/sizeof(double), qqH8_mass, qqH8_scam);
   qqH8TG_pdfp = new TGraph(sizeof(qqH8_mass)/sizeof(double), qqH8_mass, qqH8_pdfp);
   qqH8TG_pdfm = new TGraph(sizeof(qqH8_mass)/sizeof(double), qqH8_mass, qqH8_pdfm);

   //#FIXME: extrapolated from 600 to 1TeVmissing points from 650GeV to 1TeV
   double QCDScaleMass   [] = {200, 250, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900, 1000};
   double QCDScaleK0ggH0 [] = {1.15 , 1.16, 1.17, 1.20, 1.17, 1.19, 1.22, 1.24, 1.25, 1.25, 1.25, 1.25, 1.25};
   double QCDScaleK0ggH1 [] = {0.88, 0.86, 0.84, 0.83, 0.82, 0.81, 0.80, 0.78, 0.78, 0.78, 0.78, 0.78, 0.78};
   double QCDScaleK1ggH1 [] = {1.27, 1.27, 1.27, 1.27, 1.26, 1.26, 1.25, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26};
   double QCDScaleK1ggH2 [] = {0.96, 0.96, 0.95, 0.95, 0.95, 0.95, 0.95,  0.95, 0.94, 0.94, 0.94, 0.94, 0.94};
   double QCDScaleK2ggH2 [] = { 1.20, 1.17, 1.20, 1.21, 1.20, 1.20, 1.17, 1.19, 1.19, 1.19, 1.19, 1.19, 1.19};

   double UEPSf0 []         = {0.952, 0.955, 0.958, 0.964, 0.966, 0.954, 0.946, 0.931, 0.920, 0.920, 0.920, 0.920, 0.920};
   double UEPSf1 []         = {1.055, 1.058, 1.061, 1.068, 1.078, 1.092, 1.102, 1.117, 1.121, 1.121, 1.121, 1.121, 1.121};
   double UEPSf2 []         = {1.059, 0.990, 0.942, 0.889, 0.856, 0.864, 0.868, 0.861, 0.872, 0.872, 0.872, 0.872, 0.872}; 

  TG_QCDScaleK0ggH0 = new TGraph(sizeof(QCDScaleMass)/sizeof(double), QCDScaleMass, QCDScaleK0ggH0);
  TG_QCDScaleK0ggH1 = new TGraph(sizeof(QCDScaleMass)/sizeof(double), QCDScaleMass, QCDScaleK0ggH1);
  TG_QCDScaleK1ggH1 = new TGraph(sizeof(QCDScaleMass)/sizeof(double), QCDScaleMass, QCDScaleK1ggH1);
  TG_QCDScaleK1ggH2 = new TGraph(sizeof(QCDScaleMass)/sizeof(double), QCDScaleMass, QCDScaleK1ggH2);
  TG_QCDScaleK2ggH2 = new TGraph(sizeof(QCDScaleMass)/sizeof(double), QCDScaleMass, QCDScaleK2ggH2);

  TG_UEPSf0         = new TGraph(sizeof(QCDScaleMass)/sizeof(double), QCDScaleMass, UEPSf0);
  TG_UEPSf1         = new TGraph(sizeof(QCDScaleMass)/sizeof(double), QCDScaleMass, UEPSf1);
  TG_UEPSf2         = new TGraph(sizeof(QCDScaleMass)/sizeof(double), QCDScaleMass, UEPSf2);
}



        //
        // reorder the procs to get the backgrounds; total bckg, signal, data 
        //
        void AllInfo_t::sortProc(){
           std::vector<string>bckg_procs;
           std::vector<string>sign_procs;
           bool isTotal=false, isData=false;
           for(unsigned int p=0;p<sorted_procs.size();p++){
              string procName = sorted_procs[p];
              std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
              if(it==procs.end())continue;
              if(it->first=="total"){isTotal=true; continue;}
              if(it->first=="data"){isData=true; continue;}
              if(it->second.isSign)sign_procs.push_back(procName);
              if(it->second.isBckg)bckg_procs.push_back(procName);
           }
           sorted_procs.clear();
           sorted_procs.insert(sorted_procs.end(), bckg_procs.begin(), bckg_procs.end());
           if(isTotal)sorted_procs.push_back("total");
           sorted_procs.insert(sorted_procs.end(), sign_procs.begin(), sign_procs.end());
           if(isData)sorted_procs.push_back("data");
        }


        //
        // Sum up all background processes and add this as a total process
        //
        void AllInfo_t::addProc(ProcessInfo_t& dest, ProcessInfo_t& src){
           dest.xsec = src.xsec*src.br;
           for(std::map<string, ChannelInfo_t>::iterator ch = src.channels.begin(); ch!=src.channels.end(); ch++){
              if(dest.channels.find(ch->first)==dest.channels.end()){   //this channel does not exist, create it
                 dest.channels[ch->first]         = ChannelInfo_t();
                 dest.channels[ch->first].bin     = ch->second.bin;
                 dest.channels[ch->first].channel = ch->second.channel;
              }
              std::map<string, ShapeData_t>& shapesInfoDest = dest.channels[ch->first].shapes;
              std::map<string, ShapeData_t>& shapesInfoSrc  = ch->second.shapes;
              for(std::map<string, ShapeData_t>::iterator sh = shapesInfoSrc.begin(); sh!=shapesInfoSrc.end(); sh++){
                 if(shapesInfoDest.find(sh->first)==shapesInfoDest.end())shapesInfoDest[sh->first] = ShapeData_t();

                 //only care about central histogram regarding the shape... for now               
                 if(!shapesInfoDest[sh->first].histo()){ shapesInfoDest[sh->first].uncShape[""] = (TH1*) sh->second.histo()->Clone(TString(sh->second.histo()->GetName() + dest.shortName ) );
                 }else{                                  shapesInfoDest[sh->first].histo()->Add(sh->second.histo());                                                        
                 }
                 
                 //take care of the scale uncertainty 
                 for(std::map<string, double>::iterator unc = sh->second.uncScale.begin();unc!= sh->second.uncScale.end();unc++){
                    if(shapesInfoDest[sh->first].uncScale.find(unc->first)==shapesInfoDest[sh->first].uncScale.end()){
                       shapesInfoDest[sh->first].uncScale[unc->first] = unc->second;
                    }else{
                       shapesInfoDest[sh->first].uncScale[unc->first] = sqrt( pow(shapesInfoDest[sh->first].uncScale[unc->first],2) + pow(unc->second,2) );
                    }
                 }
              }             
           }
        }



        //
        // Sum up all background processes and add this as a total process
        //
        void AllInfo_t::computeTotalBackground(){
           for(std::vector<string>::iterator p=sorted_procs.begin(); p!=sorted_procs.end();p++){if((*p)=="total"){sorted_procs.erase(p);break;}}           
           sorted_procs.push_back("total");
           procs["total"] = ProcessInfo_t(); //reset
           ProcessInfo_t& procInfo_Bckgs = procs["total"];
           procInfo_Bckgs.shortName = "data";
           procInfo_Bckgs.isData = false;
           procInfo_Bckgs.isSign = false;
           procInfo_Bckgs.isBckg = true;
           procInfo_Bckgs.xsec   = 0.0;
           procInfo_Bckgs.br     = 1.0;
           for(std::map<string, ProcessInfo_t>::iterator it=procs.begin(); it!=procs.end();it++){
              if(it->first=="total" || it->second.isBckg!=true)continue;
              addProc(procInfo_Bckgs, it->second);
           }
        }

        //
        // Replace the Data process by TotalBackground
        //
        void AllInfo_t::blind(){
           if(procs.find("total")==procs.end())computeTotalBackground();
           if(procs.find("data")==procs.end())sorted_procs.push_back("data");           
           procs["data"] =  procs["total"];
        }

        //
        // Print the Yield table
        //
         void AllInfo_t::getYieldsFromShape(FILE* pFile, std::vector<TString>& selCh, string histoName){
           std::map<string, string> rows;
           std::map<string, string> rowsBin;
           string rows_header = "\\begin{tabular}{|c|";
           string rows_title  = "channel";

           //order the proc first
           sortProc();

           for(unsigned int p=0;p<sorted_procs.size();p++){
              string procName = sorted_procs[p];
              std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
              if(it==procs.end())continue;
              rows_header += "c|";
              rows_title  += "& " + it->second.shortName;
              std::map<string, double> bin_valerr;
              std::map<string, double> bin_val;
              std::map<string, double> bin_syst;
              for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
                 if(std::find(selCh.begin(), selCh.end(), ch->second.channel)==selCh.end())continue;
                 if(ch->second.shapes.find(histoName)==(ch->second.shapes).end())continue;
                 if(rows.find(ch->first)==rows.end())rows[ch->first] = string("$ ")+ch->first+" $";
                 TH1* h = ch->second.shapes[histoName].histo();
                 double valerr;
                 double val  = h->IntegralAndError(1,h->GetXaxis()->GetNbins(),valerr);
                 double syst = ch->second.shapes[histoName].getScaleUncertainty();
                 if(val<1E-6){val=0.0; valerr=0.0; syst=-1;}
                 if(it->first=="data"){valerr=-1.0; syst=-1;}
                 rows[ch->first] += "&";
                 if(it->first=="data" || it->first=="total")rows[ch->first] += "\\boldmath ";
                 if(it->first=="data"){char tmp[256];sprintf(tmp, "%.0f", val); rows[ch->first] += tmp;
                 }else{                rows[ch->first] += utils::toLatexRounded(val,valerr, syst);     }

                 bin_val   [ch->second.bin] += val;
                 bin_valerr[ch->second.bin] += pow(valerr,2);
                 bin_syst  [ch->second.bin] += syst>=0?pow(syst,2):-1;
                 if(syst<0)bin_syst  [ch->second.bin]=-1;
              }

              for(std::map<string, double>::iterator bin=bin_val.begin(); bin!=bin_val.end(); bin++){
                 if(rowsBin.find(bin->first)==rowsBin.end())rowsBin[bin->first] = string("$ ")+bin->first+" $";
                 rowsBin[bin->first] += "&";
                 if(it->first=="data" || it->first=="total")rowsBin[bin->first] += "\\boldmath ";
                 if(it->first=="data"){char tmp[256];sprintf(tmp, "%.0f", bin_val[bin->first]); rowsBin[bin->first] += tmp;
                 }else{                rowsBin[bin->first] += utils::toLatexRounded(bin_val[bin->first],sqrt(bin_valerr[bin->first]), bin_syst[bin->first]<0?-1:sqrt(bin_syst[bin->first]));   }
              }
           }

           //All Channels
           fprintf(pFile,"\\begin{sidewaystable}[htp]\n\\begin{center}\n\\caption{Event yields expected for background and signal processes and observed in data.}\n\\label{tab:table}\n");
           fprintf(pFile, "%s}\\\\\n", rows_header.c_str());
           fprintf(pFile, "%s\\\\\n", rows_title .c_str());
           for(std::map<string, string>::iterator row = rows.begin(); row!= rows.end(); row++){
              fprintf(pFile, "%s\\\\\n", row->second.c_str());
           }
           fprintf(pFile,"\\hline\n");
           fprintf(pFile,"\\end{tabular}\n\\end{center}\n\\end{sidewaystable}\n");
           fprintf(pFile,"\n\n\n\n");

           //All Bins
           fprintf(pFile,"\\begin{sidewaystable}[htp]\n\\begin{center}\n\\caption{Event yields expected for background and signal processes and observed in data.}\n\\label{tab:table}\n");
           fprintf(pFile, "%s}\\\\\n", rows_header.c_str());
           fprintf(pFile, "%s\\\\\n", rows_title .c_str());
           for(std::map<string, string>::iterator row = rowsBin.begin(); row!= rowsBin.end(); row++){
              fprintf(pFile, "%s\\\\\n", row->second.c_str());
           }
           fprintf(pFile,"\\hline\n");
           fprintf(pFile,"\\end{tabular}\n\\end{center}\n\\end{sidewaystable}\n");
           fprintf(pFile,"\n\n\n\n");
         }

        //
        // Dump efficiencies
        //
        void AllInfo_t::getEffFromShape(FILE* pFile, std::vector<TString>& selCh, string histoName)
        {
           for(unsigned int p=0;p<sorted_procs.size();p++){
              string procName = sorted_procs[p];
              std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
              if(it==procs.end())continue;
              if(!it->second.isSign)continue;
              for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
                 if(std::find(selCh.begin(), selCh.end(), ch->second.channel)==selCh.end())continue;
                 if(ch->second.shapes.find(histoName)==(ch->second.shapes).end())continue;
                 TH1* h = ch->second.shapes[histoName].histo();
                 double valerr;
                 double val  = h->IntegralAndError(1,h->GetXaxis()->GetNbins(),valerr);
                 fprintf(pFile,"%30s %30s %4.0f %6.2E %6.2E %6.2E %6.2E\n",ch->first.c_str(), it->first.c_str(), it->second.mass, it->second.xsec, it->second.br, val/(it->second.xsec*it->second.br), valerr/(it->second.xsec*it->second.br));
              }
           }
        }

        //
        // drop control channels
        //
        void AllInfo_t::dropCtrlChannels(std::vector<TString>& selCh)
        {
           for(unsigned int p=0;p<sorted_procs.size();p++){
              string procName = sorted_procs[p];
              std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
              if(it==procs.end())continue;
              for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
                 if(std::find(selCh.begin(), selCh.end(), ch->second.channel)==selCh.end()){it->second.channels.erase(ch); ch=it->second.channels.begin();}
              }
           }
        }



        //
        // drop background process that have a negligible yield
        //
        void AllInfo_t::dropSmallBckgProc(std::vector<TString>& selCh, string histoName, double threshold)
        {
           std::map<string, double> map_yields;          
           for(unsigned int p=0;p<sorted_procs.size();p++){
              string procName = sorted_procs[p];
              std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
              if(it==procs.end())continue;
              if(!it->second.isBckg)continue;
              map_yields[it->first] = 0;
              for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
                 if(std::find(selCh.begin(), selCh.end(), ch->second.channel)==selCh.end())continue;
                 if(ch->second.shapes.find(histoName)==(ch->second.shapes).end())continue;
                 map_yields[it->first] += ch->second.shapes[histoName].histo()->Integral();
              }
           }

           double total = map_yields["total"];
           for(std::map<string, double>::iterator Y=map_yields.begin();Y!=map_yields.end();Y++){
              if(Y->second/total<threshold){
                 for(std::vector<string>::iterator p=sorted_procs.begin(); p!=sorted_procs.end();p++){if((*p)==Y->first ){sorted_procs.erase(p);break;}}
                 procs.erase(procs.find(Y->first ));
              }
           }
        }

        //
        // Make a summary plot
        //
         void AllInfo_t::showShape(std::vector<TString>& selCh , TString histoName, TString SaveName)
         {
           int NLegEntry = 0;
           std::map<string, THStack*          > map_stack;
           std::map<string, TGraphErrors*     > map_unc;
           std::map<string, TH1*              > map_data;
           std::map<string, std::vector<TH1*> > map_signals;
           std::map<string, int               > map_legend;
//           TLegend* legA  = new TLegend(0.6,0.5,0.99,0.85, "NDC");
//           TLegend* legA  = new TLegend(0.03,0.00,0.97,0.70, "NDC");
           TLegend* legA  = new TLegend(0.03,0.99,0.97,0.89, "NDC");

           //order the proc first
           sortProc();

           //loop on sorted proc
           for(unsigned int p=0;p<sorted_procs.size();p++){
              string procName = sorted_procs[p];
              std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
              if(it==procs.end())continue;
              for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
                 if(std::find(selCh.begin(), selCh.end(), ch->second.channel)==selCh.end())continue;
                 if(ch->second.shapes.find(histoName.Data())==(ch->second.shapes).end())continue;
                 TH1* h = ch->second.shapes[histoName.Data()].histo();
                 
                 if(it->first=="total"){
                    double Uncertainty = ch->second.shapes[histoName.Data()].getScaleUncertainty() / h->Integral();
                    double Maximum = 0;
                    TGraphErrors* errors = new TGraphErrors(h->GetXaxis()->GetNbins());
                    errors->SetFillStyle(3427);
                    errors->SetFillColor(kGray+1);
                    errors->SetLineStyle(1);
                    errors->SetLineColor(2);
                    int icutg=0;
                    for(int ibin=1; ibin<=h->GetXaxis()->GetNbins(); ibin++){
                        if(h->GetBinContent(ibin)>0)
                        errors->SetPoint(icutg,h->GetXaxis()->GetBinCenter(ibin), h->GetBinContent(ibin));
                        errors->SetPointError(icutg,h->GetXaxis()->GetBinWidth(ibin)/2.0, sqrt(pow(h->GetBinContent(ibin)*Uncertainty,2) + pow(h->GetBinError(ibin),2) ) );
                        Maximum =  std::max(Maximum , h->GetBinContent(ibin) + errors->GetErrorYhigh(icutg));
                        icutg++;
                    }errors->Set(icutg);
                    errors->SetMaximum(Maximum);
                    map_unc[ch->first] = errors;
                    continue;//otherwise it will fill the legend
                 }else if(it->second.isBckg){                 
                    if(map_stack.find(ch->first)==map_stack.end()){
                       map_stack[ch->first] = new THStack((ch->first+"stack").c_str(),(ch->first+"stack").c_str());                    
                    }
                    map_stack   [ch->first]->Add(h,"HIST");

                 }else if(it->second.isSign){
                    map_signals [ch->first].push_back(h);
                   
                 }else if(it->first=="data"){
                    map_data[ch->first] = h;
                 }

                 if(map_legend.find(it->first)==map_legend.end()){
                    map_legend[it->first]=1;
                    legA->AddEntry(h,it->first.c_str(),it->first=="data"?"P":it->second.isSign?"L":"F");
                    NLegEntry++;
                 }
              }
           }

           int NBins = map_data.size()/selCh.size();
           TCanvas* c1 = new TCanvas("c1","c1",300*NBins,300*selCh.size());
           c1->SetTopMargin(0.00); c1->SetRightMargin(0.00); c1->SetBottomMargin(0.00);  c1->SetLeftMargin(0.00);
           TPad* t2 = new TPad("t2","t2", 0.03, 0.90, 1.00, 1.00, -1, 1);  t2->Draw();  c1->cd();
           t2->SetTopMargin(0.00); t2->SetRightMargin(0.00); t2->SetBottomMargin(0.00);  t2->SetLeftMargin(0.00);
           TPad* t1 = new TPad("t1","t1", 0.03, 0.03, 1.00, 0.90, 4, 1);  t1->Draw();  t1->cd();
           t1->SetTopMargin(0.00); t1->SetRightMargin(0.00); t1->SetBottomMargin(0.00);  t1->SetLeftMargin(0.00);
           t1->Divide(NBins, selCh.size(), 0, 0);

           int I=1;
           for(std::map<string, THStack*>::iterator p = map_stack.begin(); p!=map_stack.end(); p++){
              //init tab
              TVirtualPad* pad = t1->cd(I);
              pad->SetTopMargin(0.06); pad->SetRightMargin(0.03); pad->SetBottomMargin(0.07);  pad->SetLeftMargin(0.06);
              pad->SetLogy(true); 

              //print histograms
              TH1* axis = (TH1*)map_data[p->first]->Clone("axis");
              axis->Reset();      
              axis->GetXaxis()->SetRangeUser(0, axis->GetXaxis()->GetXmax());
              axis->SetMinimum(1E-2);
              axis->SetMaximum(1.5*std::max( map_unc[p->first]->GetMaximum(), map_data[p->first]->GetMaximum()));
              if((I-1)%NBins!=0)axis->GetYaxis()->SetTitle("");
              axis->Draw();
              p->second->Draw("same");
              map_unc [p->first]->Draw("2 same");
              map_data[p->first]->Draw("same");
              for(unsigned int i=0;i<map_signals[p->first].size();i++){
              map_signals[p->first][i]->Draw("HIST same");
              }



              //print tab channel header
              TPaveText* Label = new TPaveText(0.1,0.81,0.94,0.89, "NDC");
              Label->SetFillColor(0);  Label->SetFillStyle(0);  Label->SetLineColor(0); Label->SetBorderSize(0);  Label->SetTextAlign(31);
              TString LabelText = procs["data"].channels[p->first].channel+"  -  "+procs["data"].channels[p->first].bin;
              LabelText.ReplaceAll("mumu","#mu#mu"); LabelText.ReplaceAll("geq2jets","#geq2jets"); LabelText.ReplaceAll("eq0jets","0jet");  LabelText.ReplaceAll("eq1jets","1jet");
              Label->AddText(LabelText);  Label->Draw();
 
              I++;
           }
           //print legend
           c1->cd(0);
           legA->SetFillColor(0); legA->SetFillStyle(0); legA->SetLineColor(0);  legA->SetBorderSize(0); legA->SetHeader("");
           legA->SetNColumns((NLegEntry/2) + 1);
           legA->Draw("same");    legA->SetTextFont(42);

           //print canvas header
           t2->cd(0);
//           TPaveText* T = new TPaveText(0.1,0.995,0.84,0.95, "NDC");
           TPaveText* T = new TPaveText(0.1,0.7,0.9,1.0, "NDC");
           T->SetFillColor(0);  T->SetFillStyle(0);  T->SetLineColor(0); T->SetBorderSize(0);  T->SetTextAlign(22);
           if(systpostfix.Contains('8')){ T->AddText("CMS preliminary, #sqrt{s}=8.0 TeV");
           }else{                         T->AddText("CMS preliminary, #sqrt{s}=7.0 TeV");
           }T->Draw();

           //save canvas
           c1->SaveAs(SaveName+"_Shape.png");
           c1->SaveAs(SaveName+"_Shape.pdf");
           c1->SaveAs(SaveName+"_Shape.C");
           delete c1;
         }


        //
        // Make a summary plot
        //
         void AllInfo_t::saveHistoForLimit(string histoName, TFile* fout){
           //order the proc first
           sortProc();

           //Loop on processes and channels
           for(unsigned int p=0;p<sorted_procs.size();p++){
              string procName = sorted_procs[p];
              std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
              if(it==procs.end())continue;
              for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
                 TString chbin = ch->first;
                 if(!fout->GetDirectory(chbin)){fout->mkdir(chbin);}fout->cd(chbin);

                 if(ch->second.shapes.find(histoName)==(ch->second.shapes).end())continue;
                 ShapeData_t& shapeInfo = ch->second.shapes[histoName];      
                 TH1* h = shapeInfo.histo();
                 shapeInfo.makeStatUnc("_CMS_hzz2l2v_", (TString("_")+ch->first+"_"+it->second.shortName+systpostfix).Data());//add stat uncertainty to the uncertainty map;

                 TString proc = it->second.shortName.c_str();
                  for(std::map<string, TH1*  >::iterator unc=shapeInfo.uncShape.begin();unc!=shapeInfo.uncShape.end();unc++){
                      TString syst   = unc->first.c_str();
                      TH1*    hshape = unc->second;
                      hshape->SetDirectory(0);


                      //If cut&count keep only 1 bin in the histo
                      if(!shape){
                         hshape = hshape->Rebin(hshape->GetXaxis()->GetNbins()); 
                         //make sure to also count the underflow and overflow
                         double bin  = hshape->GetBinContent(0) + hshape->GetBinContent(1) + hshape->GetBinContent(2);
                         double bine = sqrt(hshape->GetBinError(0)*hshape->GetBinError(0) + hshape->GetBinError(1)*hshape->GetBinError(1) + hshape->GetBinError(2)*hshape->GetBinError(2));
                         hshape->SetBinContent(0,0);              hshape->SetBinError  (0,0);
                         hshape->SetBinContent(1,bin);            hshape->SetBinError  (1,bine);
                         hshape->SetBinContent(2,0);              hshape->SetBinError  (2,0);
                      }

                      if(syst==""){
                        //central shape (for data call it data_obs)
                        hshape->SetName(proc); 
                        if(it->first=="data"){
                           hshape->Write("data_obs");
                        }else{
                           hshape->Write(proc+postfix);
                        }
                      }else if(runSystematics && proc!="data" && (syst.Contains("Up") || syst.Contains("Down"))){
                        //if empty histogram --> no variation is applied
                        if(hshape->Integral()<h->Integral()*0.01 || isnan((float)hshape->Integral())){hshape->Reset(); hshape->Add(h,1); }

                        //write variation to file
                        hshape->SetName(proc+syst);
                        hshape->Write(proc+postfix+syst);
                      }else if(runSystematics){
                        //for one sided systematics the down variation mirrors the difference bin by bin
                        hshape->SetName(proc+syst);
                        hshape->Write(proc+postfix+syst+"Up");
                        TH1 *hmirrorshape=(TH1 *)hshape->Clone(proc+syst+"Down");
                        for(int ibin=1; ibin<=hmirrorshape->GetXaxis()->GetNbins(); ibin++){
                           double bin = 2*h->GetBinContent(ibin)-hmirrorshape->GetBinContent(ibin);
                           if(bin<0)bin=0;
                           hmirrorshape->SetBinContent(ibin,bin);
                        }
                        if(hmirrorshape->Integral()<=0)hmirrorshape->SetBinContent(1, 1E-10);
                        hmirrorshape->Write(proc+postfix+syst+"Down");
                      }

                      if(runSystematics && syst!=""){
                         TString systName(syst); 
                         systName.ReplaceAll("Up",""); systName.ReplaceAll("Down","");//  systName.ReplaceAll("_","");
                         if(systName.First("_")==0)systName.Remove(0,1);

                         TH1 *temp=(TH1*) hshape->Clone();
                         temp->Add(h,-1);
                         if(temp->Integral()!=0){
                            if(shape){
                               shapeInfo.uncScale[systName.Data()]=-1.0;
                            }else{
                               double Unc = fabs(temp->Integral());
                               if(shapeInfo.uncScale.find(systName.Data())==shapeInfo.uncScale.end()){
                                  shapeInfo.uncScale[systName.Data()]=Unc;
                               }else{
                                  shapeInfo.uncScale[systName.Data()]=(shapeInfo.uncScale[systName.Data()] + Unc)/2.0;
                               }
                            }
                         }
                         delete temp;
                       }else if(syst==""){
                          shapeInfo.uncScale[syst.Data()]=hshape->Integral();
                       }
                  }
               }
               fout->cd("..");
            }
         }
     
         //
         // produce the datacards 
         //
         void AllInfo_t::buildDataCards(string histoName, TString url)
         {

           //add all scale uncertainties
           for(unsigned int p=0;p<sorted_procs.size();p++){
              string procName = sorted_procs[p];
              std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
              if(it==procs.end() || it->first=="total")continue;
              for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
                 TString chbin = ch->first;
                 if(ch->second.shapes.find(histoName)==(ch->second.shapes).end())continue;
                 ShapeData_t& shapeInfo = ch->second.shapes[histoName];      
                 double integral = shapeInfo.histo()->Integral();

                 //lumi
                 if(!it->second.isData && systpostfix.Contains('8'))shapeInfo.uncScale["lumi_8TeV"] = integral*0.05;
                 if(!it->second.isData && systpostfix.Contains('7'))shapeInfo.uncScale["lumi_7TeV"] = integral*0.022;

                 //Id+Trigger efficiencies combined
                 if(!it->second.isData){
                    if(chbin.Contains("ee"  ))  shapeInfo.uncScale["CMS_eff_e"] = integral*0.03;
                    if(chbin.Contains("mumu"))  shapeInfo.uncScale["CMS_eff_m"] = integral*0.03;
                 }

                 //uncertainties to be applied only in higgs analyses
                 if(mass>0){
                    //uncertainty on Th XSec
                    //if(it->second.isSignal)shapeInfo.uncScale["theoryUncXS_HighMH"] = std::min(1.0+1.5*pow((mass/1000.0),3),2.0);

                    if(it->second.shortName.find("ggH")!=string::npos){setTGraph(it->second.shortName, systpostfix); shapeInfo.uncScale["pdf_gg"]    = integral*0.01*max(TG_pdfp->Eval(mass,NULL,"S"), TG_pdfm->Eval(mass,NULL,"S"));}
                    if(it->second.shortName.find("qqH")!=string::npos){setTGraph(it->second.shortName, systpostfix); shapeInfo.uncScale["pdf_qqbar"] = integral*0.01*max(TG_pdfp->Eval(mass,NULL,"S"), TG_pdfm->Eval(mass,NULL,"S"));}
                    if(it->second.shortName.find("zz" )!=string::npos){                                              shapeInfo.uncScale["pdf_qqbar"] = integral*(systpostfix.Contains('8')?0.0312:0.0360);}
                  //if(it->second.shortName.find("wz" )!=string::npos){setTGraph(                                    shapeInfo.uncScale["pdf_qqbar"] = integral*(systpostfix.Contains('8')?0.0455:0.0502);}

                    //underlying events
                    if(it->second.shortName.find("ggH")!=string::npos && chbin.Contains("eq0jet" )){shapeInfo.uncScale["UEPS"] = integral*(1.0-TG_UEPSf0->Eval(mass,NULL,"S"));}
                    if(it->second.shortName.find("ggH")!=string::npos && chbin.Contains("eq1jet" )){shapeInfo.uncScale["UEPS"] = integral*(1.0-TG_UEPSf1->Eval(mass,NULL,"S"));}
                    if(it->second.shortName.find("ggH")!=string::npos && chbin.Contains("eq2jet" )){shapeInfo.uncScale["UEPS"] = integral*(1.0-TG_UEPSf2->Eval(mass,NULL,"S"));}
                    if(it->second.shortName.find("ggH")!=string::npos && chbin.Contains("vbf"    )){shapeInfo.uncScale["UEPS"] = integral*(1.0-TG_UEPSf2->Eval(mass,NULL,"S"));}

                    //bin migration at th level
                    if(it->second.shortName.find("ggH")!=string::npos && chbin.Contains("eq0jet" )){shapeInfo.uncScale["QCDscale_ggH"]    = integral*(1.0-TG_QCDScaleK0ggH0->Eval(mass,NULL,"S"));}
                    if(it->second.shortName.find("ggH")!=string::npos && chbin.Contains("eq0jet" )){shapeInfo.uncScale["QCDscale_ggH1in"] = integral*(1.0-TG_QCDScaleK0ggH1->Eval(mass,NULL,"S"));}
                    if(it->second.shortName.find("ggH")!=string::npos && chbin.Contains("eq1jet" )){shapeInfo.uncScale["QCDscale_ggH1in"] = integral*(1.0-TG_QCDScaleK1ggH1->Eval(mass,NULL,"S"));}
                    if(it->second.shortName.find("ggH")!=string::npos && chbin.Contains("eq1jet" )){shapeInfo.uncScale["QCDscale_ggH2in"] = integral*(1.0-TG_QCDScaleK1ggH2->Eval(mass,NULL,"S"));}
                    if(it->second.shortName.find("ggH")!=string::npos && chbin.Contains("eq2jet" )){shapeInfo.uncScale["QCDscale_ggH2in"] = integral*(1.0-TG_QCDScaleK2ggH2->Eval(mass,NULL,"S"));}
                    if(it->second.shortName.find("ggH")!=string::npos && chbin.Contains("vbf"    )){shapeInfo.uncScale["QCDscale_ggH2in"] = integral*(1.0-TG_QCDScaleK2ggH2->Eval(mass,NULL,"S"));}
   
                    if(it->second.shortName.find("qqH")!=string::npos                             ){shapeInfo.uncScale["QCDscale_qqH"]    = integral*0.01*sqrt(pow(TG_scap->Eval(mass,NULL,"S"),2) + pow(TG_scam->Eval(mass,NULL,"S"),2));}
                 }//end of uncertainties to be applied only in higgs analyses

                 if(it->second.shortName.find("zz" )!=string::npos){shapeInfo.uncScale["QCDscale_VV"] = integral*(systpostfix.Contains('8')?0.0669:0.0700);}          //temporary removed to avoid double counts with uncertainty applied on ZZ xsection itself --> should be reintegrated for Higgs computation
                 if(it->second.shortName.find("wz" )!=string::npos){shapeInfo.uncScale["QCDscale_VV"] = integral*(systpostfix.Contains('8')?0.0767:0.0822);}

                 if(it->second.shortName.find("ww")==0){shapeInfo.uncScale["XSec_sys_WW"] = integral*(systpostfix.Contains('8')?0.097:0.097);}
                 if(it->second.shortName.find("wz")==0){shapeInfo.uncScale["XSec_sys_WZ"] = integral*(systpostfix.Contains('8')?0.056:0.056);}
              }
           }

           //Now really take care of making the datacards



           std::vector<string>clean_procs;
           std::vector<string>sign_procs;
           //make a map of all systematics considered and say if it's shape-based or not.
           std::map<string, bool> allChannels;
           std::map<string, bool> allSysts;
           for(unsigned int p=0;p<sorted_procs.size();p++){
              string procName = sorted_procs[p];
              std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
              if(it==procs.end() || it->first=="total")continue;
              if(it->second.isSign)sign_procs.push_back(procName);
              if(it->second.isBckg)clean_procs.push_back(procName);
              for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
                 TString chbin = ch->first;
                 if(ch->second.shapes.find(histoName)==(ch->second.shapes).end())continue;
                 allChannels[ch->first] = true;
                 ShapeData_t& shapeInfo = ch->second.shapes[histoName];      
                 for(std::map<string, double>::iterator unc=shapeInfo.uncScale.begin();unc!=shapeInfo.uncScale.end();unc++){
                   if(unc->first=="")continue;
                   allSysts[unc->first] = unc->second==-1?true:false;
                 }
              }
           }
           clean_procs.insert(clean_procs.begin(), sign_procs.begin(), sign_procs.end());
           int nsign = sign_procs.size();

           TString eecard = "";
           TString mumucard = "";
           TString combinedcard = "";

           for(std::map<string, bool>::iterator C=allChannels.begin(); C!=allChannels.end();C++){
              TString dcName=url;
              dcName.ReplaceAll(".root","_"+TString(C->first.c_str())+".dat");

              combinedcard += (C->first+"=").c_str()+dcName+" ";
              if(C->first.find("ee"  )!=string::npos)eecard   += (C->first+"=").c_str()+dcName+" ";
              if(C->first.find("mumu")!=string::npos)mumucard += (C->first+"=").c_str()+dcName+" ";



              FILE* pFile = fopen(dcName.Data(),"w");
              //header
              fprintf(pFile, "imax 1\n");
              fprintf(pFile, "jmax *\n");
              fprintf(pFile, "kmax *\n");
              fprintf(pFile, "-------------------------------\n");
              if(shape){
              fprintf(pFile, "shapes * * %s %s/$PROCESS %s/$PROCESS_$SYSTEMATIC\n",url.Data(), C->first.c_str(), C->first.c_str() );
              fprintf(pFile, "-------------------------------\n");
              }
              //observations
              fprintf(pFile, "bin 1\n");
              fprintf(pFile, "Observation %f\n", procs["data"].channels[C->first].shapes[histoName].histo()->Integral() );
              fprintf(pFile, "-------------------------------\n");

              //yields
              fprintf(pFile,"%55s  ", "bin");     for(unsigned int j=0; j<clean_procs.size(); j++){ fprintf(pFile,"%8i ", 1)                     ;}  fprintf(pFile,"\n");
              fprintf(pFile,"%55s  ", "process"); for(unsigned int j=0; j<clean_procs.size(); j++){ fprintf(pFile,"%8s ", procs[clean_procs[j]].shortName.c_str());}  fprintf(pFile,"\n");
              fprintf(pFile,"%55s  ", "process"); for(unsigned int j=0; j<clean_procs.size(); j++){ fprintf(pFile,"%8i ", ((int)j)-(nsign-1)    );}  fprintf(pFile,"\n");
              fprintf(pFile,"%55s  ", "rate");    for(unsigned int j=0; j<clean_procs.size(); j++){ fprintf(pFile,"%8f ", procs[clean_procs[j]].channels[C->first].shapes[histoName].histo()->Integral() );}  fprintf(pFile,"\n");
              fprintf(pFile, "-------------------------------\n");

              for(std::map<string, bool>::iterator U=allSysts.begin(); U!=allSysts.end();U++){
                 char line[2048];
                 sprintf(line,"%-45s %-10s ", U->first.c_str(), U->second?"shapeN2":"lnN");
                 bool isNonNull = false;
                 for(unsigned int j=0; j<clean_procs.size(); j++){
                    ShapeData_t& shapeInfo = procs[clean_procs[j]].channels[C->first].shapes[histoName];
                    double integral = shapeInfo.histo()->Integral();
                    if(shapeInfo.uncScale.find(U->first)!=shapeInfo.uncScale.end()){   isNonNull = true;   
                       if(U->second)                                                   sprintf(line,"%s%8s ",line,"       1");
                       else                                                            sprintf(line,"%s%8f ",line,1+(shapeInfo.uncScale[U->first]/integral));
                    }else{                                                             sprintf(line,"%s%8s ",line,"       -");        }
                 }
                 if(isNonNull)fprintf(pFile, "%s\n", line);
              }
              fclose(pFile);
           }



           FILE* pFile = fopen("combineCards.sh","w");
           fprintf(pFile,"%s;\n",(TString("combineCards.py ") + combinedcard + " > " + "card_combined.dat").Data());
           fprintf(pFile,"%s;\n",(TString("combineCards.py ") + eecard       + " > " + "card_ee.dat").Data());
           fprintf(pFile,"%s;\n",(TString("combineCards.py ") + mumucard     + " > " + "card_mumu.dat").Data());
           fclose(pFile);         

        }


         //
         // Load histograms from root file and json to memory
         //
         void AllInfo_t::getShapeFromFile(TFile* inF, std::vector<string> channelsAndShapes, int cutBin, JSONWrapper::Object &Root,  double minCut, double maxCut, bool onlyData){
           std::vector<TString> BackgroundsInSignal;

           //iterate over the processes required
           std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
           for(unsigned int i=0;i<Process.size();i++){
               TString procCtr(""); procCtr+=i;
               TString proc=Process[i].getString("tag", "noTagFound");
               TDirectory *pdir = (TDirectory *)inF->Get(proc);         
               if(!pdir)continue;

               bool isData = Process[i].getBool("isdata", false);
               if(onlyData && !isData)continue; //just here to speedup the NRB prediction
               if(proc.Contains(")cp0"))continue; // skip those samples

               bool isSignal = Process[i].getBool("issignal", false);
               if(Process[i].getBool("spimpose", false) && (proc.Contains("ggH") || proc.Contains("qqH")))isSignal=true;
               //LQ bool isInSignal = Process[i].getBool("isinsignal", false);
               int color = Process[i].getInt("color", 1);
               int lcolor = Process[i].getInt("lcolor", color);
               int mcolor = Process[i].getInt("mcolor", color);
               int lwidth = Process[i].getInt("lwidth", 1);
               int lstyle = Process[i].getInt("lstyle", 1);
               int fill   = Process[i].getInt("fill"  , 1001);
               int marker = Process[i].getInt("marker", 20);

               double procMass=0;  char procMassStr[128] = "";
               if(isSignal &&  mass>0 && proc.Contains("H(")){
                  sscanf(proc.Data()+proc.First("H(")+2,"%lf",&procMass);
                  if(!(procMass==mass || procMass==massL || procMass==massR))continue; //skip signal sample not concerned 
                  sprintf(procMassStr,"%i",(int)procMass);
               }

                    if(isSignal && mass>0 && proc.Contains("ggH") && proc.Contains("ZZ"))proc = TString("ggH")  +procMassStr;
               else if(isSignal && mass>0 && proc.Contains("qqH") && proc.Contains("ZZ"))proc = TString("qqH")  +procMassStr;
               else if(isSignal && mass>0 && proc.Contains("ggH") && proc.Contains("WW"))proc = TString("ggHWW")+procMassStr;
               else if(isSignal && mass>0 && proc.Contains("qqH") && proc.Contains("WW"))proc = TString("qqHWW")+procMassStr;

               TString shortName = proc;
               shortName.ToLower();
               shortName.ReplaceAll(procMassStr,"");
               shortName.ReplaceAll("#bar{t}","tbar");
               shortName.ReplaceAll("z-#gamma^{*}+jets#rightarrow ll","dy");
               shortName.ReplaceAll("#rightarrow","");
               shortName.ReplaceAll("(",""); shortName.ReplaceAll(")","");    shortName.ReplaceAll("+","");    shortName.ReplaceAll(" ","");   shortName.ReplaceAll("/","");  shortName.ReplaceAll("#",""); 
               shortName.ReplaceAll("=",""); shortName.ReplaceAll(".","");    shortName.ReplaceAll("^","");    shortName.ReplaceAll("}","");   shortName.ReplaceAll("{","");  shortName.ReplaceAll(",","");
               shortName.ReplaceAll("ggh", "ggH");
               shortName.ReplaceAll("qqh", "qqH");
               if(shortName.Length()>8)shortName.Resize(8);

               if(procs.find(proc.Data())==procs.end()){sorted_procs.push_back(proc.Data());}
               ProcessInfo_t& procInfo = procs[proc.Data()];
               procInfo.jsonObj = Process[i]; 
               procInfo.isData = isData;
               procInfo.isSign = isSignal;
               procInfo.isBckg = !procInfo.isData && !procInfo.isSign;
               procInfo.mass   = procMass;
               procInfo.shortName = shortName.Data();

               procInfo.xsec = procInfo.jsonObj["data"].daughters()[0].getDouble("xsec", 1);
               if(procInfo.jsonObj["data"].daughters()[0].isTag("br")){
                  std::vector<JSONWrapper::Object> BRs = procInfo.jsonObj["data"].daughters()[0]["br"].daughters();
                  double totalBR=1.0; for(size_t ipbr=0; ipbr<BRs.size(); ipbr++){totalBR*=BRs[ipbr].toDouble();}   
                  procInfo.br = totalBR;
               }               

               //Loop on all channels, bins and shape to load and store them in memory structure
               TH1* syst = (TH1*)pdir->Get("optim_systs");
               if(syst==NULL){syst=new TH1F("optim_systs","optim_systs",1,0,1);syst->GetXaxis()->SetBinLabel(1,"");}
               for(unsigned int c=0;c<channelsAndShapes.size();c++){
               TString chName    = (channelsAndShapes[c].substr(0,channelsAndShapes[c].find(";"))).c_str();
               TString binName   = (channelsAndShapes[c].substr(channelsAndShapes[c].find(";")+1, channelsAndShapes[c].rfind(";")-channelsAndShapes[c].find(";")-1)).c_str();
               TString shapeName = (channelsAndShapes[c].substr(channelsAndShapes[c].rfind(";")+1)).c_str();
               TString ch        = chName+binName;

               ChannelInfo_t& channelInfo = procInfo.channels[ch.Data()];
               channelInfo.bin        = binName.Data();
               channelInfo.channel    = chName.Data();
               ShapeData_t& shapeInfo = channelInfo.shapes[shapeName.Data()];

               for(int ivar = 1; ivar<=syst->GetNbinsX();ivar++){                 
                  TH1D* hshape   = NULL;
                  TString varName   = syst->GetXaxis()->GetBinLabel(ivar);
                  TString histoName = ch+"_"+shapeName+(isSignal?signalSufix:"")+varName ;
                  if(shapeName==histo && histoVBF!="" && ch.Contains("vbf"))histoName = ch+"_"+histoVBF+(isSignal?signalSufix:"")+varName ;
                  TH2* hshape2D = (TH2*)pdir->Get(histoName );
                  if(!hshape2D){
                     if(shapeName==histo && histoVBF!="" && ch.Contains("vbf")){   hshape2D = (TH2*)pdir->Get(histoVBF+(isSignal?signalSufix:"")+varName);
                     }else{                                                        hshape2D = (TH2*)pdir->Get(shapeName+varName);
                     }
                                
                     if(hshape2D){
                        hshape2D->Reset();
                     }else{  //if still no histo, skip this proc...
                        //printf("Histo %s does not exist for syst:%s\n", histoName.Data(), varName.Data());
                        continue;
                     }
                  }

                  //special treatment for side mass points
                  int cutBinUsed = cutBin;
                  if(shapeName == histo && !ch.Contains("vbf") && procMass==massL)cutBinUsed = indexcutL;
                  if(shapeName == histo && !ch.Contains("vbf") && procMass==massR)cutBinUsed = indexcutL;

                  
                  histoName.ReplaceAll(ch,ch+"_proj"+procCtr);
                  hshape   = hshape2D->ProjectionY(histoName,cutBinUsed,cutBinUsed);
                  if(hshape->Integral()<=0 && varName=="" && !isData){hshape->Reset(); hshape->SetBinContent(1, 1E-10);}

                  if(isnan((float)hshape->Integral())){hshape->Reset();}
                  hshape->SetDirectory(0);
                  hshape->SetTitle(proc);
                  utils::root::fixExtremities(hshape,true,true);
                  hshape->SetFillColor(color); hshape->SetLineColor(lcolor); hshape->SetMarkerColor(mcolor);
                  hshape->SetFillStyle(fill);  hshape->SetLineWidth(lwidth); hshape->SetMarkerStyle(marker); hshape->SetLineStyle(lstyle);

                  //if current shape is the one to cut on, then apply the cuts
                  if(shapeName == histo){
                     for(int x=0;x<=hshape->GetXaxis()->GetNbins()+1;x++){
                        if(hshape->GetXaxis()->GetBinCenter(x)<=minCut || hshape->GetXaxis()->GetBinCenter(x)>=maxCut){ hshape->SetBinContent(x,0); hshape->SetBinError(x,0); }
                     }
                     if(histoVBF!="" && ch.Contains("vbf")){    hshape->Rebin(5);
                     }else{                                     hshape->Rebin(rebinVal);
                     }
                     hshape->GetYaxis()->SetTitle("Entries (/25GeV)");
                  }
                  hshape->Scale(MCRescale);


                   //Do Renaming and cleaning
                   varName.ReplaceAll("down","Down");
                   varName.ReplaceAll("up","Up");
 
                         if(varName==""){//does nothing
                   }else if(varName.BeginsWith("_jes")){varName.ReplaceAll("_jes","_CMS_scale_j");
                   }else if(varName.BeginsWith("_jer")){varName.ReplaceAll("_jer","_CMS_res_j"); // continue;//skip res for now
                   }else if(varName.BeginsWith("_les")){ 
                         if(ch.Contains("ee"  ))varName.ReplaceAll("_les","_CMS_scale_e");
                         if(ch.Contains("mumu"))varName.ReplaceAll("_les","_CMS_scale_m");
                   }else if(varName.BeginsWith("_btag"  )){varName.ReplaceAll("_btag","_CMS_eff_b"); 
                   }else if(varName.BeginsWith("_pu"    )){varName.ReplaceAll("_pu", "_CMS_hzz2l2v_pu");
                   }else if(varName.BeginsWith("_ren"   )){continue;   //already accounted for in QCD scales
                   }else if(varName.BeginsWith("_fact"  )){continue; //skip this one
                   }else if(varName.BeginsWith("_interf")){varName="_CMS_hzz2l2v"+varName;
                   }else{                               varName="_CMS_hzz2l2v"+varName;
                   }

                   hshape->SetTitle(proc+varName);
                   if(shapeInfo.uncShape.find(varName.Data())==shapeInfo.uncShape.end()){
                      shapeInfo.uncShape[varName.Data()] = hshape;
                   }else{
                      shapeInfo.uncShape[varName.Data()]->Add(hshape);
                   }
                }
            }
           }
         }

         //
         // replace MC NonResonnant Backgrounds by DataDriven estimate
         //
         void AllInfo_t::doBackgroundSubtraction(FILE* pFile, std::vector<TString>& selCh,TString ctrlCh,TString mainHisto, TString sideBandHisto)
         {
              char Lcol     [1024] = "";
              char Lchan    [1024] = "";
              char Lalph1   [1024] = "";
              char Lalph2   [1024] = "";
              char Lyield   [1024] = "";
              char LyieldMC [1024] = "";

              //check that the data proc exist
              std::map<string, ProcessInfo_t>::iterator dataProcIt=procs.find("data");             
              if(dataProcIt==procs.end()){printf("The process 'data' was not found... can not do non-resonnant background prediction\n"); return;}

              //create a new proc for NRB backgrounds
              TString NRBProcName = "Top/W/WW";
              for(std::vector<string>::iterator p=sorted_procs.begin(); p!=sorted_procs.end();p++){if((*p)==NRBProcName.Data()){sorted_procs.erase(p);break;}}           
              sorted_procs.push_back(NRBProcName.Data());
              procs[NRBProcName.Data()] = ProcessInfo_t(); //reset
              ProcessInfo_t& procInfo_NRB = procs[NRBProcName.Data()];
              procInfo_NRB.shortName = "topwww";
              procInfo_NRB.isData = true;
              procInfo_NRB.isSign = false;
              procInfo_NRB.isBckg = true;
              procInfo_NRB.xsec   = 0.0;
              procInfo_NRB.br     = 1.0;

              //create an histogram containing all the MC backgrounds
              std::vector<string> toBeDelete;
              for(std::map<string, ProcessInfo_t>::iterator it=procs.begin(); it!=procs.end();it++){
                 if(!it->second.isBckg || it->second.isData)continue;
                 TString procName = it->first.c_str();
                 if(!(procName.Contains("t#bar{t}") || procName.Contains("Single top") || procName.Contains("WW") || procName.Contains("Z#rightarrow #tau#tau") || procName.Contains("W#rightarrow l#nu")) )continue;
                 addProc(procInfo_NRB, it->second);
                 for(std::vector<string>::iterator p=sorted_procs.begin(); p!=sorted_procs.end();p++){if((*p)==it->first){sorted_procs.erase(p);break;}}
                 toBeDelete.push_back(it->first);
              }
              for(std::vector<string>::iterator p=toBeDelete.begin();p!=toBeDelete.end();p++){procs.erase(procs.find((*p)));}


              for(std::map<string, ChannelInfo_t>::iterator chData = dataProcIt->second.channels.begin(); chData!=dataProcIt->second.channels.end(); chData++){            
                 if(std::find(selCh.begin(), selCh.end(), chData->second.channel)==selCh.end())continue;

                 std::map<string, ChannelInfo_t>::iterator chNRB  = procInfo_NRB.channels.find(chData->first);  
                 if(chNRB==procInfo_NRB.channels.end()){  //this channel does not exist, create it
                    procInfo_NRB.channels[chData->first] = ChannelInfo_t();     
                    chNRB                = procInfo_NRB.channels.find(chData->first);
                    chNRB->second.bin     = chData->second.bin;
                    chNRB->second.channel = chData->second.channel;
                 }

                 //load data histogram in the control channel
                 TH1* hCtrl_SB = dataProcIt->second.channels[(ctrlCh+chData->second.bin.c_str()).Data()].shapes[sideBandHisto.Data()].histo();
                 TH1* hCtrl_SI = dataProcIt->second.channels[(ctrlCh+chData->second.bin.c_str()).Data()].shapes[mainHisto    .Data()].histo();
                 TH1* hChan_SB =                    chData->second                                      .shapes[sideBandHisto.Data()].histo();
                 TH1* hNRB     =                    chNRB ->second                                      .shapes[mainHisto    .Data()].histo();

                 //compute alpha
                 double alpha=0 ,alpha_err=0;
                 double alphaUsed=0 ,alphaUsed_err=0;
                 if(hCtrl_SB->GetBinContent(5)>0){
                    alpha     = hChan_SB->GetBinContent(5) / hCtrl_SB->GetBinContent(5);
                    alpha_err = ( fabs( hChan_SB->GetBinContent(5) * hCtrl_SB->GetBinError(5) ) + fabs(hChan_SB->GetBinError(5) * hCtrl_SB->GetBinContent(5) )  ) / pow(hCtrl_SB->GetBinContent(5), 2);        
                 }
                 if(chData->second.channel.find("ee"  )==0){alphaUsed = 0.44; alphaUsed_err=0.03;}
                 if(chData->second.channel.find("mumu")==0){alphaUsed = 0.71; alphaUsed_err=0.04;}

                 double valDD, valDD_err;
                 double valMC, valMC_err;
                 valMC = hNRB->IntegralAndError(1,hNRB->GetXaxis()->GetNbins(),valMC_err);  if(valMC<1E-6){valMC=0.0; valMC_err=0.0;}

                 //for VBF stat in emu is too low, so take the shape from MC and scale it to the expected yield
                 if(chData->second.bin.find("vbf")==0){
                    hNRB->Scale(hCtrl_SI->Integral(1, hCtrl_SI->GetXaxis()->GetNbins()+1) / hNRB->Integral(1,hNRB->GetXaxis()->GetNbins()+1));
                 }else{
                    hNRB->Reset();
                    hNRB->Add(hCtrl_SI , 1.0);
                 }
                 for(int bi=1;bi<=hNRB->GetXaxis()->GetNbins()+1;bi++){
                    double val = hNRB->GetBinContent(bi);
                    double err = hNRB->GetBinError(bi);
                    double newval = val*alphaUsed;
                    double newerr = sqrt(pow(err*alphaUsed,2) + pow(val*alphaUsed_err,2));
                    hNRB->SetBinContent(bi, newval );
                    hNRB->SetBinError  (bi, newerr );
                 }
                 hNRB->Scale(DDRescale);
                 hNRB->SetTitle(NRBProcName.Data());
                 hNRB->SetFillStyle(1001);
                 hNRB->SetFillColor(592);

                 //save values for printout
                 valDD = hNRB->IntegralAndError(1,hNRB->GetXaxis()->GetNbins()+1,valDD_err); if(valDD<1E-6){valDD=0.0; valDD_err=0.0;}

                 //add syst uncertainty
                 chNRB->second.shapes[mainHisto.Data()].uncScale[string("CMS_hzz2l2v_sys_topwww") + systpostfix.Data()] = valDD*NonResonnantSyst;

                 //printout
                 sprintf(Lcol    , "%s%s"  ,Lcol,    "|c");
                 sprintf(Lchan   , "%s%25s",Lchan,   (string(" &") + chData->second.channel+string(" - ")+chData->second.bin).c_str());
                 sprintf(Lalph1  , "%s%25s",Lalph1,  (string(" &") + utils::toLatexRounded(alpha,alpha_err)).c_str());
                 sprintf(Lalph2  , "%s%25s",Lalph2,  (string(" &") + utils::toLatexRounded(alphaUsed,alphaUsed_err)).c_str());
                 sprintf(Lyield  , "%s%25s",Lyield,  (string(" &") + utils::toLatexRounded(valDD,valDD_err,valDD*NonResonnantSyst)).c_str());
                 sprintf(LyieldMC, "%s%25s",LyieldMC,(string(" &") + utils::toLatexRounded(valMC,valMC_err)).c_str());
              }

              //recompute total background
              computeTotalBackground();

              if(pFile){
                 fprintf(pFile,"\\begin{table}[htp]\n\\begin{center}\n\\caption{Non resonant background estimation.}\n\\label{tab:table}\n");
                 fprintf(pFile,"\\begin{tabular}{%s|}\\hline\n", Lcol);
                 fprintf(pFile,"channel               %s\\\\\n", Lchan);
                 fprintf(pFile,"$\\alpha$ measured    %s\\\\\n", Lalph1);
                 fprintf(pFile,"$\\alpha$ used        %s\\\\\n", Lalph2);
                 fprintf(pFile,"yield data            %s\\\\\n", Lyield);
                 fprintf(pFile,"yield mc              %s\\\\\n", LyieldMC);
                 fprintf(pFile,"\\hline\n");
                 fprintf(pFile,"\\end{tabular}\n\\end{center}\n\\end{table}\n");
              }
         }

         //
         // replace MC Z+Jets Backgrounds by DataDriven Gamma+Jets estimate
         //
         void AllInfo_t::doDYReplacement(FILE* pFile, std::vector<TString>& selCh,TString ctrlCh,TString mainHisto){
           TString DYProcName = "Z#rightarrow ll";
           TString GammaJetProcName = "Instr. background (data)";
           std::map<TString, double> LowMetIntegral;
           std::vector<string> lineprintouts;

           //open gamma+jet file
           TFile* inF = TFile::Open(DYFile);
           if( !inF || inF->IsZombie() ){ cout << "Invalid file name : " << DYFile << endl; return; }           
           TDirectory* pdir = (TDirectory *)inF->Get(GammaJetProcName);
           if(!pdir){ printf("Skip Z+Jet estimation because %s directory is missing in Gamma+Jet file\n", GammaJetProcName.Data()); return;}
           gROOT->cd(); //make sure that all histograms that will be created will be in memory and not in file

           //check that the z+jets proc exist
           std::map<string, ProcessInfo_t>::iterator mcZJetIt=procs.find(DYProcName.Data());
           if(mcZJetIt==procs.end()){printf("The process '%s' was not found... can not do DY background prediction\n",DYProcName.Data()); return;}

           //create a new proc for Z+Jets datadriven backgrounds as a copy of the MC one
           TString ZJetProcName = "Z+Jets";
           for(std::vector<string>::iterator p=sorted_procs.begin(); p!=sorted_procs.end();p++){if((*p)==ZJetProcName.Data()){sorted_procs.erase(p);break;}}
           sorted_procs.push_back(ZJetProcName.Data());
           procs[ZJetProcName.Data()] = ProcessInfo_t(); //reset
           ProcessInfo_t& procInfo_dataZJets = procs[ZJetProcName.Data()];
           procInfo_dataZJets.shortName = "zll";
           procInfo_dataZJets.isData = true;
           procInfo_dataZJets.isSign = false;
           procInfo_dataZJets.isBckg = true;
           procInfo_dataZJets.xsec   = mcZJetIt->second.xsec;
           procInfo_dataZJets.br     = mcZJetIt->second.br;
           addProc(procInfo_dataZJets, mcZJetIt->second);

           for(std::map<string, ChannelInfo_t>::iterator chMC = mcZJetIt->second.channels.begin(); chMC!=mcZJetIt->second.channels.end(); chMC++){
              if(std::find(selCh.begin(), selCh.end(), chMC->second.channel)==selCh.end())continue;

              std::map<string, ChannelInfo_t>::iterator chData  = procInfo_dataZJets.channels.find(chMC->first);
              if(chData==procInfo_dataZJets.channels.end()){  //this channel does not exist, create it
                 procInfo_dataZJets.channels[chMC->first] = ChannelInfo_t();
                 chData                 = procInfo_dataZJets.channels.find(chMC->first);
                 chData->second.bin     = chMC->second.bin;
                 chData->second.channel = chMC->second.channel;
              }

              //load G+Jets data
              int indexcut_ = indexcut; double cutMin=shapeMin; double cutMax=shapeMax;
              if(indexvbf>=0 && chMC->second.bin.find("vbf")!=string::npos){indexcut_ = indexvbf; cutMin=shapeMinVBF; cutMax=shapeMaxVBF;}

              TH2* gjets2Dshape = NULL;
              if(mainHisto==histo && histoVBF!="" && chMC->second.bin.find("vbf")!=string::npos){
                 gjets2Dshape  = (TH2*)pdir->Get(((chMC->second.channel+chMC->second.bin+"_")+histoVBF.Data()).c_str());
              }else{
                 gjets2Dshape  = (TH2*)pdir->Get(((chMC->second.channel+chMC->second.bin+"_")+mainHisto.Data()).c_str());
              }
              if(!gjets2Dshape)printf("Can't find histo: %s in g+jets template\n",((chMC->second.channel+chMC->second.bin+"_")+mainHisto.Data()).c_str());

              TH1* hMC = chMC->second.shapes[mainHisto.Data()].histo();
              TH1* hDD = gjets2Dshape->ProjectionY("tmpName",indexcut_,indexcut_);
               utils::root::fixExtremities(hDD, true, true);
               if(!(mainHisto==histo && histoVBF!="" && chMC->second.bin.find("vbf")!=string::npos)){
                  for(int x=0;x<=hDD->GetXaxis()->GetNbins()+1;x++){
                     if(hDD->GetXaxis()->GetBinCenter(x)<=cutMin || hDD->GetXaxis()->GetBinCenter(x)>=cutMax){hDD->SetBinContent(x,0); hDD->SetBinError(x,0);}
                  }
               }

              //Check the binning!!!
              if(hDD->GetXaxis()->GetXmin()!=hMC->GetXaxis()->GetXmin()){printf("Gamma+Jet templates have a different XAxis range\nStop the script here\n"); exit(0);}
              if(hDD->GetXaxis()->GetBinWidth(1)!=hMC->GetXaxis()->GetBinWidth(1)){
                 double dywidth = hDD->GetXaxis()->GetBinWidth(1);
                 printf("Gamma+Jet templates have a different bin width:");
                 double mcwidth = hMC->GetXaxis()->GetBinWidth(1);
                 if(dywidth>mcwidth){printf("bin width in Gamma+Jet templates is larger than in MC samples --> can not rebin!\nStop the script here\n"); exit(0);}
                 int rebinfactor = (int)(mcwidth/dywidth);
                 if(((int)mcwidth)%((int)dywidth)!=0){printf("bin width in Gamma+Jet templates are not multiple of the mc histograms bin width\n"); exit(0);}
                 printf("Rebinning by %i --> ", rebinfactor);
                 hDD->Rebin(rebinfactor);
                 printf("Binning DataDriven ZJets Min=%7.2f  Max=%7.2f Width=%7.2f compared to MC ZJets Min=%7.2f  Max=%7.2f Width=%7.2f\n", hDD->GetXaxis()->GetXmin(), hDD->GetXaxis()->GetXmax(), hDD->GetXaxis()->GetBinWidth(1), hMC->GetXaxis()->GetXmin(), hMC->GetXaxis()->GetXmax(), hMC->GetXaxis()->GetBinWidth(1));
              }

              //save histogram to the structure
              hDD->Scale(DDRescale);
              chData->second.shapes[mainHisto.Data()].histo()->Reset();
              chData->second.shapes[mainHisto.Data()].histo()->Add(hDD);

              //printouts
              char printout[2048];
              double valMC_err, valMC = hMC->IntegralAndError(1,hMC->GetXaxis()->GetNbins()+1,valMC_err); if(valMC<1E-6){valMC=0.0; valMC_err=0.0;}
              double valDD_err, valDD = hDD->IntegralAndError(1,hDD->GetXaxis()->GetNbins()+1,valDD_err); if(valDD<1E-6){valDD=0.0; valDD_err=0.0;}
              sprintf(printout,"%20s & %30s & %30s\\\\", chMC->first.c_str(), utils::toLatexRounded(valDD,valDD_err,valDD*GammaJetSyst).c_str(), utils::toLatexRounded(valMC,valMC_err).c_str() );
              lineprintouts.push_back(printout);

              //add syst uncertainty
              chData->second.shapes[mainHisto.Data()].uncScale[string("CMS_hzz2l2v_sys_zll") + systpostfix.Data()] = valDD*GammaJetSyst;
 
              //clean
              delete hDD;
              delete gjets2Dshape;
           }

           //all done with gamma+jet file
           inF->Close();

           //delete MC Z+Jets proc
           for(std::vector<string>::iterator p=sorted_procs.begin(); p!=sorted_procs.end();p++){if((*p)==mcZJetIt->first){sorted_procs.erase(p);break;}}
           procs.erase(mcZJetIt);

           //recompute total background
           computeTotalBackground();

           //printouts 
           if(pFile){
              fprintf(pFile,"\\begin{table}[htp]\n\\begin{center}\n\\caption{Instrumental background estimation.}\n\\label{tab:table}\n");
              fprintf(pFile,"\\begin{tabular}{|l|c|c|c|}\\hline\n");
              fprintf(pFile,"channel & rescale & yield data & yield mc\\\\\\hline\n");
              for(unsigned int i=0;i<lineprintouts.size();i++){fprintf(pFile,"%s\n",lineprintouts[i].c_str()); }
              fprintf(pFile,"\\hline\n");
              fprintf(pFile,"\\end{tabular}\n\\end{center}\n\\end{table}\n");
              fprintf(pFile,"\n\n\n\n");
           }      
         }


         //
         // Rebin histograms to make sure that high mt/met region have no empty bins
         //
         void AllInfo_t::rebinMainHisto(string histoName)
         {
           //Loop on processes and channels
           for(unsigned int p=0;p<sorted_procs.size();p++){
              string procName = sorted_procs[p];
              std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
              if(it==procs.end())continue;
              for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
                 if(ch->second.shapes.find(histoName)==(ch->second.shapes).end())continue;
                 ShapeData_t& shapeInfo = ch->second.shapes[histoName];      
                 for(std::map<string, TH1*  >::iterator unc=shapeInfo.uncShape.begin();unc!=shapeInfo.uncShape.end();unc++){
                     TH1* histo = unc->second;
                     double* xbins = NULL;  int nbins=0;
                     if(ch->first.find("vbf")!=string::npos){
                        nbins = 3;
                        xbins = new double[nbins+1];
                        xbins[0] =  0;
                        xbins[1] = 70;
                        xbins[2] = histo->GetXaxis()->GetBinLowEdge(histo->GetXaxis()->FindBin(150));
                        xbins[3] = histo->GetXaxis()->GetBinLowEdge(histo->GetNbinsX()+1);
                     }else{
                        nbins = histo->GetXaxis()->FindBin(1000)+1;
                        xbins = new double[nbins+1];
                        for(int x=0;x<=histo->GetXaxis()->FindBin(1000)+1;x++){xbins[x]=histo->GetXaxis()->GetBinLowEdge(x);}
                        xbins[nbins] = histo->GetXaxis()->GetBinLowEdge(histo->GetNbinsX()+1);
                     }
                     unc->second = histo->Rebin(nbins, histo->GetName(), (double*)xbins);
                 }
              }
           }
         }


         //
         // Interpollate the signal sample between two mass points 
         //
         void AllInfo_t::SignalInterpolation(string histoName){
           if(massL<0 || massR<0 || massL==massR)return;

           //Loop on processes and identify what is available
           std::map<string, std::pair<string, string> > procLeftRight;
           for(unsigned int p=0;p<sorted_procs.size();p++){
              string procName = sorted_procs[p];
              std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
              if(it==procs.end() || !it->second.isSign || it->second.mass<=0)continue;

              TString massStr = ""; massStr+=(int)it->second.mass;
              TString procNameWithoutMass = it->first;
              procNameWithoutMass.ReplaceAll(massStr,"");
              if(it->second.mass == massL)procLeftRight[procNameWithoutMass.Data()].first  = it->first;
              if(it->second.mass == massR)procLeftRight[procNameWithoutMass.Data()].second = it->first;
           }
         
           double Ratio = ((double)mass - massL); Ratio/=(massR - massL); 
           for(std::map<string, std::pair<string, string> >::iterator procLR = procLeftRight.begin(); procLR!=procLeftRight.end();procLR++){
              TString signProcName =  procLR->first; signProcName+=(int)mass;
              printf("Interpolate %s  based on %s and %s\n", signProcName.Data(), procLR->second.first.c_str(), procLR->second.second.c_str());

              setTGraph(signProcName.Data(), systpostfix ); //needed to eval the xsec

              ProcessInfo_t& procL = procs[ procLR->second.first ];
              ProcessInfo_t& procR = procs[ procLR->second.second];

              //create a new proc as a copy of the Left signal proc
              for(std::vector<string>::iterator p=sorted_procs.begin(); p!=sorted_procs.end();p++){if((*p)==signProcName.Data()){sorted_procs.erase(p);break;}}
              sorted_procs.push_back(signProcName.Data());
              procs[signProcName.Data()] = ProcessInfo_t(); //reset
              ProcessInfo_t& proc = procs[signProcName.Data()];
              addProc(proc, procL);
              proc.shortName = procL.shortName;
              proc.isData    = procL.isData;
              proc.isSign    = procL.isSign;
              proc.isBckg    = procL.isBckg;
              proc.mass      = mass;
              proc.xsec      = TG_xsec->Eval(proc.mass,NULL,"S");
              proc.br        = procL.br + (Ratio * (procR.br - procL.br));

              for(std::map<string, ChannelInfo_t>::iterator ch  = proc .channels.begin(); ch !=proc.channels.end(); ch++){                 
                  std::map<string, ChannelInfo_t>::iterator chL = procL.channels.find(ch->first);
                  std::map<string, ChannelInfo_t>::iterator chR = procR.channels.find(ch->first);
                  if(chL==procL.channels.end())continue; 
                  if(chR==procR.channels.end())continue;
                  if(ch ->second.shapes.find(histoName)==(ch ->second.shapes).end())continue;
                  if(chL->second.shapes.find(histoName)==(chL->second.shapes).end())continue;
                  if(chR->second.shapes.find(histoName)==(chR->second.shapes).end())continue;
                  ShapeData_t& shapeInfo  = ch ->second.shapes[histoName];
                  ShapeData_t& shapeInfoL = chL->second.shapes[histoName];
                  ShapeData_t& shapeInfoR = chR->second.shapes[histoName];

                 for(std::map<string, TH1*  >::iterator unc=shapeInfoL.uncShape.begin();unc!=shapeInfoL.uncShape.end();unc++){
                    if(shapeInfo.uncShape.find(unc->first)==shapeInfo.uncShape.end())shapeInfo.uncShape[unc->first] = (TH1*)shapeInfoL.uncShape[unc->first]->Clone(signProcName+ch->first+unc->first+"tmp");
                    TH1D* hL = (TH1D*)shapeInfoL.uncShape[unc->first];
                    TH1D* hR = (TH1D*)shapeInfoR.uncShape[unc->first];
                    TH1D* h  = (TH1D*)shapeInfo .uncShape[unc->first];
                    hL->Scale(1.0/(procL.xsec*procL.br));
                    hR->Scale(1.0/(procR.xsec*procR.br));
                    h->Reset();
                    if(hL->Integral()>0 && hR->Integral()>0)//interpolate only if the histograms are not null
                    h->Add(th1fmorph(signProcName+ch->first+unc->first,signProcName+ch->first+unc->first, hL, hR, procL.mass, procR.mass, proc.mass, (1-Ratio)*hL->Integral() + Ratio*hR->Integral(), 0), 1);
                    h->Scale(proc.xsec*proc.br);
                 }
              }
              
              //erase sideband procs
              for(std::vector<string>::iterator p=sorted_procs.begin(); p!=sorted_procs.end();p++){if((*p)==procLR->second.first ){sorted_procs.erase(p);break;}}              
              for(std::vector<string>::iterator p=sorted_procs.begin(); p!=sorted_procs.end();p++){if((*p)==procLR->second.second){sorted_procs.erase(p);break;}}       
              procs.erase(procs.find(procLR->second.first ));
              procs.erase(procs.find(procLR->second.second ));
           }
         }



         //
         // Rescale signal sample for the effect of the interference and propagate the uncertainty 
         //
         void AllInfo_t::RescaleForInterference(string histoName){
           if(mass<=0)return;

           //Loop on processes and channels
           for(unsigned int p=0;p<sorted_procs.size();p++){
              string procName = sorted_procs[p];
              std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
              if(it==procs.end() || !it->second.isSign)continue;

              double sF   = 1.0;
              double sFDn = 1.0;
              double sFUp = 1.0;

              if(doInterf && it->second.mass>400 && it->first.find("ggH")!=string::npos){
                 sF   = 0.897-0.000152*it->second.mass+7.69e-07*pow(it->second.mass,2);
                 sFDn = 0.907-2.08e-05*it->second.mass+4.63e-07*pow(it->second.mass,2);
                 sFUp = 0.889-0.000357*it->second.mass+1.21e-06*pow(it->second.mass,2);

                 if(it->second.mass>=400 && signalSufix!=""){ //scale factor for Narrow Resonnance
                    double cprime=1.0; double  brnew=0.0;
                    sscanf(signalSufix.Data(), "_cp%lf_brn%lf", &cprime, &brnew);
                    sF=1 + (sF-1)/pow(cprime,2);   sFDn = 1.0;  sFUp = 1 + (sF-1)*2;        //100% Uncertainty
                    if(sF<1){sF=1.0;  sFUp = 1 + (sF-1)*2;}
                    printf("Scale Factor for Narrow Resonnance : %f [%f,%f] applied on %s\n", sF, sFDn, sFUp, it->first.c_str());
                 }else{
                    printf("Scale Factor for Interference : %f [%f,%f] applied on %s\n",sF, sFDn, sFUp, it->first.c_str());
                 }
                 if(sFDn>sFUp){double tmp = sFUp; sFUp = sFDn; sFDn = tmp;}
              }
              printf("Total Scale Factor : %f [%f,%f] applied on %s\n",sF, sFDn, sFUp, it->first.c_str());
              it->second.xsec *= sF;
   
              for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
                 if(ch->second.shapes.find(histoName)==(ch->second.shapes).end())continue;
                 ShapeData_t& shapeInfo = ch->second.shapes[histoName];
                 TH1* histo = (TH1*)shapeInfo.histo();
                 for(std::map<string, TH1*  >::iterator unc=shapeInfo.uncShape.begin();unc!=shapeInfo.uncShape.end();unc++){
                     unc->second->Scale(sF);
                 }
                 TH1* tmp;
                 tmp = (TH1*)histo->Clone(TString("interf_ggH_") + histo->GetName() + "Down"); tmp->Scale(sFDn/sF); shapeInfo.uncShape[string("_CMS_hzz2l2v_interf_") + it->second.shortName+"Down"] = tmp;
                 tmp = (TH1*)histo->Clone(TString("interf_ggH_") + histo->GetName() + "Up"  ); tmp->Scale(sFUp/sF); shapeInfo.uncShape[string("_CMS_hzz2l2v_interf_") + it->second.shortName+"Up"  ] = tmp;
              }
           }
         }

