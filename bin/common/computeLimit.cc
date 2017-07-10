// Original Author:  Loic Quertenmont

#include <iostream>
#include <boost/shared_ptr.hpp>
#include "Math/GenVector/Boost.h"

#include "UserCode/llvv_fwk/interface/tdrstyle.h"
#include "UserCode/llvv_fwk/interface/JSONWrapper.h"
#include "UserCode/llvv_fwk/interface/RootUtils.h"
#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/HxswgUtils.h"
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
#include "TLegendEntry.h"
#include "TPaveText.h"
#include "TObjArray.h"
#include "THStack.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "Math/QuantFuncMathCore.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"

#include<iostream>
#include<fstream>
#include<map>
#include<algorithm>
#include<vector>
#include<set>
#include <regex>


using namespace std;
double NonResonnantSyst = 0.13;
double GammaJetSyst = 0.25;
double FakeLeptonDDSyst = 0.40;

TString signalSufix="";
TString histo(""), histoVBF("");
int rebinVal = 1;
double MCRescale = 1.0;
double SignalRescale = 1.0;
int mass;
bool shape = false;
TString postfix="";
TString systpostfix="";
bool runSystematics = false; 
std::vector<TString> Channels;
std::vector<string> AnalysisBins;
double DDRescale = 1.0;
TString DYFile ="";
TString FREFile="";
string signalTag="";

bool BackExtrapol  = false;
bool subNRB        = false;
bool MCclosureTest = false;
bool scaleVBF      = false;

bool mergeWWandZZ = false;
bool skipWW = true;
bool skipGGH = false;
bool skipQQH = false;
bool subDY = false;
bool subWZ = false;
bool subFake = false;
bool blindData = false;
bool blindWithSignal = false; 
TString inFileUrl(""),jsonFile("");
double shapeMin =-9999;
double shapeMax = 9999;
double shapeMinVBF =-9999;
double shapeMaxVBF = 9999;
bool doInterf = false;
double minSignalYield = 0;
float statBinByBin = -1;

bool dirtyFix1 = false;
bool dirtyFix2 = false;

std::vector<int> shapeBinToConsider;

std::vector<int> indexcutV;
std::vector<int> indexcutVL;
std::vector<int> indexcutVR;

std::map<string, int> indexcutM;
std::map<string, int> indexcutML;
std::map<string, int> indexcutMR;

std::vector<string> keywords;


int indexvbf = -1;
int massL=-1, massR=-1;


double dropBckgBelow=0.01;

bool matchKeyword(JSONWrapper::Object& process, std::vector<string>& keywords){
  if(keywords.size()<=0)return true;
  if(process.isTag("keys")){
    std::vector<JSONWrapper::Object> dsetkeywords = process["keys"].daughters();
    for(size_t ikey=0; ikey<dsetkeywords.size(); ikey++){
      for(unsigned int i=0;i<keywords.size();i++){if(std::regex_match(dsetkeywords[ikey].toString(),std::regex(keywords[i])))return true;}
    }
  }else{
    return true;
  }
  return false;
}



void filterBinContent(TH1* histo){
  if(shapeBinToConsider.size()<=0)return;
  for(int i=0;i<=histo->GetNbinsX()+1;i++){
    bool toBeConsidered=false;  for(unsigned int j=0;j<shapeBinToConsider.size();j++){if(shapeBinToConsider[j]==i){toBeConsidered=true;break;}}
    if(!toBeConsidered){histo->SetBinContent(i,0); histo->SetBinError(i,0);}
  }
}


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

  	void clearSyst(){
     	TH1* nominal = histo();
     	uncScale.clear();
     	uncShape.clear();
     	uncShape[""] = nominal;
  	}


  	void removeStatUnc(){
     	for(auto unc = uncShape.begin(); unc!= uncShape.end(); unc++){
        TString name = unc->first.c_str();
        if(name.Contains("stat") && (name.Contains("Up") || name.Contains("Down"))){
          uncShape.erase(unc);
          unc--;
        }
     	}
  	}

  	void makeStatUnc(string prefix="", string suffix="", string suffix2="", bool noBinByBin=false){
     	if(!histo() || histo()->Integral()<=0)return;
     	string delimiter = "_";
     	unsigned firstDelimiter = suffix.find(delimiter);
     	unsigned lastDelimiter = suffix.find_last_of(delimiter);
     	unsigned endPosOfFirstDelimiter = firstDelimiter + delimiter.length();
     	string channel_and_bin = suffix.substr(endPosOfFirstDelimiter, lastDelimiter-endPosOfFirstDelimiter);
     	if(suffix.find("instrmet") != std::string::npos){
        TH1* h = (TH1*) histo()->Clone("TMPFORSTAT");
        int BIN=0;
	std::vector<unsigned int > v_lowStatBin;
	v_lowStatBin.clear();
	TString InstrMET_gammaStats_Url(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/InstrMET_systematics/InstrMET_systematics_GAMMASTATS.root");
	TFile* f_InstrMET_gammaStats = TFile::Open(InstrMET_gammaStats_Url);
	TH1* h_InstrMET_Up_gammaStats = (TH1*)utils::root::GetObjectFromPath(f_InstrMET_gammaStats, (channel_and_bin+"_mt_InstrMET_absolute_shape_up").c_str() );
	TH1* h_InstrMET_Down_gammaStats = (TH1*)utils::root::GetObjectFromPath(f_InstrMET_gammaStats, (channel_and_bin+"_mt_InstrMET_absolute_shape_down").c_str() );   			
	for(int ibin=1; ibin<=h->GetXaxis()->GetNbins(); ibin++){           
          if( true /*h->GetBinContent(ibin)/h->Integral()>0.01*/){ //This condition is removed for the moment, we may put it back in the future. 

	    char ibintxt[255]; sprintf(ibintxt, "_b%i", BIN);BIN++;
	    TH1* statU=(TH1 *)h->Clone(TString(h->GetName())+"StatU"+ibintxt);//  statU->Reset();
	    TH1* statD=(TH1 *)h->Clone(TString(h->GetName())+"StatD"+ibintxt);//  statD->Reset();           

	    statU->SetBinContent(ibin, std::max(0.0, h_InstrMET_Up_gammaStats->GetBinContent(ibin))>0 ? h_InstrMET_Up_gammaStats->GetBinContent(ibin) : 0.115);
	    statD->SetBinContent(ibin, std::max(0.0, h_InstrMET_Down_gammaStats->GetBinContent(ibin))); 
	    uncShape[prefix+"stat"+suffix+ibintxt+suffix2+"Up"  ] = statU;
	    uncShape[prefix+"stat"+suffix+ibintxt+suffix2+"Down"] = statD;

					}
          else{
            v_lowStatBin.push_back(ibin);
          }
        }

     		TH1* statU=(TH1 *)h->Clone(TString(h->GetName())+"StatU");
     		TH1* statD=(TH1 *)h->Clone(TString(h->GetName())+"StatD");

				if(v_lowStatBin.size()>0){
        	for(unsigned int j=0; j < v_lowStatBin.size(); j++){
            statU->SetBinContent(v_lowStatBin[j], std::max(0.0, h_InstrMET_Up_gammaStats->GetBinContent(v_lowStatBin[j])));   
            statD->SetBinContent(v_lowStatBin[j], std::max(0.0, h_InstrMET_Down_gammaStats->GetBinContent(v_lowStatBin[j])));   
        	}
     			uncShape[prefix+"stat"+suffix+"Up"  ] = statU;
     			uncShape[prefix+"stat"+suffix+"Down"] = statD;
        }	

				//f_InstrMET_gammaStats->Close();
     		delete h; //all done with this copy

			}
     	else{
    
    TH1* h = (TH1*) histo()->Clone("TMPFORSTAT");

     	//bin by bin stat uncertainty
     	if(statBinByBin>0 && shape==true && !noBinByBin){
        int BIN=0;
        for(int ibin=1; ibin<=h->GetXaxis()->GetNbins(); ibin++){           
          if( !(h->GetBinContent(ibin)<=0 && h->GetBinError(ibin)>0) &&  (h->GetBinContent(ibin)<=0 || h->GetBinContent(ibin)/h->Integral()<0.01 || h->GetBinError(ibin)/h->GetBinContent(ibin)<statBinByBin))continue;
					//           if(h->GetBinContent(ibin)<=0)continue;
          char ibintxt[255]; sprintf(ibintxt, "_b%i", BIN);BIN++;
          TH1* statU=(TH1 *)h->Clone(TString(h->GetName())+"StatU"+ibintxt);//  statU->Reset();
          TH1* statD=(TH1 *)h->Clone(TString(h->GetName())+"StatD"+ibintxt);//  statD->Reset();           
          if(h->GetBinContent(ibin)>0){
            statU->SetBinContent(ibin,std::min(2*h->GetBinContent(ibin), std::max(0.01*h->GetBinContent(ibin), h->GetBinContent(ibin) + h->GetBinError(ibin))));   statU->SetBinError(ibin, 0.0);
            statD->SetBinContent(ibin,std::min(2*h->GetBinContent(ibin), std::max(0.01*h->GetBinContent(ibin), h->GetBinContent(ibin) - h->GetBinError(ibin))));   statD->SetBinError(ibin, 0.0);
						//            statU->SetBinContent(ibin,std::min(2*h->GetBinContent(ibin), std::max(0.0, h->GetBinContent(ibin) + h->GetBinError(ibin))));   statU->SetBinError(ibin, 0);
						//            statD->SetBinContent(ibin,std::min(2*h->GetBinContent(ibin), std::max(0.0, h->GetBinContent(ibin) - h->GetBinError(ibin))));   statD->SetBinError(ibin, 0);
          }else{
            statU->SetBinContent(ibin,              statU->GetBinContent(ibin) + statU->GetBinError(ibin));
            statD->SetBinContent(ibin,std::max(0.0, statD->GetBinContent(ibin) - statD->GetBinError(ibin)));
          }
          uncShape[prefix+"stat"+suffix+ibintxt+suffix2+"Up"  ] = statU;
          uncShape[prefix+"stat"+suffix+ibintxt+suffix2+"Down"] = statD;
          /*h->SetBinContent(ibin, 0);*/  h->SetBinError(ibin, 0);  //remove this bin from shape variation for the other ones
          //printf("%s --> %f - %f - %f\n", (prefix+"stat"+suffix+ibintxt+suffix2+"Up").c_str(), statD->Integral(), h->GetBinContent(ibin), statU->Integral() );
        }
     	}

     	//after this line, all bins with large stat uncertainty have been considered separately
     	//so now it remains to consider all the other bins for which we assume a total correlation bin by bin
     	if(h->Integral()<=0)return; //all non empty bins have already bin variated
     	TH1* statU=(TH1 *)h->Clone(TString(h->GetName())+"StatU");
     	TH1* statD=(TH1 *)h->Clone(TString(h->GetName())+"StatD");
     	for(int ibin=1; ibin<=statU->GetXaxis()->GetNbins(); ibin++){
        if(h->GetBinContent(ibin)>0){
          statU->SetBinContent(ibin,std::min(2*h->GetBinContent(ibin), std::max(0.01*h->GetBinContent(ibin), statU->GetBinContent(ibin) + statU->GetBinError(ibin))));
          statD->SetBinContent(ibin,std::min(2*h->GetBinContent(ibin), std::max(0.01*h->GetBinContent(ibin), statD->GetBinContent(ibin) - statD->GetBinError(ibin))));
        }else{
          statU->SetBinContent(ibin,              statU->GetBinContent(ibin) + statU->GetBinError(ibin));
          statD->SetBinContent(ibin,std::min(0.0, statD->GetBinContent(ibin) - statD->GetBinError(ibin)));
        }
     	}
     	uncShape[prefix+"stat"+suffix+"Up"  ] = statU;
     	uncShape[prefix+"stat"+suffix+"Down"] = statD;

     	delete h; //all done with this copy
  		}
  		}

  	double getScaleUncertainty(){
     	double Total=0;
     	for(std::map<string, double>::iterator unc=uncScale.begin();unc!=uncScale.end();unc++){
        if(unc->second<0)continue;
        Total+=pow(unc->second,2);
     	}     
     	return Total>0?sqrt(Total):-1;
  	}

		double getIntegratedShapeUncertainty(string name, string upORdown){
     	double Total=0;
     	//this = ch->second.shapes[histoName.Data()]
     	for(std::map<string, TH1*>::iterator var = uncShape.begin(); var!=uncShape.end(); var++){
       	TString systName = var->first.c_str();
       	if(var->first=="")continue; //Skip Nominal shape
       	if(!systName.Contains(upORdown))continue; //only look for syst up or down at a time (upORdown should be either "Up" or "Down"
       	TH1* hvar = (TH1*)(var->second->Clone((name+var->first).c_str()));

 			 	double varYield = hvar->Integral();
			 	TH1* h = (TH1*)(this->histo()->Clone((name+"Nominal").c_str()));
		   	double yield = h->Integral();
       	Total+=pow(varYield-yield,2); //the total shape unc is the sqrt of the quadratical sum of the difference between the nominal and the variated yields.
     	}     
     	return Total>0?sqrt(Total):-1;
  	}

		double getBinShapeUncertainty(string name, int bin, string upORdown){
     	double Total=0;
     	//this = ch->second.shapes[histoName.Data()]
     	for(std::map<string, TH1*>::iterator var = uncShape.begin(); var!=uncShape.end(); var++){
       	TString systName = var->first.c_str();
       	//if(var->first==""){
       	//}
       	if(var->first=="")continue; //Skip Nominal shape
       	if(!systName.Contains(upORdown))continue; //only look for syst up or down at a time (upORdown should be either "Up" or "Down"
       	TH1* hvar = (TH1*)(var->second->Clone((name+var->first).c_str()));

 			 	double varYield = hvar->GetBinContent(bin);
			 	TH1* h = (TH1*)(this->histo()->Clone((name+"Nominal").c_str()));
		   	double yield = h->GetBinContent(bin);
       	Total+=pow(varYield-yield,2); //the total shape unc is the sqrt of the quadratical sum of the difference between the nominal and the variated yields.
     	}     
     	return Total>0?sqrt(Total):-1;
  	}


  	void rescaleScaleUncertainties(double StartIntegral, double EndIntegral){
     	for(std::map<string, double>::iterator unc=uncScale.begin();unc!=uncScale.end();unc++){
        printf("%E/%E = %E = %E/%E\n", unc->second, StartIntegral, unc->second/StartIntegral, EndIntegral * (unc->second/StartIntegral), EndIntegral); 
        if(StartIntegral!=0){unc->second = EndIntegral * (unc->second/StartIntegral);}else{unc->second = unc->second * EndIntegral;}
     	}     
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
    void addChannel(ChannelInfo_t& dest, ChannelInfo_t& src, bool computeSyst = false);

    // Sum up all background processes and add this as a total process
    void addProc(ProcessInfo_t& dest, ProcessInfo_t& src, bool computeSyst = false);

    // Sum up all background processes and add this as a total process
    void computeTotalBackground();

    // Replace the Data process by TotalBackground
    void blind();

    // Print the Yield table
    void getYieldsFromShape(FILE* pFile, std::vector<TString>& selCh, string histoName, FILE* pFileInc=NULL);

    // Dump efficiencies
    void getEffFromShape(FILE* pFile, std::vector<TString>& selCh, string histoName);

    // drop background process that have a negligible yield
    void dropSmallBckgProc(std::vector<TString>& selCh, string histoName, double threshold);

    // drop control channels
    void dropCtrlChannels(std::vector<TString>& selCh);

    // Make a summary plot
    void showShape(std::vector<TString>& selCh , TString histoName, TString SaveName);

    // Make a summary plot of the uncertainties
    void showUncertainty(std::vector<TString>& selCh , TString histoName, TString SaveName);

    // Turn to cut&count (rebin all histo to 1 bin only)
    void turnToCC(string histoName);

    // Make a summary plot
    void saveHistoForLimit(string histoName, TFile* fout);

    // Add hardcoded uncertainties 
    void addHardCodedUncertainties(string histoName);

    // produce the datacards 
    void buildDataCards(string histoName, TString url);

    // Load histograms from root file and json to memory
    void getShapeFromFile(TFile* inF, std::vector<string> channelsAndShapes, int cutBin, JSONWrapper::Object &Root,  double minCut=0, double maxCut=9999, bool onlyData=false);

    // replace MC NonResonnant Backgrounds by DataDriven estimate
    void doBackgroundSubtraction(FILE* pFile, std::vector<TString>& selCh,TString ctrlCh,TString mainHisto, TString sideBandHisto);

    //add syst uncertainty to the instr. met background if it already exist in the file
    void addInstrMetSyst(FILE* pFile, std::vector<TString>& selCh,TString ctrlCh,TString mainHisto);
    void addInstrMetSyst_2017(std::vector<TString>& selCh,TString mainHisto); //with the method used in 2017

    // replace MC Z+Jets Backgrounds by DataDriven Gamma+Jets estimate
    void doDYReplacement(FILE* pFile, std::vector<TString>& selCh,TString ctrlCh,TString mainHisto);

    // replace MC Backgrounds with FakeLeptons by DataDriven estimate
    void doFakeLeptonEstimation(FILE* pFile, std::vector<TString>& selCh,TString ctrlCh,TString mainHisto, bool isCutAndCount);

    // Rebin histograms to make sure that high mt/met region have no empty bins
    void rebinMainHisto(string histoName);

    // Interpollate the signal sample between two mass points 
    void SignalInterpolation(string histoName);

    // scale VBF production to VBF/GGF SM prediction;
    void scaleVBF(string histoName);          

    // Rescale signal sample for the effect of the interference and propagate the uncertainty 
    void RescaleForInterference(string histoName);

    //Merge bins together
    void mergeBins(std::vector<string>& binsToMerge, string NewName);

    // Handle empty bins
    void HandleEmptyBins(string histoName);

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
  printf("--index     --> index of selection to be used (Xbin in histogram to be used); different comma separated values can be given for each analysis bin\n");
  printf("--indexL    --> index of selection to be used (Xbin in histogram to be used) used for interpolation;  different comma separated values can be given for each analysis bin\n");
  printf("--indexR    --> index of selection to be used (Xbin in histogram to be used) used for interpolation;  different comma separated values can be given for each analysis bin\n");
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
  printf("--postfix    --> use this to specify a postfix that will be added to the process names)\n");
  printf("--systpostfix    --> use this to specify a syst postfix that will be added to the process names)\n");
  printf("--MCRescale    --> use this to rescale the cross-section of all MC processes by a given factor)\n");
  printf("--signalRescale    --> use this to rescale signal cross-section by a given factor)\n");
  printf("--interf     --> use this to rescale xsection according to WW interferences)\n");
  printf("--minSignalYield   --> use this to specify the minimum Signal yield you want in each channel)\n");
  printf("--signalSufix --> use this flag to specify a suffix string that should be added to the signal 'histo' histogram\n");
  printf("--signalTag   --> use this flag to specify a tag that should be present in signal sample name\n");
  printf("--rebin         --> rebin the histogram\n");
  printf("--statBinByBin --> make bin by bin statistical uncertainty\n");
  printf("--inclusive  --> merge bins to make the analysis inclusive\n");
  printf("--dropBckgBelow --> drop all background processes that contributes for less than a threshold to the total background yields\n");
  printf("--scaleVBF    --> scale VBF signal by ggH/VBF\n");
  printf("--key        --> provide a key for sample filtering in the json\n");  
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

  //get input arguments
  for(int i=1;i<argc;i++){
    string arg(argv[i]);
    if(arg.find("--help")          !=string::npos) { printHelp(); return -1;} 
    else if(arg.find("--minSignalYield") !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%lf",&minSignalYield ); i++; printf("minSignalYield = %f\n", minSignalYield);}
    else if(arg.find("--scaleVBF") !=string::npos) { scaleVBF=true; printf("scaleVBF = True\n");}
    else if(arg.find("--subNRB")   !=string::npos) { subNRB=true; skipWW=true; printf("subNRB = True\n");}
    else if(arg.find("--subDY")    !=string::npos) { subDY=true; DYFile=argv[i+1];  i++; printf("Z+Jets will be replaced by %s\n",DYFile.Data());}
    else if(arg.find("--subFake")  !=string::npos) { subFake=true; FREFile=argv[i+1];  i++; printf("Fake lepton procs will be replaced by %s\n",FREFile.Data());}
    else if(arg.find("--subWZ")    !=string::npos) { subWZ=true; printf("WZ will be estimated from 3rd lepton SB\n");}
    else if(arg.find("--DDRescale")!=string::npos && i+1<argc)  { sscanf(argv[i+1],"%lf",&DDRescale); i++;}
    else if(arg.find("--MCRescale")!=string::npos && i+1<argc)  { sscanf(argv[i+1],"%lf",&MCRescale); i++;}
    else if(arg.find("--signalRescale")!=string::npos && i+1<argc)  { sscanf(argv[i+1],"%lf",&SignalRescale); i++;}
    else if(arg.find("--HWW")      !=string::npos) { skipWW=false; printf("HWW = True\n");}
    else if(arg.find("--skipGGH")  !=string::npos) { skipGGH=true; printf("skipGGH = True\n");}
    else if(arg.find("--skipQQH")  !=string::npos) { skipQQH=true; printf("skipQQH = True\n");}
    else if(arg.find("--blindWithSignal")    !=string::npos) { blindData=true; blindWithSignal=true; printf("blindData = True; blindWithSignal = True\n");}
    else if(arg.find("--blind")    !=string::npos) { blindData=true; printf("blindData = True\n");}
    else if(arg.find("--closure")  !=string::npos) { MCclosureTest=true; printf("MCclosureTest = True\n");}
    else if(arg.find("--shapeBinToConsider")    !=string::npos && i+1<argc)  { char* pch = strtok(argv[i+1],",");while (pch!=NULL){int C;  sscanf(pch,"%i",&C); shapeBinToConsider.push_back(C);  pch = strtok(NULL,",");} i++; printf("Only the following histo bins will be considered: "); for(unsigned int i=0;i<shapeBinToConsider.size();i++)printf(" %i ", shapeBinToConsider[i]);printf("\n");}
    else if(arg.find("--shapeMinVBF") !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%lf",&shapeMinVBF); i++; printf("Min cut on shape for VBF = %f\n", shapeMinVBF);}
    else if(arg.find("--shapeMaxVBF") !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%lf",&shapeMaxVBF); i++; printf("Max cut on shape for VBF = %f\n", shapeMaxVBF);}
    else if(arg.find("--shapeMin") !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%lf",&shapeMin); i++; printf("Min cut on shape = %f\n", shapeMin);}
    else if(arg.find("--shapeMax") !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%lf",&shapeMax); i++; printf("Max cut on shape = %f\n", shapeMax);}
    else if(arg.find("--interf")    !=string::npos) { doInterf=true; printf("doInterf = True\n");}
    else if(arg.find("--indexvbf") !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%i",&indexvbf); i++; printf("indexVBF = %i\n", indexvbf);}
    else if(arg.find("--index" )   !=string::npos && i+1<argc)   { char* pch = strtok(argv[i+1],",");while (pch!=NULL){int C;  sscanf(pch,"%i",&C); indexcutV .push_back(C);  pch = strtok(NULL,",");} i++; printf("index  = "); for(unsigned int i=0;i<indexcutV .size();i++)printf(" %i ", indexcutV [i]);printf("\n");}
    else if(arg.find("--indexL")    !=string::npos && i+1<argc)  { char* pch = strtok(argv[i+1],",");while (pch!=NULL){int C;  sscanf(pch,"%i",&C); indexcutVL.push_back(C);  pch = strtok(NULL,",");} i++; printf("indexL = "); for(unsigned int i=0;i<indexcutVL.size();i++)printf(" %i ", indexcutVL[i]);printf("\n");}
    else if(arg.find("--indexR")    !=string::npos && i+1<argc)  { char* pch = strtok(argv[i+1],",");while (pch!=NULL){int C;  sscanf(pch,"%i",&C); indexcutVR.push_back(C);  pch = strtok(NULL,",");} i++; printf("indexR = "); for(unsigned int i=0;i<indexcutVR.size();i++)printf(" %i ", indexcutVR[i]);printf("\n");}
    else if(arg.find("--in")       !=string::npos && i+1<argc)  { inFileUrl = argv[i+1];  i++;  printf("in = %s\n", inFileUrl.Data());  }
    else if(arg.find("--json")     !=string::npos && i+1<argc)  { jsonFile  = argv[i+1];  i++;  printf("json = %s\n", jsonFile.Data()); }
    else if(arg.find("--histoVBF") !=string::npos && i+1<argc)  { histoVBF  = argv[i+1];  i++;  printf("histoVBF = %s\n", histoVBF.Data()); }
    else if(arg.find("--histo")    !=string::npos && i+1<argc)  { histo     = argv[i+1];  i++;  printf("histo = %s\n", histo.Data()); }
    else if(arg.find("--mL")       !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%i",&massL ); i++; printf("massL = %i\n", massL);}
    else if(arg.find("--mR")       !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%i",&massR ); i++; printf("massR = %i\n", massR);}
    else if(arg.find("--m")        !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%i",&mass ); i++; printf("mass = %i\n", mass);}
    else if(arg.find("--bins")     !=string::npos && i+1<argc)  { char* pch = strtok(argv[i+1],",");printf("bins are : ");while (pch!=NULL){printf(" %s ",pch); AnalysisBins.push_back(pch);  pch = strtok(NULL,",");}printf("\n"); i++; }
    else if(arg.find("--channels") !=string::npos && i+1<argc)  { char* pch = strtok(argv[i+1],",");printf("channels are : ");while (pch!=NULL){printf(" %s ",pch); Channels.push_back(pch);  pch = strtok(NULL,",");}printf("\n"); i++; }
    else if(arg.find("--postfix")   !=string::npos && i+1<argc)  { postfix = argv[i+1]; systpostfix = argv[i+1]; i++;  printf("postfix '%s' will be used\n", postfix.Data());  }
    else if(arg.find("--systpostfix")   !=string::npos && i+1<argc)  { systpostfix = argv[i+1];  i++;  printf("systpostfix '%s' will be used\n", systpostfix.Data());  }
    else if(arg.find("--syst")     !=string::npos) { runSystematics=true; printf("syst = True\n");}
    else if(arg.find("--shape")    !=string::npos) { shape=true; printf("shapeBased = True\n");}
    else if(arg.find("--dirtyFix2")    !=string::npos) { dirtyFix2=true; printf("dirtyFix2 = True\n");}
    else if(arg.find("--dirtyFix1")    !=string::npos) { dirtyFix1=true; printf("dirtyFix1 = True\n");}
    else if(arg.find("--signalSufix") !=string::npos) { signalSufix = argv[i+1]; i++; printf("signalSufix '%s' will be used\n", signalSufix.Data()); }
    else if(arg.find("--signalTag") !=string::npos) { signalTag = argv[i+1]; i++; printf("signalTag '%s' will be used\n", signalTag.c_str()); }
    else if(arg.find("--rebin")    !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%i",&rebinVal); i++; printf("rebin = %i\n", rebinVal);}
    else if(arg.find("--BackExtrapol")    !=string::npos) { BackExtrapol=true; printf("BackExtrapol = True\n");}
    else if(arg.find("--statBinByBin")    !=string::npos) { sscanf(argv[i+1],"%f",&statBinByBin); i++; printf("statBinByBin = %f\n", statBinByBin);}
    else if(arg.find("--dropBckgBelow")   !=string::npos) { sscanf(argv[i+1],"%lf",&dropBckgBelow); i++; printf("dropBckgBelow = %f\n", dropBckgBelow);}
    else if(arg.find("--key"          )   !=string::npos && i+1<argc){ keywords.push_back(argv[i+1]); printf("Only samples matching this (regex) expression '%s' are processed\n", argv[i+1]); i++;  }

  }
  if(jsonFile.IsNull()) { printf("No Json file provided\nrun with '--help' for more details\n"); return -1; }
  if(inFileUrl.IsNull()){ printf("No Inputfile provided\nrun with '--help' for more details\n"); return -1; }
  if(histo.IsNull())    { printf("No Histogram provided\nrun with '--help' for more details\n"); return -1; }
  if(mass==-1)          { printf("No massPoint provided\nrun with '--help' for more details\n"); return -1; }
  if(indexcutV.size()<=0){printf("INDEX CUT SIZE IS NULL\n"); printHelp(); return -1; }
  if(AnalysisBins.size()==0)AnalysisBins.push_back("all");
  if(Channels.size()==0){Channels.push_back("ee");Channels.push_back("mumu");}

  //make sure that the index vector are well filled
  if(indexcutVL.size()==0) indexcutVL.push_back(indexcutV [0]);
  if(indexcutVR.size()==0) indexcutVR.push_back(indexcutV [0]);
  while(indexcutV .size()<AnalysisBins.size()){indexcutV .push_back(indexcutV [0]);}
  while(indexcutVL.size()<AnalysisBins.size()){indexcutVL.push_back(indexcutVL[0]);}
  while(indexcutVR.size()<AnalysisBins.size()){indexcutVR.push_back(indexcutVR[0]);}
  if(indexvbf>=0){for(unsigned int i=0;i<AnalysisBins.size();i++){if(AnalysisBins[i].find("vbf")!=string::npos){indexcutV[i]=indexvbf; indexcutVL[i]=indexvbf; indexcutVR[i]=indexvbf;} }}



  //handle merged bins
  std::vector<std::vector<string> > binsToMerge;
  for(unsigned int b=0;b<AnalysisBins.size();b++){
    if(AnalysisBins[b].find('+')!=std::string::npos){
      std::vector<string> subBins;
      char* pch = strtok(&AnalysisBins[b][0],"+"); 
      while (pch!=NULL){
        indexcutV.push_back(indexcutV[b]);
        indexcutVL.push_back(indexcutVL[b]);
        indexcutVR.push_back(indexcutVR[b]);
        AnalysisBins.push_back(pch);
        subBins.push_back(pch);
        pch = strtok(NULL,"+");
      }
      binsToMerge.push_back(subBins);
      AnalysisBins.erase(AnalysisBins.begin()+b);
      indexcutV .erase(indexcutV .begin()+b);
      indexcutVL.erase(indexcutVL.begin()+b);
      indexcutVR.erase(indexcutVR.begin()+b);
      b--;
    }
  }


  //fill the index map
  for(unsigned int i=0;i<AnalysisBins.size();i++){indexcutM[AnalysisBins[i]] = indexcutV[i]; indexcutML[AnalysisBins[i]] = indexcutVL[i]; indexcutMR[AnalysisBins[i]] = indexcutVR[i];}


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

  AllInfo_t allInfo;



  //open input file
  TFile* inF = TFile::Open(inFileUrl);
  if( !inF || inF->IsZombie() ){ printf("Invalid file name : %s\n", inFileUrl.Data());}
  gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE


  //LOAD shapes
  const size_t nsh=sh.size();
  for(size_t b=0; b<AnalysisBins.size(); b++){
    std::vector<string> channelsAndShapes;
    for(size_t i=0; i<nch; i++){
      for(size_t j=0; j<nsh; j++){
        channelsAndShapes.push_back((ch[i]+TString(";")+AnalysisBins[b]+TString(";")+sh[j]).Data());
      }
    }
    double cutMin=shapeMin; double cutMax=shapeMax;
    if((shapeMinVBF!=shapeMin || shapeMaxVBF!=shapeMax) && AnalysisBins[b].find("vbf")!=string::npos){cutMin=shapeMinVBF; cutMax=shapeMaxVBF;}
    allInfo.getShapeFromFile(inF,    channelsAndShapes, indexcutM[AnalysisBins[b]], Root, cutMin, cutMax   );     
  }


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
  if(subDY){//Hugo: This is something completely outdated... the way to choose for datadriven DY or mc-based, is just by giving the good key
  	pFile = fopen("GammaJets.tex","w");
  	allInfo.doDYReplacement(pFile, selCh,"gamma",histo);
  	fclose(pFile);
  }


  //replace Z+Jet background by Gamma+Jet estimates
  if(subFake){
  	pFile = fopen("FakeRateEstimate.tex","w");
  	allInfo.doFakeLeptonEstimation(pFile, selCh,"gamma",histo, !shape);
  	fclose(pFile);
  }

  //replace data by total MC background
  if(blindData)allInfo.blind();

  //interpollate signal sample if desired mass point is not available
  allInfo.SignalInterpolation(histo.Data());

  if(scaleVBF)  allInfo.scaleVBF(histo.Data());

  //rescale for interference
  //if(doInterf || signalSufix!="")allInfo.RescaleForInterference(histo.Data());  //disabled as the logic is moved to runXXXAnalysis code

  //extrapolate backgrounds toward higher mt/met region to make sure that there is no empty bins
  if(shape && BackExtrapol)allInfo.rebinMainHisto(histo.Data());

  // add the syst. on the Instr MET 
  allInfo.addInstrMetSyst_2017(selCh,histo);

  //drop backgrounds with rate<1%
  allInfo.dropSmallBckgProc(selCh, histo.Data(), dropBckgBelow);

  //drop control channels
  allInfo.dropCtrlChannels(selCh);

  //merge bins  
  for(unsigned int B=0;B<binsToMerge.size();B++){
    std::string NewBinName = string("["); binsToMerge[B][0];  for(unsigned int b=1;b<binsToMerge[B].size();b++){NewBinName += "+"+binsToMerge[B][b];} NewBinName+="]";
    allInfo.mergeBins(binsToMerge[B],NewBinName);
  }


  //turn to CC analysis eventually
  if(!shape)allInfo.turnToCC(histo.Data());

  allInfo.HandleEmptyBins(histo.Data()); //needed for negative bin content --> May happens due to NLO interference for instance

  if(blindData)allInfo.blind();

  //print event yields from the mt shapes
  pFile = fopen("Yields.tex","w");  FILE* pFileInc = fopen("YieldsInc.tex","w");
  allInfo.getYieldsFromShape(pFile, selCh, histo.Data(), pFileInc);
  fclose(pFile); fclose(pFileInc);

  //print signal efficiency
  pFile = fopen("Efficiency.tex","w");
  allInfo.getEffFromShape(pFile, selCh, histo.Data());
  fclose(pFile);

  //add by hand the hard coded uncertainties
  allInfo.addHardCodedUncertainties(histo.Data());

  //produce a plot
  allInfo.showShape(selCh,histo,"plot"); //this produce the final global shape

  //produce a plot
  allInfo.showUncertainty(selCh,histo,"plot"); //this produces all the plots with the syst

  //prepare the output
  string limitFile=("hzz2l2v_"+massStr+systpostfix+".root").Data();
  TFile *fout=TFile::Open(limitFile.c_str(),"recreate");

  allInfo.saveHistoForLimit(histo.Data(), fout);

  allInfo.buildDataCards(histo.Data(), limitFile);

  //all done
  fout->Close();
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
  if(isData)sorted_procs.push_back("data");
  sorted_procs.insert(sorted_procs.end(), sign_procs.begin(), sign_procs.end());
}

//
// Sum up all shapes from one src channel to a total shapes in the dest channel
//
void AllInfo_t::addChannel(ChannelInfo_t& dest, ChannelInfo_t& src, bool computeSyst){
  std::map<string, ShapeData_t>& shapesInfoDest = dest.shapes;
  std::map<string, ShapeData_t>& shapesInfoSrc  = src.shapes;

	if(!computeSyst){
	for(std::map<string, ShapeData_t>::iterator sh = shapesInfoSrc.begin(); sh!=shapesInfoSrc.end(); sh++){
    if(shapesInfoDest.find(sh->first)==shapesInfoDest.end())shapesInfoDest[sh->first] = ShapeData_t();

    //Loop on all shape systematics (including also the central value shape)
    for(std::map<string, TH1*>::iterator uncS = sh->second.uncShape.begin();uncS!= sh->second.uncShape.end();uncS++){
      if(uncS->first!="") continue; //We only take nominal shapes
      if(shapesInfoDest[sh->first].uncShape.find(uncS->first)==shapesInfoDest[sh->first].uncShape.end()){
        shapesInfoDest[sh->first].uncShape[uncS->first] = (TH1*) uncS->second->Clone(TString(uncS->second->GetName() + dest.channel + dest.bin ) );
      }else{
        shapesInfoDest[sh->first].uncShape[uncS->first]->Add(uncS->second);
      }
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
  else{



 for(std::map<string, ShapeData_t>::iterator sh = shapesInfoSrc.begin(); sh!=shapesInfoSrc.end(); sh++){
    if(shapesInfoDest.find(sh->first)==shapesInfoDest.end())shapesInfoDest[sh->first] = ShapeData_t();
 		//Loop on all shape systematics (including also the central value shape)
    for(std::map<string, TH1*>::iterator uncS = sh->second.uncShape.begin();uncS!= sh->second.uncShape.end();uncS++){
      if(uncS->first=="") continue; //We only take systematic (i.e non-nominal) shapes
      if(shapesInfoSrc[sh->first].uncShape.find("")==shapesInfoSrc[sh->first].uncShape.end()) continue;
      //1. Copy the nominal shape
      shapesInfoDest[sh->first].uncShape[uncS->first] = (TH1*) shapesInfoDest[sh->first].uncShape[""]->Clone(TString(uncS->second->GetName() + dest.channel + dest.bin ) );
      //2. we remove the nominal value of the process we are running on
      shapesInfoDest[sh->first].uncShape[uncS->first]->Add(shapesInfoSrc[sh->first].uncShape[""], -1);
      //3. and add the variation up/down
			shapesInfoDest[sh->first].uncShape[uncS->first]->Add(uncS->second); 
		}
  }



  }


}

//
// Sum up all background processes and add this as a total process
//
void AllInfo_t::addProc(ProcessInfo_t& dest, ProcessInfo_t& src, bool computeSyst){
  dest.xsec = src.xsec*src.br;
  for(std::map<string, ChannelInfo_t>::iterator ch = src.channels.begin(); ch!=src.channels.end(); ch++){
    if(dest.channels.find(ch->first)==dest.channels.end()){   //this channel does not exist, create it
      dest.channels[ch->first]         = ChannelInfo_t();
      dest.channels[ch->first].bin     = ch->second.bin;
      dest.channels[ch->first].channel = ch->second.channel;
    }

    addChannel(dest.channels[ch->first], ch->second, computeSyst);
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
  procInfo_Bckgs.shortName = "total";
  procInfo_Bckgs.isData = false;
  procInfo_Bckgs.isSign = false;
  procInfo_Bckgs.isBckg = true;
  procInfo_Bckgs.xsec   = 0.0;
  procInfo_Bckgs.br     = 1.0;
  //Compute total background nominal
  for(std::map<string, ProcessInfo_t>::iterator it=procs.begin(); it!=procs.end();it++){
    if(it->first=="total" || it->second.isBckg!=true)continue;
    addProc(procInfo_Bckgs, it->second, false);
  }
  //Compute total background systematics
  for(std::map<string, ProcessInfo_t>::iterator it=procs.begin(); it!=procs.end();it++){
    if(it->first=="total" || it->second.isBckg!=true)continue;
    addProc(procInfo_Bckgs, it->second, true);
  }

}




//
// Replace the Data process by TotalBackground
//
void AllInfo_t::blind(){
  if(procs.find("total")==procs.end())computeTotalBackground();


  if(true){ //always replace data
    //if(procs.find("data")==procs.end()){ //true only if there is no "data" samples in the json file
    sorted_procs.push_back("data");           
    procs["data"] = ProcessInfo_t(); //reset
    ProcessInfo_t& procInfo_Data = procs["data"];
    procInfo_Data.shortName = "data";
    procInfo_Data.isData = true;
    procInfo_Data.isSign = false;
    procInfo_Data.isBckg = false;
    procInfo_Data.xsec   = 0.0;
    procInfo_Data.br     = 1.0;
    for(std::map<string, ProcessInfo_t>::iterator it=procs.begin(); it!=procs.end();it++){
      if(it->first!="total")continue;
      addProc(procInfo_Data, it->second);
    }
  }
  }

  //
  // Print the Yield table
  //
	/*
     void AllInfo_t::getYieldsFromShape(FILE* pFile, std::vector<TString>& selCh, string histoName, FILE* pFileInc){
     if(!pFileInc)pFileInc=pFile;

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
	if(val<1E-5 && valerr>=10*val){val=0.0; syst=-1;}
	else if(val<1E-6){val=0.0; valerr=0.0; syst=-1;}
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
	if(it->first=="data"){char tmp[256];sprintf(tmp, "%.0f", bin_val[bin->first]); rowsBin[bin->first] += tmp;  //unblinded
	//                 if(it->first=="data"){char tmp[256];sprintf(tmp, "-"); rowsBin[bin->first] += tmp;  //blinded
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

	//All Bins
	fprintf(pFileInc,"\\begin{sidewaystable}[htp]\n\\begin{center}\n\\caption{Event yields expected for background and signal processes and observed in data.}\n\\label{tab:table}\n");
	fprintf(pFileInc, "%s}\\\\\n", rows_header.c_str());
	fprintf(pFileInc, "%s\\\\\n", rows_title .c_str());
	for(std::map<string, string>::iterator row = rowsBin.begin(); row!= rowsBin.end(); row++){
	fprintf(pFileInc, "%s\\\\\n", row->second.c_str());
	}
	fprintf(pFileInc,"\\hline\n");
	fprintf(pFileInc,"\\end{tabular}\n\\end{center}\n\\end{sidewaystable}\n");
}
*/
void AllInfo_t::getYieldsFromShape(FILE* pFile, std::vector<TString>& selCh, string histoName, FILE* pFileInc){
  if(!pFileInc)pFileInc=pFile;

  std::vector<string> VectorProc;
  std::map<string, bool> MapChannel;
  std::map<string, std::map<string, string> > MapProcChYields;         
  std::map<string, bool> MapChannelBin;
  std::map<string, std::map<string, string> > MapProcChYieldsBin;         

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
    std::map<string, double> bin_systUp;
    std::map<string, double> bin_systDown;

    VectorProc.push_back(it->first);
    for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
      if(std::find(selCh.begin(), selCh.end(), ch->second.channel)==selCh.end())continue;
      if(ch->second.shapes.find(histoName)==(ch->second.shapes).end())continue;
      TH1* h = ch->second.shapes[histoName].histo();
      double valerr;
      double val  = h->IntegralAndError(1,h->GetXaxis()->GetNbins(),valerr);
      if(procName.find("Instr. MET")!=std::string::npos) valerr =0.0; //Our systematics on the Instr. MET already includes stat unc. We would want to change that in the future
      double syst_scale = std::max(0.0, ch->second.shapes[histoName].getScaleUncertainty());
      double syst_shapeUp = std::max(0.0, ch->second.shapes[histoName].getIntegratedShapeUncertainty((it->first+ch->first).c_str(), "Up"));
      double syst_shapeDown = std::max(0.0, ch->second.shapes[histoName].getIntegratedShapeUncertainty((it->first+ch->first).c_str(), "Down"));
      double systUp = sqrt(pow(syst_scale,2)+pow(syst_shapeUp,2));
      double systDown = sqrt(pow(syst_scale,2)+pow(syst_shapeDown,2));
      systUp= (systUp >0)?systUp: -1; //Set to -1 if no syst, to be coherent with other convention in this file
      systDown= (systDown >0)?systDown: -1; //Set to -1 if no syst, to be coherent with other convention in this file
      if(val<1E-5 && valerr>=10*val && procName.find("ww")!=std::string::npos){val=0.0;}
      else if(val<1E-5 && valerr>=10*val){val=0.0; systUp=-1;systDown=-1;}
      else if(val<1E-6){val=0.0; valerr=0.0; systUp=-1;systDown=-1;}
      if(it->first=="data"){valerr=-1.0; systUp=-1;systDown=-1;}
      string YieldText = "";

      if(it->first=="data" || it->first=="total")YieldText += "\\boldmath ";
      if(it->first=="data"){char tmp[256];sprintf(tmp, "$%.0f$", val); YieldText += tmp;
      }else{                YieldText += utils::toLatexRounded(val,valerr, systUp, true, systDown);     }


      printf("%f %f %f %f --> %s\n", val, valerr, systUp, systDown, utils::toLatexRounded(val,valerr, systUp, true, systDown).c_str());

      if(rows.find(ch->first)==rows.end())rows[ch->first] = string("$ ")+ch->first+" $";
      rows[ch->first] += string("&") + YieldText;

      TString LabelText = TString("$") + ch->second.channel+ " " +ch->second.bin + TString("$");
      LabelText.ReplaceAll("eq"," ="); LabelText.ReplaceAll("g =","\\geq"); LabelText.ReplaceAll("l =","\\leq"); 
      LabelText.ReplaceAll("_OS","OS "); LabelText.ReplaceAll("el","e"); LabelText.ReplaceAll("mu","\\mu");  LabelText.ReplaceAll("ha","\\tau_{had}");

      TString BinText = TString("$") + ch->second.bin + TString("$");
      BinText.ReplaceAll("eq"," ="); BinText.ReplaceAll("g =","\\geq"); BinText.ReplaceAll("l =","\\leq");
      BinText.ReplaceAll("_OS","OS "); BinText.ReplaceAll("el","e"); BinText.ReplaceAll("mu","\\mu");  BinText.ReplaceAll("ha","\\tau_{had}");


      bin_val   [BinText.Data()] = val;
      bin_valerr[BinText.Data()] = pow(valerr,2);
      bin_systUp  [BinText.Data()] = systUp>=0?pow(systUp,2):-1;
      bin_systDown  [BinText.Data()] = systDown>=0?pow(systDown,2):-1;
      if(systUp<0)bin_systUp  [BinText.Data()]=-1;
      if(systDown<0)bin_systDown  [BinText.Data()]=-1;

      bin_val   [" Inc."] += val;
      bin_valerr[" Inc."] += pow(valerr,2);
      bin_systUp  [" Inc."] += systUp>=0?pow(systUp,2):0; //We are doing a quadratic sum here, have to add 0 if we have negative value
      bin_systDown  [" Inc."] += systDown>=0?pow(systDown,2):0; //We are doing a quadratic sum here, have to add 0 if we have negative value

      MapChannel[LabelText.Data()] = true;
      MapProcChYields[it->first][LabelText.Data()] = YieldText;
    }
    if(bin_systUp  [" Inc."] <= 0) bin_systUp  ["Inc"]=-1; //If negative value, or 0, set it to -1
    if(bin_systDown  [" Inc."] <= 0) bin_systDown  ["Inc"]=-1; //If negative value, or 0, set it to -1

    for(std::map<string, double>::iterator bin=bin_val.begin(); bin!=bin_val.end(); bin++){
      string YieldText = "";                 
      if(it->first=="data" || it->first=="total" || bin->first==" Inc.")YieldText += "\\boldmath ";
      if(it->first=="data"){char tmp[256];sprintf(tmp, "%.0f", bin_val[bin->first]); YieldText += tmp;  //unblinded
				//                 if(it->first=="data"){char tmp[256];sprintf(tmp, "-"); rowsBin[bin->first] += tmp;  //blinded
      }else{                YieldText += utils::toLatexRounded(bin_val[bin->first],sqrt(bin_valerr[bin->first]), bin_systUp[bin->first]<0?-1:sqrt(bin_systUp[bin->first]), true, bin_systDown[bin->first]<0?-1:sqrt(bin_systDown[bin->first]));   }

      if(rowsBin.find(bin->first)==rowsBin.end())rowsBin[bin->first] = string("$ ")+bin->first+" $";
      rowsBin[bin->first] += string("&") + YieldText;

      MapChannelBin[bin->first] = true;
      MapProcChYieldsBin[it->first][bin->first] = YieldText;                

      if(bin->first==" Inc."){
        MapChannel[bin->first] = true;
        MapProcChYields[it->first][bin->first] = YieldText;                
      }
      }
    }


		//           //All Channels
		//           fprintf(pFile,"\\begin{sidewaystable}[htp]\n\\begin{center}\n\\caption{Event yields expected for background and signal processes and observed in data.}\n\\label{tab:table}\n");
		//           fprintf(pFile, "%s}\\\\\n", rows_header.c_str());
		//           fprintf(pFile, "%s\\\\\n", rows_title .c_str());
		//           for(std::map<string, string>::iterator row = rows.begin(); row!= rows.end(); row++){
		//              fprintf(pFile, "%s\\\\\n", row->second.c_str());
		//           }
		//           fprintf(pFile,"\\hline\n");
		//           fprintf(pFile,"\\end{tabular}\n\\end{center}\n\\end{sidewaystable}\n");

		//           //All Bins
		//           fprintf(pFileInc,"\\begin{sidewaystable}[htp]\n\\begin{center}\n\\caption{Event yields expected for background and signal processes and observed in data.}\n\\label{tab:table}\n");
		//           fprintf(pFileInc, "%s}\\\\\n", rows_header.c_str());
		//           fprintf(pFileInc, "%s\\\\\n", rows_title .c_str());
		//           for(std::map<string, string>::iterator row = rowsBin.begin(); row!= rowsBin.end(); row++){
		//              fprintf(pFileInc, "%s\\\\\n", row->second.c_str());
		//           }
		//           fprintf(pFileInc,"\\hline\n");
		//           fprintf(pFileInc,"\\end{tabular}\n\\end{center}\n\\end{sidewaystable}\n");



    //All Channels
    fprintf(pFile,"\\begin{sidewaystable}[htp]\n\\begin{center}\n\\caption{Event yields expected for background and signal processes and observed in data.}\n\\label{tab:table}\n");
    fprintf(pFile, "\\begin{tabular}{|c|"); for(auto ch = MapChannel.begin(); ch!=MapChannel.end();ch++){ fprintf(pFile, "c|"); } fprintf(pFile, "}\\\\\n");
    fprintf(pFile, "channel");   for(auto ch = MapChannel.begin(); ch!=MapChannel.end();ch++){ fprintf(pFile, " & %s", ch->first.c_str()); } fprintf(pFile, "\\\\\\hline\n");
    for(auto proc = VectorProc.begin();proc!=VectorProc.end(); proc++){
      if(*proc=="total")fprintf(pFile, "\\hline\n");
      auto ChannelYields = MapProcChYields.find(*proc);
      if(ChannelYields == MapProcChYields.end())continue;
      fprintf(pFile, "%s ", proc->c_str()); 
      for(auto ch = MapChannel.begin(); ch!=MapChannel.end();ch++){ 
        fprintf(pFile, " & ");
        if(ChannelYields->second.find(ch->first)!=ChannelYields->second.end()){
          fprintf(pFile, " %s", (ChannelYields->second)[ch->first].c_str());
        }
      }
      fprintf(pFile, "\\\\\n");
      if(*proc=="data")fprintf(pFile, "\\hline\n");             
    }
    fprintf(pFile,"\\hline\n");
    fprintf(pFile,"\\end{tabular}\n\\end{center}\n\\end{sidewaystable}\n");

    //All Bins
    fprintf(pFileInc,"\\begin{sidewaystable}[htp]\n\\begin{center}\n\\caption{Event yields expected for background and signal processes and observed in data.}\n\\label{tab:table}\n");
    fprintf(pFileInc, "\\begin{tabular}{|c|"); for(auto ch = MapChannelBin.begin(); ch!=MapChannelBin.end();ch++){ fprintf(pFileInc, "c|"); } fprintf(pFileInc, "}\\\\\\hline\n");
    fprintf(pFileInc, "channel");   for(auto ch = MapChannelBin.begin(); ch!=MapChannelBin.end();ch++){ fprintf(pFileInc, " & %s", ch->first.c_str()); } fprintf(pFileInc, "\\\\\\hline\n");
    for(auto proc = VectorProc.begin();proc!=VectorProc.end(); proc++){
      if(*proc=="total")fprintf(pFileInc, "\\hline\n");
      auto ChannelYields = MapProcChYieldsBin.find(*proc);
      if(ChannelYields == MapProcChYieldsBin.end())continue;
      fprintf(pFileInc, "%s ", proc->c_str()); 
      for(auto ch = MapChannelBin.begin(); ch!=MapChannelBin.end();ch++){ 
        fprintf(pFileInc, " & ");
        if(ChannelYields->second.find(ch->first)!=ChannelYields->second.end()){
          fprintf(pFileInc, " %s", (ChannelYields->second)[ch->first].c_str());
        }
      }
      fprintf(pFileInc, "\\\\\n");
      if(*proc=="data")fprintf(pFileInc, "\\hline\n");             
    }
    fprintf(pFileInc,"\\hline\n");
    fprintf(pFileInc,"\\end{tabular}\n\\end{center}\n\\end{sidewaystable}\n");



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
      if(Y->first.find("FakeLep")<std::string::npos)continue;//never drop this background
      if(Y->first.find("XH")<std::string::npos)continue;//never drop this background
      if(Y->second/total<threshold){
        printf("Drop %s from the list of backgrounds because of negligible rate (%f%% of total bckq)\n", Y->first.c_str(), Y->second/total);
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
    std::map<string, TGraphAsymmErrors*     > map_unc;
    std::map<string, TH1*              > map_data;
    std::map<string, TGraphAsymmErrors*> map_dataE;
    std::map<string, std::vector<TH1*> > map_signals;
    std::map<string, int               > map_legend;
		//           TLegend* legA  = new TLegend(0.6,0.5,0.99,0.85, "NDC");
		//           TLegend* legA  = new TLegend(0.03,0.00,0.97,0.70, "NDC");
		//           TLegend* legA  = new TLegend(0.03,0.99,0.97,0.89, "NDC");
    TLegend* legA  = new TLegend(0.08,0.89,0.97,0.95, "");
    std::vector<TLegendEntry*> legEntries;  //needed to have the entry in reverse order

    //order the proc first
    sortProc();

    //loop on sorted proc
    for(unsigned int p=0;p<sorted_procs.size();p++){
      string procName = sorted_procs[p];
      std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
	    TString process(procName.c_str());
      if( process.Contains("BOnly_B") || process.Contains("SandBandInterf_SBI") ) continue;
	    if(it==procs.end())continue;
      for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
        if(std::find(selCh.begin(), selCh.end(), ch->second.channel)==selCh.end())continue;
        if(ch->second.shapes.find(histoName.Data())==(ch->second.shapes).end())continue;
        TH1* h = ch->second.shapes[histoName.Data()].histo();
        //if(process.Contains("SOnly_S") ) h->Scale(10);  
        if(it->first=="total"){
          //double Uncertainty = std::max(0.0, ch->second.shapes[histoName.Data()].getScaleUncertainty() / h->Integral() );;
          double syst_scale = std::max(0.0, ch->second.shapes[histoName.Data()].getScaleUncertainty());
          //double syst_shape = std::max(0.0, ch->second.shapes[histoName.Data()].getBinShapeUncertainty((it->first+ch->first).c_str()));
          //double syst = sqrt(pow(syst_scale,2)+pow(syst_shape,2)); //On perd de l'info ici car on considere l'ecart constant au lieu d'y aller bin par bin
          //double Uncertainty = syst / h->Integral();
          double Uncertainty_scale = syst_scale / h->Integral();

          double Maximum = 0;
          TGraphAsymmErrors* errors = new TGraphAsymmErrors(h->GetXaxis()->GetNbins());
					//                    errors->SetFillStyle(3427);
					//                    errors->SetFillColor(kGray+1);
          errors->SetFillStyle(3005);
          errors->SetFillColor(kGray+3);                    
          errors->SetLineStyle(1);
          errors->SetLineColor(1);
          int icutg=0;
          for(int ibin=1; ibin<=h->GetXaxis()->GetNbins(); ibin++){
            if(h->GetBinContent(ibin)>0)
              errors->SetPoint(icutg,h->GetXaxis()->GetBinCenter(ibin), h->GetBinContent(ibin));
            //This is the part where we define which errors will be shown on the shape plot
            double syst_shape_binUp = std::max(0.0, ch->second.shapes[histoName.Data()].getBinShapeUncertainty((it->first+ch->first).c_str(), ibin, "Up"));
            double syst_shape_binDown = std::max(0.0, ch->second.shapes[histoName.Data()].getBinShapeUncertainty((it->first+ch->first).c_str(), ibin, "Down"));
		        double syst_binUp = sqrt(pow(Uncertainty_scale*h->GetBinContent(ibin), 2) + pow(syst_shape_binUp,2));
		        double syst_binDown = sqrt(pow(Uncertainty_scale*h->GetBinContent(ibin), 2) + pow(syst_shape_binDown,2));
		        double Uncertainty_binUp = syst_binUp / h->GetBinContent(ibin);
		        double Uncertainty_binDown = syst_binDown / h->GetBinContent(ibin);


            //errors->SetPointError(icutg,h->GetXaxis()->GetBinWidth(ibin)/2.0, sqrt(pow(h->GetBinContent(ibin)*Uncertainty_bin,2) + pow(h->GetBinError(ibin),2) ) );
            errors->SetPointError(icutg,h->GetXaxis()->GetBinWidth(ibin)/2.0,h->GetXaxis()->GetBinWidth(ibin)/2.0, sqrt(pow(h->GetBinContent(ibin)*Uncertainty_binDown,2) + pow(h->GetBinError(ibin),2) ), sqrt(pow(h->GetBinContent(ibin)*Uncertainty_binUp,2) + pow(h->GetBinError(ibin),2) ) );
						//                        printf("Unc=%6.2f  X=%6.2f Y=%6.2f+-%6.2f+-%6.2f=%6.2f\n", Uncertainty, h->GetXaxis()->GetBinCenter(ibin), h->GetBinContent(ibin), h->GetBinContent(ibin)*Uncertainty, h->GetBinError(ibin), sqrt(pow(h->GetBinContent(ibin)*Uncertainty,2) + pow(h->GetBinError(ibin),2) ) );
						//                        errors->SetPointError(icutg,h->GetXaxis()->GetBinWidth(ibin)/2.0, 0 );
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
          h->SetFillStyle(0);
          h->SetFillColor(0);
          h->SetMarkerSize(0.7);
          h->SetMarkerStyle(20);
          h->SetMarkerColor(1);
          h->SetBinErrorOption(TH1::kPoisson);
          map_data[ch->first] = h;

          //poisson error bars
          const double alpha = 1 - 0.6827;
          TGraphAsymmErrors * g = new TGraphAsymmErrors(h);
          g->SetMarkerSize(0.5);
          g->SetMarkerStyle (20);

          for (int i = 0; i < g->GetN(); ++i) {
            int N = g->GetY()[i];
            double L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
            double U =  ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) ;
            g->SetPointEYlow(i, N-L);
            g->SetPointEYhigh(i, U-N);
          }
          map_dataE[ch->first] = g;
        }

        if(map_legend.find(it->first)==map_legend.end()){
          map_legend[it->first]=1;
          if(it->first=="data"){
            legA->AddEntry(h,it->first.c_str(),"PE0");
          }else if(it->second.isSign){
            legEntries.insert(legEntries.begin(), new TLegendEntry(h, it->first.c_str(), "L") );
          }else{
            legEntries.push_back(new TLegendEntry(h, it->first.c_str(), "F") );
          }
          NLegEntry++;
        }
      }
    }
    //fill the legend in reverse order
    while(!legEntries.empty()){
      legA->AddEntry(legEntries.back()->GetObject(), legEntries.back()->GetLabel(), legEntries.back()->GetOption());
      legEntries.pop_back();      
    }
    if(map_unc.begin()!=map_unc.end())legA->AddEntry(map_unc.begin()->second, "Syst. + Stat.", "F");


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
      pad->SetTopMargin(0.06); pad->SetRightMargin(0.03); pad->SetBottomMargin(I<=NBins?0.09:0.12);  pad->SetLeftMargin((I-1)%NBins!=0?0.09:0.12);
      pad->SetLogy(true); 

      //print histograms
      TH1* axis = (TH1*)map_data[p->first]->Clone("axis");
      axis->Reset();      
      axis->GetXaxis()->SetRangeUser(axis->GetXaxis()->GetXmin(), axis->GetXaxis()->GetXmax());
      double signalHeight=0; for(unsigned int s=0;s<map_signals[p->first].size();s++){signalHeight = std::max(signalHeight, map_signals[p->first][s]->GetMaximum());}
      axis->SetMaximum(1.5*std::max(signalHeight , std::max( map_unc[p->first]->GetMaximum(), map_data[p->first]->GetMaximum())));

      //hard code the range for HZZ2l2nu
      if(procs["data"].channels[p->first].bin.find("vbf")!=string::npos){
        axis->SetMinimum(1E-4);
        axis->SetMaximum(std::max(axis->GetMaximum(), 5E1));
      }else{
        axis->SetMinimum(1E-4);
        axis->SetMaximum(std::max(axis->GetMaximum(), 5E1));
      }
      if((I-1)%NBins!=0)axis->GetYaxis()->SetTitle("");
      if(I<=NBins)axis->GetXaxis()->SetTitle("");
      axis->Draw();
      p->second->Draw("same");
      map_unc [p->first]->Draw("2 same");
      for(unsigned int i=0;i<map_signals[p->first].size();i++){
        map_signals[p->first][i]->Draw("HIST same");
      }
      if(!blindData) map_dataE[p->first]->Draw("P0 same");


      bool printBinContent = false;
      if(printBinContent){
        TLatex* tex = new TLatex();
        tex->SetTextSize(0.04); tex->SetTextFont(42);
        tex->SetTextAngle(60);
        TH1* histdata = map_data[p->first];
        double Xrange = axis->GetXaxis()->GetXmax()-axis->GetXaxis()->GetXmin();
        for(int xi=1;xi<=histdata->GetNbinsX();++xi){
          double x=histdata->GetBinCenter(xi);
          double y=histdata->GetBinContent(xi);
          double yData=histdata->GetBinContent(xi);

          int graphBin=-1;  for(int k=0;k<map_unc[p->first]->GetN();k++){if(fabs(map_unc[p->first]->GetX()[k] - x)<histdata->GetBinWidth(xi)){graphBin=k;}} 
          if(graphBin<0){
            printf("MC bin not found for X=%f\n", x);
          }else{
            double yMC =  map_unc[p->first]->GetY()[graphBin];
            double yMCerr = map_unc[p->first]->GetErrorY(graphBin);
            y = std::max(y, yMC+yMCerr);
            if(yMC>=1){tex->DrawLatex(x-0.02*Xrange,y*1.15,Form("#color[4]{B=%.1f#pm%.1f}",yMC, yMCerr));
            }else{     tex->DrawLatex(x-0.02*Xrange,y*1.15,Form("#color[4]{B=%.2f#pm%.2f}",yMC, yMCerr));
            }
          }
          tex->DrawLatex(x+0.02*Xrange,y*1.15,Form("D=%.0f",yData));                 
        }
      }

      //print tab channel header
      TPaveText* Label = new TPaveText(0.1,0.81,0.94,0.89, "NDC");
      Label->SetFillColor(0);  Label->SetFillStyle(0);  Label->SetLineColor(0); Label->SetBorderSize(0);  Label->SetTextAlign(31);
      TString LabelText = procs["data"].channels[p->first].channel+"  "+procs["data"].channels[p->first].bin;
      LabelText.ReplaceAll("eq","="); LabelText.ReplaceAll("l=","#leq");LabelText.ReplaceAll("g=","#geq"); 
      LabelText.ReplaceAll("_OS","OS "); LabelText.ReplaceAll("el","e"); LabelText.ReplaceAll("mu","#mu");  LabelText.ReplaceAll("ha","#tau_{had}");
      Label->AddText(LabelText);  Label->Draw();

      gPad->RedrawAxis();

      I++;
    }
    //print legend
    c1->cd(0);
    legA->SetFillColor(0); legA->SetFillStyle(0); legA->SetLineColor(0);  legA->SetBorderSize(0); legA->SetHeader("");
    legA->SetNColumns((NLegEntry/2) + 1);
    legA->Draw("same");    legA->SetTextFont(42);

    //print canvas header
		/*
       t2->cd(0);
		//           TPaveText* T = new TPaveText(0.1,0.995,0.84,0.95, "NDC");
    TPaveText* T = new TPaveText(0.1,0.7,0.9,1.0, "NDC");
    T->SetFillColor(0);  T->SetFillStyle(0);  T->SetLineColor(0); T->SetBorderSize(0);  T->SetTextAlign(22);
    if(systpostfix.Contains('3'))      { T->AddText("CMS preliminary, #sqrt{s}=13.0 TeV");
    }else if(systpostfix.Contains('8')){ T->AddText("CMS preliminary, #sqrt{s}=8.0 TeV");
    }else{                               T->AddText("CMS preliminary, #sqrt{s}=7.0 TeV");
    }T->Draw();

*/

    c1->cd(0);
    double L=0.03, R=0.03, T=0.02, B=0.0;
    char LumiText[1024];
    if(systpostfix.Contains('3'))      { double iLumi= 35914;sprintf(LumiText, "%.1f %s^{-1} (%.0f TeV)", iLumi>100?iLumi/1000:iLumi, iLumi>100?"fb":"pb", 13.0);
    }else if(systpostfix.Contains('8')){ double iLumi=20000;sprintf(LumiText, "%.1f %s^{-1} (%.0f TeV)", iLumi>100?iLumi/1000:iLumi, iLumi>100?"fb":"pb", 8.0);
    }else{                               double iLumi= 5000;sprintf(LumiText, "%.1f %s^{-1} (%.0f TeV)", iLumi>100?iLumi/1000:iLumi, iLumi>100?"fb":"pb", 7.0); 
    }
    TPaveText* T1 = new TPaveText(1.0-R-0.50, 1.0-T-0.05, 1.02-R, 1.0-T-0.005, "NDC");
    T1->SetTextFont(43); T1->SetTextSize(23);   T1->SetTextAlign(31);
    T1->SetFillColor(0); T1->SetFillStyle(0);   T1->SetBorderSize(0);
    T1->AddText(LumiText);  T1->Draw();

    //TOP LEFT IN-FRAME
    TPaveText* T2 = new TPaveText(L+0.005, 1.0-T-0.05, L+0.20, 1.0-T-0.005, "NDC");
    T2->SetTextFont(63); T2->SetTextSize(30);   T2->SetTextAlign(11);
    T2->SetFillColor(0); T2->SetFillStyle(0);   T2->SetBorderSize(0);
    T2->AddText("CMS"); T2->Draw();

    if(true){ //Right to CMS
      TPaveText* T3 = new TPaveText(L+0.095, 1.0-T-0.05, L+0.50, 1.0-T-0.005, "NDC");
      T3->SetTextFont(53); T3->SetTextSize(23);   T3->SetTextAlign(11);
      T3->SetFillColor(0); T3->SetFillStyle(0);   T3->SetBorderSize(0);
      T3->AddText("Preliminary"); T3->Draw();
    }




    //save canvas
    c1->SaveAs(SaveName+"_Shape.png");
    c1->SaveAs(SaveName+"_Shape.pdf");
    c1->SaveAs(SaveName+"_Shape.C");
    delete c1;
  }


  //
  // Make a summary plot
  //
  void AllInfo_t::showUncertainty(std::vector<TString>& selCh , TString histoName, TString SaveName)
  {
    string UncertaintyOnYield="";  char txtBuffer[4096];

    //loop on sorted proc
    for(unsigned int p=0;p<sorted_procs.size();p++){
      int NLegEntry = 0;
      std::map<string, int               > map_legend;
      std::vector<TH1*>                    toDelete;             
      TLegend* legA  = new TLegend(0.03,0.89,0.97,0.95, "");

      string procName = sorted_procs[p];
      std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
      if(it==procs.end())continue;
      if(it->first=="total" || it->first=="data")continue;  //only do samples which have systematics

      std::map<string, bool> mapUncType;
      std::map<string, std::map< string, double> > mapYieldPerBin;
      std::map<string, std::pair< double, double> > mapYieldInc;


      int NBins = it->second.channels.size()/selCh.size();
      TCanvas* c1 = new TCanvas("c1","c1",300*NBins,300*selCh.size());
      c1->SetTopMargin(0.00); c1->SetRightMargin(0.00); c1->SetBottomMargin(0.00);  c1->SetLeftMargin(0.00);
      TPad* t2 = new TPad("t2","t2", 0.03, 0.90, 1.00, 1.00, -1, 1);  t2->Draw();  c1->cd();
      t2->SetTopMargin(0.00); t2->SetRightMargin(0.00); t2->SetBottomMargin(0.00);  t2->SetLeftMargin(0.00);
      TPad* t1 = new TPad("t1","t1", 0.03, 0.03, 1.00, 0.90, 4, 1);  t1->Draw();  t1->cd();
      t1->SetTopMargin(0.00); t1->SetRightMargin(0.00); t1->SetBottomMargin(0.00);  t1->SetLeftMargin(0.00);
      t1->Divide(NBins, selCh.size(), 0, 0);

      int I=1;
      mapYieldInc[""].first = 0;  mapYieldInc[""].second = 0;
      for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++, I++){
        if(std::find(selCh.begin(), selCh.end(), ch->second.channel)==selCh.end())continue;
        if(ch->second.shapes.find(histoName.Data())==(ch->second.shapes).end())continue;

        //add the stat uncertainty is there;
        //ch->second.shapes[histoName.Data()].makeStatUnc("_CMS_hzz2l2v_", (TString("_")+ch->first+"_"+it->second.shortName).Data(),systpostfix.Data(), false );//add stat uncertainty to the uncertainty map;

        //Li Fix
        if((it->second.shortName).find("ggH")!=std::string::npos)ch->second.shapes[histoName.Data()].makeStatUnc("_CMS_hzz2l2v_", (TString("_")+ch->first+TString("_ggH")).Data(),systpostfix.Data(), false );// attention
	else if((it->second.shortName).find("qqH")!=std::string::npos)ch->second.shapes[histoName.Data()].makeStatUnc("_CMS_hzz2l2v_", (TString("_")+ch->first+TString("_qqH")).Data(),systpostfix.Data(), false );
        else ch->second.shapes[histoName.Data()].makeStatUnc("_CMS_hzz2l2v_", (TString("_")+ch->first+"_"+it->second.shortName).Data(),systpostfix.Data(), false );
        TVirtualPad* pad = t1->cd(I); 
        pad->SetTopMargin(0.06); pad->SetRightMargin(0.03); pad->SetBottomMargin(0.07);  pad->SetLeftMargin(0.06);
				//                 pad->SetLogy(true); 

        TH1* h = (TH1*)(ch->second.shapes[histoName.Data()].histo()->Clone((it->first+ch->first+"Nominal").c_str())); 
        double yield = h->Integral();
        toDelete.push_back(h);
        mapYieldPerBin[""][ch->first] = yield;
        mapYieldInc[""].first  += yield;
        mapYieldInc[""].second = 1;

        //print histograms
        TH1* axis = (TH1*)h->Clone("axis");
        axis->Reset();      
        axis->GetXaxis()->SetRangeUser(axis->GetXaxis()->GetXmin(), axis->GetXaxis()->GetXmax());
        axis->GetYaxis()->SetRangeUser(0.5, 1.5); //100% uncertainty
        if(it->second.shortName == "instrmet") axis->GetYaxis()->SetRangeUser(-5, 5);
        if((I-1)%NBins!=0)axis->GetYaxis()->SetTitle("");
        axis->Draw();
        toDelete.push_back(axis);


        //print tab channel header
        TPaveText* Label = new TPaveText(0.1,0.81,0.94,0.89, "NDC");
        Label->SetFillColor(0);  Label->SetFillStyle(0);  Label->SetLineColor(0); Label->SetBorderSize(0);  Label->SetTextAlign(31);
        TString LabelText = ch->second.channel+"  -  "+ch->second.bin;
        LabelText.ReplaceAll("eq","="); LabelText.ReplaceAll("l=","#leq");LabelText.ReplaceAll("g=","#geq"); 
        LabelText.ReplaceAll("_OS","OS "); LabelText.ReplaceAll("el","e"); LabelText.ReplaceAll("mu","#mu");  LabelText.ReplaceAll("ha","#tau_{had}");
        Label->AddText(LabelText);  Label->Draw();

        TLine* line = new TLine(axis->GetXaxis()->GetXmin(), 1.0, axis->GetXaxis()->GetXmax(), 1.0);
        toDelete.push_back((TH1*)line);
        line->SetLineWidth(2);  line->SetLineColor(1); line->Draw("same");
        if(I==1){legA->AddEntry(line,"Nominal","L");  NLegEntry++;}

        int ColorIndex=3;

        //draw scale uncertainties
        for(std::map<string, double>::iterator var = ch->second.shapes[histoName.Data()].uncScale.begin(); var!=ch->second.shapes[histoName.Data()].uncScale.end(); var++){
          if(h->Integral()<=0)continue;
          double ScaleChange   = var->second/h->Integral();
          double ScaleUp   = 1 + ScaleChange;
          double ScaleDn   = 1 - ScaleChange;

          TString systName = var->first.c_str();
          systName.ToLower();
          systName.ReplaceAll("cms","");
          systName.ReplaceAll("hzz2l2v","");
          systName.ReplaceAll("sys","");
          systName.ReplaceAll("13tev","");
          systName.ReplaceAll("_","");
          systName.ReplaceAll("up","");
          systName.ReplaceAll("down","");

          TLine* lineUp = new TLine(axis->GetXaxis()->GetXmin(), ScaleUp, axis->GetXaxis()->GetXmax(), ScaleUp);
          TLine* lineDn = new TLine(axis->GetXaxis()->GetXmin(), ScaleDn, axis->GetXaxis()->GetXmax(), ScaleDn);
          toDelete.push_back((TH1*)lineUp);  toDelete.push_back((TH1*)lineDn);


          int color = ColorIndex;
          if(map_legend.find(systName.Data())==map_legend.end()){
            map_legend[systName.Data()]=color;
            legA->AddEntry(lineUp,systName.Data(),"L");
            NLegEntry++;
            ColorIndex++;
          }else{
            color = map_legend[systName.Data()];
          }

          if(mapYieldPerBin[systName.Data()].find(ch->first)==mapYieldPerBin[systName.Data()].end()){
            mapYieldPerBin[systName.Data()][ch->first] = fabs( ScaleChange );                       
            mapUncType[systName.Data()] = false;
          }else{
            mapYieldPerBin[systName.Data()][ch->first] = std::max( ScaleChange , mapYieldPerBin[systName.Data()][ch->first]);
          }

          if(mapYieldInc.find(systName.Data())==mapYieldInc.end()){ mapYieldInc[systName.Data()].first = 0; mapYieldInc[systName.Data()].second = 0;  }
          mapYieldInc[systName.Data()].first += ScaleChange*h->Integral();
          mapYieldInc[systName.Data()].second += h->Integral();

          lineUp->SetLineWidth(2);  lineUp->SetLineColor(color); lineUp->Draw("same");
          lineDn->SetLineWidth(2);  lineDn->SetLineColor(color); lineDn->Draw("same");
        }

        //draw shape uncertainties
        //double syst_shape = std::max(0.0, ch->second.shapes[histoName.Data()].getShapeUncertainty((it->first+ch->first).c_str()));
        for(std::map<string, TH1*>::iterator var = ch->second.shapes[histoName.Data()].uncShape.begin(); var!=ch->second.shapes[histoName.Data()].uncShape.end(); var++){
          if(var->first=="")continue;

          TH1* hvar = (TH1*)(var->second->Clone((it->first+ch->first+var->first).c_str())); 
          double varYield = hvar->Integral();
          hvar->Divide(h);
          toDelete.push_back(hvar);

          TString systName = var->first.c_str();
          systName.ToLower();
          systName.ReplaceAll("cms","");
          systName.ReplaceAll("hzz2l2v","");
          systName.ReplaceAll("sys","");
          systName.ReplaceAll("13tev","");
          systName.ReplaceAll("_","");
          systName.ReplaceAll("up","");
          systName.ReplaceAll("down","");

          int color = ColorIndex;
          if(systName.Contains("stat")){systName = "stat"; color=2;}
          if(map_legend.find(systName.Data())==map_legend.end()){
            map_legend[systName.Data()]=color;
            legA->AddEntry(hvar,systName.Data(),"L");
            NLegEntry++;
            ColorIndex++;
          }else{
            color = map_legend[systName.Data()];
          }

          if(yield>0){
            if(mapYieldPerBin[systName.Data()].find(ch->first)==mapYieldPerBin[systName.Data()].end()){
              mapYieldPerBin[systName.Data()][ch->first] = fabs( 1 - (varYield/yield));
              mapUncType[systName.Data()] = true;                        
            }else{
              mapYieldPerBin[systName.Data()][ch->first] = std::max(fabs( 1 - (varYield/yield) ), mapYieldPerBin[systName.Data()][ch->first]);
            }

            if(mapYieldInc.find(systName.Data())==mapYieldInc.end()){ mapYieldInc[systName.Data()].first = 0; mapYieldInc[systName.Data()].second = 0;  }
            mapYieldInc[systName.Data()].first +=  fabs( 1 - (varYield/yield))*yield;
            mapYieldInc[systName.Data()].second += yield;
          }

          hvar->SetFillColor(0);                  
          hvar->SetLineStyle(1);
          hvar->SetLineColor(color);
          hvar->SetLineWidth(2);
          hvar->Draw("HIST same");                   
        }
        //remove the stat uncertainty
   			ch->second.shapes[histoName.Data()].removeStatUnc(); 
      }
      //print legend
      c1->cd(0);
      legA->SetFillColor(0); legA->SetFillStyle(0); legA->SetLineColor(0);  legA->SetBorderSize(0); legA->SetHeader("");
      legA->SetNColumns((NLegEntry/2) + 1);
      legA->Draw("same");    legA->SetTextFont(42);

      //print canvas header
      t2->cd(0);
      TPaveText* T = new TPaveText(0.1,0.7,0.9,1.0, "NDC");
      T->SetFillColor(0);  T->SetFillStyle(0);  T->SetLineColor(0); T->SetBorderSize(0);  T->SetTextAlign(22);
      if(systpostfix.Contains('3'))      { T->AddText((string("CMS preliminary, #sqrt{s}=13.0 TeV,   ")+it->first).c_str());
      }else if(systpostfix.Contains('8')){ T->AddText((string("CMS preliminary, #sqrt{s}=8.0 TeV,   ")+it->first).c_str());
      }else{                               T->AddText((string("CMS preliminary, #sqrt{s}=7.0 TeV,   ")+it->first).c_str());
      }T->Draw();

      //save canvas
      c1->SaveAs(SaveName+"_Uncertainty_"+it->second.shortName+".png");
      c1->SaveAs(SaveName+"_Uncertainty_"+it->second.shortName+".pdf");
      c1->SaveAs(SaveName+"_Uncertainty_"+it->second.shortName+".C");
      delete c1;             

      for(unsigned int i=0;i<toDelete.size();i++){delete toDelete[i];} //clear the objects


      //add inclusive uncertainty as a channel
      //
      for(auto systIt=mapYieldInc.begin(); systIt!=mapYieldInc.end(); systIt++){ mapYieldPerBin[systIt->first][" Inc"] = systIt->second.first/systIt->second.second;  }
      //print uncertainty on yield            
      //
      sprintf(txtBuffer, "\\multicolumn{%i}{'c'}{\\bf{%s}}\\\\ \n", I+1, it->first.c_str());  UncertaintyOnYield+= txtBuffer;
      sprintf(txtBuffer, "%10s & %25s", "Type", "Uncertainty");
      for(auto chIt=mapYieldPerBin[""].begin();chIt!=mapYieldPerBin[""].end();chIt++){ sprintf(txtBuffer, "%s & %12s ", txtBuffer, chIt->first.c_str()); } sprintf(txtBuffer, "%s\\\\ \\hline\n", txtBuffer);  UncertaintyOnYield += txtBuffer;
      sprintf(txtBuffer, "%10s & %25s", "", "Nominal yields ");
      for(auto chIt=mapYieldPerBin[""].begin();chIt!=mapYieldPerBin[""].end();chIt++){ sprintf(txtBuffer, "%s & %12.4E ", txtBuffer, chIt->second); } sprintf(txtBuffer, "%s\\\\ \n", txtBuffer);  UncertaintyOnYield += txtBuffer;
      for(auto varIt=mapYieldPerBin.begin();varIt!=mapYieldPerBin.end();varIt++){
        if(varIt->first == "")continue;
        sprintf(txtBuffer, "%10s & %25s", mapUncType[varIt->first.c_str()]?"shape":"scale", varIt->first.c_str());
        for(auto chIt=mapYieldPerBin[""].begin();chIt!=mapYieldPerBin[""].end();chIt++){
          if(varIt->second.find(chIt->first)==varIt->second.end()){ sprintf(txtBuffer, "%s & %10s   "    , txtBuffer, "-" );                        
          }else{                                                    sprintf(txtBuffer, "%s & %+10.3f\\%% ", txtBuffer,100.0 * varIt->second[chIt->first] ); 
          }
        }sprintf(txtBuffer, "%s\\\\ \n", txtBuffer);   UncertaintyOnYield += txtBuffer;
      }sprintf(txtBuffer, "\\hline \n"); UncertaintyOnYield += txtBuffer;
    }

    FILE* pFile = fopen(SaveName+"_Uncertainty.txt", "w");
    if(pFile){ fprintf(pFile, "%s\n", UncertaintyOnYield.c_str()); fclose(pFile);}
  }




  //
  // Turn to cut&count (rebin all histo to 1 bin only)
  //
  void AllInfo_t::turnToCC(string histoName){
    //order the proc first
    sortProc();

    //Loop on processes and channels
    for(unsigned int p=0;p<sorted_procs.size();p++){
      string procName = sorted_procs[p];
      std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
      if(it==procs.end())continue;
      for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
        TString chbin = ch->first;
        if(ch->second.shapes.find(histoName)==(ch->second.shapes).end())continue;
        ShapeData_t& shapeInfo = ch->second.shapes[histoName];      
        TH1* h = shapeInfo.histo();

        TString proc = it->second.shortName.c_str();
        for(std::map<string, TH1*  >::iterator unc=shapeInfo.uncShape.begin();unc!=shapeInfo.uncShape.end();unc++){
          TString syst   = unc->first.c_str();
          TH1*    hshape = unc->second;
          hshape->SetDirectory(0);

          hshape = hshape->Rebin(hshape->GetXaxis()->GetNbins()); 
          //make sure to also count the underflow and overflow
          double bin  = hshape->GetBinContent(0) + hshape->GetBinContent(1) + hshape->GetBinContent(2);
          double bine = sqrt(hshape->GetBinError(0)*hshape->GetBinError(0) + hshape->GetBinError(1)*hshape->GetBinError(1) + hshape->GetBinError(2)*hshape->GetBinError(2));
          hshape->SetBinContent(0,0);              hshape->SetBinError  (0,0);
          hshape->SetBinContent(1,bin);            hshape->SetBinError  (1,bine);
          hshape->SetBinContent(2,0);              hshape->SetBinError  (2,0);
        }
      }
    }
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
				//                 shapeInfo.makeStatUnc("_CMS_hzz2l2v_", (TString("_")+ch->first+"_"+it->second.shortName).Data(),systpostfix.Data(), it->second.isSign );//add stat uncertainty to the uncertainty map;
        //shapeInfo.makeStatUnc("_CMS_hzz2l2v_", (TString("_")+ch->first+"_"+it->second.shortName).Data(),systpostfix.Data(), false );//add stat uncertainty to the uncertainty map;

        //Li Fix
        if((it->second.shortName).find("ggH")!=std::string::npos)shapeInfo.makeStatUnc("_CMS_hzz2l2v_", (TString("_")+ch->first+TString("_ggH")).Data(),systpostfix.Data(), false );// attention
	else if((it->second.shortName).find("qqH")!=std::string::npos)shapeInfo.makeStatUnc("_CMS_hzz2l2v_", (TString("_")+ch->first+TString("_qqH")).Data(),systpostfix.Data(), false );
        else shapeInfo.makeStatUnc("_CMS_hzz2l2v_", (TString("_")+ch->first+"_"+it->second.shortName).Data(),systpostfix.Data(), false );
				fout->cd(chbin);

        TString proc = it->second.shortName.c_str();
        for(std::map<string, TH1*  >::iterator unc=shapeInfo.uncShape.begin();unc!=shapeInfo.uncShape.end();unc++){
          TString syst   = unc->first.c_str();
          TH1*    hshape = unc->second;
          hshape->SetDirectory(0);

          if(syst==""){
            //central shape (for data call it data_obs)
            hshape->SetName(proc); 
            if(it->first=="data"){
              hshape->Write("data_obs");
            }else{
              hshape->Write(proc+postfix);
            }
          }else if(runSystematics && proc!="data" && (syst.Contains("Up") || syst.Contains("Down"))){
            //if empty histogram --> no variation is applied except for stat
            if(!syst.Contains("stat") && (hshape->Integral()<h->Integral()*0.01 || isnan((float)hshape->Integral()))){hshape->Reset(); hshape->Add(h,1); }

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
  // add hardcoded uncertainties
  //
  void AllInfo_t::addHardCodedUncertainties(string histoName){
    //7/8 TeV values
    //#FIXME: extrapolated from 600 to 1TeVmissing points from 650GeV to 1TeV
		//            double QCDScaleMass   [] = {200, 250, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900, 1000};
		//            double QCDScaleK0ggH0 [] = {1.15 , 1.16, 1.17, 1.20, 1.17, 1.19, 1.22, 1.24, 1.25, 1.25, 1.25, 1.25, 1.25};
		//            double QCDScaleK0ggH1 [] = {0.88, 0.86, 0.84, 0.83, 0.82, 0.81, 0.80, 0.78, 0.78, 0.78, 0.78, 0.78, 0.78};
		//            double QCDScaleK1ggH1 [] = {1.27, 1.27, 1.27, 1.27, 1.26, 1.26, 1.25, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26};
		//            double QCDScaleK1ggH2 [] = {0.96, 0.96, 0.95, 0.95, 0.95, 0.95, 0.95,  0.95, 0.94, 0.94, 0.94, 0.94, 0.94};
		//            double QCDScaleK2ggH2 [] = { 1.20, 1.17, 1.20, 1.21, 1.20, 1.20, 1.17, 1.19, 1.19, 1.19, 1.19, 1.19, 1.19};


    //13TeV values  
    double QCDScaleMass        [] = {  200,   300,   400,   600,   800,  1000,  1500,  2000,  2500,  3000,  9999};
    double QCDScaleggHeq0jets  [] = {2.042, 1.416, 1.283, 1.335, 1.352, 1.425, 1.542, 1.627, 1.659, 1.638, 1.638};
    double QCDScaleggHgeq1jets [] = {1.305, 1.219, 1.181, 1.168, 1.172, 1.183, 1.201, 1.219, 1.220, 1.219, 1.219};
    double QCDScaleggHvbf      [] = {1.217, 1.282, 1.200, 1.215, 1.178, 1.195, 1.201, 1.226, 1.232, 1.199, 1.199};

    double UEPSf0 []         = {0.952, 0.955, 0.958, 0.964, 0.966, 0.954, 0.946, 0.931, 0.920, 0.920, 0.920, 0.920, 0.920};
    double UEPSf1 []         = {1.055, 1.058, 1.061, 1.068, 1.078, 1.092, 1.102, 1.117, 1.121, 1.121, 1.121, 1.121, 1.121};
    double UEPSf2 []         = {1.059, 0.990, 0.942, 0.889, 0.856, 0.864, 0.868, 0.861, 0.872, 0.872, 0.872, 0.872, 0.872}; 

    TGraph* TG_QCDScaleggHeq0jets  = new TGraph(sizeof(QCDScaleMass)/sizeof(double), QCDScaleMass,  QCDScaleggHeq0jets);
    TGraph* TG_QCDScaleggHgeq1jets = new TGraph(sizeof(QCDScaleMass)/sizeof(double), QCDScaleMass, QCDScaleggHgeq1jets);
    TGraph* TG_QCDScaleggHvbf      = new TGraph(sizeof(QCDScaleMass)/sizeof(double), QCDScaleMass,      QCDScaleggHvbf);
    TGraph* TG_UEPSf0         = new TGraph(sizeof(QCDScaleMass)/sizeof(double), QCDScaleMass, UEPSf0);
    TGraph* TG_UEPSf1         = new TGraph(sizeof(QCDScaleMass)/sizeof(double), QCDScaleMass, UEPSf1);
    TGraph* TG_UEPSf2         = new TGraph(sizeof(QCDScaleMass)/sizeof(double), QCDScaleMass, UEPSf2);

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
        if(!it->second.isData && systpostfix.Contains('3'))shapeInfo.uncScale["lumi_13TeV"] = integral*0.026;
        if(!it->second.isData && systpostfix.Contains('8'))shapeInfo.uncScale["lumi_8TeV" ] = integral*0.026;
        if(!it->second.isData && systpostfix.Contains('7'))shapeInfo.uncScale["lumi_7TeV" ] = integral*0.022;

        //Id+Trigger efficiencies combined
        if(!it->second.isData){
          if(chbin.Contains("ee"  ))  shapeInfo.uncScale["CMS_eff_e"] = integral*0.072124;
          if(chbin.Contains("mumu"))  shapeInfo.uncScale["CMS_eff_m"] = integral*0.061788;
        }

        //uncertainties to be applied only in higgs analyses
        if(mass>0){
          //bin migration at th level
          if(it->second.shortName.find("ggH")!=string::npos && chbin.Contains("eq0jet" )){shapeInfo.uncScale["QCDscale_ggH"]    = integral*(TG_QCDScaleggHeq0jets->Eval(mass,NULL,"S")-1);}
          if(it->second.shortName.find("ggH")!=string::npos && chbin.Contains("eq1jet" )){shapeInfo.uncScale["QCDscale_ggH"] = integral*(TG_QCDScaleggHgeq1jets->Eval(mass,NULL,"S")-1);} 
          if(it->second.shortName.find("ggH")!=string::npos && chbin.Contains("vbf"    )){shapeInfo.uncScale["QCDscale_ggH"] = integral*(TG_QCDScaleggHvbf->Eval(mass,NULL,"S")-1);}

        }//end of uncertainties to be applied only in higgs analyses

        if(it->second.shortName.find("zz")!=string::npos && chbin.Contains("eq0jet" )){shapeInfo.uncScale["QCDscale_ZZ"]    = integral*0.063;}
        if(it->second.shortName.find("zz")!=string::npos && chbin.Contains("eq1jet" )){shapeInfo.uncScale["QCDscale_ZZ"]    = integral*0.054;}
        if(it->second.shortName.find("zz")!=string::npos && chbin.Contains("vbf" )){shapeInfo.uncScale["QCDscale_ZZ"]    = integral*0.40;}
        if(it->second.shortName.find("wz")!=string::npos && chbin.Contains("eq0jet" )){shapeInfo.uncScale["QCDscale_WZ"]    = integral*0.095;}
        if(it->second.shortName.find("wz")!=string::npos && chbin.Contains("eq1jet" )){shapeInfo.uncScale["QCDscale_WZ"]    = integral*0.051;}
        if(it->second.shortName.find("wz")!=string::npos && chbin.Contains("vbf" )){shapeInfo.uncScale["QCDscale_WZ"]    = integral*0.40;}

      }
    }
  }



  //
  // produce the datacards 
  //
  void AllInfo_t::buildDataCards(string histoName, TString url)
  {
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

      dcName.ReplaceAll("[", "");
      dcName.ReplaceAll("]", "");
      dcName.ReplaceAll("+", "");

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
      fprintf(pFile, "Observation %f\n", procs["data"].channels[C->first].shapes[histoName].histo()->Integral());
      fprintf(pFile, "-------------------------------\n");

      //yields
      fprintf(pFile,"%55s  ", "bin");     for(unsigned int j=0; j<clean_procs.size(); j++){ fprintf(pFile,"%8i ", 1)                     ;}  fprintf(pFile,"\n");
      fprintf(pFile,"%55s  ", "process"); for(unsigned int j=0; j<clean_procs.size(); j++){ fprintf(pFile,"%8s ", procs[clean_procs[j]].shortName.c_str());}  fprintf(pFile,"\n");
      fprintf(pFile,"%55s  ", "process"); for(unsigned int j=0; j<clean_procs.size(); j++){ fprintf(pFile,"%8i ", ((int)j)-(nsign-1)    );}  fprintf(pFile,"\n");
      fprintf(pFile,"%55s  ", "rate");    for(unsigned int j=0; j<clean_procs.size(); j++){ fprintf(pFile,"%8f ", procs[clean_procs[j]].channels[C->first].shapes[histoName].histo()->Integral() );}  fprintf(pFile,"\n");
      fprintf(pFile, "-------------------------------\n");

      for(std::map<string, bool>::iterator U=allSysts.begin(); U!=allSysts.end();U++){
        if(mass==125 && U->first=="CMS_hzz2l2v_lshape")continue;//skip lineshape uncertainty for 125GeV Higgs

        char line[2048];
        sprintf(line,"%-45s %-10s ", U->first.c_str(), U->second?"shapeN2":"lnN");
        bool isNonNull = false;
        for(unsigned int j=0; j<clean_procs.size(); j++){
          ShapeData_t& shapeInfo = procs[clean_procs[j]].channels[C->first].shapes[histoName];
          double integral = shapeInfo.histo()->Integral();
          if(shapeInfo.uncScale.find(U->first)!=shapeInfo.uncScale.end()){   isNonNull = true;   
            if(U->second)                                                   sprintf(line,"%s%8s ",line,"       1");
            else if(integral>0)                                             sprintf(line,"%s%8f ",line,1+(shapeInfo.uncScale[U->first]/integral));
            else                                                            sprintf(line,"%s%8f ",line,1+(shapeInfo.uncScale[U->first]));
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
      string matchingKeyword="";
      if(!utils::root::getMatchingKeyword(Process[i], keywords, matchingKeyword))continue; //only consider samples passing key filtering


      TString procCtr(""); procCtr+=i;
      TString proc=Process[i].getString("tag", "noTagFound");

      string dirName = proc.Data(); 
	    //std::<TString> keys = Process[i].getString("keys", "noKeysFound");
      if(Process[i].isTagFromKeyword(matchingKeyword, "mctruthmode") ) { char buf[255]; sprintf(buf,"_filt%d",(int)Process[i].getIntFromKeyword(matchingKeyword, "mctruthmode", 0)); dirName += buf; }
      string procSuffix = Process[i].getStringFromKeyword(matchingKeyword, "suffix", "");
      if(procSuffix!=""){dirName += "_" + procSuffix;}
      while(dirName.find("/")!=std::string::npos)dirName.replace(dirName.find("/"),1,"-");         

      TDirectory *pdir = (TDirectory *)inF->Get(dirName.c_str());         
      if(!pdir){printf("Directory (%s) for proc=%s is not in the file!\n", dirName.c_str(), proc.Data()); continue;}

      bool isData = Process[i].getBool("isdata", false);
      if(onlyData && !isData)continue; //just here to speedup the NRB prediction
      if(proc.Contains(")cp0"))continue; // skip those samples

      bool isSignal = Process[i].getBool("issignal", false);
      if(Process[i].getBool("spimpose", false) && (proc.Contains("ggH") || proc.Contains("qqH")))isSignal=true;
      //LQ bool isInSignal = Process[i].getBool("isinsignal", false);
      int color = Process[i].getInt("color", 1);
      int lcolor = Process[i].getInt("lcolor", 1);
      int mcolor = Process[i].getInt("mcolor", color);
      int lwidth = Process[i].getInt("lwidth", 1);
      int lstyle = Process[i].getInt("lstyle", 1);
      int fill   = Process[i].getInt("fill"  , 1001);
      int marker = Process[i].getInt("marker", 20);

      if(isSignal && signalTag!=""){
        if(!proc.Contains(signalTag.c_str()) )continue;
      }

      double procMass=0;  char procMassStr[128] = "";
      if(isSignal &&  mass>0 && (proc.Contains("H(") || proc.Contains("h(") || proc.Contains("A(") || proc.Contains("Rad(") || proc.Contains("RsGrav(") || proc.Contains("BulkGrav(") )){
        if(proc.Contains("H(") && proc.Contains("A(")){sscanf(proc.Data()+proc.First("A")+2,"%lf",&procMass);
        }else if(proc.Contains("H(")){sscanf(proc.Data()+proc.First("H")+2,"%lf",&procMass);
        }else if(proc.Contains("A(")){sscanf(proc.Data()+proc.First("A")+2,"%lf",&procMass);
        }else if(proc.Contains("h(")){sscanf(proc.Data()+proc.First("(")+1,"%lf",&procMass);
        }else if(proc.Contains("Rad(")){sscanf(proc.Data()+proc.First("(")+1,"%lf",&procMass);
		  	}else if(proc.Contains("RsGrav(")){sscanf(proc.Data()+proc.First("(")+1,"%lf",&procMass);
        }else if(proc.Contains("BulkGrav(")){sscanf(proc.Data()+proc.First("(")+1,"%lf",&procMass);}

        //printf("%s --> %f\n",  proc.Data(), procMass);

        //skip signal sample not needed
        if(massL!=-1 && massR!=-1){
          if(procMass!=massL && procMass!=massR)continue; 
        }else{
          if(procMass!=mass)continue;
        }
        sprintf(procMassStr,"%i",(int)procMass);
        //printf("found signal to be %s\n",  proc.Data());
      }

      if(!isSignal &&  mass>0 && proc.Contains("XH(") && proc.Contains(")#rightarrow WW")){
        sscanf(proc.Data()+proc.First("H(")+2,"%lf",&procMass);
        if(!(procMass==mass || procMass==massL || procMass==massR))continue; //skip XH-->WW background sample not concerned
      }

      TString procSave = proc;
      if(isSignal && mass>0 && proc.Contains("ggH") && proc.Contains("ZZ"))proc = TString("ggH")  +procMassStr;
      else if(isSignal && mass>0 && proc.Contains("qqH") && proc.Contains("ZZ"))proc = TString("qqH")  +procMassStr;
      else if(isSignal && mass>0 && proc.Contains("ggH") && proc.Contains("WW"))proc = TString("ggHWW")+procMassStr;
      else if(isSignal && mass>0 && proc.Contains("qqH") && proc.Contains("WW"))proc = TString("qqHWW")+procMassStr;

      if(procSave.Contains("SandBandInterf"))proc+="_SBI";
      else if(procSave.Contains("SOnly"))         proc+="_S";
      else if(procSave.Contains("BOnly"))         proc+="_B";

      if(skipGGH && isSignal && mass>0 && proc.Contains("ggH") )continue;
      if(skipQQH && isSignal && mass>0 && (proc.Contains("qqH") || proc.Contains("VBF")) )continue;

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

      if(procInfo.isSign){
        procInfo.xsec = procInfo.jsonObj["data"].daughters()[0].getDouble("xsec", 1);
        if(procInfo.jsonObj["data"].daughters()[0].isTag("br")){
          std::vector<JSONWrapper::Object> BRs = procInfo.jsonObj["data"].daughters()[0]["br"].daughters();
          double totalBR=1.0; for(size_t ipbr=0; ipbr<BRs.size(); ipbr++){totalBR*=BRs[ipbr].toDouble();}   
          procInfo.br = totalBR;
        }
      }

      //Loop on all channels, bins and shape to load and store them in memory structure
      TH1* syst = (TH1*)pdir->Get("all_optim_systs");
      if(syst==NULL){syst=new TH1F("all_optim_systs","all_optim_systs",1,0,1);syst->GetXaxis()->SetBinLabel(1,"");}
      for(unsigned int c=0;c<channelsAndShapes.size();c++){
        TString chName    = (channelsAndShapes[c].substr(0,channelsAndShapes[c].find(";"))).c_str();
        TString binName   = (channelsAndShapes[c].substr(channelsAndShapes[c].find(";")+1, channelsAndShapes[c].rfind(";")-channelsAndShapes[c].find(";")-1)).c_str();
        TString shapeName = (channelsAndShapes[c].substr(channelsAndShapes[c].rfind(";")+1)).c_str();
        TString ch        = chName+binName;

        ChannelInfo_t& channelInfo = procInfo.channels[ch.Data()];
        channelInfo.bin        = binName.Data();
        channelInfo.channel    = chName.Data();
        ShapeData_t& shapeInfo = channelInfo.shapes[shapeName.Data()];

        //printf("%s SYST SIZE=%i\n", (ch+"_"+shapeName).Data(), syst->GetNbinsX() );
        for(int ivar = 1; ivar<=syst->GetNbinsX();ivar++){                
          TH1D* hshape   = NULL;
          TString varName   = syst->GetXaxis()->GetBinLabel(ivar);
          TString histoName = ch+"_"+shapeName+(isSignal?signalSufix:"")+varName ;
          if(shapeName==histo && histoVBF!="" && ch.Contains("vbf"))histoName = ch+"_"+histoVBF+(isSignal?signalSufix:"")+varName ;
          //if(isSignal && ivar==1)printf("Syst %i = %s\n", ivar, varName.Data()); 

          TH2* hshape2D = (TH2*)pdir->Get(histoName );
          if(!hshape2D){
            if(shapeName==histo && histoVBF!="" && ch.Contains("vbf")){   hshape2D = (TH2*)pdir->Get(TString("all_")+histoVBF+(isSignal?signalSufix:"")+varName);
            }else{                                                        hshape2D = (TH2*)pdir->Get(TString("all_")+shapeName+varName);
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
          if(shapeName == histo && !ch.Contains("vbf") && procMass==massL)cutBinUsed = indexcutML[channelInfo.bin];
          if(shapeName == histo && !ch.Contains("vbf") && procMass==massR)cutBinUsed = indexcutMR[channelInfo.bin];

          histoName.ReplaceAll(ch,ch+"_proj"+procCtr);
          hshape   = hshape2D->ProjectionY(histoName,cutBinUsed,cutBinUsed);
          filterBinContent(hshape);
          //printf("%s %s %s Integral = %f\n", ch.Data(), shortName.Data(), varName.Data(), hshape->Integral() );
          //if(hshape->Integral()<=0 && varName=="" && !isData){hshape->Reset(); hshape->SetBinContent(1, 1E-10);} //TEST FOR HIGGS WIDTH MEASUREMENTS, MUST BE UNCOMMENTED ASAP

          if(isnan((float)hshape->Integral())){hshape->Reset();}
          hshape->SetDirectory(0);
          hshape->SetTitle(proc);
          utils::root::fixExtremities(hshape,false,true);
          hshape->SetFillColor(color); hshape->SetLineColor(lcolor); hshape->SetMarkerColor(mcolor);
          hshape->SetFillStyle(fill);  hshape->SetLineWidth(lwidth); hshape->SetMarkerStyle(marker); hshape->SetLineStyle(lstyle);

          //if current shape is the one to cut on, then apply the cuts
          if(shapeName == histo){
            //if(ivar==1 && isSignal)printf("A %s %s Integral = %f\n", ch.Data(), shortName.Data(), hshape->Integral() );

            for(int x=0;x<=hshape->GetXaxis()->GetNbins()+1;x++){
              if(hshape->GetXaxis()->GetBinCenter(x)<=minCut || hshape->GetXaxis()->GetBinCenter(x)>=maxCut){ hshape->SetBinContent(x,0); hshape->SetBinError(x,0); }
            }

            if(rebinVal>1){ hshape->Rebin(rebinVal); }
						//                     if(histoVBF!="" && ch.Contains("vbf")){    hshape->Rebin(30);
						//                     }else{                                     hshape->Rebin(rebinVal);
						//                     }
            hshape->GetYaxis()->SetTitle("Entries");// (/25GeV)");

          }
          hshape->Scale(MCRescale);
          if(isSignal)hshape->Scale(SignalRescale);

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
            //}else if(varName.BeginsWith("_interf")){varName="_CMS_hzz2l2v"+varName;  //commented out for the HighMass paper
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
      if(!( procName.Contains("t#bar{t}") || procName.Contains("Single top") || procName.Contains("Top") || procName.Contains("WWW") || procName.Contains("WW") || procName.Contains("WW#rightarrow 2l2#nu") ||  procName.Contains("WW#rightarrow lnu2q") || procName.Contains("W#rightarrow l#nu") || procName.Contains("W,multijets") || procName.Contains("Z#rightarrow #tau#tau") || procName.Contains("ZZ#rightarrow Z#tau#tau") ||procName.Contains("Top/W/WW") ))continue;
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
		 	int bin = 6; //bin = 5 => AllSide Region; bin = 6 => UpSide Region
      if(hCtrl_SB->GetBinContent(bin)>0){
        alpha     = hChan_SB->GetBinContent(bin) / hCtrl_SB->GetBinContent(bin);
        alpha_err = ( fabs( hChan_SB->GetBinContent(bin) * hCtrl_SB->GetBinError(bin) ) + fabs(hChan_SB->GetBinError(bin) * hCtrl_SB->GetBinContent(bin) )  ) / pow(hCtrl_SB->GetBinContent(bin), 2);        
      }
			//                 if(chData->second.channel.find("ee"  )==0){alphaUsed = 0.44; alphaUsed_err=0.03;}
			//                 if(chData->second.channel.find("mumu")==0){alphaUsed = 0.71; alphaUsed_err=0.04;}
			//                 if(chData->second.channel.find("ee"  )==0){alphaUsed = 0.47; alphaUsed_err=0.03;} //25/01/2014
			//                 if(chData->second.channel.find("mumu")==0){alphaUsed = 0.61; alphaUsed_err=0.04;}
      //if(chData->second.channel.find("ee"  )==0){alphaUsed = 0.36; alphaUsed_err=0.02;} //26/01/2016
      //if(chData->second.channel.find("mumu")==0){alphaUsed = 0.77; alphaUsed_err=0.04;}

      //if(chData->second.channel.find("ee"  )==0){alphaUsed = 0.384583; alphaUsed_err = 0.00600805;} //06/04/2017
      //if(chData->second.channel.find("mumu")==0){alphaUsed = 0.674941; alphaUsed_err = 0.00874789;}

      //if(chData->second.channel.find("ee"  )==0){alphaUsed = 0.375; alphaUsed_err = 0.006;} //09/05/2017
      //if(chData->second.channel.find("mumu")==0){alphaUsed = 0.684; alphaUsed_err = 0.005;}

      if(chData->second.channel.find("ee"  )==0){alphaUsed = 0.369; alphaUsed_err = 0.006;} //09/05/2017
      if(chData->second.channel.find("mumu")==0){alphaUsed = 0.683; alphaUsed_err = 0.0095;}
      double valDD, valDD_err;
      double valMC, valMC_err;
      valMC = hNRB->IntegralAndError(1,hNRB->GetXaxis()->GetNbins(),valMC_err);  if(valMC<1E-6){valMC=0.0; valMC_err=0.0;}

      if(hCtrl_SI->Integral(1, hCtrl_SI->GetXaxis()->GetNbins()+1)<=0){ //if no data in emu: take the shape from MC and fix integral of upper stat uncertainty to 1.8events
        double ErrInt = 0;
        for(int bi=1;bi<=hNRB->GetXaxis()->GetNbins()+1;bi++){                       
          double val = hNRB->GetBinContent(bi);
          double err = hNRB->GetBinError(bi);
          double ratio = 1.8/valMC;
          double newval = 1E-7;
          double newerr = sqrt(pow(val*ratio,2) + pow(err*ratio,2));
          hNRB->SetBinContent(bi, newval );
          hNRB->SetBinError  (bi, newerr );
          ErrInt += newerr;
        }
        printf("err Int = %f\n", ErrInt);
      }else if(chData->second.bin.find("vbf")==0){                   //for VBF stat in emu is too low, so take the shape from MC and scale it to the expected yield
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
			//                 hNRB->SetFillStyle(1001);
			//                 hNRB->SetFillColor(592);

      //save values for printout
      valDD = hNRB->IntegralAndError(1,hNRB->GetXaxis()->GetNbins()+1,valDD_err); if(valDD<1E-6){valDD=0.0; valDD_err=0.0;}

      //remove all syst uncertainty
      chNRB->second.shapes[mainHisto.Data()].clearSyst();
      //add syst uncertainty                 
      chNRB->second.shapes[mainHisto.Data()].uncScale[string("CMS_hzz2l2v_sys_topwww") + systpostfix.Data()] =valDD>=1E-4?valDD*NonResonnantSyst:1.8*valDD;

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
  // properly assign syst uncertainty to the instr Met Background 
  //
  void AllInfo_t::addInstrMetSyst(FILE* pFile, std::vector<TString>& selCh,TString ctrlCh,TString mainHisto){
    //check that the z+jets proc exist
    std::map<string, ProcessInfo_t>::iterator instrMetIt=procs.find("");
    for(std::map<string, ProcessInfo_t>::iterator it=procs.begin(); it!=procs.end();it++){
      if(!it->second.isBckg || it->second.isData)continue;
      TString procName = it->first.c_str();
      if(!( procName.Contains("Instr.") ))continue;

      it->second.isData = true;
      it->second.isSign = false;
      it->second.isBckg = true;


      for(std::map<string, ChannelInfo_t>::iterator chMC = it->second.channels.begin(); chMC!=it->second.channels.end(); chMC++){
        if(std::find(selCh.begin(), selCh.end(), chMC->second.channel)==selCh.end())continue;

        TH1* histo = chMC ->second.shapes[mainHisto.Data()].histo();
        double val_err; double val = histo->IntegralAndError(1,histo->GetXaxis()->GetNbins()+1,val_err); if(val<1E-6){val=0.0; val_err=0.0;}

        //remove all syst uncertainty
        chMC->second.shapes[mainHisto.Data()].clearSyst();
        //add syst uncertainty                 
        chMC->second.shapes[mainHisto.Data()].uncScale[string("_CMS_hzz2l2v_sys_zll") + systpostfix.Data()] = val*GammaJetSyst;
      }
    }
  }

  //
  // properly assign syst uncertainty to the instr Met Background (2017 method) 
  //
  void AllInfo_t::addInstrMetSyst_2017(std::vector<TString>& selCh,TString mainHisto){
		//check that the z+jets proc exist
    std::map<string, ProcessInfo_t>::iterator instrMetIt=procs.find("");
    for(std::map<string, ProcessInfo_t>::iterator it=procs.begin(); it!=procs.end();it++){
      if(!it->second.isBckg || it->second.isData)continue;
      TString procName = it->first.c_str();
      if(!( procName.Contains("Instr.") ))continue;

      it->second.isData = true;
      it->second.isSign = false;
      it->second.isBckg = true;

			TString InstrMET_allExceptGammaStats_Url(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/InstrMET_systematics/InstrMET_systematics_ALL_EXCEPT_GAMMASTATS.root");
			TFile* f_InstrMET_allExceptGammaStats = TFile::Open(InstrMET_allExceptGammaStats_Url);
			if(!f_InstrMET_allExceptGammaStats ){
			  std::cout<< "Missing InstrMET syst files! No syst for InstrMET!'" << std::endl;
				continue; 
			}


      for(std::map<string, ChannelInfo_t>::iterator chMC = it->second.channels.begin(); chMC!=it->second.channels.end(); chMC++){
        if(std::find(selCh.begin(), selCh.end(), chMC->second.channel)==selCh.end())continue;
				TH1* h_InstrMET_Up_allExceptGammaStats = (TH1*)utils::root::GetObjectFromPath(f_InstrMET_allExceptGammaStats, (chMC->second.channel+chMC->second.bin +"_mt_InstrMET_absolute_shape_up").c_str() );
        TH1* h_InstrMET_Down_allExceptGammaStats = (TH1*)utils::root::GetObjectFromPath(f_InstrMET_allExceptGammaStats, (chMC->second.channel+chMC->second.bin +"_mt_InstrMET_absolute_shape_down").c_str() );
        if(!h_InstrMET_Up_allExceptGammaStats || !h_InstrMET_Down_allExceptGammaStats){continue;}
  			//Careful about (re)binning...

        TH1* histo = chMC ->second.shapes[mainHisto.Data()].histo();

        //remove all syst uncertainty
        chMC->second.shapes[mainHisto.Data()].clearSyst();
        //add syst uncertainty                 
        chMC->second.shapes[mainHisto.Data()].uncShape[string("_CMS_hzz2l2v_sys_"+chMC->second.bin+"_zll") + it->second.shortName+systpostfix.Data()+"Up"] = h_InstrMET_Up_allExceptGammaStats;
        chMC->second.shapes[mainHisto.Data()].uncShape[string("_CMS_hzz2l2v_sys_"+chMC->second.bin+"_zll") + it->second.shortName+systpostfix.Data()+"Down"] =h_InstrMET_Down_allExceptGammaStats; 
     }
    }
    //Recompute the total background with correct uncertainties
    computeTotalBackground();
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
      double cutMin=shapeMin; double cutMax=shapeMax;
      if((shapeMinVBF!=shapeMin || shapeMaxVBF!=shapeMax) && chMC->second.bin.find("vbf")!=string::npos){cutMin=shapeMinVBF; cutMax=shapeMaxVBF;}

			//              int indexcut_ = indexcut;
			//              if(indexvbf>=0 && chMC->second.bin.find("vbf")!=string::npos){indexcut_ = indexvbf;}
      int indexcut_ = indexcutM[chMC->second.bin];

      TH2* gjets2Dshape = NULL;
      if(mainHisto==histo && histoVBF!="" && chMC->second.bin.find("vbf")!=string::npos){
        gjets2Dshape  = (TH2*)pdir->Get(((chMC->second.channel+chMC->second.bin+"_")+histoVBF.Data()).c_str());
      }else{
        gjets2Dshape  = (TH2*)pdir->Get(((chMC->second.channel+chMC->second.bin+"_")+mainHisto.Data()).c_str());
      }
      if(!gjets2Dshape)printf("Can't find histo: %s in g+jets template\n",((chMC->second.channel+chMC->second.bin+"_")+mainHisto.Data()).c_str());

      TH1* hMC = chMC->second.shapes[mainHisto.Data()].histo();
      TH1* hDD = gjets2Dshape->ProjectionY("tmpName",indexcut_,indexcut_);
      filterBinContent(hDD);
      utils::root::fixExtremities(hDD, false, true);
      if(!(mainHisto==histo && histoVBF!="" && chMC->second.bin.find("vbf")!=string::npos)){
        for(int x=0;x<=hDD->GetXaxis()->GetNbins()+1;x++){
          if(hDD->GetXaxis()->GetBinCenter(x)<=cutMin || hDD->GetXaxis()->GetBinCenter(x)>=cutMax){hDD->SetBinContent(x,0); hDD->SetBinError(x,0);}
        }
      }

      //Check the binning!!!
      if(hDD->GetXaxis()->GetXmin()!=hMC->GetXaxis()->GetXmin()){printf("Gamma+Jet templates have a different XAxis range\nStop the script here\n"); exit(0);}
      if(hDD->GetXaxis()->GetBinWidth(1)!=hMC->GetXaxis()->GetBinWidth(1)){
        double dywidth = hDD->GetXaxis()->GetBinWidth(1);
        printf("Gamma+Jet templates have a different bin width in %s channel:", chMC->first.c_str());
        double mcwidth = hMC->GetXaxis()->GetBinWidth(1);
        if(dywidth>mcwidth){
          printf("bin width in Gamma+Jet templates is larger than in MC samples (%f vs %f) --> can not rebin!\nStop the script here\n", dywidth,mcwidth); 
          exit(0);
        }else{
          int rebinfactor = (int)(mcwidth/dywidth);
          if(((int)mcwidth)%((int)dywidth)!=0){printf("bin width in Gamma+Jet templates are not multiple of the mc histograms bin width\n"); exit(0);}
          printf("Rebinning by %i --> ", rebinfactor);
          hDD->Rebin(rebinfactor);
          printf("Binning DataDriven ZJets Min=%7.2f  Max=%7.2f Width=%7.2f compared to MC ZJets Min=%7.2f  Max=%7.2f Width=%7.2f\n", hDD->GetXaxis()->GetXmin(), hDD->GetXaxis()->GetXmax(), hDD->GetXaxis()->GetBinWidth(1), hMC->GetXaxis()->GetXmin(), hMC->GetXaxis()->GetXmax(), hMC->GetXaxis()->GetBinWidth(1));
        }
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
      //chData->second.shapes[mainHisto.Data()].uncScale[string("CMS_hzz2l2v_sys_zll") + systpostfix.Data()] = valDD*GammaJetSyst;
 			//Hugo: We don't use this way anymore

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
  // replace MC Backgrounds with FakeLeptons by DataDriven estimate
  //
  void AllInfo_t::doFakeLeptonEstimation(FILE* pFile, std::vector<TString>& selCh,TString ctrlCh,TString mainHisto, bool isCutAndCount){
    TString DYProcName = "Z#rightarrow ll";
		//           TString GammaJetProcName = "Instr. background (data)";
    std::map<TString, double> LowMetIntegral;
    std::vector<string> lineprintouts;

    //open gamma+jet file
    TFile* inF = TFile::Open(FREFile);
    if( !inF || inF->IsZombie() ){ cout << "Invalid file name : " << FREFile << endl; return; }           
    TDirectory* pdir = (TDirectory *)inF;
		//           TDirectory* pdir = (TDirectory *)inF->Get(GammaJetProcName);
		//           if(!pdir){ printf("Skip Z+Jet estimation because %s directory is missing in Gamma+Jet file\n", GammaJetProcName.Data()); return;}
    gROOT->cd(); //make sure that all histograms that will be created will be in memory and not in file


    //check that the data proc exist
    std::map<string, ProcessInfo_t>::iterator dataProcIt=procs.find("data");             
    if(dataProcIt==procs.end()){printf("The process 'data' was not found... can not do non-resonnant background prediction\n"); return;}

    //create a new proc for Z+Jets datadriven backgrounds as a copy of the MC one
    TString DDProcName = "FakeLep";
    for(std::vector<string>::iterator p=sorted_procs.begin(); p!=sorted_procs.end();p++){if((*p)==DDProcName.Data()){sorted_procs.erase(p);break;}}           
    sorted_procs.push_back(DDProcName.Data());
    procs[DDProcName.Data()] = ProcessInfo_t(); //reset
    ProcessInfo_t& procInfo_DD = procs[DDProcName.Data()];
    procInfo_DD.shortName = "fakelep";
    procInfo_DD.isData = true;
    procInfo_DD.isSign = false;
    procInfo_DD.isBckg = true;
    procInfo_DD.xsec   = 0.0;
    procInfo_DD.br     = 1.0;

    //create an histogram containing all the MC backgrounds
    std::vector<string> toBeDelete;
    for(std::map<string, ProcessInfo_t>::iterator it=procs.begin(); it!=procs.end();it++){
      if(!it->second.isBckg || it->second.isData)continue;
      TString procName = it->first.c_str();
      if(!( procName.Contains("WZ") || procName.Contains("WW") || procName.Contains("V#gamma") || procName.Contains("tT,tTV,t,T") || procName.Contains("QCD") ||  procName.Contains("W+jets") ||  procName.Contains("Z#rightarrow ll") ) )continue;
      if(procName.Contains("ZZ#rightarrow ll#tau#tau"))continue; //do not supress ZZ to 2l2tau
      addProc(procInfo_DD, it->second);
      for(std::vector<string>::iterator p=sorted_procs.begin(); p!=sorted_procs.end();p++){if((*p)==it->first){sorted_procs.erase(p);break;}}
      toBeDelete.push_back(it->first);
    }
    for(std::vector<string>::iterator p=toBeDelete.begin();p!=toBeDelete.end();p++){procs.erase(procs.find((*p)));}


    for(std::map<string, ChannelInfo_t>::iterator chData = dataProcIt->second.channels.begin(); chData!=dataProcIt->second.channels.end(); chData++){            
      if(std::find(selCh.begin(), selCh.end(), chData->second.channel)==selCh.end())continue;

      std::map<string, ChannelInfo_t>::iterator chDD  = procInfo_DD.channels.find(chData->first);  
      if(chDD==procInfo_DD.channels.end()){  //this channel does not exist, create it
        procInfo_DD.channels[chData->first] = ChannelInfo_t();     
        chDD                = procInfo_DD.channels.find(chData->first);
        chDD->second.bin     = chData->second.bin;
        chDD->second.channel = chData->second.channel;
      }

      //load template data
      double cutMin=shapeMin; double cutMax=shapeMax;
      if((shapeMinVBF!=shapeMin || shapeMaxVBF!=shapeMax) && chData->second.bin.find("vbf")!=string::npos){cutMin=shapeMinVBF; cutMax=shapeMaxVBF;}

			//              int indexcut_ = indexcut;
			//              if(indexvbf>=0 && chData->second.bin.find("vbf")!=string::npos){indexcut_ = indexvbf;}
      int indexcut_ = indexcutM[chData->second.bin];

      TH2* h2Dshape = NULL;
      if(mainHisto==histo && histoVBF!="" && chData->second.bin.find("vbf")!=string::npos){
        h2Dshape  = (TH2*)pdir->Get(((chData->second.channel+chData->second.bin+"_")+histoVBF.Data()).c_str());
      }else{
        h2Dshape  = (TH2*)pdir->Get(((chData->second.channel+chData->second.bin+"_")+mainHisto.Data()).c_str());
      }
      if(!h2Dshape)printf("Can't find histo: %s in fake rate estimate template\n",((chData->second.channel+chData->second.bin+"_")+mainHisto.Data()).c_str());

      TH1* hMC = chData->second.shapes[mainHisto.Data()].histo();
      TH1* hDD = h2Dshape->ProjectionY("tmpName",indexcut_,indexcut_);
      double OSIntegral = hDD->GetBinContent(0); double OSIntegralError =  hDD->GetBinError(0);   hDD->SetBinContent(0, 0.0);  hDD->SetBinError(0, 0.0);  //get values from OS only
      filterBinContent(hDD);
      utils::root::fixExtremities(hDD, false, true);
      if(!(mainHisto==histo && histoVBF!="" && chData->second.bin.find("vbf")!=string::npos)){
        for(int x=0;x<=hDD->GetXaxis()->GetNbins()+1;x++){
          if(hDD->GetXaxis()->GetBinCenter(x)<=cutMin || hDD->GetXaxis()->GetBinCenter(x)>=cutMax){hDD->SetBinContent(x,0); hDD->SetBinError(x,0);}
					//                     if(hDD->GetBinContent(x)<0){hDD->SetBinContent(x,0); hDD->SetBinError(x,0);} //make sure that all bins have positive content  //DONT DO THIS AT THIS STEP, we have a dedicated function that check for negative bins
        }
      }

      //Check the binning!!!
      if(hDD->GetXaxis()->GetXmin()!=hMC->GetXaxis()->GetXmin()){printf("fake rate templates have a different XAxis range\nStop the script here\n"); exit(0);}
      if(hDD->GetXaxis()->GetBinWidth(1)!=hMC->GetXaxis()->GetBinWidth(1)){
        double dywidth = hDD->GetXaxis()->GetBinWidth(1);
        printf("fake rate templates have a different bin width in %s channel:", chData->first.c_str());
        double mcwidth = hMC->GetXaxis()->GetBinWidth(1);
        if(dywidth>mcwidth){
          printf("bin width in fake rate templates is larger than in MC samples (%f vs %f) --> can not rebin!\nStop the script here\n", dywidth,mcwidth); 
          exit(0);
        }else{
          int rebinfactor = (int)(mcwidth/dywidth);
          if(((int)mcwidth)%((int)dywidth)!=0){printf("bin width in fake rate templates are not multiple of the mc histograms bin width\n"); exit(0);}
          printf("Rebinning by %i --> ", rebinfactor);
          hDD->Rebin(rebinfactor);
          printf("Binning DataDriven fake rate Min=%7.2f  Max=%7.2f Width=%7.2f compared to MC ZJets Min=%7.2f  Max=%7.2f Width=%7.2f\n", hDD->GetXaxis()->GetXmin(), hDD->GetXaxis()->GetXmax(), hDD->GetXaxis()->GetBinWidth(1), hMC->GetXaxis()->GetXmin(), hMC->GetXaxis()->GetXmax(), hMC->GetXaxis()->GetBinWidth(1));
        }
      }

      if(isCutAndCount){  //This is ugly, but it's the best thing I've found out right now.  The integral of this thing should be the one from OS only
        hDD->SetBinContent(0, 0.0); hDD->SetBinError(0, 0.0);
        hDD->SetBinContent(hDD->GetXaxis()->GetNbins()+1, 0.0); hDD->SetBinError(hDD->GetXaxis()->GetNbins()+1, 0.0);
        for(int x=1;x<=hDD->GetXaxis()->GetNbins();x++){
          hDD->SetBinContent(x, OSIntegral/hDD->GetXaxis()->GetNbins()); hDD->SetBinError(x, OSIntegralError/sqrt(hDD->GetXaxis()->GetNbins()));
        }
      }

      //save histogram to the structure
      hDD->Scale(DDRescale);
	    //for(int x=0;x<=hDD->GetXaxis()->GetNbins()+1;x++){
	    //      cout << "DD: " << hDD->GetBinContent(x) << " in bin x " << x << endl; 
	    //      cout << "MC: " << hMC->GetBinContent(x) << " in bin x " << x << endl; 
      //      if(hDD->GetBinContent(x)==0 && hMC->GetBinContent(x)==0){
			//	      hDD->SetBinContent(x,2E-6);
			//	      hMC->SetBinContent(x,2E-6);
      //      }
		  //if(hMC->GetBinContent(x)<=1E-6) hMC->SetBinContent(x,2E-6);
	    //}
      chDD->second.shapes[mainHisto.Data()].histo()->Reset();
      chDD->second.shapes[mainHisto.Data()].histo()->Add(hDD);
      //printouts
      char printout[2048];
      double valMC_err, valMC = hMC->IntegralAndError(1,hMC->GetXaxis()->GetNbins()+1,valMC_err); if(fabs(valMC)<1E-6){valMC=0.0; valMC_err=0.0;}
      double valDD_err, valDD = hDD->IntegralAndError(1,hDD->GetXaxis()->GetNbins()+1,valDD_err); if(fabs(valDD)<1E-6){valDD=0.0; valDD_err=0.0;}
      sprintf(printout,"%20s & %30s & %30s\\\\", chData->first.c_str(), utils::toLatexRounded(valDD,valDD_err,valDD*GammaJetSyst).c_str(), utils::toLatexRounded(valMC,valMC_err).c_str() );
      lineprintouts.push_back(printout);

      //add syst uncertainty
      chDD->second.shapes[mainHisto.Data()].clearSyst();
      chDD->second.shapes[mainHisto.Data()].uncScale[string("CMS_hzz2l2v_sys_DD") + systpostfix.Data()] = valDD!=0?valDD*FakeLeptonDDSyst:FakeLeptonDDSyst;

      //clean
      delete hDD;
      delete h2Dshape;
    }

    //all done with template file
    inF->Close();

    //recompute total background
    computeTotalBackground();

    //printouts 
    if(pFile){
      fprintf(pFile,"\\begin{table}[htp]\n\\begin{center}\n\\caption{fake lepton background estimation.}\n\\label{tab:table}\n");
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
          if(!histo)continue;
	    TString jetBin = ch->second.bin.c_str();

          if(jetBin.Contains("vbf")){
            double xbins[] = {150, 225, 300, 375, 450, 600, 750, 1100, 3000};
            int nbins=sizeof(xbins)/sizeof(double);
            unc->second = histo->Rebin(nbins-1, histo->GetName(), (double*)xbins);
            utils::root::fixExtremities(unc->second, false, true);
          }else if( jetBin.Contains("eq0jets") || jetBin.Contains("geq1jets") ){
            double xbins[] = {150, 225, 300, 375, 450, 525, 600, 725, 850, 975, 1100, 1350, 1600, 2100, 3000};
            int nbins=sizeof(xbins)/sizeof(double);
            unc->second = histo->Rebin(nbins-1, histo->GetName(), (double*)xbins);
            utils::root::fixExtremities(unc->second, false, true);
          }

		  		//Old Binning
          //double xbins[] = {150, 300, 450, 600, 850, 1100, 1600, 2100, 3000}; 
          //double xbins[] = {150, 225, 300, 375, 450, 525, 600, 725, 850, 975, 1100, 1350, 1600, 1850, 2100, 2600, 3000}; 
          //int nbins=sizeof(xbins)/sizeof(double);
          //unc->second = histo->Rebin(nbins-1, histo->GetName(), (double*)xbins);
          //utils::root::fixExtremities(unc->second, true, true);
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
      TString signProcName =  procLR->first;
      signProcName.ReplaceAll("()","(");              
      signProcName+=(int)mass; signProcName+=")";
      printf("Interpolate %s  based on %s and %s\n", signProcName.Data(), procLR->second.first.c_str(), procLR->second.second.c_str());

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
      proc.xsec      = procL.xsec + (Ratio * (procR.xsec - procL.xsec)); 
      proc.br        = procL.br   + (Ratio * (procR.br   - procL.br));

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

        //morphing ignore stat unc, so add it as an uncertainty to morph it separately
        shapeInfo .makeStatUnc("", "", "", true );
        shapeInfoL.makeStatUnc("", "", "", true );
        shapeInfoR.makeStatUnc("", "", "", true );

        for(std::map<string, TH1*  >::iterator unc=shapeInfoL.uncShape.begin();unc!=shapeInfoL.uncShape.end();unc++){
          if(shapeInfo.uncShape.find(unc->first)==shapeInfo.uncShape.end())shapeInfo.uncShape[unc->first] = (TH1*)shapeInfoL.uncShape[unc->first]->Clone(signProcName+ch->first+unc->first+"tmp");
          TH1D* hL = (TH1D*)shapeInfoL.uncShape[unc->first];
          TH1D* hR = (TH1D*)shapeInfoR.uncShape[unc->first];
          TH1D* h  = (TH1D*)shapeInfo .uncShape[unc->first];
          if(!hL || !hR || !h)continue;
          hL->Scale(1.0/(procL.xsec*procL.br));
          hR->Scale(1.0/(procR.xsec*procR.br));
          h->Reset();
          if(hL->Integral()>0 && hR->Integral()>0){//interpolate only if the histograms are not null
            h->Add(th1fmorph(signProcName+ch->first+unc->first,signProcName+ch->first+unc->first, hL, hR, procL.mass, procR.mass, proc.mass, (1-Ratio)*hL->Integral() + Ratio*hR->Integral(), 0), 1);
          }
          //printf("EFF : syst%25s %f - %f -%f\n", unc->first.c_str(), hL->Integral(), h->Integral(), hR->Integral());
          h->Scale(proc.xsec*proc.br);
        }

        //incorporate back stat uncertainty into the main histo, and delete it from the list of uncertainty (will be readded later)
        TH1D* h = (TH1D*)shapeInfo.uncShape[""];  TH1D* hUp = (TH1D*)shapeInfo.uncShape["statUp"];

        if(h && hUp){
          for(unsigned int x=0;x<h->GetNbinsX();x++){
            h->SetBinError(x, fabs(hUp->GetBinContent(x)-h->GetBinContent(x)));
          }
        }
        shapeInfo.removeStatUnc();
      }

      //erase sideband procs
      for(std::vector<string>::iterator p=sorted_procs.begin(); p!=sorted_procs.end();p++){if((*p)==procLR->second.first ){sorted_procs.erase(p);break;}}              
      for(std::vector<string>::iterator p=sorted_procs.begin(); p!=sorted_procs.end();p++){if((*p)==procLR->second.second){sorted_procs.erase(p);break;}}       
      procs.erase(procs.find(procLR->second.first ));
      procs.erase(procs.find(procLR->second.second ));
    }
  }


  //
  // Rescale VBF events to match SM ggF/VBF production 
  //
  void AllInfo_t::scaleVBF(string histoName){
    if(mass<=0)return;

    //Loop on processes and channels
    for(unsigned int p=0;p<sorted_procs.size();p++){
      string procName = sorted_procs[p];
      std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
      if(it==procs.end() || !it->second.isSign)continue;              
      if(procName.find("qqH")==string::npos)continue;

      double scale = Hxswg::utils::getVBFoverGGF(systpostfix.Data())->Eval(mass);
      it->second.xsec*=scale;

      printf("scale VBF signal by %f\n", scale);
      for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
        if(ch->second.shapes.find(histoName)==(ch->second.shapes).end())continue;
        ShapeData_t& shapeInfo = ch->second.shapes[histoName];
        for(std::map<string, TH1*  >::iterator unc=shapeInfo.uncShape.begin();unc!=shapeInfo.uncShape.end();unc++){
          if(unc->second)unc->second->Scale(scale);
        }
      }
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

      double cprime=1.0; double  brnew=0.0;
      if(signalSufix!="" && signalSufix.Contains("_cpsq")){sscanf(signalSufix.Data(), "_cpsq%lf_brn%lf", &cprime, &brnew); cprime=sqrt(cprime);}
      else if(signalSufix!="" && signalSufix.Contains("_cp  ")){sscanf(signalSufix.Data(), "_cp%lf_brn%lf", &cprime, &brnew); }

      if(doInterf && it->second.mass>400){// && it->first.find("ggH")!=string::npos){
				//                 value used at run1
				//                 sF   = 0.897-0.000152*it->second.mass+7.69e-07*pow(it->second.mass,2);
				//                 sFDn = 0.907-2.08e-05*it->second.mass+4.63e-07*pow(it->second.mass,2);
				//                 sFUp = 0.889-0.000357*it->second.mass+1.21e-06*pow(it->second.mass,2);

        //value for run2
        sF   = 1.316-0.000068*it->second.mass+3.32e-06*pow(it->second.mass,2);
        sFDn = 1.590-0.000193*it->second.mass+2.26e-06*pow(it->second.mass,2);
        sFUp = 0.906+0.000162*it->second.mass+4.35e-06*pow(it->second.mass,2);


        if(it->second.mass>=400 && signalSufix!=""){ //scale factor for Narrow Resonnance
					//                    sF=1 + (sF-1)/pow(cprime,2);   sFDn = 1.0;  sFUp = 1 + (sF-1)*2;        //100% Uncertainty
          sF=1 + (sF-1)/(cprime* sqrt(1-brnew));   sFDn = 1.0;  sFUp = 1 + (sF-1)*2;        //100% Uncertainty
          if(sF<1){sF=1.0;  sFUp = 1 + (sF-1)*2;}
          printf("Scale Factor for Narrow Resonnance : %f [%f,%f] applied on %s\n", sF, sFDn, sFUp, it->first.c_str());                  
        }else{
          printf("Scale Factor for Interference : %f [%f,%f] applied on %s\n",sF, sFDn, sFUp, it->first.c_str());
        }
        if(sFDn>sFUp){double tmp = sFUp; sFUp = sFDn; sFDn = tmp;}
      }
      printf("Total Scale Factor : %f [%f,%f] applied on %s\n",sF, sFDn, sFUp, it->first.c_str());

      //rescale xsection for Narrow Resonnances
      if(signalSufix!=""){
        double xsecScale = pow(cprime,2) * (1-brnew);
        printf("Scale Factor to the narrow resonnance xsection is : %f\n", xsecScale);
        sF *= xsecScale;  sFUp *= xsecScale;  sFDn *= xsecScale;
      }
      it->second.xsec *= sF;
      fflush(stdout);

      for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
        if(ch->second.shapes.find(histoName)==(ch->second.shapes).end())continue;
        ShapeData_t& shapeInfo = ch->second.shapes[histoName];
        for(std::map<string, TH1*  >::iterator unc=shapeInfo.uncShape.begin();unc!=shapeInfo.uncShape.end();unc++){
          unc->second->Scale(sF);
        }
        TH1* histo = (TH1*)shapeInfo.histo();
        if(!histo){printf("Histo does not exit... skip it \n"); fflush(stdout); continue;}
        if(sF<=0){printf("sF has weird values : %f, set it back to one\n", sF); sF=1.0; fflush(stdout);}
        TH1* tmp;
        tmp = (TH1*)histo->Clone(TString("interf_ggH_") + histo->GetName() + "Down"); tmp->Scale(sFDn/sF); shapeInfo.uncShape[string("_CMS_hzz2l2v_interf_") + it->second.shortName+"Down"] = tmp;
        tmp = (TH1*)histo->Clone(TString("interf_ggH_") + histo->GetName() + "Up"  ); tmp->Scale(sFUp/sF); shapeInfo.uncShape[string("_CMS_hzz2l2v_interf_") + it->second.shortName+"Up"  ] = tmp;
      }
      }
    }

    //
    // merge histograms from different bins together... but keep the channel separated 
    //
    void AllInfo_t::mergeBins(std::vector<string>& binsToMerge, string NewName){
      printf("Merge the following bins of the same channel together: "); for(unsigned int i=0;i<binsToMerge.size();i++){printf("%s ", binsToMerge[i].c_str());}
      printf("The resulting bin will be called %s\n", NewName.c_str());
      for(unsigned int p=0;p<sorted_procs.size();p++){
        string procName = sorted_procs[p];
        std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
        if(it==procs.end())continue;

        for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
          if(find(binsToMerge.begin(), binsToMerge.end(), ch ->second.bin)==binsToMerge.end())continue;  //make sure this bin should be merged
          for(std::map<string, ChannelInfo_t>::iterator ch2 = ch; ch2!=it->second.channels.end(); ch2++){
            if(ch->second.channel != ch2->second.channel)continue; //make sure we merge bin in the same channel
            if(ch->second.bin     == ch2->second.bin    )continue; //make sure we do not merge with itself
            if(find(binsToMerge.begin(), binsToMerge.end(), ch2->second.bin)==binsToMerge.end())continue;  //make sure this bin should be merged
            addChannel(ch->second, ch2->second); //FIXME this only adds the nominal shapes, not also the syst
            it->second.channels.erase(ch2);  
            ch2=ch;
          }
          ch->second.bin = NewName;
        }

        //also update the map keys
        std::map<string, ChannelInfo_t> newMap;
        for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
          newMap[ch->second.channel+ch->second.bin] = ch->second;
        }
        it->second.channels = newMap;
      }
    }


    void AllInfo_t::HandleEmptyBins(string histoName){
      for(unsigned int p=0;p<sorted_procs.size();p++){
        string procName = sorted_procs[p];
        //if(procName!="FakeLep")continue; //only do this for the FakeLepbackground right now
        std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
        if(it==procs.end())continue;
        if(it->second.isData && procName!="Instr. MET")continue; //only do this for MC or InstrMET (which also contains MC and negative bins)
        for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
          ShapeData_t& shapeInfo = ch->second.shapes[histoName];
          TH1* histo = (TH1*)shapeInfo.histo();
          if(!histo){printf("Histo does not exit... skip it \n"); fflush(stdout); continue;}

          double StartIntegral = histo->Integral();
          for(int binx=1;binx<=histo->GetNbinsX();binx++){
            if(histo->GetBinContent(binx)<=0){histo->SetBinContent(binx, 1E-6); histo->SetBinError(binx, 1E-6);  }
          }
          double EndIntegral = histo->Integral();                 
          shapeInfo.rescaleScaleUncertainties(StartIntegral, EndIntegral);


          for(std::map<string, TH1*  >::iterator unc=shapeInfo.uncShape.begin();unc!=shapeInfo.uncShape.end();unc++){
            for(int binx=1;binx<=unc->second->GetNbinsX();binx++){
              if(unc->second->GetBinContent(binx)<=0){unc->second->SetBinContent(binx, 1E-6); }; //histo->SetBinError(binx, 1.8);  }
          }
        }
      }
    }

    //recompute total background
    computeTotalBackground();
  }





