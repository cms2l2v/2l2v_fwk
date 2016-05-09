#include "UserCode/llvv_fwk/interface/FRWeights.h"
#include <TFile.h>
#include <TSystem.h>
#include <math.h>
#include <vector>
#include <iostream>

using namespace std;

/*****************************************************************/
FRWeights::FRWeights()
/*****************************************************************/
{
}

/*****************************************************************/
FRWeights::~FRWeights()
/*****************************************************************/
{
}

/*****************************************************************/
bool FRWeights::init(const string& WeightsFileName)
/*****************************************************************/
{
  WeightsFile = TFile::Open(WeightsFileName.c_str());
  if(!WeightsFile) {
    cout<<"ERROR: Cannot open weights file "<<WeightsFileName<<"\n";
    return false;
  } else {

    std::vector<string> catL  = {"Fake e", "Fake #mu", "Fake #tau_{had}"};
    std::vector<string> binL  = {"Inc.", "Barrel", "Endcap", "Inc. (m_{T}>30)", "Barrel (m_{T}>30)", "Endcap (m_{T}>30)"};
    
    std::vector<string> cat   = {"FR_El", "FR_Mu", "FR_Ta"};
    std::vector<string> bin   = {"","_B","_E", "_TMCut", "_TMCut_B", "_TMCut_E"};
    std::vector<string> var   = {"", "_Id_Iso01weight", "_Id_Iso02weight", "_Id_Iso03weight", "_Id_IsoLoweight", "_Id_IsoMeweight"};
    std::vector<string> wrt   = {"_wrtJetPt", "_wrtLepPt"};
    
    for(unsigned int c=0;c<cat.size();c++){
      for(unsigned int b=0;b<bin.size();b++){
	for(unsigned int v=0;v<var.size();v++){
	  for(unsigned int w=0;w<wrt.size();w++){
	    FRWeightGraphs[cat[c]+var[v]+bin[b]+wrt[w]] = (TGraphErrors*)WeightsFile->Get((cat[c]+"FRWeights"+var[v]+bin[b]+wrt[w]).c_str());   
	  }
	}
      }
    }
    
    WeightsFile->Close();
    
    return true;
  } 
}

/*****************************************************************/
//double FRWeights::getWeight(string cat,string bin,string var,string wrt,double pT)
double FRWeights::getWeight(const std::string& cat ,const std::string& bin, const std::string& var,const std::string& wrt ,const double& pT)
/*****************************************************************/
{
  TGraphErrors* graph = FRWeightGraphs[cat+var+bin+wrt];
  double result=-1.;
  if(graph){
    result=(1 - graph->Eval(pT));
  } else {
    result = 1;
  }
  
  return result;
}
