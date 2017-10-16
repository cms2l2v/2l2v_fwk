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
    return true;
  } 
}

/*****************************************************************/
//double FRWeights::getWeight(string cat,string bin,string var,string wrt,double pT)
double FRWeights::getWeight(const std::string& cat ,const std::string& bin, const std::string& var,const std::string& wrt ,const double& pT)
/*****************************************************************/
{
  TGraphErrors* graph =  (TGraphErrors*)WeightsFile->Get((cat+"FRWeights"+var+bin+wrt).c_str());
  double result=-1.;
  if(graph){
    result=(1 - graph->Eval(pT));
  } else {
    result = 1;
  }
  
  return result;
}
