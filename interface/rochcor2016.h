#ifndef ROCHCOR2016_H
#define ROCHCOR2016_H

#include <iostream>
#include <map>
#include "TChain.h"
#include "TClonesArray.h"
#include "TString.h"

#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

#include "RoccoR.h"

class rochcor2016 {
 public:
  rochcor2016();
  rochcor2016(int seed);
  ~rochcor2016();
  
  void momcor_mc(TLorentzVector&, float, int, float&);
  void momcor_data(TLorentzVector&, float, int, float&);
  
  int aetabin(double);
  int etabin(double);
  int phibin(double);
  
 private:
  
  TRandom3 eran;
  TRandom3 sran;
    
  //  static float netabin[9] = {-2.4,-2.1,-1.4,-0.7,0.0,0.7,1.4,2.1,2.4};
  const double pi = TMath::Pi();
  static const double netabin[25];
  static const double anetabin[13];
  
  static constexpr double mu_mass = 0.105658;
  
  static constexpr double mgscl = 1.00128; //global scale in MC 
  static constexpr double dgscl = 1.00245; //global scale in data 
  static constexpr double mgscl_stat = 2.57065e-05; //stat. error of global factor for mass peak in MC 
  static constexpr double mgscl_syst = 3.81651e-03; //syst. error of global factor for mass peak in MC  
  static constexpr double dgscl_stat = 3.34256e-05; //stat. error of global factor for mass peak in data 
  static constexpr double dgscl_syst = 0.0008; //syst. error of global factor for mass peak in data 
  static constexpr double dgscl_iter = 1.00026; //one more iteration to fix offset
  static constexpr double mgscl_iter = 1.00054; //one more iteration to fix offset
  //static const double sf[3];
  //static const double sfer[3];

  //---------------------------------------------------------------------------------------------
  /*
  static const double dcor_m[16][24];  
  static const double dcor_p[16][24];
  static const double mcor_m[16][24];
  static const double mcor_p[16][24];
  static const double dcor_mer[16][24];  
  static const double dcor_per[16][24];
  static const double mcor_mer[16][24];
  static const double mcor_per[16][24];
  */

  static const double dcor_bf[16][24];  
  static const double dcor_ma[16][24];
  static const double mcor_bf[16][24];
  static const double mcor_ma[16][24];
  static const double dcor_bfer[16][24];  
  static const double dcor_maer[16][24];
  static const double mcor_bfer[16][24];
  static const double mcor_maer[16][24];
  
  //=======================================================================================================
  
  static const double dmavg[16][24];  
  static const double dpavg[16][24];  
  static const double mmavg[16][24];  
  static const double mpavg[16][24];
  static const double dmavger[16][24];  
  static const double dpavger[16][24];  
  static const double mmavger[16][24];  
  static const double mpavger[16][24];

  static const double dd[12];
  static const double dder[12];
  static const double md[12];
  static const double mder[12];  

  static const double dscl[24];
  static const double mscl[24];
  //===============================================================================================

  double mptsys_mc_dm[16][24];
  double mptsys_mc_da[16][24];
  double mptsys_da_dm[16][24];
  double mptsys_da_da[16][24];

  double gscler_mc_dev;
  double gscler_da_dev;


};
  
#endif
