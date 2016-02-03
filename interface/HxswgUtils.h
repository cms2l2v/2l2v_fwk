#ifndef HxswgUtils_h
#define HxswgUtils_h

#include<iostream>
#include<vector>

#include "TGraph.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TFile.h"
#include "TMath.h"

namespace Hxswg{

  namespace utils{

    //handy function to get a TGraph of the higgs width vs mass
    TGraph* getHWidth();

  }  
}
#endif
