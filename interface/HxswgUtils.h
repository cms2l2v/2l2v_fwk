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

    bool Equal(double x, double y, double epsilon){ double diff = std::abs(x-y); if(diff < epsilon){ return true;} else { return false; }}

    TGraph* multiplyGraph(TGraph* A, TGraph* B){   double x, y;   for(int i=0;i<A->GetN();i++){A->GetPoint(i, x, y); A->SetPoint(i, x, y*B->Eval(x));}    return A;   }
    TGraph* divideGraph  (TGraph* A, TGraph* B){   double x, y;   for(int i=0;i<A->GetN();i++){A->GetPoint(i, x, y); A->SetPoint(i, x, y/B->Eval(x));}    return A;   }    
    TGraph* makeGraphFromColXandY(std::string dataFile, int colX, int colY);
       
    TGraph* getHWidth() {  return makeGraphFromColXandY(std::string(getenv("CMSSW_BASE")) + "/src/UserCode/llvv_fwk/data/HXSWG/BrHtoGaugeBosons.dat", 0, 16);  }
    TGraph* getHWidthExtended() { TGraph* Width=getHWidth(); int I=Width->GetN();  Width->Set(I+100); for(double mH=1100;mH<3000;mH+=100){ Width->SetPoint(I,mH,mH); I++; } Width->Set(I); return Width;  } //assume WidthH = MH for MH>1TeV
   
    TGraph* getBRHtoZZ(){  return makeGraphFromColXandY(std::string(getenv("CMSSW_BASE")) + "/src/UserCode/llvv_fwk/data/HXSWG/BrHtoGaugeBosons.dat", 0, 13);  }
    TGraph* getBRRadtoZZ(){  return makeGraphFromColXandY(std::string(getenv("CMSSW_BASE")) + "/src/UserCode/llvv_fwk/data/HXSWG/BrRadtoGaugeBosons.dat", 0, 3);  }
    TGraph* getBRRsGravtoZZ(){  return makeGraphFromColXandY(std::string(getenv("CMSSW_BASE")) + "/src/UserCode/llvv_fwk/data/HXSWG/BrRsGravtoGaugeBosons.dat", 0, 2);  }
    TGraph* getBRBulkGravtoZZ(){  return makeGraphFromColXandY(std::string(getenv("CMSSW_BASE")) + "/src/UserCode/llvv_fwk/data/HXSWG/BrBulkGravtoGaugeBosons.dat", 0, 2);  }

    TGraph* get13to8ScaleGGF(){ return makeGraphFromColXandY(std::string(getenv("CMSSW_BASE")) + "/src/UserCode/llvv_fwk/data/HXSWG/pdfRatio8_13.dat", 0, 1);  }
    TGraph* get13to8ScaleVBF(){ return makeGraphFromColXandY(std::string(getenv("CMSSW_BASE")) + "/src/UserCode/llvv_fwk/data/HXSWG/pdfRatio8_13.dat", 0, 2);  }

    TGraph* getGGFXSec7TeV(){  return makeGraphFromColXandY(std::string(getenv("CMSSW_BASE")) + "/src/UserCode/llvv_fwk/data/HXSWG/XsecGGF7.dat", 0, 1);  }
    TGraph* getVBFXSec7TeV(){  return makeGraphFromColXandY(std::string(getenv("CMSSW_BASE")) + "/src/UserCode/llvv_fwk/data/HXSWG/XsecVBF7.dat", 0, 1);  }
    TGraph* getGGFXSec8TeV(){  return makeGraphFromColXandY(std::string(getenv("CMSSW_BASE")) + "/src/UserCode/llvv_fwk/data/HXSWG/XsecGGF8.dat", 0, 1);  }    
    TGraph* getVBFXSec8TeV(){  return makeGraphFromColXandY(std::string(getenv("CMSSW_BASE")) + "/src/UserCode/llvv_fwk/data/HXSWG/XsecVBF8.dat", 0, 1);  }
    TGraph* getGGFXSec13TeV(){  return makeGraphFromColXandY(std::string(getenv("CMSSW_BASE")) + "/src/UserCode/llvv_fwk/data/HXSWG/XsecGGF13.dat", 0, 1);  }
    TGraph* getVBFXSec13TeV(){  return makeGraphFromColXandY(std::string(getenv("CMSSW_BASE")) + "/src/UserCode/llvv_fwk/data/HXSWG/XsecVBF13.dat", 0, 1);  }
//    TGraph* getGGFXSec13TeV(){ TGraph* scale8To13TeV = get13to8ScaleGGF(); TGraph* graph13TeV = multiplyGraph(getGGFXSec8TeV(), scale8To13TeV); delete scale8To13TeV; return graph13TeV;  }
//    TGraph* getVBFXSec13TeV(){ TGraph* scale8To13TeV = get13to8ScaleVBF(); TGraph* graph13TeV = multiplyGraph(getVBFXSec8TeV(), scale8To13TeV); delete scale8To13TeV; return graph13TeV;  }
    TGraph*      getRadXSec13TeV(){  return makeGraphFromColXandY(std::string(getenv("CMSSW_BASE")) + "/src/UserCode/llvv_fwk/data/HXSWG/XsecRad13.dat", 0, 1);  }
    TGraph*   getRsGravXSec13TeV(){  return makeGraphFromColXandY(std::string(getenv("CMSSW_BASE")) + "/src/UserCode/llvv_fwk/data/HXSWG/XsecRsGrav13.dat", 0, 1);  }
    TGraph* getBulkGravXSec13TeV(){  return makeGraphFromColXandY(std::string(getenv("CMSSW_BASE")) + "/src/UserCode/llvv_fwk/data/HXSWG/XsecBulkGrav13.dat", 0, 1);  }
    TGraph* getXSec(std::string Name );

    TGraph* getVBFoverGGF(std::string Name);

    TGraph *getXSecMELA( float cprime);
     
  }  
}
#endif
