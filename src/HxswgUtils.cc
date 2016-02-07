#include "UserCode/llvv_fwk/interface/HxswgUtils.h"
#include "TGraphErrors.h"

namespace Hxswg{

  namespace utils{

     TGraph* makeGraphFromColXandY(std::string dataFile, int colX, int colY){
        FILE* pFile = fopen(dataFile.c_str(), "r");
        if(!pFile){  printf("Couldn't open file %s to read Higgs width values\n", dataFile.c_str()); return NULL; }
       
        TGraph* toReturn = new TGraph(9999);  int N=0;
        char line [4096];
        while(fgets(line, 4096, pFile)){
           if(std::string(line).find("//")==0)continue; //skip line starting by //
           char* pch=strtok(line,"\t"); int Arg=0; double x; double y;
           while (pch!=NULL){ 
              if(Arg==colX){        sscanf(pch, "%lf", &x); 
              }else if(Arg==colY){  sscanf(pch, "%lf", &y);
              }              
              pch=strtok(NULL,"\t");Arg++; 
           }
           toReturn->SetPoint(N, x, y); N++; 
        }fclose(pFile);
        toReturn->Set(N);
        return toReturn;
   }

    TGraph* getXSec(std::string Name){
       if(Name.find("VBF")!=std::string::npos){
          if(Name.find("13TeV")!=std::string::npos){return getVBFXSec13TeV();}
          if(Name.find("8TeV" )!=std::string::npos){return getVBFXSec8TeV();}
          if(Name.find("7TeV" )!=std::string::npos){return getVBFXSec7TeV();}
          return NULL;       
       }else{ //GGF
          if(Name.find("13TeV")!=std::string::npos){return getGGFXSec13TeV();}
          if(Name.find("8TeV" )!=std::string::npos){return getGGFXSec8TeV();}
          if(Name.find("7TeV" )!=std::string::npos){return getGGFXSec7TeV();}
          return NULL;               
       }
       return NULL;
    }

    TGraph* getVBFoverGGF(std::string Name){
       if(Name.find("13TeV")!=std::string::npos){return divideGraph(getVBFXSec13TeV(), getGGFXSec13TeV());}
       if(Name.find("8TeV" )!=std::string::npos){return divideGraph(getVBFXSec8TeV(),  getGGFXSec8TeV() );}
       if(Name.find("7TeV" )!=std::string::npos){return divideGraph(getVBFXSec7TeV(),  getGGFXSec7TeV() );}
       return NULL;
    }




  }
}
