#include "UserCode/llvv_fwk/interface/HxswgUtils.h"
#include "TGraphErrors.h"

namespace Hxswg{

  namespace utils{

     TGraph* getHWidth(){
        std::string dataFile = std::string(getenv("CMSSW_BASE")) + "/src/UserCode/llvv_fwk/data/HXSWG/BrHtoGaugeBosons.dat";
        FILE* pFile = fopen(dataFile.c_str(), "r");
        if(!pFile){  printf("Couldn't open file %s to read Higgs width values\n", dataFile.c_str()); return NULL; }
       
        TGraph* toReturn = new TGraph(9999);  int N=0;

        char line [4096];
        while(fgets(line, 4096, pFile)){
           if(std::string(line).find("//")==0)continue; //skip line starting by //
           char* pch=strtok(line,"\t"); int Arg=0; double mass; double value;
           while (pch!=NULL){ 
              if(Arg== 0){ sscanf(pch, "%lf", &mass); 
              }else{       sscanf(pch, "%lf", &value);
              }
              
              if(Arg==16){ toReturn->SetPoint(N, mass, value); N++; } //width
              pch=strtok(NULL,"\t");Arg++; 
           }
        }fclose(pFile);

        toReturn->Set(N);
        return toReturn;
   }
  }
}
