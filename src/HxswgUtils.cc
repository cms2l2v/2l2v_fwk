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
       if(Name.find("SM")!=std::string::npos){
          if(Name.find("VBF")!=std::string::npos){
	     //getHWidthExtended gives the BR values for all the mass points above 1000 GeV
             if(Name.find("13TeV")!=std::string::npos){return multiplyGraph( getVBFXSec13TeV(), getBRHtoZZ());}
             if(Name.find("8TeV" )!=std::string::npos){return multiplyGraph(  getVBFXSec8TeV(), getBRHtoZZ());}
             if(Name.find("7TeV" )!=std::string::npos){return multiplyGraph(  getVBFXSec7TeV(), getBRHtoZZ());}
             return NULL;
          }else{ //GGF
             if(Name.find("13TeV")!=std::string::npos){return multiplyGraph( getGGFXSec13TeV(), getBRHtoZZ());}
             if(Name.find("8TeV" )!=std::string::npos){return multiplyGraph(  getGGFXSec8TeV(), getBRHtoZZ());}
             if(Name.find("7TeV" )!=std::string::npos){return multiplyGraph(  getGGFXSec7TeV(), getBRHtoZZ());}
             return NULL;
          }
       }else if(Name.find("RsGrav")!=std::string::npos){
          if(Name.find("13TeV" )!=std::string::npos){return multiplyGraph(   getRsGravXSec13TeV(),  getBRRsGravtoZZ());}
          return NULL;
       }else if(Name.find("BulkGrav")!=std::string::npos){
          if(Name.find("13TeV" )!=std::string::npos){return multiplyGraph( getBulkGravXSec13TeV(), getBRBulkGravtoZZ());}
          return NULL;
       }else if(Name.find("Rad")!=std::string::npos){
          if(Name.find("13TeV" )!=std::string::npos){return multiplyGraph( getRadXSec13TeV(), getBRRadtoZZ());}
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

    TGraph *getXSecMELA( float cprime){

	double    mass[13] = {200,300,400,500,600,700,800,900,1000,1500,2000,2500,3000};
	double XSec100[13] = {0.0191670018,0.0774148727,0.0566354163,0.0324926779,0.0125720932,0.0055288361,0.0022651984,0.0011076429,0.0005393563,0.0000871796,0.0000271228,0.0000097798,0.0000076708};
	double XSec060[13] = {0.0524976539,0.2056358456,0.1650820193,0.0986099116,0.0387492768,0.0174677841,0.0073704531,0.0037865664,0.0018910718,0.0003383154,0.0000958260,0.0000303013,0.0000076708};
	double XSec030[13] = {0.2100001361,0.8076580339,0.6820127023,0.4096698103,0.1600257443,0.0734594886,0.0313789135,0.0162907538,0.0080814823,0.0014752691,0.0003679972,0.0000964785};
	double XSec010[13] = {1.9285049896,7.1367086972,6.1866400726,3.8800548591,1.4298414224,0.6505508647,0.2850172911,0.1480401099,0.0734280322,0.0134073773,0.0031436593,0.0007787744,0.0000076708};
	double XSec[13];

	if( Hxswg::utils::Equal(cprime, 1.0, 0.000001) ){
		std::cout << "Inside C'=1.0" << std::endl;
		for(unsigned int k=0;k<13;k++){ XSec[k]=XSec100[k]; }
	}else if( Hxswg::utils::Equal(cprime, 0.6, 0.000001) ){ 
		std::cout << "Inside C'=0.6" << std::endl;
		for(unsigned int k=0;k<13;k++){ XSec[k]=XSec060[k]; }
	}else if( Hxswg::utils::Equal(cprime, 0.3, 0.000001) ){
		std::cout << "Inside C'=0.3" << std::endl;
		for(unsigned int k=0;k<13;k++){ XSec[k]=XSec030[k]; }
	}else if( Hxswg::utils::Equal(cprime, 0.1, 0.000001) ){
		std::cout << "Inside C'=0.1" << std::endl;
		for(unsigned int k=0;k<13;k++){ XSec[k]=XSec010[k]; }
	}

	TGraph *XSecMELA= new TGraph( 13, mass, XSec);

	return XSecMELA;
    }


  }
}
