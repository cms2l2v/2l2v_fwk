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

	//Fixed Width
        double XSec100[13] = {0.000722821, 0.00974505, 0.0186783, 0.0170494, 0.0127224, 0.0123448, 0.00807319, 0.00694976, 0.00492894, 0.00189508, 0.000569369, 0.000196548, 5.8855e-05};
        double XSec010[13] = {0.00313343, 0.0670063, 0.218893, 0.204644, 0.145873, 0.137759, 0.086202, 0.072429, 0.0517826, 0.0192517, 0.00547137, 0.00184722, 0.000516527};
        double XSec005[13] = {0.00578392, 0.129529, 0.443913, 0.414981, 0.29471, 0.27903, 0.172004, 0.142393, 0.10418, 0.0389582, 0.0110436, 0.00363833, 0.00101596};

	//Narrow Width
	double XSec1[13] = {0.00952066, 0.0393434, 0.0357295, 0.0133214, 0.00500409, 0.00278995, 0.00112472, 0.000541774, 0.000283965, 4.16268e-05, 1.08999e-05, 4.70742e-06, 1.67446e-06};
	double XSec6[13] = {0.0260588, 0.105188, 0.103591, 0.040455, 0.0155903, 0.00888821, 0.00366068, 0.00184749, 0.000996092, 0.000161507, 3.83042e-05, 1.4616e-05, 4.60476e-06};
	double XSec3[13] = {0.103579, 0.414099, 0.426219, 0.168897, 0.0657464, 0.0379035, 0.0155764, 0.00798734, 0.00430719, 0.000703665, 0.000145865, 4.66086e-05, 1.20663e-05};
	double XSec01[13] = {0.928809, 3.63666, 3.77993, 1.58896, 0.604277, 0.352633, 0.140284, 0.0713195, 0.0401653, 0.00643292, 0.00124022, 0.00036292, 8.20836e-05};

	double XSec[13];

	//Wider Width
	if( Hxswg::utils::Equal(cprime, 100.0, 1E-6) ){
		for(unsigned int k=0;k<13;k++){ XSec[k]=XSec100[k]; }
	}else if( Hxswg::utils::Equal(cprime, 10.0, 1E-6) ){ 
		for(unsigned int k=0;k<13;k++){ XSec[k]=XSec010[k]; }
	}else if( Hxswg::utils::Equal(cprime, 5.0, 1E-6) ){	
		for(unsigned int k=0;k<13;k++){ XSec[k]=XSec005[k]; }
	}
	//Narrow Resonance
	else if( Hxswg::utils::Equal(cprime, 1.0, 1E-6) ){
		for(unsigned int k=0;k<13;k++){ XSec[k]=XSec1[k]; }
	}else if( Hxswg::utils::Equal(cprime, 0.6, 1E-6) ){ 
                for(unsigned int k=0;k<13;k++){ XSec[k]=XSec6[k]; }
        }else if( Hxswg::utils::Equal(cprime, 0.3, 1E-6) ){ 
                for(unsigned int k=0;k<13;k++){ XSec[k]=XSec3[k]; }
        }else if( Hxswg::utils::Equal(cprime, 0.1, 1E-6) ){ 
                for(unsigned int k=0;k<13;k++){ XSec[k]=XSec01[k]; }
        }

	TGraph *XSecMELA= new TGraph( 13, mass, XSec);

	return XSecMELA;
    }


  }
}
