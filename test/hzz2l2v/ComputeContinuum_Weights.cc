#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdio.h>

#include <TH1D.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TFile.h>
#include <TPaveStats.h>
#include <TPad.h>

using namespace std;

void ComputeContinuum_Weights(bool isContinuumWg, bool isSignalXSec){

	bool isContinuum(isContinuumWg);
	bool isSigCrossSection(isSignalXSec);

	string FileCont = "plotter_2017_24_04_Continuum_MCOnly.root";
	string  FileSig = "plotter_2017_01_04_full2016_MCOnly.root"; //"plotter_2017_01_24_SignalOnly.root";
	string      Dir = "/user/amagitte/HZZ_2L2Nu_Analysis_2016/LHCP2017/CMSSW_8_0_26_patch1/src/UserCode/llvv_fwk/test/hzz2l2v/";

	double    YR4[13] = {18.12,9.823,9.516,4.538,2.006,0.9235,0.4491,0.2301,0.1233,0.008913,0.001084,0.0001778,0.00003502}; 
	double   BRH2[13] = {0.255,0.307,0.269,0.261,0.272,0.285,0.296,0.304,0.311,1,1,1,1};
	string   Mass[13] = {"200","300","400","500","600","700","800","900","1000","1500","2000","2500","3000"};
	string      Cp[3] = {"100.00","10.00","5.00"}; 
	string Channel[2] = {"ggH","qqH"};
	double Luminosity = 35866.932;
	double qqZZXSec = 0.103744*6.433832215816802E-4;
	double ggZZXSec = 0.01898*2;
	double  BR2l2nu = 0.026926;
	double     BR4l = 0.0022657219;

	FILE *outfile_ContsF = fopen( "Continuum_sF.txt", "w");
	FILE *outfile_Integr = fopen( "XSection_Sig.txt", "w");

	//Loop on the channel
	for(unsigned int ch=0; ch<2; ch++){

		//Loop on the Cprime
		for(unsigned int cp=0; cp<3; cp++){

			fprintf( outfile_ContsF, "C'= %s \n", Cp[cp].c_str());
			fprintf( outfile_Integr, "C'= %s \n", Cp[cp].c_str());

			//Loop on the Mass
			for(unsigned int m=0; m<13; m++){

				double sF = 0.0;				

				string  NameHistoSig = "all_higgsMass_shape_cp"+Cp[cp]+"_brn0.00";
				string NameHistoXSec = "all_higgsMass_shape_cp"+Cp[cp]+"_brn0.00"; // "all_MELA_weights_cp"+Cp[cp]+"_brn0.00";
				string NameHistoBckg = "all_higgsMass_nnlo_raw";	
			
				//Compute the Continuum Weights
				string  DirSig = Channel[ch]+"("+Mass[m]+")_ContinuumOnly";
				string DirBckg = "ZZ_filt1113";
				
				int intMin=0; int intMax=0;

				if(isContinuum){
					TFile *RFile = new TFile((Dir+FileCont).c_str());
					if( RFile->GetDirectory(DirBckg.c_str())!=0 && RFile->GetDirectory(DirSig.c_str())!=0 ){ 
					
						RFile->cd(DirBckg.c_str());
						TH1F *HBckg = (TH1F*)gDirectory->Get(NameHistoBckg.c_str());
						RFile->cd(DirSig.c_str());
						TH1F  *HSig = (TH1F*)gDirectory->Get(NameHistoSig.c_str());
						if(ch==0){ sF = HBckg->Integral(99,310)/HSig->Integral(99,310);}
						else if(ch==1){ sF = (HBckg->Integral()/HSig->Integral())*(qqZZXSec/ggZZXSec)*(BR2l2nu/BR4l); }	
					} else { sF=1.0; }

					fprintf( outfile_ContsF, "Channel: %s Cprime: %s Mass: %s sF: %20.19f \n", Channel[ch].c_str(), Cp[cp].c_str(), Mass[m].c_str(), sF);
				}

				if(isSigCrossSection){
					//Compute the Signal Cross-Section
					double BR = BRH2[m]*0.029264;	
					string  DirSOnly = Channel[ch]+"("+Mass[m]+")_SOnly";
					TFile *SFile = new TFile((Dir+FileSig).c_str());
					SFile->cd(DirSOnly.c_str());
					TH1F  *HSigOnly = (TH1F*)gDirectory->Get(NameHistoXSec.c_str());
                                	double XSec = (HSigOnly->Integral()/Luminosity);
                                	double XSecNoBR = XSec / BR;	
					fprintf( outfile_Integr, "Channel: %s Cprime: %s Mass: %s NEvents: %10.5f XSec: %20.10f XSec(BRs Included): %10.5f XSec(YR4): %10.5f  \n", Channel[ch].c_str(), Cp[cp].c_str(), Mass[m].c_str(), HSigOnly->Integral(), XSec, XSecNoBR, YR4[m]);	
				}
			}
		}

		fprintf( outfile_ContsF, "\n");
		fprintf( outfile_Integr, "\n");
	}

	//outfile_ContsF->Close();
	//outfile_Integr->Close();
}
