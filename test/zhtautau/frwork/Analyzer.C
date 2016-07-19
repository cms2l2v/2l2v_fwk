#include "FRAnalyzer.h"
#include "TH1F.h"
#include "TH2F.h"
#include <sstream>
#include <iostream>
#include <vector>
#include "TApplication.h"
#include "TRint.h"

using namespace std;

int main(int argc, char* argv[]){
  
  if(argc<7){
    cout<<"FRAnalyzer: Systematic plots producer!"<<endl;
    cout<<"  Usage   : ./Analyser "<<endl;
    cout<<"  Needed Options : "<<endl;
    cout<<"                   --File1       [File1 Name]"<<endl;
    cout<<"                   --File2       [File2 Name]"<<endl;
    cout<<"                   --Index       [Index]"<<endl;
    cout<<"                   --Ntoy        [Ntoy]"<<endl;
    cout<<"                   --Data        [Data]"<<endl;
    cout<<"                   --Dir         [Dir]"<<endl;
    cout<<"                   --NoXserver   [Disable Xserver]"<<endl;
    cout<<"  Other Options  :"<<endl;
    //generalmente qui si mettono quelle opzioni che se non le passi vengono settate di defult dal codice
    cout<<"                   --path ]"<<endl;
    exit(0);
  };
  
  string File1, File2, Dir;
  int Index;
  int Ntoy;
  bool Data=false;
  bool NoXserver=false;

  for(int a=1; a<argc; a++){
    if(!strcmp(argv[a],"--File1")){
      File1 = argv[a+1];
    }
    else if(!strcmp(argv[a],"--File2")){
      File2 = argv[a+1];
    }
    else if(!strcmp(argv[a],"--Index")){
      string value = argv[a+1];
      istringstream(value) >> Index;
    }
    else if(!strcmp(argv[a],"--Ntoy")){
      string value = argv[a+1];
      istringstream(value) >> Ntoy;
    }
    else if(!strcmp(argv[a],"--Data")){
      Data = true;
    }
    else if(!strcmp(argv[a],"--Dir")){
      Dir = argv[a+1];
    }
    else if(!strcmp(argv[a],"--NoXserver")){
      NoXserver = true;
    }
  }
  
  /*debug*/
  cout<< " Options Passed:" <<endl;
  cout<< "           File1       = " << File1     << endl;
  cout<< "           File2       = " << File2     << endl;
  cout<< "           Index       = " << Index     << endl;
  cout<< "           Ntoy        = " << Ntoy     << endl;
  cout<< "           Data        = " << Data      << endl;
  cout<< "           Dir         = " << Dir       << endl;
  cout<< "           NoXserver   = " << NoXserver << endl;
  
#ifdef WITHRINT
  TRint *myApp;
  if(!NoXserver){
    myApp = new TRint("RootSession",&argc,argv,NULL,0);
    cout<<" Using Rint: Opening ROOT ..." << endl;
  }
#else
  TApplication *myApp;
  if(!NoXserver){
    myApp = new TApplication("myApp",0,0);
    cout<<" Using TApplication: Gui will be used ... " << endl;
  }
#endif
    
  FRAnalyzer *PM = new FRAnalyzer(File1,File2,Index,Ntoy,Data,Dir);
  cout << "Creating Dir..."<<endl;
  PM->CreateDir();
  cout << "Reading Files..."<<endl;
  PM->ReadFiles();
  cout << "Initializing..."<<endl;
  PM->Initialize();
  cout << "After initialization.."<<endl;  
  if(Data) PM->RetrieveResults();
  if(!Data) PM->RetrieveCTResults();
  if(!NoXserver)
    myApp->Run();

  return 0;
}
