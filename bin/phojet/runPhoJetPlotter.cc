// -*- C++ -*-
//
// Xin Shi <Xin.Shi@cern.ch>
// Mon Jan 5 13:10:33 EST 2015
// 
// Plotter for the photon + jet analysis 
// 


#include <iostream>
#include <string>


void print_usage() {
  printf("--help    --> print this helping text\n");
  printf("--inDir   --> path to the directory containing the .root files to process\n");
  printf("--outDir  --> path of the directory that will contains the output plots and tables\n");
  printf("--outFile --> path of the output summary .root file\n");
  printf("--json    --> containing list of process (and associated style) to process to process\n");
}


int main(int argc, char* argv[])
{

  std::string inDir   = "results/";
  std::string jsonFile = "../../data/phojet/phys14_samples.json";
  std::string outDir  = "plots/";
  std::string outFile = "plotter.root";

  // check arguments
  // if(argc<2){
  //   print_usage();
  //   return 0;
  // }
  
  for(int i=1;i<argc;i++){
    std::string arg(argv[i]);
    if (arg.find("--help") != std::string::npos){
      print_usage();
      return 0; 
    }

    if( arg.find("--inDir"  )!=std::string::npos && i+1<argc) {
      inDir = argv[i+1];  i++;  printf("inDir = %s\n", inDir.c_str());
    }

    if(arg.find("--outDir" )!=std::string::npos && i+1<argc){
      outDir = argv[i+1];  i++;  printf("outDir = %s\n", outDir.c_str());
    }

    if(arg.find("--outFile")!=std::string::npos && i+1<argc){
      outFile = argv[i+1];  i++; printf("output file = %s\n", outFile.c_str());
    }

    if(arg.find("--json"   )!=std::string::npos && i+1<argc){
      jsonFile = argv[i+1];  i++;
    }
  }

  system( (std::string("mkdir -p ") + outDir).c_str());
  
}  







