// -*- C++ -*-
//
// Xin Shi <Xin.Shi@cern.ch>
// Mon Jan 5 13:10:33 EST 2015
// 
// Plotter for the photon + jet analysis 
// 


#include <iostream>
#include <string>
#include <list>

#include "TDirectory.h"

#include "UserCode/llvv_fwk/interface/JSONWrapper.h"

struct NameAndType{
   std::string name;
   int type; 
   bool isIndexPlot;
   NameAndType(std::string name_,  int type_, bool isIndexPlot_){name = name_; type = type_; isIndexPlot = isIndexPlot_;}
   bool is1D()  {return type==1;}
   bool is2D()  {return type==2;}
   bool is3D()  {return type==3;}
   bool isTree(){return type==4;}
   bool operator==(const NameAndType& a){ return a.name == name;}
   bool operator< (const NameAndType& a){ return a.name < name;}
 };


void print_usage() {
  printf("--help    --> print this helping text\n");
  printf("--inDir   --> path to the directory containing the .root files to process\n");
  printf("--outDir  --> path of the directory that will contains the output plots and tables\n");
  printf("--outFile --> path of the output summary .root file\n");
  printf("--json    --> containing list of process (and associated style) to process to process\n");
}

void GetListOfObject(JSONWrapper::Object& Root,
		     std::string RootDir,
		     std::list<NameAndType>& histlist){

  std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
  std::cout << "Num of Process: " << Process.size() << std::endl; 
  
  // loop over all procs
  for(size_t ip=0; ip<Process.size(); ip++){
    bool isData (  Process[ip]["isdata"].toBool()  );
    std::string filtExt("");
    if(Process[ip].isTag("mctruthmode") )
      { char buf[255];
	sprintf(buf, "_filt%d", (int)Process[ip]["mctruthmode"].toInt());
	filtExt += buf;
      }
    
    std::vector<JSONWrapper::Object> Samples = (Process[ip])["data"].daughters();
    // loop over all samples 
    for(size_t id=0; id<Samples.size(); id++){
      int split = Samples[id].getInt("split", 1);
      std::cout << "Split = " << split << std::endl;
    } // end on all samples 
    
  } // end on all procs 
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

  JSONWrapper::Object Root(jsonFile, true);
  std::list<NameAndType> histlist;
  GetListOfObject(Root,inDir,histlist);
    
}  







