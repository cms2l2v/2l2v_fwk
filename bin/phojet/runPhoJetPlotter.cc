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
#include <unordered_map>

#include "TFile.h"
#include "TObject.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"


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

TObject* GetObjectFromPath(TDirectory* File, std::string Path, bool GetACopy=false)
{
  // std::cout << Path << std::endl;
  size_t pos = Path.find("/");
  // std::cout << ">>> " << pos  << std::endl;
  if(pos < 256){
    std::string firstPart = Path.substr(0,pos);
    std::string endPart   = Path.substr(pos+1,Path.length());
    TDirectory* TMP = (TDirectory*)File->Get(firstPart.c_str());
    if(TMP!=NULL){
      TObject* TMP2 =  GetObjectFromPath(TMP,endPart,GetACopy);
      return TMP2;
    }
    return NULL;
  }else{
    TObject* TMP = File->Get(Path.c_str());
    if(GetACopy){	return TMP->Clone();
    }else{            return TMP;
    }
  }
}


void GetListOfObject(JSONWrapper::Object& Root,
		     std::string RootDir,
		     std::list<NameAndType>& histlist,
		     TDirectory* dir=NULL,
		     std::string parentPath=""){
  if(parentPath=="" && !dir){
    int dataProcessed = 0;
    int signProcessed = 0;
    int bckgProcessed = 0; 
  
    std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
    std::cout << "Num of Process: " << Process.size() << std::endl; 

    std::unordered_map<std::string, bool> FileExist;

    // loop over all procs
    for(size_t ip=0; ip<Process.size(); ip++){
      bool isData (  Process[ip]["isdata"].toBool()  );
      bool isSign ( !isData &&  Process[ip].isTag("spimpose")
		    && Process[ip]["spimpose"].toBool());
      bool isMC = !isData && !isSign; 

      std::string filtExt("");
      if(Process[ip].isTag("mctruthmode") )
	{ char buf[255];
	  sprintf(buf, "_filt%d", (int)Process[ip]["mctruthmode"].toInt());
	  filtExt += buf;
	}
    
      std::vector<JSONWrapper::Object> Samples = (Process[ip])["data"].daughters();
      // loop over all samples 
      for(size_t id=0; id<Samples.size(); id++){
	int split = Samples[id].getInt("split", 1); // default 1 
	std::cout << "Split = " << split << std::endl;

	// loop over all files with CRAB convention
	for(int s=1; s<=split; s++){
	  std::string segmentExt;
	  if(split>1) {
	    char buf[255];
	    sprintf(buf,"_%i",s);
	    segmentExt += buf;
	  }
	  std::string FileName = RootDir 
	    + Samples[id].getString("dtag", "") + "/"
	    + "output" + Samples[id].getString("suffix","") + segmentExt + ".root";
	  // std::cout << FileName << std::endl;
	
	  TFile* File = new TFile(FileName.c_str());
	  bool& fileExist = FileExist[FileName];
	  if(!File || File->IsZombie() || !File->IsOpen() ||
	     File->TestBit(TFile::kRecovered) ){
	    fileExist=false;
	    continue; 
	  }else{
	    fileExist=true;
	  }
	
	  //do the following only for the first file
	  if(s>1) continue;

	  printf("Adding all objects from %25s to the list of considered objects\n",
		 FileName.c_str());

	  //just to make it faster, only consider the first 3 sample of a same kind
	  if(isData){if(dataProcessed>=1){ File->Close(); continue;}else{dataProcessed++;}}
	  if(isSign){if(signProcessed>=2){ File->Close(); continue;}else{signProcessed++;}}
	  if(isMC  ){if(bckgProcessed>0) { File->Close(); continue;}else{bckgProcessed++;}}

	  GetListOfObject(Root, RootDir, histlist, (TDirectory*)File, "" );
	  File->Close();
	} // end on all files 
      } // end on all samples 
    } // end on all procs

    // print out missing or corrupted files if any 
    for(std::unordered_map<std::string, bool>::iterator it = FileExist.begin();
	it!=FileExist.end(); it++){
      if(!it->second)
	printf("[INFO] missing or corrupted file:  %s\n", it->first.c_str());
    }

    return ; 
    
  } // end of no parentPath or no dir case
  
  if (dir==NULL) return;
  // std::cout << "dir = " << dir << std::endl;
  // std::cout << "parentPath = " << parentPath << std::endl;

  TList* list = dir->GetListOfKeys();

  // loop over list 
  for(int i=0;i<list->GetSize();i++){
    // std::cout << list->At(i)->GetName() << std::endl;
    TObject* tmp = GetObjectFromPath(dir,list->At(i)->GetName(),false);

    if(tmp->InheritsFrom("TDirectory")){
      printf("found object from TDirectory\n");
      GetListOfObject(Root, RootDir, histlist, (TDirectory*)tmp, parentPath
		      + list->At(i)->GetName() + "/" );
    }else if (tmp->InheritsFrom("TTree")){ 
      printf("found one object inheriting from a ttree\n");
      // isTree = type 4, isIndexPlot = false 
      histlist.push_back(NameAndType(parentPath+list->At(i)->GetName(), 4, false ) );
    }else if (tmp->InheritsFrom("TH1")){
      // printf("found one object inheriting from TH1\n");
      int  type = 0;
      if(tmp->InheritsFrom("TH1")) type++;
      if(tmp->InheritsFrom("TH2")) type++;
      if(tmp->InheritsFrom("TH3")) type++;
      bool hasIndex = std::string(((TH1*)tmp)->GetXaxis()->GetTitle()).find("cut index")
	<std::string::npos;
      if(hasIndex){type=1;}
      histlist.push_back(NameAndType(parentPath+list->At(i)->GetName(), type, hasIndex ) );
    }else{
      printf("The file contain an unknown object named %s\n", list->At(i)->GetName() );
    }
    delete tmp;
  } // end list loop 

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
      inDir = argv[i+1];  i++;  printf("inDir = %s\n", inDir.c_str()); }

    if(arg.find("--outDir" )!=std::string::npos && i+1<argc){
      outDir = argv[i+1];  i++;  printf("outDir = %s\n", outDir.c_str());}

    if(arg.find("--outFile")!=std::string::npos && i+1<argc){
      outFile = argv[i+1];  i++; printf("output file = %s\n", outFile.c_str());}

    if(arg.find("--json"   )!=std::string::npos && i+1<argc){
      jsonFile = argv[i+1];  i++;}
  }

  system( (std::string("mkdir -p ") + outDir).c_str());

  JSONWrapper::Object Root(jsonFile, true);
  std::list<NameAndType> histlist;
  GetListOfObject(Root, inDir, histlist);
  std::cout << "Total of hists found: " << histlist.size() << std::endl;
  histlist.sort();
  histlist.unique();   
  
}  







