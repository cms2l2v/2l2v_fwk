#include <iostream>
#include <boost/shared_ptr.hpp>

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

#include "DataFormats/FWLite/interface/LuminosityBlock.h"
#include "DataFormats/FWLite/interface/Run.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TEventList.h"
#include "TROOT.h"

using namespace std;

struct stRun {
   unsigned int runId;
   std::vector<unsigned int> lumiId;   
};

std::map<unsigned int, unsigned int> timeMap;

void GetLumiBlocks_Core(vector<string>& fileNames, std::vector<stRun*>& RunMap);
void DumpJson(const std::vector<stRun*>& RunMap, string FileName);
void DumpTime(const std::map<unsigned int, unsigned int>& TimeMap, string FileName);
void RemoveRunsAfter(unsigned int RunMax, const std::vector<stRun*>& RunMap, std::vector<stRun*>& NewRunMap);





int main(int argc, char* argv[])
{
  //##############################################
  //########    GLOBAL INITIALIZATION     ########
  //##############################################

  // check arguments
  if(argc<2){ std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl; exit(0); }
  
  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();
  
  // configure the process
  const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");

  bool isMC = runProcess.getParameter<bool>("isMC");
  if(isMC) return 0;
  std::vector<std::string> urls=runProcess.getParameter<std::vector<std::string> >("input");
  TString url = TString(urls[0]);
  TString outFileUrl(gSystem->BaseName(url));
  outFileUrl.ReplaceAll(".root","");
  TString outdir=runProcess.getParameter<std::string>("outdir");
  TString outUrl( outdir );

   std::vector<stRun*> RunMap;

   GetLumiBlocks_Core(urls, RunMap);

   DumpJson(RunMap, (outUrl + "/" + outFileUrl + "_lumi.json").Data());
   DumpTime(timeMap,(outUrl + "/" + outFileUrl + "_lumi.time").Data());
}  





void GetLumiBlocks_Core(vector<string>& fileNames, std::vector<stRun*>& RunMap)
{
   printf("Running\n");
   for(unsigned int f=0;f<fileNames.size();f++){
     cout << fileNames[f].c_str() << endl;
     //TFile file(fileNames[f].c_str() );
     TFile *file = TFile::Open(fileNames[f].c_str() );
      fwlite::LuminosityBlock ls( file);
      for(ls.toBegin(); !ls.atEnd(); ++ls){
          
        //printf("Run = %i --> Lumi =%lu\n",ls.luminosityBlockAuxiliary().run(), (unsigned long)ls.luminosityBlockAuxiliary().id().value());
        int RunIndex = -1;
        for(unsigned int r=0;r<RunMap.size();r++){
           if(RunMap[r]->runId==ls.luminosityBlockAuxiliary().run()){
              RunIndex = (int)r;
              break;
           }
        }

        if(RunIndex<0){
           stRun* tmp = new stRun();
           tmp->runId=ls.luminosityBlockAuxiliary().run();
           tmp->lumiId.push_back(ls.luminosityBlockAuxiliary().id().value());
           RunMap.push_back(tmp);           
           //std::sort(RunMap.begin(), RunMap.end(), stRunLess);

           timeMap[tmp->runId] = ls.luminosityBlockAuxiliary().beginTime().unixTime()/3600;
        }else{
            stRun* tmp = RunMap[RunIndex];
           int LumiIndex = -1;
           for(unsigned int l=0;l<tmp->lumiId.size();l++){
              //printf("%lu vs %lu\n",tmp->lumiId[l], (unsigned long) ls.luminosityBlockAuxiliary().id().value() );
              if(tmp->lumiId[l]== (unsigned int) ls.luminosityBlockAuxiliary().id().value()){
                 LumiIndex = (int)l;
                 break;
              }
           }
           if(LumiIndex<0){
               tmp->lumiId.push_back((unsigned int) ls.luminosityBlockAuxiliary().id().value());
               std::sort(tmp->lumiId.begin(), tmp->lumiId.end());
            }
        }      
      }printf("\n");
   }
}

void RemoveRunsAfter(unsigned int RunMax, const std::vector<stRun*>& RunMap, std::vector<stRun*>& NewRunMap){
   for(unsigned int r=0;r<RunMap.size();r++){
      if(RunMap[r]->runId<RunMax)NewRunMap.push_back(RunMap[r]);
   }
}


void DumpJson(const std::vector<stRun*>& RunMap, string FileName){
   FILE* json = fopen(FileName.c_str(),"w");
   if(!json)printf("Could Not open file: %s\n", FileName.c_str());
   fprintf(json,"{");
   for(unsigned int r=0;r<RunMap.size();r++){
      stRun* tmp =  RunMap[r];
      fprintf(json,"\"%i\": [",tmp->runId);
      unsigned int l=0;
      while(l<tmp->lumiId.size()){
         unsigned int FirstLumi = tmp->lumiId[l];
         unsigned Size=0; 
         for(unsigned int l2=l;l2<tmp->lumiId.size() && FirstLumi+l2-l==tmp->lumiId[l2]; l2++){Size++;}
         fprintf(json,"[%i, %i]",FirstLumi,FirstLumi+Size-1);
         l+=Size;
         if(l<tmp->lumiId.size()) fprintf(json,",");
      }
      fprintf(json,"] ");
      if(r<RunMap.size()-1)fprintf(json,",");
   }  
   fprintf(json,"}");   
   fclose(json);
}



void DumpTime(const std::map<unsigned int, unsigned int>& TimeMap, string FileName){
   FILE* pFile = fopen(FileName.c_str(),"w");
   for(std::map<unsigned int, unsigned int>::iterator it = timeMap.begin(); it!=timeMap.end(); it++){
      fprintf(pFile, "%i,%i\n", it->first, it->second);
   }
   fclose(pFile);
}



