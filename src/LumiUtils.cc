#include "UserCode/llvv_fwk/interface/LumiUtils.h"





namespace lumiUtils
{

    bool GoodLumiFilter::isGoodLumi(edm::RunNumber_t run, edm::LuminosityBlockNumber_t lumi){
       if(lumisToProcess.size()!=0){
          edm::LuminosityBlockID lumiID = edm::LuminosityBlockID(run, lumi);
          edm::LuminosityBlockRange lumiRange = edm::LuminosityBlockRange(lumiID, lumiID);
          bool(*lt)(edm::LuminosityBlockRange const&, edm::LuminosityBlockRange const&) = &edm::lessThan;
          if(!binary_search_all(lumisToProcess, lumiRange, lt))return false;//this is not a good lumiBlock
       }
       return true;
    }




   void GoodLumiFilter::FindLumiInFiles(vector<string>& fileNames)
   {
      for(unsigned int f=0;f<fileNames.size();f++){
        TFile *file = TFile::Open(fileNames[f].c_str() );
         fwlite::LuminosityBlock ls( file);
         for(ls.toBegin(); !ls.atEnd(); ++ls){
           if(!isGoodLumi(ls.luminosityBlockAuxiliary().run(),ls.luminosityBlockAuxiliary().id().value()))continue;
            
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

   void GoodLumiFilter::RemoveRunsAfter(unsigned int RunMax){
      std::vector<stRun*> NewRunMap;
      for(unsigned int r=0;r<RunMap.size();r++){
         if(RunMap[r]->runId<RunMax)NewRunMap.push_back(RunMap[r]);
      }
      RunMap = NewRunMap;
   }


   void GoodLumiFilter::DumpToJson(string FileName){
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

   void GoodLumiFilter::DumpToTime(string FileName){
      FILE* pFile = fopen(FileName.c_str(),"w");
      for(std::map<unsigned int, unsigned int>::iterator it = timeMap.begin(); it!=timeMap.end(); it++){
         fprintf(pFile, "%i,%i\n", it->first, it->second);
      }
      fclose(pFile);
   }






}
