#include <string>
#include <vector>
#include <iostream>
#include <list>
#include <iterator>
#include <algorithm>
#include <map>
#include <unordered_map>

#include "UserCode/llvv_fwk/src/JSONWrapper.cc"



JSONWrapper::Object merge_json(std::vector<std::string> jsonFiles){
   struct datasetinfo{std::string dset; std::string split; bool isdata; string br; string xsec;};
   std::map<std::string, datasetinfo> Datasets;

   for(unsigned int J=0;J<jsonFiles.size();J++){
      JSONWrapper::Object Root(jsonFiles[J], true);
      std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
      for(size_t ip=0; ip<Process.size(); ip++){
          std::vector<JSONWrapper::Object> Samples = (Process[ip])["data"].daughters();
          for(size_t id=0; id<Samples.size(); id++){
             if(!Samples[id].isTag("dtag")) continue;

               if(Datasets.find(Samples[id].getString("dtag"))!=Datasets.end()){
                  continue;//dataset already exist
                  //need to add a function to check if they are the same or not and print a warning if they are not
               }else{
                  datasetinfo info;
                  info.dset   = Samples[id].getString("dset", "");
                  info.split  = Samples[id].getString("split", "");
                  info.br     = Samples[id]["br"].DumpToString(1);
                  info.xsec   = Samples[id].getString("xsec" , "");
                  info.isdata = Process[ip].getBool  ("isdata"); 
                  Datasets[Samples[id].getString("dtag")] = info;
               }
          }
      }
   }

   JSONWrapper::Object AllJson;
   AllJson.addArray("proc");
   for(std::map<std::string, datasetinfo>::iterator it = Datasets.begin(); it!=Datasets.end();it++){
          int I = AllJson["proc"].daughters().size();
          AllJson["proc"].addList(); 
          AllJson["proc"][I].add("tag", it->first);      
          if(it->second.isdata)   AllJson["proc"][I].add("isdata","true");
          AllJson["proc"][I].addArray("data");
          AllJson["proc"][I]["data"].addList();
          if(it->first       !="")AllJson["proc"][I]["data"][0].add("dtag", it->first);      
          if(it->second.dset !="")AllJson["proc"][I]["data"][0].add("dset", it->second.dset);
          if(it->second.split!="")AllJson["proc"][I]["data"][0].add("split",it->second.split);
          if(it->second.xsec !="")AllJson["proc"][I]["data"][0].add("xsec" ,it->second.xsec);
          if(it->second.br   !="")AllJson["proc"][I]["data"][0].add("br"   ,it->second.br);
   }

   return AllJson;
}




JSONWrapper::Object make_multicrab_json(std::string jsonFile){
   struct datasetinfo{std::string dset; std::string split; bool isdata;};
   std::map<std::string, datasetinfo> Datasets;

   JSONWrapper::Object Root(jsonFile, true);
   std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
   for(size_t ip=0; ip<Process.size(); ip++){
       std::vector<JSONWrapper::Object> Samples = (Process[ip])["data"].daughters();
       for(size_t id=0; id<Samples.size(); id++){
          if(!Samples[id].isTag("dtag")) continue;

            datasetinfo info;
            info.dset   = Samples[id].getString("dset", "Unknown");
            info.split  = Samples[id].getString("split", "1");
            info.isdata = Process[ip].getBool  ("isdata"); 
            Datasets[Samples[id].getString("dtag")] = info;

//          printf("%s\n", Samples[id].getString("dtag","").c_str());
       }
   }

   JSONWrapper::Object AllJson;
   AllJson.addArray("proc");
   for(std::map<std::string, datasetinfo>::iterator it = Datasets.begin(); it!=Datasets.end();it++){
          int I = AllJson["proc"].daughters().size();
          AllJson["proc"].addList(); 
          AllJson["proc"][I].add("tag", it->first);      
          AllJson["proc"][I].addArray("data");
          AllJson["proc"][I]["data"].addList();
          AllJson["proc"][I]["data"][0].add("dtag", it->first);
          AllJson["proc"][I]["data"][0].add("dset", it->second.dset);
          AllJson["proc"][I]["data"][0].add("split",it->second.split);
          if(it->second.isdata)
          AllJson["proc"][I].add("isdata","true");
   }

   return AllJson;
}

void make_multicrab_cfg(JSONWrapper::Object& Root){
    FILE* pFile;
   /////////////  MC ///////////////////
   pFile = fopen("multicrab_mc.cfg", "w");
   fprintf(pFile, "[MULTICRAB]\n");
   fprintf(pFile, "\n");
   fprintf(pFile, "[COMMON]\n");
   fprintf(pFile, "CMSSW.pset            = runObjectProducer_mc_cfg.py\n");
   fprintf(pFile, "USER.user_remote_dir  = /store/user/quertenmont/13_08_15_2l2nu_EDMtuples/\n");
   fprintf(pFile, "\n");

   std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
   for(size_t ip=0; ip<Process.size(); ip++){
       if(Process[ip].getBool  ("isdata"))continue;          
       std::vector<JSONWrapper::Object> Samples = (Process[ip])["data"].daughters();
       for(size_t id=0; id<Samples.size(); id++){
          if(!Samples[id].isTag("dtag")) continue;

          fprintf(pFile, "[%s]\n",Samples[id].getString("dtag").c_str());
          fprintf(pFile, "CMSSW.datasetpath     = %s\n", Samples[id].getString("dset", "").c_str());
          fprintf(pFile, "\n");
       }
    }
   fclose(pFile);

   /////////////  data ///////////////////
   pFile = fopen("multicrab_data.cfg", "w");
   fprintf(pFile, "[MULTICRAB]\n");
   fprintf(pFile, "\n");
   fprintf(pFile, "[COMMON]\n");
   fprintf(pFile, "CMSSW.pset            = runObjectProducer_data_cfg.py\n");
   fprintf(pFile, "USER.user_remote_dir  = /store/user/quertenmont/13_08_15_2l2nu_EDMtuples/\n");
   fprintf(pFile, "CMSSW.lumi_mask       = /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt\n");
   fprintf(pFile, "\n");

   for(size_t ip=0; ip<Process.size(); ip++){
       if(!Process[ip].getBool  ("isdata"))continue;   
       std::vector<JSONWrapper::Object> Samples = (Process[ip])["data"].daughters();
       for(size_t id=0; id<Samples.size(); id++){
          if(!Samples[id].isTag("dtag")) continue;

          fprintf(pFile, "[%s]\n",Samples[id].getString("dtag").c_str());
          fprintf(pFile, "CMSSW.datasetpath     = %s\n", Samples[id].getString("dset", "").c_str());
          fprintf(pFile, "\n");
       }
    }
   fclose(pFile);
}


int main(int argc, char* argv[]){
   std::vector<string> args;
   for(int i=1;i<argc;i++){
     args.push_back(argv[i]);
     string arg(argv[i]);
     if(arg.find("--help")!=string::npos){
        printf("--help   --> print this helping text\n");
//        printf("--splitCanvas --> (only for 2D plots) save all the samples in separated pltos\n");
	return 0;
     }
//     if(arg.find("--iLumi"  )!=string::npos && i+1<argc){ sscanf(argv[i+1],"%lf",&iLumi); i++; printf("Lumi = %f\n", iLumi); }
   }


 
   JSONWrapper::Object AllJson = merge_json(args);


//   JSONWrapper::Object AllJson = make_multicrab_json(argv[1]);
//   make_multicrab_cfg(AllJson);

 
   FILE* pFile = fopen("tmp.json", "w"); 
   AllJson.Dump(pFile);
   fclose(pFile);



   return 0;
}
