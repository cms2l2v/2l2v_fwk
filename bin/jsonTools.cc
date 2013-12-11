#include <string>
#include <vector>
#include <iostream>
#include <list>
#include <iterator>
#include <algorithm>
#include <map>
#include <unordered_map>

#include "UserCode/llvv_fwk/interface/JSONWrapper.h"

using namespace std;


JSONWrapper::Object merge_json(std::vector<std::string> jsonFiles){
   struct datasetinfo{std::map<string, string> prop;
      datasetinfo(){};
      datasetinfo(JSONWrapper::Object& dtagObj, bool isdata=false){
         for(unsigned int i=0;i<dtagObj.key.size();i++){
            prop[dtagObj.key[i]] = dtagObj.obj[i].key.size()==0 ? dtagObj.obj[i].toString() : dtagObj.obj[i].DumpToString();
         }
         if(prop["dset"]=="")prop["dset"]="Unknown";
         if(isdata) prop["isdata"] = "true";         
      };
   }; 
   std::map<std::string, datasetinfo> Datasets;

   if(jsonFiles.size()>1)printf("Checking differences between input files and uniformize the content assuming the first files are the references:\n");

   for(unsigned int J=0;J<jsonFiles.size();J++){
      JSONWrapper::Object Root(jsonFiles[J], true);
      std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
      for(size_t ip=0; ip<Process.size(); ip++){
          std::vector<JSONWrapper::Object> Samples = (Process[ip])["data"].daughters();
          for(size_t id=0; id<Samples.size(); id++){
             if(!Samples[id].isTag("dtag")) continue;
               string dtag = Samples[id].getString("dtag");

               datasetinfo info(Samples[id], Process[ip].getBool("isdata", false) );
               
               if(Datasets.find(dtag)!=Datasets.end()){
                  //check if they are the same
                  datasetinfo& infosaved = Datasets[dtag];
                  for(std::map<string, string>::iterator it = infosaved.prop.begin(); it != infosaved.prop.end(); it++){
                     string k = it->first;
                     if(infosaved.prop[k]  != info.prop[k] ){
                        if(info.prop[k]==""){
                           printf("FIX in file=%35s dtag=%35s var=%10s : %30s changed to %30s\n", jsonFiles[J].c_str(), dtag.c_str(), k.c_str(), info.prop[k].c_str(), infosaved.prop[k].c_str());
                           infosaved.prop[k]=info.prop[k];
                        }else{
                           printf("DIF in file=%35s dtag=%35s var=%10s : %30s differs to %30s  Please fix the files\n", jsonFiles[J].c_str(), dtag.c_str(), k.c_str(), infosaved.prop[k].c_str(), info.prop[k].c_str());
                        } 
                     }
                  }
                  for(std::map<string, string>::iterator it = info.prop.begin(); it != info.prop.end(); it++){
                     string k = it->first;
                     if(infosaved.prop.find(k)!=infosaved.prop.end())continue;
                     printf("ADD in file=%35s dtag=%35s var=%10s : %30s\n", jsonFiles[J].c_str(), dtag.c_str(), k.c_str(), info.prop[k].c_str());                     
                     infosaved.prop[k]=info.prop[k];
                  }
               }else{
//                  for(unsigned int k=0;k<Samples[id].key.size();k++){printf("%30s / %30s --> %30s - %30s\n", Process[ip].getString("tag").c_str(), Samples[id].getString("dtag").c_str(), Samples[id].key[k].c_str(), Samples[id].getString(Samples[id].key[k]).c_str());}
                  Datasets[dtag] = info;
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
          if(it->second.prop["isdata"]=="true")AllJson["proc"][I].add("isdata","true");
          AllJson["proc"][I].addArray("data");
          AllJson["proc"][I]["data"].addList();
          if(it->second.prop.find("dtag" )!=it->second.prop.end())AllJson["proc"][I]["data"][0].add("dtag" , JSONWrapper::removeWhiteSpace(it->second.prop["dtag" ]));
          if(it->second.prop.find("xsec" )!=it->second.prop.end())AllJson["proc"][I]["data"][0].add("xsec" , JSONWrapper::removeWhiteSpace(it->second.prop["xsec" ]));
          if(it->second.prop.find("br"   )!=it->second.prop.end())AllJson["proc"][I]["data"][0].add("br"   , JSONWrapper::removeWhiteSpace(it->second.prop["br"   ]));
          if(it->second.prop.find("split")!=it->second.prop.end())AllJson["proc"][I]["data"][0].add("split", JSONWrapper::removeWhiteSpace(it->second.prop["split"]));
          if(it->second.prop.find("dset" )!=it->second.prop.end())AllJson["proc"][I]["data"][0].add("dset" , JSONWrapper::removeWhiteSpace(it->second.prop["dset" ]));
          for(std::map<string, string>::iterator itp = it->second.prop.begin(); itp != it->second.prop.end(); itp++){
             if(itp->first=="dtag" || itp->first=="isdata" || itp->first=="xsec" || itp->first=="br" || itp->first=="split" || itp->first=="dset")continue;
             AllJson["proc"][I]["data"][0].add(itp->first, JSONWrapper::removeWhiteSpace(itp->second));
          }
   }
   return AllJson;
}


void make_multicrab_cfg(JSONWrapper::Object& Root, bool forData){
    FILE* pFile;
   pFile = fopen("multicrab.cfg", "w");
   fprintf(pFile, "[MULTICRAB]\n");
   fprintf(pFile, "\n");
   fprintf(pFile, "[COMMON]\n");
   if(forData){
      fprintf(pFile, "CMSSW.pset            = runObjectProducer_data_cfg.py\n");
      fprintf(pFile, "USER.user_remote_dir  = /store/user/quertenmont/13_08_15_2l2nu_EDMtuples/\n");
      fprintf(pFile, "CMSSW.lumi_mask       = /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt\n");
   }else{
     fprintf(pFile, "CMSSW.pset            = runObjectProducer_mc_cfg.py\n");
     fprintf(pFile, "USER.user_remote_dir  = /store/user/quertenmont/13_08_15_2l2nu_EDMtuples/\n");
   }
   fprintf(pFile, "\n");

   std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
   for(size_t ip=0; ip<Process.size(); ip++){
       if((forData && !Process[ip].getBool  ("isdata")) or (!forData && Process[ip].getBool  ("isdata")))continue;          
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
   bool multicrab   = false;
   bool multicrabMC = false;
   std::vector<string> inJsons;
   string outJson = "";    

   std::vector<string> args;
   for(int i=1;i<argc;i++){
     args.push_back(argv[i]);
     string arg(argv[i]);
     if(arg.find("--help")!=string::npos){
        printf("--help               --> print this helping text\n");
        printf("--in  [in1,...,inN)] --> input  Json files \n");
        printf("--out [out]          --> output Json file \n");
        printf("--multicrab          --> create a multicrab.cfg file for the data samples \n");
        printf("--multicrabMC        --> create a multicrab.cfg file for the MC samples \n");
	return 0;
     }

     else if(arg.find("--in")!=string::npos){ while(i+1<argc && string(argv[i+1]).find("--")!=0){inJsons.push_back(argv[i+1]);  i++;}  printf("input  Json files:\n"); for(unsigned int j=0;j<inJsons.size();j++){printf("\t - %s\n", inJsons[j].c_str());}  }
     else if(arg.find("--out")!=string::npos && i+1<argc){outJson = argv[i+1];  i++;  printf("output Json file : %s\n", outJson.c_str());}
     else if(arg.find("--multicrabMC")!=string::npos){ multicrabMC = true; }
     else if(arg.find("--multicrab"  )!=string::npos){ multicrab   = true; }
   }


   if(multicrab || multicrabMC){
      JSONWrapper::Object AllJson = merge_json(inJsons);
      make_multicrab_cfg(AllJson, !multicrabMC);
   }


   if(outJson!=""){
      JSONWrapper::Object AllJson = merge_json(inJsons); 
      FILE* pFile = fopen("tmp.json", "w"); 
      AllJson.Dump(pFile);
      fclose(pFile);
   }



   return 0;
}
