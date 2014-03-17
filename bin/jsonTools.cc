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


struct datasetinfo{std::map<string, string> prop;
   datasetinfo(){};
   datasetinfo(JSONWrapper::Object& dtagObj, bool isdata=false){
      for(unsigned int i=0;i<dtagObj.key.size();i++){
         prop[dtagObj.key[i]] = dtagObj.obj[i].key.size()==0 ? dtagObj.obj[i].toString() : dtagObj.obj[i].DumpToString();
      }
      if(prop["dset" ]=="")prop["dset" ]="Unknown";
      if(prop["xsec" ]=="")prop["xsec" ]="1.0";
      if(prop["br"   ]=="")prop["br"   ]="[1.0]";
      if(prop["split"]=="")prop["split"]="1";
      if(isdata) prop["isdata"] = "true";         
   };
}; 


JSONWrapper::Object merge_json(std::vector<std::string> jsonFiles){
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
          if(it->second.prop.find("dtag" )!=it->second.prop.end())AllJson["proc"][I]["data"][0].add("dtag" , JSONWrapper::removeWhiteSpace(it->second.prop["dtag" ]), 50);
          if(it->second.prop.find("split")!=it->second.prop.end())AllJson["proc"][I]["data"][0].add("split", JSONWrapper::removeWhiteSpace(it->second.prop["split"]), 5);
          if(it->second.prop.find("xsec" )!=it->second.prop.end())AllJson["proc"][I]["data"][0].add("xsec" , JSONWrapper::removeWhiteSpace(it->second.prop["xsec" ]), 20);
          if(it->second.prop.find("br"   )!=it->second.prop.end())AllJson["proc"][I]["data"][0].add("br"   , JSONWrapper::removeWhiteSpace(it->second.prop["br"   ]), 30);
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
          if(!Samples[id].isTag("dtag") || Samples[id].getString("dset", "")=="EOS" || Samples[id].getString("dbsurl", "")=="local" ) continue;

          fprintf(pFile, "[%s]\n",Samples[id].getString("dtag").c_str());
          fprintf(pFile, "CMSSW.datasetpath     = %s\n", Samples[id].getString("dset", "").c_str());
          if(Samples[id].isTag("dbsurl"))
          fprintf(pFile, "CMSSW.dbs_url         = %s\n", Samples[id].getString("dbsurl", "").c_str());
          fprintf(pFile, "\n");
       }
   }
   fclose(pFile);
}



void fix_json(JSONWrapper::Object& Ref, std::string jsonFileToBeFixed){
   printf("\n");
   printf("------------------------------------------------------------------------------\n");
   printf("Fixing %50s, the file content will be overwritten\n", jsonFileToBeFixed.c_str());
   printf("------------------------------------------------------------------------------\n");

   std::map<std::string, datasetinfo> Datasets;
   //Load the datasets from the reference
   if(true){//just here to avoid variable redeclaration error
      std::vector<JSONWrapper::Object> Process = Ref["proc"].daughters();
      for(size_t ip=0; ip<Process.size(); ip++){
          std::vector<JSONWrapper::Object> Samples = (Process[ip])["data"].daughters();
          for(size_t id=0; id<Samples.size(); id++){
             if(!Samples[id].isTag("dtag")) continue;
             string dtag = Samples[id].getString("dtag");
             datasetinfo info(Samples[id], Process[ip].getBool("isdata", false) );
             Datasets[dtag] = info;                         
         }
      }
   }

   //Load the info from the file to fix
   if(true){//just here to avoid variable redeclaration error
      //prepare the output json
      JSONWrapper::Object OutJson;
      OutJson.addArray("proc");
      int I=0;

      JSONWrapper::Object Root(jsonFileToBeFixed, true);
      std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
      for(size_t ip=0; ip<Process.size(); ip++){   
          OutJson["proc"].addList();

          //copy all json properties except "data" that is taken from the reference
          for(unsigned int i=0;i<Process[ip].key.size();i++){
             if(Process[ip].key[i]=="data")continue;//data field is specialy treated afterward
             OutJson["proc"][I].add(Process[ip].key[i],Process[ip].obj[i].key.size()==0 ? Process[ip].obj[i].toString() : Process[ip].obj[i].DumpToString());
          }

          //add the data block
          OutJson["proc"][I].addArray("data");
          int J=0;

          //check the samples in the data block and replace them by the reference
          std::vector<JSONWrapper::Object> Samples = (Process[ip])["data"].daughters();
          for(size_t id=0; id<Samples.size(); id++){
             OutJson["proc"][I]["data"].addList();

             if(!Samples[id].isTag("dtag")) continue;
             string dtag = Samples[id].getString("dtag");
             if(Datasets.find(dtag)==Datasets.end()){ //make sure that the sample is in the reference, if not add it
                printf("WARNING: dtag=%s is not found in the reference files --> the block will be copied from the current file\n", dtag.c_str());
                datasetinfo info(Samples[id], Process[ip].getBool("isdata", false) );
                Datasets[dtag] = info;
             }
             datasetinfo& inforef = Datasets[dtag];

             //copy the info from the ref
             if(inforef.prop.find("dtag" )!=inforef.prop.end())OutJson["proc"][I]["data"][J].add("dtag" , JSONWrapper::removeWhiteSpace(inforef.prop["dtag" ]), 50);
             if(inforef.prop.find("split")!=inforef.prop.end())OutJson["proc"][I]["data"][J].add("split", JSONWrapper::removeWhiteSpace(inforef.prop["split"]), 5);
             if(inforef.prop.find("xsec" )!=inforef.prop.end())OutJson["proc"][I]["data"][J].add("xsec" , JSONWrapper::removeWhiteSpace(inforef.prop["xsec" ]), 20);
             if(inforef.prop.find("br"   )!=inforef.prop.end())OutJson["proc"][I]["data"][J].add("br"   , JSONWrapper::removeWhiteSpace(inforef.prop["br"   ]), 30);
             if(inforef.prop.find("dset" )!=inforef.prop.end())OutJson["proc"][I]["data"][J].add("dset" , JSONWrapper::removeWhiteSpace(inforef.prop["dset" ]));
             for(std::map<string, string>::iterator itp = inforef.prop.begin(); itp != inforef.prop.end(); itp++){
                if(itp->first=="dtag" || itp->first=="isdata" || itp->first=="xsec" || itp->first=="br" || itp->first=="split" || itp->first=="dset")continue;
                OutJson["proc"][I]["data"][J].add(itp->first, JSONWrapper::removeWhiteSpace(itp->second));
             }
             J++;
          }
          I++;
      }

      //save the output json to file
      FILE* pFile = fopen((jsonFileToBeFixed).c_str(), "w");
      OutJson.Dump(pFile);
      fclose(pFile);
   }

}






int main(int argc, char* argv[]){
   bool multicrab   = false;
   bool multicrabMC = false;
   std::vector<string> inJsons;
   std::vector<string> fixJsons;
   string outJson = "";    

   std::vector<string> args;
   for(int i=1;i<argc;i++){
     args.push_back(argv[i]);
     string arg(argv[i]);
     if(arg.find("--help")!=string::npos){
        printf("--help               --> print this helping text\n");
        printf("--in  [in1,...,inN)] --> input  Json files \n");
        printf("--fix [fx1,...,fxN)] --> Json files to be fixed.  The input files are considered as the reference and the files to be fixed will be modified according to the reference input files\n");
        printf("--out [out]          --> output Json file \n");
        printf("--multicrab          --> create a multicrab.cfg file for the data samples \n");
        printf("--multicrabMC        --> create a multicrab.cfg file for the MC samples \n");
	return 0;
     }

     else if(arg.find("--in" )!=string::npos){ while(i+1<argc && string(argv[i+1]).find("--")!=0){inJsons.push_back(argv[i+1]);  i++;}  printf("input  Json files:\n"); for(unsigned int j=0;j<inJsons.size();j++){printf("\t - %s\n", inJsons[j].c_str());}  }
     else if(arg.find("--fix")!=string::npos){ while(i+1<argc && string(argv[i+1]).find("--")!=0){fixJsons.push_back(argv[i+1]);  i++;}  printf("to be fixed Json files:\n"); for(unsigned int j=0;j<fixJsons.size();j++){printf("\t - %s\n", fixJsons[j].c_str());}  }
     else if(arg.find("--out")!=string::npos && i+1<argc){outJson = argv[i+1];  i++;  printf("output Json file : %s\n", outJson.c_str());}
     else if(arg.find("--multicrabMC")!=string::npos){ multicrabMC = true; }
     else if(arg.find("--multicrab"  )!=string::npos){ multicrab   = true; }
   }


   JSONWrapper::Object AllJson = merge_json(inJsons);

   if(multicrab || multicrabMC){
      make_multicrab_cfg(AllJson, !multicrabMC);
   }

   for(unsigned int i=0;i<fixJsons.size();i++){
      fix_json(AllJson, fixJsons[i]);
   }

   if(outJson!=""){
      FILE* pFile = fopen(outJson.c_str(), "w"); 
      AllJson.Dump(pFile);
      fclose(pFile);
   }

   return 0;
}
