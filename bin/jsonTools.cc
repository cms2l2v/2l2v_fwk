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


bool compareJSON(JSONWrapper::Object a, JSONWrapper::Object b, std::string indent=""){
   indent = indent + " - ";
   printf("%sVAL: %s vs %s", indent.c_str(),  a.val.c_str(), b.val.c_str());
   if(a.val != b.val){printf("  pasOk\n");return false;}else{ printf("\n");}
   printf("%sKEYS: %i vs %i", indent.c_str(), (int)a.key.size(), (int)b.key.size());
   if(a.key.size()!=b.key.size()){printf("  pasOk\n");return false;}else{printf("\n");}
   for(unsigned int i=0;i<a.key.size();i++){
printf("%sXXX %i\n", indent.c_str(), i);
             
       printf("%sKEY  %s vs %s", indent.c_str(), a.key[i].c_str(), b.key[i].c_str());       
       if(std::find(b.key.begin(), b.key.end(), a.key[i])==b.key.end()){printf("  pasOk\n");return false;}else{ printf("\n");}
       if(!compareJSON(a.obj[i], b[a.key[i].c_str()], indent))return false;
   }

   return true;
}


/*
JSONWrapper::Object merge_json(std::vector<std::string> jsonFiles){
   JSONWrapper::Object* Roots = new JSONWrapper::Object[jsonFiles.size()];
   for(unsigned int J=0;J<jsonFiles.size();J++){
      Roots[J] = JSONWrapper::Object(jsonFiles[J], true);

      if()

   }

  if(jsonFiles.size()>1){
     

     printf("%s and %s comparison = %i\n", jsonFiles[0].c_str(), jsonFiles[1].c_str(), compareJSON(Roots[0], Roots[1])==true?1:0);
  }


   return Roots[0];
}
*/


JSONWrapper::Object merge_json(std::vector<std::string> jsonFiles){
   struct datasetinfo{std::map<string, string> prop;
      datasetinfo(){};
      datasetinfo(JSONWrapper::Object dtagObj){
         for(unsigned int i=0;i<dtagObj.key.size();i++){
            prop[dtagObj.key[i]] = dtagObj.obj[i].key.size()==0 ? dtagObj.obj[i].toString() : dtagObj.obj[i].DumpToString();
         }
         if(prop["dset" ]=="")prop["dset" ]="Unknown";
         if(prop["split"]=="")prop["split"]="1";
         if(prop["br"   ]=="")prop["br"   ]="[1.0]";
         if(prop["xsec" ]=="")prop["xsec" ]="0.0";
      };
   }; 
   std::map<std::string, datasetinfo> Datasets;

   for(unsigned int J=0;J<jsonFiles.size();J++){
      JSONWrapper::Object Root(jsonFiles[J], true);
      std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
      for(size_t ip=0; ip<Process.size(); ip++){
          std::vector<JSONWrapper::Object> Samples = (Process[ip])["data"].daughters();
          for(size_t id=0; id<Samples.size(); id++){
             if(!Samples[id].isTag("dtag")) continue;
               string dtag = Samples[id].getString("dtag");

               datasetinfo infodef;
               datasetinfo info(Samples[id]);
               
               if(Datasets.find(dtag)!=Datasets.end()){
                  //check if they are the same
                  datasetinfo& infosaved = Datasets[dtag];
                  for(std::map<string, string>::iterator it = infosaved.prop.begin(); it != infosaved.prop.end(); it++){
                     string k = it->first;
                     if(infosaved.prop[k]  != info.prop[k] ){printf("DIFF in file=%35s dtag=%35s var=%10s : %30s vs %30s  Please fix the files\n", jsonFiles[J].c_str(), dtag.c_str(), k.c_str(), infosaved.prop[k].c_str(), info.prop[k].c_str());} 
                  }

                  continue;//dataset already exist
                  //need to add a function to check if they are the same or not and print a warning if they are not
               }else{
//                  for(unsigned int k=0;k<Samples[id].key.size();k++){printf("%30s / %30s --> %30s - %30s\n", Process[ip].getString("tag").c_str(), Samples[id].getString("dtag").c_str(), Samples[id].key[k].c_str(), Samples[id].getString(Samples[id].key[k]).c_str());}
                  Datasets[dtag] = info;
               }
          }
      }
   }


   JSONWrapper::Object AllJson;

//   AllJson.addArray("proc");
//   for(std::map<std::string, datasetinfo>::iterator it = Datasets.begin(); it!=Datasets.end();it++){
//          int I = AllJson["proc"].daughters().size();
//          AllJson["proc"].addList(); 
//          AllJson["proc"][I].add("tag", it->first);      
//          if(it->second.isdata)   AllJson["proc"][I].add("isdata","true");
//          AllJson["proc"][I].addArray("data");
//          AllJson["proc"][I]["data"].addList();
//          if(it->first       !="")AllJson["proc"][I]["data"][0].add("dtag", it->first);      
//          if(it->second.dset !="")AllJson["proc"][I]["data"][0].add("dset", it->second.dset);
//          if(it->second.split!="")AllJson["proc"][I]["data"][0].add("split",it->second.split);
//          if(it->second.xsec !="")AllJson["proc"][I]["data"][0].add("xsec" ,it->second.xsec);
//          if(it->second.br   !="")AllJson["proc"][I]["data"][0].add("br"   ,it->second.br);
//   }

   return AllJson;
}


/*
JSONWrapper::Object merge_json(std::vector<std::string> jsonFiles){
   struct datasetinfo{std::string dset; std::string split; bool isdata; string br; string xsec;
   datasetinfo(){dset="Unknown", split="1"; isdata=false; br="[1.0]"; xsec="0.0";};
   }; 
   std::map<std::string, datasetinfo> Datasets;

   for(unsigned int J=0;J<jsonFiles.size();J++){
      JSONWrapper::Object Root(jsonFiles[J], true);
      std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
      for(size_t ip=0; ip<Process.size(); ip++){
          std::vector<JSONWrapper::Object> Samples = (Process[ip])["data"].daughters();
          for(size_t id=0; id<Samples.size(); id++){
             if(!Samples[id].isTag("dtag")) continue;
               string dtag = Samples[id].getString("dtag");

               datasetinfo infodef;
               datasetinfo info;
               info.dset   = Samples[id].getString("dset", infodef.dset);
               info.split  = Samples[id].getString("split", infodef.split);                  
               info.br     = Samples[id].getFullString("br",infodef.br);
               info.xsec   = Samples[id].getString("xsec" , infodef.xsec);
               info.isdata = Process[ip].getBool  ("isdata", infodef.isdata); 
               
               if(Datasets.find(dtag)!=Datasets.end()){
                  //check if they are the same
                  datasetinfo& infosaved = Datasets[dtag];
                  if(infosaved.dset  != info.dset ){if(infosaved.dset  == infodef.dset ){infosaved.dset  = info.dset ;}else{printf("DIFF dset  in %s dtag %s : %s vs %s  Please fix the files\n", jsonFiles[J].c_str(), dtag.c_str(), infosaved.dset .c_str(), info.dset .c_str());} }
                  if(infosaved.split != info.split){if(infosaved.split == infodef.split){infosaved.split = info.split;}else{printf("DIFF split in %s dtag %s : %s vs %s  Please fix the files\n", jsonFiles[J].c_str(), dtag.c_str(), infosaved.split.c_str(), info.split.c_str());} }

                  continue;//dataset already exist
                  //need to add a function to check if they are the same or not and print a warning if they are not
               }else{
//                  for(unsigned int k=0;k<Samples[id].key.size();k++){printf("%30s / %30s --> %30s - %30s\n", Process[ip].getString("tag").c_str(), Samples[id].getString("dtag").c_str(), Samples[id].key[k].c_str(), Samples[id].getString(Samples[id].key[k]).c_str());}
                  Datasets[dtag] = info;
               }
          }
      }
   }


   JSONWrapper::Object AllJson;

//   AllJson.addArray("proc");
//   for(std::map<std::string, datasetinfo>::iterator it = Datasets.begin(); it!=Datasets.end();it++){
//          int I = AllJson["proc"].daughters().size();
//          AllJson["proc"].addList(); 
//          AllJson["proc"][I].add("tag", it->first);      
//          if(it->second.isdata)   AllJson["proc"][I].add("isdata","true");
//          AllJson["proc"][I].addArray("data");
//          AllJson["proc"][I]["data"].addList();
//          if(it->first       !="")AllJson["proc"][I]["data"][0].add("dtag", it->first);      
//          if(it->second.dset !="")AllJson["proc"][I]["data"][0].add("dset", it->second.dset);
//          if(it->second.split!="")AllJson["proc"][I]["data"][0].add("split",it->second.split);
//          if(it->second.xsec !="")AllJson["proc"][I]["data"][0].add("xsec" ,it->second.xsec);
//          if(it->second.br   !="")AllJson["proc"][I]["data"][0].add("br"   ,it->second.br);
//   }

   return AllJson;
}
*/

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
