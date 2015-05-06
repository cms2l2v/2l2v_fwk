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

#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TF1.h"
#include "TCutG.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TPaveText.h"
#include "THStack.h"
#include "TGaxis.h"


#include "UserCode/llvv_fwk/interface/JSONWrapper.h"
#include "UserCode/llvv_fwk/interface/tdrstyle.h"


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

std::string get_FileName(std::string RootDir,
			 std::vector<JSONWrapper::Object> Samples,
			 int id,
			 int s) {
  // follow CRAB convention, always start with 1 
  std::string segmentExt;
  char buf[255];
  sprintf(buf,"_%i",s);
  segmentExt += buf;

  std::string FileName = RootDir 
    + Samples[id].getString("dtag", "") + "/"
    + "output" + Samples[id].getString("suffix","") + segmentExt + ".root";
  // std::cout << FileName << std::endl;
  
  return FileName; 
}

bool isFileExist(TFile* File){
  if(!File || File->IsZombie() || !File->IsOpen() ||
     File->TestBit(TFile::kRecovered) )
    return false;

  else
    return true;
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
	// std::cout << "Split = " << split << std::endl;

	// loop over all files, follow CRAB convention start with 1 
	for(int s=0; s<split; s++){
	  // follow CRAB convention start with 1
	  std::string FileName = get_FileName(RootDir, Samples, id, s+1); 

	  TFile* File = new TFile(FileName.c_str());
	  bool& fileExist = FileExist[FileName];

	  fileExist = isFileExist(File);
	  if ( !fileExist ) continue;
 
	  //do the following only for the first file
	  if(s>0) continue;

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


void checkSumw2(TH1 *h) { if(h==0) return;  if(h->GetDefaultSumw2()) h->Sumw2();  }


TH1* set_hist_style(JSONWrapper::Object proc, TH1* hist) {
  if(proc.isTag("color" ) )hist->SetLineColor  ((int)proc[ "color"].toDouble()); else hist->SetLineColor  (1);
  if(proc.isTag("color" ) )hist->SetMarkerColor((int)proc[ "color"].toDouble()); else hist->SetMarkerColor(1);
  if(proc.isTag("color" ) )hist->SetFillColor  ((int)proc[ "color"].toDouble()); else hist->SetFillColor  (0);
  if(proc.isTag("lcolor") )hist->SetLineColor  ((int)proc["lcolor"].toDouble());
  if(proc.isTag("mcolor") )hist->SetMarkerColor((int)proc["mcolor"].toDouble()); 
  if(proc.isTag("fcolor") )hist->SetFillColor  ((int)proc["fcolor"].toDouble()); 
  if(proc.isTag("lwidth") )hist->SetLineWidth  ((int)proc["lwidth"].toDouble());// else hist->SetLineWidth  (1);
  if(proc.isTag("lstyle") )hist->SetLineStyle  ((int)proc["lstyle"].toDouble());// else hist->SetLinStyle  (1);
  if(proc.isTag("fill"  ) )hist->SetFillColor  ((int)proc["fill"  ].toDouble());
  if(proc.isTag("marker") )hist->SetMarkerStyle((int)proc["marker"].toDouble());
  return hist;
}

void Draw1DHistogram(JSONWrapper::Object& Root,
		     std::string RootDir,
		     NameAndType HistoProperties,
		     std::string outDir){

  std::vector<TObject*> ObjectToDelete;
  THStack* stack = new THStack("MC","MC");
  TCanvas* c1 = new TCanvas("c1","c1",800,800);
  TLegend* legA  = new TLegend(0.7,0.8,0.99, 0.99, "NDC"); 

  std::string SaveName = "";
   
  std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
  // loop over procs 
  for(unsigned int i=0;i<Process.size();i++){
    TH1* proc_hist = NULL;
    std::vector<JSONWrapper::Object> Samples = (Process[i])["data"].daughters();
    // loop over samples 
    for(unsigned int j=0;j<Samples.size();j++){
      double Weight = 1.0; // use 1.0 for now. 

      int split = Samples[j].getInt("split", 1);
      TH1* samp_hist = NULL;
      int NFiles=0;
      // loop over files 
      for(int s=1;s<=split;s++){
	TH1* file_hist = NULL; 
	std::string FileName = get_FileName(RootDir, Samples, j, s);

	// std::cout << FileName << std::endl;
	TFile* File = new TFile(FileName.c_str());
	if ( !isFileExist(File) ) {delete File; continue;}
	// std::cout << ">>>>>" << HistoProperties.name << std::endl;
	file_hist = (TH1*) GetObjectFromPath(File, HistoProperties.name);  
	if(!file_hist) {delete File; continue;} 

	// std::cout << "Found hist" << file_hist << std::endl;
	NFiles++;
	if(!samp_hist) {
	  gROOT->cd();
	  samp_hist = (TH1*)file_hist->Clone(file_hist->GetName());
	  checkSumw2(samp_hist);
	}else 
	  samp_hist->Add(file_hist);

	delete file_hist;
	delete File;

	// std::cout << ">>> 2 " << std::endl;
      }// end files loop 

      if(!samp_hist) continue;
      // Need to check : 
      // if(!Process[i]["isdata"].toBool())
      //  tmphist->Scale(1.0/NFiles);      

      // std::cout << ">>> 3 " << std::endl;
	
      if(!proc_hist) {
	gROOT->cd();
	proc_hist = (TH1*)samp_hist->Clone(samp_hist->GetName());
	checkSumw2(proc_hist);
	proc_hist->Scale(Weight);
      } else
	proc_hist->Add(samp_hist, Weight);
      delete samp_hist;
      
    } // end samples loop
    // std::cout << ">>> 4 " << std::endl;
    if(!proc_hist) continue;
    SaveName = proc_hist->GetName();
    ObjectToDelete.push_back(proc_hist);

    proc_hist = set_hist_style(Process[i], proc_hist); 
    // std::cout << ">>> 5 " << std::endl;
    //Add to Stack
    stack->Add(proc_hist, "HIST");   
    TLegendEntry* le = legA->AddEntry(proc_hist, Process[i]["tag"].c_str(), "F");
    le->SetTextSize(0.02); 
    
    // std::cout << ">>> 6 " << std::endl;
    
  } // end procs loop 

  // if(! stack || stack->GetStack() || stack->GetStack()->GetEntriesFast()>0)
  //   return;


  std::string SavePath = SaveName;
  while(SavePath.find("*")!=std::string::npos)SavePath.replace(SavePath.find("*"),1,"");
  while(SavePath.find("#")!=std::string::npos)SavePath.replace(SavePath.find("#"),1,"");
   while(SavePath.find("{")!=std::string::npos)SavePath.replace(SavePath.find("{"),1,"");
   while(SavePath.find("}")!=std::string::npos)SavePath.replace(SavePath.find("}"),1,"");
   while(SavePath.find("(")!=std::string::npos)SavePath.replace(SavePath.find("("),1,"");
   while(SavePath.find(")")!=std::string::npos)SavePath.replace(SavePath.find(")"),1,"");
   while(SavePath.find("^")!=std::string::npos)SavePath.replace(SavePath.find("^"),1,"");
   while(SavePath.find("/")!=std::string::npos)SavePath.replace(SavePath.find("/"),1,"-");
   if(outDir.size()) SavePath = outDir +"/"+ SavePath;
 
  if(stack && stack->GetStack() && stack->GetStack()->GetEntriesFast()>0){

    // std::cout << ">>> before draw " << std::endl;
    
    stack->Draw("");

    TH1 *hist=(TH1*)stack->GetStack()->At(0);
    if(stack->GetXaxis()) {
      stack->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());
      stack->GetYaxis()->SetTitle(hist->GetYaxis()->GetTitle());
      stack->SetMinimum(hist->GetMinimum());

      // stack->SetMaximum(maximumFound);
      // stack->SetMinimum(SignalMin);
    }

    legA->SetFillColor(0);
    legA->SetFillStyle(0);
    legA->SetLineColor(0);
    legA->SetHeader("");
    legA->Draw("same");
    legA->SetTextFont(42);
 
    ObjectToDelete.push_back(stack);
    // std::cout << "stack draw done. " << std::endl;
    // c1->SaveAs("c1.pdf");
    // std::cout << "c1 saved" << std::endl;
    system(std::string(("rm -f ") + SavePath + ".pdf").c_str());
    c1->SaveAs((SavePath + ".pdf").c_str());
    
  
    delete c1;
  }
  for(unsigned int d=0;d<ObjectToDelete.size();d++){
    delete ObjectToDelete[d];
  }
  ObjectToDelete.clear();
  
}

void SavingToFile(JSONWrapper::Object& Root,
		  std::string RootDir,
		  NameAndType HistoProperties,
		  TFile* OutputFile){
  std::vector<TObject*> ObjectToDelete;
  std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();

  using std::string; 
  std::unordered_map<std::string, bool> FileExist;

  // printf("Inside SavingToFile\n");		   
  for(unsigned int i=0;i<Process.size();i++){
      TH1* hist = NULL;
      std::vector<JSONWrapper::Object> Samples = (Process[i])["data"].daughters();
      for(unsigned int j=0;j<Samples.size();j++){
         double Weight = 1.0;
         // if(!Process[i]["isdata"].toBool() && !Process[i]["isdatadriven"].toBool() && HistoProperties.name.find("optim_")==std::string::npos )Weight*= iLumi;
	 string filtExt("");
	 if(Process[i].isTag("mctruthmode") ) { char buf[255]; sprintf(buf,"_filt%d",(int)Process[i]["mctruthmode"].toInt()); filtExt += buf; }	 

         int split = Samples[j].getInt("split", 1);
	 // printf("split = %d\n", split);
	 
         TH1* tmphist = NULL;  int NFiles=0;
         for(int s=0;s<split;s++){ 
	   //printf("Running on s = %d\n", s);
	   //string segmentExt; if(split>1){ char buf[255]; sprintf(buf,"_%i",s); segmentExt += buf;}
	   
           //string FileName = RootDir + (Samples[j])["dtag"].toString() + Samples[j].getString("suffix", "") +  segmentExt + filtExt + ".root";
	   std::string FileName = get_FileName(RootDir, Samples, j, s+1); // crab default 1 
	   //printf("FileName = %s", FileName.c_str());		   
           // if(!FileExist[FileName])continue;
	   TFile* File = new TFile(FileName.c_str());
           // if(!File || File->IsZombie() || !File->IsOpen() || File->TestBit(TFile::kRecovered) )continue;
	   if ( !isFileExist(File) ) {delete File; continue;}


           TH1* tmptmphist = (TH1*) GetObjectFromPath(File,HistoProperties.name); 
           if(!tmptmphist){delete File;continue;}            
	   // found existing histo: 
	   // std::cout << "Found the histo: " << HistoProperties.name << std::endl;
	   if(!tmphist){gROOT->cd(); tmphist = (TH1*)tmptmphist->Clone(tmptmphist->GetName()); checkSumw2(tmphist);}else{tmphist->Add(tmptmphist);}
           NFiles++;

	   // printf("NFiles = %i", NFiles);
	   
	   delete tmptmphist;
           delete File;
         }
         if(!tmphist)continue;
         if(!Process[i]["isdata"].toBool() && HistoProperties.name.find("optim_")==std::string::npos)tmphist->Scale(1.0/NFiles);
         if(!hist){gROOT->cd(); hist = (TH1*)tmphist->Clone(tmphist->GetName());checkSumw2(hist);hist->Scale(Weight);}else{hist->Add(tmphist,Weight);}
         delete tmphist;
      }   
      if(!hist)continue;      

      string dirName = Process[i]["tag"].c_str();while(dirName.find("/")!=std::string::npos)dirName.replace(dirName.find("/"),1,"-");
      OutputFile->cd();
      TDirectory* subdir = OutputFile->GetDirectory(dirName.c_str());
      if(!subdir || subdir==OutputFile) subdir = OutputFile->mkdir(dirName.c_str());
      subdir->cd();
      hist->Write();
      delete hist;
   }
}



void runPlotter(std::string inDir,
		std::string jsonFile,
		std::string outDir,
		std::string outFile)
{
  setTDRStyle();  
  // gStyle->SetPadTopMargin   (0.06);
  // gStyle->SetPadBottomMargin(0.12);
  // gStyle->SetPadRightMargin (0.16);
  // gStyle->SetPadLeftMargin  (0.14);
  gStyle->SetTitleSize(0.03, "XYZ");
  gStyle->SetLabelSize(.03, "XY");
  gStyle->SetTitleXOffset(1.1);
  gStyle->SetTitleYOffset(2.1);
  // gStyle->SetPalette(1);
  // gStyle->SetNdivisions(505);
  // TGaxis::SetMaxDigits(4);
 

  JSONWrapper::Object Root(jsonFile, true);
  std::list<NameAndType> histlist;
  GetListOfObject(Root, inDir, histlist);
  std::cout << "Total of hists found: " << histlist.size() << std::endl;
  histlist.sort();
  histlist.unique();   

  TFile* OutputFile = new TFile(outFile.c_str(),"RECREATE");
  // printf("Progressing Bar              :0%%       20%%       40%%       60%%       80%%       100%%\n");
  // printf("Processing :");
  int TreeStep = 1; // histlist.size();
  // if (TreeStep == 0) TreeStep = 1;
  
  int ictr(0);
  for(std::list<NameAndType>::iterator it= histlist.begin();
      it!= histlist.end(); it++, ictr++){
    // std::cout << "ictr = " << ictr << std::endl;
    // if (ictr > 1  ) break;
    /// std::cout << "Processing name: " << (*it).name  << std::endl;    
    // if(ictr%TreeStep==0){printf(".");fflush(stdout);}
    // if(ictr%TreeStep==0){printf("%d.", ictr);fflush(stdout);}
    if( it->is1D() ){
      // std::cout << "is 1D" << std::endl;
      // Draw1DHistogram(Root, inDir, *it, outDir);
      // std::cout << "Processing name: " << (*it).name  << std::endl;    
      SavingToFile(Root, inDir, *it, OutputFile);
    }
  }
  printf("Saved as: %s\n", outFile.c_str());
  OutputFile->Close();
  
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
  runPlotter(inDir, jsonFile, outDir, outFile); 

}  







