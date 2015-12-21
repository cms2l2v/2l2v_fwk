#include <vector>
#include <list>
#include <iterator>
#include <algorithm>
#include <unordered_map>
#include <iostream>
#include <string>
#include <regex>

 
#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TText.h"
#include "TF1.h"
#include "TCutG.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TPaveText.h"
#include "THStack.h"

#include "UserCode/llvv_fwk/interface/tdrstyle.h"
#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/JSONWrapper.h"
#include "HiggsAnalysis/CombinedLimit/interface/th1fmorph.h"

using namespace std;

int cutIndex=-1;
string cutIndexStr="";
double iLumi = 2007;
double iEcm=8;
bool showChi2 = false;
bool showUnc=false;
double baseRelUnc=0.044;
bool noLog=false;
bool isSim=false;
bool doTree = true;
bool do2D  = true;
bool do1D  = true;
bool doTex = true;
bool doPowers = true;
bool StoreInFile = true;
bool doPlot = true;
bool splitCanvas = false;
bool onlyCutIndex = false;
bool showRatioBox=true;
bool fixUnderflow=true;
bool fixOverflow=true;
double blind=-1E99;
string inDir   = "OUTNew/";
string jsonFile = "../../data/beauty-samples.json";
string outDir  = "Img/";
std::vector<std::string> plotExt;
string outFile = "/tmp/plotter.root";
string cutflowhisto = "all_cutflow";
string fileOption = "RECREATE";

//std::unordered_map<string, stSampleInfo> sampleInfoMap;
std::unordered_map<string, std::vector<string> > MissingFiles;
std::unordered_map<string, std::vector<string> > DSetFiles;


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

void checkSumw2(TH1 *h) { if(h==0) return;  if(h->GetDefaultSumw2()) h->Sumw2();  }


string dropBadCharacters(string in){
   while(in.find("*")!=std::string::npos)in.replace(in.find("*"),1,"");
   while(in.find("#")!=std::string::npos)in.replace(in.find("#"),1,"");
   while(in.find("{")!=std::string::npos)in.replace(in.find("{"),1,"");
   while(in.find("}")!=std::string::npos)in.replace(in.find("}"),1,"");
   while(in.find("(")!=std::string::npos)in.replace(in.find("("),1,"");
   while(in.find(")")!=std::string::npos)in.replace(in.find(")"),1,"");
   while(in.find("^")!=std::string::npos)in.replace(in.find("^"),1,"");
   while(in.find("/")!=std::string::npos)in.replace(in.find("/"),1,"-");
   return in;
}


void setStyle(JSONWrapper::Object& SingleProcess, TH1* hist){
   if(SingleProcess.isTag("color" ) )hist->SetLineColor  ((int)SingleProcess[ "color"].toDouble()); else hist->SetLineColor  (1);
   if(SingleProcess.isTag("color" ) )hist->SetMarkerColor((int)SingleProcess[ "color"].toDouble()); else hist->SetMarkerColor(1);
   if(SingleProcess.isTag("color" ) )hist->SetFillColor  ((int)SingleProcess[ "color"].toDouble()); else hist->SetFillColor  (0);
   if(SingleProcess.isTag("lcolor") )hist->SetLineColor  ((int)SingleProcess["lcolor"].toDouble());
   if(SingleProcess.isTag("mcolor") )hist->SetMarkerColor((int)SingleProcess["mcolor"].toDouble()); 
   if(SingleProcess.isTag("fcolor") )hist->SetFillColor  ((int)SingleProcess["fcolor"].toDouble()); 
   if(SingleProcess.isTag("lwidth") )hist->SetLineWidth  ((int)SingleProcess["lwidth"].toDouble());// else hist->SetLineWidth  (1);
   if(SingleProcess.isTag("lstyle") )hist->SetLineStyle  ((int)SingleProcess["lstyle"].toDouble());// else hist->SetLinStyle  (1);
   if(SingleProcess.isTag("fill"  ) )hist->SetFillColor  ((int)SingleProcess["fill"  ].toDouble());
   if(SingleProcess.isTag("marker") )hist->SetMarkerStyle((int)SingleProcess["marker"].toDouble());// else hist->SetMarkerStyle(1);
   if(SingleProcess.isTag("msize")  )hist->SetMarkerSize (      SingleProcess["msize"].toDouble());// else the Size is the defaoult one
}

double getXsecXbr(JSONWrapper::Object& SingleProcess){
   double toReturn = 1.0;
   if(!SingleProcess.isTag("data"))return toReturn;
   toReturn*= SingleProcess["data"][0].getDouble("xsec",1.0);
   if(SingleProcess["data"][0].isTag("br")){
      std::vector<JSONWrapper::Object> BRs = SingleProcess["data"][0]["br"].daughters();
      for(size_t ipbr=0; ipbr<BRs.size(); ipbr++){toReturn*=BRs[ipbr].toDouble();}
   }
   return toReturn;
}



// function that add the TPaveText on the current canvas with the "CMS Preliminary...." on top of the Histograms. For split Lumi
void DrawPreliminary(double L, double B, double R, double T,  bool isSim=false, bool preliminary=true){
     //TOP RIGHT OUT-FRAME
     char LumiText[1024];  sprintf(LumiText, "%.1f %s^{-1} (%.0f TeV)", iLumi>100?iLumi/1000:iLumi, iLumi>100?"fb":"pb", iEcm);       
     TPaveText* T1 = new TPaveText(1.0-R-0.50, 1.0-T, 1.02-R, 1.0-0.005, "NDC");
     T1->SetTextFont(43); T1->SetTextSize(23);   T1->SetTextAlign(31);
     T1->SetFillColor(0); T1->SetFillStyle(0);   T1->SetBorderSize(0);
     T1->AddText(LumiText);  T1->Draw();

     //TOP LEFT IN-FRAME
     TPaveText* T2 = new TPaveText(L+0.005, 1.0-T-0.05, L+0.20, 1.0-T-0.005, "NDC");
     T2->SetTextFont(63); T2->SetTextSize(30);   T2->SetTextAlign(11);
     T2->SetFillColor(0); T2->SetFillStyle(0);   T2->SetBorderSize(0);
     T2->AddText("CMS"); T2->Draw();

     if(preliminary){ //Bellow CMS
     TPaveText* T3 = new TPaveText(L+0.005, 1.0-T-0.085, L+0.20, 1.0-T-0.035, "NDC");
     T3->SetTextFont(53); T3->SetTextSize(23);   T3->SetTextAlign(11);
     T3->SetFillColor(0); T3->SetFillStyle(0);   T3->SetBorderSize(0);
     T3->AddText("Preliminary"); T3->Draw();
     }

     //if(Text!=""){  //TOP right IN-FRAME
     //TPaveText* T4 = new TPaveText(1.0-R-0.50, 1.0-T-0.06, 1.0-R, 1.0-T-0.01, "NDC");
     //T4->SetTextFont(43); T4->SetTextSize(23);   T4->SetTextAlign(32);
     //T4->SetFillColor(0); T4->SetFillStyle(0);   T4->SetBorderSize(0);
     //T4->AddText(Text.c_str());  T4->Draw();
     //}

     if(isSim){ //Right of CMS
     TPaveText* T5 = new TPaveText(L+0.085, 1.0-T-0.05, L+0.32, 1.0-T-0.005, "NDC");
     T5->SetTextFont(53); T5->SetTextSize(23);   T5->SetTextAlign(11);
     T5->SetFillColor(0); T5->SetFillStyle(0);   T5->SetBorderSize(0);
     T5->AddText("Simulation"); T5->Draw();
     }
}

void DrawPreliminary(TAttPad* pad=NULL, bool isSim=false,   bool preliminary=true){
   if     (pad )DrawPreliminary(pad->GetLeftMargin(), pad->GetBottomMargin(), pad->GetRightMargin(), pad->GetTopMargin() , isSim, preliminary);
   else if(gPad)DrawPreliminary(gPad->GetLeftMargin(),gPad->GetBottomMargin(),gPad->GetRightMargin(),gPad->GetTopMargin(), isSim, preliminary);
   else         DrawPreliminary(0.15, 0.15, 0.15, 0.15                                                                   , isSim, preliminary);
}

TObject* GetObjectFromPath(TDirectory* File, std::string Path, bool GetACopy=false)
{
   size_t pos = Path.find("/");
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

void fixExtremities(TH1* h,bool addOverflow, bool addUnderflow)
{
  if(h==0) return;

  if(addUnderflow){
      double fbin  = h->GetBinContent(0) + h->GetBinContent(1);
      double fbine = sqrt(h->GetBinError(0)*h->GetBinError(0)
                          + h->GetBinError(1)*h->GetBinError(1));
      h->SetBinContent(1,fbin);
      h->SetBinError(1,fbine);
      h->SetBinContent(0,0);
      h->SetBinError(0,0);
    }
  
  if(addOverflow){  
      int nbins = h->GetNbinsX();
      double fbin  = h->GetBinContent(nbins) + h->GetBinContent(nbins+1);
      double fbine = sqrt(h->GetBinError(nbins)*h->GetBinError(nbins) 
                          + h->GetBinError(nbins+1)*h->GetBinError(nbins+1));
      h->SetBinContent(nbins,fbin);
      h->SetBinError(nbins,fbine);
      h->SetBinContent(nbins+1,0);
      h->SetBinError(nbins+1,0);
    }
}


void GetListOfObject(JSONWrapper::Object& Root, std::string RootDir, std::list<NameAndType>& histlist, TDirectory* dir=NULL, std::string parentPath=""){
  if(parentPath=="" && !dir){
      int dataProcessed = 0;
      int signProcessed = 0;
      int bckgProcessed = 0; 

      std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
      //loop on all Proc
      for(size_t ip=0; ip<Process.size(); ip++){
          if(Process[ip].isTag("interpollation"))continue; //nothing to do with these
	  bool isData (  Process[ip]["isdata"].toBool()  );
          bool isSign ( !isData &&  Process[ip].isTag("spimpose") && Process[ip]["spimpose"].toBool());
  	  bool isMC   = !isData && !isSign; 
	  string filtExt("");
	  if(Process[ip].isTag("mctruthmode") ) { char buf[255]; sprintf(buf,"_filt%d",(int)Process[ip]["mctruthmode"].toInt()); filtExt += buf; }

	  std::vector<JSONWrapper::Object> Samples = (Process[ip])["data"].daughters();
	  for(size_t id=0; id<Samples.size(); id++){
	      int split = Samples[id].getInt("split", -1);               
              for(int s=0; s<(split>0?split:999); s++){
                 char buf[255]; sprintf(buf,"_%i",s); string segmentExt = buf;                 
                 string FileName = RootDir + Samples[id].getString("dtag", "") +  Samples[id].getString("suffix","") + segmentExt + filtExt; 
                 if(split<0){ //autosplitting --> check if there is a cfg file before checking if there is a .root file
                    FILE* pFile = fopen((FileName+"_cfg.py").c_str(), "r");
                    if(!pFile){break;}else{fclose(pFile);}
                 }
                 FileName += ".root";
	         TFile* File = new TFile(FileName.c_str());                                 
                 if(!File || File->IsZombie() || !File->IsOpen() || File->TestBit(TFile::kRecovered) ){
                    MissingFiles[Samples[id].getString("dtag")].push_back(FileName);
                    continue; 
                 }else{
                    DSetFiles[Samples[id].getString("dtag")].push_back(FileName);
                 }

                 //just to make it faster, only consider the first 3 sample of a same kind
                 if(isData){if(dataProcessed>=5){ File->Close(); continue;}else{dataProcessed++;}}
                 if(isSign){if(signProcessed>=5){ File->Close(); continue;}else{signProcessed++;}}
                 if(isMC  ){if(bckgProcessed>=5){ File->Close(); continue;}else{bckgProcessed++;}}

                 printf("Adding all objects from %25s to the list of considered objects\n",  FileName.c_str());
	         GetListOfObject(Root,RootDir,histlist,(TDirectory*)File,"" );
	         File->Close();
               }
	    }
      }

      printf("The list of missing or corrupted files, that are ignored, can be found below:\n");
      for(std::unordered_map<string, std::vector<string> >::iterator it = MissingFiles.begin(); it!=MissingFiles.end(); it++){
         if(it->second.size()<=0)continue;
         printf("Missing file in dataset %s:\n", it->first.c_str());
         for(unsigned int f=0;f<it->second.size();f++){
            printf("\t %s\n", it->second[f].c_str());
         }
      }

      //for(std::list<NameAndType>::iterator it= histlist.begin(); it!= histlist.end(); it++){printf("%s\n",it->name.c_str()); }
      return;
   }

   if(dir==NULL)return;

   TList* list = dir->GetListOfKeys();
   for(int i=0;i<list->GetSize();i++){
      TObject* tmp = GetObjectFromPath(dir,list->At(i)->GetName(),false);

      if(tmp->InheritsFrom("TDirectory")){
         GetListOfObject(Root,RootDir,histlist,(TDirectory*)tmp,parentPath+ list->At(i)->GetName()+"/" );
      }else if(tmp->InheritsFrom("TTree")){ 
        printf("found one object inheriting from a ttree\n");
        histlist.push_back(NameAndType(parentPath+list->At(i)->GetName(), 4, false ) );
      }else if(tmp->InheritsFrom("TH1")){
        int  type = 0;
        if(tmp->InheritsFrom("TH1")) type++;
        if(tmp->InheritsFrom("TH2")) type++;
        if(tmp->InheritsFrom("TH3")) type++;
        bool hasIndex = string(((TH1*)tmp)->GetXaxis()->GetTitle()).find("cut index")<string::npos;
        if(hasIndex){type=1;}
	histlist.push_back(NameAndType(parentPath+list->At(i)->GetName(), type, hasIndex ) );
      }else{
        printf("The file contain an unknown object named %s\n", list->At(i)->GetName() );
      }
      delete tmp;
   }
}



void InterpollateProcess(JSONWrapper::Object& Root, TFile* File, NameAndType& HistoProperties){
    
   std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
   for(unsigned int i=0;i<Process.size();i++){
      if(!Process[i].isTag("interpollation"))continue; //only consider intepollated processes

      string dirName = Process[i]["tag"].c_str();while(dirName.find("/")!=std::string::npos)dirName.replace(dirName.find("/"),1,"-");
      File->cd();
      TDirectory* subdir = File->GetDirectory(dirName.c_str());
      if(!subdir || subdir==File) subdir = File->mkdir(dirName.c_str());
      subdir->cd();
     
      string signal   = Process[i]["tag"].c_str();
      string signalL  = Process[i]["interpollation"][0]["tagLeft"].c_str();
      string signalR  = Process[i]["interpollation"][0]["tagRight"].c_str();
      double mass     = Process[i]["interpollation"][0]["mass"].toDouble();
      double massL    = Process[i]["interpollation"][0]["massLeft"].toDouble();
      double massR    = Process[i]["interpollation"][0]["massRight"].toDouble();
      double Ratio = ((double)mass - massL); Ratio/=(massR - massL);

      double xsecXbr  = getXsecXbr(Process[i]);
      double xsecXbrL = xsecXbr;  for(unsigned int j=0;j<Process.size();j++){if(Process[j]["tag"].c_str()==signalL){xsecXbrL = getXsecXbr(Process[j]); break;}}
      double xsecXbrR = xsecXbr;  for(unsigned int j=0;j<Process.size();j++){if(Process[j]["tag"].c_str()==signalR){xsecXbrR = getXsecXbr(Process[j]); break;}}

      while(signalL.find("/")!=std::string::npos)signalL.replace(signalL.find("/"),1,"-");
      while(signalR.find("/")!=std::string::npos)signalR.replace(signalR.find("/"),1,"-");

      TH1* objL = (TH1*)GetObjectFromPath(File,signalL + "/" + HistoProperties.name);      
      TH1* objR = (TH1*)GetObjectFromPath(File,signalR + "/" + HistoProperties.name);      
      if(!objL || !objR)continue;

      TH1* histoInterpolated = NULL;
      if(HistoProperties.isIndexPlot){
         TH2F* histo2DL = (TH2F*) objL;
         TH2F* histo2DR = (TH2F*) objR;
         if(!histo2DL or !histo2DR)continue;
         TH2F* histo2D  = (TH2F*) histo2DL->Clone(histo2DL->GetName());
         histo2D->Reset();

         for(int cutIndex=0;cutIndex<=histo2DL->GetNbinsX()+1;cutIndex++){
            TH1D* histoL = histo2DL->ProjectionY("tempL", cutIndex, cutIndex);
            TH1D* histoR = histo2DR->ProjectionY("tempR", cutIndex, cutIndex);
            if(histoL->Integral()>0 && histoR->Integral()>0){
               histoL->Scale(1.0/xsecXbrL);
               histoR->Scale(1.0/xsecXbrR);
               TH1D* histo  = (TH1D*)histoL->Clone(HistoProperties.name.c_str());  histo->Reset();
//             printf("%f - %f -%f    --> %f - %f -%f\n", massL, mass, massR, histoL->Integral(), (1-Ratio)*histoL->Integral() + Ratio*histoR->Integral(), histoR->Integral() );
               histo->Add(th1fmorph("interpolTemp","interpolTemp", histoL, histoR, massL, massR, mass, (1-Ratio)*histoL->Integral() + Ratio*histoR->Integral(), 0), 1);
               histo->Scale(xsecXbr);
               for(int y=0;y<=histo2DL->GetNbinsY()+1;y++){
                 histo2D->SetBinContent(cutIndex, y, histo->GetBinContent(y));
                 histo2D->SetBinError(cutIndex, y, histo->GetBinError(y));
               }
               delete histo;
            }
            delete histoR;
            delete histoL;
         }
         delete histo2DL;
         delete histo2DR;
         histoInterpolated = histo2D;
      }else if(HistoProperties.is1D()){
         TH1F* histoL = (TH1F*) objL;
         TH1F* histoR = (TH1F*) objR;
         if(!histoL or !histoR or histoL->Integral()<=0 or histoR->Integral()<=0)continue;
         if(histoL->Integral()>0 && histoR->Integral()>0){            
            histoL->Scale(1.0/xsecXbrL);
            histoR->Scale(1.0/xsecXbrR);
            double Integral = (1-Ratio)*histoL->Integral() + Ratio*histoR->Integral();
            TH1F* histo  = (TH1F*)histoL->Clone(HistoProperties.name.c_str());  histo->Reset();
            histo->Add(th1fmorph("interpolTemp","interpolTemp", histoL, histoR, massL, massR, mass, Integral, 0), 1.0);
//          printf("%40s - %f - %f -%f    --> %f - %f -%f\n", HistoProperties.name.c_str(), massL, mass, massR, histoL->Integral(), histo->Integral(), histoR->Integral() );
            histo->Scale(xsecXbr);          
            histoInterpolated = histo;
         }
         delete histoR;
         delete histoL;
      }

      if(histoInterpolated){
         subdir->cd();
         histoInterpolated->Write(HistoProperties.name.c_str());
         gROOT->cd();
         delete histoInterpolated;
      }
   }
}


void SavingToFile(JSONWrapper::Object& Root, std::string RootDir, TFile* OutputFile, std::list<NameAndType>& histlist){
   std::vector<TObject*> ObjectToDelete;
   std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();

   int IndexFiles = 0;
   int NFilesStep = 0;   for(std::unordered_map<string, std::vector<string> >::iterator it = DSetFiles.begin(); it!=DSetFiles.end(); it++){NFilesStep+=it->second.size();} NFilesStep=std::max(1, NFilesStep/50);
   printf("Processing input root files  :");

   for(unsigned int i=0;i<Process.size();i++){
      if(Process[i].isTag("interpollation"))continue; //treated in a specific function after the loop on Process
      string dirName = Process[i]["tag"].c_str();while(dirName.find("/")!=std::string::npos)dirName.replace(dirName.find("/"),1,"-");
      OutputFile->cd();
      TDirectory* subdir = OutputFile->GetDirectory(dirName.c_str());
      if(!subdir || subdir==OutputFile) subdir = OutputFile->mkdir(dirName.c_str());
      subdir->cd();

      float  Weight = 1.0;     
      std::vector<JSONWrapper::Object> Samples = (Process[i])["data"].daughters();
      for(unsigned int j=0;j<Samples.size();j++){
         std::vector<string>& fileList = DSetFiles[(Samples[j])["dtag"].toString()];
         if(!Process[i]["isdata"].toBool() && !Process[i]["isdatadriven"].toBool()){Weight= iLumi/fileList.size();}else{Weight=1.0;}

         for(int f=0;f<fileList.size();f++){
           if(IndexFiles%NFilesStep==0){printf(".");fflush(stdout);} IndexFiles++;
           TFile* File = new TFile(fileList[f].c_str());

           for(std::list<NameAndType>::iterator it= histlist.begin(); it!= histlist.end(); it++){
              NameAndType& HistoProperties = *it;

              TObject* inobj  = GetObjectFromPath(File,HistoProperties.name);  if(!inobj)continue;
              TObject* outobj = GetObjectFromPath(subdir,HistoProperties.name);
             
              if(HistoProperties.isTree()){
                 TTree* outtree = (TTree*)outobj;
                 if(!outtree){ //tree not yet in file, so need to add it
                    subdir->cd();
                    outtree =  ((TTree*)inobj)->CloneTree(-1, "fast");
                    outtree->SetDirectory(subdir);

                    TTree* weightTree = new TTree((HistoProperties.name+"_PWeight").c_str(),"plotterWeight");
                    weightTree->Branch("plotterWeight",&weightTree,"plotterWeight/F");
                    weightTree->SetDirectory(subdir);
                    for(unsigned int i=0;i<((TTree*)inobj)->GetEntries();i++){weightTree->Fill();}                       
                 }else{
                    outtree->CopyEntries(((TTree*)inobj), -1, "fast");
                    TTree* weightTree = (TTree*)GetObjectFromPath(subdir,HistoProperties.name+"_PWeight");
                    for(unsigned int i=0;i<((TTree*)inobj)->GetEntries();i++){weightTree->Fill();} 
                 }
              }else{        
                 if(HistoProperties.name.find("optim_")==std::string::npos) ((TH1*)inobj)->Scale(Weight);
                 TH1* outhist = (TH1*)outobj;
                 if(!outhist){ //histogram not yet in file, so need to add it
                    subdir->cd();
                    outhist = (TH1*)(inobj)->Clone(inobj->GetName());
                    setStyle(Process[i], outhist);
                    checkSumw2(outhist);
                 }else{
                    outhist->Add((TH1*)inobj);
                 }
              }
           }
           delete File;         
         }   
      }
      subdir->cd();
      subdir->Write();
   }printf("\n");


   int ictr = 0;
   int TreeStep = std::max(1,(int)(histlist.size()/50));
   printf("Sample Interpolation         :");
   for(std::list<NameAndType>::iterator it= histlist.begin(); it!= histlist.end(); it++,ictr++){
       if(ictr%TreeStep==0){printf(".");fflush(stdout);}
      NameAndType& HistoProperties = *it;
      if(!HistoProperties.isTree())InterpollateProcess(Root, OutputFile, HistoProperties);
   }printf("\n");
}


void Draw2DHistogramSplitCanvas(JSONWrapper::Object& Root, TFile* File, NameAndType& HistoProperties){
   if(HistoProperties.isIndexPlot && cutIndex<0)return;

   std::string SaveName = "";


   std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
   std::vector<TObject*> ObjectToDelete;
   for(unsigned int i=0;i<Process.size();i++){
      if(Process[i]["isinvisible"].toBool())continue;
      string filtExt("");
      if(Process[i].isTag("mctruthmode") ) { char buf[255]; sprintf(buf,"_filt%d",(int)Process[i]["mctruthmode"].toInt()); filtExt += buf; }

      TCanvas* c1 = new TCanvas("c1","c1",500,500);
      c1->SetLogz(true);

      string dirName = Process[i]["tag"].c_str();while(dirName.find("/")!=std::string::npos)dirName.replace(dirName.find("/"),1,"-");
      TH1* hist = (TH1*)GetObjectFromPath(File,dirName + "/" + HistoProperties.name);
      if(!hist)continue;
      setStyle(Process[i], hist);

      SaveName = hist->GetName();
      ObjectToDelete.push_back(hist);
      hist->SetTitle("");
      hist->SetStats(kFALSE);

      hist->Draw("COLZ");
  
      TPaveText* leg = new TPaveText(0.20,0.95,0.40,0.80, "NDC");
      leg->SetFillColor(0);
      leg->SetFillStyle(0);  leg->SetLineColor(0);
      leg->SetTextAlign(12);
      leg->AddText(Process[i]["tag"].c_str());
      leg->Draw("same");
      ObjectToDelete.push_back(leg);

      DrawPreliminary();

      string SavePath = dropBadCharacters(SaveName + "_" + (Process[i])["tag"].toString());
      if(outDir.size()) SavePath = outDir +"/"+ SavePath;
      for(auto ext = plotExt.begin(); ext != plotExt.end(); ++ext)
      {
        system(string(("rm -f ") + SavePath + *ext).c_str());
        c1->SaveAs((SavePath + *ext).c_str());
      }
      delete c1;
   }

   for(unsigned int d=0;d<ObjectToDelete.size();d++){delete ObjectToDelete[d];}ObjectToDelete.clear();
}


void Draw2DHistogram(JSONWrapper::Object& Root, TFile* File, NameAndType& HistoProperties){
   if(HistoProperties.isIndexPlot && cutIndex<0)return;

   std::string SaveName = "";

   std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
   int NSampleToDraw = 0;
   for(unsigned int i=0;i<Process.size();i++){
      if(Process[i]["isinvisible"].toBool())continue;
      NSampleToDraw++;
   }
   int CanvasX = 3;
   int CanvasY = ceil(NSampleToDraw/CanvasX);
   TCanvas* c1 = new TCanvas("c1","c1",CanvasX*350,CanvasY*350);
   c1->Divide(CanvasX,CanvasY,0,0);


   std::vector<TObject*> ObjectToDelete;
   for(unsigned int i=0;i<Process.size();i++){
      if(Process[i]["isinvisible"].toBool())continue;

      TVirtualPad* pad = c1->cd(i+1);
      pad->SetLogz(true);
      pad->SetTopMargin(0.0); pad->SetBottomMargin(0.10);  pad->SetRightMargin(0.20);

      string dirName = Process[i]["tag"].c_str();while(dirName.find("/")!=std::string::npos)dirName.replace(dirName.find("/"),1,"-");
      TH1* hist = (TH1*)GetObjectFromPath(File,dirName + "/" + HistoProperties.name);
      if(!hist)continue;
      setStyle(Process[i], hist);
  
      SaveName = hist->GetName();
      ObjectToDelete.push_back(hist);
      hist->SetTitle("");
      hist->SetStats(kFALSE);

      hist->Draw("COLZ");
  
      TPaveText* leg = new TPaveText(0.10,0.995,0.30,0.90, "NDC");
      leg->SetFillColor(0);
      leg->SetFillStyle(0);  leg->SetLineColor(0);
      leg->SetTextAlign(12);
      leg->AddText(Process[i]["tag"].c_str());
      leg->Draw("same");
      ObjectToDelete.push_back(leg);
   }
   c1->cd(0);
   DrawPreliminary();

   string SavePath = dropBadCharacters(SaveName);
   if(outDir.size()) SavePath = outDir +"/"+ SavePath; 
   for(auto ext = plotExt.begin(); ext != plotExt.end(); ++ext)
   {
     system(string(("rm -f ") + SavePath + *ext).c_str());
     c1->SaveAs((SavePath + *ext).c_str());
   }
   for(unsigned int d=0;d<ObjectToDelete.size();d++){delete ObjectToDelete[d];}ObjectToDelete.clear();
   delete c1;
}

void Draw1DHistogram(JSONWrapper::Object& Root, TFile* File, NameAndType& HistoProperties){
   
   if(HistoProperties.isIndexPlot && cutIndex<0)return;

   TCanvas* c1 = new TCanvas("c1","c1",800,800);

   TPad* t1 = new TPad("t1","t1", 0.0, 0.2, 1.0, 1.0);
   t1->SetFillColor(0);
   t1->SetBorderMode(0);
   t1->SetBorderSize(2);
   t1->SetLogy();
   t1->SetTickx(1);
   t1->SetTicky(1);
   t1->SetLeftMargin(0.10);
   t1->SetRightMargin(0.05);
   t1->SetTopMargin(0.05);
   t1->SetBottomMargin(0.10);
   t1->SetFrameFillStyle(0);
   t1->SetFrameBorderMode(0);
   t1->SetFrameFillStyle(0);
   t1->SetFrameBorderMode(0);

   t1->Draw();
   t1->cd();
   if(!noLog) t1->SetLogy(true);
   
   //TLegend* legA  = new TLegend(0.845,0.2,0.99,0.99, "NDC"); 
   //   TLegend* legA  = new TLegend(0.51,0.93,0.67,0.75, "NDC"); 
   // TLegend* legB  = new TLegend(0.67,0.93,0.83,0.75, "NDC");
   TLegend *legA = new TLegend(0.30,0.74,0.93,0.96, "NDC");
   legA->SetHeader("");
   legA->SetNColumns(3);   
   legA->SetBorderSize(0);
   legA->SetTextFont(42);   legA->SetTextSize(0.03);
   legA->SetLineColor(0);   legA->SetLineStyle(1);   legA->SetLineWidth(1);
   legA->SetFillColor(0);   legA->SetFillStyle(0);//blind>-1E99?1001:0);

   THStack* stack = new THStack("MC","MC");
   TH1*     data  = NULL;
   TH1 *     mc   = NULL;
   TH1 *     mcPlusRelUnc = NULL;
   std::vector<TH1*> spimpose;
   std::vector<TString> spimposeOpts;
   std::vector<TObject*> ObjectToDelete;
   std::string SaveName = HistoProperties.name;
   std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
   double Maximum = -1E100;
   double Minimum =  1E100;
   double SignalMin = 5E-2;
   ObjectToDelete.push_back(stack);

   for(unsigned int i=0;i<Process.size();i++){
      if(Process[i]["isinvisible"].toBool())continue;
      string dirName = Process[i]["tag"].c_str();while(dirName.find("/")!=std::string::npos)dirName.replace(dirName.find("/"),1,"-");
      TH1* hist = (TH1*)GetObjectFromPath(File,dirName + "/" + HistoProperties.name);
      if(!hist)continue;
      setStyle(Process[i], hist);
      
      fixExtremities(hist,fixOverflow,fixUnderflow); 
      if(Process[i].isTag("normto")) hist->Scale( Process[i]["normto"].toDouble()/hist->Integral() );

      if(Maximum<hist->GetMaximum())Maximum=hist->GetMaximum();
      if(Minimum>hist->GetMinimum())Minimum=hist->GetMinimum();
      ObjectToDelete.push_back(hist);

      if(Process[i]["isdata"].toBool()){
          if(!data){legA->AddEntry(hist, Process[i]["tag"].c_str(), "P E");}
          if(!data){data = (TH1D*)hist->Clone("data");checkSumw2(data);}else{data->Add(hist);}
      }else if(Process[i].isTag("spimpose") && Process[i]["spimpose"].toBool()){
          legA->AddEntry(hist, Process[i]["tag"].c_str(), "L" );
          spimposeOpts.push_back( "hist" );
          spimpose.push_back(hist);
          if(SignalMin>hist->GetMaximum()*1E-2) SignalMin=hist->GetMaximum()*1E-2;
      }else{
	 stack->Add(hist, "HIST"); //Add to Stack
         legA->AddEntry(hist, Process[i]["tag"].c_str(), "F");	 
         if(!mc){mc = (TH1D*)hist->Clone("mc");checkSumw2(mc);}else{mc->Add(hist);}      
      }
   }

   if(mc   && Maximum<mc  ->GetMaximum()) Maximum=mc  ->GetMaximum()*1.1;
   if(data && Maximum<data->GetMaximum()) Maximum=data->GetMaximum()*1.1;
   Maximum = noLog?Maximum*1.5:Maximum*2E3;
   Minimum = std::min(SignalMin, std::max(5E-2, Minimum)) * 0.9;  //important when there are negative content
 
   if(stack->GetNhists()<=0){//add a dummy histogram to make sure that we have something to draw
      TH1* dummy = NULL;
      if(!dummy && data)dummy=(TH1*)data->Clone("dummy");
      if(!dummy && spimpose.size()>0)dummy=(TH1*)spimpose[0]->Clone("dummy");
      if(!dummy) printf("Error the frame is still empty\n");
      dummy->Reset();
      stack->Add(dummy);
      ObjectToDelete.push_back(dummy);
   }

   //draw the stack which will serve as frame
   stack->Draw("");
   stack->SetTitle("");
   stack->GetXaxis()->SetTitle(((TH1*)stack->GetStack()->At(0))->GetXaxis()->GetTitle());
   stack->GetYaxis()->SetTitle(((TH1*)stack->GetStack()->At(0))->GetYaxis()->GetTitle()); 
   stack->GetXaxis()->SetLabelOffset(0.007);
   stack->GetXaxis()->SetLabelSize(0.03);
   stack->GetXaxis()->SetTitleOffset(1.0);
   stack->GetXaxis()->SetTitleFont(42);
   stack->GetXaxis()->SetTitleSize(0.035);
   stack->GetYaxis()->SetLabelFont(42);
   stack->GetYaxis()->SetLabelOffset(0.007);
   stack->GetYaxis()->SetLabelSize(0.03);
   stack->GetYaxis()->SetTitleOffset(1.3);
   stack->GetYaxis()->SetTitleFont(42);
   stack->GetYaxis()->SetTitleSize(0.035);
   stack->GetYaxis()->SetRangeUser(Minimum, Maximum);
   stack->SetMinimum(Minimum);
   stack->SetMaximum(Maximum);
   t1->Update();

   if(showUnc && mc){
      mcPlusRelUnc = (TH1 *) mc->Clone("totalmcwithunc");checkSumw2(mcPlusRelUnc); mcPlusRelUnc->SetDirectory(0); 
      for(int ibin=1; ibin<=mcPlusRelUnc->GetXaxis()->GetNbins(); ibin++){
         Double_t error=sqrt(pow(mcPlusRelUnc->GetBinError(ibin),2)+pow(mcPlusRelUnc->GetBinContent(ibin)*baseRelUnc,2));
         mcPlusRelUnc->SetBinError(ibin,error);
      }
      mcPlusRelUnc->SetFillStyle(3427);
      mcPlusRelUnc->SetFillColor(kGray+1);
      mcPlusRelUnc->SetMarkerStyle(1);
   }

   if(data){
      if(blind>-1E99){
         if(data->GetBinLowEdge(data->FindBin(blind)) != blind){printf("Warning blinding changed from %f to %f for %s inorder to stay on bin edges\n", blind, data->GetBinLowEdge(data->FindBin(blind)),  HistoProperties.name.c_str() );}
         for(unsigned int i=data->FindBin(blind);i<=data->GetNbinsX()+1; i++){   
            data->SetBinContent(i, 0); data->SetBinError(i, 0);
         }
         TH1 *hist=(TH1*)stack->GetHistogram();
         TPave* blinding_box = new TPave(data->GetBinLowEdge(data->FindBin(blind)),  hist->GetMinimum(), data->GetXaxis()->GetXmax(), hist->GetMaximum(), 0, "NB" );  
         blinding_box->SetFillColor(15);         blinding_box->SetFillStyle(3013);         blinding_box->Draw("same F");
         legA->AddEntry(blinding_box, "blinded area" , "F");
         ObjectToDelete.push_back(blinding_box);
      }

      data->Draw("E1 same");
   }

   for(size_t ip=0; ip<spimpose.size(); ip++){
     TString opt=spimposeOpts[ip];
     spimpose[ip]->Draw(opt + "same");
   }

   //compare data and MC
   if(showChi2){
       TPaveText *pave = new TPaveText(0.6,0.85,0.8,0.9,"NDC");
       pave->SetBorderSize(0);
       pave->SetFillStyle(0);
       pave->SetTextAlign(32);
       pave->SetTextFont(42);
       char buf[100];
       if(data && mc && data->Integral()>0 && mc->Integral()>0){
	   sprintf(buf,"#chi^{2}/ndof : %3.2f", data->Chi2Test(mc,"WWCHI2/NDF") );
	   pave->AddText(buf);
	   sprintf(buf,"K-S prob: %3.2f", data->KolmogorovTest(mc));
	   pave->AddText(buf);

        }else if(mc && spimpose.size()>0 && mc->Integral()>0){
	   for(size_t ip=0; ip<spimpose.size(); ip++){
	       if(spimpose[ip]->Integral()<=0) continue;
	       sprintf(buf,"#chi^{2}/ndof : %3.2f, K-S prob: %3.2f", spimpose[ip]->Chi2Test(mc,"WWCHI2/NDF"), spimpose[ip]->KolmogorovTest(mc) );
	       pave->AddText(buf);
	   }
	}
       pave->Draw();
   }

   
   legA->Draw("same");


   std::vector<TH1 *> compDists;
   if(data)                   compDists.push_back(data);
   else if(spimpose.size()>0) compDists=spimpose;
   
   if( !(mc && compDists.size() && showRatioBox) ){
       t1->SetPad(0,0,1,1);
   }else{
       c1->cd();
       TPad *t2 = new TPad("t2", "t2",0.0,0.0, 1.0,0.2);
       t2->SetFillColor(0);
       t2->SetBorderMode(0);
       t2->SetBorderSize(2);
       t2->SetGridy();
       t2->SetTickx(1);
       t2->SetTicky(1);
       t2->SetLeftMargin(0.10);
       t2->SetRightMargin(0.05);
       t2->SetTopMargin(0.0);
       t2->SetBottomMargin(0.20);
       t2->SetFrameFillStyle(0);
       t2->SetFrameBorderMode(0);
       t2->SetFrameFillStyle(0);
       t2->SetFrameBorderMode(0);
       t2->Draw();
       t2->cd();
       t2->SetGridy(true);
       t2->SetPad(0,0.0,1.0,0.2);

       //mc stats
       TH1D *denRelUncH=0;
       if(mcPlusRelUnc) denRelUncH=(TH1D *) mcPlusRelUnc->Clone("mcrelunc");
       else             denRelUncH=(TH1D *) mc->Clone("mcrelunc");
       checkSumw2(denRelUncH);
       for(int xbin=1; xbin<=denRelUncH->GetXaxis()->GetNbins(); xbin++)
	 {
	   if(denRelUncH->GetBinContent(xbin)==0) continue;
	   Double_t err=denRelUncH->GetBinError(xbin)/denRelUncH->GetBinContent(xbin);
	   denRelUncH->SetBinContent(xbin,1);
	   denRelUncH->SetBinError(xbin,err);
	 }
       TGraphErrors *denRelUnc=new TGraphErrors(denRelUncH);
       denRelUnc->SetLineColor(1);
       denRelUnc->SetFillStyle(3001);
       denRelUnc->SetFillColor(kGray);
       denRelUnc->SetMarkerColor(1);
       denRelUnc->SetMarkerStyle(1);
       denRelUncH->Reset("ICE");       
       denRelUncH->SetTitle("");
       denRelUncH->SetStats(kFALSE);
       denRelUncH->Draw();
       denRelUnc->Draw("3");
       float yscale = (1.0-0.2)/(0.2);       
       denRelUncH->GetYaxis()->SetTitle("Data/#Sigma MC");
       denRelUncH->GetXaxis()->SetTitle(""); //drop the tile to gain space
       //denRelUncH->GetYaxis()->CenterTitle(true);
       denRelUncH->SetMinimum(0.4);
       denRelUncH->SetMaximum(1.6);
       //denRelUncH->SetMinimum(0);
       //denRelUncH->SetMaximum(data->GetBinContent(data->GetMaximumBin())*1.10);

       denRelUncH->GetXaxis()->SetLabelFont(42);
       denRelUncH->GetXaxis()->SetLabelOffset(0.007);
       denRelUncH->GetXaxis()->SetLabelSize(0.03 * yscale);
       denRelUncH->GetXaxis()->SetTitleFont(42);
       denRelUncH->GetXaxis()->SetTitleSize(0.035 * yscale);
       denRelUncH->GetXaxis()->SetTitleOffset(0.8);
       denRelUncH->GetYaxis()->SetLabelFont(42);
       denRelUncH->GetYaxis()->SetLabelOffset(0.007);
       denRelUncH->GetYaxis()->SetLabelSize(0.03 * yscale);
       denRelUncH->GetYaxis()->SetTitleFont(42);
       denRelUncH->GetYaxis()->SetTitleSize(0.035 * yscale);
       denRelUncH->GetYaxis()->SetTitleOffset(0.3);
      
       //add comparisons
       for(size_t icd=0; icd<compDists.size(); icd++){
	   TString name("CompHistogram"); name+=icd;
	   TH1D *dataToObsH = (TH1D*)compDists[icd]->Clone(name);
	   checkSumw2(dataToObsH);
	   dataToObsH->Divide(mc);
	   dataToObsH->Draw("same");
       }

      if(data && blind>-1E99){
         TPave* blinding_box = new TPave(data->GetBinLowEdge(data->FindBin(blind)), 0.4, data->GetXaxis()->GetXmax(), 1.6, 0, "NB" );
         blinding_box->SetFillColor(15);         blinding_box->SetFillStyle(3013);         blinding_box->Draw("same F");
         ObjectToDelete.push_back(blinding_box);
      }


   }

   t1->cd();
   c1->cd();
   DrawPreliminary(t1);
   c1->Modified();  
   c1->Update();

   string SavePath = dropBadCharacters(SaveName);
   if(outDir.size()) SavePath = outDir +"/"+ SavePath;
   for(auto i = plotExt.begin(); i != plotExt.end(); ++i){
     //system(string(("rm -f ") + SavePath + *i).c_str());
     c1->SaveAs((SavePath + *i).c_str());
   }
   
   delete c1;
   for(unsigned int d=0;d<ObjectToDelete.size();d++){delete ObjectToDelete[d];}ObjectToDelete.clear();
   delete legA;
}


void ConvertToTex(JSONWrapper::Object& Root, TFile* File, NameAndType& HistoProperties){
   if(HistoProperties.isIndexPlot && cutIndex<0)return;

   FILE* pFile = NULL;

   std::vector<TObject*> ObjectToDelete;
   TH1* stack = NULL; 
   std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
   for(unsigned int i=0;i<Process.size();i++){
      string dirName = Process[i]["tag"].c_str();while(dirName.find("/")!=std::string::npos)dirName.replace(dirName.find("/"),1,"-");
      TH1* hist = (TH1*)GetObjectFromPath(File,dirName + "/" + HistoProperties.name);
      if(!hist)continue;

      if(!pFile){
         string SavePath = dropBadCharacters(string(hist->GetName()) + ".tex");
         SavePath = outDir + string("/") + SavePath;
         system(string(("rm -f ") + SavePath).c_str());
         pFile = fopen(SavePath.c_str(), "w");

         fprintf(pFile, "\\begin{table}[htp]\n");
         fprintf(pFile, "\\begin{center}\n");
         fprintf(pFile, "\\caption{%2.2s}\n", hist->GetName());
         fprintf(pFile, "\\label{tab:table%10s}\n", hist->GetName());

         string colfmt = "|l";  for(int b=1;b<=hist->GetXaxis()->GetNbins();b++){colfmt = colfmt + "c";} colfmt+="|";
         string colname = "";   for(int b=1;b<=hist->GetXaxis()->GetNbins();b++){
            std::string tempcolname =  hist->GetXaxis()->GetBinLabel(b);
            if(tempcolname.find("_")!=std::string::npos || tempcolname.find("^")!=std::string::npos)tempcolname = string("$") + tempcolname + "$";
	    while(tempcolname.find("#")!=std::string::npos)tempcolname.replace(tempcolname.find("#"),1,"$\\");
	    if(tempcolname.find("$\\wedge")!=std::string::npos)tempcolname.replace(tempcolname.find("$\\wedge"),8,"\\wedge");
	    if(tempcolname.find("$\\geq")!=std::string::npos)tempcolname.replace(tempcolname.find("$\\geq"),5,"$\\geq$");
	    if(tempcolname.find("{miss}")!=std::string::npos)tempcolname.replace(tempcolname.find("{miss}"),6,"{miss}$");
	    if(tempcolname.find("40 GeV$")!=std::string::npos)tempcolname.replace(tempcolname.find("40 GeV$"),7,"40 GeV"); // This is too specific. Better to change the labels from the beginning 
	    if(tempcolname.find("$\\tau")!=std::string::npos)tempcolname.replace(tempcolname.find("$\\tau"),6,"\\tau"); // This is too specific. Better to change the labels from the beginning 
	    colname = colname + "& " + tempcolname;
         }
         fprintf(pFile, "\\begin{tabular}{ %s } \\hline\n", colfmt.c_str());
         fprintf(pFile, "Process %s \\\\ \\hline\\hline\n", colname.c_str());
      }
      ObjectToDelete.push_back(hist);

      std::string CleanTag = Process[i]["tag"].c_str();
      if(CleanTag.find("#")!=std::string::npos)CleanTag = string("$") + CleanTag + "$";
      while(CleanTag.find("#")!=std::string::npos)CleanTag.replace(CleanTag.find("#"),1,"\\");

      if((!Process[i].isTag("spimpose") || !Process[i]["spimpose"].toBool()) && !Process[i]["isdata"].toBool()){
         //Add to Stack
	if(!stack){stack = (TH1*)hist->Clone("Stack");checkSumw2(stack);}else{stack->Add(hist,1.0);}

         char numberastext[2048]; numberastext[0] = '\0';
         for(int b=1;b<=hist->GetXaxis()->GetNbins();b++){sprintf(numberastext,"%s & %s",numberastext, utils::toLatexRounded(hist->GetBinContent(b), hist->GetBinError(b),-1,doPowers).c_str());}
         fprintf(pFile, "%s %s \\\\\n",CleanTag.c_str(), numberastext);         
       }else{
          //Add to Canvas   
          if(stack){
            char numberastext[2048]; numberastext[0] = '\0';
            for(int b=1;b<=hist->GetXaxis()->GetNbins();b++){sprintf(numberastext,"%s & %s",numberastext, utils::toLatexRounded(stack->GetBinContent(b), stack->GetBinError(b),-1,doPowers).c_str());}
            fprintf(pFile, "%s %s \\\\\n\\hline\n","Total expected", numberastext);
             ObjectToDelete.push_back(stack);
             stack=NULL;
          }

          if(Process[i]["isdata"].toBool()){fprintf(pFile,"\\hline\n");}

          char numberastext[2048]; numberastext[0] = '\0';
          for(int b=1;b<=hist->GetXaxis()->GetNbins();b++){sprintf(numberastext,"%s & %s",numberastext, utils::toLatexRounded(hist->GetBinContent(b), hist->GetBinError(b),-1,doPowers).c_str());}
          fprintf(pFile, "%s %s \\\\\n",CleanTag.c_str(), numberastext);
       }

   }

   if(pFile){
      fprintf(pFile,"\\hline\n");
      fprintf(pFile,"\\end{tabular}\n");
      fprintf(pFile,"\\end{center}\n");
      fprintf(pFile,"\\end{table}\n");
      fclose(pFile);
   }
   for(unsigned int d=0;d<ObjectToDelete.size();d++){delete ObjectToDelete[d];}ObjectToDelete.clear();
}






int main(int argc, char* argv[]){
   gROOT->LoadMacro("../../src/tdrstyle.C");
   setTDRStyle();  
   gStyle->SetPadTopMargin   (0.06);
   gStyle->SetPadBottomMargin(0.12);
   gStyle->SetPadRightMargin (0.16);
   gStyle->SetPadLeftMargin  (0.14);
   gStyle->SetTitleSize(0.04, "XYZ");
   gStyle->SetTitleXOffset(1.1);
   gStyle->SetTitleYOffset(1.45);
   gStyle->SetPalette(1);
   gStyle->SetNdivisions(505);

//   gErrorIgnoreLevel = kWarning; //supress info message

   std::vector<string> histoNameMask;

   for(int i=1;i<argc;i++){
     string arg(argv[i]);
     //printf("--- %i - %s\n",i,argv[i]);

     if(arg.find("--help")!=string::npos){
        printf("--help   --> print this helping text\n");

        printf("--iLumi   --> integrated luminosity to be used for the MC rescale\n");
        printf("--iEcm    --> center of mass energy (TeV) = 8 TeV by default\n");
        printf("--isSim   --> print CMS Simulation instead of the standard title\n");
        printf("--inDir   --> path to the directory containing the .root files to process\n");
        printf("--outDir  --> path of the directory that will contains the output plots and tables\n");
        printf("--outFile --> path of the output summary .root file\n");
        printf("--json    --> containing list of process (and associated style) to process to process\n");
        printf("--only    --> processing only the objects matching this regex expression\n");
        printf("--index   --> will do the projection on that index for histos of type cutIndex\n");
        printf("--chi2    --> show the data/MC chi^2\n"); 
        printf("--showUnc --> show stat uncertainty (if number is given use it as relative bin by bin uncertainty (e.g. lumi)\n"); 
	printf("--noLog   --> use linear scale\n");
        printf("--no1D   --> Skip processing of 1D objects\n");
        printf("--no2D   --> Skip processing of 2D objects\n");
        printf("--noTree --> Skip processing of Tree objects\n");
        printf("--noTex  --> Do not create latex table (when possible)\n");
        printf("--noPowers --> Do not use powers of 10 for numbers in tables\n");
        printf("--noRoot --> Do not make a summary .root file\n");
        printf("--noPlot --> Do not creates plot files (useful to speedup processing)\n");
	printf("--plotExt --> extension to save, you can specify multiple extensions by repeating this option\n");
	printf("--cutflow --> name of the histogram with the original number of events (cutflow by default)\n");
        printf("--splitCanvas --> (only for 2D plots) save all the samples in separated pltos\n");
        printf("--removeRatioPlot --> if you want to remove ratio plots between Data ad Mc\n");
        printf("--removeUnderFlow --> Remove the Underflow bin in the final plots\n");
        printf("--removeOverFlow --> Remove the Overflow bin in the final plots\n");

        printf("command line example: runPlotter --json ../data/beauty-samples.json --iLumi 2007 --inDir OUT/ --outDir OUT/plots/ --outFile plotter.root --noRoot --noPlot\n");
	return 0;
     }

     if(arg.find("--iLumi"  )!=string::npos && i+1<argc){ sscanf(argv[i+1],"%lf",&iLumi); i++; printf("Lumi = %f\n", iLumi); }
     if(arg.find("--iEcm"   )!=string::npos && i+1<argc){ sscanf(argv[i+1],"%lf",&iEcm); i++; printf("Ecm = %f TeV\n", iEcm); }
     if(arg.find("--blind"  )!=string::npos && i+1<argc){ sscanf(argv[i+1],"%lf",&blind); i++; printf("Blind above = %f\n", blind); }
     if(arg.find("--inDir"  )!=string::npos && i+1<argc){ inDir    = argv[i+1];  i++;  printf("inDir = %s\n", inDir.c_str());  }
     if(arg.find("--outDir" )!=string::npos && i+1<argc){ outDir   = argv[i+1];  i++;  printf("outDir = %s\n", outDir.c_str());  }
     if(arg.find("--outFile")!=string::npos && i+1<argc){ outFile  = argv[i+1];  i++; printf("output file = %s\n", outFile.c_str()); }
     if(arg.find("--json"   )!=string::npos && i+1<argc){ jsonFile = argv[i+1];  i++;  }
     if(arg.find("--only"   )!=string::npos && i+1<argc){ histoNameMask.push_back(argv[i+1]); printf("Only histograms matching (regex) ewpression '%s' are processed\n", argv[i+1]); i++;  }
     if(arg.find("--index"  )!=string::npos && i+1<argc){ sscanf(argv[i+1],"%d",&cutIndex); i++; onlyCutIndex=(cutIndex>=0); printf("index = %i\n", cutIndex);  }
     if(arg.find("--chi2"  )!=string::npos)             { showChi2 = true;  }
     if(arg.find("--showUnc") != string::npos) { 
       showUnc=true; 
       if(i+1<argc) { 
	 string nextArg(argv[i+1]); 
	 if(nextArg.find("--")==string::npos)
	   {
	     sscanf(argv[i+1],"%lf",&baseRelUnc);
	     i++;
	   }
       }
       printf("Uncertainty band will be included for MC with base relative uncertainty of: %3.2f",baseRelUnc);
     }
     if(arg.find("--isSim")!=string::npos){ isSim = true;    }
     if(arg.find("--noLog")!=string::npos){ noLog = true;    }
     if(arg.find("--noTree"  )!=string::npos){ doTree = false;    }
     if(arg.find("--no2D"  )!=string::npos){ do2D = false;    }
     if(arg.find("--no1D"  )!=string::npos){ do1D = false;    }
     if(arg.find("--noTex" )!=string::npos){ doTex= false;    }
     if(arg.find("--noPowers" )!=string::npos){ doPowers= false;    }
     if(arg.find("--noRoot")!=string::npos){ StoreInFile = false;    }
     if(arg.find("--noPlot")!=string::npos){ doPlot = false;    }
     if(arg.find("--removeRatioPlot")!=string::npos){ showRatioBox = false; printf("No ratio plot between Data and Mc \n"); }
     if(arg.find("--removeUnderFlow")!=string::npos){ fixUnderflow = false; printf("No UnderFlowBin \n"); }
     if(arg.find("--removeOverFlow")!=string::npos){ fixOverflow = false; printf("No OverFlowBin \n"); }
     if(arg.find("--plotExt" )!=string::npos && i+1<argc){ plotExt.push_back(argv[i+1]);  i++;  printf("saving plots as = %s\n", argv[i]);  }
     if(arg.find("--cutflow" )!=string::npos && i+1<argc){ cutflowhisto   = argv[i+1];  i++;  printf("Normalizing from 1st bin in = %s\n", cutflowhisto.c_str());  }
     if(arg.find("--splitCanvas")!=string::npos){ splitCanvas = true;    }
     if(arg.find("--fileOption" )!=string::npos && i+1<argc){ fileOption = argv[i+1];  i++;  printf("FileOption = %s\n", fileOption.c_str());  }
   } 
   system( (string("mkdir -p ") + outDir).c_str());
   if(plotExt.size() == 0)
     plotExt.push_back(".png");

   char buf[255];
   sprintf(buf, "_Index%d", cutIndex);
   cutIndexStr = buf;

   TFile* OutputFile = new TFile(outFile.c_str(),fileOption.c_str());

   JSONWrapper::Object Root(jsonFile, true);
   std::list<NameAndType> histlist;
//   if(fileOption!="READ"){
      GetListOfObject(Root,inDir,histlist);
//   }else{
//      GetListOfObject(Root,inDir,histlist, OutputFile, "");
//   }
   histlist.sort();
   histlist.unique();   

   printf("Progressing Bar              :0%%       20%%       40%%       60%%       80%%       100%%\n");
   string csvFile(outDir +"/histlist.csv");
   system(("echo \"\" > " + csvFile).c_str());

   std::list<NameAndType>::iterator it= histlist.begin();
   while(it!= histlist.end()){
       bool passMasking = (histoNameMask.size()==0);  
       for(unsigned int i=0;i<histoNameMask.size();i++){if(std::regex_match(it->name,std::regex(histoNameMask[i])))passMasking=true;}
       if(!passMasking){it=histlist.erase(it); continue; }
       if(!do2D   &&(it->is2D() || it->is3D())){it=histlist.erase(it); continue;}
       if(!do1D   && it->is1D()){it=histlist.erase(it); continue;}
       if(!doTree && it->isTree()){it=histlist.erase(it); continue;}
       it++;
   }

   if(fileOption!="READ")SavingToFile(Root,inDir,OutputFile, histlist);       

   int ictr =0;
   int TreeStep = std::max(1,(int)(histlist.size()/50));
   printf("Plotting                     :");
   for(std::list<NameAndType>::iterator it= histlist.begin(); it!= histlist.end(); it++,ictr++){
       if(ictr%TreeStep==0){printf(".");fflush(stdout);}
   
       if(doTex && (it->name.find("eventflow")!=std::string::npos || it->name.find("evtflow")!=std::string::npos) && it->name.find("optim_eventflow")==std::string::npos){    ConvertToTex(Root,OutputFile,*it); }
       if(doPlot && do2D  && it->is2D()){                      if(!splitCanvas){Draw2DHistogram(Root,OutputFile,*it); }else{Draw2DHistogramSplitCanvas(Root,OutputFile,*it);}}
       if(doPlot && do1D  && it->is1D()){ Draw1DHistogram(Root,OutputFile,*it); }
        system(("echo \"" + it->name + "\" >> " + csvFile).c_str());
   }printf("\n");
   OutputFile->Close();
   
   //system(("python ${CMSSW_BASE}/src/CMGTools/HtoZZ2l2nu/data/html/generateJSONplotterFromList.py -i " + csvFile + " -o "+outDir+"/plotter.json").c_str());
   system(("rm " + csvFile).c_str());
   //system(("cp ${CMSSW_BASE}/src/CMGTools/HtoZZ2l2nu/data/html/index.html " + outDir).c_str());
   //printf("You can browse the results using %s/index.html\n",outDir.c_str());
}

