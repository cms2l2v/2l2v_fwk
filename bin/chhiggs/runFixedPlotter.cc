#include <string>
#include <vector>
#include <iostream>
#include <list>
#include <iterator>
#include <algorithm>
#include <unordered_map>

#include "TSystem.h"
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
#include "TF1.h"
#include "TCutG.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TPaveText.h"
#include "THStack.h"
#include "TRandom.h"

#include "UserCode/llvv_fwk/interface/tdrstyle.h"
#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/JSONWrapper.h"

using namespace std;

bool debug=false;

int cutIndex=-1;
string cutIndexStr="";
double iLumi = 2007;
double iEcm=8;
bool showChi2 = false;
bool showUnc=false;
double baseRelUnc=0.027;
bool noLog=false;
bool logX=false;
bool isSim=false;
bool do2D  = true;
bool do1D  = true;
bool doTex = true;
bool doPowers = true;
bool StoreInFile = true;
bool doPlot = true;
bool splitCanvas = false;
bool onlyCutIndex = false;
string inDir   = "OUTNew/";
string jsonFile = "../../data/beauty-samples.json";
string outDir  = "Img/";
string plotExt = ".png";
string outFile = "plotter.root";
string cutflowhisto = "all_cutflow";
bool forceMerge=false;
bool useMerged=false;
bool jodorStyle=false;
bool generatePseudoData=false;
TH1* myPseudoData=NULL;
TString myPseudoDataName="";

struct stSampleInfo{ double PURescale_up; double PURescale_down; double initialNumberOfEvents;};
std::unordered_map<string, stSampleInfo> sampleInfoMap;
std::unordered_map<string, bool> FileExist;

struct NameAndType{
   std::string name;
   bool type; 
   bool isIndexPlot;
   NameAndType(std::string name_,  bool type_, bool isIndexPlot_){name = name_; type = type_; isIndexPlot = isIndexPlot_;}

   bool operator==(const NameAndType& a){ return a.name == name;}
   bool operator< (const NameAndType& a){ return a.name < name;}
 };

void checkSumw2(TH1 *h) { if(h==0) return;  if(h->GetDefaultSumw2()) h->Sumw2();  }
void MergeSplittedSamples(JSONWrapper::Object& Root, std::string RootDir);


TObject* GetObjectFromPath(TDirectory* File, std::string Path, bool GetACopy=false)
{
   size_t pos = Path.find("/");
   if(pos < 256){
      std::string firstPart = Path.substr(0,pos);
      std::string endPart   = Path.substr(pos+1,Path.length());
      TDirectory* TMP = (TDirectory*)File->Get(firstPart.c_str());
//      if(TMP!=NULL)return GetObjectFromPath(TMP,endPart,GetACopy);
      if(TMP!=NULL){
         TObject* TMP2 =  GetObjectFromPath(TMP,endPart,GetACopy);
//         if(!TMP2)printf("BUG: %s\n",Path.c_str());
         return TMP2;
      }

//      printf("BUG: %s\n",Path.c_str());
      return NULL;
   }else{
      TObject* TMP = File->Get(Path.c_str());
//      if(!TMP)printf("BUG: %s\n",Path.c_str());

      if(GetACopy){
	return TMP->Clone();
      }else{
         return TMP;
      }
   }

}

void GetListOfObject(JSONWrapper::Object& Root, std::string RootDir, std::list<NameAndType>& histlist, TDirectory* dir=NULL, std::string parentPath=""){

  if(parentPath=="" && !dir)
    {
      int dataProcessed = 0;
      int signProcessed = 0;
      int bckgProcessed = 0; 

      std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
      //loop on all Proc
      for(size_t ip=0; ip<Process.size(); ip++){
          if(Process[ip]["isinvisible"].toBool())continue;
	  bool isData (  Process[ip]["isdata"].toBool()  );
          bool isSign ( !isData &&  ((Process[ip].isTag("spimpose") && Process[ip]["spimpose"].toBool()) || (Process[ip].isTag("issignal") && Process[ip]["issignal"].toBool()) ) );
  	  bool isMC   = !isData && !isSign; 
	  string filtExt("");
	  // TEST // if(Process[ip].isTag("mctruthmode") ) { char buf[255]; sprintf(buf,"_filt%d",(int)Process[ip]["mctruthmode"].toInt()); filtExt += buf; }
          //just to make it faster, only consider the first 3 sample of a same kind
          if(isData){if(dataProcessed>=1){continue;}else{dataProcessed++;}}
          if(isSign){if(signProcessed>=2){continue;}else{signProcessed++;}}
          if(isMC  ){if(bckgProcessed>3) {continue;}else{bckgProcessed++;}}

	  std::vector<JSONWrapper::Object> Samples = (Process[ip])["data"].daughters();
          //to make it faster only consider the first samples 
	  //for(size_t id=0; id<Samples.size()&&id<2; id++){
	  for(size_t id=0; id<Samples.size(); id++){
	      int split = 1;
	      if(!useMerged && Samples[id].isTag("split"))split = Samples[id]["split"].toInt();

	      string segmentExt; if(split>1) { char buf[255]; sprintf(buf,"_%i",0); segmentExt += buf; }
              // TEST // string secondFiltExt=filtExt;
              string secondFiltExt("");
              // TEST // if(TString((Samples[id])["dtag"].toString()).Contains("filt5")) secondFiltExt="_filt5";
              // TEST // if(TString((Samples[id])["dtag"].toString()).Contains("filt6")) secondFiltExt="_filt6";
              string FileName = RootDir + (Samples[id])["dtag"].toString() +  (Samples[id].isTag("suffix")?(Samples[id])["suffix"].toString():string("")) + filtExt + segmentExt + secondFiltExt + ".root";
              printf("Adding all objects from %25s to the list of considered objects\n",  FileName.c_str());
	      TFile* file = new TFile(FileName.c_str());
	      if(file->IsZombie())
		{
		  if(isData) dataProcessed--;
		  if(isSign) signProcessed--;
		  if(isMC) bckgProcessed--;
		}
	      GetListOfObject(Root,RootDir,histlist,(TDirectory*)file,"" );
	      file->Close();
	    }
	}

    //for(std::list<NameAndType>::iterator it= histlist.begin(); it!= histlist.end(); it++){printf("%s\n",it->name.c_str()); }

      return;
   }

   if(dir==NULL)return;

   TList* list = dir->GetListOfKeys();
   for(int i=0;i<list->GetSize();i++){
      TObject* tmp = GetObjectFromPath(dir,list->At(i)->GetName(),false);
      if(tmp->InheritsFrom("TTree")) continue;
//      if(tmp->InheritsFrom("TH1")){
//         histlist.push_back(NameAndType(parentPath+list->At(i)->GetName(), !(tmp->InheritsFrom("TH2") || tmp->InheritsFrom("TH3")) ) );
//      }else if(tmp->InheritsFrom("TDirectory")){
//         GetListOfObject(Root,RootDir,histlist,(TDirectory*)tmp,parentPath+ list->At(i)->GetName()+"/" );
//      }


      if(tmp->InheritsFrom("TDirectory")){
         GetListOfObject(Root,RootDir,histlist,(TDirectory*)tmp,parentPath+ list->At(i)->GetName()+"/" );
      }else if(tmp->InheritsFrom("TH1")){
        bool isTH1 = !(tmp->InheritsFrom("TH2") || tmp->InheritsFrom("TH3"));
        bool hasIndex = string(((TH1*)tmp)->GetXaxis()->GetTitle()).find("cut index")<string::npos;
        if(hasIndex){isTH1=true;}
	histlist.push_back(NameAndType(parentPath+list->At(i)->GetName(), isTH1, hasIndex ) );
      }

      delete tmp;
   }
   

}

void MergeSplittedSamples(JSONWrapper::Object& Root, std::string RootDir)
{
  std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
  for(unsigned int i=0;i<Process.size();i++){
    string filtExt("");
    // TEST // if(Process[i].isTag("mctruthmode") ) { char buf[255]; sprintf(buf,"_filt%d",(int)Process[i]["mctruthmode"].toInt()); filtExt += buf; }
    std::vector<JSONWrapper::Object> Samples = (Process[i])["data"].daughters();
    for(unsigned int j=0;j<Samples.size();j++){
      int split = 1;
      if(Samples[j].isTag("split"))split = Samples[j]["split"].toInt();
      if(split==1) continue;
      // TEST // string secondFiltExt=filtExt;
      string secondFiltExt("");
      // TEST // if(TString((Samples[j])["dtag"].toString()).Contains("filt5")) secondFiltExt="_filt5";
      // TEST // if(TString((Samples[j])["dtag"].toString()).Contains("filt6")) secondFiltExt="_filt6";

      TString toHadd("");
      for(int s=0;s<split;s++){
	string segmentExt; if(split>1) { char buf[255]; sprintf(buf,"_%i",s); segmentExt += buf; }
	string FileName = RootDir + (Samples[j])["dtag"].toString() + ((Samples[j].isTag("suffix"))?(Samples[j])["suffix"].toString():string("")) + filtExt + segmentExt + secondFiltExt + ".root";
	toHadd += " " + FileName;
      }
      TString outName(RootDir + (Samples[j])["dtag"].toString() + ((Samples[j].isTag("suffix"))?(Samples[j])["suffix"].toString():string("")) + filtExt + filtExt + ".root");
      //cout << "hadd -f -k "+ outName + toHadd << endl;
      gSystem->Exec("hadd -f -k "+ outName + toHadd);
    }
  }
}

void GetInitialNumberOfEvents(JSONWrapper::Object& Root, std::string RootDir, NameAndType HistoProperties){
   std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
   for(unsigned int i=0;i<Process.size();i++){
      if(!StoreInFile && Process[i]["isinvisible"].toBool())continue;
      string filtExt("");
      // TEST // if(Process[i].isTag("mctruthmode") ) { char buf[255]; sprintf(buf,"_filt%d",(int)Process[i]["mctruthmode"].toInt()); filtExt += buf; }
      std::vector<JSONWrapper::Object> Samples = (Process[i])["data"].daughters();
      for(unsigned int j=0;j<Samples.size();j++){
	 int split = 1;
	 if(!useMerged && Samples[j].isTag("split"))split = Samples[j]["split"].toInt();

         TH1* tmphist = NULL;
         // TEST // string secondFiltExt=filtExt;
         string secondFiltExt("");
         // TEST // if(TString((Samples[j])["dtag"].toString()).Contains("filt5")) secondFiltExt="_filt5";
         // TEST // if(TString((Samples[j])["dtag"].toString()).Contains("filt6")) secondFiltExt="_filt6";
         for(int s=0;s<split;s++){
	   string segmentExt; if(split>1) { char buf[255]; sprintf(buf,"_%i",s); segmentExt += buf; }
            string FileName = RootDir + (Samples[j])["dtag"].toString() + ((Samples[j].isTag("suffix"))?(Samples[j])["suffix"].toString():string("")) + filtExt + segmentExt + secondFiltExt + ".root";
            TFile* File = TFile::Open(FileName.c_str());
            bool& fileExist = FileExist[FileName];
	    if(!File || File->IsZombie() || !File->IsOpen() || File->TestBit(TFile::kRecovered) ){fileExist=false;  continue; }else{fileExist=true;}

            TH1* tmptmphist = (TH1*) GetObjectFromPath(File,HistoProperties.name); 
	    if(tmptmphist)
	      {
		if(!tmphist){gROOT->cd(); tmphist = (TH1*)tmptmphist->Clone(tmptmphist->GetName()); checkSumw2(tmphist); }else{tmphist->Add(tmptmphist);}
		delete tmptmphist;
	      }
            delete File;
         }
         if(!tmphist)continue;
         
	 bool isMC( !Process[i]["isdata"].toBool() && !Process[i]["isdatadriven"].toBool() );

         stSampleInfo& sampleInfo = sampleInfoMap[(Samples[j])["dtag"].toString()];

         double PUCentralnnorm =  1; if(tmphist->GetBinContent(3)>0)PUCentralnnorm = tmphist->GetBinContent(2) / tmphist->GetBinContent(3);
         double PUDownnorm     =  1; if(tmphist->GetBinContent(4)>0)PUDownnorm     = tmphist->GetBinContent(3) / tmphist->GetBinContent(4);
         double PUUpnorm       =  1; if(tmphist->GetBinContent(5)>0)PUUpnorm       = tmphist->GetBinContent(3) / tmphist->GetBinContent(5);
         sampleInfo.PURescale_down = PUDownnorm;
         sampleInfo.PURescale_up   = PUUpnorm;
         if(isMC)printf("PU Renormalization %25s Shift Down --> %6.2f  Central = %6.2f  Up Down --> %6.2f\n",(Samples[j])["dtag"].toString().c_str(),PUDownnorm, PUCentralnnorm, PUUpnorm);	


         double cnorm = 1.0;
         if(tmphist)cnorm = tmphist->GetBinContent(1);
         if(cnorm<=0 || !isMC)cnorm = 1.0;         
         if(cnorm==1 && isMC)printf("is there a problem with %s ? cnorm = %f\n",(Samples[j])["dtag"].toString().c_str(), cnorm);          
         if(!isMC)PUCentralnnorm = 1;
	 
	 double VBFMCRescale = tmphist->GetXaxis()->GetNbins()>5 ? tmphist->GetBinContent(6) / tmphist->GetBinContent(2) : 1.0;
	 if(VBFMCRescale!=0)          cnorm *= VBFMCRescale;
	 //printf("VBFMCRescale for sample %s is %f\n", (Samples[j])["dtag"].toString().c_str(), VBFMCRescale );
         sampleInfo.initialNumberOfEvents = cnorm; /// FIXME TO BE REACTIVATED ONCE FIXED THE initNorm HISTOGRAM FILLING / PUCentralnnorm;
         if(debug) cout << "cnorm: " << cnorm << ", PUCentralnnorm: " << PUCentralnnorm << endl;
         delete tmphist;
      }   
   }
}



void SavingPseudoDataToFile(JSONWrapper::Object& Root, std::string RootDir, NameAndType HistoProperties, TFile* OutputFile){
  
  TH1* hist = (TH1*) myPseudoData->Clone(myPseudoDataName);
  string dirName = "data";
  OutputFile->cd();
  TDirectory* subdir = OutputFile->GetDirectory(dirName.c_str());
  if(!subdir || subdir==OutputFile) subdir = OutputFile->mkdir(dirName.c_str());
  subdir->cd();
  hist->Write();
  delete hist;
}


void SavingToFile(JSONWrapper::Object& Root, std::string RootDir, NameAndType HistoProperties, TFile* OutputFile){
   std::vector<TObject*> ObjectToDelete;

   std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
   for(unsigned int i=0;i<Process.size();i++){
      TH1* hist = NULL;
      std::vector<JSONWrapper::Object> Samples = (Process[i])["data"].daughters();
      for(unsigned int j=0;j<Samples.size();j++){
         double Weight = 1.0;
         if(!Process[i]["isdata"].toBool() && !Process[i]["isdatadriven"].toBool())Weight*= iLumi;
         if(Samples[j].isTag("xsec")     )Weight*= Samples[j]["xsec"].toDouble();
	 string filtExt("");
	 // TEST // if(Process[i].isTag("mctruthmode") ) { char buf[255]; sprintf(buf,"_filt%d",(int)Process[i]["mctruthmode"].toInt()); filtExt += buf; }	 

	 std::vector<JSONWrapper::Object> BR = Samples[j]["br"].daughters();for(unsigned int b=0;b<BR.size();b++){Weight*=BR[b].toDouble();}
         stSampleInfo& sampleInfo = sampleInfoMap[(Samples[j])["dtag"].toString()];
         Weight /= sampleInfo.initialNumberOfEvents;

         if(HistoProperties.name.find("puup"  )!=string::npos){Weight *= sampleInfo.PURescale_up;}
         if(HistoProperties.name.find("pudown")!=string::npos){Weight *= sampleInfo.PURescale_down;}

         if(HistoProperties.name.find("optim_cut")!=string::npos){Weight=1.0;}

         int split = 1;
	 if(!useMerged && Samples[j].isTag("split"))split = Samples[j]["split"].toInt();
         // TEST // string secondFiltExt=filtExt;
         string secondFiltExt("");
         // TEST // if(TString((Samples[j])["dtag"].toString()).Contains("filt5")) secondFiltExt="_filt5";
         // TEST // if(TString((Samples[j])["dtag"].toString()).Contains("filt6")) secondFiltExt="_filt6";

         TH1* tmphist = NULL;
         for(int s=0;s<split;s++){
	   string segmentExt; if(split>1){ char buf[255]; sprintf(buf,"_%i",s); segmentExt += buf;}
	   
            string FileName = RootDir + (Samples[j])["dtag"].toString() + ((Samples[j].isTag("suffix"))?(Samples[j])["suffix"].toString():string("")) + filtExt + segmentExt + secondFiltExt + ".root";
            if(!FileExist[FileName])continue;
            TFile* File = new TFile(FileName.c_str());
            if(!File || File->IsZombie() || !File->IsOpen() || File->TestBit(TFile::kRecovered) )continue;
            TH1* tmptmphist = (TH1*) GetObjectFromPath(File,HistoProperties.name); 
            if(!tmptmphist){delete File;continue;}
            if(!tmphist){gROOT->cd(); tmphist = (TH1*)tmptmphist->Clone(tmptmphist->GetName()); checkSumw2(tmphist);}else{tmphist->Add(tmptmphist);}
            delete tmptmphist;
            delete File;
         }
         if(!tmphist)continue;
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

void GeneratePseudoData(TH1* pseudoData, TH1* mc){
  //pseudoData = (TH1*) mc->Clone("pseudoData");
  
  pseudoData->Reset("ICES"); // Reset Integral, Contents , Errors and Statistics (maintains Maximum ("M") and so on)

  // Add pseudo-data (poisson smearing of the total SM background)
  for(int xbin=1; xbin<=mc->GetXaxis()->GetNbins(); xbin++){
    pseudoData->SetBinContent(xbin, gRandom->Poisson(mc->GetBinContent(xbin)));
    pseudoData->SetBinError(xbin, sqrt(pseudoData->GetBinContent(xbin)));
  }
  //for(int i=1; i<=gRandom->Poisson(nexp); i++) hpseudo->Fill(h->GetRandom()); ?
  //Pietro Vischia
  //  Hm this fills a random number of times the histogram with a random x value, no? I think I need something like filling each bin with the bin content of the base histo smeared randomly along a poissonian using setbincontwnt
  //But I am not sure
  //Pedro Silva
  //  for(int xbin=1; xbin<=h->GetBins(); xbin++) for(int i=1; i<=gRandom->Poisson(h->GetBinContent(xbin)); i++) hpseudo->Fill( h->GetXaxis()->GetBinCenter(xbin));
  //can compare both
  //and see how different they are
  ///  or simply for(int xbin=1; xbin<=h->GetBins(); xbin++) hpseudo->SetBinContent(xbin, gRandom->Poisson(h->GetBinContent(xbin)) );
  //as you say
}




void Draw2DHistogramSplitCanvas(JSONWrapper::Object& Root, std::string RootDir, NameAndType HistoProperties){
   if(HistoProperties.isIndexPlot && cutIndex<0)return;

   std::string SaveName = "";

   TPaveText* T = new TPaveText(0.10,0.995,0.85,0.945, "NDC");
   T->SetFillColor(0);
   T->SetFillStyle(0);  T->SetLineColor(0);
   T->SetTextAlign(12);
   char Buffer[1024]; 
   if(isSim) sprintf(Buffer, "CMS simulation, #sqrt{s}=%.0f TeV, #scale[0.5]{#int} L=%.1f fb^{-1}", iEcm, iLumi/1000);
   else      sprintf(Buffer, "CMS preliminary, #sqrt{s}=%.0f TeV, #scale[0.5]{#int} L=%.1f fb^{-1}", iEcm, iLumi/1000);
   T->AddText(Buffer);

   std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
   std::vector<TObject*> ObjectToDelete;
   for(unsigned int i=0;i<Process.size();i++){
      if(Process[i]["isinvisible"].toBool())continue;
      string filtExt("");
      // TEST // if(Process[i].isTag("mctruthmode") ) { char buf[255]; sprintf(buf,"_filt%d",(int)Process[i]["mctruthmode"].toInt()); filtExt += buf; }

      TCanvas* c1 = new TCanvas("c1","c1",500,500);
      c1->SetLogz(true);

      TH1* hist = NULL;
      std::vector<JSONWrapper::Object> Samples = (Process[i])["data"].daughters();
      for(unsigned int j=0;j<Samples.size();j++){
         double Weight = 1.0;
         if(!Process[i]["isdata"].toBool()  && !Process[i]["isdatadriven"].toBool())Weight*= iLumi;
         if(Samples[j].isTag("xsec")     )Weight*= Samples[j]["xsec"].toDouble();
         std::vector<JSONWrapper::Object> BR = Samples[j]["br"].daughters();for(unsigned int b=0;b<BR.size();b++){Weight*=BR[b].toDouble();}
         stSampleInfo& sampleInfo = sampleInfoMap[(Samples[j])["dtag"].toString()];
         Weight /= sampleInfo.initialNumberOfEvents;

         if(HistoProperties.name.find("puup"  )!=string::npos){Weight *= sampleInfo.PURescale_up;}
         if(HistoProperties.name.find("pudown")!=string::npos){Weight *= sampleInfo.PURescale_down;}

         int split = 1;
	 if(!useMerged && Samples[j].isTag("split"))split = Samples[j]["split"].toInt();

         TH1* tmphist = NULL;
         for(int s=0;s<split;s++){
	   string segmentExt; if(split>1) { char buf[255]; sprintf(buf,"_%i",s); segmentExt += buf; }
	    
            string FileName = RootDir + (Samples[j])["dtag"].toString() + ((Samples[j].isTag("suffix"))?(Samples[j])["suffix"].toString():string("")) + filtExt + segmentExt + filtExt + ".root";
            if(!FileExist[FileName])continue;
            TFile* File = new TFile(FileName.c_str());
            if(!File || File->IsZombie() || !File->IsOpen() || File->TestBit(TFile::kRecovered) )continue;
            TH1* tmptmphist = (TH1*) GetObjectFromPath(File,HistoProperties.name); 
            if(!tmptmphist){delete File;continue;}
            if(!tmphist){gROOT->cd(); tmphist = (TH1*)tmptmphist->Clone(tmptmphist->GetName());checkSumw2(tmphist);}else{tmphist->Add(tmptmphist);}
            //if(Process[i]["isdata"].toBool())printf("%s --> %f*%f(%f)\n",(Samples[j])["dtag"].toString().c_str(), tmptmphist->Integral(),Weight, initialNumberOfEvents[(Samples[j])["dtag"].toString()]);
            delete tmptmphist;
            delete File;
         }
         if(!tmphist)continue;
         if(!hist){gROOT->cd(); hist = (TH1*)tmphist->Clone(tmphist->GetName());checkSumw2(hist);hist->Scale(Weight);}else{hist->Add(tmphist,Weight);}
         delete tmphist;
      }   
      if(!hist)continue;

      SaveName = hist->GetName();
      ObjectToDelete.push_back(hist);
      hist->SetTitle("");
      hist->SetStats(kFALSE);

      hist->Draw("COLZ");
  
      TPaveText* leg = new TPaveText(0.20,0.95,0.40,0.80, "NDC");
      leg->SetFillColor(0);
      leg->SetFillStyle(0);  leg->SetLineColor(0);
      leg->SetTextAlign(12);
      if(Process[i].isTag("label")) leg->AddText(Process[i]["label"].c_str());
      else                          leg->AddText(Process[i]["tag"].c_str());
      leg->Draw("same");
      ObjectToDelete.push_back(leg);
//      delete leg;
//      delete hist;

      T->Draw("same");

      string SavePath = SaveName + "_" + (Process[i])["tag"].toString() + plotExt;
      while(SavePath.find("*")!=std::string::npos)SavePath.replace(SavePath.find("*"),1,"");
      while(SavePath.find("#")!=std::string::npos)SavePath.replace(SavePath.find("#"),1,"");
      while(SavePath.find("{")!=std::string::npos)SavePath.replace(SavePath.find("{"),1,"");
      while(SavePath.find("}")!=std::string::npos)SavePath.replace(SavePath.find("}"),1,"");
      while(SavePath.find("(")!=std::string::npos)SavePath.replace(SavePath.find("("),1,"");
      while(SavePath.find(")")!=std::string::npos)SavePath.replace(SavePath.find(")"),1,"");
      while(SavePath.find("^")!=std::string::npos)SavePath.replace(SavePath.find("^"),1,"");
      while(SavePath.find("/")!=std::string::npos)SavePath.replace(SavePath.find("/"),1,"-");
      if(outDir.size()) SavePath = outDir +"/"+ SavePath;
      system(string(("rm -f ") + SavePath).c_str());
      c1->SaveAs(SavePath.c_str());
      delete c1;
   }

   for(unsigned int d=0;d<ObjectToDelete.size();d++){delete ObjectToDelete[d];}ObjectToDelete.clear();
   delete T;
}


void Draw2DHistogram(JSONWrapper::Object& Root, std::string RootDir, NameAndType HistoProperties){
   if(HistoProperties.isIndexPlot && cutIndex<0)return;

   std::string SaveName = "";

   TPaveText* T = new TPaveText(0.10,0.995,0.85,0.945, "NDC");
   T->SetFillColor(0);
   T->SetFillStyle(0);  T->SetLineColor(0);
   T->SetTextAlign(12);
   char Buffer[1024]; sprintf(Buffer, "CMS preliminary, #sqrt{s}=%.0f TeV, #scale[0.5]{#int} L=%.1f fb^{-1}", iEcm, iLumi/1000);
   T->AddText(Buffer);

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

      TH1* hist = NULL;
      std::vector<JSONWrapper::Object> Samples = (Process[i])["data"].daughters();
      for(unsigned int j=0;j<Samples.size();j++){
         double Weight = 1.0;
         if(!Process[i]["isdata"].toBool()  && !Process[i]["isdatadriven"].toBool())Weight*= iLumi;
         string filtExt("");
         // TEST // if(Process[i].isTag("mctruthmode") ) { char buf[255]; sprintf(buf,"_filt%d",(int)Process[i]["mctruthmode"].toInt()); filtExt += buf; }
         if(Samples[j].isTag("xsec")     )Weight*= Samples[j]["xsec"].toDouble();
         std::vector<JSONWrapper::Object> BR = Samples[j]["br"].daughters();for(unsigned int b=0;b<BR.size();b++){Weight*=BR[b].toDouble();}
         stSampleInfo& sampleInfo = sampleInfoMap[(Samples[j])["dtag"].toString()];
         Weight /= sampleInfo.initialNumberOfEvents;

         if(HistoProperties.name.find("puup"  )!=string::npos){Weight *= sampleInfo.PURescale_up;}
         if(HistoProperties.name.find("pudown")!=string::npos){Weight *= sampleInfo.PURescale_down;}         

         int split = 1;
	 if(!useMerged && Samples[j].isTag("split"))split = Samples[j]["split"].toInt();

         TH1* tmphist = NULL;
         for(int s=0;s<split;s++){
	   string segmentExt; if(split>1) { char buf[255]; sprintf(buf,"_%i",s); segmentExt += buf; } 
	    
            string FileName = RootDir + (Samples[j])["dtag"].toString() + ((Samples[j].isTag("suffix"))?(Samples[j])["suffix"].toString():string("")) + filtExt + segmentExt + filtExt + ".root";
            if(!FileExist[FileName])continue;
            TFile* File = new TFile(FileName.c_str());
            if(!File || File->IsZombie() || !File->IsOpen() || File->TestBit(TFile::kRecovered) )continue;
            TH1* tmptmphist = (TH1*) GetObjectFromPath(File,HistoProperties.name); 
            if(!tmptmphist){delete File;continue;}
            if(!tmphist){gROOT->cd(); tmphist = (TH1*)tmptmphist->Clone(tmptmphist->GetName());checkSumw2(tmphist);}else{tmphist->Add(tmptmphist);}
            //if(Process[i]["isdata"].toBool())printf("%s --> %f\n",(Samples[j])["dtag"].toString().c_str(), tmptmphist->Integral());
            delete tmptmphist;
            delete File;
         }
         if(!tmphist)continue;
         if(!hist){gROOT->cd(); hist = (TH1*)tmphist->Clone(tmphist->GetName());checkSumw2(hist);hist->Scale(Weight);}else{hist->Add(tmphist,Weight);}
         delete tmphist;
      }   
      if(!hist)continue;

      SaveName = hist->GetName();
      ObjectToDelete.push_back(hist);
      hist->SetTitle("");
      hist->SetStats(kFALSE);

      hist->Draw("COLZ");
  
      TPaveText* leg = new TPaveText(0.10,0.995,0.30,0.90, "NDC");
      leg->SetFillColor(0);
      leg->SetFillStyle(0);  leg->SetLineColor(0);
      leg->SetTextAlign(12);
      if(Process[i].isTag("label")) leg->AddText(Process[i]["label"].c_str());
      else                          leg->AddText(Process[i]["tag"].c_str());
      leg->Draw("same");
      ObjectToDelete.push_back(leg);
//      delete leg;
//      delete hist;
   }
   c1->cd(0);
   //T->Draw("same");
   string SavePath = SaveName + plotExt;
   while(SavePath.find("*")!=std::string::npos)SavePath.replace(SavePath.find("*"),1,"");
   while(SavePath.find("#")!=std::string::npos)SavePath.replace(SavePath.find("#"),1,"");
   while(SavePath.find("{")!=std::string::npos)SavePath.replace(SavePath.find("{"),1,"");
   while(SavePath.find("}")!=std::string::npos)SavePath.replace(SavePath.find("}"),1,"");
   while(SavePath.find("(")!=std::string::npos)SavePath.replace(SavePath.find("("),1,"");
   while(SavePath.find(")")!=std::string::npos)SavePath.replace(SavePath.find(")"),1,"");
   while(SavePath.find("^")!=std::string::npos)SavePath.replace(SavePath.find("^"),1,"");
   while(SavePath.find("/")!=std::string::npos)SavePath.replace(SavePath.find("/"),1,"-");
   if(outDir.size()) SavePath = outDir +"/"+ SavePath; 
   system(string(("rm -f ") + SavePath).c_str());
   c1->SaveAs(SavePath.c_str());
   for(unsigned int d=0;d<ObjectToDelete.size();d++){delete ObjectToDelete[d];}ObjectToDelete.clear();
   delete c1;
   delete T;
}


void Draw1DHistogram(JSONWrapper::Object& Root, std::string RootDir, NameAndType HistoProperties){
   if(HistoProperties.isIndexPlot && cutIndex<0)return;

   TCanvas* c1 = new TCanvas("c1","c1",800,800);
   TPad* t1 = new TPad("t1","t1", 0.0, 0.20, 1.0, 1.0);
   t1->Draw();
   t1->cd();
   if(!noLog) t1->SetLogy(true);
   if(logX)   t1->SetLogx(true);
   float maximumFound(false);//noLog);

   TLegend* legA  = new TLegend(0.845,0.2,0.99,0.99, "NDC"); 
   //   TLegend* legA  = new TLegend(0.51,0.93,0.67,0.75, "NDC"); 
   // TLegend* legB  = new TLegend(0.67,0.93,0.83,0.75, "NDC");
   THStack* stack = new THStack("MC","MC");
   TH1 *     mc   = NULL;
   TH1 *     mcPlusRelUnc = NULL;
   std::vector<TH1 *> spimpose;
   std::vector<TString> spimposeOpts;
   TH1*     data = NULL;
   std::vector<TObject*> ObjectToDelete;
   std::string SaveName = "";
   std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
   for(unsigned int i=0;i<Process.size();i++){
      if(Process[i]["isinvisible"].toBool())continue;
      TH1* hist = NULL;
      std::vector<JSONWrapper::Object> Samples = (Process[i])["data"].daughters();
      for(unsigned int j=0;j<Samples.size();j++){
         double Weight = 1.0;
         if(!Process[i]["isdata"].toBool() && !Process[i]["isdatadriven"].toBool() ){Weight*= iLumi; if(debug) cout << "Lumi is : " << iLumi << endl; }
         string filtExt("");
         // TEST // if(Process[i].isTag("mctruthmode") ) { char buf[255]; sprintf(buf,"_filt%d",(int)Process[i]["mctruthmode"].toInt()); filtExt += buf; }
         if(Samples[j].isTag("xsec")     )Weight*= Samples[j]["xsec"].toDouble();
         std::vector<JSONWrapper::Object> BR = Samples[j]["br"].daughters();for(unsigned int b=0;b<BR.size();b++){Weight*=BR[b].toDouble();}
         stSampleInfo& sampleInfo = sampleInfoMap[(Samples[j])["dtag"].toString()];
         Weight /= sampleInfo.initialNumberOfEvents;
	 //	 cout << (Samples[j])["dtag"].toString() << " " << (iLumi/1000)*(1/Weight) <<  " /fb" << endl;
         if(HistoProperties.name.find("puup"  )!=string::npos){Weight *= sampleInfo.PURescale_up  ;}
         if(HistoProperties.name.find("pudown")!=string::npos){Weight *= sampleInfo.PURescale_down;}
         if(HistoProperties.name.find("optim_cut")!=string::npos){Weight=1.0;}
         // TEST // string secondFiltExt=filtExt;
         string secondFiltExt("");
         // TEST // if(TString((Samples[j])["dtag"].toString()).Contains("filt5")) secondFiltExt="_filt5";
         // TEST // if(TString((Samples[j])["dtag"].toString()).Contains("filt6")) secondFiltExt="_filt6";

         int split = 1; 
	 if(!useMerged && Samples[j].isTag("split"))split = Samples[j]["split"].toInt();

         TH1* tmphist = NULL;
         for(int s=0;s<split;s++){
	   string segmentExt; if(split>1) { char buf[255]; sprintf(buf,"_%i",s); segmentExt += buf; }

	    string FileName = RootDir + (Samples[j])["dtag"].toString() + ((Samples[j].isTag("suffix"))?(Samples[j])["suffix"].toString():string("")) + filtExt + segmentExt + secondFiltExt + ".root";
            if(!FileExist[FileName]){
              //cout << FileName << " does not exist " << endl;
              continue;}
            TFile* File = new TFile(FileName.c_str());
            if(!File || File->IsZombie() || !File->IsOpen() || File->TestBit(TFile::kRecovered) ){ cout << FileName << " not recovered" << endl; continue;}
            TH1* tmptmphist = NULL; 
            if(HistoProperties.isIndexPlot && cutIndex>=0){
               TH2* tmp2D = (TH2*) GetObjectFromPath(File,HistoProperties.name);
               if(tmp2D){tmptmphist = tmp2D->ProjectionY((string(tmp2D->GetName())+cutIndexStr).c_str(),cutIndex,cutIndex); delete tmp2D;}
//            }else if(HistoProperties.name.find("mt_shape")){
//               TH2* tmp2D = (TH2*) GetObjectFromPath(File,HistoProperties.name);
//               if(tmp2D){tmp2D->GetXaxis()->SetTitle("#events"); tmptmphist = tmp2D->ProjectionX((string(tmp2D->GetName())+cutIndexStr).c_str()); delete tmp2D;}
            }else if(!HistoProperties.isIndexPlot){
               tmptmphist = (TH1*) GetObjectFromPath(File,HistoProperties.name);
	    }
	    if(!tmptmphist){
              //cout << FileName << " histogram " << HistoProperties.name << " does not exist " << endl; 
              delete File;continue;}
            if(!tmphist){gROOT->cd(); tmphist = (TH1*)tmptmphist->Clone(tmptmphist->GetName());checkSumw2(tmphist);}else{tmphist->Add(tmptmphist);}
            //if(Process[i]["isdata"].toBool())printf("%s --> %f\n",(Samples[j])["dtag"].toString().c_str(), tmptmphist->Integral());
            delete tmptmphist;
            delete File;
         }
         if(!tmphist)continue;
         if(debug) cout << "Histo " << tmphist->GetName() << " has integral " << tmphist->Integral() << ", entries " << tmphist->GetEntries() << " and will be weighted " << Weight << ", because lumi is " << iLumi << ", xsec is " << Samples[j]["xsec"].toDouble() << " and init number of events is " << sampleInfo.initialNumberOfEvents << ", yielding " << iLumi*Samples[j]["xsec"].toDouble()/sampleInfo.initialNumberOfEvents << endl;
         if(!hist){gROOT->cd(); hist = (TH1*)tmphist->Clone(tmphist->GetName());checkSumw2(hist);hist->Scale(Weight);}else{hist->Add(tmphist,Weight);}
         delete tmphist;
      }   
      if(!hist)continue;

      SaveName = hist->GetName();
      myPseudoDataName=hist->GetName();
      TString postfix(""); postfix+=i;
      hist->SetName(SaveName+postfix);
      if(Process[i].isTag("color" ) )hist->SetLineColor  ((int)Process[i][ "color"].toDouble()); else hist->SetLineColor  (1);
      if(Process[i].isTag("color" ) )hist->SetMarkerColor((int)Process[i][ "color"].toDouble()); else hist->SetMarkerColor(1);
      if(Process[i].isTag("color" ) )hist->SetFillColor  ((int)Process[i][ "color"].toDouble()); else hist->SetFillColor  (0);
      if(Process[i].isTag("lcolor") )hist->SetLineColor  ((int)Process[i]["lcolor"].toDouble());
      if(Process[i].isTag("mcolor") )hist->SetMarkerColor((int)Process[i]["mcolor"].toDouble()); 
      if(Process[i].isTag("fcolor") )hist->SetFillColor  ((int)Process[i]["fcolor"].toDouble()); 
      if(Process[i].isTag("lwidth") )hist->SetLineWidth  ((int)Process[i]["lwidth"].toDouble());// else hist->SetLineWidth  (1);
      if(Process[i].isTag("lstyle") )hist->SetLineStyle  ((int)Process[i]["lstyle"].toDouble());// else hist->SetLinStyle  (1);
      if(Process[i].isTag("fill"  ) )hist->SetFillColor  ((int)Process[i]["fill"  ].toDouble());
      if(Process[i].isTag("marker") )hist->SetMarkerStyle((int)Process[i]["marker"].toDouble());// else hist->SetMarkerStyle(1);
      
      if(std::string(hist->GetName()).find("eventflow") == std::string::npos )
        fixExtremities(hist,true,true);
      hist->SetTitle("");
      hist->SetStats(kFALSE);
      if( jodorStyle){
	//hist->SetMinimum(5e-2);
	//if(!Process[i]["isdata"].toBool()) hist->SetMinimum(hist->GetBinContent(hist->GetMinimumBin())*0.9);
	hist->SetMinimum(hist->GetBinContent(hist->GetMinimumBin())*0.90> 5e-2 ? hist->GetBinContent(hist->GetMinimumBin())*0.90 : 5e-2 );
      }
      else
	hist->SetMinimum(5e-2);      
      //hist->SetMaximum(1E6);
      hist->SetMaximum(hist->GetBinContent(hist->GetMaximumBin())*1.10);
      ObjectToDelete.push_back(hist);
      if(Process[i].isTag("normto")) hist->Scale( Process[i]["normto"].toDouble()/hist->Integral() );
      
      
      if((!Process[i].isTag("spimpose") || !Process[i]["spimpose"].toBool()) && !Process[i]["isdata"].toBool()){
         //Add to Stack
	stack->Add(hist, "HIST");               
	if(Process[i].isTag("label")) legA->AddEntry(hist, Process[i]["label"].c_str(), "F");	 
	else                          legA->AddEntry(hist, Process[i]["tag"].c_str(), "F");	 
	if(!mc){mc = (TH1D*)hist->Clone("mc");checkSumw2(mc);}else{mc->Add(hist);}
      }
      else if(Process[i].isTag("spimpose") && Process[i]["spimpose"].toBool())
	{
   	  //legB->AddEntry(hist, Process[i]["tag"].c_str(), "L");
	  if(Process[i].isTag("label")) legA->AddEntry(hist, Process[i]["label"].c_str(), Process[i]["isdata"].toBool() ? "P" : "L" );
	  else                          legA->AddEntry(hist, Process[i]["tag"].c_str(), Process[i]["isdata"].toBool() ? "P" : "L" );
	  spimposeOpts.push_back( Process[i]["isdata"].toBool() ? "e1" : "hist" );
	  spimpose.push_back(hist);
	  if(maximumFound<hist->GetMaximum()) maximumFound=hist->GetMaximum()*1.1;
	}
      else{
	if(Process[i]["isdata"].toBool()){
	  if(!data){
	    data = hist; 
	    if(Process[i].isTag("label")) legA->AddEntry(hist, Process[i]["label"].c_str(), "P"); 
	    else                          legA->AddEntry(hist, Process[i]["tag"].c_str(), "P"); 
	  }
	  else data->Add(hist);
	}
      }
   }
   if(mc &&  maximumFound<mc->GetMaximum()) maximumFound=mc->GetMaximum()*1.1;
   if(data && maximumFound<data->GetMaximum()) maximumFound=data->GetMaximum()*1.1;

   bool canvasIsFilled(false);
   if(stack && stack->GetStack() && stack->GetStack()->GetEntriesFast()>0){
     stack->Draw("");
     TH1 *hist=(TH1*)stack->GetStack()->At(0);
     if(stack->GetXaxis())
       {
	 stack->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());
	 stack->GetYaxis()->SetTitle(hist->GetYaxis()->GetTitle());
	 stack->SetMinimum(hist->GetMinimum());
	 stack->SetMaximum(maximumFound);
	 if(noLog)
	   {
	     stack->SetMaximum(maximumFound);
	     //stack->GetXaxis()->SetRangeUser(hist->GetMinimum(),maximumFound);
	   }
       }
     ObjectToDelete.push_back(stack);
     canvasIsFilled=true;

     if(showUnc && mc)
       {
	 mcPlusRelUnc = (TH1 *) mc->Clone("totalmcwithunc");checkSumw2(mcPlusRelUnc); mcPlusRelUnc->SetDirectory(0); 
	 for(int ibin=1; ibin<=mcPlusRelUnc->GetXaxis()->GetNbins(); ibin++)
	   {
	     Double_t error=sqrt(pow(mcPlusRelUnc->GetBinError(ibin),2)+pow(mcPlusRelUnc->GetBinContent(ibin)*baseRelUnc,2));
	     mcPlusRelUnc->SetBinError(ibin,error);
	   }
	 mcPlusRelUnc->SetFillStyle(3427);
	 mcPlusRelUnc->SetFillColor(kGray+1);
	 mcPlusRelUnc->SetMarkerStyle(1);
       }
   }
   if(data){
       data->Draw(canvasIsFilled ? "E1 same" : "E1");
       canvasIsFilled=true;
   }
   
   TH1* pseudoData = (TH1*) mc->Clone("pseudoData");
   pseudoData->Sumw2();
   if(generatePseudoData){
     GeneratePseudoData(pseudoData,mc);
     legA->AddEntry(pseudoData, "pseudo-data", "P");	 
     if(pseudoData){
       //pseudoData->SetLineColor(1);
       pseudoData->SetFillStyle  (0);
       //pseudoData->SetLineStyle (0);
       pseudoData->SetMarkerColor(1);
       pseudoData->SetMarkerStyle(20);
       pseudoData->Sumw2();
       pseudoData->Draw(canvasIsFilled ? "PE1 same" : "PE1");
       canvasIsFilled=true;
       myPseudoData= (TH1*) pseudoData->Clone("globalPseudoData"+myPseudoDataName);
     }
   }
   
   for(size_t ip=0; ip<spimpose.size(); ip++){
     TString opt=spimposeOpts[ip];
    spimpose[ip]->Draw(opt + (canvasIsFilled ? "same": "") );
     canvasIsFilled=true;
   }

   //compare data and MC
   if(showChi2)
     {
       TPaveText *pave = new TPaveText(0.6,0.85,0.8,0.9,"NDC");
       pave->SetBorderSize(0);
       pave->SetFillStyle(0);
       pave->SetTextAlign(32);
       pave->SetTextFont(42);
       char buf[100];
       if(data && mc && data->Integral()>0 && mc->Integral()>0)
	 {
	   sprintf(buf,"#chi^{2}/ndof : %3.2f", data->Chi2Test(mc,"WWCHI2/NDF") );
	   pave->AddText(buf);
	   sprintf(buf,"K-S prob: %3.2f", data->KolmogorovTest(mc));
	   pave->AddText(buf);

	 }
       else if(mc && spimpose.size()>0 && mc->Integral()>0)
	 {
	   for(size_t ip=0; ip<spimpose.size(); ip++)
	     {
	       if(spimpose[ip]->Integral()<=0) continue;
	       sprintf(buf,"#chi^{2}/ndof : %3.2f, K-S prob: %3.2f", spimpose[ip]->Chi2Test(mc,"WWCHI2/NDF"), spimpose[ip]->KolmogorovTest(mc) );
	       pave->AddText(buf);
	     }
	 }
       pave->Draw();
     }
   
   TPaveText* T = new TPaveText(0.1,0.995,0.84,0.95, "NDC");
   T->SetFillColor(0);
   T->SetFillStyle(0);  T->SetLineColor(0);
   T->SetTextAlign(12);
   char Buffer[1024]; 
   //if(isSim) sprintf(Buffer, "CMS simulation, #sqrt{s}=%.1f TeV, #scale[0.5]{#int} L=%.1f fb^{-1}", iEcm, iLumi/1000);
   if(isSim) sprintf(Buffer, "CMS simulation, #sqrt{s}=%.0f TeV", iEcm);
   else      sprintf(Buffer, "CMS preliminary, #sqrt{s}=%.0f TeV, #scale[0.5]{#int} L=%.1f fb^{-1}", iEcm, iLumi/1000);
   T->AddText(Buffer);
   T->Draw("same");
   T->SetBorderSize(0);

   
   legA->SetFillColor(0); legA->SetFillStyle(0); legA->SetLineColor(0);
   legA->SetHeader("");
   legA->Draw("same");
   legA->SetTextFont(42);
   //    legB->SetFillColor(0); legB->SetFillStyle(0); legB->SetLineColor(0);
   //    legB->SetHeader("");
   //    legB->Draw("same");
   //   legB->SetTextFont(42);
   
   //
   std::vector<TH1 *> compDists;
   if(data)                   compDists.push_back(data);
   else if(spimpose.size()>0) compDists=spimpose;
   if(mc && compDists.size())
     {
       c1->cd();
       TPad* t2 = new TPad("t2","t2", 0.0, 0.0, 1.0, 0.2);
       t2->Draw();
       t2->cd();
       t2->SetGridy(true);
       t2->SetPad(0,0.0,1.0,0.2);
       t2->SetTopMargin(0);
       t2->SetBottomMargin(0.5);
       if(logX)   t2->SetLogx(true);

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
       denRelUncH->Draw();
       denRelUnc->Draw("3");
       float yscale = (1.0-0.2)/(0.18-0);       
       denRelUncH->GetYaxis()->SetTitle("Data/#Sigma MC");
       denRelUncH->SetMinimum(0.4);
       denRelUncH->SetMaximum(1.6);
       denRelUncH->GetXaxis()->SetTitle("");
       //denRelUncH->SetMinimum(0);
       //denRelUncH->SetMaximum(data->GetBinContent(data->GetMaximumBin())*1.10);
       if(jodorStyle){
	 denRelUncH->GetXaxis()->SetTitleOffset(1.3);
	 denRelUncH->GetXaxis()->SetLabelOffset(0.05);
	 denRelUncH->GetXaxis()->SetLabelSize(0.050*yscale);
	 denRelUncH->GetXaxis()->SetTitleSize(0.042*yscale);
	 denRelUncH->GetXaxis()->SetTickLength(0.03*yscale);
       }
       else{
	 denRelUncH->GetXaxis()->SetTitleOffset(1.3);
	 denRelUncH->GetXaxis()->SetLabelSize(0.033*yscale);
	 denRelUncH->GetXaxis()->SetTitleSize(0.036*yscale);
	 denRelUncH->GetXaxis()->SetTickLength(0.03*yscale);
       }
       denRelUncH->GetYaxis()->SetTitleOffset(0.3);
       denRelUncH->GetYaxis()->SetNdivisions(5);
       denRelUncH->GetYaxis()->SetLabelSize(0.033*yscale);
       denRelUncH->GetYaxis()->SetTitleSize(0.036*yscale);
       
       //add comparisons
       for(size_t icd=0; icd<compDists.size(); icd++)
	 {
	   TString name("CompHistogram"); name+=icd;
	   TH1D *dataToObsH = (TH1D*)compDists[icd]->Clone(name);
	   checkSumw2(dataToObsH);
	   dataToObsH->Divide(mc);
	   dataToObsH->Draw("same");
	 }
     }
   else
     {
       //if not comparison resize the canvas
       c1->SetWindowSize(600,400);
       c1->SetCanvasSize(600,400);
       t1->SetPad(0,0,1,1);
     }

   c1->Modified();
   c1->Update();
   c1->cd();
   string SavePath = SaveName + plotExt;
   while(SavePath.find("*")!=std::string::npos)SavePath.replace(SavePath.find("*"),1,"");
   while(SavePath.find("#")!=std::string::npos)SavePath.replace(SavePath.find("#"),1,"");
   while(SavePath.find("{")!=std::string::npos)SavePath.replace(SavePath.find("{"),1,"");
   while(SavePath.find("}")!=std::string::npos)SavePath.replace(SavePath.find("}"),1,"");
   while(SavePath.find("(")!=std::string::npos)SavePath.replace(SavePath.find("("),1,"");
   while(SavePath.find(")")!=std::string::npos)SavePath.replace(SavePath.find(")"),1,"");
   while(SavePath.find("^")!=std::string::npos)SavePath.replace(SavePath.find("^"),1,"");
   while(SavePath.find("/")!=std::string::npos)SavePath.replace(SavePath.find("/"),1,"-");
   if(outDir.size()) SavePath = outDir +"/"+ SavePath;
   system(string(("rm -f ") + SavePath).c_str());
   c1->SaveAs(SavePath.c_str());
   delete c1;
   for(unsigned int d=0;d<ObjectToDelete.size();d++){delete ObjectToDelete[d];}ObjectToDelete.clear();
   delete legA;
   //   delete legB;
   delete T;
}


void ConvertToTex(JSONWrapper::Object& Root, std::string RootDir, NameAndType HistoProperties){
   if(HistoProperties.isIndexPlot && cutIndex<0)return;

   FILE* pFile = NULL;

   std::vector<TObject*> ObjectToDelete;
   TH1* stack = NULL; 
   std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
   for(unsigned int i=0;i<Process.size();i++){
      TH1* hist = NULL;
      std::vector<JSONWrapper::Object> Samples = (Process[i])["data"].daughters();




      for(unsigned int j=0;j<Samples.size();j++){
         double Weight = 1.0;
         if(!Process[i]["isdata"].toBool()  && !Process[i]["isdatadriven"].toBool())Weight*= iLumi;
         string filtExt("");
         // TEST // if(Process[i].isTag("mctruthmode") ) { char buf[255]; sprintf(buf,"_filt%d",(int)Process[i]["mctruthmode"].toInt()); filtExt += buf; }
         if(Samples[j].isTag("xsec")     )Weight*= Samples[j]["xsec"].toDouble();
         std::vector<JSONWrapper::Object> BR = Samples[j]["br"].daughters();for(unsigned int b=0;b<BR.size();b++){Weight*=BR[b].toDouble();}
         stSampleInfo& sampleInfo = sampleInfoMap[(Samples[j])["dtag"].toString()];
         Weight /= sampleInfo.initialNumberOfEvents;

         if(HistoProperties.name.find("puup"  )!=string::npos){Weight *= sampleInfo.PURescale_up  ;}
         if(HistoProperties.name.find("pudown")!=string::npos){Weight *= sampleInfo.PURescale_down;}


         int split = 1;
	 if(!useMerged && Samples[j].isTag("split"))split = Samples[j]["split"].toInt();

         TH1* tmphist = NULL;
         for(int s=0;s<split;s++){
	   string segmentExt; if(split>1) { char buf[255]; sprintf(buf,"_%i",s); segmentExt += buf; }

            string FileName = RootDir + (Samples[j])["dtag"].toString() + ((Samples[j].isTag("suffix"))?(Samples[j])["suffix"].toString():string("")) + filtExt + segmentExt + filtExt + ".root";
            if(!FileExist[FileName])continue;
            TFile* File = new TFile(FileName.c_str());
            if(!File || File->IsZombie() || !File->IsOpen() || File->TestBit(TFile::kRecovered) )continue;
            TH1* tmptmphist = NULL;
            if(HistoProperties.isIndexPlot && cutIndex>=0){
               TH2* tmp2D = (TH2*) GetObjectFromPath(File,HistoProperties.name);
               if(tmp2D){tmptmphist = tmp2D->ProjectionY("_py",cutIndex,cutIndex); delete tmp2D;}
            }else if(!HistoProperties.isIndexPlot){
               tmptmphist = (TH1*) GetObjectFromPath(File,HistoProperties.name);
            }
            if(!tmptmphist){delete File;continue;}
            if(!tmphist){gROOT->cd(); tmphist = (TH1*)tmptmphist->Clone(tmptmphist->GetName());checkSumw2(tmphist);}else{tmphist->Add(tmptmphist);}
            delete tmptmphist;
            delete File;
         }
         if(!tmphist)continue;
         if(!hist){gROOT->cd(); hist = (TH1*)tmphist->Clone(tmphist->GetName());checkSumw2(hist);hist->Scale(Weight);}else{hist->Add(tmphist,Weight);}
         delete tmphist;
      }
      if(!hist)continue;

      if(!pFile){
         string SavePath = string(hist->GetName()) + ".tex";
         while(SavePath.find("*")!=std::string::npos)SavePath.replace(SavePath.find("*"),1,"");
         while(SavePath.find("#")!=std::string::npos)SavePath.replace(SavePath.find("#"),1,"");
         while(SavePath.find("{")!=std::string::npos)SavePath.replace(SavePath.find("{"),1,"");
         while(SavePath.find("}")!=std::string::npos)SavePath.replace(SavePath.find("}"),1,"");
         while(SavePath.find("(")!=std::string::npos)SavePath.replace(SavePath.find("("),1,"");
         while(SavePath.find(")")!=std::string::npos)SavePath.replace(SavePath.find(")"),1,"");
         while(SavePath.find("^")!=std::string::npos)SavePath.replace(SavePath.find("^"),1,"");
         while(SavePath.find("/")!=std::string::npos)SavePath.replace(SavePath.find("/"),1,"-");
         SavePath = outDir + SavePath;
         system(string(("rm -f ") + SavePath).c_str());
         pFile = fopen(SavePath.c_str(), "w");

         fprintf(pFile, "\\begin{table}[htp]\n");
         fprintf(pFile, "\\begin{center}\n");
         fprintf(pFile, "\\caption{}\n");
         fprintf(pFile, "\\label{tab:table}\n");

         string colfmt = "|l";  for(int b=1;b<=hist->GetXaxis()->GetNbins();b++){colfmt = colfmt + "c";} colfmt+="|";
         string colname = "";   for(int b=1;b<=hist->GetXaxis()->GetNbins();b++){
            std::string tempcolname =  hist->GetXaxis()->GetBinLabel(b);
            if(tempcolname.find("_")!=std::string::npos || tempcolname.find("^")!=std::string::npos)tempcolname = string("$") + tempcolname + "$";
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


void clean_json(JSONWrapper::Object& Root, std::vector<string>& dtagstokeep){
   if(dtagstokeep.size()<=0)return;

  //loop on the json list of dtags and remove the one that are not in the list
   std::vector<JSONWrapper::Object>& Process = Root["proc"].daughters();
   for(size_t ip=0; ip<Process.size(); ip++){
       string tag = Process[ip].getString("tag");
       //search if this tag is in the vector of object to keep
       bool drop=true;
       for(unsigned int s=0;s<dtagstokeep.size();s++){ if(tag.find(dtagstokeep[s])!=std::string::npos){drop=false; break;} }
       if(drop){
          Root["proc"].obj.erase(Root["proc"].obj.begin()+ip); 
          Root["proc"].key.erase(Root["proc"].key.begin()+ip);
          ip--;
       }else{
          printf("Keep the sample %s\n", tag.c_str());
       }
   }
}


int main(int argc, char* argv[]){
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

   std::vector<string> histoNameMask;
   std::vector<string> histoNameMaskStart;
   std::vector<string> dtagstokeep;

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
        printf("--onlyContain    --> processing only the objects containing the following argument in their name\n");
        printf("--onlyStartWith  --> processing only the objects starting with the following argument in their name\n");
        printf("--index   --> will do the projection on that index for histos of type cutIndex\n");
        printf("--chi2    --> show the data/MC chi^2\n"); 
        printf("--showUnc --> show stat uncertainty (if number is given use it as relative bin by bin uncertainty (e.g. lumi)\n"); 
	printf("--noLog   --> use linear scale\n");
	printf("--logX    --> use log scale on X\n");
        printf("--no1D   --> Skip processing of 1D objects\n");
        printf("--no2D   --> Skip processing of 2D objects\n");
        printf("--noTex  --> Do not create latex table (when possible)\n");
        printf("--noPowers --> Do not use powers of 10 for numbers in tables\n");
        printf("--noRoot --> Do not make a summary .root file\n");
        printf("--noPlot --> Do not creates plot files (useful to speedup processing)\n");
	printf("--plotExt --> extension to save\n");
	printf("--cutflow --> name of the histogram with the original number of events (cutflow by default)\n");
        printf("--splitCanvas --> (only for 2D plots) save all the samples in separated pltos\n");
        printf("--forceMerge --> merge splitted samples\n");
	printf("--useMerged --> use merged splitted samples\n");
	printf("--jodorStyle --> use plotting style requested from Jodor");
        printf("--generatePseudoData --> generate pseudo-data by poisson-smearing the total SM MC distribution");
        
        printf("command line example: runPlotter --json ../data/beauty-samples.json --iLumi 2007 --inDir OUT/ --outDir OUT/plots/ --outFile plotter.root --noRoot --noPlot\n");
	return 0;
     }

     if(arg.find("--tag"   )!=string::npos && i+1<argc){ dtagstokeep.push_back(argv[i+1]); i++; printf("Tags containing '%s' will be kept but not the others\n", argv[i]);    }
     if(arg.find("--iLumi"  )!=string::npos && i+1<argc){ sscanf(argv[i+1],"%lf",&iLumi); i++; printf("Lumi = %f\n", iLumi); }
     if(arg.find("--iEcm"  )!=string::npos && i+1<argc){ sscanf(argv[i+1],"%lf",&iEcm); i++; printf("Ecm = %f TeV\n", iEcm); }

     if(arg.find("--inDir"  )!=string::npos && i+1<argc){ inDir    = argv[i+1];  i++;  printf("inDir = %s\n", inDir.c_str());  }
     if(arg.find("--outDir" )!=string::npos && i+1<argc){ outDir   = argv[i+1];  i++;  printf("outDir = %s\n", outDir.c_str());  }
     if(arg.find("--outFile")!=string::npos && i+1<argc){ outFile  = argv[i+1];  i++; printf("output file = %s\n", outFile.c_str()); }
     if(arg.find("--json"   )!=string::npos && i+1<argc){ jsonFile = argv[i+1];  i++;  }
     if(arg.find("--onlyStartWith"   )!=string::npos && i+1<argc){ histoNameMaskStart.push_back(argv[i+1]); printf("Only process Histo starting with '%s'\n", argv[i+1]); i++;  }
     if(arg.find("--onlyContain"   )!=string::npos && i+1<argc)         { histoNameMask.push_back(argv[i+1]); printf("Only process Histo containing '%s'\n", argv[i+1]); i++;  }
     if(arg.find("--index"  )!=string::npos && i+1<argc)         { sscanf(argv[i+1],"%d",&cutIndex); i++; onlyCutIndex=(cutIndex>=0); printf("index = %i\n", cutIndex);  }
     if(arg.find("--chi2"  )!=string::npos)                      { showChi2 = true;  }
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
     if(arg.find("--forceMerge")!=string::npos){ forceMerge = true;  useMerged=true;  }
     if(arg.find("--useMerged")!=string::npos){ useMerged = true;    }
     if(arg.find("--jodorStyle")!=string::npos){ jodorStyle = true;  }
     if(arg.find("--generatePseudoData")!=string::npos){ generatePseudoData = true; }
     if(arg.find("--noLog")!=string::npos){ noLog = true;    }
     if(arg.find("--logX")!=string::npos){ logX = true;    }
     if(arg.find("--no2D"  )!=string::npos){ do2D = false;    }
     if(arg.find("--no1D"  )!=string::npos){ do1D = false;    }
     if(arg.find("--noTex" )!=string::npos){ doTex= false;    }
     if(arg.find("--noPowers" )!=string::npos){ doPowers= false;    }
     if(arg.find("--noRoot")!=string::npos){ StoreInFile = false;    }
     if(arg.find("--noPlot")!=string::npos){ doPlot = false;    }
     if(arg.find("--plotExt" )!=string::npos && i+1<argc){ plotExt   = argv[i+1];  i++;  printf("saving plots as = %s\n", plotExt.c_str());  }
     if(arg.find("--cutflow" )!=string::npos && i+1<argc){ cutflowhisto   = argv[i+1];  i++;  printf("Normalizing from 1st bin in = %s\n", cutflowhisto.c_str());  }
     if(arg.find("--splitCanvas")!=string::npos){ splitCanvas = true;    }
   } 
   system( (string("mkdir -p ") + outDir).c_str());

   char buf[255];
   sprintf(buf, "_Index%d", cutIndex);
   cutIndexStr = buf;

   JSONWrapper::Object Root(jsonFile, true);
   if(forceMerge) MergeSplittedSamples(Root,inDir);
   GetInitialNumberOfEvents(Root,inDir,NameAndType(cutflowhisto,true, false));  //Used to get the rescale factor based on the total number of events geenrated
   std::list<NameAndType> histlist;
   GetListOfObject(Root,inDir,histlist);
   histlist.sort();
   histlist.unique();   
   clean_json(Root, dtagstokeep);

   TFile* OutputFile = NULL;
   if(StoreInFile) OutputFile = new TFile(outFile.c_str(),"RECREATE");
   printf("Progressing Bar              :0%%       20%%       40%%       60%%       80%%       100%%\n");
   printf("                             :");
   int TreeStep = histlist.size()/50;if(TreeStep==0)TreeStep=1;
   string csvFile(outDir +"/histlist.csv");
   system(("echo \"\" > " + csvFile).c_str());
   int ictr(0);
   for(std::list<NameAndType>::iterator it= histlist.begin(); it!= histlist.end(); it++,ictr++)
     {
       if(ictr%TreeStep==0){printf(".");fflush(stdout);}
       bool passMasking = false;
       for(unsigned int i=0;i<histoNameMask.size();i++){if(it->name.find(histoNameMask[i])!=std::string::npos)passMasking=true;}
       for(unsigned int i=0;i<histoNameMaskStart.size();i++){if(it->name.find(histoNameMaskStart[i])==0)passMasking=true;}
       if(histoNameMask.size()==0 && histoNameMaskStart.size()==0)passMasking = true;
       if(!passMasking)continue;

       system(("echo \"" + it->name + "\" >> " + csvFile).c_str());

       if(doTex && (it->name.find("eventflow")!=std::string::npos || it->name.find("evtflow")!=std::string::npos) && it->name.find("optim_eventflow")==std::string::npos){    ConvertToTex(Root,inDir,*it); }
       if(doPlot && do2D  && !it->type){                      if(!splitCanvas){Draw2DHistogram(Root,inDir,*it); }else{Draw2DHistogramSplitCanvas(Root,inDir,*it);}}
       if(doPlot && do1D  &&  it->type){                                       Draw1DHistogram(Root,inDir,*it); }
      
       if(StoreInFile && do2D  && !it->type){                                  SavingToFile(Root,inDir,*it, OutputFile); }
       if(StoreInFile && do1D  &&  it->type){                                  SavingToFile(Root,inDir,*it, OutputFile); }
       if(generatePseudoData)                                           SavingPseudoDataToFile(Root,inDir,*it, OutputFile);

     }printf("\n");
   if(StoreInFile) OutputFile->Close();
   
   //system(("python ${CMSSW_BASE}/src/CMGTools/HtoZZ2l2nu/data/html/generateJSONplotterFromList.py -i " + csvFile + " -o "+outDir+"/plotter.json").c_str());
   system(("rm " + csvFile).c_str());
   //system(("cp ${CMSSW_BASE}/src/CMGTools/HtoZZ2l2nu/data/html/index.html " + outDir).c_str());
   //printf("You can browse the results using %s/index.html\n",outDir.c_str());
   printf("runPlotter : All done\n");
}

