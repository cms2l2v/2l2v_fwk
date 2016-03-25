#include <vector>
#include <list>
#include <iterator>
#include <algorithm>
#include <unordered_map>
#include <iostream>
#include <string>
#include <regex>

#include "UserCode/llvv_fwk/interface/tdrstyle.h"
#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/RootUtils.h"
#include "UserCode/llvv_fwk/interface/JSONWrapper.h"

using namespace std;


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

   //FIXME these are currently hard-coded, but we should pass them via parameters using the code just bellow
   bool isData=true;
   int  rbin=1;
   bool asym=true;
   string DataFilePath  = "plotter_2016_01_18_forPhotonWeights.root";
   string outDir    = "photonWeights/";
   string outFile   = outDir + "photonWeights.root";

   for(int i=1;i<argc;i++){
     string arg(argv[i]);
     if(arg.find("--help")!=string::npos){
        printf("--help   --> print this helping text\n");
        printf("--inFile  --> path of the plotter .root file with inputs\n");
        printf("--outDir  --> path of the output directory to save the plots\n");
        printf("--outFile --> path of the output summary .root file\n");
	return 0;
     }
     if(arg.find("--inFile" )!=string::npos && i+1<argc){ DataFilePath = argv[i+1];  i++; printf("input file = %s\n", DataFilePath.c_str()); }
     if(arg.find("--outFile")!=string::npos && i+1<argc){ outFile  = argv[i+1];  i++; printf("output file = %s\n", outFile.c_str()); }
     if(arg.find("--outDir" )!=string::npos && i+1<argc){ outDir   = argv[i+1];  i++; printf("outDir = %s\n", outDir.c_str());  }
   }
   system( (string("mkdir -p ") + outDir).c_str());


  std::vector<string> fDataDir;  //fake data
  std::vector<string> tDataDir;  //true lepton data (contamination to be subtracted)
  if(isData){
     fDataDir.push_back("data");
     tDataDir.push_back("ZH#rightarrow ll#tau#tau");
     tDataDir.push_back("ZZ#rightarrow 4l");
     tDataDir.push_back("WZ#rightarrow 2l+X");
     tDataDir.push_back("VVV");
  }else{
     fDataDir.push_back("ZH#rightarrow ll#tau#tau");
     fDataDir.push_back("ZZ#rightarrow 4l");
     fDataDir.push_back("ZZ#rightarrow 2l+X");
     fDataDir.push_back("WW#rightarrow l#nu+X");
     fDataDir.push_back("WZ#rightarrow 2l+X");
     fDataDir.push_back("VVV");
     fDataDir.push_back("Top");
     fDataDir.push_back("W#rightarrow l#nu");
     fDataDir.push_back("Z#rightarrow ll");
     tDataDir.push_back("ZH#rightarrow ll#tau#tau");
     tDataDir.push_back("ZZ#rightarrow 4l");
     tDataDir.push_back("WZ#rightarrow 2l+X");
     tDataDir.push_back("VVV");
     tDataDir.push_back("Z#rightarrow ll");
  }


  std::vector<int> catColor = {4     , 2   , 8  };
  std::vector<string> catL  = {"Fake e", "Fake #mu", "Fake #tau_{had}"};
  std::vector<string> binL  = {"Inc.", "Barrel", "Endcap", "Inc. (m_{T}>30)", "Barrel (m_{T}>30)", "Endcap (m_{T}>30)"};

  std::vector<string> cat   = {"FR_El", "FR_Mu", "FR_Ta"};
  std::vector<string> bin   = {"","_B","_E", "_TMCut", "_TMCut_B", "_TMCut_E"};
  std::vector<string> var   = {"", "_Id_Iso01", "_Id_Iso02", "_Id_Iso03", "_Id_IsoLo", "_Id_IsoMe"};
  std::vector<string> wrt   = {"_wrtJetPt", "_wrtLepPt"};



  std::map<string, TH1D*> DataHistos;
  for(unsigned int c=0;c<cat.size();c++){
  for(unsigned int b=0;b<bin.size();b++){
  for(unsigned int v=0;v<var.size();v++){
  for(unsigned int w=0;w<wrt.size();w++){
     string& DataFile             = DataFilePath;
     std::vector<string>& DataDir = fDataDir;



     TFile* File = new TFile(DataFile.c_str(), "READ");

     if(!File || File->IsZombie() || !File->IsOpen() || File->TestBit(TFile::kRecovered) ){
        printf("File %s coudn't be opened\n", DataFile.c_str());
     }

     for(unsigned int d=0;d<DataDir.size();d++){
        TH1D* hist = (TH1D*) File->Get((DataDir[d]+"/"+cat[c]+var[v]+bin[b]+wrt[w]).c_str());        
        if(!hist)continue;



        hist->GetSumw2();
        if(rbin>1){ hist->Rebin(rbin);  hist->Scale(1.0, "width"); }
        hist->GetXaxis()->SetTitleSize(.055);
        hist->GetYaxis()->SetTitleSize(.055);
        hist->GetXaxis()->SetLabelSize(.05);
        hist->GetYaxis()->SetLabelSize(.05);
        hist->SetFillColor(0);
        hist->SetLineColor(catColor[c]);
        hist->SetMarkerColor(catColor[c]);
        hist->SetStats(kFALSE);



        if(DataHistos.find(cat[c]+var[v]+bin[b])==DataHistos.end()){
           gROOT->cd(); //make sure that the file is saved in memory and not in file
           DataHistos[cat[c]+var[v]+bin[b]+wrt[w]] = (TH1D*)hist->Clone();; //create a new histo, since it's not found
        }else{
           DataHistos[cat[c]+var[v]+bin[b]+wrt[w]]->Add(hist); //add to existing histogram
        }



     }
     File->Close();
  }}}}//all histos are now loaded


  //compute the FR
  for(unsigned int w=0;w<wrt.size();w++){
  for(unsigned int v=0;v<var.size();v++){
     if(var[v]=="" || var[v].find("weight")!=std::string::npos)continue;
     for(unsigned int b=0;b<bin.size();b++){     
        for(unsigned int c=0;c<cat.size();c++){           
           if(DataHistos.find(cat[c]+""+bin[b]+wrt[w])==DataHistos.end()){printf("Histo missing for %s\n", (cat[c]+""+bin[b]+wrt[w]).c_str()); continue;}
           TH1D* histDen = DataHistos[cat[c]+""+bin[b]+wrt[w] ];

           if(DataHistos.find(cat[c]+var[v]+bin[b]+wrt[w])==DataHistos.end()){printf("Histo missing for %s\n", (cat[c]+var[v]+bin[b]+wrt[w]).c_str()); continue;}
           TH1D* histNum = DataHistos[cat[c]+var[v]+bin[b]+wrt[w]];
           histNum = (TH1D*)histNum->Clone((cat[c]+var[v]+"weight"+bin[b]+wrt[w]).c_str());
           histNum->Divide(histDen);
           DataHistos[cat[c]+var[v]+"weight"+bin[b]+wrt[w]] = histNum;
        }
     }
     printf("add variable %s\n", (var[v]+"weight").c_str());
     var.push_back(var[v]+"weight");
  }}


  //make the plots
  for(unsigned int w=0;w<wrt.size();w++){
  for(unsigned int v=0;v<var.size();v++){
     double xmin,xmax;
     double ymin=0.5, ymax=1E6;
     xmin=1.; xmax=1000.0;   
     if(var[v].find("weight")!=std::string::npos){ymin=0;  ymax=1.0;}

//     TCanvas* c1 =new TCanvas("c1","c1",500*bin.size(), 500);
//     c1->Divide(bin.size(),1);

     for(unsigned int b=0;b<bin.size();b++){
//        c1->cd(1+b);
        TCanvas* c1 =new TCanvas("c1","c1",500, 500);       
        if(var[v].find("weight")==std::string::npos) gPad->SetLogy();
        TLegend* leg = new TLegend(0.52,0.67,0.92,0.90);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.04);
        leg->SetTextFont(42);
        leg->SetHeader(binL[b].c_str());

        bool axisDrawn=false;
        for(unsigned int c=0;c<cat.size();c++){

           if(DataHistos.find(cat[c]+var[v]+bin[b]+wrt[w])==DataHistos.end()){printf("Histo missing for %s\n", (cat[c]+var[v]+bin[b]+wrt[w]).c_str()); continue;}
           TH1D* hist = DataHistos[cat[c]+var[v]+bin[b]+wrt[w]];
           hist->GetXaxis()->SetRangeUser(xmin,xmax);
           hist->GetYaxis()->SetRangeUser(ymin,ymax);

           if(!axisDrawn){hist->Draw("AXIS"); axisDrawn=true;}
           hist->Draw("HIST E1 same");           
           leg->AddEntry(hist, catL[c].c_str(), "LP"); 
        }
        leg->Draw("SAME");
        c1->SaveAs((outDir + "/LeptonFR"+var[v]+bin[b]+wrt[w]+".png").c_str());
     }
  }}

  //SAVE THE material to the weight file
  TFile* OutputFile = new TFile(outFile.c_str(),"RECREATE");
  OutputFile->cd();
  for(unsigned int c=0;c<cat.size();c++){
     for(unsigned int b=0;b<bin.size();b++){
     for(unsigned int v=0;v<var.size();v++){
     if(var[v].find("weight")==std::string::npos)continue;     
     for(unsigned int w=0;w<wrt.size();w++){
        if(DataHistos.find(cat[c]+var[v]+bin[b]+wrt[w])==DataHistos.end()){printf("Histo missing for %s\n", (cat[c]+var[v]+bin[b]+wrt[w]).c_str()); continue;}
        TH1D* hist = DataHistos[cat[c]+var[v]+bin[b]+wrt[w]];

        TGraphErrors* graph = new TGraphErrors(hist);
        graph->Write((cat[c]+"FRWeights"+bin[b]+wrt[w]).c_str());
        }}
     }
  }
  OutputFile->Write(); OutputFile->Close();
}

