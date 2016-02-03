
#include <string>
#include <vector>

#include "TROOT.h"
#include "TStyle.h"
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
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TPaveText.h"

#include "UserCode/llvv_fwk/interface/tdrstyle.h"
#include "UserCode/llvv_fwk/src/tdrstyle.C"
#include "UserCode/llvv_fwk/interface/RootUtils.h"


std::map<double, double> xsecXbr;
void getXSecXBR(TString inputXSec){
   if(inputXSec=="")return;
   FILE* pFile = fopen(inputXSec.Data(),"r");
   if(!pFile){printf("Can't open %s\n",inputXSec.Data()); return;}
//   while(true){
    for(int n=0;n<1000;n++){
      float mass; double xsec, br;
      if(fscanf(pFile,"%f %lf %lf\n",&mass, &xsec, &br)==EOF)break;
      //printf("%f %f %f\n", mass, xsec, br);
      xsecXbr[mass] = xsec*br;
   }
   fclose(pFile);
}

//tree variables
double Tmh, Tlimit, TlimitErr; float TquantExp;

void fillLimitArray(TTree* tree, float Quantil, double* limitArray, double* massArray=NULL){
  int i=0;
  for(int ientry=0;ientry<tree->GetEntriesFast();ientry++){
     tree->GetEntry(ientry);
     if(TquantExp==Quantil){
        //printf("Quantil = %f - mH=%f --> %f\n",TquantExp,Tmh,Tlimit);
        limitArray[i]=Tlimit;
        if(massArray){massArray[i]=Tmh;} 
        i++;
     }
  }
}


TCutG* GetErrorBand(string name, int N, double* Mass, double* Low, double* High)
{
   TCutG* cutg = new TCutG(name.c_str(),2*N+2);
   cutg->SetFillColor(kGreen-7);
   cutg->SetLineStyle(0);
   cutg->SetLineColor(0);
   int I = 0;
   for(int i=0;i<N;i++){cutg->SetPoint(I,Mass[i]    , Low[i]     );I++; }
                        cutg->SetPoint(I,Mass[N-1]  , Low[N-1]   );I++;
                        cutg->SetPoint(I,Mass[N-1]  , High[N-1]  );I++;
   for(int i=0;i<N;i++){cutg->SetPoint(I,Mass[N-1-i], High[N-1-i]);I++;}
   return cutg;
}

void scaleLimitsByXSecXBR(int N, double* Mass, double* Limit){
   for(int i=0;i<N;i++){
      if(xsecXbr.find(Mass[i])==xsecXbr.end())continue; //skip points missing
      Limit[i] = Limit[i]*xsecXbr[Mass[i]];
   }   
}

void printLimits(FILE* pFile, TGraph* graph, double Mmin=200, double Mmax=600){
   double previous = graph->Eval(Mmin);
   fprintf(pFile, "Exclude ");
   bool opened = false;
   for(double M=Mmin;M<=Mmax;M+=1.0){
      double NEW = graph->Eval(M);
      //printf("%04f - %6.3f\n", M, NEW);
      if(previous>1 && NEW<=1){fprintf(pFile,"[%f,",M); opened=true;}
      if(previous<1 && NEW>=1){fprintf(pFile,"%f]",M-1); opened=false;}
      if(M==Mmax    && opened){fprintf(pFile,">%f]",Mmax); opened=false;}
      previous = NEW;
   }fprintf(pFile,"\n");
}



void plotLimit(string outputDir="./", TString inputs="", TString inputXSec="", bool strengthLimit=true, bool blind=false, double energy=7, double luminosity=5.035, TString legendName="ee and #mu#mu channels")
{
   setTDRStyle();  
   gStyle->SetPadTopMargin   (0.05);
   gStyle->SetPadBottomMargin(0.12);
   gStyle->SetPadRightMargin (0.16);
   gStyle->SetPadLeftMargin  (0.14);
   gStyle->SetTitleSize(0.04, "XYZ");
   gStyle->SetTitleXOffset(1.1);
   gStyle->SetTitleYOffset(1.45);
   gStyle->SetPalette(1);
   gStyle->SetNdivisions(505);

  string suffix = outputDir;
  double cprime=1.0; double  brnew=0.0;
  double XSecScaleFactor = 1.0;
  if(suffix.find("_cp")!=string::npos){
     sscanf(suffix.c_str()+suffix.find("_cp"), "_cp%lf_brn%lf", &cprime, &brnew);
//     XSecScaleFactor = pow(cprime,2) * (1-brnew);
  }
  //XSecScaleFactor = 1.0/1000; //pb to fb

 
  //get xsec * br from summary file
  getXSecXBR(inputXSec); 
  //get the limits from the tree
  TFile* file = TFile::Open(inputs);
  printf("Looping on %s\n",inputs.Data());
  if(!file) return;
  if(file->IsZombie()) return;
  TTree* tree = (TTree*)file->Get("limit");
  tree->GetBranch("mh"              )->SetAddress(&Tmh      );
  tree->GetBranch("limit"           )->SetAddress(&Tlimit   );
  tree->GetBranch("limitErr"        )->SetAddress(&TlimitErr);
  tree->GetBranch("quantileExpected")->SetAddress(&TquantExp);
  unsigned int N = tree->GetEntriesFast() / 6 ;// 6Limits per mass point (observed, meand , +-1sigma, +-2sigma)
  double* MassAxis   = new double[N];
  double* ObsLimit   = new double[N];  fillLimitArray(tree,-1   ,ObsLimit,MassAxis);  if(!strengthLimit)scaleLimitsByXSecXBR(N, MassAxis, ObsLimit);
  double* ExpLimitm2 = new double[N];  fillLimitArray(tree,0.025,ExpLimitm2);         if(!strengthLimit)scaleLimitsByXSecXBR(N, MassAxis, ExpLimitm2);
  double* ExpLimitm1 = new double[N];  fillLimitArray(tree,0.160,ExpLimitm1);         if(!strengthLimit)scaleLimitsByXSecXBR(N, MassAxis, ExpLimitm1);
  double* ExpLimit   = new double[N];  fillLimitArray(tree,0.500,ExpLimit  );         if(!strengthLimit)scaleLimitsByXSecXBR(N, MassAxis, ExpLimit);
  double* ExpLimitp1 = new double[N];  fillLimitArray(tree,0.840,ExpLimitp1);         if(!strengthLimit)scaleLimitsByXSecXBR(N, MassAxis, ExpLimitp1);
  double* ExpLimitp2 = new double[N];  fillLimitArray(tree,0.975,ExpLimitp2);         if(!strengthLimit)scaleLimitsByXSecXBR(N, MassAxis, ExpLimitp2);
  file->Close();



  //get the pValue
  inputs = inputs.ReplaceAll("/LimitTree.root", "/PValueTree.root");
  file = TFile::Open(inputs);
  printf("Looping on %s\n",inputs.Data());
  if(!file) return;
  if(file->IsZombie()) return;
  tree = (TTree*)file->Get("limit");
  tree->GetBranch("limit"           )->SetAddress(&Tlimit   );
  double* pValue     = new double[N];  fillLimitArray(tree,-1   ,pValue);
  file->Close();


  //make TH Cross-sections
  double* ThXSec   = new double[N]; for(unsigned int i=0;i<N;i++){ThXSec[i] = xsecXbr[MassAxis[i]];} 

  //scale TH cross-section and limits according to scale factor 
  //this only apply to NarrowResonnance case
  for(unsigned int i=0;i<N;i++){
    if(strengthLimit){
       ObsLimit[i]  /= XSecScaleFactor;
       ExpLimitm2[i]/= XSecScaleFactor;
       ExpLimitm1[i]/= XSecScaleFactor;
       ExpLimit  [i]/= XSecScaleFactor;
       ExpLimitp1[i]/= XSecScaleFactor;
       ExpLimitp2[i]/= XSecScaleFactor;
    }
    ThXSec[i]    *= XSecScaleFactor;
  }

  //FIXME
  strengthLimit=false;

  string prod = "pp";
  if(outputDir.find("GGF")!=std::string::npos)prod="gg";
  if(outputDir.find("VBF")!=std::string::npos)prod="qq";

  //limits in terms of signal strength
  TCanvas* c = new TCanvas("c", "c",800,800);
  TH1F* framework = new TH1F("Graph","Graph",1,strengthLimit?150:150,1050);
  framework->SetStats(false);
  framework->SetTitle("");
  framework->GetXaxis()->SetTitle("Higgs boson mass (GeV)");
  framework->GetYaxis()->SetTitleOffset(1.70);
  if(strengthLimit){
  framework->GetYaxis()->SetTitle("#mu = #sigma_{95%} / #sigma_{th}");
  framework->GetYaxis()->SetRangeUser(1E-1,1E4);
  c->SetLogy(true);
  }else{
  framework->GetYaxis()->SetTitle((string("#sigma_{95%} (") + prod +" #rightarrow H #rightarrow ZZ) (fb)").c_str());
  framework->GetYaxis()->SetRangeUser(1E1,1E4);
  c->SetLogy(true);
  }
  framework->GetXaxis()->SetLabelOffset(0.007);
  framework->GetXaxis()->SetLabelSize(0.03);
  framework->GetXaxis()->SetTitleOffset(1.0);
  framework->GetXaxis()->SetTitleFont(42);
  framework->GetXaxis()->SetTitleSize(0.035);
  framework->GetYaxis()->SetLabelFont(42);
  framework->GetYaxis()->SetLabelOffset(0.007);
  framework->GetYaxis()->SetLabelSize(0.03);
  framework->GetYaxis()->SetTitleOffset(1.3);
  framework->GetYaxis()->SetTitleFont(42);
  framework->GetYaxis()->SetTitleSize(0.035);
  framework->Draw();

  TGraph* TGObsLimit   = new TGraph(N,MassAxis,ObsLimit);  TGObsLimit->SetLineWidth(2);
  TGraph* TGExpLimit   = new TGraph(N,MassAxis,ExpLimit);  TGExpLimit->SetLineWidth(2); TGExpLimit->SetLineStyle(2);
  TCutG* TGExpLimit1S  = GetErrorBand("1S", N, MassAxis, ExpLimitm1, ExpLimitp1);  
  TCutG* TGExpLimit2S  = GetErrorBand("2S", N, MassAxis, ExpLimitm2, ExpLimitp2);  TGExpLimit2S->SetFillColor(5);
  TGraph* THXSec        = new TGraph(N,MassAxis,ThXSec); THXSec->SetLineWidth(2); THXSec->SetLineStyle(2); THXSec->SetLineColor(4);


  TGExpLimit->SetLineColor(1);  TGExpLimit->SetLineStyle(2);
  TGObsLimit->SetLineWidth(2);  TGObsLimit->SetMarkerStyle(20);
  TGExpLimit2S->Draw("fc same");
  TGExpLimit1S->Draw("fc same");
  if(!blind) TGObsLimit->Draw("same CP");
  TGExpLimit->Draw("same C");

  if(strengthLimit){
     TLine* SMLine = new TLine(framework->GetXaxis()->GetXmin(),1.0,framework->GetXaxis()->GetXmax(),1.0);
     SMLine->SetLineWidth(2); SMLine->SetLineStyle(1); SMLine->SetLineColor(4);      
     SMLine->Draw("same C");
  }else{
//     THXSec->Draw("same C");
  }

  utils::root::DrawPreliminary(luminosity, energy, c);

  TLegend* LEG = new TLegend(0.55,0.75,0.85,0.95);
  LEG->SetHeader("");
  LEG->SetFillColor(0);
  LEG->SetFillStyle(0);
  LEG->SetTextFont(42);
  LEG->SetBorderSize(0);
//  if(!strengthLimit)LEG->AddEntry(THXSec  , "TH prediction"  ,"L");
  LEG->AddEntry(TGExpLimit  , "median expected"  ,"L");
  LEG->AddEntry(TGExpLimit1S  , "expected #pm 1#sigma"  ,"F");
  LEG->AddEntry(TGExpLimit2S  , "expected #pm 2#sigma"  ,"F");
  if(!blind) LEG->AddEntry(TGObsLimit  , "observed"  ,"LP");
  LEG->Draw();
  c->RedrawAxis();
  c->SaveAs((outputDir+"Limit.png").c_str());
  c->SaveAs((outputDir+"Limit.C").c_str());
  c->SaveAs((outputDir+"Limit.pdf").c_str()); 


  //save a summary of the limits
  FILE* pFileSum = fopen((outputDir+"LimitSummary").c_str(),"w");
  for(unsigned int i=0;i<N;i++){
    fprintf(pFileSum, "$%8.6E$ & $%8.6E$ & $[%8.6E,%8.6E]$ & $[%8.6E,%8.6E]$ & $%8.6E$ & Th=$%8.6E$ & pValue=$%8.6E$\\\\\\hline\n",MassAxis[i], ExpLimit[i], ExpLimitm1[i], ExpLimitp1[i], ExpLimitm2[i],  ExpLimitp2[i], ObsLimit[i], ThXSec[i], pValue[i]);
    if(int(MassAxis[i])%50!=0)continue; printf("%f ",ObsLimit[i]);
  }printf("\n");
  fclose(pFileSum);

  pFileSum = fopen((outputDir+"LimitRange").c_str(),"w");
  fprintf(pFileSum, "EXPECTED LIMIT --> ");                   printLimits(pFileSum,TGExpLimit, MassAxis[0], MassAxis[N-1]);
  if(!blind) fprintf(pFileSum, "OBSERVED LIMIT --> ");        printLimits(pFileSum,TGObsLimit, MassAxis[0], MassAxis[N-1]);
  fprintf(pFileSum, "Exp Limits for Model are: ");              for(unsigned int i=0;i<N;i++){if(int(MassAxis[i])%50==0) fprintf(pFileSum, "%f+-%f ",ExpLimit[i], (ExpLimitp1[i]-ExpLimitm1[i])/2.0);}fprintf(pFileSum,"\n");
  if(!blind) { fprintf(pFileSum, "Obs Limits for Model are: "); for(unsigned int i=0;i<N;i++){if(int(MassAxis[i])%50==0) fprintf(pFileSum, "%f ",ObsLimit[i]);}fprintf(pFileSum,"\n"); }
  fclose(pFileSum);
}




