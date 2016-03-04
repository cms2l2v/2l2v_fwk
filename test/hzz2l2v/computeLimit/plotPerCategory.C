
#include <string>
#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraph2D.h"
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
#include "TStyle.h"
#include "TMarker.h"
//#include "tdrstyle.C"
#include "UserCode/llvv_fwk/interface/HxswgUtils.h"
#include "UserCode/llvv_fwk/src/HxswgUtils.cc"
#include "UserCode/llvv_fwk/interface/RootUtils.h"


TGraph* getGraph(string name, int color, int width, int style, TLegend* LEG, TGraph* Ref, int type, string filePath){
//   filePath+="/LimitSummary";
   FILE* pFile = fopen(filePath.c_str(),"r");
   if(!pFile){printf("Can't open %s\n",filePath.c_str()); exit(0);}
   double mass, th, exp, obs, unused;// char buffer[1024];

   TGraph* graph = new TGraph(250);
   int N=0;
   while(fscanf(pFile,"$%le$ & $%le$ & $[%le,%le]$ & $[%le,%le]$ & $%le$ & Th=$%le$ & pValue=$%le$\\\\\\hline\n",&mass, &exp, &unused, &unused, &unused,&unused,&obs, &th, &unused) != EOF){
//      printf("%i %f - %f - %f\n",N,mass,exp, th);

      double value = exp;
      if(abs(type)==0) value = th;
      else if(abs(type)==1) value = exp;
      else value = obs;

//      if(type<0 && filePath.find("cp0.80")!=string::npos) value *= pow(0.8, 2);
//      if(type<0 && filePath.find("cp0.60")!=string::npos) value *= pow(0.6, 2);
//      if(type<0 && filePath.find("cp0.30")!=string::npos) value *= pow(0.3, 2);
//      if(type<0 && filePath.find("cp0.10")!=string::npos) value *= pow(0.1, 2);

      if(Ref){value/=Ref->Eval(mass);}
      graph->SetPoint(N, mass, value);N++;
      if(N>100)break;
   }
   graph->Set(N);

   graph->SetName(name.c_str());
   graph->SetLineColor(color);
   graph->SetLineWidth(width);
   graph->SetLineStyle(style);
   if(LEG)LEG->AddEntry(graph, name.c_str()      ,"L");
   return graph;
}

void scaleGraph(TGraph* Limit, double scale){
   for(int i=0;i<Limit->GetN();i++){
      Limit->SetPoint(i, Limit->GetX()[i], Limit->GetY()[i]*scale);
   }   
}

void plotPerCategory(){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //   setTDRStyle();
   gStyle->SetPadTopMargin   (0.04);
   gStyle->SetPadBottomMargin(0.12);
   gStyle->SetPadRightMargin (0.05);
   gStyle->SetPadLeftMargin  (0.12);
   gStyle->SetTitleSize(0.04, "XYZ");
   gStyle->SetTitleXOffset(1.1);
   gStyle->SetTitleYOffset(1.45);
   gStyle->SetPalette(1);
   gStyle->SetNdivisions(505);

   TCanvas* c1;
   TH1F* framework;
   TH2F* framework2d;
   TLegend* LEG, *LEGTH;
   TGraph* Ref;

   string Directories[]={"cards_SB13TeV_cp1.00_brn0.00", "cards_SB13TeV_GGF_cp1.00_brn0.00", "cards_SB13TeV_VBF_cp1.00_brn0.00"};
   for(unsigned int D=0;D<sizeof(Directories)/sizeof(string);D++){
      string Dir = Directories[D];

     string prod = "pp";
     if(Dir.find("GGF")!=std::string::npos)prod="gg";
     if(Dir.find("VBF")!=std::string::npos)prod="qq";
     bool strengthLimit = false;
     if(prod=="pp")strengthLimit=true;

      c1 = new TCanvas("c", "c",600,600);
      c1->SetLogy(true);
      framework = new TH1F("Graph","Graph",1,190,1510);
      framework->SetStats(false);
      framework->SetTitle("");
      framework->GetXaxis()->SetTitle("Higgs boson mass [GeV]");
      if(strengthLimit){
         framework->GetYaxis()->SetTitle("#mu = #sigma_{95%} / #sigma_{th}");
         framework->GetYaxis()->SetRangeUser(1E-2,1E3);
      }else{
         framework->GetYaxis()->SetTitle((string("#sigma_{95%} (") + prod + " #rightarrow H #rightarrow ZZ) (fb)").c_str());        
         framework->GetYaxis()->SetRangeUser(1E1,1E5);
      }
      framework->GetYaxis()->SetTitleOffset(1.40);
      framework->Draw();

      LEG = new TLegend(0.70,0.70,0.95,0.94);
      LEG->SetFillStyle(0);
      LEG->SetBorderSize(0);
      LEG->SetHeader("Observed:");

      TLegend* LEGExp = NULL;
      LEGExp = new TLegend(0.45,0.70,0.70,0.94);
      LEGExp->SetFillStyle(0);
      LEGExp->SetBorderSize(0);
      LEGExp->SetHeader("Expected:");


      getGraph("=0 Jet"                       , 2, 2, 1, LEG  , NULL, 2, Dir+ "_eq0jets/Stength_LimitSummary")->Draw("C same");
      getGraph("#geq1 Jets"                   , 4, 2, 1, LEG  , NULL, 2, Dir+"_geq1jets/Stength_LimitSummary")->Draw("C same");
      getGraph("VBF"                          , 6, 2, 1, LEG  , NULL, 2, Dir+     "_vbf/Stength_LimitSummary")->Draw("C same");
      getGraph("Combined"                     , 1, 2, 1, LEG  , NULL, 2, Dir+         "/Stength_LimitSummary")->Draw("C same");

      getGraph("=0 Jet"                       , 2, 2, 2, LEGExp  , NULL, 1, Dir+ "_eq0jets/Stength_LimitSummary")->Draw("C same");
      getGraph("#geq1 Jets"                   , 4, 2, 2, LEGExp  , NULL, 1, Dir+"_geq1jets/Stength_LimitSummary")->Draw("C same");
      getGraph("VBF"                          , 6, 2, 2, LEGExp  , NULL, 1, Dir+     "_vbf/Stength_LimitSummary")->Draw("C same");
      getGraph("Combined"                     , 1, 2, 2, LEGExp  , NULL, 1, Dir+         "/Stength_LimitSummary")->Draw("C same");


   //   LEGTH->Draw("same");
      LEGExp  ->Draw("same");
      LEG  ->Draw("same");


      /*
      char LumiLabel[1024];
      sprintf(LumiLabel,"CMS preliminary,  #sqrt{s}=%.0f TeV #scale[0.5]{#int} L=%6.1ffb^{-1}",13.0,2.3);
      TPaveText *pave = new TPaveText(0.1,0.96,0.94,0.99,"NDC");
      pave->SetBorderSize(0);
      pave->SetFillStyle(0);
      pave->SetTextAlign(32);
      pave->SetTextFont(42);
      pave->AddText(LumiLabel);
      pave->Draw("same");*/

      utils::root::DrawPreliminary(2268.759, 13, gPad->GetLeftMargin(),gPad->GetBottomMargin(),gPad->GetRightMargin(),gPad->GetTopMargin()+0.025);     
      if(strengthLimit){
         TLine* SMLine = new TLine(framework->GetXaxis()->GetXmin(),1.0,framework->GetXaxis()->GetXmax(),1.0);
         SMLine->SetLineWidth(2); SMLine->SetLineStyle(1); SMLine->SetLineColor(4);      
         SMLine->Draw("same C");
      }else{
         TGraph* THXSec   = Hxswg::utils::getXSec(Dir); 
         THXSec->SetLineWidth(2); THXSec->SetLineStyle(1); THXSec->SetLineColor(4);
         scaleGraph(THXSec, 1000);  //convert cross-section to fb
         THXSec->Draw("same C");
      }

      c1->SaveAs((Dir+"/perCat_FinalPlot.png").c_str());
      c1->SaveAs((Dir+"/perCat_FinalPlot.pdf").c_str());
      c1->SaveAs((Dir+"/perCat_FinalPlot.C"  ).c_str());
   }



/*
   for(unsigned int D=0;D<sizeof(Directories)/sizeof(string);D++){
      string Dir = Directories[D];

      c1 = new TCanvas("c", "c",600,600);
      c1->SetLogy(true);
      framework = new TH1F("Graph","Graph",1,150,1050);
      framework->SetStats(false);
      framework->SetTitle("");
      framework->GetXaxis()->SetTitle("Higgs boson mass [GeV]");
//      framework->GetYaxis()->SetTitle("#mu = #sigma_{95%} / #sigma_{th}");
//      framework->GetYaxis()->SetTitle("#sigma_{95%} (fb)");
      framework->GetYaxis()->SetTitle("#sigma_{95%} (pp #rightarrow H #rightarrow ZZ) (fb)");        
      framework->GetYaxis()->SetTitleOffset(1.40);
      framework->GetYaxis()->SetRangeUser(1E1,1E4);
      framework->Draw();

      LEG = new TLegend(0.70,0.70,0.95,0.94);
      LEG->SetFillStyle(0);
      LEG->SetBorderSize(0);
//      LEG->SetHeader("Expected @95% C.L.");

      LEGTH = new TLegend(0.45,0.70,0.70,0.94);
      LEGTH->SetFillStyle(0);
      LEGTH->SetBorderSize(0);
      LEGTH->SetHeader("Theoretical");

//      getGraph("SM-like"                     , 1, 2, 1, LEG  , NULL, 1, Dir+               "/Stength_LimitSummary")->Draw("C same");
      getGraph("C'=1.0"                      , 2, 2, 2, LEG  , NULL, 1, Dir+"_cp1.00_brn0.00/Stength_LimitSummary")->Draw("C same");
      getGraph("C'=0.8"                      , 4, 2, 2, LEG  , NULL, 1, Dir+"_cp0.80_brn0.00/Stength_LimitSummary")->Draw("C same");
      getGraph("C'=0.6"                      , 6, 2, 2, LEG  , NULL, 1, Dir+"_cp0.60_brn0.00/Stength_LimitSummary")->Draw("C same");
      getGraph("C'=0.4"                      , 7, 2, 2, LEG  , NULL, 1, Dir+"_cp0.40_brn0.00/Stength_LimitSummary")->Draw("C same");
      getGraph("C'=0.2"                      , 8, 2, 2, LEG  , NULL, 1, Dir+"_cp0.20_brn0.00/Stength_LimitSummary")->Draw("C same");

   //   LEGTH->Draw("same");
      LEG  ->Draw("same");

      char LumiLabel[1024];
      sprintf(LumiLabel,"CMS preliminary,  #sqrt{s}=%.0f TeV #scale[0.5]{#int} L=%6.1ffb^{-1}",13.0,2.2);
      TPaveText *pave = new TPaveText(0.1,0.96,0.94,0.99,"NDC");
      pave->SetBorderSize(0);
      pave->SetFillStyle(0);
      pave->SetTextAlign(32);
      pave->SetTextFont(42);
      pave->AddText(LumiLabel);
      pave->Draw("same");

      TLine* SMLine = new TLine(framework->GetXaxis()->GetXmin(),1.0,framework->GetXaxis()->GetXmax(),1.0);
      SMLine->SetLineWidth(2); SMLine->SetLineStyle(1); SMLine->SetLineColor(4);      
//      SMLine->Draw("same C");


      c1->SaveAs((Dir+"/perC_FinalPlot.png").c_str());
      c1->SaveAs((Dir+"/perC_FinalPlot.pdf").c_str());
      c1->SaveAs((Dir+"/perC_FinalPlot.C"  ).c_str());
   }
*/
}






