
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
#include "UserCode/llvv_fwk/interface/HxswgUtils.h"
#include "UserCode/llvv_fwk/src/HxswgUtils.cc"

//tree variables
double Tmh, Tlimit, TlimitErr; float TquantExp;


TGraph* getLimitGraph(TTree* tree, float Quantil){
  TGraph* toReturn = new TGraph(999);
  int i=0; 
  for(int ientry=0;ientry<tree->GetEntriesFast();ientry++){
     tree->GetEntry(ientry);
     
     if(TquantExp==Quantil){
        //printf("Quantil = %f - mH=%f --> %f\n",TquantExp,Tmh,Tlimit);
        //
  //      if(Tmh>1000)continue;
        toReturn->SetPoint(i, Tmh, Tlimit);
        i++;
     }
  }toReturn->Set(i);
  return toReturn;

}


TCutG* GetErrorBand(string name, TGraph* Low, TGraph* High)
{
   TCutG* cutg = new TCutG(name.c_str(),Low->GetN()+High->GetN()+2);
   cutg->SetFillColor(kGreen-7);
   cutg->SetLineStyle(0);
   cutg->SetLineColor(0);
   int I = 0;
   for(int i=0;i<Low->GetN();i++){  cutg->SetPoint(I,Low ->GetX()[i]               , Low ->GetY()[i]               );I++; }
                                    cutg->SetPoint(I,Low ->GetX()[Low ->GetN()-1]  , Low ->GetY()[Low ->GetN()-1]  );I++;
                                    cutg->SetPoint(I,High->GetX()[High->GetN()-1]  , High->GetY()[High->GetN()-1]  );I++;
   for(int i=0;i<High->GetN() ;i++){cutg->SetPoint(I,High->GetX()[High->GetN()-1-i], High->GetY()[High->GetN()-1-i]);I++;}
   return cutg;
}

void scaleGraph(TGraph* Limit, double scale){ 
   for(int i=0;i<Limit->GetN();i++){ 
      Limit->SetPoint(i, Limit->GetX()[i], Limit->GetY()[i]*scale);
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
  TGraph* ObsLimit   = getLimitGraph(tree,-1   ); 
  TGraph* ExpLimitm2 = getLimitGraph(tree,0.025);
  TGraph* ExpLimitm1 = getLimitGraph(tree,0.160);
  TGraph* ExpLimit   = getLimitGraph(tree,0.500);
  TGraph* ExpLimitp1 = getLimitGraph(tree,0.840);
  TGraph* ExpLimitp2 = getLimitGraph(tree,0.975);
  file->Close(); 

  FILE* pFileSStrenght = fopen((outputDir+"SignalStrenght").c_str(),"w");
  std::cout << "Printing Signal Strenght" << std::endl;
  for(int i=0;i<ExpLimit->GetN();i++){
     double M = ExpLimit->GetX()[i];
     std::cout << "Mass: " << M << "; ExpLimit: " << ExpLimit->Eval(M) << std::endl; 
     printf("$%8.6E$ & $%8.6E$ & $[%8.6E,%8.6E]$ & $[%8.6E,%8.6E]$ \\\\\\hline\n",M, ExpLimit->Eval(M), ExpLimitm1->Eval(M), ExpLimitp1->Eval(M), ExpLimitm2->Eval(M),  ExpLimitp2->Eval(M));
     fprintf(pFileSStrenght, "$%8.6E$ & $%8.6E$ & $[%8.6E,%8.6E]$ & $[%8.6E,%8.6E]$ & $%8.6E$ \\\\\\hline\n",M, ExpLimit->Eval(M), ExpLimitm1->Eval(M), ExpLimitp1->Eval(M), ExpLimitm2->Eval(M),  ExpLimitp2->Eval(M), ObsLimit->Eval(M));
    if(int(ExpLimit->GetX()[i])%50!=0)continue; //printf("%f ",ObsLimit->Eval(M));
  }printf("\n");
  fclose(pFileSStrenght); 
 

  //get the pValue
  inputs = inputs.ReplaceAll("/LimitTree.root", "/PValueTree.root");
  file = TFile::Open(inputs);
  
  printf("Looping on %s\n",inputs.Data());
  if(!file) return;
  if(file->IsZombie()) return;
  
  tree = (TTree*)file->Get("limit");
  
  tree->GetBranch("limit"           )->SetAddress(&Tlimit   );
  
  TGraph* pValue     = getLimitGraph(tree,-1);
  
  file->Close();

  
  //make TH Cross-sections
   string suffix = outputDir;
   TGraph* THXSec   = Hxswg::utils::getXSec(outputDir); 
   scaleGraph(THXSec, 1000);  //convert cross-section to fb
   double cprime=1.0; double  brnew=0.0;
   double XSecScaleFactor = 1.0;
   if(suffix.find("_cp")!=string::npos){
     sscanf(suffix.c_str()+suffix.find("_cp"), "_cp%lf_brn%lf", &cprime, &brnew);
     XSecScaleFactor = pow(cprime,2) * (1-brnew);
   }
  //XSecScaleFactor = 0.001; //pb to fb
  scaleGraph(THXSec, XSecScaleFactor);


  string prod = "pp";
  if(outputDir.find("GGF")!=std::string::npos)prod="gg";
  if(outputDir.find("VBF")!=std::string::npos)prod="qq";

  
  strengthLimit = false;
  if(prod=="pp")strengthLimit=true;
 
  //TGraph *XSecMELA = Hxswg::utils::getXSecMELA(cprime);

  //Hxswg::utils::multiplyGraph(   ObsLimit, XSecMELA);
  //Hxswg::utils::multiplyGraph( ExpLimitm2, XSecMELA);
  //Hxswg::utils::multiplyGraph( ExpLimitm1, XSecMELA);
  //Hxswg::utils::multiplyGraph(   ExpLimit, XSecMELA);
  //Hxswg::utils::multiplyGraph( ExpLimitp1, XSecMELA);
  //Hxswg::utils::multiplyGraph( ExpLimitp2, XSecMELA);
 
  //Scale exclusion XSec in fb
  scaleGraph(ObsLimit  , 0.001); //pb to fb
  scaleGraph(ExpLimitm2, 0.001); //pb to fb
  scaleGraph(ExpLimitm1, 0.001); //pb to fb
  scaleGraph(ExpLimit  , 0.001); //pb to fb
  scaleGraph(ExpLimitp1, 0.001); //pb to fb
  scaleGraph(ExpLimitp2, 0.001); //pb to fb

  //scal eTH cross-section and limits according to scale factor 
  //this only apply to NarrowResonnance case
  if(strengthLimit){
     Hxswg::utils::divideGraph(ObsLimit   , THXSec);
     Hxswg::utils::divideGraph(ExpLimitm2 , THXSec);
     Hxswg::utils::divideGraph(ExpLimitm1 , THXSec);
     Hxswg::utils::divideGraph(ExpLimit   , THXSec);
     Hxswg::utils::divideGraph(ExpLimitp1 , THXSec);
     Hxswg::utils::divideGraph(ExpLimitp2 , THXSec);
     Hxswg::utils::divideGraph(THXSec     , THXSec);
  }


  //limits in terms of signal strength
  TCanvas* c = new TCanvas("c", "c",800,800);
  c->SetGridx();
  c->SetGridy();
  TH1F* framework = new TH1F("Graph","Graph",1,strengthLimit?199:199,2500); //3000);
  framework->SetStats(false);
  framework->SetTitle("");
  framework->GetXaxis()->SetTitle("M_{H} [GeV]");
  framework->GetYaxis()->SetTitleOffset(1.70);
  if(strengthLimit){
  framework->GetYaxis()->SetTitle("#mu = #sigma_{95%} / #sigma_{th}");
  framework->GetYaxis()->SetRangeUser(1E-4,1E3);
  c->SetLogy(true);
  }else{
  framework->GetYaxis()->SetTitle((string("#sigma_{95%} (") + prod +" #rightarrow H #rightarrow ZZ) (pb)").c_str());
  framework->GetYaxis()->SetRangeUser(1E-3,1E3);
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

  
  TGraph* TGObsLimit   = ObsLimit;  TGObsLimit->SetLineWidth(2);
  TGraph* TGExpLimit   = ExpLimit;  TGExpLimit->SetLineWidth(2); TGExpLimit->SetLineStyle(2);
  TCutG* TGExpLimit1S  = GetErrorBand("1S", ExpLimitm1, ExpLimitp1);  
  TCutG* TGExpLimit2S  = GetErrorBand("2S", ExpLimitm2, ExpLimitp2);  TGExpLimit2S->SetFillColor(5);
  THXSec->SetLineWidth(2); THXSec->SetLineStyle(1); THXSec->SetLineColor(4);

  TGExpLimit->SetLineColor(1);  TGExpLimit->SetLineStyle(2);
  TGObsLimit->SetLineWidth(2);  TGObsLimit->SetMarkerStyle(20);
  TGExpLimit2S->Draw("fc same");
  TGExpLimit1S->Draw("fc same");
  if(!blind) TGObsLimit->Draw("same P");
  TGExpLimit->Draw("same c");

  
  /*if(strengthLimit){
     TLine* SMLine = new TLine(framework->GetXaxis()->GetXmin(),1.0,framework->GetXaxis()->GetXmax(),1.0);
     SMLine->SetLineWidth(2); SMLine->SetLineStyle(1); SMLine->SetLineColor(4);      
     SMLine->Draw("same C");
  }else{
     THXSec->Draw("same C");
  }*/

  utils::root::DrawPreliminary(luminosity, energy, c);

  
  TLegend* LEG = new TLegend(0.55,0.75,0.85,0.95);
  LEG->SetHeader("");
  LEG->SetFillColor(0);
  LEG->SetFillStyle(0);
  LEG->SetTextFont(42);
  LEG->SetBorderSize(0);
  //LEG->AddEntry(THXSec  , "Th prediction"  ,"L");
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
  for(int i=0;i<TGExpLimit->GetN();i++){
     double M = ExpLimit->GetX()[i];
     fprintf(pFileSum, "$%8.6E$ & $%8.6E$ & $[%8.6E,%8.6E]$ & $[%8.6E,%8.6E]$ & $%8.6E$ & Th=$%8.6E$ & pValue=$%8.6E$\\\\\\hline\n",M, ExpLimit->Eval(M), ExpLimitm1->Eval(M), ExpLimitp1->Eval(M), ExpLimitm2->Eval(M),  ExpLimitp2->Eval(M), ObsLimit->Eval(M), (THXSec!=NULL)?THXSec->Eval(M):-1, pValue->Eval(M));
    if(int(ExpLimit->GetX()[i])%50!=0)continue; printf("%f ",ObsLimit->Eval(M));
  }printf("\n");
  fclose(pFileSum);

  pFileSum = fopen((outputDir+"LimitRange").c_str(),"w");
  fprintf(pFileSum, "EXPECTED LIMIT --> ");                   printLimits(pFileSum,TGExpLimit, TGExpLimit->GetX()[0], TGExpLimit->GetX()[TGExpLimit->GetN()-1]);
  if(!blind) fprintf(pFileSum, "OBSERVED LIMIT --> ");        printLimits(pFileSum,TGObsLimit, TGObsLimit->GetX()[0], TGObsLimit->GetX()[TGObsLimit->GetN()-1]);
  fprintf(pFileSum, "Exp Limits for Model are: ");              for(int i=0;i<TGExpLimit->GetN();i++){if(int(TGExpLimit->GetX()[i])%50==0) fprintf(pFileSum, "%f+-%f ",TGExpLimit->GetY()[i], (ExpLimitp1->GetY()[i]-ExpLimitm1->GetY()[i])/2.0);}fprintf(pFileSum,"\n");
  if(!blind) { fprintf(pFileSum, "Obs Limits for Model are: "); for(int i=0;i<TGObsLimit->GetN();i++){if(int(TGObsLimit->GetX()[i])%50==0) fprintf(pFileSum, "%f ",TGObsLimit->GetY()[i]);}fprintf(pFileSum,"\n"); }
  fclose(pFileSum); 
}




