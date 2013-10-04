#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLegend.h"

void compareLineShapes(){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  
  //TFile *fIn=TFile::Open("/afs/cern.ch/user/p/psilva/work/CMSSW_5_3_11/src/UserCode/llvv_fwk/data/weights/GGtoHtoZZLineShapes.root");
  //TString proc("gg#rightarrow H #rightarrow ZZ");

  TFile *fIn=TFile::Open("/afs/cern.ch/user/p/psilva/work/CMSSW_5_3_11/src/UserCode/llvv_fwk/data/weights/VBFtoHtoZZLineShapes.root");
  TString proc("qq#rightarrow qqH #rightarrow qqZZ");
  
  TString masses[]={"H200","H300","H400","H500","H600","H700","H800","H900","H1000"};
  Int_t nmasses=sizeof(masses)/sizeof(TString);
  for(Int_t i=0; i<nmasses; i++)
    {
      
      TCanvas *c=new TCanvas("c","c",600,600);
      TPad *p=new TPad("t1","t1",0,0.4,1.0,1.0);
      p->Draw();
      p->cd();
      TH1 *gen=(TH1 *)fIn->Get(masses[i]+"/gen");               gen->SetLineColor(1);     gen->SetLineWidth(1);     gen->SetTitle("Powheg");      gen->SetDirectory(0);
      TH1 *cps=(TH1 *)fIn->Get(masses[i]+"/cps_shape");         cps->SetLineColor(4);                               cps->SetTitle("CPS");         cps->SetDirectory(0);
      TH1 *cpspint=(TH1 *)fIn->Get(masses[i]+"/nominal_shape"); cpspint->SetLineColor(2); cpspint->SetLineWidth(2); cpspint->SetTitle("CPS+Int"); cpspint->SetDirectory(0);
      gen->Draw("hist");
      cps->Draw("histsame");
      cpspint->Draw("histsame");
      gen->GetYaxis()->SetTitle("Events (a.u.)");
      gen->GetXaxis()->SetTitle("Higgs mass [GeV]");
      gen->GetYaxis()->SetRangeUser(1e-4,1.0);
      gen->GetXaxis()->SetRangeUser(0,1500);
      p->SetLogy();
      
      TLegend *leg=new TLegend(0.15,0.8,0.9,0.9);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->AddEntry(gen,gen->GetTitle(),"f");
      leg->AddEntry(cps,cps->GetTitle(),"f");
      leg->AddEntry(cpspint,cpspint->GetTitle(),"f");
      leg->SetTextSize(0.05);
      leg->SetNColumns(3);
      leg->Draw();

      TPaveText *pt=new TPaveText(0.1,0.95,0.9,0.99,"brNDC");
      pt->SetBorderSize(0);
      pt->SetFillStyle(0);
      pt->SetTextAlign(12);
      pt->AddText("CMS simulation, #sqrt{s}=8 TeV, " + proc);
      pt->Draw();

      c->cd();
      p=new TPad("t2","t2",0.0,0.2,1.0,0.4);
      p->Draw();
      p->cd();
      TGraph *gen2cps=(TGraph *)fIn->Get(masses[i]+"/cps");
      gen2cps->Draw("al");
      gen2cps->GetYaxis()->SetTitle("CPS/Powheg");
      gen2cps->GetXaxis()->SetRangeUser(0,1500);
      gen2cps->GetYaxis()->SetRangeUser(0,2);
      gen2cps->GetYaxis()->SetTitleOffset(0.5);
      gen2cps->GetYaxis()->SetTitleSize(0.09);
      gen2cps->GetYaxis()->SetLabelSize(0.08);
      gen2cps->GetXaxis()->SetLabelSize(0.08);


      c->cd();
      p=new TPad("t3","t3",0.0,0.0,1.0,0.2);
      p->Draw();
      p->cd();
      TGraph *cps2nominal=(TGraph *)fIn->Get(masses[i]+"/nominal");
      cps2nominal->Draw("al");
      cps2nominal->GetYaxis()->SetTitle("CPS+Int/CPS");
      cps2nominal->GetXaxis()->SetRangeUser(0,1500);      
      cps2nominal->GetYaxis()->SetRangeUser(0,2);
      cps2nominal->GetYaxis()->SetTitleOffset(0.5);
      cps2nominal->GetYaxis()->SetTitleSize(0.09);
      cps2nominal->GetYaxis()->SetLabelSize(0.08);
      cps2nominal->GetXaxis()->SetLabelSize(0.08);

      c->cd();
      c->Modified();
      c->Update();
      c->SaveAs(masses[i]+"_lineshape.png");
    }



}
