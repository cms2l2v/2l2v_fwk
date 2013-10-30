#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TString.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLine.h"

TH1F *generateBtagsFor(float R=1.0, TString url="~/work/top_5311/plotter_syst.root");
void generateRpdfs(TString url="~/work/top_5311/plotter_syst.root");

//
TH1F *generateBtagsFor(float R, TString url)
{

  TFile *inF=TFile::Open(url);
  TH1F *r1=(TH1F*)inF->Get("t#bar{t}systr1/csvLbtagsextended")->Clone("r1");       r1->SetDirectory(0);   r1->Scale(R*R/r1->Integral());
  TH1F *r05a=(TH1F*)inF->Get("t#bar{t}systr05a/csvLbtagsextended")->Clone("r05a"); r05a->SetDirectory(0); r05a->Scale(R*(1-R)/r05a->Integral());
  TH1F *r05b=(TH1F*)inF->Get("t#bar{t}systr05b/csvLbtagsextended")->Clone("r05b"); r05b->SetDirectory(0); r05b->Scale(R*(1-R)/r05b->Integral());
  TH1F *r05=(TH1F *)r05a->Clone("r05"); r05->Add(r05b); r05->SetDirectory(0);
  TH1F *r0=(TH1F*)inF->Get("t#bar{t}systr0/csvLbtagsextended")->Clone("r0");       r0->SetDirectory(0);   r0->Scale((1-R)*(1-R)/r0->Integral());
  inF->Close();

  r1->SetFillStyle(1001); r1->SetFillColor(614);
  r05->SetFillStyle(1001); r05->SetFillColor(824);
  r0->SetFillStyle(1001); r0->SetFillColor(592);


  THStack *total=new THStack;
  total->Add(r0,"hist");
  total->Add(r05,"hist");
  total->Add(r1,"hist");

  TCanvas *c=new TCanvas("c","c",600,600);
  total->Draw();
  total->GetXaxis()->SetTitle(r1->GetXaxis()->GetTitle());
  total->GetYaxis()->SetTitle("Events (a.u.)");
  total->GetYaxis()->SetRangeUser(0,0.2);


  for(int i=0; i<3; i++)
    {
      TString ch("ee");
      if(i==1) ch="#mu#mu";
      if(i==2) ch="e#mu";
      
      TPaveText *pt=new TPaveText(3*5*i+7,0.12,3*5*i+12,0.15);
      pt->SetFillStyle(0);
      pt->SetTextAlign(12);
      pt->SetBorderSize(0);
      pt->SetTextFont(52);
      pt->AddText(ch);
      pt->SetTextSize(0.025);
      pt->Draw();
      
      TLine *l=new TLine(3*5*(1+i),0,3*5*(1+i),0.15);
      l->SetLineStyle(7);
      l->Draw();
    }



  TPaveText *pt=new TPaveText(0.1,0.95,0.6,0.99,"brNDC");
  pt->SetFillStyle(0);
  pt->SetTextAlign(12);
  pt->SetBorderSize(0);
  pt->SetTextSize(0.03);
  char buf[50];
  sprintf(buf,"|V_{tb}|=%3.2f",sqrt(R));
  pt->AddText("CMS simulation, #sqrt{s}=8 TeV, " + TString(buf));
  pt->Draw();

  TLegend *leg=new TLegend(0.65,0.95,0.96,0.99);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->AddEntry(r0,"R=0","f");
  leg->AddEntry(r05,"R=0.5","f");
  leg->AddEntry(r1,"R=1.0","f");
  leg->SetNColumns(3);
  leg->Draw();

  TString name("csvLbtagsextended_r");
  name += (int)(R*100);

  c->Modified();
  c->Update();
  c->SaveAs(name+".pdf");



  TH1F *h=(TH1F *)total->GetStack()->At( total->GetStack()->GetEntriesFast()-1 )->Clone(name);
  h->SetDirectory(0);
  return h;
}

//
void generateRpdfs(TString url)
{
  TObjArray pdfs;
  for(float r=0; r<=1.05; r+=0.05) pdfs.Add( generateBtagsFor(r,url) );

  TFile *outF=TFile::Open("RclosurePDFs.root","RECREATE");
  for(size_t i=0; i<pdfs.GetEntriesFast(); i++) pdfs.At(i)->Write();
  outF->Close();
}


