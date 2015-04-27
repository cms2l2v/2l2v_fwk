// 
// Author: Xin Shi <Xin.Shi@cern.ch> 
// Created: 2015.04.27
//
// Draw efficiencies 

#include <iostream>

#include <TApplication.h> 
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TStyle.h> 
#include <TROOT.h> 
#include <TLegend.h> 
#include <TEfficiency.h> 

void set_root_style(int stat=1110, int grid=0){
  gROOT->Reset();

  gStyle->SetTitleFillColor(0) ; 
  gStyle->SetTitleBorderSize(0); 
    
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasDefX(0); 
  gStyle->SetCanvasDefY(0); 
  gStyle->SetFrameBorderMode(0); 
  gStyle->SetFrameBorderSize(1); 
  gStyle->SetFrameFillColor(0); 
  gStyle->SetFrameFillStyle(0); 
  gStyle->SetFrameLineColor(1); 
  gStyle->SetFrameLineStyle(1); 
  gStyle->SetFrameLineWidth(1); 

  // gStyle->SetPadTopMargin(PadTopMargin);  
  gStyle->SetPadLeftMargin(0.15);  
  gStyle->SetPadRightMargin(0.05);  

  gStyle->SetLabelSize(0.03, "XYZ");  
  gStyle->SetTitleSize(0.04, "XYZ");  
  gStyle->SetTitleOffset(1.2, "Y");  

  gStyle->SetPadBorderMode(0);  
  gStyle->SetPadColor(0);  
  gStyle->SetPadTickX(1); 
  gStyle->SetPadTickY(1); 
  gStyle->SetPadGridX(grid); 
  gStyle->SetPadGridY(grid); 

  gStyle->SetOptStat(stat); 
  gStyle->SetStatColor(0); 
  gStyle->SetStatBorderSize(1); 
}

TH1F * get_hist(TString dir, TString hname,
		std::vector<TString> inputFiles) {
 
  TH1F *hall = NULL;
  for (std::vector<int>:: size_type i = 0; i != inputFiles.size(); i++) {    
    TFile *f = new TFile(inputFiles[i]);
    TString histName = Form("%s/%s", dir.Data(), hname.Data()); 
    TH1F *h = (TH1F*)f->Get(histName);
    if (!h) {
      std::cout << "Not able to find hist: " << histName << std::endl; 
      return NULL; 
    }
    if (!hall) hall = h;
    else hall->Add(h); 
  }
  return hall; 
}

TCanvas* drawEff(TString label, std::vector<TString> inputFiles){
  int ww(800), wh(800);
  set_root_style();
  TCanvas *c = new TCanvas("c", "Efficiency", ww, wh);

  TString dir = "#gamma+jets_pT-15to3000";

  TH1F *h_total = get_hist(dir, "slim_phopt", inputFiles);
  if (!h_total) return NULL;

  TH1F *h_pass = get_hist(dir, "trg_phopt", inputFiles);
  if (!h_pass) return NULL;

  h_total->Draw(); 
  h_pass->Draw("same");

  int xlow = 0;
  int xup = 1000; 
  int nbins = h_total->GetSize()-2;
  // printf("Total bins = %d", nbins); 

  TEfficiency* eff = new TEfficiency("eff", ";p_{T} [GeV];Efficiency", nbins, xlow, xup);

  for (int bin=1; bin<=nbins; bin++) {
    int iTotal = h_total->GetBinContent(bin); 
    int iPass = h_pass->GetBinContent(bin); 

    if (iPass  == 0 ) continue; 

    eff->SetTotalEvents(bin, iTotal); 
    eff->SetPassedEvents(bin, iPass); 
  }

  eff->Draw();
  // eff->GetYaxis()->SetRangeUser(0., 1.0);

  // c->SetLogx();
  // c->SetLogy();
  c->Update(); 
  return c; 
}


void print_usage(){
  printf("NAME\n\tdrawEff - draw efficiencies\n");
  printf("\nSYNOPSIS\n\tdrawPt LABEL input.root\n "); 
  printf("\nLABEL\n");
  printf("\t%-5s  %-40s\n", "photrg", "efficiencies after photon triggers");
  printf("\nAUTHOR\n\tXin Shi <Xin.Shi@cern.ch>\n");
}


int main(int argc, char** argv) {
  if (argc < 2) {
    print_usage() ;  
    return -1; 
  }

  TString label(argv[1]); 
  std::vector<TString> inputFiles(argv+2, argv+argc);

  TApplication theApp("App", 0, 0);
  theApp.SetReturnFromRun(true);
  drawEff(label, inputFiles);
  theApp.Run();
}


