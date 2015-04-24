// 
// Author: Xin Shi <Xin.Shi@cern.ch> 
// Created: 2015.04.21
// 

#include <iostream>

#include <TApplication.h> 
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TStyle.h> 
#include <TROOT.h> 
#include <TLegend.h> 

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


TCanvas* drawPt(TString label, std::vector<TString> inputFiles){
  if (label == "pho2") printf("label = %s \n", label.Data());
  int ww(800), wh(800);
  set_root_style();
  TCanvas *c = new TCanvas("c", "Transverse momentum", ww, wh);

  TH1F *h = NULL;
  TH1F *h2 = NULL; 

  TLegend *leg = new TLegend(0.2, 0.6, 0.5, 0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetNColumns(1);
  leg->SetTextSize(0.02);
  leg->SetTextSizePixels(25);


  for (std::vector<int>:: size_type i = 0; i != inputFiles.size(); i++) {    
    TFile *f = new TFile(inputFiles[i]);
    TString dir = "#gamma+jets_pT-15to3000"; 
    TString histName = Form("%s/qt", dir.Data()); 
    h = (TH1F*)f->Get(histName);
    if (!h) {
      std::cout << "Not able to find hist: " << histName << std::endl; 
      return NULL; 
    }
    
    if (label == "pho2") {
      TString histName2 = Form("%s/trg_phopt", dir.Data()); 
      h2 = (TH1F*)f->Get(histName2);
      if (!h2) {
	std::cout << "Not able to find hist: " << histName2 << std::endl; 
	return NULL; 
      }

      leg->AddEntry(h2, "triggerd");
      h2->SetMarkerColor(2);
      h2->SetTitle(""); 
      h2->Draw();
    }

    h->GetXaxis()->SetRangeUser(10,1000);
    if (label == "pho2") {
      leg->AddEntry(h, "selected"); 
      h->Draw("same");
      leg->Draw(); 
    }
    else {h->Draw();}
  }

  c->SetLogx();
  c->SetLogy();
  c->Update(); 
  return c; 
}


void print_usage(){
  printf("NAME\n\tdrawPt - draw transverse momentum spectrum\n");
  printf("\nSYNOPSIS\n\tdrawPt label input.root\n "); 
  // printf("\nSYNOPSIS\n\tdrawPt [-t hist-type ] [-opt draw-option]\n "); 
  // printf("\t[-h hist-name ] [-vmax max-value] [-npad num-pad] [-b] input1 input2 ...\n");
  printf("\nOPTIONS\n");
  printf("\t%-5s  %-40s\n", "label", "options: pho, pho2");
  // printf("\t%-5s  %-40s\n", "-t", "hist type [TH1D, TH2D]");
  // printf("\n\t%-5s  %-40s\n", "-opt", "draw option for histgram [colz, surf2]");
  // printf("\n\t%-5s  %-40s\n", "-h", "histogram name");
  // printf("\n\t%-5s  %-40s\n", "-vmax", "limit the histogram within vmax-value");
  // printf("\n\t%-5s  %-40s\n", "-vmin", "limit the histogram within vmin-value");
  // printf("\n\t%-5s  %-40s\n", "-p", "print pixel values to stdout");
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
  drawPt(label, inputFiles);
  theApp.Run();
}

