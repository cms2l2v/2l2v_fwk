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

TH1F * get_hist(TString dir, std::vector<TString> hnames,
		std::vector<TString> inputFiles) {
 
  TH1F *hall = NULL;
  for (std::vector<int>:: size_type i = 0; i != inputFiles.size(); i++) {    
    TFile *f = new TFile(inputFiles[i]);
    for (std::vector<TString>:: size_type j=0; j!=hnames.size(); j++) {    
      TString hname = hnames[j]; 
      TString histName = Form("%s/%s", dir.Data(), hname.Data()); 
      TH1F *h = (TH1F*)f->Get(histName);
      if (!h) {
	std::cout << "Not able to find hist: " << histName << std::endl; 
	return NULL; 
      }
      if (!hall) hall = h;
      else hall->Add(h); 
    }
  }
  
  return hall; 
}


TCanvas* drawPt_phojn(std::vector<TString> inputFiles){
  int ww(1200), wh(400);
  set_root_style();
  TCanvas *c = new TCanvas("c", "Transverse momentum", ww, wh);
  c->Divide(3, 1); 
  c->cd(1); 
  TString dir = "#gamma+jets_pT-15to3000";

  std::vector<TString> hnames;
  hnames.push_back("eeeq0jets_qt");
  hnames.push_back("mumueq0jets_qt");
  TH1F *heq0jets = get_hist(dir, hnames, inputFiles);
  if (!heq0jets) return NULL;
  heq0jets->SetTitle("eq0 jets"); 
  heq0jets->Draw();
  TPad *c_1 = (TPad*) c->GetListOfPrimitives()->FindObject("c_1");
  c_1->SetLogx();
  c_1->SetLogy();

  c->cd(2); 
  hnames.clear();
  hnames.push_back("eegeq1jets_qt");
  hnames.push_back("mumugeq1jets_qt");
  TH1F *hgeq1jets = get_hist(dir, hnames, inputFiles);
  if (!hgeq1jets) return NULL; 
  hgeq1jets->SetTitle("geq1 jets"); 
  hgeq1jets->Draw();
  TPad *c_2 = (TPad*) c->GetListOfPrimitives()->FindObject("c_2");
  c_2->SetLogx();
  c_2->SetLogy();

  
  c->cd(3); 
  hnames.clear();
  hnames.push_back("eevbf_qt");
  hnames.push_back("mumuvbf_qt");
  TH1F *hvbf = get_hist(dir, hnames, inputFiles);
  if (!hvbf) return NULL; 
  hvbf->SetTitle("VBF"); 
  hvbf->Draw();
  TPad *c_3 = (TPad*) c->GetListOfPrimitives()->FindObject("c_3");
  c_3->SetLogx();
  c_3->SetLogy();

  c->Update(); 
  return c; 
}



TCanvas* drawPt(TString label, std::vector<TString> inputFiles){
  if (label == "phojn") {
    printf("label = %s \n", label.Data());
    return drawPt_phojn(inputFiles); 
  }
  
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
  printf("\nSYNOPSIS\n\tdrawPt LABEL input.root\n "); 
  // printf("\nSYNOPSIS\n\tdrawPt [-t hist-type ] [-opt draw-option]\n "); 
  // printf("\t[-h hist-name ] [-vmax max-value] [-npad num-pad] [-b] input1 input2 ...\n");
  printf("\nLABEL\n");
  printf("\t%-5s  %-40s\n", "pho", "photon pT");
  printf("\t%-5s  %-40s\n", "pho2", "triggered and selected photon pT");
  printf("\t%-5s  %-40s\n", "phojn", "photon pT with jet multiplicity");
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

