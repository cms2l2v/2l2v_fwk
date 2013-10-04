#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TObjArray.h"
#include "TMath.h"
#include "TGraphSmooth.h"

#include <iostream>
#include <vector>
#include <map>

using namespace std;

enum AnalysisMode {EWKVJJ, HZZ};
enum MCMode      { ALL,   PUREG };
void FitQtSpectrum(TString url="plotter.root", TString gUrl="plotter_gamma.root", int mcMode=PUREG, int anMode=EWKVJJ);
void printCMSheader(bool isSim);
void printCategoryToPlot(TString title);
TGraph *computeWeights(TH1F *target, TH1F *ctrl,TString name);
void formatQtHisto(TH1 *h,TString name, TString title,Int_t color, Int_t marker);


//
void formatQtHisto(TH1 *h,TString name, TString title,Int_t color, Int_t marker)
{
  if(h==0) return;
  h->SetDirectory(0); 
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetFillStyle(0);
  h->SetMarkerStyle(marker);
  h->SetName(name);   
  h->SetTitle(title);
  h->Rebin(4);
  h->GetXaxis()->SetTitle("Transverse momentum [GeV]");
  h->GetXaxis()->SetRangeUser(50,h->GetXaxis()->GetXmax());
  h->GetYaxis()->SetRangeUser(1e-1,1e5);
  h->GetYaxis()->SetTitle("Events");
}


//
TGraph *computeWeights(TH1F *target, TH1F *ctrl,TString name)
{
  if(target==0 || ctrl==0) return 0;

  //divide to compute the weights
  TH1F *ratio=(TH1F *) target->Clone(name+"raw");
  ratio->Divide(ctrl);
  TGraph *ratioGr=new TGraph(ratio);

  //smooth weights
  TGraphSmooth *gs = new TGraphSmooth(name+"smooth");
  //  TGraph *smoothWgtGr=gs->SmoothSuper(ratioGr,"",3);
  //  TGraph *smoothWgtGr=gs->SmoothLowess(ratioGr,"normal");
  //TGraph *smoothWgtGr=gs->SmoothLowess(ratioGr,"normal");

  Double_t xout[]={50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200,250,300};
  Int_t nout=sizeof(xout)/sizeof(Double_t);
  Double_t yout[nout];
  for(Int_t ip=0; ip<nout; ip++){
    yout[ip]=ratioGr->Eval(xout[ip]);
  }
  TGraph *smoothWgtGr=gs->Approx(ratioGr,"linear", nout, xout, 0, yout[nout-1],1);
  smoothWgtGr->SetName(name);
  smoothWgtGr->SetTitle(name);
  smoothWgtGr->SetFillColor(target->GetFillColor());
  smoothWgtGr->SetFillStyle(target->GetFillStyle());
  smoothWgtGr->SetLineColor(target->GetLineColor());
  smoothWgtGr->SetLineStyle(target->GetLineStyle());
  smoothWgtGr->SetMarkerColor(target->GetMarkerColor());
  smoothWgtGr->SetMarkerStyle(target->GetMarkerStyle());
  
  //return result
  delete ratio;
  delete ratioGr; 
  return smoothWgtGr;
}

//
void printCategoryToPlot(TString title)
{
  TPaveText *pave = new TPaveText(0.7,0.84,0.95,0.94,"brNDC");
  pave->SetFillStyle(0);
  pave->SetBorderSize(0);
  pave->AddText(title);
  pave->SetTextFont(42);
  pave->SetTextAlign(12);
  pave->Draw("same");
}

//
void printCMSheader(bool isSim)
{
  TPaveText *pave = new TPaveText(0.1,0.95,0.9,0.99,"brNDC");
  pave->SetFillStyle(0);
  pave->SetBorderSize(0);
  pave->SetTextAlign(12);
  //  pave->SetTextSize(0.08);
  if(isSim) pave->AddText("CMS simulation, #sqrt{s}=8 TeV");
  else      pave->AddText("CMS preliminary, #sqrt{s}=8 TeV");
  pave->Draw("same");
}

//
void FitQtSpectrum(TString url, TString gUrl, int mcMode, int catMode)
{
  std::vector<int> colors, markers;
  std::vector<TString> categs,titles,mcg;
  if(catMode==EWKVJJ){
    cout << "[FitQtSpectrum] Adding mjj categories" << endl;
    categs.push_back("mjjq016"); titles.push_back("M_{jj}<250");      colors.push_back(30); markers.push_back(20);
    categs.push_back("mjjq033"); titles.push_back("250<M_{jj}<350");  colors.push_back(32); markers.push_back(24);
    categs.push_back("mjjq049"); titles.push_back("350<M_{jj}<450");  colors.push_back(30); markers.push_back(24);
    categs.push_back("mjjq066"); titles.push_back("450<M_{jj}<550");  colors.push_back(32); markers.push_back(20);
    categs.push_back("mjjq083"); titles.push_back("550<M_{jj}<750");  colors.push_back(49); markers.push_back(20);
    categs.push_back("mjjq092"); titles.push_back("750<M_{jj}<1000"); colors.push_back(46); markers.push_back(24);
    categs.push_back("mjjq100"); titles.push_back("M_{jj}>1000");     colors.push_back(38); markers.push_back(24);
  }
  else{
    cout << "[FitQtSpectrum] Adding jet-bin categories" << endl;
    categs.push_back("eq0jets");  titles.push_back("=0 jets");      colors.push_back(30); markers.push_back(20);
    categs.push_back("geq1jets");  titles.push_back("#geq1 jets");  colors.push_back(32); markers.push_back(24);
    categs.push_back("vbf");      titles.push_back("VBF");          colors.push_back(30); markers.push_back(24); 
  }
  const size_t ncategs=categs.size();
  
  //mc for closure
  if(mcMode!=PUREG)
    {
      mcg.push_back("EWK #gammajj");
      mcg.push_back("V#gamma");
      mcg.push_back("Multijets");
    }
  mcg.push_back("#gamma+jets");
  TString mcdy("Z#rightarrow ll");

  //data to use
  TString llDir("data");
  if(url==gUrl) llDir="data (ll)";

  //variable to weight on
  TString llVar("qt");
  
  //plots to save in file
  TObjArray toSave;
  
  //prepare visualization
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);
  TCanvas *c= new TCanvas("c","c",500,1000);          c->Divide(1,3);
  TCanvas *mcc= new TCanvas("mcc","mcc",500,1000);    mcc->Divide(1,3);
  int ncdivs(ncategs/2);
  if(ncategs%2==1) ncdivs++;
  TCanvas *wc=new TCanvas("wc","wc",500,1000);        wc->Divide(2,ncdivs);
  TCanvas *mcwc=new TCanvas("mcwc","mcwc",1000,1000);  mcwc->Divide(2,ncdivs);


  //get histograms from files
  for(size_t icat=0; icat<ncategs; icat++)
    {
      //get histos from file
      TFile *fIn=TFile::Open(url);
      TFile *gIn=TFile::Open(gUrl);

      //mass distribution
      TH1F *eemass=(TH1F *) fIn->Get(llDir+"/ee"+categs[icat]+"_qmass");      
      if(eemass){
	eemass->SetName("ee"+categs[icat]+"_zmass");
	eemass->SetDirectory(0);
	toSave.Add(eemass);
      }

      //mass distributions
      TH1F *mmmass=(TH1F *) fIn->Get(llDir+"/mumu"+categs[icat]+"_qmass");      
      if(mmmass){
	mmmass->SetName("mumu"+categs[icat]+"_zmass");
	mmmass->SetDirectory(0);
	toSave.Add(mmmass);
      }

      //qt distributions
      TH1F *mmqt = (TH1F *)fIn->Get(llDir+"/mumu"+categs[icat]+"_"+llVar);
      TH1F *eeqt = (TH1F *)fIn->Get(llDir+"/ee"+categs[icat]+"_"+llVar);  
      TH1F *gqt  = (TH1F *)gIn->Get("data (#gamma)/mumu"+categs[icat]+"_qt");
      //if(categs[icat]=="mjjq092")
      //	{
      //	  if(mmqt) mmqt->Add( (TH1F *)fIn->Get(llDir+"/mumumjjq100_"+llVar) );
      //	  if(eeqt) eeqt->Add( (TH1F *)fIn->Get(llDir+"/eemjjq100_"+llVar) );
      //	  if(gqt) gqt->Add(  (TH1F *)gIn->Get("data (#gamma)/mumumjjq100_qt") );
      //	}
      formatQtHisto(mmqt,"mm"+categs[icat],   titles[icat], colors[icat],markers[icat]);
      formatQtHisto(eeqt,"ee"+categs[icat],   titles[icat], colors[icat],markers[icat]);
      formatQtHisto(gqt, "gamma"+categs[icat],titles[icat], colors[icat],markers[icat]);

     
      //qt distributions in MC
      TH1F *mcmmqt=(TH1F *)fIn->Get(mcdy+"/mumu"+categs[icat]+"_"+llVar);
      //if(mcmmqt)
      //	{
      //	  if(categs[icat]=="mjjq092") mcmmqt->Add( (TH1F *)fIn->Get(mcdy+"/mumumjjq100_"+llVar) );
      //	  mcmmqt->SetDirectory(0);
      //	}

      TH1F *mceeqt=(TH1F *)fIn->Get(mcdy+"/ee"+categs[icat]+"_"+llVar);
      //      if(mceeqt!=0)
      //	{
      //	  if(categs[icat]=="mjjq092") mceeqt->Add( (TH1F *)fIn->Get(mcdy+"/eemjjq100_"+llVar) );
      //	  mceeqt->SetDirectory(0);
      //	}
      
      TH1F *mcgqt=0;
      for(size_t imcg=0; imcg<mcg.size(); imcg++)
	{
	  TH1F *ih  = (TH1F *)gIn->Get(mcg[imcg]+"/mumu"+categs[icat]+"_qt");
	  if(ih!=0)
	    {
	      if(mcgqt==0)
		{
		  mcgqt=ih;
		  mcgqt->SetDirectory(0);  
		}
	      else
		mcgqt->Add(ih);
	      //if(categs[icat]=="mjjq092")
	      //	mcgqt->Add(  (TH1F *)gIn->Get(mcg[imcg]+"/mumumjjq100_qt") );
	    }
	}
      formatQtHisto(mcmmqt,"mcmm"+categs[icat],   titles[icat], colors[icat],markers[icat]);
      formatQtHisto(mceeqt,"mcee"+categs[icat],   titles[icat], colors[icat],markers[icat]);
      formatQtHisto(mcgqt, "mcgamma"+categs[icat],titles[icat], colors[icat],markers[icat]);

      fIn->Close();
      gIn->Close();
      
      //
      // PHOTON DATA -> REWEIGHT TO DILEPTON DATA
      //
      if(gqt)
	{      
	  if(mmqt) {c->cd(1); mmqt->Draw(icat==0 ? "e1" : "e1same"); }
	  if(eeqt) {c->cd(2); eeqt->Draw(icat==0 ? "e1" : "e1same"); }
	  if(gqt)  {c->cd(3); gqt ->Draw(icat==0 ? "e1" : "e1same"); }
	  
	  TGraph *eewgtGr=0,*mmwgtGr=0;
	  if(eeqt) eewgtGr = computeWeights(eeqt,gqt,"ee"  +categs[icat]+"_qt_datafitwgts");
	  if(mmqt) mmwgtGr = computeWeights(mmqt,gqt,"mumu"+categs[icat]+"_qt_datafitwgts");
	  TH1 *eeratio=(TH1 *) eeqt->Clone("ee"+categs[icat]+"ratio"); eeratio->Divide(gqt);
	  TH1 *mmratio=(TH1 *) mmqt->Clone("mm"+categs[icat]+"ratio"); mmratio->Divide(gqt);


	  wc->cd();
	  TPad *p=(TPad *) wc->cd(icat+1); 
	  p->SetLogx();
	  p->SetLogy();
	  bool fill(false);
	  //TGraph *frame=mmwgtGr;
	  TH1 *frame=mmratio;
	  if(eewgtGr)    {eeratio->Draw("e2"); eewgtGr->Draw("l");       fill=true; eewgtGr->SetLineColor(1); eewgtGr->SetLineWidth(2); frame=eeratio; toSave.Add(eewgtGr); eeratio->SetMarkerStyle(24); }
	  if(mmwgtGr)    {mmratio->Draw(fill ? "e2same" : "e2");  mmwgtGr->Draw("l");  fill=true; mmwgtGr->SetLineColor(1); mmwgtGr->SetLineWidth(1);       toSave.Add(mmwgtGr); mmratio->SetMarkerStyle(20); }
	
	  if(fill)
	    {
	      frame->GetXaxis()->SetTitle("Transverse momentum [GeV]");
	      frame->GetYaxis()->SetTitle("Weight");
	      frame->GetXaxis()->SetLabelSize(0.04);
	      frame->GetXaxis()->SetTitleSize(0.05);
	      frame->GetYaxis()->SetLabelSize(0.04);
	      frame->GetYaxis()->SetTitleSize(0.05);
	      frame->GetXaxis()->SetRangeUser(50,1000);
	      frame->GetYaxis()->SetRangeUser(1e-4,10.);
	      printCategoryToPlot(titles[icat]);
	      if(icat==0) {
		printCMSheader(false);
		TLegend *leg = new TLegend(0.2,0.84,0.5,0.94,"","brNDC");
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextFont(42);
		leg->SetNColumns(2);
		if(eewgtGr) leg->AddEntry(eewgtGr,"ee","f");
		if(mmwgtGr) leg->AddEntry(mmwgtGr,"#mu#mu","f");
		leg->Draw();
	      }
	    }
	}
      
      //
      // MC PHOTON SAMPLES -> REWEIGHT TO DATA
      //
      if(mcgqt)
	{
	  if(mcmmqt) {mcc->cd(1); mcmmqt->Draw(icat==0 ? "e1" : "e1same"); }
	  if(mceeqt) {mcc->cd(2); mceeqt->Draw(icat==0 ? "e1" : "e1same"); }
	  if(mcgqt)  {mcc->cd(3); mcgqt ->Draw(icat==0 ? "e1" : "e1same"); }
	  
	  TGraph *mceewgtGr=0,*mcmmwgtGr=0;
	  if(mceeqt) mceewgtGr = computeWeights(mceeqt,mcgqt,"ee"  +categs[icat]+"_qt_mcfitwgts");
	  if(mcmmqt) mcmmwgtGr = computeWeights(mcmmqt,mcgqt,"mumu"+categs[icat]+"_qt_mcfitwgts");
	  
	  
	  mcwc->cd();
	  TPad *p=(TPad *) mcwc->cd(icat+1); 
	  p->SetLogy();
	  bool fill(false);
	  TGraph *frame=mcmmwgtGr;
	  if(mceewgtGr)    { mceewgtGr->Draw("al");                fill=true; frame=mceewgtGr; mceewgtGr->SetLineWidth(2); toSave.Add(mceewgtGr); }
	  if(mcmmwgtGr)    { mcmmwgtGr->Draw(fill ? "l" : "al");   fill=true;                  mcmmwgtGr->SetLineWidth(1); toSave.Add(mcmmwgtGr); }
	  if(fill)
	    {
	      frame->GetXaxis()->SetTitle("Transverse momentum [GeV]");
	      frame->GetYaxis()->SetTitle("Weight");
	      frame->GetXaxis()->SetLabelSize(0.04);
	      frame->GetXaxis()->SetTitleSize(0.05);
	      frame->GetYaxis()->SetLabelSize(0.04);
	      frame->GetYaxis()->SetTitleSize(0.05);
	      frame->GetXaxis()->SetRangeUser(50,1000);
	      frame->GetYaxis()->SetRangeUser(1e-4,10.);
	      printCategoryToPlot(titles[icat]);
	      if(icat==0) {
		printCMSheader(true);
		TLegend *leg = new TLegend(0.2,0.84,0.5,0.94,"","brNDC");
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextFont(42);
		leg->SetNColumns(2);
		if(mceewgtGr) leg->AddEntry(mceewgtGr,"ee","f");
		if(mcmmwgtGr) leg->AddEntry(mcmmwgtGr,"#mu#mu","f");
		leg->Draw();
	      }
	    }
	}
    }


  TString channels[]={"#mu#mu events","ee events","#gamma events"};
  for(size_t i=1; i<=3; i++)
    {
      TPad *p=(TPad *)c->cd(i);
      p->SetLogy(); 
      p->SetLogx();
      printCategoryToPlot(channels[i-1]); 
      if(i==1)  { 
	p->SetTopMargin(0.05);
	printCMSheader(false);
	if(p->GetListOfPrimitives()->GetSize()>2){
	  TLegend *leg = p->BuildLegend(0.7,0.5,0.95,0.84);
	  leg->SetFillStyle(0);
	  leg->SetBorderSize(0);
	  leg->SetTextFont(42);
	}
      }

      p=(TPad *)mcc->cd(i);
      p->SetLogy();
      p->SetLogx();
      printCategoryToPlot(channels[i-1]);
      if(i==1) {
	printCMSheader(true);
	p->SetTopMargin(0.05);
	if(p->GetListOfPrimitives()->GetSize()>2){
	  TLegend *leg = p->BuildLegend(0.7,0.5,0.95,0.84);
	  leg->SetFillStyle(0);
	  leg->SetBorderSize(0);
	  leg->SetTextFont(42);
	}
      }
    }

  //save plots      
  c->Modified();   c->Update();   c->SaveAs("qtFit.png");      c->SaveAs("qtFit.pdf");       c->SaveAs("qtFit.C");
  mcc->Modified(); mcc->Update(); mcc->SaveAs("qtFit_mc.png"); mcc->SaveAs("qtFit_mc.pdf");  mcc->SaveAs("qtFit_mc.C");
  wc->SaveAs("qtWeights.png");       wc->SaveAs("qtWeights.pdf");       wc->SaveAs("qtWeights.C");
  mcwc->SaveAs("qtWeights_mc.png");  mcwc->SaveAs("qtWeights_mc.pdf");  mcwc->SaveAs("qtWeights_mc.C");

  //save results
  TFile *fOut=TFile::Open("gammawgts.root","RECREATE");
  fOut->cd();
  for(int i=0; i<toSave.GetEntriesFast(); i++)
    {
      if(toSave.At(i)==0) continue;
      toSave.At(i)->Write();
    }
  fOut->Close();
}
