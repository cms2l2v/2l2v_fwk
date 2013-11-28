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
TGraph *computeWeights(TH1F *target, TH1F *ctrl,TString name, bool smooth=false);
void formatHisto(TH1 *h,TString name, TString title,Int_t color, Int_t marker,Int_t rebin);


//
void formatHisto(TH1 *h,TString name, TString title,Int_t color, Int_t marker, Int_t rebin)
{
  if(h==0) return;
  h->SetDirectory(0); 
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetFillStyle(0);
  h->SetMarkerStyle(marker);
  h->SetName(name);   
  h->SetTitle(title);
  if(rebin>0) h->Rebin(rebin);
  h->GetXaxis()->SetTitle("Transverse momentum [GeV]");
  h->GetXaxis()->SetRangeUser(50,h->GetXaxis()->GetXmax());
  h->GetYaxis()->SetRangeUser(1e-1,1e8);
  h->GetYaxis()->SetTitle("Events");
}


//
TGraph *computeWeights(TH1F *target, TH1F *ctrl,TString name,bool smooth)
{
  if(target==0 || ctrl==0) return 0;

  //divide to compute the weights
  TH1F *ratio=(TH1F *) target->Clone(name+"raw");
  ratio->Divide(ctrl);
  TGraph *ratioGr=new TGraph(ratio);

  if(!smooth) {   ratioGr->SetName(name); return ratioGr; }

  //smooth weights
  TGraphSmooth *gs = new TGraphSmooth(name);
  TGraph *smoothWgtGr=gs->SmoothSuper(ratioGr,"",3);
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
  TPaveText *pave = new TPaveText(0.65,0.9,0.9,0.94,"brNDC");
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
  std::vector<TString> categs,titles, photonMC;

  TString kinVar("qt"),kinVarTitle("Transverse momentum [GeV]");
  TString ewkSuppVar("mjj"),ewkSuppVarTitle("E_{T}^{miss}/q_{T}");
  TString ewkPhotonContribution("EWK #gammajj");
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
    kinVar="qt";                       kinVarTitle="Transverse momentum [GeV]";
    ewkSuppVar="balance";              ewkSuppVarTitle="E_{T}^{miss}/q_{T}";
    ewkPhotonContribution="V#gamma";
    categs.push_back("eq0jets");  titles.push_back("=0 jets");      colors.push_back(1); markers.push_back(20);
    categs.push_back("geq1jets"); titles.push_back("#geq1 jets");   colors.push_back(2); markers.push_back(24);
    categs.push_back("vbf");      titles.push_back("VBF");          colors.push_back(3); markers.push_back(20); 
  }
  const size_t ncategs=categs.size();
  
  //mc for closure
  if(mcMode!=PUREG)
    {
      photonMC.push_back("EWK #gammajj");
      photonMC.push_back("V#gamma");
      photonMC.push_back("Multijets");
    }
  photonMC.push_back("#gamma+jets");
  TString mcdy("Z#rightarrow ll");

  //data to use
  TString llDir("data");
  if(url==gUrl) llDir="data (ll)";
 
  //plots to save in file
  TObjArray toSave;
  
  //prepare visualization
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);

  //get histograms from files
  int ncdivs(ncategs/2);
  if(ncategs%2==1) ncdivs++;
  TString wgtVars[]     = {kinVar,      ewkSuppVar      };
  TString wgtTitles[]   = {kinVarTitle, ewkSuppVarTitle };
  Int_t   rebinFactor[] = {2, 0};
  TString channels[]={"#mu#mu events","ee events","#gamma events"};
  for(size_t iw=0; iw<2; iw++)
    {
      TString var=wgtVars[iw];
      TString varTitle=wgtTitles[iw];
      Int_t rebin=rebinFactor[iw];
      TString pf("");
      if(iw==1) pf="_ewksup";

      TCanvas *c= new TCanvas("c","c",500,1000);          c->Divide(1,3);
      TCanvas *mcc= new TCanvas("mcc","mcc",500,1000);    mcc->Divide(1,3);
      
      TCanvas *wc=new TCanvas("wc","wc",500,1000);        wc->Divide(2,ncdivs);
      TCanvas *ewkwc=new TCanvas("ewkwc","ewkwc",500,1000); ewkwc->Divide(2,ncdivs);
      TCanvas *mcwc=new TCanvas("mcwc","mcwc",1000,1000);  mcwc->Divide(2,ncdivs);
      
      for(size_t icat=0; icat<ncategs; icat++)
	{
	  //get histos from file
	  TFile *fIn=TFile::Open(url);
	  TFile *gIn=TFile::Open(gUrl);

	  //mass distributions
	  if(iw==0)
	    {
	      TH1F *eemass=(TH1F *) fIn->Get(llDir+"/ee"+categs[icat]+"_qmass");      
	      if(eemass){
		eemass->SetName("ee"+categs[icat]+"_zmass");
		eemass->SetDirectory(0);
		toSave.Add(eemass);
	      }
	      TH1F *mmmass=(TH1F *) fIn->Get(llDir+"/mumu"+categs[icat]+"_qmass");      
	      if(mmmass){
		mmmass->SetName("mumu"+categs[icat]+"_zmass");
		mmmass->SetDirectory(0);
		toSave.Add(mmmass);
	      }
	      TH1F *eemcmass=(TH1F *) fIn->Get(mcdy+"/ee"+categs[icat]+"_qmass");      
	      if(eemcmass){
		eemcmass->SetName("ee"+categs[icat]+"_mczmass");
		eemcmass->SetDirectory(0);
		toSave.Add(eemcmass);
	      }
	      TH1F *mmmcmass=(TH1F *) fIn->Get(mcdy+"/mumu"+categs[icat]+"_qmass");      
	      if(mmmcmass){
		mmmcmass->SetName("mumu"+categs[icat]+"_mczmass");
		mmmcmass->SetDirectory(0);
		toSave.Add(mmmcmass);
	      }
	    }
	  
	  //kinematics distribution in data
	  TH1F *mm = (TH1F *)fIn->Get(llDir+"/mumu"+categs[icat]+"_"+var);
	  TH1F *ee = (TH1F *)fIn->Get(llDir+"/ee"+categs[icat]+"_"+var);  
	  TH1F *g  = (TH1F *)gIn->Get("data (#gamma)/mumu"+categs[icat]+"_"+var);
	  formatHisto(mm,"mm"+categs[icat]+pf,   titles[icat], colors[icat],markers[icat],rebin);
	  formatHisto(ee,"ee"+categs[icat]+pf,   titles[icat], colors[icat],markers[icat],rebin);
	  formatHisto(g, "gamma"+categs[icat]+pf,titles[icat], colors[icat],markers[icat],rebin);
     
	  //kinematics distribution in MC
	  TH1F *mcmm=(TH1F *)fIn->Get(mcdy+"/mumu"+categs[icat]+"_"+var);
	  TH1F *mcee=(TH1F *)fIn->Get(mcdy+"/ee"+categs[icat]+"_"+var);
	  TH1F *mcg=0;
	  for(size_t imcg=0; imcg<photonMC.size(); imcg++)
	    {
	      TH1F *ih  = (TH1F *)gIn->Get(photonMC[imcg]+"/mumu"+categs[icat]+"_"+var);
	      if(ih!=0)
		{
		  if(mcg==0)
		    {
		      mcg=ih;
		      mcg->SetDirectory(0);  
		    }
		  else
		    mcg->Add(ih);
		}
	    }
	  TH1F *mcewkg=(TH1F *)gIn->Get(ewkPhotonContribution+"/mumu"+categs[icat]+"_"+var);
	  if(mcewkg) {
	    mcewkg->SetDirectory(0);
	    formatHisto(mcewkg, "mcewkg"+categs[icat]+pf,titles[icat], colors[icat],markers[icat],rebin);
	    mcewkg->SetLineStyle(7);
	    if(iw) mcewkg->Smooth();
	  }
	  formatHisto(mcmm,"mcmm"+categs[icat]+pf,   titles[icat], colors[icat],markers[icat],rebin);
	  formatHisto(mcee,"mcee"+categs[icat]+pf,   titles[icat], colors[icat],markers[icat],rebin);
	  formatHisto(mcg, "mcgamma"+categs[icat]+pf,titles[icat], colors[icat],markers[icat],rebin);

	  //all done with the files
	  fIn->Close();
	  gIn->Close();


	  //special procedure for EWK ad-hoc suppression
	  if(mcewkg && iw==1)
	    {
	      mm=(TH1F *) g->Clone("mumu"+categs[icat]+"_"+var);
	      mm->SetDirectory(0);
	      ee=(TH1F *) g->Clone("ee"+categs[icat]+"_"+var);
	      ee->SetDirectory(0);
	      cout <<  "Dilepton histograms to be replaced with EWK subtracted photon data plots for " << var << endl;
 
	      //re-scale from sideband
	      if(catMode==HZZ)
		{
		  Float_t minVar(2.0);
		  if(icat==0) minVar=1.5;
		  Int_t ibin=mcewkg->GetXaxis()->FindBin(minVar);
		  Int_t ebin=mcewkg->GetXaxis()->GetNbins();
		  Float_t mcExpected=mcewkg->Integral(ibin,ebin+1);
		  Float_t dataObserved=g->Integral(ibin,ebin+1);
		  Float_t ewksf(1.0);
		  if(mcExpected>0 && dataObserved>0) ewksf=dataObserved/mcExpected;
		  mm->Add(mcewkg,-ewksf);
		  ee->Add(mcewkg,-ewksf);
		  cout << "...applying EWK scale factor of " << ewksf << endl;
		}
	      else
		{
		  mm->Add(mcewkg,-1);
		  ee->Add(mcewkg,-1);
		}

	      //no negative entries
	      for(int ibin=1; ibin<mm->GetXaxis()->GetNbins(); ibin++)
		{
		  Float_t y=mm->GetBinContent(ibin);
		  if(y<0) mm->SetBinContent(ibin,0.);
		  y=ee->GetBinContent(ibin);
		  if(y<0) ee->SetBinContent(ibin,0.);

		}
	    }
	  
	  
	  // Derive weights to apply to data
	  if(g && (mm || ee))
	    {      
	      if(mm) {c->cd(1); mm->Draw(icat==0 ? "e1" : "e1same"); }
	      if(ee) {c->cd(2); ee->Draw(icat==0 ? "e1" : "e1same"); }
	      if(g)  {c->cd(3); g ->Draw(icat==0 ? "e1" : "e1same"); if(iw==1) mcewkg->Draw("histsame"); }
	  
	      TGraph *eewgtGr=0,*mmwgtGr=0;
	      if(ee) eewgtGr = computeWeights(ee,g,"ee"  +categs[icat]+"_"+var+"_datafitwgts", (iw==1));
	      if(mm) mmwgtGr = computeWeights(mm,g,"mumu"+categs[icat]+"_"+var+"_datafitwgts", (iw==1));
	      TH1 *eeratio=(TH1 *) ee->Clone("ee"+categs[icat]+"ratio"+pf); eeratio->Divide(g);
	      TH1 *mmratio=(TH1 *) mm->Clone("mm"+categs[icat]+"ratio"+pf); mmratio->Divide(g);

	      wc->cd();
	      TPad *p=(TPad *) wc->cd(icat+1); 
	      if(iw==0) 
		{
		  p->SetLogx();
		  p->SetLogy();
		}
	      bool fill(false);
	      TH1 *frame=mmratio;
	      if(eewgtGr)    {eeratio->Draw("e2");                    eewgtGr->Draw("l");  fill=true; eewgtGr->SetLineColor(1); eewgtGr->SetLineWidth(2); frame=eeratio; toSave.Add(eewgtGr); eeratio->SetMarkerStyle(24); }
	      if(mmwgtGr)    {mmratio->Draw(fill ? "e2same" : "e2");  mmwgtGr->Draw("l");  fill=true; mmwgtGr->SetLineColor(1); mmwgtGr->SetLineWidth(1);                toSave.Add(mmwgtGr); mmratio->SetMarkerStyle(20); }
	      if(fill)
		{
		  frame->GetXaxis()->SetTitle(varTitle);
		  frame->GetYaxis()->SetTitle("Weight");
		  frame->GetXaxis()->SetLabelSize(0.04);
		  frame->GetXaxis()->SetTitleSize(0.05);
		  frame->GetYaxis()->SetLabelSize(0.04);
		  frame->GetYaxis()->SetTitleSize(0.05);
		  frame->GetXaxis()->SetRangeUser(50,1000);
		  frame->GetYaxis()->SetRangeUser(1e-4,1.);
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
	  if(mcg)
	    {
	      if(mcmm) {mcc->cd(1); mcmm->Draw(icat==0 ? "e1" : "e1same"); }
	      if(mcee) {mcc->cd(2); mcee->Draw(icat==0 ? "e1" : "e1same"); }
	      if(mcg)  {mcc->cd(3); mcg ->Draw(icat==0 ? "e1" : "e1same"); }
	        
	      TGraph *mceewgtGr=0,*mcmmwgtGr=0;
	      if(mcee) mceewgtGr = computeWeights(mcee,mcg,"ee"  +categs[icat]+"_"+var+"_mcfitwgts", (iw==1));
	      if(mcmm) mcmmwgtGr = computeWeights(mcmm,mcg,"mumu"+categs[icat]+"_"+var+"_mcfitwgts", (iw==1));
	  	  
	      mcwc->cd();
	      TPad *p=(TPad *) mcwc->cd(icat+1); 
	      if(iw==0) p->SetLogx();
	      //p->SetLogy();
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
		  frame->GetYaxis()->SetRangeUser(1e-4,1.);
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
    

      //finish with the label and legends (not sure how this cam here)
      for(size_t i=1; i<=3; i++)
	{
	  TPad *p=(TPad *)c->cd(i);
	  if(iw==0) p->SetLogx();
	  p->SetLogy(); 
	  printCategoryToPlot(channels[i-1]); 
	  if(i==1)  { 
	    p->SetTopMargin(0.05);
	    printCMSheader(false);
	    if(p->GetListOfPrimitives()->GetSize()>2){
	      TLegend *leg = p->BuildLegend(0.65,0.6,0.9,0.88);
	      leg->SetFillStyle(0);
	      leg->SetBorderSize(0);
	      leg->SetTextFont(42);
	    }
	  }
	  
	  p=(TPad *)mcc->cd(i);
	  if(iw==0) p->SetLogx();
	  p->SetLogy();
	  printCategoryToPlot(channels[i-1]);
	  if(i==1) {
	    printCMSheader(true);
	    p->SetTopMargin(0.05);
	    if(p->GetListOfPrimitives()->GetSize()>2){
	      TLegend *leg = p->BuildLegend(0.65,0.6,0.9,0.88);
	      leg->SetFillStyle(0);
	      leg->SetBorderSize(0);
	      leg->SetTextFont(42);
	    }
	  }
	}

      //save plots      
      c->Modified();    c->Update();    c->SaveAs(var+"Distribution.png");       c->SaveAs(var+"Distribution.pdf");       c->SaveAs(var+"Distribution.C");
      mcc->Modified();  mcc->Update();  mcc->SaveAs(var+"Distribution_mc.png");  mcc->SaveAs(var+"Distribution_mc.pdf");  mcc->SaveAs(var+"Distribution_mc.C");
      wc->Modified();   wc->Update();   wc->SaveAs(var+"Weights.png");           wc->SaveAs(var+"Weights.pdf");           wc->SaveAs(var+"Weights.C");
      mcwc->Modified(); mcwc->Update(); mcwc->SaveAs(var+"Weights_mc.png");      mcwc->SaveAs(var+"Weights_mc.pdf");      mcwc->SaveAs(var+"Weights_mc.C");
    }
  
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
