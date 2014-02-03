#include "TCanvas.h"
#include "TString.h"
#include "TH1.h"
#include "TFile.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TObjArray.h"
#include "TGraphSmooth.h"
#include "TF1.h"
#include "TSystem.h"
#include "TMath.h"

#include <iostream>
#include <map>
#include <vector>

using namespace std;

TObjArray toSave;
std::map<TString,TString> systForClosure;

void closureTest(TFile *llfile, TFile *gfile,TString distr,TString ch, TString cat,bool purePhoton);
std::vector<TH1F *> getRatioOnly(TFile *llF,TFile *gF,TString distr,TString ch, TString cat, bool purePhoton);
void closureDY(TFile *llfile,TFile *gfile,TString distr,bool purePhoton);
void runVBFZClosure(TFile *llfile,TFile *gfile, TString outfile,bool purePhoton);
void runFinalClosure();

bool smoothFakesHisto=false;
TString dilCh="ll";
TFile *llInF=0;
std::vector<TFile *> gInF;

//
void runFinalHZZClosure()
{
  systForClosure.clear();

  dilCh="ll";
  //dilCh="ee";
  //dilCh="mumu";


  //open the files with the input plots
  TString llfile="~/work/hzz_5311/plotter_dy_closure.root";
  TFile *llInF=TFile::Open(llfile);

  TString gfile="~/work/hzz_5311/plotter_dy_closure_g_qt_pure.root";
  TFile *gInF=TFile::Open(gfile);

  TString distr[]={"met","mindphijmet"}; //"mt","mtNM1","axialmet","mindphijmet","mindphijmetNM1","balance"};
  TString cat[]={"","eq0jets","geq1jets","vbf"};
  const size_t ncat=sizeof(cat)/sizeof(TString);
  for(size_t icat=0; icat<ncat; icat++)
    {
      for(size_t ich=0; ich<sizeof(distr)/sizeof(TString); ich++) 
	{
	  closureTest(llInF,gInF,distr[ich],dilCh,cat[icat],true);
	}
    }
  
  gSystem->Exec("mkdir -p closure");
  gSystem->Exec("mv "+dilCh+"*closure*.* closure");
  
  //close all opened files
  llInF->Close();
  gInF->Close();
  toSave.Clear();
}

//
void runFinalClosure()
{
  systForClosure.clear();

  dilCh="ll";
  //dilCh="ee";
  //dilCh="mumu";


  //open the files with the input plots
  TString llfile="/afs/cern.ch/user/p/psilva/work/ewkzp2j_5311/plotter_dy_closure.root";
  llInF=TFile::Open(llfile);

  TString gfiles[]={
    "~/work/ewkzp2j_5311/plotter_dy_closure_g_qt_loose.root",
    "~/work/ewkzp2j_5311/plotter_dy_closure_g_qt_pure.root",
    "~/work/ewkzp2j_5311/plotter_dy_closure_g_qt_tight.root"
  };
  for(size_t i=0; i<3; i++) gInF.push_back( TFile::Open(gfiles[i]) );

  //run the closure test
  for(size_t i=0; i<3; i++) 
    {
      TString outfile(gSystem->BaseName(gfiles[i]));
      outfile.ReplaceAll("plotter","closure");

      bool purePhoton(false);
      TString outDir("tight");
      if(outfile.Contains("loose")) outDir="loose";
      if(outfile.Contains("pure")) { outDir="pureg"; purePhoton=true; }
      
      systForClosure[outDir]=outfile;
      runVBFZClosure(llInF,gInF[i],outfile,purePhoton);
      
      gSystem->Exec("mkdir -p "+outDir);
      gSystem->Exec("mv "+dilCh+"*closure*.* " + outDir);
    }

  //close all opened files
  llInF->Close();
  for(size_t i=0; i<3; i++) gInF[i]->Close();
  gInF.clear();
}

//
void runVBFZClosure(TFile *llfile,TFile *gfile, TString outfile, bool purePhoton)
{
  TString distr[]={
    "qt",                                                                                                     //boson qT
    "vbfcandjet1eta", "vbfcandjet2eta", "vbfcandjet1pt",     "vbfcandjet2pt", "vbfqgmva1", "vbfqgmva2",       //tag jets
    "vbfcandjetdeta", "vbfcandjetseta", "vbfcandjetetaprod", "vbfdphijj",     "vbfmjj",    "vbfspt", "BDTD", "MLP", //dijet 
    "vbfystar", 	     "vbfhardpt",                                                                            //dijet+Z
    "met",            "metL",                                                                                 //met
    "vbfcjv15",  "vbfhtcjv15",  "vbfmaxcjvjpt", "vbfystar3" //, "vbfcjv20", "vbfhtcjv20", "vbfcjv", "vbfhtcjv",            //central jet activity
  };
  for(size_t ich=0; ich<sizeof(distr)/sizeof(TString); ich++) closureDY(llfile,gfile,distr[ich],purePhoton);
  
  
  //save all the ratios to a file
  TFile *fOut=TFile::Open(outfile,"RECREATE");
  fOut->cd();
  for(int i=0; i<toSave.GetEntriesFast(); i++) toSave.At(i)->Write();
  fOut->Close();
  toSave.Clear();
}


//
void closureDY(TFile *llfile,TFile *gfile, TString distr,bool purePhoton)
{
  TString cat[]={"mjjq016","mjjq033","mjjq049","mjjq066","mjjq083","mjjq092","mjjq100","highmjj","mjjgt092","lowhardpt","highhardpt"};
  const size_t ncat=sizeof(cat)/sizeof(TString);
  for(size_t icat=0; icat<ncat; icat++)
    {
      if(distr=="Fisher" && icat==0) continue;
      closureTest(llfile,gfile,distr,dilCh,cat[icat],purePhoton);

      //
      //PDF variation stability for selected distributions
      //
      if(distr!="vbfcandjetdeta" && distr!="vbfcandjet1pt" && distr!="vbfcandjet2pt" && distr!="vbfcandjet1eta" && distr!="vbfcandjet2eta") continue;
      
      std::map<std::pair<int,int>, TGraph *> ratios;
      std::vector<TH1F *>baseLine;
      for(int ivar=0; ivar<=44; ivar++)
	{
	  TString var("_"); var+= ivar;
	  std::vector<TH1F *>result=getRatioOnly(llfile,gfile,distr+var,dilCh,cat[icat],purePhoton);
	  if(result.size()<3) continue;
	  if(baseLine.size()==0) { baseLine=result; } 
	  else{
	    for(size_t i=0; i<result.size(); i++){
	      result[i]->Add(baseLine[i],-1);
	      if(i<2) result[i]->Divide(baseLine[i]);
	      TGraph *gr=new TGraph(result[i]);
	      gr->SetName(result[i]->GetName());
	      gr->SetFillStyle(0);
	      gr->SetLineColor(kGray);
	      gr->SetMarkerColor(kGray);
	      gr->SetMarkerStyle(1);
	      std::pair<int,int> key(i,ivar);
	      ratios[key]=gr;
	    }
	  }
	}
      
      //sum in quadrature
      for(size_t i=0; i<baseLine.size(); i++)
	{
	  TString pf("dy"); if(i==1) pf="g"; if(i==2) pf="";
	  TGraph *totalUp=new TGraph;   totalUp->SetFillStyle(0);   totalUp->SetMarkerStyle(1);   totalUp->SetName(distr+pf+"Up");     totalUp->SetTitle("+1-#sigma");
	  TGraph *totalDown=new TGraph; totalDown->SetFillStyle(0); totalDown->SetMarkerStyle(1); totalDown->SetName(distr+pf+"Down"); totalDown->SetTitle("-1-#sigma");
	  Double_t c90(1.64485); //needed to normalize CTEQ variations to 68% CL
	  
	  std::pair<int,int> firstKey(i,1);
	  for(int ip=1; ip<=baseLine[i]->GetXaxis()->GetNbins(); ip++)
	    {
	      Double_t varUp(0),varDown(0),x(0);
	      
	      //for(std::map<int, TGraph *>::iterator it=ratios.begin(); it!=ratios.end(); it++)
	      for(int ivar=1; ivar<=44; ivar+=2)
		{
		  std::pair<int,int> key1(i,ivar),key2(i,ivar+1);
		  if(ratios.find(key1)==ratios.end() || ratios.find(key2)==ratios.end()) continue;
		  
		  Double_t y1,y2;
		  ratios[key1]->GetPoint(ip-1,x,y1);
		  ratios[key2]->GetPoint(ip-1,x,y2);
		  
		  varUp   += TMath::Power(TMath::Max( y1, y2), 2);
		  varDown += TMath::Power(TMath::Max(-y1,-y2), 2);
		} 
	      totalUp  ->SetPoint(ip-1,x,sqrt(varUp)/c90);
	      totalDown->SetPoint(ip-1,x,-sqrt(varDown)/c90);
	    }
	  
	  //show the results
	  
	  TCanvas *c=new TCanvas("c"+pf,"c"+pf,600,600*0.3);
	  totalUp->Draw("al");
	  totalUp->GetXaxis()->SetTitle( baseLine[i]->GetXaxis()->GetTitle());
	  TString title("#Delta QCD Z+2j");
	  if(i==1) title="#Delta QCD #gamma+2j";
	  if(i==2) title="#Delta closure test";
	  totalUp->GetYaxis()->SetTitle(title);
	  totalUp->GetYaxis()->SetRangeUser(-0.5,0.5);
	  totalUp->GetYaxis()->SetTitleOffset(0.5);
	  totalUp->GetXaxis()->SetTitleOffset(0.5);
	  totalUp->GetYaxis()->SetTitleSize(0.08);
	  totalUp->GetXaxis()->SetTitleSize(0.08);
	  totalUp->GetYaxis()->SetLabelSize(0.07);
	  totalUp->GetXaxis()->SetLabelSize(0.07);
	  totalDown->Draw("l");
	  for(int ivar=1; ivar<=44; ivar+=2)
	    {
	      std::pair<int,int> key(i,ivar);
	      ratios[key]->Draw("l");
	    }
	  /*
	    TLegend *leg=new TLegend(0.7,0.95,0.95,0.99);
	    leg->AddEntry(totalUp,totalUp->GetTitle(),"F");
	    leg->AddEntry(totalDown,totalDown->GetTitle(),"F");
	    leg->SetBorderSize(0);
	    leg->SetFillStyle(0);
	    leg->SetTextFont(42);
	    leg->SetTextAlign(11);
	    leg->SetTextSize(0.05);
	    leg->Draw("same");
	    leg->SetNColumns(2);
	  */
	  
	  /*
	    TPaveText *pave = new TPaveText(0.5,0.8,0.9,0.9,"brNDC");
	    pave->SetBorderSize(0);
	    pave->SetFillStyle(0);
	    pave->SetTextAlign(12);
	    pave->SetTextSize(0.08);
	    TString mjjCat("M_{jj}>1000");
	    if(cat[icat].Contains("mjjq016")) mjjCat="M_{jj}<250";
	    if(cat[icat].Contains("mjjq033")) mjjCat="250<M_{jj}<350";
	    if(cat[icat].Contains("mjjq049")) mjjCat="350<M_{jj}<450";
	    if(cat[icat].Contains("mjjq066")) mjjCat="450<M_{jj}<550";
	    if(cat[icat].Contains("mjjq083")) mjjCat="550<M_{jj}<750";
	    if(cat[icat].Contains("mjjq092")) mjjCat="750<M_{jj}<1000";
	    if(cat[icat].Contains("mjjgt092")) mjjCat="M_{jj}750";
	    if(cat[icat].Contains("highmjj")) mjjCat="M_{jj}>1250";
	    pave->AddText("CMS simulation, #sqrt{s}=8 TeV, "+mjjCat);
	    pave->Draw();
	  */
	  
	  c->SaveAs(dilCh+cat[icat]+"_"+pf+distr+"_pdf_closure.png");
	  c->SaveAs(dilCh+cat[icat]+"_"+pf+distr+"_pdf_closure.C");
	  c->SaveAs(dilCh+cat[icat]+"_"+pf+distr+"_pdf_closure.pdf");
	}
    }
}


//
std::vector<TH1F *> getRatioOnly(TFile *llF,TFile *gF,TString distr,TString ch, TString cat, bool purePhoton)
{
  std::vector<TH1F *> toReturn;

  bool rebin(distr.Contains("jetdeta") || distr.Contains("spt"));

  //
  //GET HISTOS FROM FILES
  //
  TString mcdy("Z#rightarrow ll");
  TH1 *hdy = 0;
  if(ch=="ll")
    {
      hdy=(TH1 *) llF->Get(mcdy+"/ee"+cat+"_"+distr);
      if(hdy)
	{
	  hdy= (TH1 *)hdy->Clone();
	  hdy->Add( (TH1 *) llF->Get(mcdy+"/mumu"+cat+"_"+distr) );
	}
    }
  else
    {
      hdy=(TH1 *) llF->Get(mcdy+"/"+ch+cat+"_"+distr);
      if(hdy) hdy=(TH1 *)hdy->Clone();
    }
  if(hdy==0) return toReturn;
  hdy->SetDirectory(0);
  if(rebin)hdy->Rebin();


  std::vector<TString> mcg;
  mcg.push_back("#gamma+jets");
  if(!purePhoton) mcg.push_back("Multijets");
  TH1 *hg=0;
  for(size_t ig=0; ig<mcg.size(); ig++)
    {
      if(ch=="ll")
	{
	  if(hg) {
	    hg->Add( (TH1 *)gF->Get(mcg[ig]+"/mumu"+cat+"_"+distr) );
	  }
	  else{
	    hg=(TH1 *)gF->Get(mcg[ig]+"/mumu"+cat+"_"+distr);
	    hg=(TH1 *)hg->Clone("mcg_"+cat+"_"+distr);
	  }
	  hg->Add( (TH1 *)gF->Get(mcg[ig]+"/ee"+cat+"_"+distr) );
	  
	}
      else
	{
	  if(hg){
	    hg->Add( (TH1 *) gF->Get(mcg[ig]+"/"+ch+cat+"_"+distr) );
	  }
	  else{
	    hg=(TH1 *) gF->Get(mcg[ig]+"/"+ch+cat+"_"+distr);
	    hg=(TH1 *)hg->Clone("mcg_"+ch+cat+"_"+distr);
	  }
	}
    }
  hg->SetDirectory(0);
  if(rebin)hg->Rebin();
  
  //compute the ratio
  TH1F *hratio=(TH1F *)hdy->Clone(distr+"_ratio");
  hratio->Divide(hg);
  toReturn.push_back((TH1F *)hdy);
  toReturn.push_back((TH1F *)hg);
  toReturn.push_back((TH1F *)hratio);
  return toReturn;
}


//
void closureTest(TFile *llF,TFile *gF,TString distr,TString ch, TString cat, bool purePhoton)
{

  bool rebin(distr.Contains("jetdeta") || distr.Contains("spt") || distr.Contains("qgmva"));

  //
  //GET HISTOS FROM FILES
  //
  TString mcdy("Z#rightarrow ll");
  TH1 *hdy = 0;
  if(ch=="ll")
    {
      hdy=(TH1 *) llF->Get(mcdy+"/ee"+cat+"_"+distr);	       
      hdy=(TH1 *) hdy->Clone("mcdy_"+cat+"_"+distr);
      hdy->Add((TH1 *) llF->Get(mcdy+"/mumu"+cat+"_"+distr) );
    }
  else
    {
      hdy=(TH1 *) llF->Get(mcdy+"/"+ch+cat+"_"+distr);
      hdy=(TH1 *) hdy->Clone("mcdy_"+ch+cat+"_"+distr); 
    }
  if(hdy==0) return;
  hdy->SetDirectory(0);
  if(rebin) hdy->Rebin();
  
  std::vector<TString> mcg;
  mcg.push_back("#gamma+jets");
  if(!purePhoton) mcg.push_back("Multijets");
  TH1 *hg=0, *hpureg=0;
  for(size_t ig=0; ig<mcg.size(); ig++)
    {
      if(ch=="ll")
	{
	  if(hg) {
	    hg->Add( (TH1 *)gF->Get(mcg[ig]+"/mumu"+cat+"_"+distr) );
	  }
	  else{
	    hg=(TH1 *)gF->Get(mcg[ig]+"/mumu"+cat+"_"+distr);
	    hg=(TH1 *)hg->Clone("mcg_"+cat+"_"+distr);
	  }
	  hg->Add( (TH1 *)gF->Get(mcg[ig]+"/ee"+cat+"_"+distr) );
	
	}
      else
	{
	  if(hg){
	    hg->Add( (TH1 *) gF->Get(mcg[ig]+"/"+ch+cat+"_"+distr) );
	  }
	  else{
	    hg=(TH1 *) gF->Get(mcg[ig]+"/"+ch+cat+"_"+distr);
	    hg=(TH1 *)hg->Clone("mcg_"+ch+cat+"_"+distr);
	  }
	}
      
      if(ig==0){
	TString pureName(hg->GetName());
	pureName.ReplaceAll("mcg","mcpureg");
	hpureg=(TH1 *)hg->Clone(pureName);
      }
    }
  if(hg==0 || hpureg==0) return;
  hg->SetDirectory(0);
  hpureg->SetDirectory(0);
  if(rebin) {
    hg->Rebin();
    hpureg->Rebin();
  }
  //if(distr=="qt") { hg->Rebin(4); hdy->Rebin(4); hpureg->Rebin(4); }
  //  else if(!distr.Contains("cjv")) { hg->Rebin(2); hdy->Rebin(2); hpureg->Rebin(2); }

  //
  //SHOW CLOSURE TEST
  //
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TCanvas *c = new TCanvas("c","c");
  c->SetWindowSize(500,500);
  c->cd();
  
  //distributions
  TPad *t1 = new TPad("p1","p1",0,0.3,1.0,1.0);
  t1->Draw();
  t1->cd();
  t1->SetTopMargin(0.08);
  t1->SetBottomMargin(0);
  t1->SetRightMargin(0.05);

  //find limits
  hdy->Scale(1./hdy->Integral());
  Double_t xmin(hdy->GetXaxis()->GetXmin());
  Double_t xmax(hdy->GetXaxis()->GetXmax());
  if(distr.Contains("met") && !distr.Contains("axial")) {xmin=0;   xmax=250;}
  if(distr.Contains("mt"))  {xmin=120; xmax=600;}
  float ymin(3e-4),ymax(hdy->GetMaximum()*1.2);

  //draw
  float gscale(hg->Integral());
  hg->Scale(1./gscale);
  hpureg->Scale(1./gscale);
  TH1 *hfakes=(TH1 *)hg->Clone("fakes");
  hfakes->SetDirectory(0);
  if(smoothFakesHisto) hfakes->Smooth();

  hg->SetTitle("QCD #gamma jj");
  hg->GetYaxis()->SetLabelSize(0.06);
  hg->GetYaxis()->SetTitleSize(0.05);      
  hg->GetYaxis()->SetTitleOffset(1.0);
  hg->GetYaxis()->SetLabelSize(0.05);
  hg->GetXaxis()->SetRangeUser(xmin,xmax); 
  if(distr.Contains("qt"))
    {
      ymin=3e-6;
      t1->SetLogy();
      t1->SetLogx();
      hg->GetXaxis()->SetRangeUser(50,1e3);
    }
  hg->GetYaxis()->SetRangeUser(ymin,ymax);
  hg->GetYaxis()->SetTitle("Events (a.u.)");
  hg->SetLineColor(1);
  hg->SetMarkerColor(1);
  hg->SetMarkerStyle(1);
  hg->SetFillStyle(1001);
  hg->SetFillColor(800);
  hg->Draw("hist");


  if(!purePhoton){
    hfakes->Add(hpureg,-1);
    hfakes->SetTitle("Fakes");
    hfakes->SetLineColor(1);
    hfakes->SetMarkerColor(1);
    hfakes->SetMarkerStyle(1);
    hfakes->SetFillStyle(3001);
    hfakes->SetFillColor(kGray);
    hfakes->Draw("histsame");   
  }


  hdy->SetTitle("QCD Z jj");
  hdy->Draw("e1same");

  bool setLogY(false);
  TString mjjCat("M_{jj}>1000");
  if(cat.Contains("mjjq016")) mjjCat="M_{jj}<250";
  if(cat.Contains("mjjq033")) mjjCat="250<M_{jj}<350";
  if(cat.Contains("mjjq049")) mjjCat="350<M_{jj}<450";
  if(cat.Contains("mjjq066")) mjjCat="450<M_{jj}<550";
  if(cat.Contains("mjjq083")) mjjCat="550<M_{jj}<750";
  if(cat.Contains("mjjq092")) mjjCat="750<M_{jj}<1000";
  if(cat.Contains("mjjgt092")) mjjCat="M_{jj}>750";
  if(cat.Contains("highmjj")) mjjCat="M_{jj}>1250";
  if(cat.Contains("highhardpt")) mjjCat="Hard p_{T}>50";
  if(cat.Contains("lowhardpt")) mjjCat="Hard p_{T}<50";
  if(cat=="") { mjjCat="inclusive"; }
  if(cat.Contains("eq0jets"))  {mjjCat="=0 jets"; setLogY=true;}
  if(cat.Contains("geq1jets")) {mjjCat="#geq1 jets"; setLogY=true;}
  if(cat.Contains("vbf"))      {mjjCat="VBF"; setLogY=true;}
  if(setLogY)t1->SetLogy();

  TPaveText *pave = new TPaveText(0.7,0.85,0.95,0.9,"brNDC");
  pave->SetBorderSize(0);
  pave->SetFillStyle(0);
  pave->SetTextAlign(32);
  pave->SetTextFont(42);
  pave->SetTextSize(0.05);
  pave->SetTextColor(kBlue);
  pave->AddText("["+mjjCat+"]");
  pave->Draw();

  TLegend *leg=new TLegend(0.6,0.95,0.95,0.98);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextAlign(12);
  leg->SetTextSize(0.05);
  //leg->SetNColumns(3);
  leg->SetNColumns(2);
  if(!setLogY)
    {
      leg->AddEntry(hdy,hdy->GetTitle(),"P");
      leg->AddEntry(hg,hg->GetTitle(),"F");
      if(!purePhoton) {
	leg->AddEntry(hfakes,hfakes->GetTitle(),"F");
	leg->SetNColumns(3);
      }
    }
  else
    {
      leg->AddEntry(hdy,"Z","P");
      leg->AddEntry(hg,"#gamma","F");
    }
  leg->Draw("same");


  pave = new TPaveText(0.94,0.4,1.0,0.83,"brNDC");
  pave->SetBorderSize(0);
  pave->SetFillStyle(0);
  pave->SetTextAlign(21);
  pave->SetTextFont(42);
  pave->SetTextColor(kGray+2);
  pave->SetTextSize(0.05);
  char buf[1000];  
  sprintf(buf,"#chi^{2}/ndof : %3.2f , K-S prob : %3.2f",hdy->Chi2Test(hg,"WWCHI2/NDF"),hdy->KolmogorovTest(hg,"") );
  pave->AddText(buf)->SetTextAngle(-90);
  pave->Draw();

  pave = new TPaveText(0.09,0.95,0.35,0.98,"NDC");
  pave->SetBorderSize(0);
  pave->SetFillStyle(0);
  pave->SetTextAlign(12);
  pave->SetTextSize(0.05);
  pave->AddText("CMS simulation, #sqrt{s}=8 TeV");
  pave->Draw();
  
  //closure
  c->cd();
  TPad *t2 = new TPad("p2","p2",0,0.0,1.0,0.3);
  t2->SetTopMargin(0);
  t2->SetBottomMargin(0.25);
  t2->SetRightMargin(0.05);
  t2->Draw();
  t2->cd();

  leg = new TLegend(0.15,0.82,0.9,0.92,"","brNDC");
  leg->SetBorderSize(0);
  leg->SetFillStyle(3001);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.11);
  leg->SetTextAlign(12);


  //mc stats
  TH1F *denRelUncH=(TH1F *) hg->Clone("mcrelunc");
  for(int xbin=1; xbin<=denRelUncH->GetXaxis()->GetNbins(); xbin++)
    {
      if(denRelUncH->GetBinContent(xbin)==0) continue;
      Double_t err=denRelUncH->GetBinError(xbin)/denRelUncH->GetBinContent(xbin);
      denRelUncH->SetBinContent(xbin,1);
      denRelUncH->SetBinError(xbin,err);
    }
  TGraphErrors *denRelUnc=new TGraphErrors(denRelUncH);
  denRelUnc->SetLineColor(1);
  denRelUnc->SetFillStyle(3001);
  denRelUnc->SetFillColor(kGray);
  denRelUnc->SetMarkerColor(1);
  denRelUnc->SetMarkerStyle(1);
  denRelUncH->Reset("ICE");       
  denRelUncH->Draw();
  denRelUnc->Draw("3");
  leg->AddEntry(denRelUnc,"stat unc.","f");
  denRelUncH->GetYaxis()->SetRangeUser(0.2,1.74);
  denRelUncH->GetXaxis()->SetTitle(hdy->GetXaxis()->GetTitle());
  denRelUncH->GetXaxis()->SetLabelSize(0.12);
  denRelUncH->GetXaxis()->SetTitleSize(0.14);
  denRelUncH->GetXaxis()->SetTitleOffset(0.8);
  denRelUncH->GetYaxis()->SetLabelSize(0.12);
  denRelUncH->GetYaxis()->SetNdivisions(5);
  denRelUncH->GetYaxis()->SetTitleSize(0.12);
  //gr->GetYaxis()->SetTitle("Pred. rel. bias");
  denRelUncH->GetYaxis()->SetTitle("Ratio");
  denRelUncH->GetYaxis()->SetTitleOffset(0.4);
  denRelUncH->SetMarkerStyle(20);
  denRelUncH->SetMarkerColor(1);
  denRelUncH->SetLineColor(1);
  if(distr.Contains("qt")) { t2->SetLogx(); denRelUncH->GetXaxis()->SetRangeUser(50,1000); }

  //now the actual ratio 
  TH1 *gr = (TH1 *) hdy->Clone("closuregr");
  gr->SetDirectory(0);
  gr->Divide(hg);
 
  //smooth curve and symmetrize it
  TGraphErrors *uncGr=new TGraphErrors(gr);   
  TGraph *uncGrUp=new TGraph(gr);     uncGrUp->SetName(ch+cat+"_"+distr+"_unc");  uncGrUp->SetLineWidth(2);
  TGraph *uncGrDown=new TGraph;                                               uncGrDown->SetLineWidth(2);
  for(int ip=0; ip<uncGrUp->GetN(); ip++)
    {
      Double_t x,y;
      uncGrUp->GetPoint(ip,x,y);
      uncGr->SetPoint(ip,x,1+TMath::Abs(y-1));
      uncGrUp->SetPoint(ip,x,1+TMath::Abs(y-1));
      uncGrDown->SetPoint(ip,x,1-TMath::Abs(y-1));
    }
  uncGrUp->SetLineWidth(2);
  uncGrDown->SetLineWidth(2);
  if(setLogY)   leg->AddEntry(uncGrUp,"Z/#gamma","l");
  else          leg->AddEntry(uncGrUp,"QCD Z jj/QCD #gamma jj","l");
  uncGrUp->Draw("l");
  uncGrDown->Draw("l");
  toSave.Add(uncGrUp);
  
  //add L-T and Pure-T systs
  if(systForClosure.find("loose")!=systForClosure.end()
     && systForClosure.find("tight")!=systForClosure.end()
     && systForClosure.find("pureg")!=systForClosure.end()
     )
    {
      TFile *looseF=TFile::Open(systForClosure["loose"]);
      TGraph *looseGr=(TGraph *)looseF->Get(ch+cat+"_"+distr+"_unc");
      looseF->Close();

      TFile *puregF=TFile::Open(systForClosure["pureg"]);
      TGraph *puregGr=(TGraph *)puregF->Get(ch+cat+"_"+distr+"_unc");
      puregF->Close();
      
      TGraph *tightGr=uncGrUp;
      
      TGraph *totalSystUp         = new TGraph;
      TGraph *loose2tightSystUp   = new TGraph;   loose2tightSystUp->SetLineColor(kGreen-3);    
      TGraph *loose2tightSystDown = new TGraph;   loose2tightSystDown->SetLineColor(kGreen-3);  
      TGraph *pureg2tightSystUp   = new TGraph;   pureg2tightSystUp->SetLineColor(kRed-3);    
      TGraph *pureg2tightSystDown = new TGraph;   pureg2tightSystDown->SetLineColor(kRed-3);  
      for(int ip=0; ip<looseGr->GetN(); ip++)
	{
	  Double_t x(0),yl(0),yt(0),yg(0);
	  looseGr->GetPoint(ip,x,yl);
	  puregGr->GetPoint(ip,x,yg);
	  tightGr->GetPoint(ip,x,yt);
 
	  Double_t iDiff(TMath::Abs(yl-yt));
	  Double_t diffUp  (1+iDiff);
	  Double_t diffDown(1-iDiff);
	  loose2tightSystUp->SetPoint(ip,x,diffUp);
	  loose2tightSystDown->SetPoint(ip,x,diffDown);
	  
	  iDiff=TMath::Abs(yg-yt);
	  diffUp  =(1+iDiff);
	  diffDown=(1-iDiff);
	  pureg2tightSystUp->SetPoint(ip,x,diffUp);
	  pureg2tightSystDown->SetPoint(ip,x,diffDown);

	  //final uncertainty
	  // 2% PDF envelope
	  // loose2tight taken from data
	  // pure2tight is taken from MC 	  
	  Double_t totalDiff(pow(0.02,2)); 
	  totalDiff += pow(yg-yt, 2);
	  totalDiff += pow(yl-yt, 2);
	  totalDiff += 1+sqrt(totalDiff);
	  totalSystUp->SetPoint(ip,x,totalDiff);
	}

      loose2tightSystUp->SetLineWidth(2);
      loose2tightSystDown->SetLineWidth(2);
      leg->AddEntry(loose2tightSystUp,"loose-tight","l");
      loose2tightSystUp->Draw("l");
      loose2tightSystDown->Draw("l");

      pureg2tightSystUp->SetLineWidth(2);
      pureg2tightSystDown->SetLineWidth(2);
      leg->AddEntry(pureg2tightSystUp,"pure-tight","l");
      pureg2tightSystUp->Draw("l");
      pureg2tightSystDown->Draw("l");
      
      //only need to save one of the variations
      loose2tightSystUp->SetName(ch+cat+"_"+distr+"_ltunc");
      pureg2tightSystUp->SetName(ch+cat+"_"+distr+"_ptunc");
      totalSystUp->SetName(ch+cat+"_"+distr+"_totalunc");

      TF1 *ffunc=new TF1(ch+"_"+distr+"_func","max([0]*((1.+[1]+x)/(1.+[2]*x)),0.)",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());
      TF1 *ltfunc=(TF1 *) ffunc->Clone(ch+"_"+distr+"_ltfunc");
      TF1 *ptfunc=(TF1 *) ffunc->Clone(ch+"_"+distr+"_ptfunc");

      loose2tightSystUp->Fit(ltfunc,"MQE0N+");
      pureg2tightSystUp->Fit(ptfunc,"MQE0N+");
      totalSystUp->Fit(ffunc,"MQE0N+");
      
      toSave.Add(totalSystUp);
      toSave.Add(ffunc);
      toSave.Add(ltfunc);
      toSave.Add(loose2tightSystUp);
      toSave.Add(ptfunc);
      toSave.Add(pureg2tightSystUp);
    }

  //  TLine *l = new TLine(xmin,0.5,xmax,0.5);
  //l->SetLineColor(kRed);
  //l->Draw();
  //l = new TLine(xmin,1.5,xmax,1.5);
  //l->SetLineColor(kRed);
  //l->Draw();
  
  leg->Draw();
  leg->SetNColumns(4);

  c->Modified();
  c->Update();

  c->SaveAs(ch+cat+"_"+distr+"_closure.png");
  c->SaveAs(ch+cat+"_"+distr+"_closure.C");
  c->SaveAs(ch+cat+"_"+distr+"_closure.pdf");
}
