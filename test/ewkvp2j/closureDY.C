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

#include <iostream>
#include <map>

using namespace std;

TObjArray toSave;
std::map<TString,TString> systForClosure;

void closureTest(TString url,TString gUrl,TString distr,TString ch, TString cat,bool purePhoton);
TH1F *getRatioOnly(TString url,TString gUrl,TString distr,TString ch, TString cat, bool purePhoton);
void closureDY(TString url,TString gUrl,TString distr,bool purePhoton);
void runVBFZClosure(TString llfile,TString gfile, TString outfile,bool purePhoton);
void runFinalClosure();


//
void runFinalClosure()
{
  systForClosure.clear();

  TString llfile="/afs/cern.ch/user/p/psilva/work/ewkzp2j_539/plotter_dy_closure.root";
  
  TString gfiles[]={
    "~/work/ewkzp2j_539/plotter_dy_closure_g_qt_loose.root",
    "~/work/ewkzp2j_539/plotter_dy_closure_g_qt_pure.root",
    "~/work/ewkzp2j_539/plotter_dy_closure_g_qt_tight.root"
  };
  
  for(size_t i=0; i<3; i++)
    {
      TString outfile(gSystem->BaseName(gfiles[i]));
      outfile.ReplaceAll("plotter","closure");

      bool purePhoton(false);
      TString outDir("tight");
      if(outfile.Contains("loose")) outDir="loose";
      if(outfile.Contains("pure")) { outDir="pureg"; purePhoton=true; }
      
      systForClosure[outDir]=outfile;
      runVBFZClosure(llfile,gfiles[i],outfile,purePhoton);
      
      gSystem->Exec("mkdir -p "+outDir);
      gSystem->Exec("mv *.{png,pdf} " + outDir);
    }
}

//
void runVBFZClosure(TString llfile,TString gfile, TString outfile, bool purePhoton)
{

  
  TString distr[]={//"Fisher", "LikelihoodD", "BDTD", 
		   "vbfcandjetdeta", "vbfcandjetetaprod", "vbfcandjet1eta",   "vbfcandjet2eta", "vbfcandjet1pt","vbfcandjet2pt","vbfdphijj",  
		   "vbfspt", "vbfhardpt",  "met", "metL",           
		   "vbfcjv15",       "vbfhtcjv15",   "vbfcjv20",     "vbfhtcjv20", "vbfcjv", "vbfhtcjv", "vbfmaxcjvjpt"};
  for(size_t ich=0; ich<sizeof(distr)/sizeof(TString); ich++) closureDY(llfile,gfile,distr[ich],purePhoton);
  
  
  //save all the ratios to a file
  TFile *fOut=TFile::Open(outfile,"RECREATE");
  fOut->cd();
  for(int i=0; i<toSave.GetEntriesFast(); i++) toSave.At(i)->Write();
  fOut->Close();
  toSave.Clear();
}

//
void closureDY(TString url,TString gUrl, TString distr,bool purePhoton)
{
  //TString dilCh="ll";
  TString dilCh="ee";
  //  TString dilCh="mumu";
  TString cat[]={"mjjq016","mjjq033","mjjq049","mjjq066","mjjq083","mjjq092","mjjq100"};
  const size_t ncat=sizeof(cat)/sizeof(TString);
  for(size_t icat=0; icat<ncat; icat++)
    {
      closureTest(url,gUrl,distr,dilCh,cat[icat],purePhoton);

      //variation stability
      if(distr!="vbfcandjetdeta") continue;
      std::map<int, TGraph *> ratios;
      TH1F *baseLine=0;
      for(int ivar=0; ivar<=44; ivar++)
	{
	  TString var("_"); var+= ivar;
	  TH1F *iratio=getRatioOnly(url,gUrl,distr+var,dilCh,cat[icat],purePhoton);
	  if(iratio==0) continue;
	  if(ivar==0) baseLine=iratio;
	  else if(baseLine){
	    iratio->Add(baseLine,-1);
	    ratios[ivar]=new TGraph(iratio);
	    ratios[ivar]->SetName(iratio->GetName());
	    ratios[ivar]->SetFillStyle(0);
	    ratios[ivar]->SetLineColor(kGray);
	    ratios[ivar]->SetMarkerColor(kGray);
	    ratios[ivar]->SetMarkerStyle(1);
	  }
	}

      //sum in quadrature
      if(baseLine!=0){
	TGraph *totalUp=new TGraph;   totalUp->SetFillStyle(0);   totalUp->SetMarkerStyle(1);   totalUp->SetName(distr+"Up");     totalUp->SetTitle("+1-#sigma");
	TGraph *totalDown=new TGraph; totalDown->SetFillStyle(0); totalDown->SetMarkerStyle(1); totalDown->SetName(distr+"Down"); totalDown->SetTitle("-1-#sigma");
	Double_t c90(1.64485); //needed to normalize CTEQ variations to 68% CL
	for(int ip=1; ip<=baseLine->GetXaxis()->GetNbins(); ip++)
	  {
	    Double_t varUp(0),varDown(0),x(0);
	    
	    //for(std::map<int, TGraph *>::iterator it=ratios.begin(); it!=ratios.end(); it++)
	    for(int ivar=1; ivar<=(int)ratios.size(); ivar+=2)
	      {
		if(ratios.find(ivar)==ratios.end() || ratios.find(ivar+1)==ratios.end()) continue;
		
		Double_t y1,y2;
		ratios[ivar]->GetPoint(ip-1,x,y1);
		ratios[ivar+1]->GetPoint(ip-1,x,y2);
		
		varUp   += pow(max( y1, y2), 2);
		varDown += pow(max(-y1,-y2), 2);
	      } 
	    totalUp  ->SetPoint(ip-1,x,sqrt(varUp)/c90);
	    totalDown->SetPoint(ip-1,x,-sqrt(varDown)/c90);
	  }
      
	//show the results
	
	TCanvas *c=new TCanvas("c","c",600,600);
	totalUp->Draw("al");
	totalUp->GetXaxis()->SetTitle( baseLine->GetXaxis()->GetTitle());
	totalUp->GetYaxis()->SetTitle( "Difference");
	totalUp->GetYaxis()->SetRangeUser(-0.5,0.5);
	totalDown->Draw("l");
	for(std::map<int, TGraph *>::iterator it=ratios.begin(); it!=ratios.end(); it++) it->second->Draw("l");
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
	
	TPaveText *pave = new TPaveText(0.5,0.8,0.8,0.9,"brNDC");
	pave->SetBorderSize(0);
	pave->SetFillStyle(0);
	pave->SetTextAlign(11);
	pave->SetTextFont(42);
	pave->SetTextSize(0.06);
	TString mjjCat("M_{jj}>1000");
	if(cat[icat].Contains("mjjq016")) mjjCat="M_{jj}<250";
	if(cat[icat].Contains("mjjq033")) mjjCat="250<M_{jj}<350";
	if(cat[icat].Contains("mjjq049")) mjjCat="350<M_{jj}<450";
	if(cat[icat].Contains("mjjq066")) mjjCat="450<M_{jj}<550";
	if(cat[icat].Contains("mjjq083")) mjjCat="550<M_{jj}<750";
	if(cat[icat].Contains("mjjq092")) mjjCat="750<M_{jj}<1000";
	pave->AddText(mjjCat);
	pave->Draw();

	pave = new TPaveText(0.1,0.95,0.6,0.99,"NDC");
	pave->SetBorderSize(0);
	pave->SetFillStyle(0);
	pave->SetTextAlign(11);
	pave->SetTextSize(0.04);
	pave->AddText("CMS simulation, #sqrt{s}=8 TeV");
	pave->Draw();
	
	c->SaveAs(dilCh+cat[icat]+"_"+distr+"_pdf_closure.png");
	c->SaveAs(dilCh+cat[icat]+"_"+distr+"_pdf_closure.pdf");
      }
    }
}


//
TH1F *getRatioOnly(TString url,TString gUrl,TString distr,TString ch, TString cat, bool purePhoton)
{
  //
  //GET HISTOS FROM FILES
  //
  TString mcdy("Z#rightarrow ll");
  TFile *llF = TFile::Open(url);
  TH1 *hdy = 0;
  if(ch=="ll")
    {
      hdy=(TH1 *) llF->Get(mcdy+"/ee"+cat+"_"+distr);
      if(hdy)
	{
	  hdy->Add( (TH1 *) llF->Get(mcdy+"/mumu"+cat+"_"+distr) );
	}
    }
  else
    {
      hdy=(TH1 *) llF->Get(mcdy+"/"+ch+cat+"_"+distr);
    }
  if(hdy==0) return 0;
  hdy->SetDirectory(0);
  llF->Close();

  std::vector<TString> mcg;
  if(!purePhoton) {
    mcg.push_back("V#gamma");
    mcg.push_back("Multijets");
  }
  mcg.push_back("#gamma+jets");
  TFile *gF = TFile::Open(gUrl);
  TH1 *hg=0;
  for(size_t ig=0; ig<mcg.size(); ig++)
    {
      TH1 *h=0;
      if(ch=="ll")
	{
	  h=(TH1 *) gF->Get(mcg[ig]+"/ee"+cat+"_"+distr);
	  if(h) h->Add( (TH1 *) gF->Get(mcg[ig]+"/mumu"+cat+"_"+distr) );
	}
      else
	{
	  h=(TH1 *) gF->Get(mcg[ig]+"/"+ch+cat+"_"+distr);
	}
      if(h==0) continue;
      if(hg==0) { hg= h;      hg->SetDirectory(0); }
      else      { hg->Add(h);                      }
    }
  if(hg==0) return 0;
  gF->Close();

  hdy->Divide(hg);
  hdy->SetName(distr+"_ratio");
  return (TH1F *)hdy;
}


//
void closureTest(TString url,TString gUrl,TString distr,TString ch, TString cat, bool purePhoton)
{
  //
  //GET HISTOS FROM FILES
  //
  TString mcdy("Z#rightarrow ll");
  TFile *llF = TFile::Open(url);
  TH1 *hdy = 0;
  if(ch=="ll")
    {
      hdy=(TH1 *) llF->Get(mcdy+"/ee"+cat+"_"+distr);	       
      hdy->Add((TH1 *) llF->Get(mcdy+"/mumu"+cat+"_"+distr) );
      hdy->SetName("ll"+cat+"_"+distr);
    }
  else
    {
      hdy=(TH1 *) llF->Get(mcdy+"/"+ch+cat+"_"+distr);
    }
  if(hdy==0) { llF->Close(); return; }
  hdy->SetDirectory(0);
  llF->Close();

  std::vector<TString> mcg;
  if(!purePhoton) {
    mcg.push_back("V#gamma");
    mcg.push_back("Multijets");
  }
  mcg.push_back("#gamma+jets");
  TFile *gF = TFile::Open(gUrl);
  TH1 *hg=0,*hpureg=0;
  for(size_t ig=0; ig<mcg.size(); ig++)
    {
      TH1 *h=0;
      if(ch=="ll")
	{
	  h=(TH1 *) gF->Get(mcg[ig]+"/mumu"+cat+"_"+distr);
	  if(h){
	    h->Add( (TH1 *) gF->Get(mcg[ig]+"/ee"+cat+"_"+distr) );
	    h->SetName("ll"+cat+"_"+distr);
	  }
	}
      else
	{
	  h=(TH1 *) gF->Get(mcg[ig]+"/"+ch+cat+"_"+distr);
	}
      if(h==0) continue;
      if(hg==0) { hg= h;      hg->SetDirectory(0); }
      else      { hg->Add(h);                      }
      if(mcg[ig]!="#gamma+jets") continue;
      hpureg=(TH1 *) h->Clone(ch+cat+"_"+distr+"_pureg");
      hpureg->SetDirectory(0);
    }
  if(hg==0) return;
  gF->Close();

  if(distr=="qt") { hg->Rebin(4); hdy->Rebin(4); hpureg->Rebin(4); }
  //else if(!distr.Contains("cjv")) { hg->Rebin(4); hdy->Rebin(4); hpureg->Rebin(4); }
  else if(!distr.Contains("cjv")) { hg->Rebin(2); hdy->Rebin(2); hpureg->Rebin(2); }

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
  t1->SetBottomMargin(0);
  
  //find limits
  hdy->Scale(1./hdy->Integral());
  Double_t xmin(hdy->GetXaxis()->GetXmin());
  Double_t xmax(hdy->GetXaxis()->GetXmax());
  if(distr.Contains("met")) {xmin=0;   xmax=250;}
  if(distr.Contains("mt"))  {xmin=120; xmax=600;}
  float ymin(3e-4),ymax(hdy->GetMaximum()*1.2);

  //draw
  float gscale(hg->Integral());
  hg->Scale(1./gscale);
  hg->Draw("hist");   
  hg->SetTitle("All #gamma");
  hg->GetYaxis()->SetRangeUser(ymin,ymax);
  hg->GetXaxis()->SetRangeUser(xmin,xmax); 
  hg->GetYaxis()->SetLabelSize(0.06);
  hg->GetYaxis()->SetTitleSize(0.05);      
  hg->GetYaxis()->SetTitleOffset(1.2);
  hg->GetYaxis()->SetLabelSize(0.05);
  hg->SetLineColor(1);
  hg->SetMarkerColor(1);
  hg->SetMarkerStyle(1);
  hg->SetFillStyle(3001);
  hg->SetFillColor(kGray);
  if(distr.Contains("qt"))
    {
      ymin=3e-6;
      t1->SetLogy();
      t1->SetLogx();
      hg->GetXaxis()->SetRangeUser(50,1e3);
    }
  hg->GetYaxis()->SetRangeUser(ymin,ymax);

  hpureg->Scale(1./gscale);
  hpureg->SetTitle("Prompt #gamma");
  hpureg->GetYaxis()->SetRangeUser(ymin,ymax);
  hpureg->SetLineColor(1);
  hpureg->SetMarkerColor(1);
  hpureg->SetMarkerStyle(1);
  hpureg->SetFillStyle(1001);
  hpureg->SetFillColor(800);
  hpureg->Draw("histsame");

  hdy->SetTitle("Z#rightarrow ll");
  hdy->Draw("e1same");

  TLegend *leg=new TLegend(0.5,0.95,0.95,0.99);
  leg->AddEntry(hg,hg->GetTitle(),"F");
  leg->AddEntry(hpureg,hpureg->GetTitle(),"F");
  leg->AddEntry(hdy,hdy->GetTitle(),"P");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextAlign(11);
  leg->SetTextSize(0.05);
  leg->Draw("same");
  leg->SetNColumns(3);
  
  TPaveText *pave = new TPaveText(0.65,0.7,0.85,0.9,"brNDC");
  pave->SetBorderSize(0);
  pave->SetFillStyle(0);
  pave->SetTextAlign(11);
  pave->SetTextFont(42);
  pave->SetTextSize(0.06);
  TString mjjCat("M_{jj}>1000");
  if(cat.Contains("mjjq016")) mjjCat="M_{jj}<250";
  if(cat.Contains("mjjq033")) mjjCat="250<M_{jj}<350";
  if(cat.Contains("mjjq049")) mjjCat="350<M_{jj}<450";
  if(cat.Contains("mjjq066")) mjjCat="450<M_{jj}<550";
  if(cat.Contains("mjjq083")) mjjCat="550<M_{jj}<750";
  if(cat.Contains("mjjq092")) mjjCat="750<M_{jj}<1000";
  pave->AddText(mjjCat);
  char buf[100];
  sprintf(buf,"#chi^{2}/ndof : %3.2f", hdy->Chi2Test(hg,"WWCHI2/NDF") );
  pave->AddText(buf);
  sprintf(buf,"K-S prob: %3.2f", hdy->KolmogorovTest(hg,"") );
  pave->AddText(buf);
  //  sprintf(buf,"Prob: %3.3f", hdy->Chi2Test(hg,"WW") );
  //pave->AddText(buf);
  pave->Draw();

  pave = new TPaveText(0.1,0.95,0.6,0.99,"NDC");
  pave->SetBorderSize(0);
  pave->SetFillStyle(0);
  pave->SetTextAlign(11);
  pave->SetTextSize(0.04);
  pave->AddText("CMS simulation, #sqrt{s}=8 TeV");
  pave->Draw();
  
  //closure
  c->cd();
  TPad *t2 = new TPad("p2","p2",0,0.0,1.0,0.3);
  t2->SetTopMargin(0);
  t2->SetBottomMargin(0.25);
  t2->Draw();
  t2->cd();

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
  denRelUncH->GetYaxis()->SetRangeUser(0.2,1.74);
  denRelUncH->GetXaxis()->SetTitle(hdy->GetXaxis()->GetTitle());
  denRelUncH->GetXaxis()->SetLabelSize(0.1);
  denRelUncH->GetXaxis()->SetTitleSize(0.12);
  denRelUncH->GetXaxis()->SetTitleOffset(0.8);
  denRelUncH->GetYaxis()->SetLabelSize(0.1);
  denRelUncH->GetYaxis()->SetTitleSize(0.12);
  denRelUncH->GetYaxis()->SetTitleOffset(0.5);
  //gr->GetYaxis()->SetTitle("Pred. rel. bias");
  denRelUncH->GetYaxis()->SetTitle("Ratio");
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

      loose2tightSystUp->Draw("l");
      loose2tightSystDown->Draw("l");

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

  TLine *l = new TLine(xmin,0.5,xmax,0.5);
  l->SetLineColor(kRed);
  l->Draw();
  l = new TLine(xmin,1.5,xmax,1.5);
  l->SetLineColor(kRed);
  l->Draw();
  

  c->Modified();
  c->Update();

  c->SaveAs(ch+cat+"_"+distr+"_closure.png");
  c->SaveAs(ch+cat+"_"+distr+"_closure.pdf");
}
