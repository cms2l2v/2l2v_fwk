#include "TFile.h"
#include "TString.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TSystem.h"
#include "TF1.h"

#include <iostream>

using namespace std;

TString outDir("");

std::pair<float,int> computeChiSquareFor(TGraphErrors *dataGr,TGraphErrors *mcGr)
{
  std::pair<float,int> toRet(0,0);
  if(dataGr==0||mcGr==0) return toRet;
  
  for(int ip=0; ip<dataGr->GetN(); ip++)
    {
      Double_t x,y;
      dataGr->GetPoint(ip,x,y);
      Double_t ymc(-1),yerrmc(-1);
      for(int jp=0; jp<mcGr->GetN(); jp++)
	{
	  Double_t ix,iy;
	  mcGr->GetPoint(jp,ix,iy);
	  if(ix!=x) continue;
	  ymc=iy;
	  yerrmc=mcGr->GetErrorY(jp);
	  break;
	}
      if(ymc<0 || yerrmc<=0) continue;
      toRet.first += pow((y-ymc)/yerrmc,2);
      toRet.second++;
    }

  return toRet;
}


void profileDistributions(TString dist);
void profileUE(TString outUrl="");
void drawCMSHeader(TString ch,float txtSize=0)
{
  TPaveText *pt=new TPaveText(0.1,0.95,0.9,0.99,"brNDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->SetTextAlign(12);
  if(txtSize!=0) pt->SetTextSize(txtSize);
  TString label("CMS preliminary, #sqrt{s}=8 TeV, #scale[0.5]{#int}L=19.7 fb^{-1}, " + ch + " events");
  label.ReplaceAll("mu","#mu");
  pt->AddText(label);
  pt->Draw();
}

//
void profileUE(TString outUrl)
{
  //prepare output
  if(outUrl!=""){
    gSystem->ExpandPathName(outUrl);
    if(!outUrl.EndsWith("/")) outUrl = outUrl +"/";
    gSystem->Exec("mkdir -p " + outUrl);
    outDir=outUrl;
  }

  profileDistributions("nchprofpt");
  profileDistributions("ptfluxprofpt");
  profileDistributions("avgptfluxprofpt");
  profileDistributions("nchprofmt");
  profileDistributions("ptfluxprofmt");
  profileDistributions("avgptfluxprofmt");
  profileDistributions("ptfluxprofphi");
}

//
void profileDistributions(TString dist)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(1);

  TFile *inF=TFile::Open("~/work/top_539/plotter.root");
  TFile *inSystF=TFile::Open("~/work/top_539/plotter_syst.root");

  bool isPhi(dist.Contains("phi"));
  TString profs[]={dist,
		   dist+"toward",
		   dist+"transverse",
		   dist+"away"};
  const size_t nprofs= isPhi ? 1 : sizeof(profs)/sizeof(TString);


  Int_t colors[]={1,kBlue,kRed,kGreen-3};

  TString cats[]={"inclusive","away","transverse","toward"};

  /*
  TFile *inNomF=inF;
  TString nomTTbar="t#bar{t}";
  TString nomTTbarTitle="MG+PY";
  
  //TString systList[]={"t#bar{t}systmcatnlo"};
  //TString systLabels[]={"MC@NLO+HW"};
  
  //TString systList[]={"t#bar{t}systq2down","t#bar{t}systq2up"};
  //TString systLabels[]={"-Q^{2}","+Q^{2}"};

  //TString systList[]={"t#bar{t}systmepsdown","t#bar{t}systmepsup"};
  // TString systLabels[]={"-ME-PS","+ME-PS"};
  */

  TFile *inNomF=inSystF;
  TString nomTTbar="t#bar{t}systtunep11";
  TString nomTTbarTitle="P11";
  TString systList[]={"t#bar{t}systtunep11nocr","t#bar{t}systtunep11tev"};
  TString systLabels[]={"P11-noCR","P11TeV"};


  Int_t systColors[]={1,kAzure,kRed-6,kGreen+3};

  const size_t nSysts=sizeof(systList)/sizeof(TString);

  TString ch[]={"emu"};
  const size_t nch=sizeof(ch)/sizeof(TString);
  
  for(size_t ich=0; ich<nch; ich++)
    {
      TCanvas *c=new TCanvas("c"+ch[ich]+dist,"c"+ch[ich]+dist,800,800); c->SetTopMargin(0); c->SetBottomMargin(0);
      c->Divide(1,nprofs);
      
      TCanvas *cratios=new TCanvas("cratios"+ch[ich]+dist,"cratios"+ch[ich]+dist,800,800); cratios->SetTopMargin(0); cratios->SetBottomMargin(0);
      cratios->Divide(1,nprofs);

      TCanvas *cproj=new TCanvas("cproj"+ch[ich]+dist,"cproj"+ch[ich]+dist,600,600);
      TLegend *dataprojLeg=new TLegend(0.5,0.9,0.6,0.7,"Data","brNDC");
      TLegend *mcprojLeg=new TLegend(0.7,0.9,0.92,0.7,"MC","brNDC");
      
      for(size_t i=0; i<nprofs; i++)
	{
	  TH2 *bckgMC   = (TH2 *) inF->Get("Z#rightarrow ll/"+ch[ich]+"_"+profs[i]); 
	  bckgMC->Add( (TH2 *) inF->Get("VV/"+ch[ich]+"_"+profs[i]) );
	  bckgMC->Add( (TH2 *) inF->Get("Single top/"+ch[ich]+"_"+profs[i]) );
	  bckgMC->Add( (TH2 *) inF->Get("W+jets-multijets/"+ch[ich]+"_"+profs[i]) );
	  
	  //build the total MC with signal alternatives
	  TH2 *MC = (TH2 *) inNomF->Get(nomTTbar+"/"+ch[ich]+"_"+profs[i]) ; 
	  Double_t totalTTbar(MC->Integral());
	  MC->Add(bckgMC);
	  MC->SetDirectory(0);
	  MC->Sumw2();	  
	  
	  std::vector<TH2 *> systMC;
	  for(size_t isyst=0; isyst<nSysts; isyst++)
	    {
	      TH2F *h=(TH2F *)inSystF->Get(systList[isyst]+"/"+ch[ich]+"_"+profs[i]) ;
	      h->Scale(totalTTbar/h->Integral());
	      h->Add(bckgMC);
	      h->Sumw2();
	      h->SetDirectory(0);
	      systMC.push_back(h);
	    }
		

	  //get the data
	  TH2 *Data = (TH2 *) inF->Get("data/"+ch[ich]+"_"+profs[i]);            
	  Data->SetDirectory(0);
	  Data->Sumw2();

	  TGraphErrors *MCProf   = new TGraphErrors(MC->ProfileX());    MCProf->SetMarkerStyle(24);   MCProf->SetFillStyle(0); MCProf->SetName(ch[ich]+profs[i]+"mc");
	  TGraphErrors *DataProf = new TGraphErrors(Data->ProfileX());  DataProf->SetMarkerStyle(20); DataProf->SetFillStyle(0); DataProf->SetName(ch[ich]+profs[i]+"data");
	  
	  //build data/MC scale factors
	  std::vector<TGraphErrors *> data2mcProfs;
	  for(size_t isyst=0; isyst<=systMC.size(); isyst++)
	    {
	      TGraphErrors *prof= (isyst==0 ? MCProf : new TGraphErrors(systMC[isyst-1]->ProfileX()));
	      TString baseName(ch[ich]+profs[i]);
	      if(isyst) baseName += systList[isyst];
	      prof = (TGraphErrors *) prof->Clone(baseName+"data2mc");
	      for(int ip=0; ip<DataProf->GetN(); ip++)
		{
		  Double_t x,y,ydata,y_err,ydata_err;
		  prof->GetPoint(ip,x,y);         y_err=prof->GetErrorY(ip);
		  DataProf->GetPoint(ip,x,ydata); ydata_err=DataProf->GetErrorY(ip);
		  if(y<=0) continue;
		  prof->SetPoint(ip,x,ydata/y);
		  prof->SetPointError(ip,0,sqrt(pow(ydata*y_err,2)+pow(ydata_err*y,2))/pow(y,2));
		}
	      prof->SetFillColor(systColors[isyst]);
	      prof->SetFillStyle(3001+isyst%2);
	      prof->SetTitle( isyst==0 ? nomTTbarTitle : systLabels[isyst-1] );
	      //prof->SetMarkerStyle(24);   prof->SetFillStyle(3001); prof->SetMarkerColor(1+isyst); prof->SetMarkerColor(1+isyst); prof->SetLineColor(1+isyst);
	      data2mcProfs.push_back(prof);
	    }

	  
	  TH1D *MCProjY=MC->ProjectionY();
	  MCProjY->Scale(1./MCProjY->Integral()); 
	  TGraphErrors *MCProj   = new TGraphErrors(MCProjY);    MCProj->SetMarkerStyle(24); MCProj->SetFillStyle(0);  MCProj->SetName(ch[ich]+profs[i]+"projmc");
	  TH1D *DataProjY=Data->ProjectionY();
	  DataProjY->Scale(1./DataProjY->Integral());
	  TGraphErrors *DataProj = new TGraphErrors(DataProjY);  DataProj->SetMarkerStyle(20); DataProj->SetFillStyle(0); DataProj->SetName(ch[ich]+profs[i]+"projdata");
	  MCProj->SetLineColor(colors[i]);   MCProj->SetMarkerColor(colors[i]);   MCProj->SetFillColor(colors[i]); MCProj->SetFillStyle(1001);
	  DataProj->SetLineColor(colors[i]); DataProj->SetMarkerColor(colors[i]);   DataProj->SetFillColor(colors[i]);
	  
	  TPad *p=(TPad *)cproj->cd();
	  p->SetLeftMargin(0.15);
	  p->SetRightMargin(0.02);
	  p->SetTopMargin(0.05);
	  p->SetLogy();
	  MCProj->SetFillStyle(0);
	  MCProj->Draw(i==0 ? "al3" : "l3");
	  MCProj->GetYaxis()->SetRangeUser(1e-5,1.0);
	  MCProj->GetXaxis()->SetTitle( MC->GetYaxis()->GetTitle() );
	  MCProj->GetYaxis()->SetTitle( "1/N dN/dx" );
	  MCProj->GetYaxis()->SetTitleOffset(1.8);
	  DataProj->Draw("p");
	  std::pair<float,int> chi2=computeChiSquareFor(DataProj,MCProj);
	  char buf[200];
	  sprintf(buf,"#scale[0.7]{#chi^{2}/ndof=%3.1f}", chi2.first/chi2.second );
	  dataprojLeg->AddEntry(DataProj,buf,"p");
	  mcprojLeg->AddEntry(MCProj,cats[i],"l");
	  if(i==0) drawCMSHeader(ch[ich]);
	  

	  p=(TPad *)c->cd(i+1);
	  if(i<nprofs-1) p->SetBottomMargin(0.01);
	  if(i>0) p->SetTopMargin(0);
	  if(i==0)p->SetTopMargin(0.1);
	  if(i==nprofs-1) p->SetBottomMargin(0.15);
	  TGraphErrors *frame=DataProf;
	  frame->Draw("ap");
	  //frame->GetYaxis()->SetRangeUser(0.54,1.46);
	  frame->GetXaxis()->SetTitle(MC->GetXaxis()->GetTitle());
	  //frame->GetXaxis()->SetRangeUser(0,4.75);
	  TString yTit("<"); yTit+=MC->GetYaxis()->GetTitle(); yTit +=">";
	  frame->GetYaxis()->SetTitle(yTit);
	  frame->GetYaxis()->SetLabelSize(0.07);
	  frame->GetYaxis()->SetTitleSize(0.09);
	  frame->GetYaxis()->SetTitleOffset(0.5);
	  frame->GetXaxis()->SetLabelSize(0.07);
	  frame->GetXaxis()->SetTitleSize(0.09);
	  frame->GetXaxis()->SetTitleOffset(0.7);
	  //p->SetGridy();
	  //DataProf->Draw("p");
	  MCProf->Draw("p");
	  if(i==0)
	    {
	      drawCMSHeader(ch[ich],0.08);

	      TLegend *leg=new TLegend(0.6,0.95,1.0,0.99);
	      leg->SetBorderSize(0);
	      leg->SetFillStyle(0);
	      leg->SetTextFont(42);
	      leg->SetTextSize(0.09);
	      leg->AddEntry(DataProf,"data","p");
	      leg->AddEntry(MCProf,"simulation","p");
	      leg->SetNColumns(2);
	      leg->Draw();
	    }

	  TPaveText *pt=new TPaveText(0.15,0.5,0.8,0.85,"brNDC");
	  pt->SetBorderSize(0);
	  pt->SetFillStyle(0);
	  pt->SetTextAlign(13);
	  pt->SetTextFont(42);
	  pt->SetTextColor(kBlue);
	  pt->SetTextSize(0.08);
	  pt->AddText("[ "+cats[i]+" ]");
	  if(i==0)
	    {
	      pt->AddText("p_{T}>0.5 GeV");
	      pt->AddText("|#eta|<2.1");
	    }
	  pt->Draw();


	  p=(TPad *) cratios->cd(i+1);
	  if(i<nprofs-1) p->SetBottomMargin(0.01);
	  if(i>0) p->SetTopMargin(0);
	  if(i==0)p->SetTopMargin(0.1);
	  if(i==nprofs-1) p->SetBottomMargin(0.15);
	  frame=data2mcProfs[0];
	  frame->Draw("a3");
	  TLine *l=new TLine(frame->GetXaxis()->GetXmin(),1,frame->GetXaxis()->GetXmax(),1);
	  l->Draw();
	  frame->GetYaxis()->SetRangeUser(0.54,1.46);
	  frame->GetXaxis()->SetTitle(MC->GetXaxis()->GetTitle());
	  frame->GetYaxis()->SetTitle("Data/Simulation");
	  frame->GetYaxis()->SetLabelSize(0.07);
	  frame->GetYaxis()->SetTitleSize(0.09);
	  frame->GetYaxis()->SetTitleOffset(0.5);
	  frame->GetXaxis()->SetLabelSize(0.07);
	  frame->GetXaxis()->SetTitleSize(0.09);
	  frame->GetXaxis()->SetTitleOffset(0.7);
	  //p->SetGridy();
	  for(size_t ip=1; ip<data2mcProfs.size(); ip++) data2mcProfs[ip]->Draw("3");
	  if(i==0)
	    {
	      drawCMSHeader(ch[ich],0.08);

	      TLegend *leg=new TLegend(0.6,0.94,1.0,0.98);
	      leg->SetBorderSize(0);
	      leg->SetFillStyle(0);
	      leg->SetTextFont(42);
	      leg->SetTextSize(0.09);
	      leg->SetNColumns(data2mcProfs.size());
	      for(size_t ip=0; ip<data2mcProfs.size(); ip++)
		leg->AddEntry(data2mcProfs[ip],data2mcProfs[ip]->GetTitle(),"f");
	      leg->Draw();
	    }
	  
	  pt=new TPaveText(0.12,0.65,0.8,0.9,"brNDC");
	  pt->SetBorderSize(0);
	  pt->SetFillStyle(0);
	  pt->SetTextAlign(13);
	  pt->SetTextFont(42);
	  pt->SetTextColor(kBlue);
	  pt->SetTextSize(0.08);
	  pt->AddText("[ "+cats[i]+" ]");
	  if(i==0)
	    {
	      pt->AddText("p_{T}>0.5 GeV");
	      pt->AddText("|#eta|<2.1");
	    }
	  pt->Draw();
	}

      cproj->cd();
      dataprojLeg->SetFillStyle(0);
      dataprojLeg->SetBorderSize(0);
      dataprojLeg->SetTextFont(42);
      dataprojLeg->SetTextSize(0.03);
      dataprojLeg->Draw();
      mcprojLeg->SetFillStyle(0);
      mcprojLeg->SetBorderSize(0);
      mcprojLeg->SetTextFont(42);
      mcprojLeg->SetTextSize(0.03);
      mcprojLeg->Draw();


      TString pf(".png");
      c->SaveAs(outDir+c->GetName()+pf);
      cproj->SaveAs(outDir+cproj->GetName()+pf);
      cratios->SaveAs(outDir+cratios->GetName()+pf);

    }

  inF->Close();
  inSystF->Close();
}

void showThrustProfile(TString est="raw")
{


  TFile *inF=TFile::Open("~/work/top_539/plotter_syst.root");
  TH2 *ptResponse=(TH2 *)inF->Get("t#bar{t}172.5/emu_thrustptresponse_"+est);
  ptResponse->SetDirectory(0);
  TH2 *phiResponse=(TH2 *)inF->Get("t#bar{t}172.5/emu_thrustphiresponse_"+est);
  phiResponse->SetDirectory(0);
  inF->Close();
  
  gStyle->SetPalette(55);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TCanvas *c=new TCanvas("c","c",1000,500);
  c->Divide(2,1);
  c->cd(1);
  ptResponse->Draw("colz");
  ptResponse->GetYaxis()->SetTitle("#Delta p_{T} [GeV]");
  drawCMSHeader("emu");
  c->cd(2);
  phiResponse->Draw("colz");

  TCanvas *cprof=new TCanvas("cprof","cprof",500,500);
  TGraphErrors *gr=new TGraphErrors; gr->SetMarkerStyle(20); gr->SetName("pt");
  for(int xbin=1; xbin<=ptResponse->GetXaxis()->GetNbins(); xbin++)
    {
      TH1 *py=ptResponse->ProjectionY("py",xbin,xbin);
      Int_t ip=gr->GetN();
      float x=ptResponse->GetXaxis()->GetBinCenter(xbin);
      float ex=ptResponse->GetXaxis()->GetBinWidth(xbin)/2;
      gr->SetPoint(ip,x,py->GetMean());
      gr->SetPointError(ip,ex,py->GetRMS());
      py->Delete();
    }
  gr->Draw("ap");
  gr->GetYaxis()->SetTitle("#Delta p_{T} [GeV]");
  gr->GetXaxis()->SetTitle(ptResponse->GetXaxis()->GetTitle());
  drawCMSHeader("emu");
}
