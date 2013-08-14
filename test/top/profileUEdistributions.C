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
#include "TMath.h"

#include <iostream>
#include <vector>

using namespace std;

struct ComparisonScheme_t
{
  TString title;
  bool nomFromSyst;
  TString              nom,      nomLabel;
  std::vector<TString> systList, systLabels;
  std::vector<int>     systColors;
};

typedef std::vector<ComparisonScheme_t> ComparisonSchemeCollection_t;

TString outDir="539/ueprofiles";
TString stdPlotterUrl="539/plotter_ue.root";
TString systPlotterUrl="539/plotter_syst_ue.root";
ComparisonSchemeCollection_t m_compSchemes;


std::pair<float,int> computeChiSquareFor(TGraphErrors *dataGr,TGraphErrors *mcGr);
void drawCMSHeader(TString ch="",float txtSize=0, bool isSim=false);
void profile(TString dist);
void profileUEdistributions(TString type="mgpy");
TGraphErrors *getEvolutionWithErrors(TH1 *h);
TGraphErrors *getRatio(TGraphErrors *cdfgr,TGraphErrors *datacdf);

//
TGraphErrors *getEvolutionWithErrors(TH1 *h)
{
  TGraphErrors *gr=new TGraphErrors;
  gr->SetName(h->GetName()+TString("_evol"));
 
  TString htype(h->ClassName());
  Int_t nbins=h->GetXaxis()->GetNbins();
  for(Int_t ibin=1; ibin<=nbins; ibin++)
    {
      Double_t x=h->GetXaxis()->GetBinCenter(ibin);
      Double_t xerr=h->GetXaxis()->GetBinWidth(ibin)/2;
      Double_t val(0),val_err(0);
	  
      //for TH1 we compute the gap fraction
      if(htype.Contains("TH1"))
	{
	  Double_t total=h->Integral(1,nbins);
	  val=h->IntegralAndError(1,ibin,val_err)/total;
	  val_err/=total;
	}
    
      //for TH2 we project the average y as function of x
      if(htype.Contains("TH2"))
	{
	  TH2* h2d=(TH2 *)h;
	  TH1 *py=h2d->ProjectionY("py",ibin,ibin);
	  val=py->GetMean();
	  val_err=py->GetMeanError();
	  delete py;
	}

      if(val==0) continue;
      Int_t ip=gr->GetN();
      gr->SetPoint(ip,x,val);
      gr->SetPointError(ip,xerr,val_err);
    }

  return gr;
}

//
TGraphErrors *getRatio(TGraphErrors *ngr,TGraphErrors *dgr)
{
  TGraphErrors *gr=new TGraphErrors;
  gr->SetName( ngr->GetName()+TString("_")+dgr->GetName() );
  for(int ip=0; ip<dgr->GetN(); ip++)
    {
      Double_t x,xerr,yn,ynerr,yd,yderr;
      ngr->GetPoint(ip,x,yn); xerr=ngr->GetErrorX(ip); ynerr=ngr->GetErrorY(ip);
      dgr->GetPoint(ip,x,yd);                          yderr=dgr->GetErrorY(ip); 
      
      gr->SetPoint(ip,x,yd==0 ? 1 : yn/yd);
      gr->SetPointError(ip,xerr,yd==0 ? 1 : sqrt(pow(ynerr*yd,2)+pow(yderr*yn,2))/(yd*yd));
    }
  return gr;
}



//
std::pair<float,int> computeChiSquareFor(TGraphErrors *dataGr,TGraphErrors *mcGr)
{
  std::pair<float,int> toRet(0,0);
  if(dataGr==0||mcGr==0) return toRet;
  
  for(int ip=0; ip<dataGr->GetN(); ip++)
    {
      Double_t x,y,yerrdata(-1);
      dataGr->GetPoint(ip,x,y);
      yerrdata=dataGr->GetErrorY(ip);

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
      Double_t yerr=sqrt(pow(yerrmc,2)+pow(yerrdata,2));
      toRet.first += TMath::Power((y-ymc)/yerr,2);
      toRet.second++;
    }

  return toRet;
}


//
void drawCMSHeader(TString ch,float txtSize,bool isSim)
{
  TPaveText *pt=new TPaveText(0.05,0.93,0.9,0.97,"brNDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->SetTextAlign(12);
  if(txtSize!=0) pt->SetTextSize(txtSize);
  TString label("CMS preliminary, #sqrt{s}=8 TeV, #scale[0.5]{#int}L=19.7 fb^{-1}");
  if(isSim) label="CMS simulation, #sqrt{s}=8 TeV";
  pt->AddText(label);
  pt->Draw();
}

//
void profile(TString dist)
{
  bool isSimple(dist.Contains("phi") || dist.Contains("nvtx"));
  bool isPhi(dist.Contains("phi"));

  TString ch[]={"emu","emueq0j","emueq1j","emugeq2j",
		"ll", "lleq0j", "lleq1j", "llgeq2j"};
  const size_t nch=sizeof(ch)/sizeof(TString);
  
  TString profs[]={dist, dist+"away", dist+"transverse", dist+"toward"};
  const size_t nprofs= isSimple ? 1 : sizeof(profs)/sizeof(TString);

  TString cats[]={"inclusive","away","transverse","toward"};
  Int_t colors[]={1,kBlue,kRed,kGreen-3};

  //run all comparison schemes possible
  for(ComparisonSchemeCollection_t::iterator it=m_compSchemes.begin(); it != m_compSchemes.end(); it++)
    {  
      ComparisonScheme_t &comp=*it;
      cout << "[ " << comp.title << " ]" << endl;

      //prepare the output
      TString localOutDir=outDir+"/"+comp.title+"/";
      gSystem->ExpandPathName(localOutDir);
      gSystem->Exec("mkdir -p " + localOutDir);
      
      //retrieve the plots
      TFile *inF     = TFile::Open(stdPlotterUrl);
      TFile *inSystF = TFile::Open(systPlotterUrl);
      TFile *inNomF  = comp.nomFromSyst ?  inSystF : inF;
      for(size_t ich=0; ich<nch; ich++)
	{
	  //if(ch[ich].Contains("ll") && it!=m_compSchemes.begin()) continue;
	  if(isSimple){
	    if(ch[ich].Contains("eq0j")) cats[0]="=0 jets";
	    if(ch[ich].Contains("eq1j")) cats[0]="=1 jet";
	    if(ch[ich].Contains("geq2j")) cats[0]="#geq2 jets";
	  }

	  TCanvas *c=isSimple?
	    new TCanvas("c"+ch[ich]+dist,"c"+ch[ich]+dist,600,600):
	    new TCanvas("c"+ch[ich]+dist,"c"+ch[ich]+dist,800,1200);
	  c->SetTopMargin(0); 
	  c->SetBottomMargin(0);
	  c->Divide(1,nprofs);
	  float yprofmin(0),yprofmax(40);
	  
	  TCanvas *cratios=isSimple ? 
	    new TCanvas("cratios"+ch[ich]+dist,"cratios"+ch[ich]+dist,600,600):
	    new TCanvas("cratios"+ch[ich]+dist,"cratios"+ch[ich]+dist,600,600); 
	  cratios->SetTopMargin(0); cratios->SetBottomMargin(0);
	  cratios->Divide(1,nprofs);
	  
	  TCanvas *cproj= isSimple? 
	    new TCanvas("cproj"+ch[ich]+dist,"cproj"+ch[ich]+dist,600,600) : 
	    new TCanvas("cproj"+ch[ich]+dist,"cproj"+ch[ich]+dist,600,600);
	  TLegend *dataprojLeg=new TLegend(0.5,0.9,0.6,0.7,"Data","brNDC");
	  TLegend *mcprojLeg=new TLegend(0.7,0.9,0.92,0.7,"MC","brNDC");
	  
	  for(size_t i=0; i<nprofs; i++)
	    {
	      TH2 *bckgMC = (TH2 *) inF->Get("Single top/"      + ch[ich]+"_"+profs[i]); 
	      bckgMC->Add(  (TH2 *) inF->Get("VV/"              + ch[ich]+"_"+profs[i]) );
	      bckgMC->Add(  (TH2 *) inF->Get("Z#rightarrow ll/" + ch[ich]+"_"+profs[i]) );
	      bckgMC->Add(  (TH2 *) inF->Get("W,multijets/"     + ch[ich]+"_"+profs[i]) );
	      bckgMC->Add(  (TH2 *) inF->Get("other t#bar{t}/"  + ch[ich]+"_"+profs[i]) );
	     
	      //build the total MC with signal alternatives
	      TH2 *MC = (TH2 *) inNomF->Get(comp.nom+"/"   +ch[ich]+"_"+profs[i]) ; 
	      Double_t totalTTbar(MC->Integral());
	      MC->Add(bckgMC);
	      MC->SetDirectory(0);
	      MC->Sumw2();	  

	      //alternative MCs
	      std::vector<TH2 *> systMC;
	      for(size_t isyst=0; isyst<comp.systList.size(); isyst++)
		{
		  TH2F *h=(TH2F *)inSystF->Get(comp.systList[isyst]+"/"+ch[ich]+"_"+profs[i]) ;
		  h->Scale(totalTTbar/h->Integral());
		  h->Add(bckgMC);
		  h->Sumw2();
		  h->SetDirectory(0);
		  systMC.push_back(h);
		} 
	     
	      //data
	      TH2 *Data = (TH2 *) inF->Get("data/"+ch[ich]+"_"+profs[i]);            
	      Data->SetDirectory(0);
	      Data->Sumw2();
	      
	      TGraphErrors *MCProf   = new TGraphErrors(MC->ProfileX());    MCProf->SetMarkerStyle(24);   MCProf->SetFillStyle(0);   MCProf->SetName(ch[ich]+profs[i]+"mc");
	      TGraphErrors *DataProf = new TGraphErrors(Data->ProfileX());  DataProf->SetMarkerStyle(20); DataProf->SetFillStyle(0); DataProf->SetName(ch[ich]+profs[i]+"data");
	      
	      //build data/MC scale factors
	      std::vector<TGraphErrors *> data2mcProfs;
	      for(size_t isyst=0; isyst<=systMC.size(); isyst++)
		{
		  TGraphErrors *prof= (isyst==0 ? MCProf : new TGraphErrors(systMC[isyst-1]->ProfileX()));
		  TString baseName(ch[ich]+profs[i]);
		  if(isyst>0) baseName += comp.systList[isyst];
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
		  prof->SetFillStyle(3001+isyst%2);
		  prof->SetFillColor(isyst==0 ? 1 : comp.systColors[isyst-1]);
		  prof->SetTitle( isyst==0 ? comp.nomLabel : comp.systLabels[isyst-1] );
		  data2mcProfs.push_back(prof);
		}
	      
	      //TH1D *MCProjY=MC->ProjectionY();
	      TH1D *MCProjY=0;
	      for(int xbin=1; xbin<=MC->GetXaxis()->GetNbins(); xbin++){
		TH1D *py=MC->ProjectionY("imcpy",xbin,xbin);
		if(MCProjY==0) { MCProjY=(TH1D *)py->Clone("mcpy"); MCProjY->Reset("ICE"); }
		MCProjY->Add( py, MC->GetXaxis()->GetBinWidth(xbin) );
	      }
	      MCProjY->Scale(1./MCProjY->Integral()); 
	      TGraphErrors *MCProj   = new TGraphErrors(MCProjY);    MCProj->SetMarkerStyle(24); MCProj->SetFillStyle(0);  MCProj->SetName(ch[ich]+profs[i]+"projmc");
	      //TH1D *DataProjY=Data->ProjectionY();
	      TH1D *DataProjY=0;
	      for(int xbin=1; xbin<=Data->GetXaxis()->GetNbins(); xbin++){
		TH1D *py=Data->ProjectionY("imcpy",xbin,xbin);
		if(DataProjY==0) { DataProjY=(TH1D *)py->Clone("mcpy"); DataProjY->Reset("ICE"); }
		DataProjY->Add( py, Data->GetXaxis()->GetBinWidth(xbin) );
	      }
	      DataProjY->Scale(1./DataProjY->Integral());
	      TGraphErrors *DataProj = new TGraphErrors(DataProjY);  DataProj->SetMarkerStyle(20); DataProj->SetFillStyle(0); DataProj->SetName(ch[ich]+profs[i]+"projdata");
	      MCProj->SetLineColor(colors[i]);   MCProj->SetMarkerColor(colors[i]);   MCProj->SetFillColor(colors[i]); MCProj->SetFillStyle(1001);
	      DataProj->SetLineColor(colors[i]); DataProj->SetMarkerColor(colors[i]);   DataProj->SetFillColor(colors[i]);

	      TPad *p=(TPad *)cproj->cd();
	      p->SetLeftMargin(0.15);
	      p->SetRightMargin(0.05);
	      p->SetTopMargin(0.1);
	      p->SetLogy();
	      MCProj->SetFillStyle(0);
	      MCProj->Draw(i==0 ? "al3" : "l3");
	      MCProj->GetYaxis()->SetRangeUser(1e-5,1.0);
	      MCProj->GetXaxis()->SetTitle( MC->GetYaxis()->GetTitle() );
	      MCProj->GetYaxis()->SetTitle( TString("#scale[0.5]{#int}d#Delta p_{t#bar{t}}") + MC->GetZaxis()->GetTitle() );
	      MCProj->GetYaxis()->SetTitleOffset(1.8);
	      DataProj->Draw("p");
	      std::pair<float,int> chi2=computeChiSquareFor(DataProj,MCProj);
	      char buf[200];
	      sprintf(buf,"#scale[0.7]{#chi^{2}/ndof=%3.1f}", chi2.first/chi2.second );
	      dataprojLeg->AddEntry(DataProj,buf,"p");
	      mcprojLeg->AddEntry(MCProj,cats[i],"l");
	      if(i==0) {
		drawCMSHeader(ch[ich],0.04);
	      }    

	      p=(TPad *)c->cd(i+1);
	      if(dist.EndsWith("pt")) p->SetLogx();
	      if(i<nprofs-1) p->SetBottomMargin(0.01);
	      if(i>0) p->SetTopMargin(0);
	      if(i==0)p->SetTopMargin(0.1);
	      if(i==nprofs-1) p->SetBottomMargin(0.15);
	      p->SetRightMargin(0.05);
	      TGraphErrors *frame=DataProf;
	      frame->Draw("ap");
	      if(i==0) {
		yprofmin=0;//frame->GetYaxis()->GetXmin()*0.95;
		yprofmax=frame->GetYaxis()->GetXmax()*1.2;
	      }
	      frame->GetYaxis()->SetRangeUser(yprofmin,yprofmax);
	      frame->GetXaxis()->SetTitle(MC->GetXaxis()->GetTitle());
	      TString yTit("<"); yTit+=MC->GetYaxis()->GetTitle(); yTit +=">";
	      frame->GetYaxis()->SetTitle(yTit);
	      frame->GetYaxis()->SetLabelSize(isSimple ? 0.04 : 0.07);
	      frame->GetYaxis()->SetTitleSize(isSimple ? 0.05 : 0.09);
	      frame->GetYaxis()->SetTitleOffset(0.5);
	      frame->GetXaxis()->SetLabelSize(isSimple ? 0.04 : 0.07);
	      frame->GetXaxis()->SetTitleSize(isSimple ? 0.05 : 0.09);
	      frame->GetXaxis()->SetTitleOffset(0.7);
	      if(!isSimple){
		if(i<nprofs-1)
		  {
		    frame->GetXaxis()->SetLabelSize(0);
		    frame->GetXaxis()->SetTitleSize(0);
		  }
		else
		  {
		    frame->GetXaxis()->SetMoreLogLabels();
		  }
	      }
	      else{
		p->SetBottomMargin(0.15);
		frame->GetYaxis()->SetTitleOffset(1.0);
	      }
	      MCProf->Draw("p");
	      if(i==0)
		{
		  drawCMSHeader(ch[ich],isSimple ? 0.04 : 0.08);

		  if(isPhi){
		    TLine *l=new TLine(60,0,60,yprofmax); l->SetLineColor(15); l->Draw();
		    l=new TLine(120,0,120,yprofmax); l->SetLineColor(15); l->Draw();
		  }	  

		  TLegend *leg= isSimple ?
		    new TLegend(0.12,0.85,0.8,0.90):
		    new TLegend(0.6,0.95,1.0,0.99);
		  leg->SetBorderSize(0);
		  leg->SetFillStyle(0);
		  leg->SetTextFont(42);
		  leg->SetTextSize(isSimple ? 0.05 : 0.09);
		  leg->AddEntry(DataProf,"data","p");
		  leg->AddEntry(MCProf,"simulation","p");
		  leg->SetNColumns(2);
		  leg->Draw();
		}
	      
	      TPaveText *pt=isSimple ? 
		new TPaveText(0.15,0.75,0.8,0.85,"brNDC"):
		new TPaveText(0.15,0.5,0.8,0.85,"brNDC");
	      pt->SetBorderSize(0);
	      pt->SetFillStyle(0);
	      pt->SetTextAlign(13);
	      pt->SetTextFont(42);
	      pt->SetTextColor(kBlue);
	      pt->SetTextSize(isSimple ? 0.05 : 0.08);
	      pt->AddText("[ "+cats[i]+" ]");
	      if(i==0 && !isSimple)
		{
		  pt->AddText("p_{T}>0.5 GeV");
		  pt->AddText("|#eta|<2.1");
		}
	      pt->Draw();
	      
	      
	      p=(TPad *) cratios->cd(i+1);
	      if(dist.EndsWith("pt")) p->SetLogx();
	      if(i<nprofs-1) p->SetBottomMargin(0.01);
	      if(i>0) p->SetTopMargin(0);
	      if(i==0)p->SetTopMargin(0.1);
	      if(i==nprofs-1) p->SetBottomMargin(0.15);
	      p->SetRightMargin(0.05);
	      frame=data2mcProfs[0];
	      frame->Draw("a3");
	      TLine *l=new TLine(frame->GetXaxis()->GetXmin(),1,frame->GetXaxis()->GetXmax(),1);
	      l->Draw();
	      frame->GetYaxis()->SetRangeUser(0.54,1.46);
	      frame->GetXaxis()->SetTitle(MC->GetXaxis()->GetTitle());
	      TString yTitVar=MC->GetYaxis()->GetTitle();
	      frame->GetYaxis()->SetTitle( "#Kappa_{data/sim}^{"+yTitVar+"}");
	      frame->GetYaxis()->SetLabelSize(isSimple ? 0.04 : 0.07);
	      frame->GetYaxis()->SetTitleSize(isSimple ? 0.05 : 0.09);
	      frame->GetYaxis()->SetTitleOffset(0.5);
	      frame->GetXaxis()->SetLabelSize(isSimple ? 0.04 : 0.07);
	      frame->GetXaxis()->SetTitleSize(isSimple ? 0.05 : 0.09);
	      frame->GetXaxis()->SetTitleOffset(0.7);
	      if(!isSimple){
		if(i<nprofs-1)
		  {
		    frame->GetXaxis()->SetLabelSize(0);
		    frame->GetXaxis()->SetTitleSize(0);
		  }
		else
		  {
		    frame->GetXaxis()->SetMoreLogLabels();
		  }
	      }
	      else{
		p->SetBottomMargin(0.15);
		p->SetLeftMargin(0.15);
		frame->GetYaxis()->SetTitleOffset(1.1);
	      }

	      for(size_t ip=1; ip<data2mcProfs.size(); ip++) data2mcProfs[ip]->Draw("3");
	      if(i==0)
		{
		  drawCMSHeader(ch[ich],isSimple ? 0.05 : 0.08);
		  
		  TLegend *leg= isSimple ?
		    new TLegend(0.16,0.85,0.8,0.90):
		    new TLegend(0.6,0.94,1.0,0.98);

		  leg->SetBorderSize(0);
		  leg->SetFillStyle(0);
		  leg->SetTextFont(42);
		  leg->SetTextSize(isSimple ? 0.05 : 0.09);
		  leg->SetNColumns(data2mcProfs.size());
		  for(size_t ip=0; ip<data2mcProfs.size(); ip++)
		    leg->AddEntry(data2mcProfs[ip],data2mcProfs[ip]->GetTitle(),"f");
		  leg->Draw();
		}
	      
	      pt=isSimple?
		new TPaveText(0.16,0.75,0.8,0.85,"brNDC"):
		new TPaveText(0.12,0.65,0.8,0.9,"brNDC");
	      pt->SetBorderSize(0);
	      pt->SetFillStyle(0);
	      pt->SetTextAlign(13);
	      pt->SetTextFont(42);
	      pt->SetTextColor(kBlue);
	      pt->SetTextSize(isSimple ? 0.05 : 0.08);
	      pt->AddText("[ "+cats[i]+" ]");
	      if(i==0 && !isSimple)
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
	 	 
	  TString pf[]={".png",".pdf",".C"};
	  for(size_t ipf=0; ipf<sizeof(pf)/sizeof(TString); ipf++){
	    c->SaveAs(localOutDir+c->GetName()+pf[ipf]);
	    cproj->SaveAs(localOutDir+cproj->GetName()+pf[ipf]);
	    cratios->SaveAs(localOutDir+cratios->GetName()+pf[ipf]);
	  }
	}

      inF->Close();
      inSystF->Close();
    }
}


//
void profileUEdistributions(TString type)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(1);
   
  m_compSchemes.clear();

  //define the comparisons to perform
  ComparisonScheme_t comp;
  comp.title=type;
  if(type=="mgpy")
    {
      comp.nomFromSyst=false;
      comp.nom      = "t#bar{t}";
      comp.nomLabel = "Z2^{*}";
      comp.systList.push_back("t#bar{t}systtunep11"); comp.systLabels.push_back("P11");      comp.systColors.push_back(kAzure);
    }
  else if (type=="p11")
    {
      comp.nomFromSyst=true;
      comp.nom = "t#bar{t}systtunep11";
      comp.nomLabel="P11";
      comp.systList.clear(); comp.systLabels.clear(); comp.systColors.clear();
      comp.systList.push_back("t#bar{t}systtunep11mpihi"); comp.systLabels.push_back("MPI hi");   comp.systColors.push_back(kAzure);
      comp.systList.push_back("t#bar{t}systtunep11tev");   comp.systLabels.push_back("TEV");      comp.systColors.push_back(kRed-6);
      comp.systList.push_back("t#bar{t}systtunep11nocr");  comp.systLabels.push_back("No CR");    comp.systColors.push_back(kGreen+3);
    }
  else if (type=="mgpyq2scale")
    {
      comp.nomFromSyst=false;
      comp.nom      = "t#bar{t}";
      comp.nomLabel = "Z2^{*}";
      comp.systList.clear(); comp.systLabels.clear(); comp.systColors.clear();
      comp.systList.push_back("t#bar{t}systq2down"); comp.systLabels.push_back("(Q/2)^{2}");  comp.systColors.push_back(kAzure);
      comp.systList.push_back("t#bar{t}systq2up");   comp.systLabels.push_back("(2Q)^{2}");   comp.systColors.push_back(kRed-6);
    }
  else if (type=="mgpymepsscale")
    {
      comp.nomFromSyst=false;
      comp.nom      = "t#bar{t}";
      comp.nomLabel = "Z2^{*}";
      comp.systList.clear(); comp.systLabels.clear();  comp.systColors.clear();
      comp.systList.push_back("t#bar{t}systmepsdown"); comp.systLabels.push_back("-ME-PS"); comp.systColors.push_back(kAzure);
      comp.systList.push_back("t#bar{t}systmepsup");   comp.systLabels.push_back("+ME-PS"); comp.systColors.push_back(kRed-6);
    }
  /*
  else if (type=="py2hw")
    {
      comp.nomFromSyst=true;
      comp.nom      = "t#bar{t}systpowhegpy";
      comp.nomLabel = "Z2^{*}";
      comp.systList.clear(); comp.systLabels.clear(); comp.systColors.clear();
      comp.systList.push_back("t#bar{t}systpowheghw"); comp.systLabels.push_back("Herwig");    comp.systColors.push_back(kAzure);
    }
  */
  else if (type=="py2hw")
    {
      comp.nomFromSyst=false;
      comp.nom      = "t#bar{t}";
      comp.nomLabel = "Z2^{*}";
      comp.systList.clear(); comp.systLabels.clear(); comp.systColors.clear();
      comp.systList.push_back("t#bar{t}systmcatnlo"); comp.systLabels.push_back("Herwig");    comp.systColors.push_back(kAzure);
    }
  m_compSchemes.push_back( comp );
  
  //run the comparison for these distributions
  profile("nchprofpt");
  profile("nchprofavgptflux");
  profile("ptfluxprofpt");
  profile("avgptfluxprofpt");
  profile("nchprofphi");
  profile("ptfluxprofphi");
  profile("avgptfluxprofphi");
  profile("nchprofnvtx");
}



void showThrustProfile()
{
  //  TString proc("t#bar{t}systtunep11");
  TString proc("t#bar{t}systtunep11tev");

  //get histos from file
  TFile *inF=TFile::Open("plotter_syst_ue.root");
  TH2 *ptResponse=(TH2 *)inF->Get(proc+"/emu_ptresponse");
  ptResponse->Add((TH2 *)inF->Get(proc+"/ll_ptresponse"));
  ptResponse->SetDirectory(0);
  TH2 *phiResponse=(TH2 *)inF->Get(proc+"/emu_phiresponse");
  phiResponse->Add((TH2 *)inF->Get(proc+"/ll_phiresponse"));
  phiResponse->SetDirectory(0);
  TH2 *phivsptResponse=(TH2 *)inF->Get(proc+"/emu_phivspt");
  phivsptResponse->Add((TH2 *)inF->Get(proc+"/ll_phivspt"));
  phivsptResponse->SetDirectory(0);
  inF->Close();
  
  //show correlation plots
  gStyle->SetPalette(55);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);
  TCanvas *c=new TCanvas("c","c",1500,500);
  c->Divide(3,1);
  TPad *p=(TPad *)c->cd(1);
  ptResponse->Draw("colz");
  ptResponse->GetYaxis()->SetTitle("#Delta p_{T} [GeV]");
  ptResponse->GetXaxis()->SetMoreLogLabels();
  p->SetLogx();
  p->SetLogz();
  p->SetRightMargin(0.15);
  p->SetLeftMargin(0.15);
  p->SetTopMargin(0.05);
  p->SetBottomMargin(0.1);
  drawCMSHeader("",0.04,true);
  p=(TPad *)c->cd(2);
  phiResponse->Draw("colz");
  p->SetRightMargin(0.15);
  p->SetLeftMargin(0.15);
  p->SetTopMargin(0.05);
  p->SetBottomMargin(0.1);
  p->SetLogz();
  p=(TPad *)c->cd(3);
  phivsptResponse->Draw("colz");
  phivsptResponse->GetXaxis()->SetMoreLogLabels();
  p->SetRightMargin(0.15);
  p->SetLeftMargin(0.15);
  p->SetTopMargin(0.05);
  p->SetBottomMargin(0.1);
  p->SetLogx();
  p->SetLogz();

  TF1 *fitFunc = new TF1("ffunc","gaus");
  for(size_t i=0; i<3; i++){
    TString grName("gr"); grName+=i;
    TGraphErrors *gr=new TGraphErrors; gr->SetMarkerStyle(20); gr->SetName(grName);
    TH2 *h=ptResponse;
    if(i==1) h=phiResponse;
    if(i==2) h=phivsptResponse;

    for(int xbin=1; xbin<=h->GetXaxis()->GetNbins(); xbin++)
      {
	TH1 *py=h->ProjectionY(grName+"py",xbin,xbin);

	fitFunc->SetRange(h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax());
	fitFunc->SetParLimits(1,0,h->GetXaxis()->GetXmax());
	py->Fit(fitFunc,"Q0R");
	Int_t ip=gr->GetN();
	float x=h->GetXaxis()->GetBinCenter(xbin);
	float ex=h->GetXaxis()->GetBinWidth(xbin)/2;
	gr->SetPoint(ip,x,fitFunc->GetParameter(1));
	gr->SetPointError(ip,ex,fitFunc->GetParameter(2));
	py->Delete();
      }
    c->cd(i+1);
    gr->SetMarkerStyle(20);
    gr->SetLineWidth(2);
    gr->Draw("p");
  }
}


void drawGapFor(TString plot="_softleadpt")
{
  std::map<TString, TH1 *> histos;

  //get nominal prediction and subtract backgrounds from data
  TString dataProc="data";
  TString signalProc="t#bar{t}";
  TString bckgProcs[]={"VV","W,multijets","other t#bar{t}","Single top","Z#rightarrow ll"};
  TFile *inF=TFile::Open("plotter_ue.root");
  histos["data"]     = (TH1 *) inF->Get(dataProc+"/"+plot)->Clone("data");        histos["data"]->SetDirectory(0);
  for(size_t i=0; i<sizeof(bckgProcs)/sizeof(TString); i++)
    histos["data"]->Add( (TH1 *) inF->Get(bckgProcs[i]+"/"+plot), -1 );
  histos["signal"]   = (TH1 *) inF->Get(signalProc+"/"+plot)->Clone("signal");    histos["signal"]->SetDirectory(0);
  inF->Close();

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(1);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.1);
  TCanvas *c=new TCanvas("c","c",600,600);
 
  TPad *p=new TPad("p1","p1",0,0.53,1,1);
  p->Draw();
  p->cd();
  p->SetBottomMargin(0.01);
  TGraphErrors *sigcdf= getEvolutionWithErrors(histos["signal"]);
  sigcdf->SetMarkerStyle(24);
  sigcdf->Draw("ap");
  sigcdf->GetYaxis()->SetRangeUser(0.05,4);
  sigcdf->GetXaxis()->SetRangeUser(0,5);
  TString htype(histos["signal"]->ClassName());
  if(htype.Contains("TH1"))  sigcdf->GetYaxis()->SetTitle("Gap fraction");
  else sigcdf->GetYaxis()->SetTitle( "<" + TString( histos["signal"]->GetYaxis()->GetTitle() ) +">" );
  sigcdf->GetXaxis()->SetTitleSize(0);
  sigcdf->GetXaxis()->SetLabelSize(0);
  sigcdf->GetYaxis()->SetTitleSize(0.06);
  sigcdf->GetYaxis()->SetLabelSize(0.04);
  sigcdf->GetYaxis()->SetTitleOffset(0.75);
  
  TGraphErrors *datacdf= getEvolutionWithErrors(histos["data"]);
  datacdf->Draw("p");
  datacdf->SetMarkerStyle(20);
  drawCMSHeader("",0.04);

  TLegend *leg=new TLegend(0.6,0.75,0.9,0.85);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.09);
  leg->AddEntry(datacdf,"data","p");
  leg->AddEntry(sigcdf,"simulation","p");
  leg->SetNColumns(2);
  leg->Draw();


  //systematic variations
  inF=TFile::Open("plotter_syst_ue.root");
  for(int isyst=0; isyst<2; isyst++)
    {
      c->cd();
      if(isyst==0) { p=new TPad("p2","p2",0,0.53,1,0.28); p->SetBottomMargin(0.008); }
      else         { p=new TPad("p3","p3",0,0.28,1,0);   p->SetBottomMargin(0.2); }
      p->Draw();
      p->cd();
      p->SetTopMargin(0);
      ComparisonScheme_t comp;
      if(isyst==0){
	comp.systList.push_back("t#bar{t}systtunep11");      comp.systLabels.push_back("P11");      comp.systColors.push_back(1);
	comp.systList.push_back("t#bar{t}systtunep11mpihi"); comp.systLabels.push_back("MPI hi");   comp.systColors.push_back(kAzure);
	comp.systList.push_back("t#bar{t}systtunep11tev");     comp.systLabels.push_back("TEV");      comp.systColors.push_back(kRed-6);
	comp.systList.push_back("t#bar{t}systtunep11nocr");    comp.systLabels.push_back("No CR");    comp.systColors.push_back(kGreen+3);
      }
      else{
	 comp.systList.push_back("t#bar{t}");           comp.systLabels.push_back("Z2^{*}");     comp.systColors.push_back(1 );
	 comp.systList.push_back("t#bar{t}systq2down"); comp.systLabels.push_back("(Q/2)^{2}");  comp.systColors.push_back(kAzure);
	 comp.systList.push_back("t#bar{t}systq2up");   comp.systLabels.push_back("(2Q)^{2}");   comp.systColors.push_back(kRed-6);
      }
      
      TLegend *leg=new TLegend(0.6,0.84,0.9,0.9);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.09);
      leg->SetNColumns(comp.systList.size());
      for(size_t i=0; i<comp.systList.size(); i++)
	{
	  TString procName( comp.systList[i] );
	  
	  TGraphErrors *cdfgr( procName=="t#bar{t}" ? sigcdf : getEvolutionWithErrors( (TH1 *)inF->Get(procName+"/"+plot) ) );
	  TGraphErrors *gr( getRatio(cdfgr,datacdf) );
	  gr->SetFillColor(comp.systColors[i]);
	  gr->SetFillStyle(3002);
	  gr->SetLineColor(comp.systColors[i]);
	  gr->SetMarkerStyle(20);
	  gr->Draw(i==0? "ap" : "3");
	  leg->AddEntry(gr,comp.systLabels[i],i==0 ? "p" : "f");
	  if(i==0){
	    gr->GetYaxis()->SetTitle("(Data-#Sigma Bkg)/Signal");
	    gr->GetXaxis()->SetTitle( histos["signal"]->GetXaxis()->GetTitle() );
	    gr->GetXaxis()->SetTitleSize(isyst==0 ? 0 : 0.09);
	    gr->GetXaxis()->SetLabelSize(isyst==0 ? 0 : 0.07);
	    gr->GetXaxis()->SetTitleOffset(1);
	    gr->GetYaxis()->SetTitleOffset(0.5);
	    gr->GetYaxis()->SetTitleSize(0.09);
	    gr->GetYaxis()->SetLabelSize(0.07);
	    gr->GetYaxis()->SetNdivisions(5);
	    gr->GetYaxis()->SetRangeUser(0.76,1.24);
	    gr->GetXaxis()->SetRangeUser(0,5);
	  }
	}
      leg->Draw();
    }
  inF->Close();


}
