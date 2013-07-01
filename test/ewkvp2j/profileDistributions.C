#include "TFile.h"
#include "TString.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "THStack.h"
#include "TLine.h"

#include <iostream>

using namespace std;

//
void showCMSHeader(TObject *data=0, TObject *mc=0,bool isSim=false)
{
  TPaveText *pt=new TPaveText(0.1,0.91,0.8,0.99,"brNDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->SetTextAlign(12);
  //  pt->SetTextSize(0.08);
  if(!isSim) pt->AddText("CMS preliminary, #sqrt{s}=8 TeV, #scale[0.5]{#int}L=19.7 fb^{-1}");
  else       pt->AddText("CMS simulation, #sqrt{s}=8 TeV");
  pt->Draw();
  
  if(data==0 && mc==0) return;
  TLegend *leg=new TLegend(0.15,0.15,0.7,0.2);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetNColumns(2);
  //leg->SetTextSize(0.08);
  if(data) leg->AddEntry(data,"data","p");
  if(mc) leg->AddEntry(mc,"simulation","p");
  leg->Draw();
}

//
void drawDistributionWithVariations(TString dist="mumu_jetpt")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(1);
  TString ch[]={"#alpha^{4}_{EW}-ll+2j","VV","QCD-W+jets","Top","Z#rightarrow ll"};
  int colors[]={804,592,809,824,831};
  const size_t nch(sizeof(ch)/sizeof(TString));
  

  THStack *stack=new THStack("MC","MC");
  TH1 *totalMC=0,*totalMCup=0,*totalMCdown=0;
  TLegend *leg=new TLegend(0.7,0.7,0.9,0.9,"","brNDC"); leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42);
  TFile *llF=TFile::Open("~/work/ewkzp2j_539/plotter.root");
  for(size_t i=0; i<nch; i++)
    {
      TH1 *h=(TH1 *)llF->Get(ch[i]+"/"+dist);            
      h->SetTitle(ch[i]);
      h->SetFillColor(colors[i]);
      h->SetFillStyle(1001);
      //      cout << ch[i]+"/"+dist << " " << h << endl;
      h->SetDirectory(0);
      TH1 *hup=(TH1 *)llF->Get(ch[i]+"/"+dist+"_jesup");        hup->SetDirectory(0);
      TH1 *hdown=(TH1 *)llF->Get(ch[i]+"/"+dist+"_jesdown");    hdown->SetDirectory(0);
      if(totalMC==0){
	totalMC=(TH1 *)h->Clone("totalmc");         totalMC->Reset("ICE");     totalMC->SetDirectory(0);
	totalMCup=(TH1 *)h->Clone("totalmcup");     totalMCup->Reset("ICE");   totalMCup->SetDirectory(0);
	totalMCdown=(TH1 *)h->Clone("totalmcdown"); totalMCdown->Reset("ICE"); totalMCdown->SetDirectory(0);
      }
      stack->Add(h,"HIST");
      totalMC->Add(h);
      totalMCup->Add(hup);
      totalMCdown->Add(hdown);
      leg->AddEntry(h,h->GetTitle(),"F");
    }
  TH1 *data=(TH1*)llF->Get("data/"+dist); data->SetDirectory(0);
  leg->AddEntry(data,"data","p");
  llF->Close();

  TCanvas *c=new TCanvas("c","c",600,600);
  TPad* t1 = new TPad("t1","t1", 0.0, 0.20, 1.0, 1.0);
  t1->Draw();
  t1->SetLogy();
  t1->cd();
  stack->Draw("");
  stack->GetXaxis()->SetTitle(data->GetXaxis()->GetTitle());
  stack->GetYaxis()->SetTitle(data->GetYaxis()->GetTitle());
  stack->GetYaxis()->SetRangeUser(0.5,data->GetMaximum()*1.1);
  data->Draw("e1same");
  leg->Draw();
  showCMSHeader();

  c->cd();
  TPad* t2 = new TPad("t2","t2", 0.0, 0.0, 1.0, 0.2);
  t2->Draw();
  t2->cd();
  t2->SetGridy(true);
  t2->SetPad(0,0.0,1.0,0.2);
  t2->SetTopMargin(0);
  t2->SetBottomMargin(0.5);
  TH1 *nomRatio=(TH1 *)data->Clone("nomratio"); 
  nomRatio->Divide(totalMC); 
  nomRatio->SetDirectory(0); 
  nomRatio->SetMarkerStyle(20); 
  nomRatio->SetMarkerColor(1); 
  nomRatio->SetLineColor(1); 
  nomRatio->Draw("e1");
  nomRatio->GetYaxis()->SetTitle("Data/#Sigma MC");
  nomRatio->GetXaxis()->SetTitle("");
  nomRatio->GetYaxis()->SetRangeUser(0.46,1.54);
  nomRatio->GetYaxis()->SetTitleSize(0.12);
  nomRatio->GetYaxis()->SetTitleOffset(0.5);
  nomRatio->GetYaxis()->SetLabelSize(0.1);
  nomRatio->GetXaxis()->SetTitleSize(0.12);
  nomRatio->GetXaxis()->SetLabelSize(0.1);
  TH1 *upRatio=(TH1 *)data->Clone("upratio"); 
  upRatio->Divide(totalMCup); 
  upRatio->SetLineColor(kRed);
  upRatio->Draw("histsame");
  TH1 *downRatio=(TH1 *)data->Clone("downratio"); 
  downRatio->Divide(totalMCdown); 
  downRatio->SetLineColor(kBlue);
  downRatio->Draw("histsame");

  c->cd();
  c->Modified();
  c->Update();
}

//project pu jet kinematics
void projectPUjetKinematics()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(55);
  
  TFile *inF=TFile::Open("~/work/ewkzp2j_539/plotter.root");
  TH2 *incj=(TH2 *) inF->Get("Z#rightarrow ll/mumu_ptllvsdphi");        incj->SetDirectory(0);
  TH2 *pu=(TH2 *) inF->Get("Z#rightarrow ll/mumu_ptllvsdphipu");        pu->SetDirectory(0); 
  TH2 *truej=(TH2 *) inF->Get("Z#rightarrow ll/mumu_ptllvsdphitrue");   truej->SetDirectory(0); 
  inF->Close();

  TH2 *puPur=(TH2 *)pu->Clone("pupur");        puPur->Divide(incj);   puPur->GetZaxis()->SetTitle("Purity");   puPur->SetDirectory(0);
  TH2 *truePur=(TH2 *)truej->Clone("truepur"); truePur->Divide(incj); truePur->GetZaxis()->SetTitle("Purity"); truePur->SetDirectory(0);

  //draw
  TCanvas *c=new TCanvas("c","c",1000,1000);    
  c->Divide(2,2);

  TPad *p=(TPad *)c->cd(1);  p->SetLogz(); 
  pu->Scale(1./pu->Integral()); pu->Draw("colz");
  showCMSHeader(0,0,true); 
  pu->GetZaxis()->SetRangeUser(1e-3,1);
  TPaveText *pt=new TPaveText(0.6,0.3,0.9,0.4,"brNDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->SetTextFont(42);
  pt->AddText("Pileup jets");
  pt->Draw();

  p=(TPad *)c->cd(2);  p->SetLogz(); 
  truej->Scale(1./truej->Integral());  truej->Draw("colz");
  truej->GetZaxis()->SetRangeUser(1e-3,1);
  pt=new TPaveText(0.6,0.3,0.9,0.4,"brNDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->SetTextFont(42);
  pt->AddText("True jets");
  pt->Draw();


  p=(TPad *)c->cd(3); puPur->Draw("colz");
  TLine *l=new TLine(0,1,15,1); l->SetLineColor(kRed); l->SetLineWidth(2); l->Draw();
  l=new TLine(15,1,15,0);       l->SetLineColor(kRed); l->SetLineWidth(2); l->Draw();
  pt=new TPaveText(0.6,0.3,0.9,0.4,"brNDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->SetTextFont(42);
  pt->AddText("Pileup jets purity");
  pt->Draw();

  p=(TPad *)c->cd(4); truePur->Draw("colz");
  l=new TLine(50,3.2,50,2.7);   l->SetLineColor(kRed); l->SetLineWidth(2); l->Draw();
  l=new TLine(50,2.7,200,2.7);  l->SetLineColor(kRed); l->SetLineWidth(2); l->Draw();
   pt=new TPaveText(0.6,0.3,0.9,0.4,"brNDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->SetTextFont(42);
  pt->AddText("True jets purity");
  pt->Draw();

  
}


//
void measureJetIdEfficiency()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(1);
  
  
  //  TString profs[]= {"recoilbalanceidlt15collinear", "recoilbalanceidgt50back2back"};
  //TString cats[] = {"Pileup control",               "True jets control"};
  TString profs[]= {"recoilbalanceidgt50back2back"};
  TString cats[] = {"True jets control"};
  const size_t nprofs=sizeof(profs)/sizeof(TString);

  TCanvas *c=new TCanvas("c","c",600,600);       c->Divide(1,nprofs);   c->cd(1)->SetTopMargin(0.1);
  TCanvas *csf=new TCanvas("csf","csf",600,600); csf->Divide(1,nprofs); csf->cd(1)->SetTopMargin(0.1);

  TFile *llF=TFile::Open("~/work/ewkzp2j_539/plotter.root");
  TLegend *dataprojLeg = new TLegend(0.2,0.2,0.3,0.5,"Data","brNDC"); dataprojLeg->SetBorderSize(0); dataprojLeg->SetFillStyle(0); dataprojLeg->SetTextFont(42); dataprojLeg->SetTextSize(0.05);
  TLegend *mcprojLeg   = new TLegend(0.3,0.2,0.5,0.5,"MC","brNDC");   mcprojLeg->SetBorderSize(0);   mcprojLeg->SetFillStyle(0);   mcprojLeg->SetTextFont(42);   mcprojLeg->SetTextSize(0.05); 
  TLegend *data2mcLeg  = (TLegend *) dataprojLeg->Clone();
  for(size_t i=0; i<nprofs; i++)
    {
      TH2 *dyMC   = (TH2 *) llF->Get("Z#rightarrow ll/ee_"+profs[i]); 
      TH2 *dyData = (TH2 *) llF->Get("data/ee_"+profs[i]);      
  
      TH1 *allJetsMC=dyMC->ProjectionX("alljetsmc",1,1);    
      TH1 *allJetsData=dyData->ProjectionX("alljetsdata",1,1); 
      for(int ybin=2; ybin<=4+0*dyMC->GetYaxis()->GetNbins(); ybin+=2)
	{
	  TString label(dyMC->GetYaxis()->GetBinLabel(ybin));
	  TString pf(""); pf+=ybin; pf+=i;

	  TH1 *passJetsMC=dyMC->ProjectionX("passjetsmc",ybin,ybin); 
	  TH1 *passJetsData=dyData->ProjectionX("passjetsdata",ybin,ybin);
	  
	  TGraphAsymmErrors *mcGr=new TGraphAsymmErrors; mcGr->Divide(passJetsMC,allJetsMC,"n");   mcGr->SetMarkerStyle(24+ybin-2);      mcGr->SetFillStyle(0);      mcGr->SetName("mceff"+pf);        mcGr->SetTitle(label);
	  TGraphAsymmErrors *dataGr=new TGraphAsymmErrors; dataGr->Divide(passJetsData,allJetsData,"cp"); dataGr->SetMarkerStyle(20+ybin-2);    dataGr->SetFillStyle(0);    dataGr->SetName("dataeff"+pf);    dataGr->SetTitle("");
	  TGraphAsymmErrors *data2mcGr=(TGraphAsymmErrors *) dataGr->Clone("data2mc"+pf);
	  for(int ip=0; ip<mcGr->GetN(); ip++)
	    {
	      Double_t x,y1,ey1,y2,ey2;
	      mcGr->GetPoint(ip,x,y1);   ey1=mcGr->GetErrorY(ip);
	      dataGr->GetPoint(ip,x,y2); ey2=dataGr->GetErrorY(ip);
	      Double_t sf=y2/y1;  
	      Double_t sf_err=sqrt(pow(y2*ey1,2)+pow(y1*ey2,2))/pow(y1,2);
	      data2mcGr->SetPoint(ip,x,sf);
	      data2mcGr->SetPointError(ip,0,0,sf_err,sf_err);
	    }

	  if(i==0){
	    dataprojLeg->AddEntry(dataGr,"","pf");
	    mcprojLeg->AddEntry(mcGr,label,"pf");
	    data2mcLeg->AddEntry(data2mcGr,label,"pf");
	  }

	  c->cd(i+1);
	  mcGr->Draw(ybin==2 ? "ap" : "p");
	  if(ybin==2) {
	    mcGr->GetYaxis()->SetRangeUser(0,1.2); mcGr->GetYaxis()->SetTitle("Efficiency"); mcGr->GetXaxis()->SetTitle("Pseudo-rapidity"); 
	    mcGr->GetYaxis()->SetTitleSize(0.07);
	    mcGr->GetYaxis()->SetTitleOffset(0.9);
	    mcGr->GetYaxis()->SetLabelSize(0.06);
	    mcGr->GetXaxis()->SetTitleSize(0.07);
	    mcGr->GetXaxis()->SetLabelSize(0.06);
	    mcGr->GetXaxis()->SetTitleOffset(0.9);

	    TPaveText *pt=new TPaveText(0.15,0.5,0.8,0.85,"brNDC");
	    pt->SetBorderSize(0);
	    pt->SetFillStyle(0);
	    pt->SetTextAlign(13);
	    pt->SetTextFont(42);
	    pt->SetTextColor(kBlue);
	    pt->SetTextSize(0.08);
	    pt->AddText("[ "+cats[i]+" ]");
	    pt->Draw();
	  }
	  dataGr->Draw("p");
	  if(i==0)  showCMSHeader();
  
	  csf->cd(i+1);
	  data2mcGr->Draw(ybin==2 ? "ap" : "p");
	  if(ybin==2) { 
	    data2mcGr->GetYaxis()->SetRangeUser(0.74,1.26); data2mcGr->GetYaxis()->SetTitle("Data/MC scale factor"); data2mcGr->GetXaxis()->SetTitle("Pseudo-rapidity"); 
	    data2mcGr->GetYaxis()->SetTitleSize(0.07);
	    data2mcGr->GetYaxis()->SetTitleOffset(0.9);
	    data2mcGr->GetYaxis()->SetLabelSize(0.06);
	    data2mcGr->GetXaxis()->SetTitleSize(0.07);
	    data2mcGr->GetXaxis()->SetLabelSize(0.06);
	    data2mcGr->GetXaxis()->SetTitleOffset(0.9);

	    TPaveText *pt=new TPaveText(0.15,0.5,0.8,0.85,"brNDC");
	    pt->SetBorderSize(0);
	    pt->SetFillStyle(0);
	    pt->SetTextAlign(13);
	    pt->SetTextFont(42);
	    pt->SetTextColor(kBlue);
	    pt->SetTextSize(0.08);
	    pt->AddText("[ "+cats[i]+" ]");
	    pt->Draw();
	  }
	  if(i==0) showCMSHeader();
	  
	  //passJetsMC->Delete();
	  //passJetsData->Delete();
	  //data2MC->Delete();
	}
    }      
  llF->Close();

  //update the canvases
  c->cd(1); dataprojLeg->Draw();  mcprojLeg->Draw(); c->cd();   c->Modified();   c->Update();
  csf->cd(1); data2mcLeg->Draw();                    csf->cd(); csf->Modified(); csf->Update();
}


//
void profileResidualJEC()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(1);

  TFile *llF=TFile::Open("~/work/ewkzp2j_539/plotter.root");
  //  TFile *gF=TFile::Open("~/work/ewkzp2j_539/plotter_g_raw.root");
  
  TString profs[]={//"recoilbalancevseta",
    //"recoilbalancevseta30to50",
		   "recoilbalancevseta50toInf"};
  TString cats[]={//"inclusive",
    // "30<p_{T}/GeV<50",
		  "p_{T}>50 GeV"};
  const size_t nprofs=sizeof(profs)/sizeof(TString);

  TString ch[]={"ee","mumu"};
  const size_t nch=sizeof(ch)/sizeof(TString);
  
  TGraph *frame=new TGraph; 
  frame->SetName("frame");
  frame->SetPoint(0,0,0.4);
  frame->SetPoint(1,5,1.8);
  frame->SetMarkerStyle(1);
  TGraph *sfFrame=0;
  TObjArray toSave;
  for(size_t ich=0; ich<nch; ich++)
    {
      TCanvas *c=new TCanvas("c"+ch[ich],"c"+ch[ich],600,600); c->Divide(1,nprofs);
      TCanvas *cg=new TCanvas("cg","cg",600,600);              cg->Divide(1,nprofs);

      TCanvas *cr=new TCanvas("cr"+ch[ich],"cr"+ch[ich],600,600); cr->Divide(1,nprofs);
      TCanvas *cgr=new TCanvas("cgr","cgr",600,600);              cgr->Divide(1,nprofs);

      for(size_t i=0; i<nprofs; i++)
	{
	  TH2 *dyMC   = (TH2 *) llF->Get("Z#rightarrow ll/"+ch[ich]+"_"+profs[i]); 
	  TH2 *dyData = (TH2 *) llF->Get("data/"+ch[ich]+"_"+profs[i]);            
	  dyMC->SetDirectory(0);
	  dyData->SetDirectory(0);
	  
	  TH1D *dyMCProfH   = dyMC->ProfileX("dymcprof");
	  TGraphErrors *dyMCProf   = new TGraphErrors(dyMCProfH);
	  dyMCProf->SetMarkerStyle(24);  
	  dyMCProf->SetName(ch[ich]+profs[i]+"dymc");
	  
	  TH1D *dyDataProfH = dyData->ProfileX("dydataprof");	  
	  TGraphErrors *dyDataProf = new TGraphErrors(dyDataProfH);  
	  dyDataProf->SetMarkerStyle(20); 
	  dyDataProf->SetName(ch[ich]+profs[i]+"dydata");
	  
	  dyDataProfH->Divide(dyMCProfH);
	  TGraph *dyRatio = new TGraph(dyDataProfH); 
	  dyRatio->SetName(ch[ich]+profs[i]+"dydata2mc");

	  /*
	  TH2 *gMC   = (TH2 *) gF->Get("#gamma+jets/"+ch[ich]+"_"+profs[i]); 
	  gMC->Add(    (TH2 *) gF->Get("Multijets/"+ch[ich]+"_"+profs[i]) ); 
	  TH2 *gData = (TH2 *) gF->Get("data (#gamma)/"+ch[ich]+"_"+profs[i]);            
	  gMC->SetDirectory(0);
	  gData->SetDirectory(0);

	  TH1D *gMCProfH = gMC->ProfileX("gmcprof");
	  TGraphErrors *gMCProf    = new TGraphErrors(gMCProfH);    
	  gMCProf->SetMarkerStyle(24);  
	  gMCProf->SetName(ch[ich]+profs[i]+"gmc");
	 
	  TH1D *gDataProfH = gData->ProfileX("gdataprof");
	  TGraphErrors *gDataProf  = new TGraphErrors(gDataProfH);   
	  gDataProf->SetMarkerStyle(20); 
	  gDataProf->SetName(ch[ich]+profs[i]+"gdata");

	  gDataProfH->Divide(gMCProfH);
	  TGraph *gRatio = new TGraph(gDataProfH);
	  gRatio->SetName(ch[ich]+profs[i]+"gdata2mc");
	  */
	  toSave.Add(dyRatio);
	  
	  /*
	  toSave.Add(gRatio);
	  */
	  //
	  TPad *p=(TPad *)c->cd(i+1);
	  if(i<nprofs-1) p->SetBottomMargin(0.01);
	  if(i>0) p->SetTopMargin(0);
	  frame->Draw("ap");
	  frame->GetYaxis()->SetRangeUser(0.34,1.66);
	  frame->GetYaxis()->SetTitle("Jet response");
	  frame->GetXaxis()->SetRangeUser(0,4.75);
	  frame->GetXaxis()->SetTitle("Pseudo-rapidity");
	  frame->GetYaxis()->SetLabelSize(0.06);
	  frame->GetYaxis()->SetTitleSize(0.08);
	  frame->GetYaxis()->SetTitleOffset(0.5);
	  frame->GetXaxis()->SetLabelSize(0.06);
	  frame->GetXaxis()->SetTitleSize(0.08);
	  frame->GetXaxis()->SetTitleOffset(0.6);
	  dyDataProf->Draw("p");
	  dyMCProf->Draw("p");
	  if(i==0) showCMSHeader(dyDataProf,dyMCProf);

	  TPaveText *pt=new TPaveText(0.1,0.8,0.8,0.85,"brNDC");
	  pt->SetBorderSize(0);
	  pt->SetFillStyle(0);
	  pt->SetTextAlign(12);
	  pt->SetTextFont(42);
	  pt->SetTextSize(0.08);
	  TString label("Z#rightarrow "+ch[ich]+ ", " + cats[i]);
	  label.ReplaceAll("mumu","#mu#mu");
	  pt->AddText(label);
	  pt->Draw();


	  //
	  p=(TPad *)cr->cd(i+1);
	  if(i<nprofs-1) p->SetBottomMargin(0.01);
	  if(i>0) p->SetTopMargin(0);
	  if(sfFrame==0) {
	    sfFrame =(TGraph *) frame->Clone("sfframe");
	    sfFrame->GetYaxis()->SetTitle("Data/MC jet response");
	    sfFrame->GetYaxis()->SetRangeUser(0.86,1.14);
	  }
	  sfFrame->Draw("ap");
	  dyRatio->Draw("l");
	  if(i==0) showCMSHeader();
	  pt->Clone()->Draw();

	  //
	  /*
	  p=(TPad *)cg->cd(i+1);
	  if(i<nprofs-1) p->SetBottomMargin(0.01);
	  if(i>0) p->SetTopMargin(0);
	  frame->Draw("ap");
	  gDataProf->Draw("p");
	  gMCProf->Draw("p");

	  if(i==0) showCMSHeader(gDataProf,gMCProf);

	  pt=new TPaveText(0.1,0.8,0.8,0.85,"brNDC");
	  pt->SetBorderSize(0);
	  pt->SetFillStyle(0);
	  pt->SetTextAlign(12);
	  pt->SetTextFont(42);
	  pt->SetTextSize(0.08);
	  pt->AddText("#gamma, " + cats[i]);
	  pt->Draw();


	  //
	  p=(TPad *)cgr->cd(i+1);
	  if(i<nprofs-1) p->SetBottomMargin(0.01);
	  if(i>0) p->SetTopMargin(0);
	  sfFrame->Draw("ap");
	  gRatio->Draw("l");
	  if(i==0) showCMSHeader();
	  pt->Clone()->Draw();
	  */
	}

      c->SaveAs(ch[ich]+"jetresponse.pdf");
      cg->SaveAs("gjetresponse.pdf");
      cr->SaveAs(ch[ich]+"jetresponseSF.pdf");
      cgr->SaveAs("gjetresponseSF.pdf");
    }
  llF->Close();
  // gF->Close();
  

  //save the results
  TFile *fOut = TFile::Open("jecResponse.root","RECREATE");
  for(int i=0; i<toSave.GetEntriesFast(); i++)
      toSave.At(i)->Write();
  fOut->Close();
}
