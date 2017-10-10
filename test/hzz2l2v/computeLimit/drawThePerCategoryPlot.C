#include "UserCode/llvv_fwk/test/hzz2l2v/computeLimit/tdrstyle.C"
#include "UserCode/llvv_fwk/test/hzz2l2v/computeLimit/CMS_lumi.C"
#include <map>
#include "TH1D.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "TPad.h"
#include "THStack.h"
#include "TFile.h"
#include "TMath.h"
#include "TGaxis.h"


TGraphAsymmErrors *convertTheHistoInTGraph(TH1D *theInputHisto){
  TGraphAsymmErrors *theGraph = new TGraphAsymmErrors(theInputHisto);
  const double alpha = 1 - 0.6827;
  for (int i = 0; i < theGraph->GetN(); ++i) {
    int N = theGraph->GetY()[i]; //Nominal
    double L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.)); //Low
    double U =  ROOT::Math::gamma_quantile_c(alpha/2,N+1,1); //Up
    theGraph->SetPointEYlow(i, N-L);
    theGraph->SetPointEYhigh(i, U-N);
  }
  return theGraph;
}

void drawPlots(TGraphAsymmErrors* theGraph, THStack* hs, TH1D* ggH_sigInterf, TH1D* qqH_sigInterf, TGraphAsymmErrors* gr_statUnc, TGraphAsymmErrors* gr_systUnc, std::vector<TString> v_listMCProcesses, TString theLeptonCategoryText, TString theJetCategoryText, TString categoryName, std::map<TString, TH1D*> map_MCProcesses, std::vector<TString> v_legendNames, bool StatThenSyst, bool SystNormThenSystShapeThenStat, TH1D *h_MC_shape_up, TH1D * h_MC_shape_down, bool mergeSyst, TFile* systFile, bool doErrorPreFit, double sigma_signal){

  int W = 800;
  int H = 600;
  int H_ref = 600;
  int W_ref = 800;

  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.12*H_ref;
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;

  //Compute TGraph in the case of SystNormThenSystShapeThenStat = true 
  TGraphAsymmErrors* gr_systNandSUnc = (TGraphAsymmErrors*) gr_statUnc->Clone();

  for (int i = 0; i < gr_systNandSUnc->GetN(); ++i) {
    gr_systNandSUnc->SetPointEYlow(i, h_MC_shape_up->GetBinContent(i+1));
    gr_systNandSUnc->SetPointEYhigh(i, h_MC_shape_down->GetBinContent(i+1));
  }


  TCanvas *c0 = new TCanvas("c0","coucou",50,50,W,H);
  c0->SetFillColor(0);
  c0->SetBorderMode(0);
  c0->SetFrameFillStyle(0);
  c0->SetFrameBorderMode(0);
  c0->SetLeftMargin( L/W );
  c0->SetRightMargin( R/W );
  c0->SetTopMargin( T/H );
  c0->SetBottomMargin( B/H );
  c0->SetTickx(1);
  c0->SetTicky(1);
  c0->SetLogy();

	//Fix the size of the canvas: xmin,ymin,xmax,ymax
	gPad->DrawFrame(0.0, 0.001, 3000.0, 1000);
	gPad->Update();

  hs->GetXaxis()->SetTitle("M_{T} [GeV]");
  hs->GetYaxis()->SetTitleOffset(1);
  hs->GetYaxis()->SetTitle("Events");
  hs->Draw("same:hist");
  ggH_sigInterf->SetLineWidth(2);
  ggH_sigInterf->SetLineColor(619);
  ggH_sigInterf->Draw("same:hist");
  qqH_sigInterf->SetLineWidth(2);
  qqH_sigInterf->SetLineColor(619);
  qqH_sigInterf->SetLineStyle(2);
  qqH_sigInterf->Draw("same:hist");
  theGraph->Draw("same:P0");
  gr_statUnc->SetFillStyle(3005);
  gr_statUnc->SetFillColor(kGray+3);
  gr_statUnc->SetLineStyle(1);
  gr_statUnc->SetLineColor(1);
  gr_statUnc->Draw("2 same");
  if(!mergeSyst){
  	gr_systUnc->SetFillStyle(3004);
  	if(SystNormThenSystShapeThenStat) gr_systUnc->SetFillStyle(3016);
  	//gr_systUnc->SetFillColor(kGray+1);
  	gr_systUnc->SetFillColor(kBlue+3);
  	gr_systUnc->SetLineStyle(1);
  	gr_systUnc->SetLineColor(1);
		gr_systUnc->Draw("2 same");
	}
	if(SystNormThenSystShapeThenStat){
  	gr_systNandSUnc->SetFillStyle(3004);
  	gr_systNandSUnc->SetFillColor(kMagenta);
  	gr_systNandSUnc->SetLineStyle(1);
  	gr_systNandSUnc->SetLineColor(1);
  	gr_systNandSUnc->Draw("2 same");
  }
  gPad->RedrawAxis();

  CMS_lumi( c0, 16, 0, false);

  TLegend *t = new TLegend(0.38,0.62,0.95,0.85);
  t->SetLineColor(0);
  t->SetBorderSize(1);
  t->SetNColumns(3);
  t->AddEntry(theGraph, "data","ELP");
  for (unsigned int i= 0; i < v_listMCProcesses.size(); i++) t->AddEntry(map_MCProcesses[v_listMCProcesses[i]],v_legendNames[i],"f");
  t->AddEntry(ggH_sigInterf,"ggF X(M=800;W=100;#sigma=50fb)","l");
  t->AddEntry(qqH_sigInterf,"VBF X(M=800;W=100;#sigma=30fb)","l");
  if(StatThenSyst){
    t->AddEntry(gr_statUnc,"Stat.","F");
    t->AddEntry(gr_systUnc,"Stat. + Syst. (Norm)","F");
  }
  else if(SystNormThenSystShapeThenStat){
    if(!mergeSyst) t->AddEntry(gr_systUnc,"Syst. (Norm)","F");
    if(mergeSyst) t->AddEntry(gr_systNandSUnc,"Syst.","F");
    else t->AddEntry(gr_systNandSUnc,"Syst. (Norm + Shape)","F");
    t->AddEntry(gr_statUnc,"Syst. + Stat.","F");
  }
  else{
    t->AddEntry(gr_systUnc,"Syst. (Norm)","F");
    t->AddEntry(gr_statUnc,"Syst. (Norm) + Stat.","F");
  }
  t->Draw();


  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);

  latex.SetTextFont(42);
  latex.SetTextAlign(31);
  latex.SetTextSize(0.05);
  latex.DrawLatex(0.94,0.86,theLeptonCategoryText+ "   " + theJetCategoryText);

  c0->Print(categoryName+".png");
  c0->Print(categoryName+".pdf");

}

void drawRatioPlots(TGraphAsymmErrors* theGraph, THStack* hs, TH1D* ggH_sigInterf, TH1D* qqH_sigInterf, TGraphAsymmErrors* gr_statUnc, TGraphAsymmErrors* gr_systUnc, std::vector<TString> v_listMCProcesses, TString theLeptonCategoryText, TString theJetCategoryText, TString categoryName, std::map<TString, TH1D*> map_MCProcesses, std::vector<TString> v_legendNames, TH1D* h_MC_shape_up, TH1D* h_MC_shape_down, bool showSystShape, bool StatThenSyst, bool SystNormThenSystShapeThenStat, bool mergeSyst, TFile* systFile, bool doErrorPreFit, TString errorFitType, double sigma_signal){

  //////////////////////////////////////////
  /////////// Configuration part ///////////
  //////////////////////////////////////////

 	int fontType = 63; //precision 3 font, so label size is expressed in pixel now
  int pixelFontSize = 20;

  int W = 800;
  int H = 600;
  int H_ref = 600;
  int W_ref = 800;



  //references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.30*H_ref;
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;

	//////////////////////////////////////////
  //////// End of Configuration part ///////
  //////////////////////////////////////////

	//Compute TGraph in the case of SystNormThenSystShapeThenStat = true 
  TGraphAsymmErrors* gr_systNandSUnc = (TGraphAsymmErrors*) gr_statUnc->Clone();

  for (int i = 0; i < gr_systNandSUnc->GetN(); ++i) {
    gr_systNandSUnc->SetPointEYlow(i, h_MC_shape_up->GetBinContent(i+1));
    gr_systNandSUnc->SetPointEYhigh(i, h_MC_shape_down->GetBinContent(i+1));
  }

  //Compute pre-fit syst
  TH1D *totalBackground_fromFit = (TH1D*) ((TH1D*) systFile->Get(errorFitType+"/"+categoryName+"/total_background"))->Clone();
  TGraphAsymmErrors* gr_preFitSyst = (TGraphAsymmErrors*) gr_statUnc->Clone();
  TH1D* nominalHisto = (TH1D*) hs->GetStack()->Last()->Clone();
  for (int i = 0; i < gr_preFitSyst->GetN(); ++i) {
    gr_preFitSyst->SetPointEYlow(i, totalBackground_fromFit->GetBinError(i+1));
    gr_preFitSyst->SetPointEYhigh(i, totalBackground_fromFit->GetBinError(i+1));
	}



  TCanvas *c0 = new TCanvas("c0",categoryName,50,50,W,H);

	//Fix the size of the canvas: xmin,ymin,xmax,ymax
  TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
  pad1->SetLeftMargin( L/W );
  pad1->SetRightMargin( R/W );
  pad1->SetTopMargin( T/H );
  pad1->SetBottomMargin(0);
	pad1->Draw();
  pad1->cd();

  TH1F *hFrame1 = pad1->DrawFrame(0.0, 0.00101, 3000.0, 1300, ";M_{T} [GeV];Events");
  pad1->SetLogy();
  pad1->Modified();
  pad1->Update();

  hFrame1 = ((TH1F *)(gPad->FindObject("hframe")));
  hFrame1->GetYaxis()->SetLabelFont(fontType);
  hFrame1->GetYaxis()->SetLabelSize(pixelFontSize);
  hFrame1->GetYaxis()->SetNdivisions(505);
  hFrame1->GetYaxis()->SetTitleOffset(1.4);
  hFrame1->GetYaxis()->SetTitleFont(fontType);
  hFrame1->GetYaxis()->SetTitleSize(pixelFontSize);

  hs->Draw("same:hist");
  ggH_sigInterf->SetLineWidth(2);
  ggH_sigInterf->SetLineColor(619);
  ggH_sigInterf->Draw("same:hist");
  qqH_sigInterf->SetLineWidth(2);
  qqH_sigInterf->SetLineColor(619);
  qqH_sigInterf->SetLineStyle(2);
  qqH_sigInterf->Draw("same:hist");
  theGraph->Draw("same:P0");
  if(doErrorPreFit){
  	gr_preFitSyst->SetFillStyle(3005);
  	gr_preFitSyst->SetFillColor(kGray+3);
  	gr_preFitSyst->SetLineStyle(1);
  	gr_preFitSyst->SetLineColor(1);
  	gr_preFitSyst->Draw("2 same");
  }
  else{
 		gr_statUnc->SetFillStyle(3005);
  	gr_statUnc->SetFillColor(kGray+3);
  	gr_statUnc->SetLineStyle(1);
  	gr_statUnc->SetLineColor(1);
  	gr_statUnc->Draw("2 same");

		if(!mergeSyst){
  		gr_systUnc->SetFillStyle(3004);
  		if(SystNormThenSystShapeThenStat) gr_systUnc->SetFillStyle(3016);
  		//gr_systUnc->SetFillColor(kGray+1);
  		gr_systUnc->SetFillColor(kBlue+3);
  		gr_systUnc->SetLineStyle(1);
  		gr_systUnc->SetLineColor(1);
  		gr_systUnc->Draw("2 same");
  	}
  	if(SystNormThenSystShapeThenStat){
  		gr_systNandSUnc->SetFillStyle(3004);
  		gr_systNandSUnc->SetFillColor(kMagenta);
  		gr_systNandSUnc->SetLineStyle(1);
  		gr_systNandSUnc->SetLineColor(1);
  		gr_systNandSUnc->Draw("2 same");
  	}
  }

  gPad->RedrawAxis();

  CMS_lumi( pad1, 16, 0, false);

	double legend_posXlow = 0.38;
	double legend_posXhigh = 0.95;
	double legend_posYlow = 0.68;
	double legend_posYhigh = 0.90;

 	TLegend *t = new TLegend(legend_posXlow, legend_posYlow, legend_posXhigh, legend_posYhigh);
  t->SetLineColor(0);
  t->SetBorderSize(1);
  t->SetNColumns(3);
  t->AddEntry(theGraph, "data","ELP");
  for (unsigned int i= 0; i < v_listMCProcesses.size(); i++) t->AddEntry(map_MCProcesses[v_listMCProcesses[i]],v_legendNames[i],"f");
  if(doErrorPreFit){
    t->AddEntry(gr_preFitSyst,"Bkg. Unc.","F");
  }
  else{
  	if(StatThenSyst){
    	t->AddEntry(gr_statUnc,"Stat.","F");
    	t->AddEntry(gr_systUnc,"Stat. + Syst. (Norm)","F");
  	}
  	else if(SystNormThenSystShapeThenStat){
    	if(!mergeSyst) t->AddEntry(gr_systUnc,"Syst. (Norm)","F");
    	if(mergeSyst) t->AddEntry(gr_systNandSUnc,"Syst.","F");
    	else t->AddEntry(gr_systNandSUnc,"Syst. (Norm + Shape)","F");
    	t->AddEntry(gr_statUnc,"Syst. + Stat.","F");
  	}  
  	else{
    	t->AddEntry(gr_systUnc,"Syst. (Norm)","F");
    	t->AddEntry(gr_statUnc,"Syst. (Norm) + Stat.","F");
  	}
  }
  t->Draw();

  //TLegend *signal_text = new TLegend(legend_posXlow, legend_posYlow-0.1, legend_posXlow+(legend_posXhigh-legend_posXlow)*4./5., legend_posYlow);
  //TLegend *signal_text = new TLegend(legend_posXlow-0.155, legend_posYlow-0.1, legend_posXhigh-0.07, legend_posYlow);
  //TLegend *signal_text = new TLegend(legend_posXlow-0.155, legend_posYlow-0.065, legend_posXhigh-0.07, legend_posYlow-0.01);
  TLegend *signal_text = new TLegend(legend_posXlow-0.155, legend_posYlow-0.12, legend_posXhigh-0.07, legend_posYlow-0.01);
  signal_text->SetFillStyle(0);
  signal_text->SetLineColor(0);
  signal_text->SetBorderSize(0);
  //signal_text->AddEntry((TObject*)0, "Signal + Interference: (M_{X}, #Gamma_{X}) = (800,100) GeV", "");
  signal_text->AddEntry((TObject*)0, "signal + interference ", "");
  signal_text->AddEntry((TObject*)0, "(M_{X}, #Gamma_{X})=(800,100) GeV", "");

	signal_text->Draw();

  /*TLegend *signal_legend = new TLegend(legend_posXlow+0.3, legend_posYlow-0.15, legend_posXhigh, legend_posYlow-0.1);
  signal_legend->SetLineColor(0);
  signal_legend->SetBorderSize(0);
  signal_legend->SetNColumns(3);
  signal_legend->AddEntry(ggH_sigInterf,"ggF","l");
  signal_legend->AddEntry((TObject*)0, "", "");
  signal_legend->AddEntry(qqH_sigInterf,"VBF","l");

	signal_legend->Draw();
*/
  TLegend *signal_legend = new TLegend(legend_posXlow+0.32, legend_posYlow-0.09, legend_posXhigh-0.05, legend_posYlow-0.03);
  signal_legend->SetLineColor(0);
  signal_legend->SetBorderSize(0);
  signal_legend->SetNColumns(2);
  signal_legend->SetFillStyle(0);
  //signal_legend->AddEntry((TObject*)0, "#bf{(M_{X}, #Gamma_{X})=(800,100) GeV}", "");
  //signal_legend->AddEntry((TObject*)0, "(M_{X}, #Gamma_{X})=(800,100) GeV", "");
  signal_legend->AddEntry(ggH_sigInterf,"ggF ","l");
  signal_legend->AddEntry(qqH_sigInterf,"VBF","l");

	signal_legend->Draw();

  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);

  latex.SetTextFont(42);
  //latex.SetTextAlign(31);
  latex.SetTextSize(0.07*6/5.);
  //latex.DrawLatex(0.15,0.85,"#bf{"+theCategoryText+"}");
  latex.DrawLatex(0.155,0.83,theLeptonCategoryText);

  TLatex latex2;
  latex2.SetNDC();
  latex2.SetTextAngle(0);
  latex2.SetTextColor(kBlack);

  //latex.SetTextFont(42);
  //latex.SetTextAlign(31);
  latex2.SetTextSize(0.07);
  //latex.DrawLatex(0.15,0.85,"#bf{"+theCategoryText+"}");
  if(theJetCategoryText.Contains("VBF")) latex2.DrawLatex(0.155,0.75,theJetCategoryText);
  else latex2.DrawLatex(0.215,0.829,theJetCategoryText);

  pad1->Modified();
  pad1->Update();

	//Pad 2
	c0->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.25);
  pad2->SetLeftMargin( L/W );
  pad2->SetRightMargin( R/W );
	pad2->SetTopMargin(0);
  pad2->SetBottomMargin( B/H );
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();

	TH1F *hFrame2 = pad2->DrawFrame(0.0, 0, 3000.0, 2.15, ";M_{T} [GeV];#frac{Data}{#Sigma Bkg.}");
	pad2->Modified();
	pad2->Update();

	hFrame2 = ((TH1F *)(gPad->FindObject("hframe")));
  hFrame2->GetYaxis()->SetLabelFont(fontType);
  hFrame2->GetYaxis()->SetLabelSize(pixelFontSize);
  hFrame2->GetYaxis()->SetNdivisions(505);
  hFrame2->GetYaxis()->CenterTitle();
  hFrame2->GetYaxis()->SetTitleFont(fontType);
  hFrame2->GetYaxis()->SetTitleSize(pixelFontSize);
  hFrame2->GetYaxis()->SetTitleOffset(1.4);
  hFrame2->GetXaxis()->SetLabelFont(fontType);
  hFrame2->GetXaxis()->SetLabelSize(pixelFontSize);
  hFrame2->GetXaxis()->SetTitleFont(fontType);
  hFrame2->GetXaxis()->SetTitleSize(pixelFontSize);
  hFrame2->GetXaxis()->SetTitleOffset(4);

	//data over MC points
  TGraphAsymmErrors* numerator = (TGraphAsymmErrors*) theGraph->Clone();
  TH1D* denominator = (TH1D*) hs->GetStack()->Last()->Clone();

	double num =0., num_up=0., num_down=0., den=0.;
	for (int i = 0; i < numerator->GetN(); ++i) {
    num = numerator->GetY()[i];
    num_up = numerator->GetEYhigh()[i];
    num_down = numerator->GetEYlow()[i];
    den = denominator->GetBinContent(i+1);

    numerator->SetPoint(i, numerator->GetX()[i], num/den);
    numerator->SetPointEYlow(i, num_down/den);
    numerator->SetPointEYhigh(i, num_up/den);
	}

  numerator->Draw("same:P0");

  //If errors are coming from the prefit
	if(doErrorPreFit){
    TGraphAsymmErrors* numerator_preFit = (TGraphAsymmErrors*) gr_preFitSyst->Clone();

  	double num_preFit =0., num_preFit_up=0., num_preFit_down=0., den_preFit=0.;
  	for (int i = 0; i < numerator_preFit->GetN(); ++i) {
    	num_preFit = numerator_preFit->GetY()[i];
    	num_preFit_up = numerator_preFit->GetEYhigh()[i];
    	num_preFit_down = numerator_preFit->GetEYlow()[i];
    	den_preFit = numerator_preFit->GetY()[i];;

    	numerator_preFit->SetPoint(i, numerator_preFit->GetX()[i], num_preFit/den_preFit);
    	numerator_preFit->SetPointEYlow(i, num_preFit_down/den_preFit);
    	numerator_preFit->SetPointEYhigh(i, num_preFit_up/den_preFit);
  	}

 	  numerator_preFit->Draw("2 same");


	}
	else{
		//Uncertainty band: stat
		TGraphAsymmErrors* numerator_stat = (TGraphAsymmErrors*) gr_statUnc->Clone();

  	double num_stat =0., num_stat_up=0., num_stat_down=0., den_stat=0.;
  	for (int i = 0; i < numerator_stat->GetN(); ++i) {
    	num_stat = numerator_stat->GetY()[i];
    	num_stat_up = numerator_stat->GetEYhigh()[i];
    	num_stat_down = numerator_stat->GetEYlow()[i];
    	den_stat = numerator_stat->GetY()[i];;

    	numerator_stat->SetPoint(i, numerator_stat->GetX()[i], num_stat/den_stat);
    	numerator_stat->SetPointEYlow(i, num_stat_down/den_stat);
    	numerator_stat->SetPointEYhigh(i, num_stat_up/den_stat);
  	}

		numerator_stat->Draw("2 same");

		//Uncertainty band: syst
		TGraphAsymmErrors* numerator_syst = (TGraphAsymmErrors*) gr_systUnc->Clone();

  	double num_syst =0., num_syst_up=0., num_syst_down=0., den_syst=0.;
  	for (int i = 0; i < numerator_syst->GetN(); ++i) {
    	num_syst = numerator_syst->GetY()[i];
    	num_syst_up = numerator_syst->GetEYhigh()[i];
    	num_syst_down = numerator_syst->GetEYlow()[i];
    	den_syst = numerator_syst->GetY()[i];;

    	numerator_syst->SetPoint(i, numerator_syst->GetX()[i], num_syst/den_syst);
    	numerator_syst->SetPointEYlow(i, num_syst_down/den_syst);
    	numerator_syst->SetPointEYhigh(i, num_syst_up/den_syst);
  	}

		if(!mergeSyst) numerator_syst->Draw("2 same");

		//Uncertainty shape: syst
		if(showSystShape && !SystNormThenSystShapeThenStat){
			TGraphAsymmErrors* shape_syst = (TGraphAsymmErrors*) gr_systUnc->Clone();
    	TH1D* denominator_shape = (TH1D*) hs->GetStack()->Last()->Clone();

  		double num_shape_up=0., num_shape_down=0., den_shape=0.;
  		for (int i = 0; i < shape_syst->GetN(); ++i) {
    		num_shape_up = h_MC_shape_up->GetBinContent(i+1);
    		num_shape_down = h_MC_shape_down->GetBinContent(i+1);
    		den_shape = denominator_shape->GetBinContent(i+1);;

				//std::cout<< "Setting point " << i << " to value 1 with errors +" <<num_shape_up/den_shape<<"/-"<<num_shape_down/den_shape<<std::endl;
				//std::cout<<"Value num_up "<<num_shape_up<<" and down "<<num_shape_down<<" over den "<<den_shape<<std::endl;
    		shape_syst->SetPoint(i, shape_syst->GetX()[i], 1);
    		//shape_syst->SetPoint(i, gr_systUnc->GetX()[i], 1);
    		shape_syst->SetPointEYlow(i, num_shape_down/den_shape); //have to look at the difference between nominal (1) and error
    		shape_syst->SetPointEYhigh(i, num_shape_up/den_shape); //I could use fabs here... However for the graph below we don't want difference but the actual value... So I don't substract anything here!
  		}
    	shape_syst->SetFillStyle(3005);
    	shape_syst->SetFillColor(kBlue);
    	shape_syst->SetLineStyle(1);
    	shape_syst->SetLineColor(1);
			//shape_syst->Draw("C same");
			//for (int i = 0; i < shape_syst->GetN(); ++i) std::cout<< shape_syst->GetX()[i] << " and " << shape_syst->GetEYhigh()[i] << std::endl;
    	TGraph *gr_shape_syst_up = new TGraph(shape_syst->GetN(), shape_syst->GetX(), shape_syst->GetEYhigh());
    	TGraph *gr_shape_syst_down = new TGraph(shape_syst->GetN(), shape_syst->GetX(), shape_syst->GetEYlow());
    	gr_shape_syst_up->SetLineColor(kRed);
    	gr_shape_syst_down->SetLineColor(kRed);
    	gr_shape_syst_up->Draw("L same"); //C, L, or CF
    	gr_shape_syst_down->Draw("L same");

 			TLegend *t_fit = new TLegend(0.78,0.78,0.95,0.95);
  		t_fit->SetLineColor(0);
  		t_fit->SetBorderSize(1);
  		t_fit->AddEntry(gr_shape_syst_up,"Syst. (Shape)","l");
  		t_fit->Draw();


		}
		else if(showSystShape && SystNormThenSystShapeThenStat){
  		TGraphAsymmErrors* numerator_systNandS = (TGraphAsymmErrors*) gr_systNandSUnc->Clone();

  		double num_systNandS =0., num_systNandS_up=0., num_systNandS_down=0., den_systNandS=0.;
  		for (int i = 0; i < numerator_systNandS->GetN(); ++i) {
    		num_systNandS = numerator_systNandS->GetY()[i];
    		num_systNandS_up = numerator_systNandS->GetEYhigh()[i];
    		num_systNandS_down = numerator_systNandS->GetEYlow()[i];
    		den_systNandS = numerator_systNandS->GetY()[i];;

    		numerator_systNandS->SetPoint(i, numerator_systNandS->GetX()[i], num_systNandS/den_systNandS);
    		numerator_systNandS->SetPointEYlow(i, num_systNandS_down/den_systNandS);
    		numerator_systNandS->SetPointEYhigh(i, num_systNandS_up/den_systNandS);
  		}

  		numerator_systNandS->Draw("2 same");


		}
	}
	pad2->Modified();
	pad2->Update();

	c0->Print(categoryName+".png");
	c0->Print(categoryName+".pdf");
	c0->Print(categoryName+".root");

}


void makeThePlotForOneCat(TString categoryName, TFile *theFile, TString theLeptonCategoryText, TString theJetCategoryText, bool doRatio, bool showSystShape, bool doFit, bool produceControlPlots, bool removeInstrMETStatUnc, bool showSystOnly, bool StatThenSyst, bool SystNormThenSystShapeThenStat, bool mergeSyst, TFile* systFile, bool doErrorPreFit, TString errorFitType, bool takePostNormalization, TString normaFitType, double sigma_signal){
	//Global syst from datacards (BY HAND ! In a perfect world should be done by reading a datacard... World is not perfect.)
	double CMS_eff_e = 1.072124; //For ee
	double CMS_eff_m = 1.061788; //For mumu
	double CMS_eff = 1.;
	if(categoryName.BeginsWith("ee"))        CMS_eff = CMS_eff_e;
	else if(categoryName.BeginsWith("mumu")) CMS_eff = CMS_eff_m;

  double CMS_hzz2l2v_sys_topwww_13TeV = 1.13;

  double QCDscale_WZ_0 = 1.095; //For 0jets
  double QCDscale_WZ_1 = 1.051; //For 1jets
  double QCDscale_WZ_VBF = 1.4; //For VBF
  double QCDscale_WZ = 1.;
  if(categoryName.EndsWith("eq0jets"))       QCDscale_WZ = QCDscale_WZ_0;
  else if(categoryName.EndsWith("geq1jets")) QCDscale_WZ = QCDscale_WZ_1;
  else if(categoryName.EndsWith("vbf"))      QCDscale_WZ = QCDscale_WZ_VBF;

	double QCDscale_ZZ_0 = 1.063; //For 0jets
  double QCDscale_ZZ_1 = 1.054; //For 1jets
  double QCDscale_ZZ_VBF = 1.4; //For VBF
  double QCDscale_ZZ = 1.;
  if(categoryName.EndsWith("eq0jets"))       QCDscale_ZZ = QCDscale_ZZ_0;
  else if(categoryName.EndsWith("geq1jets")) QCDscale_ZZ = QCDscale_ZZ_1;
  else if(categoryName.EndsWith("vbf"))      QCDscale_ZZ = QCDscale_ZZ_VBF;

	double QCDscale_ggZZ = 1.1;
  double kfactor_ggZZ = 1.1;

  double lumi_13TeV = 1.026;

	//All initialisation
	std::vector<TString> v_listMCProcesses   = {"ggH_bonl",          "zz", "wz", "zvv",   "instrmet",   "topwww"}; //The process names in the root file
	std::vector<int>     v_colorsMCProcesses = { 596,                 595,  594,   869,          831,          8}; //The colors associated to the above processes
	std::vector<TString> v_legendNames       = {"gg #rightarrow ( H #rightarrow ) ZZ", "qq #rightarrow ZZ", "WZ", "ZVV", "Instr. MET", "Top/W/WW"}; //The legend names associated to above
	std::map<TString, std::vector<double> >  map_globalUnc;
	map_globalUnc["zz"] = (std::vector<double>) {CMS_eff, QCDscale_ZZ, lumi_13TeV}; 
	map_globalUnc["ggH_bonl"] = (std::vector<double>) {CMS_eff, QCDscale_ggZZ, kfactor_ggZZ, lumi_13TeV}; 
	map_globalUnc["wz"] = (std::vector<double>) {CMS_eff, QCDscale_WZ, lumi_13TeV}; 
	map_globalUnc["zvv"] = (std::vector<double>) {CMS_eff, lumi_13TeV}; 
	map_globalUnc["instrmet"] = (std::vector<double>) {1.0}; 
	map_globalUnc[""] = (std::vector<double>) {CMS_hzz2l2v_sys_topwww_13TeV}; 

  //Get all the MC histos for nominal values
	std::map<TString, TH1D*> map_MCProcesses;
 	for(unsigned int i =0; i < v_listMCProcesses.size(); i++) map_MCProcesses[v_listMCProcesses[i]] = (TH1D*) ((TH1D*) theFile->Get(categoryName+"/"+v_listMCProcesses[i]))->Clone();
  	//Get ggH and qqH signal+interference for nominal values
  	TH1D *ggH_bonl = (TH1D*) ((TH1D*) theFile->Get(categoryName+"/ggH_bonl"))->Clone();
  	TH1D *qqH_bonl = (TH1D*) ((TH1D*) theFile->Get(categoryName+"/qqH_bonl"))->Clone();

  	TH1D *ggH_sand = (TH1D*) ((TH1D*) theFile->Get(categoryName+"/ggH_sand"))->Clone();
  	TH1D *qqH_sand = (TH1D*) ((TH1D*) theFile->Get(categoryName+"/qqH_sand"))->Clone();

  	TH1D *ggH_sonl = (TH1D*) ((TH1D*) theFile->Get(categoryName+"/ggH_sonl"))->Clone();
  	TH1D *qqH_sonl = (TH1D*) ((TH1D*) theFile->Get(categoryName+"/qqH_sonl"))->Clone();


		double mu_ggH = sigma_signal*0.0673*0.2*2/(1000.*0.00723416);
		double mu_qqH = sigma_signal*0.0673*0.2*2/(1000.*0.00253221);

		TH1D *ggH_Interf = (TH1D*) ggH_sand->Clone("ggH_Interf");
  	ggH_Interf->Add(ggH_sonl,-1);
  	ggH_Interf->Add(ggH_bonl,-1);
  	TH1D *ggH_sigInterf = (TH1D*) ggH_sonl->Clone("ggH_sigInterf");
  	ggH_sigInterf->Scale(mu_ggH);
    ggH_sigInterf->Add(ggH_Interf, -1*sqrt(mu_ggH));

		TH1D *qqH_Interf = (TH1D*) qqH_sand->Clone("qqH_Interf");
  	qqH_Interf->Add(qqH_sonl,-1);
  	qqH_Interf->Add(qqH_bonl,-1);
  	TH1D *qqH_sigInterf = (TH1D*) qqH_sonl->Clone("qqH_sigInterf");
  	qqH_sigInterf->Scale(mu_qqH);
    qqH_sigInterf->Add(qqH_Interf, -1*sqrt(mu_qqH));


		//Get the data for nominal value
  	TH1D *data_obs = (TH1D*) ((TH1D*) theFile->Get(categoryName+"/data_obs"))->Clone();
    TGraphAsymmErrors *theGraph = convertTheHistoInTGraph(data_obs);

		if(takePostNormalization){
    //Since we need to keep the bin width, the easiest is to use the histo defined above but to replace the bin content with the one from syst
    for(unsigned int i =0; i < v_listMCProcesses.size(); i++){
      for(unsigned int j = 1; j <=  map_MCProcesses[v_listMCProcesses[i]]->GetNbinsX(); j++){
    	map_MCProcesses[v_listMCProcesses[i]]->SetBinContent(j, ((TH1D*) systFile->Get(normaFitType+"/"+categoryName+"/"+v_listMCProcesses[i]))->GetBinContent(j));
			}
		}
  	//Get ggH and qqH signal+interference for nominal values
/*
		ggH_bonl = (TH1D*) ((TH1D*) systFile->Get(normaFitType+"/"+categoryName+"/ggH_bonl"))->Clone();
  	qqH_bonl = (TH1D*) ((TH1D*) systFile->Get(normaFitType+"/"+categoryName+"/qqH_bonl"))->Clone();

  	ggH_sand = (TH1D*) ((TH1D*) systFile->Get(normaFitType+"/"+categoryName+"/ggH_sand"))->Clone();
  	qqH_sand = (TH1D*) ((TH1D*) systFile->Get(normaFitType+"/"+categoryName+"/qqH_sand"))->Clone();

  	ggH_sigInterf = (TH1D*) ggH_sand->Clone("ggH_sigInterf");
  	*ggH_sigInterf = *ggH_sand;
  	ggH_sigInterf->Add(ggH_bonl,-1);
  	qqH_sigInterf = (TH1D*) qqH_sand->Clone("ggH_sigInterf");
  	*qqH_sigInterf = *qqH_sand;
  	qqH_sigInterf->Add(qqH_bonl,-1);

		//Get the data for nominal value
  	theGraph = (TGraphAsymmErrors*) ((TGraphAsymmErrors*) systFile->Get(normaFitType+"/"+categoryName+"/data"))->Clone();
*/
		}

  //Stack the histos for nominal value
  THStack *hs = new THStack("hs","Stacked of the backgrounds");

	for(unsigned int i =0; i < v_listMCProcesses.size(); i++){
    map_MCProcesses[v_listMCProcesses[i]]->SetFillColor(v_colorsMCProcesses[i]);
	  hs->Add(map_MCProcesses[v_listMCProcesses[i]]);
	}

  //Take all stat histo for uncertainty bands
	std::map<TString, std::vector<TH1D*> > map_MC_stats_up, map_MC_stats_down;
  theFile->cd(categoryName);
  TIter nextkey(gDirectory->GetListOfKeys());
  TKey *key;
  while ((key=(TKey*)nextkey())) {
    TObject *obj = key->ReadObj();
		for (unsigned int i= 0; i < v_listMCProcesses.size(); i++){
			if( ((TString) key->GetName()).BeginsWith(v_listMCProcesses[i]+"_CMS_hzz2l2v_stat_")){
				TH1D *h_tmp = (TH1D*) ((TH1D*) theFile->Get(categoryName+"/"+key->GetName()))->Clone();
				h_tmp->Add(map_MCProcesses[v_listMCProcesses[i]],-1);
				if( ((TString) key->GetName()).EndsWith("Up")) map_MC_stats_up[v_listMCProcesses[i]].push_back((TH1D*) h_tmp->Clone());
				if( ((TString) key->GetName()).EndsWith("Down")) map_MC_stats_down[v_listMCProcesses[i]].push_back((TH1D*) h_tmp->Clone());
			}
		}
    delete obj;      
  }
  //Create the TGraphAsymmErrors coming from stat
	std::map<TString, TH1D*> map_MC_TotalStat_up, map_MC_TotalStat_down; 
  TH1D* h_MC_stat_down = (TH1D*) ((TH1D*) theFile->Get(categoryName+"/total"))->Clone();
  TH1D* h_MC_stat_up = (TH1D*) ((TH1D*) theFile->Get(categoryName+"/total"))->Clone();
  for(unsigned int j =0; j < v_listMCProcesses.size(); j++){
    for( unsigned int i = 0; i < map_MC_stats_up[v_listMCProcesses[j]].size(); i++){
    	if(map_MC_TotalStat_up[v_listMCProcesses[j]]) map_MC_TotalStat_up[v_listMCProcesses[j]]->Add(map_MC_stats_up[v_listMCProcesses[j]].at(i));
    	else map_MC_TotalStat_up[v_listMCProcesses[j]] = (TH1D*) map_MC_stats_up[v_listMCProcesses[j]].at(i)->Clone();
    }
    for( unsigned int i = 0; i < map_MC_stats_down[v_listMCProcesses[j]].size(); i++){
    	if(map_MC_TotalStat_down[v_listMCProcesses[j]]) map_MC_TotalStat_down[v_listMCProcesses[j]]->Add(map_MC_stats_down[v_listMCProcesses[j]].at(i));
      else map_MC_TotalStat_down[v_listMCProcesses[j]] = (TH1D*) map_MC_stats_down[v_listMCProcesses[j]].at(i)->Clone();
		}
	}


 	double temp_binContent_up = 0;
  double temp_binContent_down = 0;
  for( unsigned int i = 1; i <= ((TH1D*) hs->GetStack()->Last())->GetNbinsX(); i++){
    temp_binContent_up = 0; 
    temp_binContent_down = 0; 
    for(unsigned int j =0; j < v_listMCProcesses.size(); j++){
    	if(showSystOnly) continue;
		  if(v_listMCProcesses[j] == "instrmet" && removeInstrMETStatUnc) continue;
		  //if(v_listMCProcesses[j] != "zvv") continue;
		  //if(v_listMCProcesses[j] != "instrmet") continue;
      temp_binContent_up += map_MC_TotalStat_up[v_listMCProcesses[j]]->GetBinContent(i)*map_MC_TotalStat_up[v_listMCProcesses[j]]->GetBinContent(i);
      temp_binContent_down += map_MC_TotalStat_down[v_listMCProcesses[j]]->GetBinContent(i)*map_MC_TotalStat_down[v_listMCProcesses[j]]->GetBinContent(i);
		}
    h_MC_stat_up->SetBinContent(i, sqrt(temp_binContent_up));
    h_MC_stat_down->SetBinContent(i, sqrt(temp_binContent_down));
	}
	std::vector<double> vx, vy, vexl, vexh, veyl, veyh;
  for(unsigned int i = 1; i <= ((TH1D*) hs->GetStack()->Last())->GetNbinsX(); i++){
    vx.push_back(((TH1D*) hs->GetStack()->Last())->GetBinCenter(i));
    vy.push_back(((TH1D*) hs->GetStack()->Last())->GetBinContent(i));
    vexl.push_back(vx[i-1]-((TH1D*) hs->GetStack()->Last())->GetBinLowEdge(i));
    vexh.push_back(((TH1D*) hs->GetStack()->Last())->GetBinLowEdge(i)+ ((TH1D*) hs->GetStack()->Last())->GetBinWidth(i)-vx[i-1]);
    veyl.push_back(h_MC_stat_down->GetBinContent(i));
    veyh.push_back(h_MC_stat_up->GetBinContent(i));
    //Warning: TGraphAsymmErrors need DEVIATION from central values (i.e. errors) not the lower and upper values
  }






	//Syst
	std::vector<double> systx = (std::vector<double>) vx;
	std::vector<double> systy = (std::vector<double>) vy;
	std::vector<double> systexl = (std::vector<double>) vexl;
	std::vector<double> systexh = (std::vector<double>) vexh;
	std::vector<double> systeyl = (std::vector<double>) veyl;
	std::vector<double> systeyh = (std::vector<double>) veyh;

	double temp_globalSystMC_square = 0.;
	for(unsigned int i = 1; i <= ((TH1D*) hs->GetStack()->Last())->GetNbinsX(); i++){
	  temp_globalSystMC_square = 0.;

		for(unsigned int j = 0; j < v_listMCProcesses.size(); j++){
			for(unsigned int k =0; k < map_globalUnc[v_listMCProcesses[j]].size(); k++){
				temp_globalSystMC_square += (map_globalUnc[v_listMCProcesses[j]][k]-1)*map_MCProcesses[v_listMCProcesses[j]]->GetBinContent(i) * (map_globalUnc[v_listMCProcesses[j]][k]-1)*map_MCProcesses[v_listMCProcesses[j]]->GetBinContent(i);
			}
		}

    if(StatThenSyst && !SystNormThenSystShapeThenStat){
		  systeyl[i-1] = sqrt(systeyl[i-1]*systeyl[i-1]+temp_globalSystMC_square);
	    systeyh[i-1] = sqrt(systeyh[i-1]*systeyh[i-1]+temp_globalSystMC_square);
		}
		else{
		  systeyl[i-1] = sqrt(temp_globalSystMC_square); //take syst only
	    systeyh[i-1] = sqrt(temp_globalSystMC_square);
		  veyl[i-1] = sqrt(veyl[i-1]*veyl[i-1]+temp_globalSystMC_square); //take stat like syst+stat
		  veyh[i-1] = sqrt(veyh[i-1]*veyh[i-1]+temp_globalSystMC_square);
		}
	}


	//Handle shape syst
  std::map<TString, std::vector<TH1D*> > map_MC_shape_up, map_MC_shape_down;
  std::map<TString, std::vector<TH1D*> > map_MC_shape_up_fit, map_MC_shape_down_fit;
	std::map<TString, std::vector<TString> > map_shapeLegendName_up, map_shapeLegendName_down;
  theFile->cd(categoryName);
  TIter nextkey2(gDirectory->GetListOfKeys());
  TKey *key2;
  while ((key2=(TKey*)nextkey2())) {
    TObject *obj = key2->ReadObj();
    for (unsigned int i= 0; i < v_listMCProcesses.size(); i++){
      //if( v_listMCProcesses[i] == "ggH_bonl" ) continue;
      //if( v_listMCProcesses[i] == "wz" ) continue;
      //if( v_listMCProcesses[i] == "zz" ) continue;
      //if( v_listMCProcesses[i] == "zvv" ) continue;
      if( ((TString) key2->GetName()).BeginsWith(v_listMCProcesses[i]+"_CMS_hzz2l2v_")){
      	if( ((TString) key2->GetName()).BeginsWith(v_listMCProcesses[i]+"_CMS_hzz2l2v_stat_")) continue; //loop on syst and not stat
        TH1D *h_tmp = (TH1D*) ((TH1D*) theFile->Get(categoryName+"/"+key2->GetName()))->Clone();
        //if(((TString) key2->GetName()).Contains("th_pdf")) continue;
        //if(((TString) key2->GetName()).Contains("th_ewk")) continue;
        //if(!((TString) key2->GetName()).Contains("zllinstrmet")) continue;
        //if(((TString) key2->GetName()).Contains("scale_j")) continue;
        //if(((TString) key2->GetName()).Contains("th_alphas")) continue;
        if(doFit){
        	if( ((TString) key2->GetName()).EndsWith("Up")) map_shapeLegendName_up[v_listMCProcesses[i]].push_back(key2->GetName());
        	if( ((TString) key2->GetName()).EndsWith("Down")) map_shapeLegendName_down[v_listMCProcesses[i]].push_back(key2->GetName());
        	h_tmp->Scale(1.*map_MCProcesses[v_listMCProcesses[i]]->Integral()/h_tmp->Integral());
          if( ((TString) key2->GetName()).EndsWith("Up")) map_MC_shape_up_fit[v_listMCProcesses[i]].push_back((TH1D*) h_tmp->Clone());
          if( ((TString) key2->GetName()).EndsWith("Down")) map_MC_shape_down_fit[v_listMCProcesses[i]].push_back((TH1D*) h_tmp->Clone());
        }
        h_tmp->Add(map_MCProcesses[v_listMCProcesses[i]],-1);
        if( ((TString) key2->GetName()).EndsWith("Up")) map_MC_shape_up[v_listMCProcesses[i]].push_back((TH1D*) h_tmp->Clone());
        if( ((TString) key2->GetName()).EndsWith("Down")) map_MC_shape_down[v_listMCProcesses[i]].push_back((TH1D*) h_tmp->Clone());
      }
    }
    delete obj;
  }
  //Create the TH1 coming from shape syst
  TH1D* h_MC_shape_down = (TH1D*) ((TH1D*) theFile->Get(categoryName+"/total"))->Clone();
  TH1D* h_MC_shape_up = (TH1D*) ((TH1D*) theFile->Get(categoryName+"/total"))->Clone();
  temp_binContent_up = 0;
  temp_binContent_down = 0;
  for( unsigned int i = 1; i <= ((TH1D*) hs->GetStack()->Last())->GetNbinsX(); i++){
    temp_binContent_up = 0;
    temp_binContent_down = 0;
    for(unsigned int j =0; j < v_listMCProcesses.size(); j++){
      //if(v_listMCProcesses[j] != "instrmet") continue;
      for(unsigned int k =0; k < map_MC_shape_up[v_listMCProcesses[j]].size(); k++){
        temp_binContent_up += map_MC_shape_up[v_listMCProcesses[j]][k]->GetBinContent(i)*map_MC_shape_up[v_listMCProcesses[j]][k]->GetBinContent(i);
        temp_binContent_down += map_MC_shape_down[v_listMCProcesses[j]][k]->GetBinContent(i)*map_MC_shape_down[v_listMCProcesses[j]][k]->GetBinContent(i);
      }
    }

		if(SystNormThenSystShapeThenStat){
      //Syst Norm already defined above
      h_MC_shape_down->SetBinContent(i, sqrt(systeyl[i-1]*systeyl[i-1]+temp_binContent_down)); //Take syst Norm + shape, and this time only take the events due to shape diff, not nominal + diff
      h_MC_shape_up->SetBinContent(i, sqrt(systeyh[i-1]*systeyh[i-1]+temp_binContent_up)); 
      veyl[i-1] = sqrt(veyl[i-1]*veyl[i-1]+temp_binContent_down); //take stat like systN+systS+stat
      veyh[i-1] = sqrt(veyh[i-1]*veyh[i-1]+temp_binContent_up);
		}
		else{
			h_MC_shape_up->SetBinContent(i, ((TH1D*) hs->GetStack()->Last())->GetBinContent(i)+sqrt(temp_binContent_up));
      h_MC_shape_down->SetBinContent(i, ((TH1D*) hs->GetStack()->Last())->GetBinContent(i)-sqrt(temp_binContent_down));
 		}
 	}


  if(doFit && produceControlPlots){ //produce control plots
  	int W = 800;
  	int H = 600;
  	int H_ref = 600;
  	int W_ref = 800;

  	// references for T, B, L, R
  	float T = 0.08*H_ref;
  	float B = 0.12*H_ref;
  	float L = 0.12*W_ref;
  	float R = 0.04*W_ref;

    std::map<TString, TCanvas*> map_c_fit;
    for(unsigned int j =0; j < v_listMCProcesses.size(); j++){
  		map_c_fit[v_listMCProcesses[j]] = new TCanvas("c_fit",categoryName+"_"+v_listMCProcesses[j]+"_fit",50,50,W,H);
  		map_c_fit[v_listMCProcesses[j]]->SetFillColor(0);
  		map_c_fit[v_listMCProcesses[j]]->SetBorderMode(0);
  		map_c_fit[v_listMCProcesses[j]]->SetFrameFillStyle(0);
  		map_c_fit[v_listMCProcesses[j]]->SetFrameBorderMode(0);
  		map_c_fit[v_listMCProcesses[j]]->SetLeftMargin( L/W );
  		map_c_fit[v_listMCProcesses[j]]->SetRightMargin( R/W );
  		map_c_fit[v_listMCProcesses[j]]->SetTopMargin( T/H );
  		map_c_fit[v_listMCProcesses[j]]->SetBottomMargin( B/H );
  		map_c_fit[v_listMCProcesses[j]]->SetTickx(1);
  		map_c_fit[v_listMCProcesses[j]]->SetTicky(1);
  		map_c_fit[v_listMCProcesses[j]]->SetLogy();  

  		//Fix the size of the canvas: xmin,ymin,xmax,ymax
  		gPad->DrawFrame(0.0, 0.001, 3000.0, 1000);
  		gPad->Update();

      map_MCProcesses[v_listMCProcesses[j]]->Draw("hist same");
      //Double_t colors[12] = {632, 416, 600, 400, 616, 432, 800, 820, 840, 860, 880, 900};
      Double_t colors[13] = {2,1,3,4,5,6,7,8,9,808,606,921, 402};
  		for(unsigned int k =0; k < map_MC_shape_up_fit[v_listMCProcesses[j]].size(); k++){
        map_MC_shape_up_fit[v_listMCProcesses[j]][k]->SetLineColor(colors[k]);
        map_MC_shape_down_fit[v_listMCProcesses[j]][k]->SetLineColor(colors[k]);
        //map_MC_shape_up_fit[v_listMCProcesses[j]][k]->SetLineColor(k+1);
        //map_MC_shape_down_fit[v_listMCProcesses[j]][k]->SetLineColor(k+1);
        map_MC_shape_up_fit[v_listMCProcesses[j]][k]->Draw("same hist");
				map_MC_shape_down_fit[v_listMCProcesses[j]][k]->Draw("same hist");

  		}
  		TLegend *t_fit = new TLegend(0.38,0.58,0.95,0.85);
  		t_fit->SetLineColor(0);
  		t_fit->SetBorderSize(1);
  		t_fit->SetNColumns(3);
  		for(unsigned int k =0; k < map_MC_shape_up_fit[v_listMCProcesses[j]].size(); k++){
  			t_fit->AddEntry(map_MC_shape_up_fit[v_listMCProcesses[j]][k],(TString) ((TString) map_shapeLegendName_up[v_listMCProcesses[j]][k].Remove(map_shapeLegendName_up[v_listMCProcesses[j]][k].Length()-2,map_shapeLegendName_up[v_listMCProcesses[j]][k].Length()))(map_shapeLegendName_up[v_listMCProcesses[j]][k].Index("2l2v_")+5, map_shapeLegendName_up[v_listMCProcesses[j]][k].Length()),"l");
  		}
  		t_fit->AddEntry(map_MCProcesses[v_listMCProcesses[j]],v_legendNames[j],"f");
  		t_fit->Draw();
  		gPad->RedrawAxis();

  		map_c_fit[v_listMCProcesses[j]]->Print("shape/"+categoryName+"_"+v_listMCProcesses[j]+"_fit.pdf");
  		map_c_fit[v_listMCProcesses[j]]->Print("shape/"+categoryName+"_"+v_listMCProcesses[j]+"_fit.png");

  	}

	}

	TGraphAsymmErrors* gr_statUnc = new TGraphAsymmErrors(vx.size(), &(vx[0]), &(vy[0]), &(vexl[0]), &(vexh[0]), &(veyl[0]), &(veyh[0]));
	TGraphAsymmErrors* gr_systUnc = new TGraphAsymmErrors(systx.size(), &(systx[0]), &(systy[0]), &(systexl[0]), &(systexh[0]), &(systeyl[0]), &(systeyh[0]));



  ///////////////////////////////////////////////////
  ///////////////   DRAWING SECTION   ///////////////
  ///////////////////////////////////////////////////
  if(doRatio) drawRatioPlots(theGraph, hs, ggH_sigInterf, qqH_sigInterf, gr_statUnc, gr_systUnc, v_listMCProcesses, theLeptonCategoryText, theJetCategoryText, categoryName, map_MCProcesses, v_legendNames, h_MC_shape_up, h_MC_shape_down, showSystShape, StatThenSyst, SystNormThenSystShapeThenStat, mergeSyst, systFile, doErrorPreFit, errorFitType, sigma_signal);
  else drawPlots(theGraph, hs, ggH_sigInterf, qqH_sigInterf, gr_statUnc, gr_systUnc, v_listMCProcesses, theLeptonCategoryText, theJetCategoryText, categoryName, map_MCProcesses, v_legendNames, StatThenSyst, SystNormThenSystShapeThenStat, h_MC_shape_up, h_MC_shape_down, mergeSyst, systFile, doErrorPreFit, sigma_signal);

}



void drawThePerCategoryPlot(){
  gROOT->SetBatch(); 
	gROOT->ForceStyle();
  setTDRStyle();
  TFile *histoFile = new TFile("hzz2l2v_800_13TeV.root");
  TFile *systFile = new TFile("mlFitForPapers/mlfit_allSyst.root");

	bool doErrorPreFit = false; //This option erases most of options below as it's taking errors from the pre/post fit
  //TString errorFitType = "shapes_prefit"; //Please select one of the three types of error to show
  TString errorFitType = "shapes_fit_b";
  //TString errorFitType = "shapes_fit_s"; //What is this shape exactly?
  bool takePostNormalization = false;
  //TString normaFitType = "shapes_prefit"; //Please select one of the three types of error to show
  TString normaFitType = "shapes_fit_b";
  //TString normaFitType = "shapes_fit_s"; //What is this shape exactly?
  //if(takePostNormalization) doErrorPreFit = true; //Only use errors from fit if normalization from fit
  
  //Do not touch - begin
  bool doRatio = true;
  bool showSystShape = true;
  bool doFit = true;
  bool produceControlPlots = false;
  //if(!doRatio) showSystShape = false;
  if(!showSystShape) doFit = false;
  bool removeInstrMETStatUnc = false;
  bool showSystOnly = false;
  bool StatThenSyst = false; //show stat+syst if true and syst+stat if false
  bool SystNormThenSystShapeThenStat = true;
	if(StatThenSyst) SystNormThenSystShapeThenStat = false;
	bool mergeSyst = true;
	if(!SystNormThenSystShapeThenStat) mergeSyst = false;
  //Do not touch - end

  double sigma_signal = 50;

	makeThePlotForOneCat("eeeq0jets",histoFile, "ee", "0 jet", doRatio, showSystShape, doFit, produceControlPlots, removeInstrMETStatUnc, showSystOnly, StatThenSyst, SystNormThenSystShapeThenStat, mergeSyst, systFile, doErrorPreFit, errorFitType, takePostNormalization, normaFitType, sigma_signal); //First argument is super important!
  makeThePlotForOneCat("eegeq1jets",histoFile, "ee", "#geq1 jet", doRatio, showSystShape, doFit, produceControlPlots, removeInstrMETStatUnc, showSystOnly, StatThenSyst, SystNormThenSystShapeThenStat, mergeSyst, systFile, doErrorPreFit, errorFitType, takePostNormalization, normaFitType, sigma_signal);
  makeThePlotForOneCat("eevbf",histoFile, "ee", "VBF-tagged", doRatio, showSystShape, doFit, produceControlPlots, removeInstrMETStatUnc, showSystOnly, StatThenSyst, SystNormThenSystShapeThenStat, mergeSyst, systFile, doErrorPreFit, errorFitType, takePostNormalization, normaFitType, sigma_signal);
  makeThePlotForOneCat("mumueq0jets",histoFile, "#mu#mu", "0 jet", doRatio, showSystShape, doFit, produceControlPlots, removeInstrMETStatUnc, showSystOnly, StatThenSyst, SystNormThenSystShapeThenStat, mergeSyst, systFile, doErrorPreFit, errorFitType, takePostNormalization, normaFitType, sigma_signal);
  makeThePlotForOneCat("mumugeq1jets",histoFile, "#mu#mu", "#geq1 jet", doRatio, showSystShape, doFit, produceControlPlots, removeInstrMETStatUnc, showSystOnly, StatThenSyst, SystNormThenSystShapeThenStat, mergeSyst, systFile, doErrorPreFit, errorFitType, takePostNormalization, normaFitType, sigma_signal);
  makeThePlotForOneCat("mumuvbf",histoFile, "#mu#mu", "VBF-tagged", doRatio, showSystShape, doFit, produceControlPlots, removeInstrMETStatUnc, showSystOnly, StatThenSyst, SystNormThenSystShapeThenStat, mergeSyst, systFile, doErrorPreFit, errorFitType, takePostNormalization, normaFitType, sigma_signal);
}
