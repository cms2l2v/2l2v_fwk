//To execute : .X MakeBeauty_custom_forMT.C

#include <TH1F.h>
#include <TH2F.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TRandom.h>
#include <TMath.h>
#include <math.h>
#include <vector>
#include <TFile.h>
#include <iostream>
#include "TROOT.h"
#include <TGaxis.h>
#include <TLine.h>
#include <algorithm> 
#include <map>

//You can tune easily the options here
#define VERBOSE false
#define HIDE_WARNING true
//Choose which syst you want in the end. ALL is a global switch that overwrites the others if set to TRUE.
#define ALL true
#define WGAMMA_ONLY true
#define ZGAMMA_ONLY true
#define WJETS_ONLY true
#define METHOD_ONLY true
#define GAMMASTATS_ONLY true
//Set to false if InstrMET not in root file
#define InstrMET_control_plots true

void MakeSyst_custom(TFile* f_input, TFile* f_output, TString SystType, int projectionBin_arg){
	if(VERBOSE)std::cout<< "Starting macro" << std::endl;
  gROOT->SetBatch();  

	if(HIDE_WARNING) gErrorIgnoreLevel=kError;

	int projectionBin;
	//bool doProjection=false; //tell if the histo you are looking at needs a projection or not
	//if(SystType == "met" || SystType == "mt") doProjection=true;
	
	//Don't touch this, because the name of the final plot produced is metfinal or mtfinal... and those are with 125GeV MET cut!
	if(SystType == "met") projectionBin = projectionBin_arg; //set a cut at 1 for a 50GeV MET cut
	if(SystType == "mt") projectionBin = projectionBin_arg; //set a cut at 17 for 125GeV MET cut


	TString XaxisName;
	if(SystType == "met") XaxisName = "E_{T}^{miss} [GeV]";
	if(SystType == "mt") XaxisName = "M_{T} [GeV]";

	TString globalDirectoryName= SystType;
  TDirectory* globalDirectory = f_output->mkdir(globalDirectoryName+"/");
  globalDirectory->cd();

	//Def of the histo of interest
	std::vector<TString> chTags;
	chTags.push_back("ee"); 	
  chTags.push_back("mumu"); 
  //chTags.push_back("ll"); 
  //chTags.push_back("emu");
	//chTags.push_back("gamma"); 

	std::vector<TString> evCat;
	evCat.push_back("eq0jets");	
	evCat.push_back("vbf");
	evCat.push_back("geq1jets");
	evCat.push_back(""); //all cat

  //Definition of my histograms
	std::vector<TH2F*> met_shapes, met_shapesUp, met_shapesDown; 
	std::vector<TH2F*> met_shapes_fromTheory, met_shapesUp_fromTheory, met_shapesDown_fromTheory; 

  //My canvas for the output
  std::vector<TCanvas*> c;
  std::vector<TCanvas*> c_fromTheory;
  std::vector<TCanvas*> c_wrapUp;
  std::vector<TCanvas*> c_final_wrapUp;
  std::vector<TCanvas*> c_InstrMET_genuineMET;

	//My different types of correction
	std::vector< TString> correctionsUp;
	correctionsUp.push_back("_th_factup");
	correctionsUp.push_back("_res_jup"); 
	correctionsUp.push_back("_scale_jup"); 
	correctionsUp.push_back("_scale_mup"); 
	correctionsUp.push_back("_scale_eup"); 
	correctionsUp.push_back("_puup"); 
	correctionsUp.push_back("_scale_umetup"); 
	correctionsUp.push_back("_eff_bup"); 

	std::vector< TString> correctionsDown;
	correctionsDown.push_back("_th_factdown");
	correctionsDown.push_back("_res_jdown"); 
	correctionsDown.push_back("_scale_jdown"); 
	correctionsDown.push_back("_scale_mdown"); 
	correctionsDown.push_back("_scale_edown"); 
	correctionsDown.push_back("_pudown"); 
	correctionsDown.push_back("_scale_umetdown"); 
	correctionsDown.push_back("_eff_bdown"); 

	std::vector< TString> samples;
	samples.push_back("W#gamma #rightarrow l#nu#gamma_reweighted");
	samples.push_back("Z#gamma #rightarrow #nu#nu#gamma_reweighted"); 
	samples.push_back("W#rightarrow l#nu_reweighted");


	// create a subdirectory in this file
	std::vector< TDirectory*> directory;
	for(unsigned int s =0; s < samples.size(); s++){  
		for ( unsigned int ich = 0; ich < chTags.size() ; ich++){
			globalDirectory->cd();
   		directory.push_back(f_output->mkdir(globalDirectoryName+"/"+samples[s]+"_"+chTags[ich]));
		}
	}
	int position =0;
	int cpos=0;
	int emptyPosition=0;
	int position_canvas_reset =0;
	int directory_pos=0;
	if(VERBOSE)std::cout<< "Let s go in the loops" << std::endl;



	std::vector<Double_t> v_custom_axis;
	if(SystType == "met") v_custom_axis ={0,5,10,15,20,25,30,35,40,45,50,55,60,70,80,90,100,125,150,175,200,250,300,400,500,600,700,800,900,1000};
	//else if(SystType == "mt") v_custom_axis ={150, 225, 300, 375, 450, 525, 600, 725, 850, 975, 1100, 1350, 1600, 1850, 2100, 2600, 3000}; //from https://github.com/cms2l2v/2l2v_fwk/blob/master/bin/common/computeLimit.cc#L2587
	else if(SystType == "mt") v_custom_axis ={100,120,140,160,180,200,220,240,260,280,300,325,350,375,400,450,500,600,700,800,900,1000,1500,2000}; //from https://github.com/cms2l2v/2l2v_fwk/blob/master/bin/hzz2l2v/runHZZ2l2vAnalysis.cc#L564
	//else if(SystType == "mt") v_custom_axis = {150, 225, 300, 375, 450, 525, 600, 725, 850, 975, 1100, 1350, 1600, 1850, 2100, 2600, 3000}; //from limits: https://github.com/cms2l2v/2l2v_fwk/blob/master/bin/common/computeLimit.cc#L2587
	else{
		std::cout<<"Error: change the mode" << std::endl;
	}
	Double_t* metaxis = new Double_t[v_custom_axis.size()];
	std::copy(v_custom_axis.begin(), v_custom_axis.end(), metaxis);
	Int_t nmetAxis=v_custom_axis.size();

	//The same axis than the one from limits: https://github.com/cms2l2v/2l2v_fwk/blob/master/bin/common/computeLimit.cc#L2587
	Double_t eq0jets_axis_limits[] = {150, 225, 300, 375, 450, 525, 600, 725, 850, 975, 1100, 1350, 1600, 2100, 3000};
	Double_t geq1jets_axis_limits[] = {150, 225, 300, 375, 450, 525, 600, 725, 850, 975, 1100, 1350, 1600, 2100, 3000};
	Double_t vbf_axis_limits[] = {150, 225, 300, 375, 450, 525, 600, 725,  1100, 3000};
	Double_t all_axis_limits[] = {150, 225, 300, 375, 450, 525, 600, 725, 850, 975, 1100, 1350, 1600, 2100, 3000};
	std::map<TString, Double_t* > metaxis_limits;
	metaxis_limits["eq0jets"] = eq0jets_axis_limits;
	metaxis_limits["geq1jets"] = geq1jets_axis_limits;
	metaxis_limits["vbf"] = vbf_axis_limits;
	metaxis_limits[""] = all_axis_limits;
	std::map<TString, Int_t> nmetAxis_limits = {
		        { "eq0jets", sizeof(eq0jets_axis_limits)/sizeof(Double_t) },
		        { "geq1jets", sizeof(geq1jets_axis_limits)/sizeof(Double_t) },
		        { "vbf", sizeof(vbf_axis_limits)/sizeof(Double_t) },
		        { "", sizeof(all_axis_limits)/sizeof(Double_t) } };

	//Initialize a place to save histo:
	std::vector<std::vector<std::vector<std::vector<TH1D*> > >> savedHisto(samples.size(), std::vector<std::vector<std::vector<TH1D*> > >(chTags.size(), std::vector<std::vector<TH1D*> >(evCat.size()/*, std::vector<TH1D*>*/))) ;
	std::vector<std::vector<std::vector<std::vector<TH1D*> > >> savedHisto_Up(samples.size(), std::vector<std::vector<std::vector<TH1D*> > >(chTags.size(), std::vector<std::vector<TH1D*> >(evCat.size()/*, std::vector<TH1D*>*/))) ;
	std::vector<std::vector<std::vector<std::vector<TH1D*> > >> savedHisto_Down(samples.size(), std::vector<std::vector<std::vector<TH1D*> > >(chTags.size(), std::vector<std::vector<TH1D*> >(evCat.size()/*, std::vector<TH1D*>*/))) ;

	//Here are the histo with the wrap up shapes (i.e. final shape up and down for all syst+stat on the three main genuineMET processes)
	std::vector<std::vector<std::vector<std::vector<TH1D*> > >> wrapUp_shape(samples.size(), std::vector<std::vector<std::vector<TH1D*> > >(chTags.size(), std::vector<std::vector<TH1D*> >(evCat.size()/*, std::vector<TH1D*>*/))) ;
	std::vector<std::vector<std::vector<std::vector<TH1D*> > >> wrapUp_shape_Up(samples.size(), std::vector<std::vector<std::vector<TH1D*> > >(chTags.size(), std::vector<std::vector<TH1D*> >(evCat.size()/*, std::vector<TH1D*>*/))) ;
	std::vector<std::vector<std::vector<std::vector<TH1D*> > >> wrapUp_shape_Down(samples.size(), std::vector<std::vector<std::vector<TH1D*> > >(chTags.size(), std::vector<std::vector<TH1D*> >(evCat.size()/*, std::vector<TH1D*>*/))) ;

	//The list of histo of shape uncertainty on the Instr.MET from genuineMET
	std::vector<std::vector<std::vector<std::vector<TH1D*> > >> savedInstrMET_shape(samples.size(), std::vector<std::vector<std::vector<TH1D*> > >(chTags.size(), std::vector<std::vector<TH1D*> >(evCat.size()/*, std::vector<TH1D*>*/))) ;
	std::vector<std::vector<std::vector<std::vector<TH1D*> > >> savedInstrMET_shape_Up(samples.size(), std::vector<std::vector<std::vector<TH1D*> > >(chTags.size(), std::vector<std::vector<TH1D*> >(evCat.size()/*, std::vector<TH1D*>*/))) ;
	std::vector<std::vector<std::vector<std::vector<TH1D*> > >> savedInstrMET_shape_Down(samples.size(), std::vector<std::vector<std::vector<TH1D*> > >(chTags.size(), std::vector<std::vector<TH1D*> >(evCat.size()/*, std::vector<TH1D*>*/))) ;

	//List of histo saved for instr MET global unc
	std::vector<std::vector<std::vector<TH1D*> > > saved_cst_method_InstrMET_shape_Up(chTags.size(), std::vector<std::vector<TH1D*> >(evCat.size()/*, std::vector<TH1D*>*/)) ;
	std::vector<std::vector<std::vector<TH1D*> > > saved_cst_method_InstrMET_shape_Down(chTags.size(), std::vector<std::vector<TH1D*> >(evCat.size()/*, std::vector<TH1D*>*/)) ;
	std::vector<std::vector<std::vector<TH1D*> > > saved_cst_stat_InstrMET_shape_Up(chTags.size(), std::vector<std::vector<TH1D*> >(evCat.size()/*, std::vector<TH1D*>*/)) ;
	std::vector<std::vector<std::vector<TH1D*> > > saved_cst_stat_InstrMET_shape_Down(chTags.size(), std::vector<std::vector<TH1D*> >(evCat.size()/*, std::vector<TH1D*>*/)) ;

	for(unsigned int s =0; s < samples.size(); s++){
		for ( unsigned int i = 0; i < correctionsUp.size(); i ++){
			for ( unsigned int ich = 0; ich < chTags.size() ; ich++){
				directory_pos = s*chTags.size() + ich;
				if(correctionsUp[i] == "_scale_mup" && chTags[ich] != "mumu"){
					emptyPosition += evCat.size();
					continue;
				}
				if(correctionsUp[i] == "_scale_eup" && chTags[ich] != "ee"){
					emptyPosition += evCat.size();
					continue;
				}
				for ( unsigned ev = 0; ev < evCat.size() ; ev++){
					TString tags_full=chTags[ich]+evCat[ev]; 
					position = s*correctionsUp.size()*chTags.size()*evCat.size() + i*chTags.size()*evCat.size() + ich*evCat.size() + ev - emptyPosition;
					TString num ="_";
					num += correctionsUp[i];

					if(VERBOSE)std::cout<< "define canvas" << std::endl;
					//Canvas
					c.push_back( new TCanvas(tags_full+"_met_shapes_corrections"+num, tags_full+"_met_shapes_corrections"+num, 800, 800));
					if(VERBOSE)std::cout<< "take histo" << std::endl;
					//Take the histos of interest
					if(VERBOSE)std::cout << "Sample: " << samples[s] << std::endl;
					f_input->cd(samples[s]);		//This argument gives the file on which we should take the shape. We have to loop on my samples
					//Nominal
					met_shapes.push_back( (TH2F*) gDirectory->Get(tags_full+"_"+SystType+"_shapes"));
					//Up
					met_shapesUp.push_back( (TH2F*) gDirectory->Get(tags_full+"_"+SystType+"_shapes"+correctionsUp[i]));
					//Down
					met_shapesDown.push_back( (TH2F*) gDirectory->Get(tags_full+"_"+SystType+"_shapes"+correctionsDown[i]));
					if(VERBOSE)std::cout<< "projection... " << tags_full << " and algo " << i << std::endl;
					if( !met_shapes[position]){
						if(VERBOSE)std::cout<< "Empty histogram, skip it" << std:: endl;
						continue;
					}

					if( !met_shapesUp[position] || !met_shapesDown[position]){
						if(VERBOSE)std::cout << "No corrections found, skip it" << std::endl;
						continue;
					}
					//Projection
					if(VERBOSE)std::cout << "Number of events = " << met_shapes[position]->GetEntries() << std::endl;
      		TH1D* Nominal_met = met_shapes[position]->ProjectionY("Nominal_met", projectionBin, projectionBin); //We do a projection of slice 17 (from slice 17 to 17 actually). This is related to a MET cut. The met cut is 50+(5*sliceNumber-2). Slice 1 corresponds to no MET cut,slide 17 to 125GeV MET cut.
      		TH1D* Up_met = met_shapesUp[position]->ProjectionY("Up_met", projectionBin, projectionBin);
      		TH1D* Down_met = met_shapesDown[position]->ProjectionY("Down_met", projectionBin, projectionBin);
					if(VERBOSE)std::cout<< "projection met done" << std::endl;


					//Save the values in a vector so we can loop on it at the end:
					if(VERBOSE)std::cout << "Number of events = " << Nominal_met->GetEntries() << std::endl;
					if(Nominal_met->GetEntries() ==0) continue;
					savedHisto[s][ich][ev].push_back((TH1D*) Nominal_met->Clone());
					savedHisto_Up[s][ich][ev].push_back((TH1D*) Up_met->Clone());
					savedHisto_Down[s][ich][ev].push_back((TH1D*) Down_met->Clone());

					TH1D *Nominal_met_new = (TH1D*) Nominal_met->Rebin(nmetAxis-1, "Nominal_met_new", metaxis);
					Nominal_met_new->Scale(1.0, "width");
					TH1D *Up_met_new = (TH1D*) Up_met->Rebin(nmetAxis-1, "Up_met_new", metaxis);
					Up_met_new->Scale(1.0, "width");
					TH1D *Down_met_new = (TH1D*) Down_met->Rebin(nmetAxis-1, "Down_met_new", metaxis);
					Down_met_new->Scale(1.0, "width");


					//Style part
					//Nominal in black
					Nominal_met_new->SetLineColor(kBlack);

					//Up in red
        	Up_met_new->SetLineColor(kRed);

					//Down in green
        	Down_met_new->SetLineColor(kGreen);

					//No stats
 					Nominal_met_new->SetStats(kFALSE);


					//Make the legends

        	TLegend *leg_met = new TLegend(0.6,0.7,0.89,0.89);
        	leg_met->AddEntry(Nominal_met_new, "Nominal", "l");
        	leg_met->AddEntry(Up_met_new, correctionsUp[i], "l");
        	leg_met->AddEntry(Down_met_new, correctionsDown[i], "l");






					if(VERBOSE)std::cout<< "ready to fill" << std::endl;
		  		cpos=position;
		  		c[cpos]->cd();
        	TPad *pad_met1 = new TPad("pad_met1","pad_met1",0,0.3,1,1);
        	TPad *pad_met2 = new TPad("pad_met2","pad_met2",0,0,1,0.3);
        	pad_met1->SetTopMargin(0.05);
        	pad_met1->SetBottomMargin(0);
        	pad_met1->Draw();
        	pad_met2->SetTopMargin(0);
        	pad_met2->SetBottomMargin(0.25);
        	pad_met2->Draw();
        	pad_met1->cd();
        	TH1 *h_nominal_met = Nominal_met_new->DrawCopy(); 
        	h_nominal_met->SetTitle("");
        	TH1 *h_up_met = Up_met_new->DrawCopy("same"); 
        	TH1 *h_down_met = Down_met_new->DrawCopy("same"); 
        	leg_met->Draw(); 

        	// Do not draw the Y axis label on the upper plot and redraw a small
		    	// axis instead, in order to avoid the first label (0) to be clipped.
   				//h_nominal_met->GetYaxis()->SetLabelSize(0.);
   				//TGaxis *axis = new TGaxis( -5, -5, 5, -5, 0.01,200,510,"");
   				//axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   				//axis->SetLabelSize(15);
   				//axis->Draw();

        	// Y axis ratio plot settings
        	h_nominal_met->GetYaxis()->SetTitle("Events (a.u.)");
        	h_nominal_met->GetYaxis()->SetNdivisions(505);
        	h_nominal_met->GetYaxis()->SetTitleSize(20);
        	h_nominal_met->GetYaxis()->SetTitleFont(43);
        	h_nominal_met->GetYaxis()->SetTitleOffset(1.55);
        	h_nominal_met->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
        	h_nominal_met->GetYaxis()->SetLabelSize(15); 

        	pad_met2->cd();
        	TH1D *h_nominal_met2 = (TH1D*)Up_met_new->Clone("h_nominal_met2");//Be careful: non intuitive naming conventions to earn time. 
        	TH1D *h_nominal_met3 = (TH1D*)Down_met_new->Clone("h_nominal_met3");
        	TH1D *h_up_met2 = (TH1D*)Nominal_met_new->Clone("h_up_met2");
        	TH1D *h_down_met2 = (TH1D*)Nominal_met_new->Clone("h_down_met2");
        	h_nominal_met2->Sumw2();
        	h_nominal_met2->SetStats(0);
        	h_nominal_met2->Divide(h_up_met2);
        	h_nominal_met2->SetMarkerStyle(6);
        	//h_nominal_met2->SetMarkerSize(20);
        	h_nominal_met2->SetMarkerColor(kRed);
        	h_nominal_met2->Draw("hist p e1");
        	h_nominal_met3->Sumw2();
        	h_nominal_met3->SetStats(0);
        	h_nominal_met3->Divide(h_down_met2);
        	h_nominal_met3->SetMarkerStyle(6);
        	h_nominal_met3->SetMarkerColor(kGreen);
        	h_nominal_met3->Draw("hist p e1 same");
        	//Nominal_met->Draw("ep same");

  				// Ratio plot (h3) settings
   				h_nominal_met2->SetTitle(""); // Remove the ratio title

   				// Y axis ratio plot settings
   				h_nominal_met2->GetYaxis()->SetTitle("Up(Down)/Nominal ");
   				h_nominal_met2->GetYaxis()->SetNdivisions(505);
   				h_nominal_met2->GetYaxis()->SetTitleSize(20);
   				h_nominal_met2->GetYaxis()->SetTitleFont(43);
   				h_nominal_met2->GetYaxis()->SetTitleOffset(1.55);
   				h_nominal_met2->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   				h_nominal_met2->GetYaxis()->SetLabelSize(15);

   				//X axis ratio plot settings
   				h_nominal_met2->GetXaxis()->SetTitle(XaxisName);
   				h_nominal_met2->GetXaxis()->SetTitleSize(20);
   				h_nominal_met2->GetXaxis()->SetTitleFont(43);
   				h_nominal_met2->GetXaxis()->SetTitleOffset(4.);
   				h_nominal_met2->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   				h_nominal_met2->GetXaxis()->SetLabelSize(15);

					//Draw a line at 1
					TLine *line = new TLine(h_nominal_met2->GetXaxis()->GetXmin(),1.,h_nominal_met2->GetXaxis()->GetXmax() ,1.);
					line->SetLineColor(kBlack);
					line->SetLineWidth(1);
					line->SetLineStyle(3);
					line->Draw();

				}
				if(VERBOSE)std::cout << __LINE__ << std::endl;
				f_output->cd(globalDirectoryName+"/"+samples[s]+"_"+chTags[ich]);
				for(unsigned int k=position_canvas_reset; k < c.size(); k++) c[k]->Write();
				position_canvas_reset = c.size();
				if(VERBOSE)std::cout << __LINE__ << std::endl;
			}

		}


	}





	//
	//
	//Let's add specifc corrections per sample:
	//
	//
	if(VERBOSE)std::cout << "IN THE NEW PART" << std::endl;
	int position_fromTheory =0;
	int cpos_fromTheory=0;
	int emptyPosition_fromTheory=0;
	int position_canvas_reset_fromTheory =0;
	int canvasOldSize = 0;

	for(unsigned int s =0; s < samples.size(); s++){
		for ( unsigned int ich = 0; ich < chTags.size() ; ich++){
			directory_pos = s*chTags.size() + ich;
			for ( unsigned ev = 0; ev < evCat.size() ; ev++){
				canvasOldSize = c_fromTheory.size();

				TString tags_full=chTags[ich]+evCat[ev]; 
				position_fromTheory = s*chTags.size()*evCat.size() + ich*evCat.size() + ev;
				if(VERBOSE)std::cout<< "POSITON: " << position_fromTheory << std::endl;
				//Canvas
				if(VERBOSE)std::cout<< "take histo" << std::endl;
				//Take the histos of interest
				if(VERBOSE)std::cout << "Sample: " << samples[s] << std::endl;
				f_input->cd(samples[s]);		//This argument gives the file on which we should take the shape. We have to loop on my samples
				//Nominal
				met_shapes_fromTheory.push_back( (TH2F*) gDirectory->Get(tags_full+"_"+SystType+"_shapes"));

				if( !met_shapes_fromTheory[position_fromTheory]){
					if(VERBOSE)std::cout<< "Empty histogram, skip it" << std:: endl;
					//emptyPosition++;
					continue;
				}



 				TH1D* Nominal_met_fromTheory = met_shapes_fromTheory[position_fromTheory]->ProjectionY("Nominal_met_fromTheory", projectionBin, projectionBin); 

				//1. PDF+scale uncertainty
				if(VERBOSE)std::cout<< "In PDF+scale unc" << std::endl;
				c_fromTheory.push_back( new TCanvas(tags_full+"_met_shapes_corrections_PDFandScale", tags_full+"_met_shapes_corrections_PDFandScale", 800, 800));
				TH1D *Nominal_met_PDFandScale_up = (TH1D*) Nominal_met_fromTheory->Clone();
				TH1D *Nominal_met_PDFandScale_down = (TH1D*) Nominal_met_fromTheory->Clone();
				double binContent =0;
				double PDFandScale_unc =0;
				if(samples[s].Contains("W#gamma")) PDFandScale_unc = 0.050; // W#gamma #rightarrow l#nu#gamma : from http://journals.aps.org/prd/pdf/10.1103/PhysRevD.89.092005
				if(samples[s].Contains("Z#gamma")) PDFandScale_unc = 0.50;//0.057; //Due to QCD NLO corrections this is set to 50% (see: http://link.springer.com/article/10.1007%2FJHEP02%282016%29057) // Z#gamma #rightarrow #nu#nu#gamma : from http://journals.aps.org/prd/pdf/10.1103/PhysRevD.89.092005
				if(samples[s].Contains("W#rightarrow")) PDFandScale_unc = 0.035; // W#rightarrow l#nu : from CMS-PAS-SMP-15-004 -- https://cds.cern.ch/record/2093537
				for(unsigned int bin=1 ; bin < Nominal_met_fromTheory->GetNbinsX(); bin++){
					binContent = Nominal_met_fromTheory->GetBinContent(bin);
					if(binContent){
						Nominal_met_PDFandScale_up->SetBinContent(bin, binContent+ binContent*PDFandScale_unc); 
						Nominal_met_PDFandScale_down->SetBinContent(bin, binContent- binContent*PDFandScale_unc); 
					}
				}

				//2. Electroweak corrections
				if(VERBOSE)std::cout<< "In Ewk corr" << std::endl;
				c_fromTheory.push_back( new TCanvas(tags_full+"_met_shapes_corrections_EwkCorrections", tags_full+"_met_shapes_corrections_EwkCorrections", 800, 800));
				TH1D *Nominal_met_EwkCorrections_up = (TH1D*) Nominal_met_fromTheory->Clone();
				TH1D *Nominal_met_EwkCorrections_down = (TH1D*) Nominal_met_fromTheory->Clone();
				binContent =0;
				double EwkCorrections_unc =0;
				if(samples[s].Contains("W#gamma")) EwkCorrections_unc = 0.070; // W#gamma #rightarrow l#nu#gamma; from arXiv:1412.7421
				if(samples[s].Contains("Z#gamma")) EwkCorrections_unc = 0.110; // Z#gamma #rightarrow #nu#nu#gamma; from arXiv:1510.08742
				if(samples[s].Contains("W#rightarrow")) EwkCorrections_unc = 0.000; // W#rightarrow l#nu 
				for(unsigned int bin=1 ; bin < Nominal_met_fromTheory->GetNbinsX(); bin++){
					binContent = Nominal_met_fromTheory->GetBinContent(bin);
					if(binContent){
						Nominal_met_EwkCorrections_up->SetBinContent(bin, binContent+ binContent*EwkCorrections_unc); 
						Nominal_met_EwkCorrections_down->SetBinContent(bin, binContent- binContent*EwkCorrections_unc); 
					}
				}

				//3. Statistical fluctuations on MC sample
				if(VERBOSE)std::cout<< "In stats" << std::endl;
				c_fromTheory.push_back( new TCanvas(tags_full+"_met_shapes_corrections_stats", tags_full+"_met_shapes_corrections_stats", 800, 800));
				TH1D *Nominal_met_stats_up = (TH1D*) Nominal_met_fromTheory->Clone();
				TH1D *Nominal_met_stats_down = (TH1D*) Nominal_met_fromTheory->Clone();
				binContent =0;
				for(unsigned int bin=1 ; bin < Nominal_met_fromTheory->GetNbinsX(); bin++){
          binContent = Nominal_met_fromTheory->GetBinContent(bin);
          if(binContent){ //Empty bins should also have a statistical error I guess...
            Nominal_met_stats_up->SetBinContent(bin, binContent+ Nominal_met_fromTheory->GetBinError(bin));
            Nominal_met_stats_down->SetBinContent(bin, binContent- Nominal_met_fromTheory->GetBinError(bin));
          }
       	}


				//Style part
				std::vector< TH1D*> Nominal_met_new, Up_met_new, Down_met_new;
				std::vector< TH1D*> Nominal_met_new2, Up_met_new2, Down_met_new2;

				//non-rebinned plots that have to be passed
				Nominal_met_new2.push_back((TH1D*) Nominal_met_fromTheory->Clone());
				Nominal_met_new2.push_back((TH1D*) Nominal_met_fromTheory->Clone());
				Nominal_met_new2.push_back((TH1D*) Nominal_met_fromTheory->Clone());

        Up_met_new2.push_back((TH1D*) Nominal_met_PDFandScale_up->Clone());
        Up_met_new2.push_back((TH1D*) Nominal_met_EwkCorrections_up->Clone());
        Up_met_new2.push_back((TH1D*) Nominal_met_stats_up->Clone());

        Down_met_new2.push_back((TH1D*) Nominal_met_PDFandScale_down->Clone());
        Down_met_new2.push_back((TH1D*) Nominal_met_EwkCorrections_down->Clone());
        Down_met_new2.push_back((TH1D*) Nominal_met_stats_down->Clone());

				//rebin just to plot things
				TH1D *Nominal_met_fromTheory_rebin = (TH1D*) Nominal_met_fromTheory->Rebin(nmetAxis-1, "Nominal_met_fromTheory_rebin", metaxis);
        Nominal_met_fromTheory_rebin->Scale(1.0, "width");
        Nominal_met_fromTheory_rebin->Sumw2();

				TH1D *Nominal_met_PDFandScale_up_rebin = (TH1D*) Nominal_met_PDFandScale_up->Rebin(nmetAxis-1, "Nominal_met_PDFandScale_up_rebin", metaxis);
        Nominal_met_PDFandScale_up_rebin->Scale(1.0, "width");
        Nominal_met_PDFandScale_up_rebin->Sumw2();

				TH1D *Nominal_met_EwkCorrections_up_rebin = (TH1D*) Nominal_met_EwkCorrections_up->Rebin(nmetAxis-1, "Nominal_met_EwkCorrections_up_rebin", metaxis);
        Nominal_met_EwkCorrections_up_rebin->Scale(1.0, "width");
        Nominal_met_EwkCorrections_up_rebin->Sumw2();

				TH1D *Nominal_met_stats_up_rebin = (TH1D*) Nominal_met_stats_up->Rebin(nmetAxis-1, "Nominal_met_stats_up_rebin", metaxis);
        Nominal_met_stats_up_rebin->Scale(1.0, "width");
        Nominal_met_stats_up_rebin->Sumw2();


				TH1D *Nominal_met_PDFandScale_down_rebin = (TH1D*) Nominal_met_PDFandScale_down->Rebin(nmetAxis-1, "Nominal_met_PDFandScale_down_rebin", metaxis);
        Nominal_met_PDFandScale_down_rebin->Scale(1.0, "width");
        Nominal_met_PDFandScale_down_rebin->Sumw2();

				TH1D *Nominal_met_EwkCorrections_down_rebin = (TH1D*) Nominal_met_EwkCorrections_down->Rebin(nmetAxis-1, "Nominal_met_EwkCorrections_down_rebin", metaxis);
        Nominal_met_EwkCorrections_down_rebin->Scale(1.0, "width");
        Nominal_met_EwkCorrections_down_rebin->Sumw2();

				TH1D *Nominal_met_stats_down_rebin = (TH1D*) Nominal_met_stats_down->Rebin(nmetAxis-1, "Nominal_met_stats_down_rebin", metaxis);
        Nominal_met_stats_down_rebin->Scale(1.0, "width");
        Nominal_met_stats_down_rebin->Sumw2();

				Nominal_met_new.push_back((TH1D*) Nominal_met_fromTheory_rebin->Clone());
				Nominal_met_new.push_back((TH1D*) Nominal_met_fromTheory_rebin->Clone());
				Nominal_met_new.push_back((TH1D*) Nominal_met_fromTheory_rebin->Clone());

        Up_met_new.push_back((TH1D*) Nominal_met_PDFandScale_up_rebin->Clone());
        Up_met_new.push_back((TH1D*) Nominal_met_EwkCorrections_up_rebin->Clone());
        Up_met_new.push_back((TH1D*) Nominal_met_stats_up_rebin->Clone());

        Down_met_new.push_back((TH1D*) Nominal_met_PDFandScale_down_rebin->Clone());
        Down_met_new.push_back((TH1D*) Nominal_met_EwkCorrections_down_rebin->Clone());
        Down_met_new.push_back((TH1D*) Nominal_met_stats_down_rebin->Clone());





        TLegend *leg_met = new TLegend(0.6,0.7,0.89,0.89);
				std::vector<TString > corrections;
				corrections.push_back("PDF+Scale");
				corrections.push_back("Electroweak corrections");
				corrections.push_back("Statistical fluctuations");

				for(unsigned int index=0; index < Nominal_met_new.size(); index++){
					if(VERBOSE)std::cout << "In my custom loop" << std::endl;

					if(VERBOSE)std::cout << "Number of events = " << Nominal_met_new2[index]->GetEntries() << std::endl;
					if(Nominal_met_new[index]->GetEntries() ==0) continue;
					savedHisto[s][ich][ev].push_back((TH1D*) Nominal_met_new2[index]->Clone());
					savedHisto_Up[s][ich][ev].push_back((TH1D*) Up_met_new2[index]->Clone());
					savedHisto_Down[s][ich][ev].push_back((TH1D*) Down_met_new2[index]->Clone());



					//Nominal in black
					Nominal_met_new[index]->SetLineColor(kBlack);

					//Up in red
        	Up_met_new[index]->SetLineColor(kRed);

					//Down in green
        	Down_met_new[index]->SetLineColor(kGreen);

					//No stats
 					Nominal_met_new[index]->SetStats(kFALSE);


					//Make the legends
					leg_met->Clear();
        	leg_met->AddEntry(Nominal_met_new[index], "Nominal", "l");
        	leg_met->AddEntry(Up_met_new[index], corrections[index]+" up", "l");
        	leg_met->AddEntry(Down_met_new[index], corrections[index]+" down", "l");






					if(VERBOSE)std::cout<< "ready to fill" << std::endl;
		  		//For met_shapes
		  		//cpos++;
		  		cpos_fromTheory=canvasOldSize+index;
		  		c_fromTheory[cpos_fromTheory]->cd();
        	//c[cpos]->SetBottomMargin(6);
        	TPad *pad_met1 = new TPad("pad_met1","pad_met1",0,0.3,1,1);
        	TPad *pad_met2 = new TPad("pad_met2","pad_met2",0,0,1,0.3);
        	pad_met1->SetTopMargin(0.05);
        	pad_met1->SetBottomMargin(0);
        	pad_met1->Draw();
        	pad_met2->SetTopMargin(0);
        	pad_met2->SetBottomMargin(0.25);
        	pad_met2->Draw();
        	pad_met1->cd();
        	TH1 *h_nominal_met = Nominal_met_new[index]->DrawCopy(); 
        	h_nominal_met->SetTitle("");
        	TH1 *h_up_met = Up_met_new[index]->DrawCopy("same"); 
        	TH1 *h_down_met = Down_met_new[index]->DrawCopy("same"); 

        	leg_met->Draw(); //same not needed?

        	// Y axis ratio plot settings
        	h_nominal_met->GetYaxis()->SetTitle("Events (a.u.)");
        	h_nominal_met->GetYaxis()->SetNdivisions(505);
        	h_nominal_met->GetYaxis()->SetTitleSize(20);
        	h_nominal_met->GetYaxis()->SetTitleFont(43);
        	h_nominal_met->GetYaxis()->SetTitleOffset(1.55);
        	h_nominal_met->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
        	h_nominal_met->GetYaxis()->SetLabelSize(15); 

        	pad_met2->cd();
        	TH1D *h_nominal_met2 = (TH1D*)Up_met_new[index]->Clone("h_nominal_met2"); //Be careful: non intuitive naming conventions to earn time.  
        	TH1D *h_nominal_met3 = (TH1D*)Down_met_new[index]->Clone("h_nominal_met3");
        	TH1D *h_up_met2 = (TH1D*)Nominal_met_new[index]->Clone("h_up_met2");
        	TH1D *h_down_met2 = (TH1D*)Nominal_met_new[index]->Clone("h_down_met2");
        	h_nominal_met2->Sumw2();
        	h_nominal_met2->SetStats(0);
        	h_nominal_met2->Divide(h_up_met2);
        	h_nominal_met2->SetMarkerStyle(6);
        	//h_nominal_met2->SetMarkerSize(20);
        	h_nominal_met2->SetMarkerColor(kRed);
        	h_nominal_met2->Draw("hist p e1");
        	h_nominal_met3->Sumw2();
        	h_nominal_met3->SetStats(0);
        	h_nominal_met3->Divide(h_down_met2);
        	h_nominal_met3->SetMarkerStyle(6);
        	h_nominal_met3->SetMarkerColor(kGreen);
        	h_nominal_met3->Draw("hist p e1 same");
        	//Nominal_met->Draw("ep same");

  				// Ratio plot (h3) settings
   				h_nominal_met2->SetTitle(""); // Remove the ratio title

   				// Y axis ratio plot settings
   				h_nominal_met2->GetYaxis()->SetTitle("Up(Down)/Nominal ");
   				h_nominal_met2->GetYaxis()->SetNdivisions(505);
   				h_nominal_met2->GetYaxis()->SetTitleSize(20);
   				h_nominal_met2->GetYaxis()->SetTitleFont(43);
   				h_nominal_met2->GetYaxis()->SetTitleOffset(1.55);
   				h_nominal_met2->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   				h_nominal_met2->GetYaxis()->SetLabelSize(15);

   				//X axis ratio plot settings
   				h_nominal_met2->GetXaxis()->SetTitle(XaxisName);
   				h_nominal_met2->GetXaxis()->SetTitleSize(20);
   				h_nominal_met2->GetXaxis()->SetTitleFont(43);
   				h_nominal_met2->GetXaxis()->SetTitleOffset(4.);
   				h_nominal_met2->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   				h_nominal_met2->GetXaxis()->SetLabelSize(15);

					//Draw a line at 1
					//c1.Range( 0., -10., 1., 10. );
					TLine *line = new TLine(h_nominal_met2->GetXaxis()->GetXmin(),1.,h_nominal_met2->GetXaxis()->GetXmax() ,1.);
					line->SetLineColor(kBlack);
					line->SetLineWidth(1);
					line->SetLineStyle(3);
					line->Draw();



					if(VERBOSE)std::cout << "In my custom loop: end!" << std::endl;




					f_output->cd(globalDirectoryName+"/"+samples[s]+"_"+chTags[ich]);
					c_fromTheory[cpos_fromTheory]->Write();

				}
			}
		}
	}



	//
	//Now, let's make the final shapes up and down for the mains processes:
	//
	if(VERBOSE)std::cout<< "In the wrapping up loop" << std::endl;
	double square_sum_up=0;
	double square_sum_down=0;
	double syst_binContent_nominal = 0;
	double syst_binContent_up = 0;
	double syst_binContent_down = 0;
	int num_canvas =0;
	double updatedBinContent_up = 0;
	double updatedBinContent_down = 0;

	for(unsigned int s =0; s < samples.size(); s++){
    for ( unsigned int ich = 0; ich < chTags.size() ; ich++){
      directory_pos = s*chTags.size() + ich;
      for ( unsigned ev = 0; ev < evCat.size() ; ev++){
				TString tags_full=chTags[ich]+evCat[ev];
				if(savedHisto[s][ich][ev].size() ==0) continue;	
				TH1D *h_nominal_sum = (TH1D*)savedHisto[s][ich][ev][0]->Clone("h_nominal_sum");
				TH1D *h_nominal_sum_up = (TH1D*)savedHisto[s][ich][ev][0]->Clone("h_nominal_sum_up");
				TH1D *h_nominal_sum_down = (TH1D*)savedHisto[s][ich][ev][0]->Clone("h_nominal_sum_down");
				for(unsigned int bin=1; bin < savedHisto[s][ich][ev][0]->GetNbinsX(); bin++){
					for(unsigned int v =0; v < savedHisto[s][ich][ev].size() ; v++){

						syst_binContent_nominal = savedHisto[s][ich][ev][v]->GetBinContent(bin); 
						syst_binContent_up = savedHisto_Up[s][ich][ev][v]->GetBinContent(bin); 
						syst_binContent_down = savedHisto_Down[s][ich][ev][v]->GetBinContent(bin); 

						if(VERBOSE)std::cout << "v = " << v << " et syst_binContent_nominal = " << syst_binContent_nominal << " et syst_binContent_up = " << syst_binContent_up << " et syst_binContent_down = " << syst_binContent_down << " et square_sum_up = " << square_sum_up << " et square_sum_down = " << square_sum_down << std::endl;
						square_sum_up += 1.0*(syst_binContent_up - syst_binContent_nominal)*(syst_binContent_up - syst_binContent_nominal);
						square_sum_down += 1.0*(syst_binContent_down - syst_binContent_nominal)*(syst_binContent_down - syst_binContent_nominal);

					}
					if(VERBOSE)std::cout<< "Bin = " << bin << " et sqrt(square_sum_up) = " << sqrt(square_sum_up) << std::endl;
					if(h_nominal_sum_up->GetBinContent(bin) !=0 && h_nominal_sum_down->GetBinContent(bin) !=0){
						updatedBinContent_up = h_nominal_sum_up->GetBinContent(bin) + sqrt(square_sum_up);
						updatedBinContent_down = h_nominal_sum_down->GetBinContent(bin) - sqrt(square_sum_down);
						if(VERBOSE)std::cout<< "UpdatedBinContent_up = " << updatedBinContent_up << " and updatedBinContent_down = " << updatedBinContent_down << std::endl;	
						h_nominal_sum_up->SetBinContent(bin, updatedBinContent_up);
						h_nominal_sum_down->SetBinContent(bin, updatedBinContent_down);
					}
					square_sum_up=0;
					square_sum_down=0;

				}


				//Style part

				TH1D *Nominal_met_new2 = (TH1D*) h_nominal_sum->Clone();
        TH1D *Up_met_new2 = (TH1D*) h_nominal_sum_up->Clone();
        TH1D *Down_met_new2 = (TH1D*) h_nominal_sum_down->Clone();

				//Rebin to plot
				TH1D *Nominal_met_new = (TH1D*) h_nominal_sum->Rebin(nmetAxis-1, "Nominal_met_new", metaxis);
        Nominal_met_new->Scale(1.0, "width");
        Nominal_met_new->Sumw2();
				TH1D *Up_met_new = (TH1D*) h_nominal_sum_up->Rebin(nmetAxis-1, "Up_met_new", metaxis);
        Up_met_new->Scale(1.0, "width");
        Up_met_new->Sumw2();
				TH1D *Down_met_new = (TH1D*) h_nominal_sum_down->Rebin(nmetAxis-1, "Down_met_new", metaxis);
        Down_met_new->Scale(1.0, "width");
        Down_met_new->Sumw2();



				//Save the final wrap up shapes up and down for each category
				if(VERBOSE)std::cout << "Number of events = " << Nominal_met_new->GetEntries() << std::endl;
        if(Nominal_met_new->GetEntries() ==0) continue;
        wrapUp_shape[s][ich][ev].push_back((TH1D*) Nominal_met_new2->Clone());
        wrapUp_shape_Up[s][ich][ev].push_back((TH1D*) Up_met_new2->Clone());
        wrapUp_shape_Down[s][ich][ev].push_back((TH1D*) Down_met_new2->Clone());



				//Nominal in black
				Nominal_met_new->SetLineColor(kBlack);

				//Up in red
        Up_met_new->SetLineColor(kRed);

				//Down in green
        Down_met_new->SetLineColor(kGreen);

				//No stats
 				Nominal_met_new->SetStats(kFALSE);


				//Make the legends

        TLegend *leg_met = new TLegend(0.6,0.7,0.89,0.89);
        leg_met->AddEntry(Nominal_met_new, "Nominal", "l");
        leg_met->AddEntry(Up_met_new, "full corrections up", "l");
        leg_met->AddEntry(Down_met_new, "full corrections down", "l");






				if(VERBOSE)std::cout<< "ready to fill" << std::endl;

				c_wrapUp.push_back( new TCanvas(tags_full+"_met_shapes_FINAL_WRAP_UP", tags_full+"_met_shapes_corrections_FINAL_WRAP_UP", 800, 800));
        c_wrapUp[num_canvas]->cd();
				//c[cpos]->SetBottomMargin(6);
        TPad *pad_met1 = new TPad("pad_met1","pad_met1",0,0.3,1,1);
        TPad *pad_met2 = new TPad("pad_met2","pad_met2",0,0,1,0.3);
        pad_met1->SetTopMargin(0.05);
        pad_met1->SetBottomMargin(0);
        pad_met1->Draw();
        pad_met2->SetTopMargin(0);
        pad_met2->SetBottomMargin(0.25);
        pad_met2->Draw();
        pad_met1->cd();
        TH1 *h_nominal_met = Nominal_met_new->DrawCopy(); 
        h_nominal_met->SetTitle("");
        TH1 *h_up_met = Up_met_new->DrawCopy("same"); 
        TH1 *h_down_met = Down_met_new->DrawCopy("same"); 

        leg_met->Draw(); //same not needed?

        // Y axis ratio plot settings
        h_nominal_met->GetYaxis()->SetTitle("Events (a.u.)");
        h_nominal_met->GetYaxis()->SetNdivisions(505);
        h_nominal_met->GetYaxis()->SetTitleSize(20);
        h_nominal_met->GetYaxis()->SetTitleFont(43);
        h_nominal_met->GetYaxis()->SetTitleOffset(1.55);
        h_nominal_met->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
        h_nominal_met->GetYaxis()->SetLabelSize(15); 

        pad_met2->cd();
        TH1D *h_nominal_met2 = (TH1D*)Up_met_new->Clone("h_nominal_met2"); //Be careful: non intuitive naming conventions to earn time. 
        TH1D *h_nominal_met3 = (TH1D*)Down_met_new->Clone("h_nominal_met3");
        TH1D *h_up_met2 = (TH1D*)Nominal_met_new->Clone("h_up_met2");
        TH1D *h_down_met2 = (TH1D*)Nominal_met_new->Clone("h_down_met2");
        h_nominal_met2->Sumw2();
        h_nominal_met2->SetStats(0);
        h_nominal_met2->Divide(h_up_met2);
        h_nominal_met2->SetMarkerStyle(6);
        //h_nominal_met2->SetMarkerSize(20);
        h_nominal_met2->SetMarkerColor(kRed);
        h_nominal_met2->Draw("hist p e1");
        h_nominal_met3->Sumw2();
        h_nominal_met3->SetStats(0);
        h_nominal_met3->Divide(h_down_met2);
        h_nominal_met3->SetMarkerStyle(6);
        h_nominal_met3->SetMarkerColor(kGreen);
        h_nominal_met3->Draw("hist p e1 same");
        //Nominal_met->Draw("ep same");

  			// Ratio plot (h3) settings
   			h_nominal_met2->SetTitle(""); // Remove the ratio title

   			// Y axis ratio plot settings
   			h_nominal_met2->GetYaxis()->SetTitle("Up(Down)/Nominal ");
   			h_nominal_met2->GetYaxis()->SetNdivisions(505);
   			h_nominal_met2->GetYaxis()->SetTitleSize(20);
   			h_nominal_met2->GetYaxis()->SetTitleFont(43);
   			h_nominal_met2->GetYaxis()->SetTitleOffset(1.55);
   			h_nominal_met2->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   			h_nominal_met2->GetYaxis()->SetLabelSize(15);

   			//X axis ratio plot settings
   			h_nominal_met2->GetXaxis()->SetTitle(XaxisName);
   			h_nominal_met2->GetXaxis()->SetTitleSize(20);
   			h_nominal_met2->GetXaxis()->SetTitleFont(43);
   			h_nominal_met2->GetXaxis()->SetTitleOffset(4.);
   			h_nominal_met2->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   			h_nominal_met2->GetXaxis()->SetLabelSize(15);

				//Draw a line at 1
				//c1.Range( 0., -10., 1., 10. );
				TLine *line = new TLine(h_nominal_met2->GetXaxis()->GetXmin(),1.,h_nominal_met2->GetXaxis()->GetXmax() ,1.);
				line->SetLineColor(kBlack);
				line->SetLineWidth(1);
				line->SetLineStyle(3);
				line->Draw();


				if(VERBOSE)std::cout << __LINE__ << std::endl;




				f_output->cd(globalDirectoryName+"/"+samples[s]+"_"+chTags[ich]);
				c_wrapUp[num_canvas]->Write();
				num_canvas++;



			}
		}
	}









	//
	//Now let's construct the instrumental MET shape up and down for the genuine MET
	//
	if(VERBOSE)std::cout<< "Instrumental MET construction" << std::endl;

	std::vector<TString> genuineMET;
	genuineMET.push_back("W#gamma #rightarrow l#nu#gamma_reweighted"); 
	genuineMET.push_back("Z#gamma #rightarrow ll#gamma_reweighted");
	genuineMET.push_back("Z#rightarrow #nu#nu_reweighted");
	genuineMET.push_back("W#rightarrow l#nu_reweighted"); 
	genuineMET.push_back("Top_reweighted");
	genuineMET.push_back("WZ_reweighted");
	genuineMET.push_back("WW_reweighted");
	genuineMET.push_back("Z#rightarrow #tau#tau_filt15_reweighted");
	genuineMET.push_back("Z#gamma #rightarrow #nu#nu#gamma_reweighted"); 
	genuineMET.push_back("Top+#gamma_reweighted");
	genuineMET.push_back("Z#rightarrow ee-#mu#mu_filt1113_reweighted");


  std::vector<std::vector<std::vector<TH1D*> > > v_genuineMET(chTags.size(), std::vector<std::vector<TH1D*> >(evCat.size()));
  std::vector<std::vector<std::vector<TH1D*> > > v_gammaData(chTags.size(), std::vector<std::vector<TH1D*> >(evCat.size()));
  std::vector<std::vector<std::vector<TH1D*> > > v_genuineMET_nominal(chTags.size(), std::vector<std::vector<TH1D*> >(evCat.size()));
  std::vector<std::vector<std::vector<TH1D*> > > v_genuineMET_nominal2(chTags.size(), std::vector<std::vector<TH1D*> >(evCat.size()));
  std::vector<std::vector<std::vector<TH1D*> > > v_gammaData_nominal(chTags.size(), std::vector<std::vector<TH1D*> >(evCat.size()));

	std::vector<TH2F*> met_shape_tmp, met_shape_tmp_genuineMET, met_shape_tmp_genuineMET2;
  std::vector<std::vector<std::vector<TH1D*> > > Instr_MET_nominal(chTags.size(), std::vector<std::vector<TH1D*> >(evCat.size()));
	std::vector<std::vector<std::vector<std::vector<TH1D*> > > > Instr_MET_up(samples.size(), std::vector<std::vector<std::vector<TH1D*> > >(chTags.size(), std::vector<std::vector<TH1D*> >(evCat.size())));
	std::vector<std::vector<std::vector<std::vector<TH1D*> > > > Instr_MET_down(samples.size(), std::vector<std::vector<std::vector<TH1D*> > >(chTags.size(), std::vector<std::vector<TH1D*> >(evCat.size())));
  std::vector<std::vector<std::vector<TH1D*> > > Instr_MET_nominal_rebin(chTags.size(), std::vector<std::vector<TH1D*> >(evCat.size()));
	std::vector<std::vector<std::vector<std::vector<TH1D*> > > > Instr_MET_up_rebin(samples.size(), std::vector<std::vector<std::vector<TH1D*> > >(chTags.size(), std::vector<std::vector<TH1D*> >(evCat.size())));
	std::vector<std::vector<std::vector<std::vector<TH1D*> > > > Instr_MET_down_rebin(samples.size(), std::vector<std::vector<std::vector<TH1D*> > >(chTags.size(), std::vector<std::vector<TH1D*> >(evCat.size())));

	//Let's take all the histo of interest for the genuineMET

	int counter_instrMET=-1;
	int counter_genuineMET =-1;
	int current_position=0;
	int current_position_genuineMET=-1;
	globalDirectory->cd();
	TDirectory* dir_InstrMET = 	f_output->mkdir(globalDirectoryName+"/"+"Instr_MET_nominal");

  std::vector< TDirectory*> directory_Instr_MET;
  for(unsigned int s =0; s < samples.size(); s++){
    for ( unsigned int ich = 0; ich < chTags.size() ; ich++){
			globalDirectory->cd();
      directory_Instr_MET.push_back(f_output->mkdir(globalDirectoryName+"/"+"Instr_MET_"+samples[s]+"_"+chTags[ich]));
    }
  } 	

 	for ( unsigned int ich = 0; ich < chTags.size() ; ich++){
    for ( unsigned ev = 0; ev < evCat.size() ; ev++){
      TString tags_full=chTags[ich]+evCat[ev];
			if(VERBOSE)std::cout<< "Tag: " << tags_full << std::endl;

			current_position = ich*evCat.size() +ev;
			if(VERBOSE)std::cout<< "current_position: " << current_position << std::endl;
			//and for the gamma data reweighted
			f_input->cd("#gamma data_reweighted");
  		met_shape_tmp.push_back((TH2F*) gDirectory->Get(tags_full+"_"+SystType+"_shapes"));	
			if(VERBOSE)std::cout << "th2f taken" << std::endl;
      if(met_shape_tmp[current_position]){
				//Initialize Instr.MET with gamma data
				Instr_MET_nominal[ich][ev].push_back(met_shape_tmp[current_position]->ProjectionY(tags_full+"_InstrMET_nominal", projectionBin, projectionBin));
				v_gammaData_nominal[ich][ev].push_back((TH1D*) Instr_MET_nominal[ich][ev][0]->Clone());

				for(unsigned int s =0; s < samples.size(); s++){
					Instr_MET_up[s][ich][ev].push_back((TH1D*) Instr_MET_nominal[ich][ev][0]->Clone());
					Instr_MET_down[s][ich][ev].push_back((TH1D*) Instr_MET_nominal[ich][ev][0]->Clone());
				}

			}
			else{
				continue;
			}

			if(Instr_MET_nominal[ich][ev][0]->GetEntries() ==0) continue;
			if(VERBOSE)std::cout<< "Gamma data found" << std::endl;
			counter_genuineMET=-1;
			for(unsigned int index=0; index < genuineMET.size(); index++){
				f_input->cd(genuineMET[index]); 
  			met_shape_tmp_genuineMET.push_back((TH2F*) gDirectory->Get(tags_full+"_"+SystType+"_shapes"));	
  			current_position_genuineMET++;
				if(VERBOSE)std::cout <<"current_position_genuineMET" << current_position_genuineMET << std::endl;
				if(met_shape_tmp_genuineMET[current_position_genuineMET]){
					if(VERBOSE)std::cout<< "genuineMET found" << std::endl;
      		v_genuineMET_nominal[ich][ev].push_back( met_shape_tmp_genuineMET[current_position_genuineMET]->ProjectionY(tags_full+genuineMET[index], projectionBin, projectionBin));

					if(VERBOSE)std::cout<< __LINE__ << std::endl;
					counter_genuineMET++;

					TH1D *nominal_tmp = (TH1D*) v_genuineMET_nominal[ich][ev][counter_genuineMET]->Clone();
					if(nominal_tmp->GetEntries() == 0) continue;
					Instr_MET_nominal[ich][ev][0]->Add(nominal_tmp,-1);
					Instr_MET_nominal[ich][ev][0]->Sumw2();

					if(VERBOSE)std::cout<< __LINE__ << std::endl;

					for(unsigned int s =0; s < samples.size(); s++){
						if(genuineMET[index] == samples[s]) continue;

						TH1D *up_tmp = (TH1D*) v_genuineMET_nominal[ich][ev][counter_genuineMET]->Clone();
				  	Instr_MET_up[s][ich][ev][0]->Add(up_tmp,-1);
						Instr_MET_up[s][ich][ev][0]->Sumw2();

						TH1D *down_tmp = (TH1D*) v_genuineMET_nominal[ich][ev][counter_genuineMET]->Clone();
				  	Instr_MET_down[s][ich][ev][0]->Add(down_tmp,-1);
						Instr_MET_down[s][ich][ev][0]->Sumw2();
					}
				}
			}

			if(VERBOSE)std::cout<< "Writing" << std::endl;
			//write histo to debug
			f_output->cd(globalDirectoryName+"/"+"Instr_MET_nominal");
			//Rebin to plot

			Instr_MET_nominal_rebin[ich][ev].push_back( (TH1D*) Instr_MET_nominal[ich][ev][0]->Rebin(nmetAxis-1, "Instr_MET_nominal_rebin_"+chTags[ich]+evCat[ev], metaxis));
			Instr_MET_nominal_rebin[ich][ev][0]->Scale(1.0, "width");
					if(VERBOSE)std::cout<< __LINE__ << std::endl;

			for(unsigned int s =0; s < samples.size(); s++){
				Instr_MET_up_rebin[s][ich][ev].push_back( (TH1D*) Instr_MET_up[s][ich][ev][0]->Rebin(nmetAxis-1, "Instr_MET_up_rebin", metaxis));
				Instr_MET_up_rebin[s][ich][ev][0]->Scale(1.0, "width");
				Instr_MET_down_rebin[s][ich][ev].push_back( (TH1D*) Instr_MET_down[s][ich][ev][0]->Rebin(nmetAxis-1, "Instr_MET_down_rebin", metaxis));
				Instr_MET_down_rebin[s][ich][ev][0]->Scale(1.0, "width");
			}
			Instr_MET_nominal_rebin[ich][ev][0]->Write();
					if(VERBOSE)std::cout<< __LINE__ << std::endl;

			//Instr MET from the plooter
			if(InstrMET_control_plots){
			  f_input->cd("Instr. MET");
			  TH2F* met_shape_true_InstrMET = (TH2F*) gDirectory->Get(tags_full+"_"+SystType+"_shapes");
			  if(VERBOSE)std::cout<< "Projecting " << tags_full <<"_"<<SystType<<"_shapes" << std::endl;
			  TH1D* h_met_shape_true_InstrMET = (TH1D*) met_shape_true_InstrMET->ProjectionY(tags_full+"_InstrMET_true", projectionBin, projectionBin);
			  		if(VERBOSE)std::cout<< __LINE__ << std::endl;
			  f_output->cd(globalDirectoryName+"/"+"Instr_MET_nominal");
			  h_met_shape_true_InstrMET = (TH1D*) h_met_shape_true_InstrMET->Rebin(nmetAxis-1, "", metaxis);
			  h_met_shape_true_InstrMET->Scale(1.0, "width");
			  h_met_shape_true_InstrMET->Write();
			  if(VERBOSE)std::cout<< "Writing done" << std::endl;
			}
		}
	}



	//We loop on the three main sources of genuine MET (the samples std::vector)
	//At this state, the instrMET up and down is without the sample of interest up/down so we add it
	num_canvas=0;
	for(unsigned int s =0; s < samples.size(); s++){


 		for ( unsigned int ich = 0; ich < chTags.size() ; ich++){
 			directory_pos = s*chTags.size() + ich;

 			for ( unsigned ev = 0; ev < evCat.size() ; ev++){
      	TString tags_full=chTags[ich]+evCat[ev];
				if(VERBOSE)std::cout<< "Tag: " << tags_full << std::endl;

				if(wrapUp_shape_Up[s][ich][ev].size() ==0 && wrapUp_shape_Down[s][ich][ev].size()== 0 ) continue;
				TH1D* h_sample_up = (TH1D*) wrapUp_shape_Up[s][ich][ev][0]->Clone();
				TH1D* h_sample_down = (TH1D*) wrapUp_shape_Down[s][ich][ev][0]->Clone();

				if(VERBOSE)std::cout << __LINE__ << std::endl;

				//This removes the samples s from the corresponding Instr.MET shape
       	Instr_MET_up[s][ich][ev][0]->Add(h_sample_down,-1);
       	Instr_MET_down[s][ich][ev][0]->Add(h_sample_up,-1);

				if(VERBOSE)std::cout << __LINE__ << std::endl;


				//Style part

				TH1D *Nominal_met_new2 = (TH1D*) Instr_MET_nominal[ich][ev][0]->Clone();
        TH1D *Up_met_new2 = (TH1D*) Instr_MET_up[s][ich][ev][0]->Clone();
        TH1D *Down_met_new2 = (TH1D*) Instr_MET_down[s][ich][ev][0]->Clone();

				//Rebin to plot
				TH1D *Nominal_met_new = (TH1D*) Instr_MET_nominal[ich][ev][0]->Rebin(nmetAxis-1, "Nominal_met_new", metaxis);
        Nominal_met_new->Scale(1.0, "width");
        Nominal_met_new->Sumw2();
				TH1D *Up_met_new = (TH1D*) Instr_MET_up[s][ich][ev][0]->Rebin(nmetAxis-1, "Up_met_new", metaxis);
        Up_met_new->Scale(1.0, "width");
        Up_met_new->Sumw2();
				TH1D *Down_met_new = (TH1D*) Instr_MET_down[s][ich][ev][0]->Rebin(nmetAxis-1, "Down_met_new", metaxis);
        Down_met_new->Scale(1.0, "width");
        Down_met_new->Sumw2();




				//Save the final wrap up shapes up and down for each category
				if(VERBOSE)std::cout << "Number of events = " << Nominal_met_new->GetEntries() << std::endl;
        if(Nominal_met_new->GetEntries() ==0) continue;
        savedInstrMET_shape[s][ich][ev].push_back((TH1D*) Nominal_met_new2->Clone());
        savedInstrMET_shape_Up[s][ich][ev].push_back((TH1D*) Up_met_new2->Clone());
        savedInstrMET_shape_Down[s][ich][ev].push_back((TH1D*) Down_met_new2->Clone());



				//Nominal in black
				Nominal_met_new->SetLineColor(kBlack);

				//Up in red
        Up_met_new->SetLineColor(kRed);

				//Down in green
        Down_met_new->SetLineColor(kGreen);

				//No stats
 				Nominal_met_new->SetStats(kFALSE);


				//Make the legends

        TLegend *leg_met = new TLegend(0.6,0.7,0.89,0.89);
        leg_met->AddEntry(Nominal_met_new, "Nominal", "l");
        leg_met->AddEntry(Up_met_new, "full corrections up", "l");
        leg_met->AddEntry(Down_met_new, "full corrections down", "l");






				if(VERBOSE)std::cout<< "ready to fill" << std::endl;

				c_InstrMET_genuineMET.push_back( new TCanvas(tags_full+"_met_shapes_FINAL_genuineMET_contribution", tags_full+"_met_shapes_corrections_FINAL_genuineMET_contribution", 800, 800));
        c_InstrMET_genuineMET[num_canvas]->cd();
				//c[cpos]->SetBottomMargin(6);
        TPad *pad_met1 = new TPad("pad_met1","pad_met1",0,0.3,1,1);
        TPad *pad_met2 = new TPad("pad_met2","pad_met2",0,0,1,0.3);
        pad_met1->SetTopMargin(0.05);
        pad_met1->SetBottomMargin(0);
        pad_met1->Draw();
        pad_met2->SetTopMargin(0);
        pad_met2->SetBottomMargin(0.25);
        pad_met2->Draw();
        pad_met1->cd();
        TH1 *h_nominal_met = Nominal_met_new->DrawCopy(); 
        h_nominal_met->SetTitle("");
        TH1 *h_up_met = Up_met_new->DrawCopy("same"); 
        TH1 *h_down_met = Down_met_new->DrawCopy("same"); 

        leg_met->Draw(); //same not needed?

        // Do not draw the Y axis label on the upper plot and redraw a small
		    // axis instead, in order to avoid the first label (0) to be clipped.
   			//h_nominal_met->GetYaxis()->SetLabelSize(0.);
   			//TGaxis *axis = new TGaxis( -5, -5, 5, -5, 0.01,200,510,"");
   			//axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   			//axis->SetLabelSize(15);
   			//axis->Draw();

        // Y axis ratio plot settings
        h_nominal_met->GetYaxis()->SetTitle("Events (a.u.)");
        h_nominal_met->GetYaxis()->SetNdivisions(505);
        h_nominal_met->GetYaxis()->SetTitleSize(20);
        h_nominal_met->GetYaxis()->SetTitleFont(43);
        h_nominal_met->GetYaxis()->SetTitleOffset(1.55);
        h_nominal_met->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
        h_nominal_met->GetYaxis()->SetLabelSize(15); 

        pad_met2->cd();
        TH1D *h_nominal_met2 = (TH1D*)Up_met_new->Clone("h_nominal_met2"); //Be careful: non intuitive naming conventions to earn time.  
        TH1D *h_nominal_met3 = (TH1D*)Down_met_new->Clone("h_nominal_met3");
        TH1D *h_up_met2 = (TH1D*)Nominal_met_new->Clone("h_up_met2");
        TH1D *h_down_met2 = (TH1D*)Nominal_met_new->Clone("h_down_met2");
        h_nominal_met2->Sumw2();
        h_nominal_met2->SetStats(0);
        h_nominal_met2->Divide(h_up_met2);
        h_nominal_met2->SetMarkerStyle(6);
        //h_nominal_met2->SetMarkerSize(20);
        h_nominal_met2->SetMarkerColor(kRed);
        h_nominal_met2->Draw("hist p e1");
        h_nominal_met3->Sumw2();
        h_nominal_met3->SetStats(0);
        h_nominal_met3->Divide(h_down_met2);
        h_nominal_met3->SetMarkerStyle(6);
        h_nominal_met3->SetMarkerColor(kGreen);
        h_nominal_met3->Draw("hist p e1 same");
        //Nominal_met->Draw("ep same");

  			// Ratio plot (h3) settings
   			h_nominal_met2->SetTitle(""); // Remove the ratio title

   			// Y axis ratio plot settings
   			h_nominal_met2->GetYaxis()->SetTitle("Up(Down)/Nominal ");
   			h_nominal_met2->GetYaxis()->SetNdivisions(505);
   			h_nominal_met2->GetYaxis()->SetTitleSize(20);
   			h_nominal_met2->GetYaxis()->SetTitleFont(43);
   			h_nominal_met2->GetYaxis()->SetTitleOffset(1.55);
   			h_nominal_met2->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   			h_nominal_met2->GetYaxis()->SetLabelSize(15);

   			//X axis ratio plot settings
   			h_nominal_met2->GetXaxis()->SetTitle(XaxisName);
   			h_nominal_met2->GetXaxis()->SetTitleSize(20);
   			h_nominal_met2->GetXaxis()->SetTitleFont(43);
   			h_nominal_met2->GetXaxis()->SetTitleOffset(4.);
   			h_nominal_met2->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   			h_nominal_met2->GetXaxis()->SetLabelSize(15);

				//Draw a line at 1
				//c1.Range( 0., -10., 1., 10. );
				TLine *line = new TLine(h_nominal_met2->GetXaxis()->GetXmin(),1.,h_nominal_met2->GetXaxis()->GetXmax() ,1.);
				line->SetLineColor(kBlack);
				line->SetLineWidth(1);
				line->SetLineStyle(3);
				line->Draw();




				//Now make ratio:
   			//h1->GetXaxis()->SetLabelFont(63); //font in pixels
   			//h1->GetXaxis()->SetLabelSize(16); //in pixels
   			//h1->GetYaxis()->SetLabelFont(63); //font in pixels
   			//h1->GetYaxis()->SetLabelSize(16); //in pixels
				if(VERBOSE)std::cout << __LINE__ << std::endl;




				f_output->cd(globalDirectoryName+"/"+"Instr_MET_"+samples[s]+"_"+chTags[ich]);
				c_InstrMET_genuineMET[num_canvas]->Write();
				num_canvas++;






			}
		}


	}





	std::vector< TCanvas*> c_InstrMET_cst;
	double binContent=0;
	current_position_genuineMET=-1;
	num_canvas=0;
	//
	//Now we add a constant uncertainty for the method (for each category + the stat uncertainty)
	//

  std::vector< TDirectory*> directory_Instr_MET_cst;
  for ( unsigned int ich = 0; ich < chTags.size() ; ich++){
		globalDirectory->cd();
    directory_Instr_MET_cst.push_back(f_output->mkdir(globalDirectoryName+"/"+"Instr_MET_cst_"+chTags[ich]));
  }

  for ( unsigned int ich = 0; ich < chTags.size() ; ich++){

    for ( unsigned ev = 0; ev < evCat.size() ; ev++){
      TString tags_full=chTags[ich]+evCat[ev];
      if(VERBOSE)std::cout<< "Tag: " << tags_full << std::endl;
			canvasOldSize = c_InstrMET_cst.size();

			//1. Compute error from method (here 10% everywhere on the Instr MET
      if(VERBOSE)std::cout<< "In InstrMET method corr" << std::endl;
      c_InstrMET_cst.push_back( new TCanvas(tags_full+"_met_shapes_InstrMET_method", tags_full+"_met_shapes_InstrMET_method", 800, 800));

      TH1D *Nominal_met_method_up = (TH1D*) Instr_MET_nominal[ich][ev][0]->Clone();
      TH1D *Nominal_met_method_down = (TH1D*) Instr_MET_nominal[ich][ev][0]->Clone();
      binContent =0;
      double method_unc =0.10; //For the moment we take a conservative value of 10% (extrapolated from the closure test). As already discussed, we take the same value for all lepton and jet categories
      for(unsigned int bin=1 ; bin < Instr_MET_nominal[ich][ev][0]->GetNbinsX(); bin++){
        binContent = Instr_MET_nominal[ich][ev][0]->GetBinContent(bin);

				/*if(SystType == "met"){
          if( Instr_MET_nominal[ich][ev][0]->GetBinLowEdge(bin) > 75) method_unc = 0.20 ;
          else if( Instr_MET_nominal[ich][ev][0]->GetBinLowEdge(bin) > 55) method_unc = 0.15 ;
          else if( Instr_MET_nominal[ich][ev][0]->GetBinLowEdge(bin) > 40) method_unc = 0.10 ;
          else method_unc = 0.05 ;
        }
        else if(SystType == "mt"){
          *//*if( Instr_MET_nominal[ich][ev][0]->GetBinLowEdge(bin) > 350) method_unc = 0.20 ;
          else*//* if( Instr_MET_nominal[ich][ev][0]->GetBinLowEdge(bin) > 300) method_unc = 0.05 ;
          else if( Instr_MET_nominal[ich][ev][0]->GetBinLowEdge(bin) > 275) method_unc = 0.15 ;
          else if( Instr_MET_nominal[ich][ev][0]->GetBinLowEdge(bin) > 250) method_unc = 0.10 ;
          else if( Instr_MET_nominal[ich][ev][0]->GetBinLowEdge(bin) > 175) method_unc = 0.05 ;
          else method_unc = 0.10 ;
        }*/
        method_unc = 0.10;
        if(binContent){
          if(SystType == "met"){
            Nominal_met_method_up->SetBinContent(bin, binContent); //Asymetric limit: goes only down
            Nominal_met_method_down->SetBinContent(bin, binContent - binContent*method_unc);
          }
          if(SystType == "mt"){
            Nominal_met_method_up->SetBinContent(bin, binContent + binContent*method_unc); //Asymetric limit: goes only up
            Nominal_met_method_down->SetBinContent(bin, binContent);
          }

        
        }
      }

      //2.a Statistical fluctuations on gamma+jet sample
      if(VERBOSE)std::cout<< "In gamma+jet stats" << std::endl;
      c_InstrMET_cst.push_back( new TCanvas(tags_full+"_met_shapes_InstrMET_stats", tags_full+"_met_shapes_InstrMET_stats", 800, 800));
      c_InstrMET_cst.push_back( new TCanvas(tags_full+"_met_shapes_gammaData_stats", tags_full+"_met_shapes_gammaData_stats", 800, 800));
      TH1D *Nominal_met_gammaOnly_stats = (TH1D*) v_gammaData_nominal[ich][ev][0]->Clone();
      TH1D *Nominal_met_gamma_stats_up = (TH1D*) v_gammaData_nominal[ich][ev][0]->Clone();
      TH1D *Nominal_met_gammaOnly_stats_up = (TH1D*) v_gammaData_nominal[ich][ev][0]->Clone();
      TH1D *Nominal_met_gamma_stats_down = (TH1D*) v_gammaData_nominal[ich][ev][0]->Clone();
      TH1D *Nominal_met_gammaOnly_stats_down = (TH1D*) v_gammaData_nominal[ich][ev][0]->Clone();
      binContent =0;
      for(unsigned int bin=1 ; bin < v_gammaData_nominal[ich][ev][0]->GetNbinsX(); bin++){
        binContent = v_gammaData_nominal[ich][ev][0]->GetBinContent(bin);
        if(binContent){ //Empty bins should also have a statistical error I guess...
					if(VERBOSE)std::cout<< "INFO : bin = " << bin << " content = " << binContent << " error = " << v_gammaData_nominal[ich][ev][0]->GetBinError(bin) << std::endl;
					if(VERBOSE)std::cout << "INFObis : down -- " << binContent- v_gammaData_nominal[ich][ev][0]->GetBinError(bin) << " --nominal-- " << binContent << " --up-- " << binContent+ v_gammaData_nominal[ich][ev][0]->GetBinError(bin) <<std::endl;
					Nominal_met_gamma_stats_up->SetBinContent(bin, binContent+ v_gammaData_nominal[ich][ev][0]->GetBinError(bin));
          Nominal_met_gammaOnly_stats_up->SetBinContent(bin, binContent+ v_gammaData_nominal[ich][ev][0]->GetBinError(bin));
          Nominal_met_gamma_stats_down->SetBinContent(bin, binContent- v_gammaData_nominal[ich][ev][0]->GetBinError(bin));
          Nominal_met_gammaOnly_stats_down->SetBinContent(bin, binContent- v_gammaData_nominal[ich][ev][0]->GetBinError(bin));
        }
      }

			//2.b Now that we have gamma up and down, compute the corresponding InstrMET up and down
			//

			if(v_gammaData_nominal[ich][ev][0]->GetEntries() ==0) continue;
			if(VERBOSE)std::cout<< "Gamma data found" << std::endl;
			counter_genuineMET=-1;
			for(unsigned int index=0; index < genuineMET.size(); index++){
				f_input->cd(genuineMET[index]); 
  			met_shape_tmp_genuineMET2.push_back((TH2F*) gDirectory->Get(tags_full+"_"+SystType+"_shapes"));	
  			current_position_genuineMET++;
				if(VERBOSE)std::cout <<"current_position_genuineMET" << current_position_genuineMET << std::endl;
				if(met_shape_tmp_genuineMET2[current_position_genuineMET]){
					if(VERBOSE)std::cout<< "genuineMET found" << std::endl;
      		v_genuineMET_nominal2[ich][ev].push_back( met_shape_tmp_genuineMET2[current_position_genuineMET]->ProjectionY(tags_full+genuineMET[index], projectionBin, projectionBin));

					if(VERBOSE)std::cout<< __LINE__ << std::endl;
					//remove genuineMET from Instr.MET
					counter_genuineMET++;

					TH1D *nominal_tmp = (TH1D*) v_genuineMET_nominal2[ich][ev][counter_genuineMET]->Clone();
					TH1D *nominal_tmp2 = (TH1D*) v_genuineMET_nominal2[ich][ev][counter_genuineMET]->Clone();
					if(nominal_tmp->GetEntries() == 0) continue;
					Nominal_met_gamma_stats_up->Add(nominal_tmp,-1);

					Nominal_met_gamma_stats_down->Add(nominal_tmp2,-1);

					if(VERBOSE)std::cout<< __LINE__ << std::endl;
				}
			}


			TH1D *Nominal_met_method_up2 = (TH1D*) Nominal_met_method_up->Clone();
			TH1D *Nominal_met_method_down2 = (TH1D*) Nominal_met_method_down->Clone();
			TH1D *Nominal_met_gamma_stats_up2 = (TH1D*) Nominal_met_gamma_stats_up->Clone();
			TH1D *Nominal_met_gamma_stats_down2 = (TH1D*) Nominal_met_gamma_stats_down->Clone();


			//Rebin to plot
			Nominal_met_gammaOnly_stats = (TH1D*) Nominal_met_gammaOnly_stats->Rebin(nmetAxis-1, "", metaxis);
			Nominal_met_gammaOnly_stats->Scale(1.0, "width");
			Nominal_met_gammaOnly_stats->Sumw2();
			Nominal_met_gamma_stats_up = (TH1D*) Nominal_met_gamma_stats_up->Rebin(nmetAxis-1, "", metaxis);
			Nominal_met_gamma_stats_up->Scale(1.0, "width");
			Nominal_met_gamma_stats_up->Sumw2();
			Nominal_met_gamma_stats_down = (TH1D*) Nominal_met_gamma_stats_down->Rebin(nmetAxis-1, "", metaxis);
			Nominal_met_gamma_stats_down->Scale(1.0, "width");
			Nominal_met_gamma_stats_down->Sumw2();
			Nominal_met_gammaOnly_stats_up = (TH1D*) Nominal_met_gammaOnly_stats_up->Rebin(nmetAxis-1, "", metaxis);
			Nominal_met_gammaOnly_stats_up->Scale(1.0, "width");
			Nominal_met_gammaOnly_stats_up->Sumw2();
			Nominal_met_gammaOnly_stats_down = (TH1D*) Nominal_met_gammaOnly_stats_down->Rebin(nmetAxis-1, "", metaxis);
			Nominal_met_gammaOnly_stats_down->Scale(1.0, "width");
			Nominal_met_gammaOnly_stats_down->Sumw2();
			Nominal_met_method_up = (TH1D*) Nominal_met_method_up->Rebin(nmetAxis-1, "", metaxis);
			Nominal_met_method_up->Scale(1.0, "width");
			Nominal_met_method_up->Sumw2();
			Nominal_met_method_down = (TH1D*) Nominal_met_method_down->Rebin(nmetAxis-1, "", metaxis);
			Nominal_met_method_down->Scale(1.0, "width");
			Nominal_met_method_down->Sumw2();



			//Save the shape up and down for method
      if(Nominal_met_method_up->GetEntries() ==0) continue;
      saved_cst_method_InstrMET_shape_Up[ich][ev].push_back((TH1D*) Nominal_met_method_up2->Clone());
      saved_cst_method_InstrMET_shape_Down[ich][ev].push_back((TH1D*) Nominal_met_method_down2->Clone());



			//Save the shape up and down for gamma stats 
      if(Nominal_met_gamma_stats_up->GetEntries() ==0) continue;
      saved_cst_stat_InstrMET_shape_Up[ich][ev].push_back((TH1D*) Nominal_met_gamma_stats_up2->Clone());
      saved_cst_stat_InstrMET_shape_Down[ich][ev].push_back((TH1D*) Nominal_met_gamma_stats_down2->Clone());



			//Style part
			std::vector< TH1D*> Nominal_met_new, Up_met_new, Down_met_new;

			Nominal_met_new.push_back((TH1D*) Instr_MET_nominal_rebin[ich][ev][0]->Clone());
			Nominal_met_new.push_back((TH1D*) Instr_MET_nominal_rebin[ich][ev][0]->Clone());
			Nominal_met_new.push_back((TH1D*) Nominal_met_gammaOnly_stats->Clone());

      Up_met_new.push_back((TH1D*) Nominal_met_method_up->Clone());
      Up_met_new.push_back((TH1D*) Nominal_met_gamma_stats_up->Clone());
      Up_met_new.push_back((TH1D*) Nominal_met_gammaOnly_stats_up->Clone());

      Down_met_new.push_back((TH1D*) Nominal_met_method_down->Clone());
      Down_met_new.push_back((TH1D*) Nominal_met_gamma_stats_down->Clone());
      Down_met_new.push_back((TH1D*) Nominal_met_gammaOnly_stats_down->Clone());

      TLegend *leg_met = new TLegend(0.6,0.7,0.89,0.89);
			std::vector<TString > corrections;
			corrections.push_back("Method");
			corrections.push_back("Gamma statistical fluctuations");
			corrections.push_back("Gamma sample only");

			if(VERBOSE)std::cout<< "Nominat_met size: " << Nominal_met_new.size() << std::endl;
			for(unsigned int index=0; index < Nominal_met_new.size(); index++){
				if(VERBOSE)std::cout << "In my custom loop" << std::endl;

				if(VERBOSE)std::cout << "Number of events = " << Nominal_met_new[index]->GetEntries() << std::endl;
				if(Nominal_met_new[index]->GetEntries() ==0) continue;

				//Nominal in black
				Nominal_met_new[index]->SetLineColor(kBlack);

				//Up in red
        Up_met_new[index]->SetLineColor(kRed);

				//Down in green
        Down_met_new[index]->SetLineColor(kGreen);

				//No stats
 				Nominal_met_new[index]->SetStats(kFALSE);


				//Make the legends

				leg_met->Clear();
        leg_met->AddEntry(Nominal_met_new[index], "Nominal", "l");
        leg_met->AddEntry(Up_met_new[index], corrections[index]+" up", "l");
        leg_met->AddEntry(Down_met_new[index], corrections[index]+" down", "l");






				cpos_fromTheory=canvasOldSize+index;
				if(VERBOSE)std::cout<< "position in Instr MET syst: " << cpos_fromTheory << std::endl;
		  	c_InstrMET_cst[cpos_fromTheory]->cd();
        //c[cpos]->SetBottomMargin(6);
        TPad *pad_met1 = new TPad("pad_met1","pad_met1",0,0.3,1,1);
        TPad *pad_met2 = new TPad("pad_met2","pad_met2",0,0,1,0.3);
        pad_met1->SetTopMargin(0.05);
        pad_met1->SetBottomMargin(0);
        pad_met1->Draw();
        pad_met2->SetTopMargin(0);
        pad_met2->SetBottomMargin(0.25);
        pad_met2->Draw();
        pad_met1->cd();
        TH1 *h_nominal_met = Nominal_met_new[index]->DrawCopy(); 
        h_nominal_met->SetTitle("");
        TH1 *h_up_met = Up_met_new[index]->DrawCopy("same"); 
        TH1 *h_down_met = Down_met_new[index]->DrawCopy("same"); 
        leg_met->Draw(); //same not needed?

				// Y axis ratio plot settings
        h_nominal_met->GetYaxis()->SetTitle("Events (a.u.)");
        h_nominal_met->GetYaxis()->SetNdivisions(505);
        h_nominal_met->GetYaxis()->SetTitleSize(20);
        h_nominal_met->GetYaxis()->SetTitleFont(43);
        h_nominal_met->GetYaxis()->SetTitleOffset(1.55);
        h_nominal_met->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
        h_nominal_met->GetYaxis()->SetLabelSize(15); 

        pad_met2->cd();
        TH1D *h_nominal_met2 = (TH1D*)Up_met_new[index]->Clone("h_nominal_met2"); //Be careful: non intuitive naming conventions to earn time.  
        TH1D *h_nominal_met3 = (TH1D*)Down_met_new[index]->Clone("h_nominal_met3");
        TH1D *h_up_met2 = (TH1D*)Nominal_met_new[index]->Clone("h_up_met2");
        TH1D *h_down_met2 = (TH1D*)Nominal_met_new[index]->Clone("h_down_met2");
        h_nominal_met2->Sumw2();
        h_nominal_met2->SetStats(0);
        h_nominal_met2->Divide(h_up_met2);
        h_nominal_met2->SetMarkerStyle(6);
        //h_nominal_met2->SetMarkerSize(20);
        h_nominal_met2->SetMarkerColor(kRed);
        h_nominal_met2->Draw("hist p e1");
        h_nominal_met3->Sumw2();
        h_nominal_met3->SetStats(0);
        h_nominal_met3->Divide(h_down_met2);
        h_nominal_met3->SetMarkerStyle(6);
        h_nominal_met3->SetMarkerColor(kGreen);
        h_nominal_met3->Draw("hist p e1 same");
        //Nominal_met->Draw("ep same");

  			// Ratio plot (h3) settings
   			h_nominal_met2->SetTitle(""); // Remove the ratio title

   			// Y axis ratio plot settings
   			h_nominal_met2->GetYaxis()->SetTitle("Up(Down)/Nominal ");
   			h_nominal_met2->GetYaxis()->SetNdivisions(505);
   			h_nominal_met2->GetYaxis()->SetTitleSize(20);
   			h_nominal_met2->GetYaxis()->SetTitleFont(43);
   			h_nominal_met2->GetYaxis()->SetTitleOffset(1.55);
   			h_nominal_met2->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   			h_nominal_met2->GetYaxis()->SetLabelSize(15);

   			//X axis ratio plot settings
   			h_nominal_met2->GetXaxis()->SetTitle(XaxisName);
   			h_nominal_met2->GetXaxis()->SetTitleSize(20);
   			h_nominal_met2->GetXaxis()->SetTitleFont(43);
   			h_nominal_met2->GetXaxis()->SetTitleOffset(4.);
   			h_nominal_met2->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   			h_nominal_met2->GetXaxis()->SetLabelSize(15);

				//Draw a line at 1
				//c1.Range( 0., -10., 1., 10. );
				TLine *line = new TLine(h_nominal_met2->GetXaxis()->GetXmin(),1.,h_nominal_met2->GetXaxis()->GetXmax() ,1.);
				line->SetLineColor(kBlack);
				line->SetLineWidth(1);
				line->SetLineStyle(3);
				line->Draw();


				if(VERBOSE)std::cout << "In my custom loop: end!" << std::endl;



				f_output->cd(globalDirectoryName+"/"+"Instr_MET_cst_"+chTags[ich]);
				c_InstrMET_cst[cpos_fromTheory]->Write();

			}


		}
	}

	//
	//Final wrap up of the uncertainties for Instr. MET !
	//
	if(VERBOSE)std::cout<< "In the very last wrapping up loop" << std::endl;
	square_sum_up=0;
	square_sum_down=0;
	syst_binContent_nominal = 0;
	std::vector<double> syst_binContent_genuineMET_up(samples.size(), 0);
	double syst_binContent_method_up = 0;
	double syst_binContent_stat_up = 0;
	std::vector<double> syst_binContent_genuineMET_down(samples.size(), 0);
	double syst_binContent_method_down = 0;
	double syst_binContent_stat_down = 0;
	num_canvas =0;
	updatedBinContent_up = 0;
	updatedBinContent_down = 0;

  std::vector< TDirectory*> directory_Instr_MET_FINAL;
  for ( unsigned int ich = 0; ich < chTags.size() ; ich++){
		globalDirectory->cd();
    directory_Instr_MET_FINAL.push_back(f_output->mkdir(globalDirectoryName+"/"+"Instr_MET_FINAL_"+chTags[ich]));
  }
  if(VERBOSE)std::cout<< "Just before the loop" << std::endl;
  for ( unsigned int ich = 0; ich < chTags.size() ; ich++){
    for ( unsigned ev = 0; ev < evCat.size() ; ev++){
			TString tags_full=chTags[ich]+evCat[ev];

			TH1D *h_nominal_finalSum = (TH1D*) Instr_MET_nominal[ich][ev][0]->Clone();	
			TH1D *h_nominal_finalSum_up = (TH1D*) Instr_MET_nominal[ich][ev][0]->Clone();	
			TH1D *h_nominal_finalSum_down = (TH1D*) Instr_MET_nominal[ich][ev][0]->Clone();	

			if(VERBOSE)std::cout<< "A bit further" << std::endl;

			for(unsigned int bin=1; bin < Instr_MET_nominal[ich][ev][0]->GetNbinsX(); bin++){
				if(VERBOSE)std::cout<< __LINE__ << std::endl;
				syst_binContent_nominal = Instr_MET_nominal[ich][ev][0]->GetBinContent(bin); 
				if(VERBOSE)std::cout<< __LINE__ << std::endl;
				//genuineMET
				for(unsigned int s =0; s < samples.size(); s++){
				if(VERBOSE)std::cout<< __LINE__ << std::endl;
				if(VERBOSE)std::cout<< samples[s] << std::endl;
					syst_binContent_genuineMET_up[s] = savedInstrMET_shape_Up[s][ich][ev][0]->GetBinContent(bin); 
				if(VERBOSE)std::cout<< __LINE__ << std::endl;
					syst_binContent_genuineMET_down[s] = savedInstrMET_shape_Down[s][ich][ev][0]->GetBinContent(bin); 
				}
				if(VERBOSE)std::cout<< __LINE__ << std::endl;
				//method
				syst_binContent_method_up = saved_cst_method_InstrMET_shape_Up[ich][ev][0]->GetBinContent(bin); 
				syst_binContent_method_down = saved_cst_method_InstrMET_shape_Down[ich][ev][0]->GetBinContent(bin); 
				//gammaData stat
				syst_binContent_stat_up = saved_cst_stat_InstrMET_shape_Up[ich][ev][0]->GetBinContent(bin); 
				syst_binContent_stat_down = saved_cst_stat_InstrMET_shape_Down[ich][ev][0]->GetBinContent(bin); 
			
				if(VERBOSE)std::cout<< __LINE__ << std::endl;

				for(unsigned int s = 0; s < samples.size(); s++){
					if(ALL || (samples[s].Contains("W#gamma") && WGAMMA_ONLY) || (samples[s].Contains("Z#gamma") && ZGAMMA_ONLY) || (samples[s].Contains("W#rightarrow") && WJETS_ONLY) )square_sum_up += 1.0*(syst_binContent_genuineMET_up[s] - syst_binContent_nominal)*(syst_binContent_genuineMET_up[s] - syst_binContent_nominal);
				}
				if(ALL || METHOD_ONLY) square_sum_up += 1.0*(syst_binContent_method_up - syst_binContent_nominal)*(syst_binContent_method_up - syst_binContent_nominal);
				if(ALL || GAMMASTATS_ONLY) square_sum_up += 1.0*(syst_binContent_stat_up - syst_binContent_nominal)*(syst_binContent_stat_up - syst_binContent_nominal);
				if(VERBOSE)std::cout<< __LINE__ << std::endl;
				for(unsigned int s = 0; s < samples.size(); s++){
				  if(ALL || (samples[s].Contains("W#gamma") && WGAMMA_ONLY) || (samples[s].Contains("Z#gamma") && ZGAMMA_ONLY) || (samples[s].Contains("W#rightarrow") && WJETS_ONLY) )	square_sum_down += 1.0*(syst_binContent_genuineMET_down[s] - syst_binContent_nominal)*(syst_binContent_genuineMET_down[s] - syst_binContent_nominal);
        }
        if(ALL || METHOD_ONLY) square_sum_down += 1.0*(syst_binContent_method_down - syst_binContent_nominal)*(syst_binContent_method_down - syst_binContent_nominal);
        if(ALL || GAMMASTATS_ONLY) square_sum_down += 1.0*(syst_binContent_stat_down - syst_binContent_nominal)*(syst_binContent_stat_down - syst_binContent_nominal);

				if(VERBOSE)std::cout<< "Bin = " << bin << " et sqrt(square_sum_up) = " << sqrt(square_sum_up) << std::endl;

				if(h_nominal_finalSum_up->GetBinContent(bin) !=0 && h_nominal_finalSum_down->GetBinContent(bin) !=0){
					updatedBinContent_up = h_nominal_finalSum_up->GetBinContent(bin) + sqrt(square_sum_up);
          updatedBinContent_down = h_nominal_finalSum_down->GetBinContent(bin) - sqrt(square_sum_down);
          if(VERBOSE)std::cout<< "Case: "<<tags_full<< " UpdatedBinContent_up = " << updatedBinContent_up << " and updatedBinContent_down = " << updatedBinContent_down << std::endl;
          h_nominal_finalSum_up->SetBinContent(bin, updatedBinContent_up);
          h_nominal_finalSum_down->SetBinContent(bin, updatedBinContent_down);
        }
        square_sum_up=0;
        square_sum_down=0;

			}



			//Style part

			TH1D *Nominal_met_new = (TH1D*) h_nominal_finalSum->Rebin(nmetAxis-1, "Nominal_met_new", metaxis);
			Nominal_met_new->Scale(1.0, "width");
			TH1D *Up_met_new = (TH1D*) h_nominal_finalSum_up->Rebin(nmetAxis-1, "Up_met_new", metaxis);
			Up_met_new->Scale(1.0, "width");
			TH1D *Down_met_new = (TH1D*) h_nominal_finalSum_down->Rebin(nmetAxis-1, "Down_met_new", metaxis);
			Down_met_new->Scale(1.0, "width");


			if(VERBOSE)std::cout << "Number of events = " << Nominal_met_new->GetEntries() << std::endl;
      if(Nominal_met_new->GetEntries() ==0) continue;



			//Nominal in black
			Nominal_met_new->SetLineColor(kBlack);

			//Up in red
      Up_met_new->SetLineColor(kRed);

			//Down in green
      Down_met_new->SetLineColor(kGreen);

			//No stats
 			Nominal_met_new->SetStats(kFALSE);


			//Make the legends

      TLegend *leg_met = new TLegend(0.6,0.7,0.89,0.89);
      leg_met->AddEntry(Nominal_met_new, "Nominal", "l");
      leg_met->AddEntry(Up_met_new, "full corrections up", "l");
      leg_met->AddEntry(Down_met_new, "full corrections down", "l");






			if(VERBOSE)std::cout<< "ready to fill" << std::endl;

			c_final_wrapUp.push_back( new TCanvas(tags_full+"_FINAL_RESULT_met_shapes", tags_full+"_FINAL_RESULT_met_shapes", 800, 800));
      c_final_wrapUp[num_canvas]->cd();
			//c[cpos]->SetBottomMargin(6);
      TPad *pad_met1 = new TPad("pad_met1","pad_met1",0,0.3,1,1);
      TPad *pad_met2 = new TPad("pad_met2","pad_met2",0,0,1,0.3);
      pad_met1->SetTopMargin(0.05);
      pad_met1->SetBottomMargin(0);
      pad_met1->Draw();
      pad_met2->SetTopMargin(0);
      pad_met2->SetBottomMargin(0.25);
      pad_met2->Draw();
      pad_met1->cd();
      TH1 *h_nominal_met = Nominal_met_new->DrawCopy(); 
      h_nominal_met->SetTitle("");
      TH1 *h_up_met = Up_met_new->DrawCopy("same"); 
      TH1 *h_down_met = Down_met_new->DrawCopy("same"); 

      leg_met->Draw(); //same not needed?

      // Y axis ratio plot settings
      h_nominal_met->GetYaxis()->SetTitle("Events (a.u.)");
      h_nominal_met->GetYaxis()->SetNdivisions(505);
      h_nominal_met->GetYaxis()->SetTitleSize(20);
      h_nominal_met->GetYaxis()->SetTitleFont(43);
      h_nominal_met->GetYaxis()->SetTitleOffset(1.55);
      h_nominal_met->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
      h_nominal_met->GetYaxis()->SetLabelSize(15); 

      pad_met2->cd();
      TH1D *h_nominal_met2 = (TH1D*)Up_met_new->Clone("h_nominal_met2"); //Be careful: non intuitive naming conventions to earn time.  
      TH1D *h_nominal_met3 = (TH1D*)Down_met_new->Clone("h_nominal_met3");
      TH1D *h_up_met2 = (TH1D*)Nominal_met_new->Clone("h_up_met2");
      TH1D *h_down_met2 = (TH1D*)Nominal_met_new->Clone("h_down_met2");
      h_nominal_met2->Sumw2();
      h_nominal_met2->SetStats(0);
      h_nominal_met2->Divide(h_up_met2);
      h_nominal_met2->SetMarkerStyle(6);
      //h_nominal_met2->SetMarkerSize(20);
      h_nominal_met2->SetMarkerColor(kRed);
      h_nominal_met2->Draw("hist p e1");
      h_nominal_met3->Sumw2();
      h_nominal_met3->SetStats(0);
      h_nominal_met3->Divide(h_down_met2);
      h_nominal_met3->SetMarkerStyle(6);
      h_nominal_met3->SetMarkerColor(kGreen);
      h_nominal_met3->Draw("hist p e1 same");
      //Nominal_met->Draw("ep same");

  		// Ratio plot (h3) settings
   		h_nominal_met2->SetTitle(""); // Remove the ratio title

   		// Y axis ratio plot settings
   		h_nominal_met2->GetYaxis()->SetTitle("Up(Down)/Nominal ");
   		h_nominal_met2->GetYaxis()->SetNdivisions(505);
   		h_nominal_met2->GetYaxis()->SetTitleSize(20);
   		h_nominal_met2->GetYaxis()->SetTitleFont(43);
   		h_nominal_met2->GetYaxis()->SetTitleOffset(1.55);
   		h_nominal_met2->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   		h_nominal_met2->GetYaxis()->SetLabelSize(15);

   		//X axis ratio plot settings
   		h_nominal_met2->GetXaxis()->SetTitle(XaxisName);
   		h_nominal_met2->GetXaxis()->SetTitleSize(20);
   		h_nominal_met2->GetXaxis()->SetTitleFont(43);
   		h_nominal_met2->GetXaxis()->SetTitleOffset(4.);
   		h_nominal_met2->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   		h_nominal_met2->GetXaxis()->SetLabelSize(15);

			//Draw a line at 1
			//c1.Range( 0., -10., 1., 10. );
			TLine *line = new TLine(h_nominal_met2->GetXaxis()->GetXmin(),1.,h_nominal_met2->GetXaxis()->GetXmax() ,1.);
			line->SetLineColor(kBlack);
			line->SetLineWidth(1);
			line->SetLineStyle(3);
			line->Draw();




			if(VERBOSE)std::cout << __LINE__ << std::endl;




			f_output->cd(globalDirectoryName+"/"+"Instr_MET_FINAL_"+chTags[ich]);
			c_final_wrapUp[num_canvas]->Write();
			num_canvas++;
			//write at the root the shapes up and down with convention of the runPlotter
			//Convention is, for example, eevbf_metfinal_InstrMET_up : WARNING, we look only at metfinal and mtfinal since we have a cut at 125 GeV
			f_output->cd();
			TString naming = "ERROR";
			if(projectionBin == 17) naming="final_";
			if(projectionBin == 1) naming="_";
			TH1D* h_temp_up = (TH1D*) Up_met_new->Clone(tags_full+"_"+SystType+naming+"InstrMET_shape_up");
			TH1D* h_temp_down = (TH1D*) Down_met_new->Clone(tags_full+"_"+SystType+naming+"InstrMET_shape_down");
			h_temp_up->Write();
			h_temp_down->Write();

			//Also write the shapes with a total number of events as a Y-axis and not just a #events/GeV with the convention of the computeLimits code
			if(SystType == "mt" && projectionBin == 17){
				f_output->cd();
				TH1D* h_temp_up_absolute = (TH1D*) h_nominal_finalSum_up->Clone(tags_full+"_"+SystType+"_InstrMET_absolute_shape_up");
				TH1D* h_temp_down_absolute = (TH1D*) h_nominal_finalSum_down->Clone(tags_full+"_"+SystType+"_InstrMET_absolute_shape_down");
				//Rebin like computeLimits expect (be careful, here we want a total number of events and not just a #events/GeV)
				h_temp_up_absolute = (TH1D*) h_temp_up_absolute->Rebin(nmetAxis_limits.at(evCat[ev])-1, "", metaxis_limits.at(evCat[ev]));
				h_temp_down_absolute = (TH1D*) h_temp_down_absolute->Rebin(nmetAxis_limits.at(evCat[ev])-1, "", metaxis_limits.at(evCat[ev]));
				h_temp_up_absolute->Write();
    		h_temp_down_absolute->Write();
			}









		}
	}
}










void MakeSyst_custom_forMTandMET(){
	TString path = std::string(getenv("CMSSW_BASE")) + "/src/UserCode/llvv_fwk/test/hzz2l2v/";
  TFile *f_input = TFile::Open(path+"plotter.root");
	TFile *f_output = new TFile(path+"/macro_Intr_MET_Syst/InstrMET_systematics.root","RECREATE");
	MakeSyst_custom(f_input, f_output, "mt", 17); //produce up and down variations for the mt_shapes (for limits and plotter) -- this is for the final version, with a MET cut of 125GeV
	//MakeSyst_custom(f_input, f_output, "met", 17);//produce up and down variations for the met_shapes (for the plotter) -- this is for the final version, with a MET cut of 125GeV
	//MakeSyst_custom(f_input, f_output, "mt", 1); //produce up and down variations for the mt_shapes (for the plotter) -- this is for the final version, with no MET cut
	//MakeSyst_custom(f_input, f_output, "met", 1);//produce up and down variations for the met_shapes (for the plotter) -- this is for the final version, with no MET cut

	
	
	//MakeSyst_custom(f_input, f_output, "mt_Inbveto50"); //produce up and down variation for Inbveto50 control plots
	//MakeSyst_custom(f_input, f_output, "mt_Inbveto80"); //produce up and down variation for Inbveto80 control plots
	//MakeSyst_custom(f_input, f_output, "mt_Inbveto125");//produce up and down variation for Inbveto125 control plots
}
