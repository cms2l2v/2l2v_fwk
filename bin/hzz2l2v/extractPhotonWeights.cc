#include <vector>
#include <list>
#include <iterator>
#include <algorithm>
#include <unordered_map>
#include <iostream>
#include <string>
#include <regex>

#include "UserCode/llvv_fwk/interface/tdrstyle.h"
#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/RootUtils.h"
#include "UserCode/llvv_fwk/interface/JSONWrapper.h"

using namespace std;

std::string function_used = "([0] + [1]*x + [2]*x^2 + [3]*x^3 + [4]*x^4 )";

TF1* FitFunction(TH1D *hf, bool cutOnBinContent){

        TH1D *loghf = new TH1D("","", 1500, 0, 1500);
        int  ncells = hf->GetSize();

        double min = 55;
        double max = 0;

        for(int bin=0;bin<(ncells-1);bin++){
           if(hf->GetBinContent(bin)>0){ 
               loghf->SetBinContent( bin, log10(hf->GetBinContent(bin)));
               loghf->SetBinError( bin, (hf->GetBinError(bin)/hf->GetBinContent(bin)));	
               if( hf->GetBinContent(bin)>1 && cutOnBinContent ){ max = hf->GetXaxis()->GetBinCenter(bin); }
	       else if( hf->GetBinContent(bin)>15 && !cutOnBinContent ){ max = hf->GetXaxis()->GetBinCenter(bin); }
           }
        }
        TF1 *f1 = new TF1( "fit1", function_used.c_str(), min, max);
        loghf->Fit( f1, "R");
        TF1 *fit = loghf->GetFunction("fit1");

	TF1 *ffit1 = new TF1( "ffit1", function_used.c_str(), min, max);
        ffit1->SetParameters( fit->GetParameter(0), fit->GetParameter(1), fit->GetParameter(2), fit->GetParameter(3), fit->GetParameter(4));

	return ffit1;   
}

TH1D* SetLogHisto(TH1D *hf){

        TH1D *loghf = new TH1D("","", 1500, 0, 1500);
        int  ncells = hf->GetSize();

        for(int bin=0;bin<(ncells-1);bin++){
           if(hf->GetBinContent(bin)>0){
               loghf->SetBinContent( bin, log10(hf->GetBinContent(bin)));
               loghf->SetBinError( bin, (hf->GetBinError(bin)/hf->GetBinContent(bin)));
           }
        }
    
        return loghf;
}

TH1D* FillHistoWgts(TF1* flep, TF1* fphot, bool cutOnBinContent){
    TH1D *hweights = new TH1D( "", "", 450, 50,500);
    int n=0;
    for(int k=5;k<450;k++){
	double mom = (50 + k);
        if( std::pow(10, flep->Eval(mom))<1 && cutOnBinContent ){ n++; }
        else if( std::pow(10, flep->Eval(mom))<15 && !cutOnBinContent ){ n++; }
	if(n>0) continue;
        hweights->SetBinContent(k, std::pow(10,flep->Eval(mom))/std::pow(10,fphot->Eval(mom)));
	hweights->SetBinError(k,0);
    }
    return hweights;
}  

int main(int argc, char* argv[]){
   setTDRStyle();  
   gStyle->SetPadTopMargin   (0.06);
   gStyle->SetPadBottomMargin(0.12);
   gStyle->SetPadRightMargin (0.16);
   gStyle->SetPadLeftMargin  (0.14);
   gStyle->SetTitleSize(0.04, "XYZ");
   gStyle->SetTitleXOffset(1.1);
   gStyle->SetTitleYOffset(1.45);
   gStyle->SetPalette(1);
   gStyle->SetNdivisions(505);

   //FIXME these are currently hard-coded, but we should pass them via parameters using the code just bellow
   bool isData=true;
   bool bfit=false;
   int  rbin=1;
   bool asym=false;
   int rebin=1;
   string gDataFile = "plotter_2016_01_18_forPhotonWeights.root";
   string zDataFile = "plotter_2016_01_18_forPhotonWeights.root";
   string outDir    = "photonWeights/";
   string outFile   = outDir + "photonWeights.root";

   for(int i=1;i<argc;i++){
     string arg(argv[i]);
     if(arg.find("--help")!=string::npos){
        printf("--help    --> print this helping text\n");
        printf("--inFile  --> path of the plotter .root file with inputs\n");
        printf("--outDir  --> path of the output directory to save the plots\n");
        printf("--outFile --> path of the output summary .root file\n");
        printf("--fitf    --> use fitting functions to compute weights\n");
	return 0;
     }
     if(arg.find("--inFile" )!=string::npos && i+1<argc){ zDataFile= argv[i+1];  gDataFile= argv[i+1]; i++; printf("input file = %s\n", zDataFile.c_str()); }
     if(arg.find("--outFile")!=string::npos && i+1<argc){ outFile  = argv[i+1];  i++; printf("output file = %s\n", outFile.c_str()); }
     if(arg.find("--outDir" )!=string::npos && i+1<argc){ outDir   = argv[i+1];  i++; printf("outDir = %s\n", outDir.c_str());  }
     if(arg.find("--fitf")!=string::npos && i+1<argc){ bfit = argv[i+1];  i++; }
   }
   system( (string("mkdir -p ") + outDir).c_str());


  std::vector<string> gDataDir;
  std::vector<string> zDataDir;
  if(isData){
     gDataDir.push_back("#gamma data_reweighted");
     zDataDir.push_back("data");
  }else{
     gDataDir.push_back("GJets_HT-40to100");
     gDataDir.push_back("GJets_HT-100to200");
     gDataDir.push_back("GJets_HT-200to400");
     gDataDir.push_back("GJets_HT-400to600");
     gDataDir.push_back("GJets_HT-600toInf");
     zDataDir.push_back("Z#rightarrow ll");
  }

  std::vector<int> catColor = {4     , 2   , 8   , 1    };
  std::vector<string> catL  = {"#mu#mu", "ee", "ll", "#gamma"};
  std::vector<string> binL  = {"=0jet", "#geq1jets", "vbf"};

  std::vector<string> cat   = {"mumu", "ee", "ll", "gamma"};
  std::vector<string> bin   = {"eq0jets","geq1jets","vbf"};
  std::vector<string> var   = {"_qt", "_qmass", "_met", "_mt"};

  std::map<string, TH1D*> DataHistos;
  std::map<string, TF1* > FitFunctionMap;
  std::map<string, TH1D*> WeightsFitFunction;
 
  for(unsigned int c=0;c<cat.size();c++){
  for(unsigned int b=0;b<bin.size();b++){
  for(unsigned int v=0;v<var.size();v++){
     string& DataFile             = (cat[c]=="gamma")?gDataFile:zDataFile;
     std::vector<string>& DataDir = (cat[c]=="gamma")?gDataDir :zDataDir;

     TFile* File = new TFile(DataFile.c_str(), "READ");
     if(!File || File->IsZombie() || !File->IsOpen() || File->TestBit(TFile::kRecovered) ){
        printf("File %s coudn't be opened\n", DataFile.c_str());
     }

     for(unsigned int d=0;d<DataDir.size();d++){
        TH1D* hist = NULL;
//        if(cat[c]=="photon"){ //read photon info form ee channel
//           hist = (TH1D*) File->Get((DataDir[d]+"/"+"ee"+bin[b]+var[v]).c_str());        
//        }else{
           hist = (TH1D*) File->Get((DataDir[d]+"/"+cat[c]+bin[b]+var[v]).c_str());        
//        }
        if(!hist)continue;

        hist->GetSumw2();
        if(rbin>1){ hist->Rebin(rbin);  hist->Scale(1.0, "width"); }
        hist->GetXaxis()->SetTitleSize(.055);
        hist->GetYaxis()->SetTitleSize(.055);
        hist->GetXaxis()->SetLabelSize(.05);
        hist->GetYaxis()->SetLabelSize(.05);
        hist->SetFillColor(0);
        hist->SetLineColor(catColor[c]);
        hist->SetMarkerColor(catColor[c]);
        hist->SetStats(kFALSE);

        if(DataHistos.find(cat[c]+bin[b]+var[v])==DataHistos.end()){
           gROOT->cd(); //make sure that the file is saved in memory and not in file
           DataHistos[cat[c]+bin[b]+var[v]] = (TH1D*)hist->Clone(); //create a new histo, since it's not found
        }else{
           DataHistos[cat[c]+bin[b]+var[v]]->Add(hist); //add to existing histogram
        }
     
        //Compute the fitting function for evey histogram in eq0jets and geq1jets, in ee-mumu-gamma channel. For VBF, we used geq1jets fuctions reweighted in such a way the integral is still ok.
        if( var[v]!="_qt" ) continue; //Considering only qt distribution
	if( bin[b] == "vbf" && ( cat[c] == "ee" || cat[c] == "mumu" || cat[c] == "ll" )){
		TH1D *hvbf = (TH1D*) File->Get((DataDir[d]+"/"+cat[c]+"vbf"+var[v]).c_str());
		hist = (TH1D*) File->Get((DataDir[d]+"/"+cat[c]+"geq1jets"+var[v]).c_str());	
		hist->Scale(hvbf->Integral()/hist->Integral());
		FitFunctionMap[cat[c]+"vbf"+var[v]] = FitFunction( hist, true);
	}else if( bin[b] !="vbf" || cat[c] == "gamma" ){
                FitFunctionMap[cat[c]+bin[b]+var[v]] =  FitFunction( hist, false);
	}
     } 

     File->Close();
  }}}//all histos are now loaded

  
  //non fixed-width rebins
  if(asym){
     double xbins[] = {55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110, 120, 130, 140, 150, 165, 250, 300, 400, 500, 700, 1000, 1500};    
     for(unsigned int v=0;v<var.size();v++){
        if(var[v]!="_qt")continue;
        for(unsigned int b=0;b<bin.size();b++){
           for(unsigned int c=0;c<cat.size();c++){
              if(DataHistos.find(cat[c]+bin[b]+var[v])==DataHistos.end()){continue;}
              TH1D* hist = DataHistos[cat[c]+bin[b]+var[v]];        
              hist = (TH1D*)hist->Clone((cat[c]+bin[b]+var[v]+"rebin").c_str());
              hist = (TH1D*)hist->Rebin((sizeof(xbins)/sizeof(double))-1, hist->GetName(), xbins);             
              hist->Scale(1.0, "width");
              DataHistos[cat[c]+bin[b]+var[v]+"rebin"] = hist;
           }
        }
        printf("add variable %s\n", (var[v]+"rebin").c_str());
        var.push_back(var[v]+"rebin");
     }
  }
 
  //compute the weight
  for(unsigned int v=0;v<var.size();v++){
     if(var[v]!="_qt" && var[v]!="_qtrebin")continue;
     for(unsigned int b=0;b<bin.size();b++){     
	if(DataHistos.find(string("gamma")+bin[b]+var[v])==DataHistos.end()){printf("Histo missing for %s\n", (string("gamma")+bin[b]+var[v]).c_str()); continue;}
        TH1D* histgamma = DataHistos[string("gamma")+bin[b]+var[v]];
	if(FitFunctionMap.find(string("gamma")+bin[b]+var[v])==FitFunctionMap.end()){printf("Function missing for %s\n", (string("gamma")+bin[b]+var[v]).c_str()); continue;}
        TF1 *fphot = FitFunctionMap[string("gamma")+bin[b]+var[v]];

        for(unsigned int c=0;c<cat.size();c++){
           if(cat[c]=="gamma")continue;
	    
           if(bfit){
	       if(FitFunctionMap.find(cat[c]+bin[b]+var[v])==FitFunctionMap.end()){printf("Function missing for %s\n", (cat[c]+bin[b]+var[v]).c_str()); continue;}
	       TF1  *flep = FitFunctionMap[cat[c]+bin[b]+var[v]]; 
	       if( bin[b] == "vbf" && ( cat[c] == "ee" || cat[c] == "mumu" || cat[c] == "ll" )){
                   WeightsFitFunction[ cat[c]+bin[b]+var[v]+"weight"]= FillHistoWgts(flep,fphot, true);
	       } else if( bin[b] !="vbf" || cat[c] == "gamma" ){
                   WeightsFitFunction[ cat[c]+bin[b]+var[v]+"weight"]= FillHistoWgts(flep,fphot, false);
	       }
	        
	   } else if(!bfit){
               if(DataHistos.find(cat[c]+bin[b]+var[v])==DataHistos.end()){printf("Histo missing for %s\n", (cat[c]+bin[b]+var[v]).c_str()); continue;}
               TH1D* hist = DataHistos[cat[c]+bin[b]+var[v]];
               hist = (TH1D*)hist->Clone((cat[c]+bin[b]+var[v]+"weight").c_str());
               hist->Divide(histgamma);
               DataHistos[cat[c]+bin[b]+var[v]+"weight"] = hist;
           }
        }
     }
     printf("add variable %s\n", (var[v]+"weight").c_str());
     var.push_back(var[v]+"weight");
  }

  //make the plots
  if(!bfit){
  for(unsigned int v=0;v<var.size();v++){
     double xmin,xmax;
     double ymin=0.5, ymax=1E6;
     if (var[v]=="_qt" || var[v]=="_qtrebin") {                    xmin=55.00; xmax=1000.00;
     } else if (var[v]=="_qtweight" || var[v]=="_qtrebinweight") { xmin=55.00; xmax=1000.00; ymin=0.0001;  ymax=0.5;
     } else if( (var[v]=="_met") || var[v]=="_mt") {               xmin=1.; xmax=1000.0;
     } else{                                                       xmin=1.; xmax=1000.0;    
     }

     TCanvas* c1 =new TCanvas("c1","c1",500*bin.size(), 500);
     c1->Divide(bin.size(),1);

     for(unsigned int b=0;b<bin.size();b++){
        c1->cd(1+b);
        if(var[v]!="_qtweight" || var[v]!="_qtrebinweight") gPad->SetLogy(); if(var[v]=="_qt" || var[v]=="_qtrebin")gPad->SetLogx();
        TLegend* leg = new TLegend(0.52,0.67,0.92,0.90);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.04);
        leg->SetTextFont(42);
        leg->SetHeader(binL[b].c_str());

        bool axisDrawn=false;
        for(unsigned int c=0;c<cat.size();c++){
           if(cat[c]=="ll")continue;  //do not ll on the plots
           if(DataHistos.find(cat[c]+bin[b]+var[v])==DataHistos.end()){printf("Histo missing for %s\n", (cat[c]+bin[b]+var[v]).c_str()); continue;}
           TH1D* hist = DataHistos[cat[c]+bin[b]+var[v]];
           hist->GetXaxis()->SetRangeUser(xmin,xmax);
           hist->GetYaxis()->SetRangeUser(ymin,ymax);

           if(!axisDrawn){hist->Draw("AXIS"); axisDrawn=true;}
           hist->Draw("HIST E1 same");           
           leg->AddEntry(hist, catL[c].c_str(), "LP"); 
       }
       leg->Draw("SAME");
    }
    c1->SaveAs((outDir + "/PhotonWeight"+var[v]+".png").c_str());
   }
   }

  //Save Fit
  if(bfit){ 
    for(unsigned int v=0;v<var.size();v++){
       if(var[v]=="_qtrebin" || var[v]=="_qtweight" || var[v]=="_qtrebinweight" ) continue;
       for(unsigned int c=0;c<cat.size();c++){
	  string canvas_name;
          TCanvas* c1 =new TCanvas("c1","c1",500*bin.size(), 500);
          c1->Divide(bin.size(),1);
          for(unsigned int b=0;b<bin.size();b++){
	     c1->cd(1+b);
             TLegend* leg = new TLegend(0.52,0.67,0.92,0.90);
             leg->SetFillStyle(0);
             leg->SetBorderSize(0);
             leg->SetTextSize(0.04);
             leg->SetTextFont(42);
             leg->SetHeader(binL[b].c_str());
	     TH1D* hist = NULL;
             if(var[v]=="_qt"){
		TF1 *f1 = NULL;
		hist = SetLogHisto(DataHistos[cat[c]+bin[b]+var[v]]); 
		f1 = FitFunctionMap[cat[c]+bin[b]+var[v]];
		f1->SetParameter(0, rebin*f1->GetParameter(0));
                f1->SetLineColor(kRed);
                f1->SetLineWidth(2);
	        hist->SetStats(false);
		hist->Rebin(rebin);
		hist->Draw("AP");
		f1->Draw("same");
        	hist->GetXaxis()->SetTitle("q_{T} [GeV]");
		hist->GetXaxis()->SetRangeUser(55,500); 
		hist->GetXaxis()->SetTitleOffset(0.8);
        	hist->GetYaxis()->SetTitle("Events");
		hist->GetYaxis()->SetTitleOffset(0.8); 
		canvas_name = "/FittingFunction_"+cat[c]+var[v];
	     } else {
		hist = DataHistos[cat[c]+bin[b]+var[v]];
                hist->GetXaxis()->SetRangeUser(  1., 1000);
                hist->GetYaxis()->SetRangeUser( 0.5,  1E6);
		hist->Draw("HIST E1 same`");
		canvas_name = "/PhotonDistribution_"+cat[c]+var[v];
	     }
             leg->AddEntry(hist, catL[c].c_str(), "LP");
	     leg->Draw("SAME");
	  } 
          c1->SaveAs((outDir + canvas_name +".png").c_str());
       }
    }
  }

  //SAVE THE material to the weight file
  TFile* OutputFile = new TFile(outFile.c_str(),"RECREATE");
  OutputFile->cd();
  for(unsigned int c=0;c<cat.size();c++){
     if(cat[c]=="gamma")continue;
     for(unsigned int b=0;b<bin.size();b++){
     for(unsigned int v=0;v<var.size();v++){
             if(!bfit){ if(DataHistos.find(cat[c]+bin[b]+var[v])==DataHistos.end()){printf("Histo missing for %s\n", (cat[c]+bin[b]+var[v]).c_str()); continue;} }

             TH1D* hist = DataHistos[cat[c]+bin[b]+var[v]];
             if(bfit){
                 if(var[v]=="_qtweight"){ 
			   TH1D* histweight = WeightsFitFunction[cat[c]+bin[b]+"_qtweight"]; 
			   TGraphErrors* graph = new TGraphErrors(histweight);
			   graph->Write((cat[c]+bin[b]+"_qt_datafitfunctionwgts").c_str()); 
                 }
             } else if(!bfit) { 
             	 if(var[v]=="_qtweight"){
                	   TGraphErrors* graph = new TGraphErrors(hist);
                	   graph->Write((cat[c]+bin[b]+"_qt_datafitwgts").c_str());
            	 }else if(var[v]=="_qtrebinweight"){
                 	   TGraphErrors* graph = new TGraphErrors(hist);
                	   graph->Write((cat[c]+bin[b]+"_qt_datafitwgtsrebin").c_str()); 
             	}
             }

             if(var[v]=="_qmass"){
                     //replace VBF by geq1jets has there is no stat in VBF for Z candidates
                     if(bin[b]=="vbf")continue;  
                     if(bin[b]=="geq1jets"){
                          hist->Write((cat[c]+"vbf"+"_zmass").c_str());  
                     }
                     hist->Write((cat[c]+bin[b]+"_zmass").c_str());
             }
        }
     }
  }
  OutputFile->Write(); OutputFile->Close();
}

