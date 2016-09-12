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

//std::string function_used = "([0] + [1]*x + [2]*pow(x,[3]) + [4]*pow(x,[5]) )";
//std::string function_used = "[0] + [1]*TMath::Exp(-(log(x + [2])-[3])^2/[4])/[4]*x";
std::string function_used = "[0]*pow(x,[1]) + [2]*pow(x,[3]) + [4]";

TF1* FitFunction(TH1D *hf, double SF=1.0, std::string cat=" ", std::string bin=" "){

  if(SF!=1.0)printf("SF = %f\n", SF);        


        TH1D* loghf = (TH1D*)hf->Clone("TmpForFitting");
        loghf->Reset();
        int  ncells = hf->GetSize();

        double par[5] = {0.,0.,0.,0.,0.};
	if( cat == "eq0jets" && bin == "all" ){ par[0]=2.5; par[1]=2.5; par[2]=2.5; par[3]=1.5; par[4]=2.5; }
	else if( cat == "eq0jets"  && bin == "ee"    ){ par[0]=-2.5; par[1]=-1.5;  par[2]=2.5; par[3]=2.5; par[4]=2.5; }
	else if( cat == "eq0jets"  && bin == "mumu"  ){  par[0]=2.5;  par[1]=1.5;  par[2]=2.5; par[3]=3.5; par[4]=1.5; }
        
	else if( cat == "geq1jets" && bin == "all"   ){ par[0]=-2.5; par[1]=-3.5;  par[2]=2.5; par[3]=1.5; par[4]=4.5; }
        else if( cat == "geq1jets" && bin == "ee"    ){ par[0]=-1.5; par[1]=-1.5; par[2]=1.5; par[3]=2.5; par[4]=2.5; }
        else if( cat == "geq1jets" && bin == "mumu"  ){ par[0]=-1.5; par[1]=-1.5;  par[2]=1.5; par[3]=2.5; par[4]=2.5; }

        double min = 55;
        double max = 0;

        for(int nbin=0;nbin<(ncells-1);nbin++){
           int nbin_content_limit = 0;
           if( bin == "ee" ){ nbin_content_limit=4; } else { nbin_content_limit=1; } 
           if(hf->GetBinContent(nbin)>nbin_content_limit ){ 
               loghf->SetBinContent( nbin, log10(hf->GetBinContent(nbin)*SF ));
               loghf->SetBinError( nbin, log10(hf->GetBinError(nbin)*SF ) );	

               if(hf->GetBinContent(nbin)>3*SF){ max = hf->GetXaxis()->GetBinCenter(nbin); }  //only consider bins with >=1 entries
           }
        }
        TF1 *f1 = new TF1( "fit1", function_used.c_str(), min, max);
        f1->SetParameters( par[0], par[1], par[2], par[3], par[4]);
        /*if( cat == "eq0jets" && bin == "ee" ){ 
	   f1->SetParLimits(  0,  1e-08, 5e-01);
           f1->SetParLimits(  1,     0,     10);
           f1->SetParLimits(  2,  1e-08, 5e-01);
           f1->SetParLimits(  3,     0,     10);
           f1->SetParLimits(  4,     0,     10);
        }
        else if( cat == "eq0jets" && bin == "mumu" ){
           f1->SetParLimits(  0, -2e+04, -2e+02);
           f1->SetParLimits(  1, -1e-02, -1e-08);
           f1->SetParLimits(  2, -1e-05, -1e-15);
           f1->SetParLimits(  3,      0,     10);
           f1->SetParLimits(  4,  2e+02,  2e+04);
        }
        else if( cat == "geq1jets" && bin == "mumu" ){
           f1->SetParLimits(  0, -5e-04, -1e-01);
           f1->SetParLimits(  1,  1e-04,      1);
           f1->SetParLimits(  2,  1e-08,  1e-02);
           f1->SetParLimits(  3,      0,     10);
           f1->SetParLimits(  4,      0,     10);
        }
        else if( cat == "geq1jets" && bin == "ee" ){
           f1->SetParLimits(  0, -5e-03,  -1e-01);
           f1->SetParLimits(  1,  1e-03,       1);
           f1->SetParLimits(  2,  1e-09,   1e-05);
           f1->SetParLimits(  3,      0,      10);
           f1->SetParLimits(  4,      0,      10);
        }*/
        loghf->Fit( f1, "R");
        return f1;   
}

TH1D* SetLogHisto(TH1D *hf){

        TH1D *loghf = new TH1D("","", 1500, 0, 1500);
        int  ncells = hf->GetSize();

        for(int bin=0;bin<(ncells-1);bin++){
           if(hf->GetBinContent(bin)>0){
               loghf->SetBinContent( bin, log10(hf->GetBinContent(bin)));
               loghf->SetBinError( bin, log10(hf->GetBinError(bin)));
           }
        }
    
        return loghf;
}

TGraphErrors* FillHistoWgts(TF1* flep, TF1* fphot, bool cutOnBinContent, TH1D* hlep, TH1D* hphot){
    TGraphErrors *hweights = new TGraphErrors(3000);
    int n=0;
    int I=0;

    double LepMin = flep->GetMinimumX(0, 2000);
    double PhotpMin = hphot->GetMinimumBin();
   
    for(double bin=56; bin<((hphot->GetNbinsX())-1); bin++){
	double mom = hphot->GetBinCenter(bin); 
        double ratio = 0;
        if( hphot->GetBinContent(bin) > 0 ) ratio = pow(10,flep->Eval(std::min(LepMin,mom)))/hphot->GetBinContent(bin);
        if(pow(10,flep->Eval(mom))<0.1 || isnan((float)ratio) || mom>LepMin){
           hweights->SetPoint(I, mom, 0 ); 
   	   hweights->SetPointError(I, 1, 0);
        }else{
           //printf("debug mom=%f lep=%f pot=%f ratio=%f\n", mom, flep->Eval(mom), fphot->Eval(mom), ratio);
           hweights->SetPoint(I, mom, ratio ); 
   	   hweights->SetPointError(I, 1, 0);
        }
        I++;
    }hweights->Set(I);

/*    
    //CHECK Intergral of reweighted photon
    double reweightedIntegral=0;
    for(int bin=1;bin<=hphot->GetNbinsX();bin++){
       double weight=hweights->Eval(hphot->GetBinCenter(bin));
       if(weight>10)continue; //should be less than 1 in principle
       reweightedIntegral += hphot->GetBinContent(bin) * hweights->Eval(hphot->GetBinCenter(bin));
       printf("bin=%i content=%f weigh=%f integral=%f\n", bin, hphot->GetBinContent(bin), hweights->Eval(hphot->GetBinCenter(bin)), reweightedIntegral);
    }
    double scaleFactor = hlep->Integral()/reweightedIntegral;
    printf("ReweightedIntegral=%f vs lep integral=%f --> scale weights by %f  (debug%f)\n", reweightedIntegral, hlep->Integral(), scaleFactor, hphot->Integral());

   //rescale weight graph
   for(int i=0;i<hweights->GetN();i++){  hweights->SetPoint(i,hweights->GetX()[i], hweights->GetY()[i]*scaleFactor);}
*/

    return hweights;
}

TGraphErrors* FillHistWgtsVbf( TGraphErrors* Weights, TH1D* hphot_geq1jet, TH1D* hphot_vbf, TH1D* hqt_geq1jet, TH1D* hqt_vbf){

	TH1D*   geq1jets = (TH1D*)hphot_geq1jet->Clone("photon_geq1jets");
	TH1D*        vbf = (TH1D*)hphot_vbf->Clone("photon_vbf");
	TH1D* qtgeq1jets = (TH1D*)hqt_geq1jet->Clone("qt_geq1jets");
	TH1D*      qtvbf = (TH1D*)hqt_vbf->Clone("qt_vbf");
        TGraphErrors* weights = (TGraphErrors*)Weights->Clone("weights");
        geq1jets->Divide( vbf );
	qtvbf->Divide( qtgeq1jets );
	for(int i=0; i<weights->GetN(); i++){ 
		double qt = weights->GetX()[i];
		for(double bin=55; bin<geq1jets->GetNbinsX(); bin++){
			if( bin == qt ){ 
				double  beta = geq1jets->GetBinContent( bin );
				double alpha = qtvbf->GetBinContent( bin );
				weights->GetY()[i] *= beta*alpha; 
			}
		}
	}

	return weights;	
} 

TGraphErrors* NormalizeWeight( TGraphErrors* Weights, TH1D* hlep, TH1D* hphot){

        TH1D* phot = (TH1D*)hphot->Clone("PhotonRescale");
	TH1D*  lep = (TH1D*)hlep->Clone("LepRescale");
        TGraphErrors* Weights_Nr = (TGraphErrors*)Weights->Clone("TmpGraph");
        int total_bins = hphot->GetSize();
	for( int bin=1; bin<total_bins; bin++){
		phot->SetBinContent( bin, phot->GetBinContent(bin)*Weights->Eval(phot->GetBinCenter(bin)));
	}
	double norm_fact = 0;
	norm_fact = lep->Integral()/phot->Integral();
	for (int i=0;i< Weights_Nr->GetN();i++){ Weights_Nr->GetY()[i] *= norm_fact;}

	return Weights_Nr;	
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
  std::map<string, TGraphErrors*> WeightsFitFunction;
  std::map<string, TGraphErrors*> WeightsFitFunctionNorm;
 
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
     } 

     File->Close();
  }}}//all histos are now loaded


  //Adding ee and mumu histos togheter
  for(unsigned int b=0;b<bin.size();b++){
  for(unsigned int v=0;v<var.size();v++){
     TH1D *hist_mumu = DataHistos["mumu"+bin[b]+var[v]];
     TH1D *hist_ee = DataHistos["ee"+bin[b]+var[v]];
     TH1D *hist = new TH1D( "", "", (hist_mumu->GetSize()-2), 0, 1.25*(hist_mumu->GetBinContent(1)+hist_ee->GetBinContent(1)));

     for(int k=0; k<hist_mumu->GetSize(); k++){
	hist->SetBinContent( k, hist_mumu->GetBinContent(k)+hist_ee->GetBinContent(k));
     }

     if(DataHistos.find("all"+bin[b]+var[v])==DataHistos.end()){
        gROOT->cd(); //make sure that the file is saved in memory and not in file
        DataHistos["all"+bin[b]+var[v]] = (TH1D*)hist->Clone(); //create a new histo, since it's not found
     }else{
        DataHistos["all"+bin[b]+var[v]]->Add(hist); //add to existing histogram
     }
  }}

  cat.push_back("all");
  catColor.push_back(6);
  catL.push_back("all");

  //Compute the fitting function for evey histogram in eq0jets and geq1jets, in ee-mumu-gamma channel. For VBF, we used geq1jets fuctions reweighted in such a way the integral is still ok.
  for(unsigned int c=0;c<cat.size();c++){
  for(unsigned int b=0;b<bin.size();b++){
  for(unsigned int v=0;v<var.size();v++){
     if( var[v]!="_qt" ) continue; //Considering only qt distribution
     if( cat[c] == "ll" ) continue;
     if( bin[b] == "vbf" ){
        FitFunctionMap[cat[c]+bin[b]+var[v]] =  FitFunction(DataHistos[cat[c]+"geq1jets"+var[v]], (DataHistos[cat[c]+"vbf"+var[v]]->Integral()/DataHistos[cat[c]+"geq1jets"+var[v]]->Integral()), bin[b], cat[c]);
     }else if( bin[b] !="vbf" ){
        FitFunctionMap[cat[c]+bin[b]+var[v]] =  FitFunction(DataHistos[cat[c]+bin[b]+var[v]], 1.0, bin[b], cat[c]);
     }    
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
        TH1D* hphot = DataHistos[string("gamma")+bin[b]+var[v]];

        for(unsigned int c=0;c<cat.size();c++){
           if(cat[c]=="gamma" || cat[c] == "ll" )continue;
	    
           if(bfit){
	       if(FitFunctionMap.find(cat[c]+bin[b]+var[v])==FitFunctionMap.end()){printf("Function missing for %s\n", (cat[c]+bin[b]+var[v]).c_str()); continue;}
	       TF1  *flep = FitFunctionMap[cat[c]+bin[b]+var[v]];
               TH1D* hlep = DataHistos[cat[c]+bin[b]+var[v]];
	       if( bin[b] !="vbf" ){ WeightsFitFunction[ cat[c]+bin[b]+var[v]+"weight"]= FillHistoWgts(flep,fphot, false, hlep, hphot);}
	        
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

   //The weights of VBF are computed startig from the weights geq1jets
   for(unsigned int c=0;c<cat.size();c++){
	if( cat[c] == "gamma" || cat[c] == "ll" ) continue;
   	WeightsFitFunction[ cat[c]+"vbf_qtweight"] = FillHistWgtsVbf( WeightsFitFunction[ cat[c]+"geq1jets_qtweight"], DataHistos["gammageq1jets_qt"], DataHistos["gammavbf_qt"], DataHistos[cat[c]+"geq1jets_qt"], DataHistos[cat[c]+"vbf_qt"]);
   }

  //Normalize the weights to the integral.
  for(unsigned int v=0;v<var.size();v++){
     for(unsigned int b=0;b<bin.size();b++){
        for(unsigned int c=0;c<cat.size();c++){
		if( var[v] != "_qt" ) continue;
		if( cat[c] == "gamma" || cat[c] == "ll" ) continue;
		TH1D*  hlep = DataHistos[string("gamma")+bin[b]+var[v]];
		TH1D* hphot = DataHistos[cat[c]+bin[b]+var[v]];
		
		WeightsFitFunctionNorm[cat[c]+bin[b]+var[v]+"weight_nr"] = NormalizeWeight( WeightsFitFunction[ cat[c]+bin[b]+var[v]+"weight"], hlep, hphot); 
  }}}

  var.push_back("_qtweight_nr");
  

  //make the plots
  if(!bfit){
  for(unsigned int v=0;v<var.size();v++){
     double xmin,xmax;
     double ymin=0.5, ymax=1E6;
     if (var[v]=="_qt" || var[v]=="_qtrebin") {                    xmin=55.00; xmax=1000.00;
     } else if (var[v]=="_qtweight" || var[v]=="_qtrebinweight" || var[v]=="_qtweight_nr") { xmin=55.00; xmax=1000.00; ymin=0.0001;  ymax=0.5;
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
	
       if( var[v]=="_qtrebin" || var[v]=="_qtrebinweight" ) continue;
       for(unsigned int c=0;c<cat.size();c++){
	  string canvas_name;
          TCanvas* c1 = new TCanvas("c1","c1",500*bin.size(), 500);
          c1->Divide(bin.size(),2);
 
          for(unsigned int b=0;b<bin.size();b++){
	     c1->cd(1+b);//->SetLogy(true);
            
             TLegend* leg = new TLegend(0.52,0.67,0.92,0.90);
             leg->SetFillStyle(0);
             leg->SetBorderSize(0);
             leg->SetTextSize(0.04);
             leg->SetTextFont(42);
             leg->SetHeader(binL[b].c_str());
	     TH1D* hist = NULL;

             if(var[v]=="_qt"){
          	if( cat[c] == "ll" || cat[c] == "gamma" ) continue;
                
		TF1 *f1 = NULL;
		hist = SetLogHisto(DataHistos[cat[c]+bin[b]+var[v]]); 
//		hist = DataHistos[cat[c]+bin[b]+var[v]]; 
		f1 = FitFunctionMap[cat[c]+bin[b]+var[v]];
//		f1->SetParameter(0, rebin*f1->GetParameter(0));
                f1->SetLineColor(kRed);
                f1->SetLineWidth(2);
	        hist->SetStats(false);
//		hist->Rebin(rebin);
        	hist->GetXaxis()->SetTitle("q_{T} [GeV]");
		hist->GetXaxis()->SetRangeUser(55,500); 
		hist->GetXaxis()->SetTitleOffset(0.8);
        	hist->GetYaxis()->SetTitle("Events");
		hist->GetYaxis()->SetTitleOffset(0.8); 
 		hist->Draw("HIST P same");
		f1->Draw("same");              
		canvas_name = "/FittingFunction_"+cat[c]+var[v];

                c1->cd(1+b+bin.size());//->SetLogy(true);
                f1->Draw();

             } else if(var[v]=="_qtweight_nr" ){
                if( cat[c] == "ll" || cat[c] == "gamma" || cat[c] == "all" ) continue;
		
                TGraphErrors* graph = NULL;
                graph = WeightsFitFunctionNorm[cat[c]+bin[b]+"_qtweight_nr"];
                graph->SetMarkerColor(kRed);
                graph->GetXaxis()->SetRangeUser(  0., 1000);
                graph->GetXaxis()->SetTitleSize(.055);
                graph->GetYaxis()->SetTitleSize(.055);
                graph->GetXaxis()->SetLabelSize(.05);
                graph->GetYaxis()->SetLabelSize(.05);
                graph->GetXaxis()->SetTitle("Z Boson pT [GeV]");
                graph->GetYaxis()->SetTitle("Weights");
                graph->Draw();
                leg->AddEntry(graph, catL[c].c_str(), "LP");
                canvas_name = "/Weights_"+cat[c]+var[v];

	     } else {
                if( cat[c] == "all" || cat[c] == "ll" || var[v] == "_qtweight") continue;
                
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
          c1->SaveAs((outDir + canvas_name +".root").c_str());


       }
    }
  }

  //SAVE THE material to the weight file
  TFile* OutputFile = new TFile(outFile.c_str(),"RECREATE");
  OutputFile->cd();
  for(unsigned int c=0;c<cat.size();c++){
     if( cat[c] == "ll" || cat[c] == "gamma" || cat[c] == "all" )continue;
     for(unsigned int b=0;b<bin.size();b++){
     for(unsigned int v=0;v<var.size();v++){
             if(!bfit){ if(DataHistos.find(cat[c]+bin[b]+var[v])==DataHistos.end()){printf("Histo missing for %s\n", (cat[c]+bin[b]+var[v]).c_str()); continue;} }

             TH1D* hist = DataHistos[cat[c]+bin[b]+var[v]];
             if(bfit){
                 if(var[v]=="_qtweight_nr"){ 		
			   
			   TGraphErrors* graph = WeightsFitFunctionNorm[cat[c]+bin[b]+"_qtweight_nr"];
			   graph->Write((cat[c]+bin[b]+"_qt_datafitfunctionwgts_norm").c_str()); 
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
                          hist->Write((cat[c]+"vbf"+"_qmass").c_str());  
                     }
                     hist->Write((cat[c]+bin[b]+"_qmass").c_str());
             }
        }
     }
  }
  OutputFile->Write(); OutputFile->Close();
}

