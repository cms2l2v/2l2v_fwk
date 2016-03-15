
#include <string>
#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TF1.h"
#include "TCutG.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TMarker.h"
#include "TColor.h"

#include "UserCode/llvv_fwk/interface/tdrstyle.h"
#include "UserCode/llvv_fwk/src/tdrstyle.C"
#include "UserCode/llvv_fwk/interface/RootUtils.h"

#include "UserCode/llvv_fwk/interface/HxswgUtils.h"
#include "UserCode/llvv_fwk/src/HxswgUtils.cc"


TGraph2D* get2HDMLimitsTgVsAmB(string FilePath, TGraph2D* g2dMassVsWidth, TGraph* SMhiggsWidth, double OnlyCosAmB, int limitType=4){
   if(FilePath=="")FilePath="Weight_2HDM_Model.txt";
   FILE* pFile = fopen(FilePath.c_str(),"r");
   if(!pFile){printf("Can't open %s\n",FilePath.c_str()); exit(0);}

   TGraph2D* graph = new TGraph2D(9999);    int N=0;
   char line[4096];
   double tgB, cosAmB=OnlyCosAmB, mH2, width, BRtoZZ, XSec, unused;
   while(fgets(line, 4096, pFile)){
      if(line[0]=='\\' && line[1]=='\\')continue;
//      sscanf(line,"tgBeta=%lf cosAmB=%lf MH2=%lf Width=%lf BrH2toZZ=%lf XSec=%lf\n",&tgB, &cosAmB, &mH2, &width, &BRtoZZ, &XSec);
         // tanBeta mH  xs_bbH   width_H   mH    xs_ggH  br_HW+H- br_Htautau br_Hss   br_Hgg   br_Hbb  br_HZgamma br_Htt br_HZA   br_Hmumu br_Hgammagamma br_Hee   br_HWW    br_HZZ  br_HZh   br_Hhh   br_Hcc
       sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &tgB, &mH2, &unused, &width, &unused, &XSec, &unused, &unused,   &unused, &unused, &unused, &unused, &unused, &unused, &unused, &unused,      &unused, &unused, &BRtoZZ, &unused, &unused, &unused);     
       if(mH2<200 || int(mH2)%20!=0)continue;
       //printf("line=%s\n--> log10(tan)=%f mh=%f w=%f xsec=%f br=%f\n", line, log10(tgB), mH2, width, XSec, BRtoZZ);
//      if(fabs(cosAmB-OnlyCosAmB)<0.05){

        double widthEquivalent = std::min(std::max(width/SMhiggsWidth->Eval(mH2), 0.05), 1.0);;  //take Width=0.01 if < 0.01
      
        if      (limitType==0){ graph->SetPoint(N, mH2, tgB, 1000*XSec);  N++;
        }else if(limitType==1){ graph->SetPoint(N, mH2, tgB,    BRtoZZ);  N++;          
        }else if(limitType==2){ graph->SetPoint(N, mH2, tgB, 1000*XSec*BRtoZZ);  N++;          
        }else if(limitType==3){ graph->SetPoint(N, mH2, tgB, width/SMhiggsWidth->Eval(mH2));  N++;          
        }else if(limitType==4){ graph->SetPoint(N, mH2, tgB, g2dMassVsWidth->Interpolate(mH2, widthEquivalent)  );  N++;          
//             printf("Signal Strength limit = %f = %f/%f\n",  g2dMassVsWidth->Interpolate(mH2, width/SMhiggsWidth->Eval(mH2)) / (XSec*BRtoZZ*1000),  g2dMassVsWidth->Interpolate(mH2, width/SMhiggsWidth->Eval(mH2)) , (XSec*BRtoZZ*1000) );
         }
//      }
//if(N>1500)break;
   }fclose(pFile);
   graph->Set(N);
   return graph;
}
  
void scaleGraph(TGraph* Limit, double scale){
   for(int i=0;i<Limit->GetN();i++){
      Limit->SetPoint(i, Limit->GetX()[i], Limit->GetY()[i]*scale);
   }   
}

void drawPointOnGrid(TGraph2D* graph){
        TH2D* h2d = graph->GetHistogram();              
	double* Xs = graph->GetX();
	double* Ys = graph->GetY();
	for(int p=0;p<graph->GetN();p++){
           if(Xs[p]<h2d->GetXaxis()->GetXmin())continue;
           if(Xs[p]>h2d->GetXaxis()->GetXmax())continue;
           if(Ys[p]<h2d->GetYaxis()->GetXmin())continue;
           if(Ys[p]>h2d->GetYaxis()->GetXmax())continue;
           TMarker* mark = new TMarker(Xs[p], Ys[p], 20);
	   mark->SetMarkerSize(0.3);
	   mark->Draw("same");
           //printf("draw a point at %f - %f\n", Xs[p], Ys[p]);
	}
}

void Histo2DToExp(TH2D* histo){
              for(int x=1;x<=histo->GetNbinsX();x++){
                 for(int y=1;y<=histo->GetNbinsY();y++){
                    histo->SetBinContent(x,y,pow(10,histo->GetBinContent(x,y)));
              }}
}

void HideUnexplored(TH2D* histo, TH2D* width){
              for(int x=1;x<=histo->GetNbinsX();x++){
                 for(int y=1;y<=histo->GetNbinsY();y++){
                    if(width->GetBinContent(x,y)>1) histo->SetBinContent(x,y,0);
              }}
}


TGraph* getGraph(string name, int color, int width, int style, TLegend* LEG, TGraph* Ref, int type, string filePath){
//   filePath+="/LimitSummary";
   FILE* pFile = fopen(filePath.c_str(),"r");
   if(!pFile){printf("Can't open %s\n",filePath.c_str()); exit(0);}
   double mass=0, th, exp=0, obs=0, unused=0;// char buffer[1024];
   double explow2=0, explow1=0, expup1=0, expup2=0;

   TGraph* graph = new TGraph(100);
   int N=0;
   while(fscanf(pFile,"$%le$ & $%le$ & $[%le,%le]$ & $[%le,%le]$ & $%le$ & Th=$%le$ & pValue=$%le$\\\\\\hline\n",&mass, &exp, &explow1, &expup1, &explow2,&expup2,&obs, &th, &unused) != EOF){

//      printf("%i %f - %f - %f\n",N,mass,exp, th);

      double value = exp;
      if(abs(type)==0) value = th;
      else if(abs(type)==1) value = exp;
      else if(abs(type)==2) value = obs;
      else if(abs(type)==3) value = explow1;
      else if(abs(type)==4) value = expup1;
      else if(abs(type)==5) value = explow2;
      else if(abs(type)==6) value = expup2;
      else value = obs;


//      if(type<0 && filePath.find("cp0.80")!=string::npos) value *= pow(0.8, 2);
//      if(type<0 && filePath.find("cp0.60")!=string::npos) value *= pow(0.6, 2);
//      if(type<0 && filePath.find("cp0.30")!=string::npos) value *= pow(0.3, 2);
//      if(type<0 && filePath.find("cp0.10")!=string::npos) value *= pow(0.1, 2);

      if(Ref){value/=Ref->Eval(mass);}
      graph->SetPoint(N, mass, value);N++;
   }
   fclose(pFile);
   graph->Set(N);

   graph->SetName(name.c_str());
   graph->SetLineColor(color);
   graph->SetLineWidth(width);
   graph->SetLineStyle(style);
   if(LEG)LEG->AddEntry(graph, name.c_str()      ,"L");
   return graph;
}


TGraph** getGraphs(string name, int color, int width, TLegend* LEG, TGraph* Ref, string filePath){
   TGraph** graphs = new TGraph*[7];
   for(int i=0;i<7;i++){  
      char nameBuf[255];sprintf(nameBuf,"%s_%i",name.c_str(),i);  
      graphs[i] = getGraph(nameBuf, color, width, i<=1?2:1, NULL, Ref, i, filePath);
      graphs[i]->SetTitle(name.c_str());
   }
   if(LEG)LEG->AddEntry(graphs[2], name.c_str()      ,"L");
   return graphs;
}


TCutG* GetErrorBand(string name, int N, double* Mass, double* Low, double* High, double ymin=0.0, double ymax=1.0)
{
   TCutG* cutg = new TCutG(name.c_str(),2*N+2);
   cutg->SetFillColor(kGreen-7);
   cutg->SetLineStyle(0);
   cutg->SetLineColor(0);
   int I = 0;
   for(int i=0;i<N;i++){cutg->SetPoint(I,Mass[i]    , std::max(ymin,std::min(ymax,Low[i]     )));I++; }
                        cutg->SetPoint(I,Mass[N-1]  , std::max(ymin,std::min(ymax,Low[N-1]   )));I++;
                        cutg->SetPoint(I,Mass[N-1]  , std::max(ymin,std::min(ymax,High[N-1]  )));I++;
   for(int i=0;i<N;i++){cutg->SetPoint(I,Mass[N-1-i], std::max(ymin,std::min(ymax,High[N-1-i])));I++;}
   return cutg;
}


TGraph* getContour(TH2D* inputHisto, TCanvas* goodCanvas, int Width, int Style, int Color){
	TCanvas* c1 = new TCanvas("temp", "temp",600,600);

	TH2D* histo = (TH2D*)inputHisto->Clone("temp");
	for(int x=0;x<=histo->GetNbinsX();x++){
		for(int y=0;y<=histo->GetNbinsY()+1;y++){
			if(histo->GetBinContent(x,y)<=0 || histo->GetBinContent(x,y)>=1E10)histo->SetBinContent(x,y, 1E9);
		}
	}

	double levels[] = {0.9, 1.0, 1.1};
	histo->SetContour(3, levels);
	histo->Draw("CONT LIST");
	c1->Update();
	TObjArray* contours = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
	Int_t ncontours     = contours->GetSize();
	TList *list         = (TList*)contours->At(1);
	delete c1;

	goodCanvas->cd();
	printf("list size = %i\n", (int)list->GetSize());
	if(list->GetSize()<=0)return new TGraph(0);

	TGraph* EXCLUSION = new TGraph(0);
	for(int i=0;i<list->GetSize();i++){               
		EXCLUSION   = (TGraph*)(list->At(i)->Clone("copy"));
		EXCLUSION->SetLineColor(Color);
		EXCLUSION->SetLineWidth(Width);
		EXCLUSION->SetLineStyle(Style);


                printf("EXCLUSION->GetMean = %f\n", EXCLUSION->GetMean(2));
                if(EXCLUSION->GetMean(2)<=0.15)continue; //BUGFIX

		EXCLUSION->Draw("CL same");
	}

	return EXCLUSION;
}




int mapIndex(double CP, double BR){return CP<=0?-1:(int)(CP*100+BR*10);}

void finalLimitPlot(){
   setTDRStyle();  
   gStyle->SetPadTopMargin   (0.06);
   gStyle->SetPadBottomMargin(0.12);
   gStyle->SetPadRightMargin (0.16);
   gStyle->SetPadLeftMargin  (0.14);
   gStyle->SetTitleSize(0.04, "XYZ");
   gStyle->SetTitleXOffset(1.1);
   gStyle->SetTitleYOffset(1.45);
//   gStyle->SetPalette(1);
   gStyle->SetNdivisions(505);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);



     //inverted deep blue color palette
//      Double_t stops[9] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000};
      Double_t stops[2] = { 0.0000, 1.0000};
      Double_t red[2]   = {180./255.,0./255.};
      Double_t green[2] = {245./255.,0./255.};
      Double_t blue[2]  = {245./255.,255./255.};
      TColor::CreateGradientColorTable(2, stops, red, green, blue, 255);

//   gStyle->SetPalette(51);




   TGraph* higgsWidth = Hxswg::utils::getHWidthExtended();


   TCanvas* c1;
   TH1F* framework;
   TH2F* framework2d;
   TLegend* LEG, *LEGTH;
   TGraph* Ref;

   std::map<int, int> colorMap;
   colorMap[mapIndex(1.0, 0.0)] = 1;
   colorMap[mapIndex(0.9, 0.0)] = 1;
   colorMap[mapIndex(0.8, 0.0)] = 6;
   colorMap[mapIndex(0.7, 0.0)] = 6;
   colorMap[mapIndex(0.6, 0.0)] = 4;
   colorMap[mapIndex(0.5, 0.0)] = 4;
   colorMap[mapIndex(0.4, 0.0)] = 3;
   colorMap[mapIndex(0.3, 0.0)] = 3;
   colorMap[mapIndex(0.2, 0.0)] = 2;
   colorMap[mapIndex(0.1, 0.0)] = 2;

   std::vector<double> CPs = {1.0,0.6,0.3,0.1,0.0};
   std::vector<double> BRs = {0.0};

  //LIMIT ON SIGNAL STRENGTH
   string Directories2[]={"cards_SB13TeV", "cards_SB13TeV_GGF", "cards_SB13TeV_VBF"};  //DEBUG
   for(unsigned int D=0;D<sizeof(Directories2)/sizeof(string);D++){
      string Dir     = Directories2[D];
      char tmp[1024];sprintf(tmp, "%s_cp%4.2f_brn%4.2f/", Dir.c_str(), 1.0, 0.0);
      string SaveDir = tmp;

      string prod = "pp";
      if(Dir.find("GGF")!=std::string::npos)prod="gg";
      if(Dir.find("VBF")!=std::string::npos)prod="qq";
      bool strengthLimit = false;
      if(prod=="pp")strengthLimit=true;

  
      std::map<int, TGraph**> gCPBR;
      for(unsigned int CPi = 0; CPi<CPs.size();CPi++){
      for(unsigned int BRi = 0; BRi<BRs.size();BRi++){
         double CP = CPs[CPi];  double BR = BRs[BRi]; 
         if ( CP*CP / (1-BR)>1.0 )continue; //skip point leading to a width larger than SM
         char limitpath[1024]; sprintf(limitpath, "%s_cp%4.2f_brn%4.2f/Stength_LimitSummary", Dir.c_str(), CP<=0.1?0.1:CP, BR);
//         char limitlegend[1024]; sprintf(limitlegend, "C'=%3.1f BRnew=%3.1f", CP, BR);
         char limitlegend[1024]; sprintf(limitlegend, "C'=%3.1f", CP);
         std::cout<<"LIMITPATH="<<limitpath<< "\n";         
         gCPBR[ mapIndex(CP, BR) ] = getGraphs(limitlegend, colorMap[mapIndex(CP, BR)], 2, NULL  , NULL, limitpath);

         //overwrite the theory graph with a higher granularity one
         TGraph* THXSec   = Hxswg::utils::getXSec(Dir); 
         scaleGraph(THXSec, 1000);  //convert cross-section to fb
         double XSecScaleFactor = pow(CP,2) * (1-BR);
         scaleGraph(THXSec, XSecScaleFactor);
         gCPBR[ mapIndex(CP, BR) ][0] = THXSec;
      }}

      // build 2D graphs of Mass vs C' limits
      TGraph2D* g2dMassVsCp[7];    
      for(unsigned int i=0;i<7;i++){
         g2dMassVsCp[i]    = new TGraph2D( 9999 ); 
         int Ig2d=0;
         for(unsigned int CPi = 0; CPi<CPs.size();CPi++){
            double CP = CPs[CPi];  double BR = 0.0;             
            TGraph* graph = gCPBR[ mapIndex(CP, BR) ][i];          
            for(int p=0;p<graph->GetN ();p++){
               double mass, limit;
               graph->GetPoint(p, mass, limit);
               g2dMassVsCp[i]   ->SetPoint(Ig2d, mass, CP   , limit);
               Ig2d++;              
            }
         }
         g2dMassVsCp   [i]->Set(Ig2d);
//
//         char buffer[255]; sprintf(buffer, "%s%i", "gridMvsC",i);
//         TH2D* histGrid  = new TH2D(buffer  , buffer  , 26, 200, 1500, 20, 0.1, 1.0);
//         histGrid->SetDirectory(0);
//         g2dMassVsCp   [i]->SetHistogram(histGrid);
      } 

      TGraph2D* g2dMassVsWidth[7];  
      std::vector<  double> WidthVect = {0.000001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

      for(unsigned int i=0;i<7;i++){
         g2dMassVsWidth[i] = new TGraph2D( 9999 );        
         int Ig2d=0;

         for(unsigned int W=0;W<WidthVect.size();W++){
         for(double Mass = 200; Mass<=1500;Mass+=25){
               double Width = WidthVect[W];
               double limit= g2dMassVsCp[i]->Interpolate(Mass, sqrt(Width));
               g2dMassVsWidth[i]->SetPoint(Ig2d, Mass, Width, limit);
               Ig2d++;              
         }}
         g2dMassVsWidth[i]->Set(Ig2d);        
      } 


      ///////////////////////////////////////////////
      //limit Versus mass for different C' values
      ///////////////////////////////////////////////

         c1 = new TCanvas("c", "c",600,600);
         c1->SetLogy(true);
         framework = new TH1F("Graph","Graph",1,190,1510);
         framework->SetStats(false);
         framework->SetTitle("");
         framework->GetXaxis()->SetTitle("M_{H} [GeV]");
         if(strengthLimit){
            framework->GetYaxis()->SetTitle("#mu = #sigma_{95%} / #sigma_{th}");
            framework->GetYaxis()->SetRangeUser(1E-2,1E3);
         }else{
            framework->GetYaxis()->SetTitle((string("#sigma_{95%} (") + prod +" #rightarrow H #rightarrow ZZ) (fb)").c_str());
            framework->GetYaxis()->SetRangeUser(1E1,1E5);
         }
         framework->GetYaxis()->SetTitleOffset(1.40);
         framework->Draw();

         if(strengthLimit){
            TLine* SMLine = new TLine(framework->GetXaxis()->GetXmin(),1.0,framework->GetXaxis()->GetXmax(),1.0);
            SMLine->SetLineWidth(2); SMLine->SetLineStyle(1); SMLine->SetLineColor(4);      
            SMLine->Draw("same C");
         }


         LEG = new TLegend(0.70,0.70,0.85,0.94);
         LEG->SetTextFont(43); LEG->SetTextSize(15);   LEG->SetTextAlign(12);
         LEG->SetHeader("Observed");
         LEG->SetFillStyle(0);
         LEG->SetBorderSize(0);

         TLegend* LEGExp = new TLegend(0.55,0.70,0.70,0.94);
         LEGExp->SetTextFont(43); LEGExp->SetTextSize(15);   LEGExp->SetTextAlign(12);         
         LEGExp->SetHeader("Expected");
         LEGExp->SetFillStyle(0);
         LEGExp->SetBorderSize(0);

         TLegend* LEGTh = NULL;
         if(!strengthLimit){
            LEGTh = new TLegend(0.40,0.70,0.55,0.94);
            LEGTh->SetTextFont(43); LEGTh->SetTextSize(15);   LEGTh->SetTextAlign(12);         
            LEGTh->SetHeader("Predicted");
            LEGTh->SetFillStyle(0);
            LEGTh->SetBorderSize(0);
         }
        
        
         for(unsigned int CPi = 0; CPi<CPs.size();CPi++){
         for(unsigned int BRi = 0; BRi<BRs.size();BRi++){
            double CP = CPs[CPi];  double BR = BRs[BRi]; 
            if(CP<0.1)continue;


            if(BR!=0.0)continue;
            if(!gCPBR[ mapIndex(CP, BR) ])continue;
            TGraph* g = (gCPBR[ mapIndex(CP, BR) ])[1+1];
            if(!g) continue;
            g->SetLineStyle(1);
            g->Draw("C same");
            LEG->AddEntry(g, g->GetTitle(), "L");            


            TGraph* gExp = (gCPBR[ mapIndex(CP, BR) ])[1];
            if(!gExp) continue;
            gExp->SetLineStyle(2);
            gExp->Draw("C same");
            LEGExp->AddEntry(gExp, g->GetTitle(), "L");            


           if(!strengthLimit){
              TGraph* THXSec = (gCPBR[ mapIndex(CP, BR) ])[0];
                THXSec->SetLineStyle(3);  THXSec->SetLineWidth(1);  THXSec->SetLineColor(g->GetLineColor());
                THXSec->Draw("same C");
                LEGTh->AddEntry(THXSec, g->GetTitle(), "L");
           }
         }}
         LEG  ->Draw("same");
         LEGExp  ->Draw("same");
         if(LEGTh)LEGTh  ->Draw("same");


         utils::root::DrawPreliminary(2268.759, 13, gPad->GetLeftMargin(),gPad->GetBottomMargin(),gPad->GetRightMargin(),gPad->GetTopMargin()+0.03);
//         if(observed==0){

            gPad->RedrawAxis();
            c1->SaveAs((SaveDir+"/Stength_FinalPlot.png").c_str());
            c1->SaveAs((SaveDir+"/Stength_FinalPlot.pdf").c_str());
            c1->SaveAs((SaveDir+"/Stength_FinalPlot.C"  ).c_str());
//         }else{
//            c1->SaveAs((SaveDir+"/Stength_FinalPlot_Obs.png").c_str());
//            c1->SaveAs((SaveDir+"/Stength_FinalPlot_Obs.pdf").c_str());
//            c1->SaveAs((SaveDir+"/Stength_FinalPlot_Obs.C"  ).c_str());
//         }

      ///////////////////////////////////////////////
     //Mass Versus Cprime limits  2D
     ///////////////////////////////////////////////
     for(int mode=0; mode<=1; mode++){
     //mode=0 --> Mass versus C'     interpolated from 2D limit(CPrime, BRNew)
     //mode=1 --> Mass versus Gamma  interpolated from 2D limit(CPrime, BRNew)
     
        TH2D* grid   = new TH2D("grid"  , "grid"  , 300, 0, 1500, 500, 0, 1);
        TGraph2D*   graph[] = {g2dMassVsCp   [1], g2dMassVsCp   [2], g2dMassVsCp   [3], g2dMassVsCp   [4]};
        if(mode==1){ graph[0] = g2dMassVsWidth[1]; graph[1] = g2dMassVsWidth[2]; graph[2] = g2dMassVsWidth[3]; graph[3] = g2dMassVsWidth[4];}
        graph[0]->SetHistogram((TH2D*)grid->Clone("GRIDg2dMassVsCpA"));
        graph[1]->SetHistogram((TH2D*)grid->Clone("GRIDg2dMassVsCpB"));
        graph[2]->SetHistogram((TH2D*)grid->Clone("GRIDg2dMassVsCpC"));
        graph[3]->SetHistogram((TH2D*)grid->Clone("GRIDg2dMassVsCpD"));


        //build th cross-section 2D plane
        TGraph* THXSec   = Hxswg::utils::getXSec(Dir); 
        scaleGraph(THXSec, 1000);  //convert cross-section to fb
        TH2D* THXSec2D   = (TH2D*)grid->Clone("ThXSec");
        for(int x=1;x<=THXSec2D->GetNbinsX();x++){
           double mass = THXSec2D->GetXaxis()->GetBinCenter(x);
           for(int y=1;y<=THXSec2D->GetNbinsY();y++){
              double CP = THXSec2D->GetYaxis()->GetBinCenter(y);  double BR=0.0;
              if(mode==1){CP = sqrt(CP);}  //go from width=C'² --> to C'
              double XSecScaleFactor = pow(CP,2) * (1-BR);
              THXSec2D->SetBinContent(x,y,THXSec->Eval(mass)*XSecScaleFactor);
        }}
        TH2D* signalStrength[] = { (TH2D*) (graph[0]->GetHistogram()->Clone("signalStrengthExp")),  (TH2D*) (graph[1]->GetHistogram()->Clone("signalStrengthObs")),  (TH2D*) (graph[2]->GetHistogram()->Clone("signalStrengthExpM1")), (TH2D*) (graph[3]->GetHistogram()->Clone("signalStrengthExpP1")) };
        if(!strengthLimit){
           signalStrength[0]->Divide(THXSec2D);
           signalStrength[1]->Divide(THXSec2D);
           signalStrength[2]->Divide(THXSec2D);
           signalStrength[3]->Divide(THXSec2D);
        }

           
       for(int observed=0;observed<=1;observed++){
           TH2D* h2d = graph[observed]->GetHistogram();
           

           c1 = new TCanvas("c", "c",600,600);
           c1->SetLogz(true);      
           c1->SetRightMargin(0.17);
           framework2d = new TH2F("Graph","Graph",1,250, 1500, 1,0,1.0);
           framework2d->SetStats(false);
           framework2d->SetTitle("");
           framework2d->GetXaxis()->SetTitle("M_{H} [GeV]");
           if(mode==0)framework2d->GetYaxis()->SetTitle("C'");
           if(mode==1)framework2d->GetYaxis()->SetTitle("#Gamma/#Gamma_{SM}");          
           framework2d->GetYaxis()->SetTitleOffset(1.40);
           framework2d->Draw("");
           h2d->GetZaxis()->SetTitleOffset(1.60);
           if(strengthLimit){
              h2d->GetZaxis()->SetTitle((string(observed==0?"Expected":"Observed") + " #mu = #sigma_{95%} / #sigma_{th}").c_str());
              h2d->GetZaxis()->SetRangeUser(1E-1,1E2);
           }else{
              h2d->GetZaxis()->SetTitle((string(observed==0?"Expected":"Observed") + string(" #sigma_{95%} (") + prod +" #rightarrow H #rightarrow ZZ) (fb)").c_str());
              h2d->GetZaxis()->SetRangeUser(1E1,1E3);
           }
           h2d->Draw("COLZ same");

           //strength contour
           getContour(signalStrength[0], c1, 2, 7, 1); //exp
           getContour(signalStrength[2], c1, 1, 2, 1); //expM1
           getContour(signalStrength[3], c1, 1, 2, 1); //expP1
           if(observed==1) getContour(signalStrength[observed], c1, 3, 1, 1); //obs
          
           if(strengthLimit){
              TPaveText* pave1 = NULL;
              TGraph* indirectLimit = new TGraph(2);
              indirectLimit->SetLineColor(17);
              indirectLimit->SetFillColor(17);
              indirectLimit->SetLineWidth(502);
              indirectLimit->SetFillStyle(3005);

              //µ_h125 = 1.00 +- 0.14 --> C'²<0.28 at 95% C.L.
              //µ_h125 = 1.09 +- 0.1  --> C'²<13 at 95% C.L.
              if(mode==0){
                 indirectLimit->SetPoint(0, framework2d->GetXaxis()->GetXmin(), sqrt(0.13) );
                 indirectLimit->SetPoint(1, framework2d->GetXaxis()->GetXmax(), sqrt(0.13) );

                  pave1 = new TPaveText(framework2d->GetXaxis()->GetXmax()*0.7, sqrt(0.13)+0.10, framework2d->GetXaxis()->GetXmax(),sqrt(0.13));
              }else{
                 indirectLimit->SetPoint(0, framework2d->GetXaxis()->GetXmin(), 0.13 );
                 indirectLimit->SetPoint(1, framework2d->GetXaxis()->GetXmax(), 0.13 );

                 pave1 = new TPaveText(framework2d->GetXaxis()->GetXmax()*0.7, 0.13+0.10, framework2d->GetXaxis()->GetXmax(),0.13);
              }
              indirectLimit->Draw("same");

              pave1->SetBorderSize(0);
              pave1->SetFillStyle(0);
              pave1->SetTextAlign(11);
              pave1->SetTextFont(62);
              pave1->SetTextSize(0.03);
              pave1->SetTextColor(17);
              pave1->AddText("#mu_{h(125)} = 1.09 #pm 0.11");
              pave1->Draw("same");
           }


            utils::root::DrawPreliminary(2268.759, 13, gPad->GetLeftMargin(),gPad->GetBottomMargin(),gPad->GetRightMargin(),gPad->GetTopMargin()+0.03);

            char massStr[512];
            if(mode==0)sprintf(massStr, "MassVsCp%s", observed==0?"":"_Obs");
            if(mode==1)sprintf(massStr, "MassVsWidth%s", observed==0?"":"_Obs");          

            gPad->RedrawAxis();
            c1->SaveAs((SaveDir+"/Stength_FinalPlot2D_"+massStr+".png").c_str());
            c1->SaveAs((SaveDir+"/Stength_FinalPlot2D_"+massStr+".pdf").c_str());
            c1->SaveAs((SaveDir+"/Stength_FinalPlot2D_"+massStr+".C").c_str());
         }//end of observed loop

      }//end of mode loop


      ///////////////////////////////////////////////
     //Limits on 2HDM
     ///////////////////////////////////////////////

      if(prod=="gg"){
         TH2D* grid   = new TH2D("grid2HDM"  , "grid"  , 650, 200, 1500, 1000, 0.5, 60);         
         for(int type=1;type<=2;type++){
             for(double cosBmA=0.1; cosBmA<=0.1; cosBmA+=0.1){

             TGraph2D* graph[] = {
                get2HDMLimitsTgVsAmB(type==1?"2HDM_type1_H.txt":"2HDM_type2_H.txt", NULL                      , higgsWidth, cosBmA, 0),  //0 --> sigma(gg-->H2)
                get2HDMLimitsTgVsAmB(type==1?"2HDM_type1_H.txt":"2HDM_type2_H.txt", NULL                      , higgsWidth, cosBmA, 1),  //1 --> BR(H2 -->ZZ) 
                get2HDMLimitsTgVsAmB(type==1?"2HDM_type1_H.txt":"2HDM_type2_H.txt", NULL                      , higgsWidth, cosBmA, 2),  //2 --> sigma(gg-->H2-->ZZ)
                get2HDMLimitsTgVsAmB(type==1?"2HDM_type1_H.txt":"2HDM_type2_H.txt", NULL                      , higgsWidth, cosBmA, 3),  //3 --> Gamma/Gamma_SM
                get2HDMLimitsTgVsAmB(type==1?"2HDM_type1_H.txt":"2HDM_type2_H.txt", g2dMassVsWidth[1+0]       , higgsWidth, cosBmA, 4),  //4 --> Exp Limit on sigma(gg-->H2-->ZZ)
                get2HDMLimitsTgVsAmB(type==1?"2HDM_type1_H.txt":"2HDM_type2_H.txt", g2dMassVsWidth[1+1]       , higgsWidth, cosBmA, 4),  //5 --> Obs Limit on sigma(gg-->H2-->ZZ)
                get2HDMLimitsTgVsAmB(type==1?"2HDM_type1_H.txt":"2HDM_type2_H.txt", g2dMassVsWidth[1+2]       , higgsWidth, cosBmA, 4),  //6 --> Exp M1 Limit on sigma(gg-->H2-->ZZ)
                get2HDMLimitsTgVsAmB(type==1?"2HDM_type1_H.txt":"2HDM_type2_H.txt", g2dMassVsWidth[1+3]       , higgsWidth, cosBmA, 4)   //7 --> Exp P1 Limit on sigma(gg-->H2-->ZZ)
             };

             graph[0]->SetHistogram((TH2D*)grid->Clone("GRID2HDMmH0"));
             graph[1]->SetHistogram((TH2D*)grid->Clone("GRID2HDMmH0"));
             graph[2]->SetHistogram((TH2D*)grid->Clone("GRID2HDMmH0"));
             graph[3]->SetHistogram((TH2D*)grid->Clone("GRID2HDMmH0"));
             graph[4]->SetHistogram((TH2D*)grid->Clone("GRID2HDMmH0"));
             graph[5]->SetHistogram((TH2D*)grid->Clone("GRID2HDMmH0"));
             graph[6]->SetHistogram((TH2D*)grid->Clone("GRID2HDMmH0"));
             graph[7]->SetHistogram((TH2D*)grid->Clone("GRID2HDMmH0"));

             TH2D* h2ds[] = {
                graph[0]->GetHistogram(),
                graph[1]->GetHistogram(),
                graph[2]->GetHistogram(),
                graph[3]->GetHistogram(),
                graph[4]->GetHistogram(),
                graph[5]->GetHistogram(),
                (TH2D*)(graph[4]->GetHistogram()->Clone("Strength")),
                (TH2D*)(graph[5]->GetHistogram()->Clone("StrengthObs")),
                graph[6]->GetHistogram(),
                graph[7]->GetHistogram()

             };
             h2ds[6]->Divide(h2ds[2]); //signal strength
             h2ds[7]->Divide(h2ds[2]); //signal strength
             h2ds[8]->Divide(h2ds[2]); //signal strength
             h2ds[9]->Divide(h2ds[2]); //signal strength


             printf("DEBUG 2HDM %f %f %f\n", graph[6]->Interpolate(300, 1.0), graph[4]->Interpolate(300, 1.0), graph[7]->Interpolate(300, 1.0));

             for(int Limittype=0;Limittype<8;Limittype++){
                c1 = new TCanvas("c", "c",600,600);              
                c1->SetLogz(true);      
                c1->SetLogy(true);      
                c1->SetRightMargin(0.17);
                framework2d = new TH2F("Graph","Graph",1,200,  600, 1,0.5, 10);
                framework2d->SetStats(false);
                framework2d->SetTitle("");
                framework2d->GetXaxis()->SetTitle("M_{H2} [GeV]");          
                framework2d->GetYaxis()->SetTitle("tan( #beta )");
                framework2d->GetYaxis()->SetTitleOffset(1.60);
                framework2d->Draw("");

                TH2D* h2d = h2ds[Limittype];
                if(Limittype>=4)HideUnexplored(h2d, h2ds[3]);

                double zmin=1E-1; double zmax=1E3;
                if(Limittype==0){ h2d->GetZaxis()->SetTitle("#sigma(gg #rightarrow H2) (fb)" );  zmin=1E-2; zmax=1E5; }
                if(Limittype==1){ h2d->GetZaxis()->SetTitle("BR(H2 #rightarrow ZZ)" );  zmin=1E-4; zmax=1;  }
                if(Limittype==2){ h2d->GetZaxis()->SetTitle("#sigma(gg #rightarrow H2 #rightarrow ZZ) (fb)");  zmin=1E-3; zmax=1E5;  }
                if(Limittype==3){ h2d->GetZaxis()->SetTitle("#Gamma_{H2} / #Gamma_{SM}" );  zmin=5E-3; zmax=2;  }
                if(Limittype==4){ h2d->GetZaxis()->SetTitle("Expected #sigma_{95%} (gg #rightarrow H2 #rightarrow ZZ) (fb)"); }
                if(Limittype==5){ h2d->GetZaxis()->SetTitle("Observed #sigma_{95%} (gg #rightarrow H2 #rightarrow ZZ) (fb)"); }
                if(Limittype==6){ h2d->GetZaxis()->SetTitle("Expected #sigma_{95%} (gg #rightarrow H2 #rightarrow ZZ) / #sigma_{th}"); }
                if(Limittype==7){ h2d->GetZaxis()->SetTitle("Observed #sigma_{95%} (gg #rightarrow H2 #rightarrow ZZ) / #sigma_{th}"); }
                h2d->GetZaxis()->SetTitleOffset(1.55);
                h2d->Draw("COLZ same");
                if(zmin!=zmax)h2d->GetZaxis()->SetRangeUser(zmin, zmax);

                TPaveText* pave1 = new TPaveText(0.5, 0.82,0.65,0.76,"NDC");
                pave1->SetBorderSize(0);
                pave1->SetFillStyle(0);
                pave1->SetTextAlign(22);
                pave1->SetTextFont(62);
                pave1->SetTextSize(0.04);
                pave1->AddText(type==1?"2HDM Type-1":"2HDM Type-2");
                char textBuffer[512]; sprintf(textBuffer, "cos(#beta - #alpha) = %6.2f", cosBmA);
                pave1->AddText(textBuffer);
                pave1->Draw("same");


                if(Limittype>=4 && Limittype<=7){ getContour(h2ds[6], c1, 2, 7, 1);  getContour(h2ds[8], c1, 1, 2, 1);  getContour(h2ds[9], c1, 1, 2, 1);    }//expected contour
                if(Limittype==5 || Limittype==7){ getContour(h2ds[7], c1, 3, 1, 1); }
              
               //drawPointOnGrid(graph);           
               utils::root::DrawPreliminary(2268.759, 13, gPad->GetLeftMargin(),gPad->GetBottomMargin(),gPad->GetRightMargin(),gPad->GetTopMargin()+0.03);


               char Str[512];
               sprintf(Str, "type_%i_cos%3.1f", type, cosBmA);
               if(Limittype==0)sprintf(Str, "%s_Thxsec", Str);
               if(Limittype==1)sprintf(Str, "%s_Thbr", Str);
               if(Limittype==2)sprintf(Str, "%s_ThxsecXbr", Str);            
               if(Limittype==3)sprintf(Str, "%s_Thwidth", Str);            
               if(Limittype==4)sprintf(Str, "%s", Str);            
               if(Limittype==5)sprintf(Str, "%s_Obs", Str);            
               if(Limittype==6)sprintf(Str, "%s_Strength", Str);            
               if(Limittype==7)sprintf(Str, "%s_Obs_Strength", Str);            

               gPad->RedrawAxis();
               c1->SaveAs((SaveDir+"/Stength_2HDM_"+Str+".png").c_str());
               c1->SaveAs((SaveDir+"/Stength_2HDM_"+Str+".pdf").c_str());
               //c1->SaveAs((SaveDir+"/Stength_2HDM_"+Str+".C").c_str());
             }//end of limit type
           }//end of cos loop
         }//end of type loop
      }//gg production



   }//end of directory loop
}






