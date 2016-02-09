
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

#include "UserCode/llvv_fwk/interface/tdrstyle.h"
#include "UserCode/llvv_fwk/src/tdrstyle.C"
#include "UserCode/llvv_fwk/interface/RootUtils.h"

#include "UserCode/llvv_fwk/interface/HxswgUtils.h"
#include "UserCode/llvv_fwk/src/HxswgUtils.cc"


TGraph2D* get2HDMLimitsTgVsAmB(TGraph2D* g2dMassVsWidth, TGraph* SMhiggsWidth, double OnlyCosAmB, bool getThXSec=false){
   FILE* pFile = fopen("Weight_2HDM_Model.txt","r");
   if(!pFile){printf("Can't open %s\n","Weight_2HDM_Model.txt"); exit(0);}

   TGraph2D* graph = new TGraph2D(9999);    int N=0;
   double tgB, cosAmB, mH2, width, BRtoZZ, XSec;
   while(fscanf(pFile,"tgBeta=%lf cosAmB=%lf MH2=%lf Width=%lf BrH2toZZ=%lf XSec=%lf\n",&tgB, &cosAmB, &mH2, &width, &BRtoZZ, &XSec) != EOF){
      if(fabs(cosAmB-OnlyCosAmB)<0.05){
         if(getThXSec){
            graph->SetPoint(N, mH2, tgB, 1000*XSec*BRtoZZ);  N++;
         }else{
            graph->SetPoint(N, mH2, tgB, g2dMassVsWidth->Interpolate(mH2, width/SMhiggsWidth->Eval(mH2)) / (XSec*BRtoZZ*1000)   );  N++;
//            printf("Signal Strength limit = %f = %f/%f\n",  g2dMassVsWidth->Interpolate(mH2, width/SMhiggsWidth->Eval(mH2)) / (XSec*BRtoZZ*1000),  g2dMassVsWidth->Interpolate(mH2, width/SMhiggsWidth->Eval(mH2)) , (XSec*BRtoZZ*1000) );
         }
      }
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
	double* Xs = graph->GetX();
	double* Ys = graph->GetY();
	for(int p=0;p<graph->GetN();p++){
           TMarker* mark = new TMarker(Xs[p], Ys[p], 20);
	   mark->SetMarkerSize(0.3);
	   mark->Draw("same");
           //printf("draw a point at %f - %f\n", Xs[p], Ys[p]);
	}
}


TGraph* getGraph(string name, int color, int width, int style, TLegend* LEG, TGraph* Ref, int type, string filePath){
//   filePath+="/LimitSummary";
   FILE* pFile = fopen(filePath.c_str(),"r");
   if(!pFile){printf("Can't open %s\n",filePath.c_str()); exit(0);}
   double mass, th, exp, obs, unused;// char buffer[1024];
   double explow2, explow1, expup1, expup2;

   TGraph* graph = new TGraph(100);
   int N=0;
   while(fscanf(pFile,"$%le$ & $%le$ & $[%le,%le]$ & $[%le,%le]$ & $%le$ & Th=$%le$ & pValue=$%le$\\\\\\hline\n",&mass, &exp, &unused, &unused, &unused,&unused,&obs, &th, &unused) != EOF){

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
	for(unsigned int i=0;i<list->GetSize();i++){
		EXCLUSION   = (TGraph*)(list->At(i)->Clone("copy"));
		EXCLUSION->SetLineColor(Color);
		EXCLUSION->SetLineWidth(Width);
		EXCLUSION->SetLineStyle(Style);
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
   gStyle->SetPalette(1);
   gStyle->SetNdivisions(505);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);

   std::cout<<"A\n";
   TGraph* higgsWidth = Hxswg::utils::getHWidth();
   std::cout<<"B\n";


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

   double CPs[] = {1.0,0.6,0.3,0.1};
   double BRs[] = {0.0};

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
      for(int CPi = 0; CPi<(sizeof(CPs)/sizeof(double));CPi++){
      for(int BRi = 0; BRi<(sizeof(BRs)/sizeof(double));BRi++){
         double CP = CPs[CPi];  double BR = BRs[BRi]; 
         if ( CP*CP / (1-BR)>1.0 )continue; //skip point leading to a width larger than SM
         char limitpath[1024]; sprintf(limitpath, "%s_cp%4.2f_brn%4.2f/Stength_LimitSummary", Dir.c_str(), CP, BR);
//         char limitlegend[1024]; sprintf(limitlegend, "C'=%3.1f BRnew=%3.1f", CP, BR);
         char limitlegend[1024]; sprintf(limitlegend, "C'=%3.1f", CP);
         std::cout<<limitpath<< "\n";         
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
      TGraph2D* g2dMassVsWidth[7];    
      for(int i=0;i<7;i++){
         g2dMassVsCp[i]    = new TGraph2D( 9999 ); 
         g2dMassVsWidth[i] = new TGraph2D( 9999 );        
         int Ig2d=0;
         for(int CPi = 0; CPi<(sizeof(CPs)/sizeof(double));CPi++){
            double CP = CPs[CPi];  double BR = 0.0;             
            TGraph* graph = gCPBR[ mapIndex(CP, BR) ][i];           
            for(int p=0;p<graph->GetN ();p++){
               double mass, limit;
               graph->GetPoint(p, mass, limit);
               g2dMassVsCp[i]   ->SetPoint(Ig2d, mass, CP   , limit);
               g2dMassVsWidth[i]->SetPoint(Ig2d, mass, CP*CP, limit);
               Ig2d++;              
            }
         }
         g2dMassVsCp   [i]->Set(Ig2d);
         g2dMassVsWidth[i]->Set(Ig2d);        
      } 


      ///////////////////////////////////////////////
      //limit Versus mass for different C' values
      ///////////////////////////////////////////////

      for(int observed=0;observed<=1;observed++){
         c1 = new TCanvas("c", "c",600,600);
         c1->SetLogy(true);
         framework = new TH1F("Graph","Graph",1,390,1010);
         framework->SetStats(false);
         framework->SetTitle("");
         framework->GetXaxis()->SetTitle("M_{H} [GeV/c^{2}]");
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


         LEG = new TLegend(0.50,0.70,0.75,0.94);
         LEG->SetFillStyle(0);
         LEG->SetBorderSize(0);
         for(int CPi = 0; CPi<(sizeof(CPs)/sizeof(double));CPi++){
         for(int BRi = 0; BRi<(sizeof(BRs)/sizeof(double));BRi++){
            double CP = CPs[CPi];  double BR = BRs[BRi]; 
            if(BR!=0.0)continue;
            if(!gCPBR[ mapIndex(CP, BR) ])continue;
            TGraph* g = (gCPBR[ mapIndex(CP, BR) ])[1+observed];
            if(!g) continue;
            g->SetLineStyle(2);
            g->Draw("C same");
            LEG->AddEntry(g, g->GetTitle(), "L");            


           if(!strengthLimit){
              TGraph* THXSec = (gCPBR[ mapIndex(CP, BR) ])[0];
                THXSec->SetLineStyle(1);  THXSec->SetLineWidth(1);  THXSec->SetLineColor(g->GetLineColor());
                THXSec->Draw("same C");
           }
         }}
         LEG->SetHeader(observed==0?"Expected @95% CL":"Observed @95% CL");
         LEG  ->Draw("same");

         utils::root::DrawPreliminary(2215, 13);
         if(observed==0){
            c1->SaveAs((SaveDir+"/Stength_FinalPlot.png").c_str());
            c1->SaveAs((SaveDir+"/Stength_FinalPlot.pdf").c_str());
            c1->SaveAs((SaveDir+"/Stength_FinalPlot.C"  ).c_str());
         }else{
            c1->SaveAs((SaveDir+"/Stength_FinalPlot_Obs.png").c_str());
            c1->SaveAs((SaveDir+"/Stength_FinalPlot_Obs.pdf").c_str());
            c1->SaveAs((SaveDir+"/Stength_FinalPlot_Obs.C"  ).c_str());
         }
     }


      ///////////////////////////////////////////////
     //Mass Versus Cprime limits  2D
     ///////////////////////////////////////////////

      for(int observed=0;observed<=1;observed++){
         if(observed>0)continue;//DEBUG
         for(int mode=0; mode<=1; mode++){
           //mode=0 --> Mass versus C'     interpolated from 2D limit(CPrime, BRNew)
           //mode=1 --> Mass versus Gamma  interpolated from 2D limit(CPrime, BRNew)
 
           TH2D* grid   = new TH2D("grid"  , "grid"  , 200, 0, 1006, 200, 0, 1);
           TGraph2D*   graph = g2dMassVsCp   [1+observed];
           if(mode==1) graph = g2dMassVsWidth[1+observed];
           graph->SetHistogram((TH2D*)grid->Clone("GRIDg2dMassVsCp"));
           TH2D* h2d = graph->GetHistogram();

           c1 = new TCanvas("c", "c",600,600);
           c1->SetLogz(true);      
           c1->SetRightMargin(0.17);
           framework2d = new TH2F("Graph","Graph",1,400, 1000, 1,0,1.0);
           framework2d->SetStats(false);
           framework2d->SetTitle("");
           framework2d->GetXaxis()->SetTitle("H mass (GeV)");
           if(mode==0)framework2d->GetYaxis()->SetTitle("C'");
           if(mode==1)framework2d->GetYaxis()->SetTitle("#Gamma/#Gamma_{SM}");          
           framework2d->GetYaxis()->SetTitleOffset(1.40);
           framework2d->Draw("");
           h2d->GetZaxis()->SetTitleOffset(1.33);
           if(strengthLimit){
              h2d->GetZaxis()->SetTitle((string(observed==0?"Expected":"Observed") + " #mu = #sigma_{95%} / #sigma_{th}").c_str());
              h2d->GetZaxis()->SetRangeUser(1E-1,1E2);
           }else{
              h2d->GetZaxis()->SetTitle((string(observed==0?"Expected":"Observed") + string(" #sigma_{95%} (") + prod +" #rightarrow H #rightarrow ZZ) (fb)").c_str());
              h2d->GetZaxis()->SetRangeUser(1E1,3E3);
           }
           h2d->Draw("COLZ same");


           if(strengthLimit){
              getContour(h2d, c1, 3, 1, 1);
           }else{
              //build th cross-section 2D plane
              TGraph* THXSec   = Hxswg::utils::getXSec(Dir); 
              scaleGraph(THXSec, 1000);  //convert cross-section to fb
              TH2D* THXSec2D   = (TH2D*)grid->Clone("ThXSec");
              for(unsigned int x=1;x<=THXSec2D->GetNbinsX();x++){
                 double mass = THXSec2D->GetXaxis()->GetBinCenter(x);
                 for(unsigned int y=1;y<=THXSec2D->GetNbinsY();y++){
                    double CP = THXSec2D->GetYaxis()->GetBinCenter(y);  double BR=0.0;
                    if(mode==1){CP = sqrt(CP);}  //go from width=C'² --> to C'
                    double XSecScaleFactor = pow(CP,2) * (1-BR);
                    THXSec2D->SetBinContent(x,y,THXSec->Eval(mass)*XSecScaleFactor);
              }}

              TH2D* signalStrength = (TH2D*) h2d->Clone("signalStrength");
              signalStrength->Divide(THXSec2D);
              getContour(signalStrength, c1, 3, 1, 1);
           }


            utils::root::DrawPreliminary(2215, 13);         
            char massStr[512];
            if(mode==0)sprintf(massStr, "MassVsCp");
            if(mode==1)sprintf(massStr, "MassVsWidth");          
            c1->SaveAs((SaveDir+"/Stength_FinalPlot2D_"+massStr+".png").c_str());
            c1->SaveAs((SaveDir+"/Stength_FinalPlot2D_"+massStr+".pdf").c_str());
            c1->SaveAs((SaveDir+"/Stength_FinalPlot2D_"+massStr+".C").c_str());
         }//end of mode loop
      }//end of observed loop

 

      ///////////////////////////////////////////////
     //Limits on 2HDM
     ///////////////////////////////////////////////

      if(prod!="gg")continue;

      for(int observed=0;observed<=1;observed++){
         if(observed>0)continue;//DEBUG
         for(double cosBmA=0.1; cosBmA<=0.1; cosBmA+=0.1){
           TH2D* grid   = new TH2D("grid2HDM"  , "grid"  , 1000, 400, 1000, 1000, 0.5, 60);
           TGraph2D* graph = get2HDMLimitsTgVsAmB(g2dMassVsWidth[1+observed], higgsWidth, cosBmA);
           graph->SetHistogram((TH2D*)grid->Clone("GRID2HDMmH"));
           TH2D* h2d = graph->GetHistogram();

           c1 = new TCanvas("c", "c",600,600);
           c1->SetLogz(true);      
           c1->SetLogy(true);      
           c1->SetRightMargin(0.17);
           framework2d = new TH2F("Graph","Graph",1,400, 1000, 1,0.5,60.0);
           framework2d->SetStats(false);
           framework2d->SetTitle("");
           framework2d->GetYaxis()->SetTitle("mH2 [GeV]");          
           framework2d->GetYaxis()->SetTitle("tan( #beta )");
           framework2d->GetYaxis()->SetTitleOffset(1.60);
           framework2d->Draw("");

//           h2d->GetZaxis()->SetTitle((string(observed==0?"Expected":"Observed") + " #sigma_{95%} (" + prod +" #rightarrow H #rightarrow ZZ) (fb)").c_str() );
           h2d->GetZaxis()->SetTitle((string(observed==0?"Expected":"Observed") + " #sigma_{95%} / #sigma_{TH}").c_str() );
           h2d->GetZaxis()->SetTitleOffset(1.33);
           h2d->Draw("COLZ same");

 	   TPaveText* pave1 = new TPaveText(0.5, 0.82,0.65,0.76,"NDC");
  	   pave1->SetBorderSize(0);
 	   pave1->SetFillStyle(0);
	   pave1->SetTextAlign(22);
  	   pave1->SetTextFont(62);
	   pave1->SetTextSize(0.04);
 	   pave1->AddText("2HDM Type-2");
           char textBuffer[512]; sprintf(textBuffer, "cos(#beta - #alpha) = %6.2f", cosBmA);
	   pave1->AddText(textBuffer);
           pave1->Draw("same");

           drawPointOnGrid(graph);           
            utils::root::DrawPreliminary(2215, 13);         
            char Str[512];
            sprintf(Str, "cos%3.1f", cosBmA);
            c1->SaveAs((SaveDir+"/Stength_2HDM_"+Str+".png").c_str());
            c1->SaveAs((SaveDir+"/Stength_2HDM_"+Str+".pdf").c_str());
            c1->SaveAs((SaveDir+"/Stength_2HDM_"+Str+".C").c_str());
         }//end of mode loop
      }//end of observed loop


   }//end of directory loop
}






