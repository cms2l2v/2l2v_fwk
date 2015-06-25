
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
//#include "tdrstyle.C"


TGraph* getGraph(string name, int color, int width, int style, TLegend* LEG, TGraph* Ref, int type, string filePath){
//   filePath+="/LimitSummary";
   FILE* pFile = fopen(filePath.c_str(),"r");
   if(!pFile){printf("Can't open %s\n",filePath.c_str()); exit(0);}
   double mass, th, exp, obs, unused;// char buffer[1024];
   double explow2, explow1, expup1, expup2;

   TGraph* graph = new TGraph(100);
   int N=0;
   while(fscanf(pFile,"$%le$ & $%le$ & $[%le,%le]$ & $[%le,%le]$ & $%le$ & Th=$%le$\\\\\\hline\n",&mass, &exp, &explow1, &expup1, &explow2,&expup2,&obs, &th) != EOF){
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

int mapIndex(double CP, double BR){return CP<=0?-1:(int)(CP*100+BR*10);}

void finalLimitPlot(){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //   setTDRStyle();
   gStyle->SetPadTopMargin   (0.04);
   gStyle->SetPadBottomMargin(0.12);
   gStyle->SetPadRightMargin (0.05);
   gStyle->SetPadLeftMargin  (0.12);
   gStyle->SetTitleSize(0.04, "XYZ");
   gStyle->SetTitleXOffset(1.1);
   gStyle->SetTitleYOffset(1.45);
   gStyle->SetPalette(1);
   gStyle->SetNdivisions(505);

   TCanvas* c1;
   TH1F* framework;
   TH2F* framework2d;
   TLegend* LEG, *LEGTH;
   TGraph* Ref;

   double CPs[] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
   double BRs[] = {0.0,0.1,0.2,0.3,0.4,0.5};

  //LIMIT ON SIGNAL STRENGTH
   string Directories2[]={"cards_SB8TeV", "cards_SB8TeV_GG", "cards_SB8TeV_QQ"};
   for(unsigned int D=0;D<sizeof(Directories2)/sizeof(string);D++){
      string Dir = Directories2[D];

      std::map<int, int> colorMap;
      colorMap[mapIndex(1.0, 0.0)] = 2;
      colorMap[mapIndex(0.9, 0.0)] = 2;
      colorMap[mapIndex(0.8, 0.0)] = 3;
      colorMap[mapIndex(0.7, 0.0)] = 3;
      colorMap[mapIndex(0.6, 0.0)] = 4;
      colorMap[mapIndex(0.5, 0.0)] = 4;
      colorMap[mapIndex(0.4, 0.0)] = 6;
      colorMap[mapIndex(0.3, 0.0)] = 6;
      colorMap[mapIndex(0.2, 0.0)] = 8;
      colorMap[mapIndex(0.1, 0.0)] = 8;


      std::map<int, TGraph**> gCPBR;
      gCPBR[ mapIndex(-1.0, 0.0) ]=getGraphs("SM"                         , 1, 2, NULL  , NULL, Dir+               "/Stength_LimitSummary");      
      for(int CPi = 0; CPi<(sizeof(CPs)/sizeof(double));CPi++){
      for(int BRi = 0; BRi<(sizeof(BRs)/sizeof(double));BRi++){
         double CP = CPs[CPi];  double BR = BRs[BRi]; 
         if ( (CP / sqrt(1-BR))>1.0 )continue; //skip point leading to a width larger than SM
         char limitpath[1024]; sprintf(limitpath, "%s_cp%4.2f_brn%4.2f/Stength_LimitSummary", Dir.c_str(), CP, BR);
         char limitlegend[1024]; sprintf(limitlegend, "C'=%3.1f BRnew=%3.1f", CP, BR);
         gCPBR[ mapIndex(CP, BR) ] = getGraphs(limitlegend, colorMap[mapIndex(CP, BR)], 2, NULL  , NULL, limitpath);
      }}


      char LumiLabel[1024];
      if(Dir.find("7TeV")!=string::npos)sprintf(LumiLabel,"CMS preliminary,  #sqrt{s}=7 TeV, #int L=%6.1ffb^{-1}",5.035);
      if(Dir.find("8TeV")!=string::npos)sprintf(LumiLabel,"CMS preliminary,  #sqrt{s}=8 TeV, #int L=%6.1ffb^{-1}",19.6);
      if(Dir.find("Comb")!=string::npos)sprintf(LumiLabel,"CMS preliminary,  #sqrt{s}=%.0f TeV #scale[0.5]{#int} L=%6.1ffb^{-1}, #sqrt{s}=%.0f TeV #scale[0.5]{#int} L=%6.1ffb^{-1}",7.0,5.0,8.0,19.7);
      TPaveText *pave = new TPaveText(0.1,0.96,0.94,0.99,"NDC");
      pave->SetBorderSize(0);
      pave->SetFillStyle(0);
      pave->SetTextAlign(32);
      pave->SetTextFont(42);
      pave->AddText(LumiLabel);
      pave->Draw("same");

          ///////////////////////////////////////////////
          //limit Versus mass for different C' values
          ///////////////////////////////////////////////

      for(int observed=0;observed<=1;observed++){
         c1 = new TCanvas("c", "c",600,600);
         c1->SetLogy(true);
         framework = new TH1F("Graph","Graph",1,150,1050);
         framework->SetStats(false);
         framework->SetTitle("");
         framework->GetXaxis()->SetTitle("M_{H} [GeV/c^{2}]");
         framework->GetYaxis()->SetTitle("#sigma_{95%} / #sigma_{TH}");
         framework->GetYaxis()->SetTitleOffset(1.40);
         framework->GetYaxis()->SetRangeUser(0.1,500);
         framework->Draw();

         LEG = new TLegend(0.70,0.70,0.95,0.94);
         LEG->SetFillStyle(0);
         LEG->SetBorderSize(0);
         for(double cp=0.0;cp<=1.0;cp+=0.2){
            TGraph* g = (gCPBR[ mapIndex(cp, 0.0) ])[1+observed];
            g->SetLineStyle(1);
            g->Draw("C same");
            LEG->AddEntry(g, g->GetTitle(), "L");
         }
         LEG->SetHeader(observed==0?"Expected @95% CL":"Observed @95% CL");
         LEG  ->Draw("same");

         TLine* SMLine = new TLine(framework->GetXaxis()->GetXmin(),1.0,framework->GetXaxis()->GetXmax(),1.0);
         SMLine->SetLineWidth(2); SMLine->SetLineStyle(2); SMLine->SetLineColor(1);
         SMLine->Draw("same C");

         //system("mkdir -p LimitPlots");
         if(observed==0){
            c1->SaveAs((Dir+"/Stength_FinalPlot.png").c_str());
            c1->SaveAs((Dir+"/Stength_FinalPlot.pdf").c_str());
            c1->SaveAs((Dir+"/Stength_FinalPlot.C"  ).c_str());
         }else{
            c1->SaveAs((Dir+"/Stength_FinalPlot_Obs.png").c_str());
            c1->SaveAs((Dir+"/Stength_FinalPlot_Obs.pdf").c_str());
            c1->SaveAs((Dir+"/Stength_FinalPlot_Obs.C"  ).c_str());
         }
     }

     ///////////////////////////////////////////////
     //Mass Versus Cprime limits
     ///////////////////////////////////////////////

     TGraph* gMvsCp[7];
     for(int type=1;type<7;type++){
        if(type==0)continue;
        double* Masses = (gCPBR[ mapIndex(1.0, 0.0) ])[0]->GetX();  
        int NMasses = (gCPBR[ mapIndex(1.0, 0.0) ])[0]->GetN();
           
        gMvsCp[type] = new TGraph(NMasses);
        for(int Mi=0;Mi<NMasses;Mi++){
           //create a 1d graph per mass point having C' on x axis and limit on Y axis, such that we can find the C' value for which we have the intersection of the limit with "1" 
           TGraph* g1d = new TGraph( (sizeof(CPs)/sizeof(double)) );
           for(unsigned int CPi=0;CPi<(sizeof(CPs)/sizeof(double));CPi++){
              g1d->SetPoint(CPi, CPs[CPi], (gCPBR[ mapIndex(CPs[CPi], 0.0) ])[type]->Eval(Masses[Mi]));
           }

           //find C' value that is excluded
           double cpExcluded = 2.0;  double cpPrev=-1;  double limitPrev=-1;
           for(double cp=1.0; cp>0.0;cp-=0.01){
              if(limitPrev!=-1 && limitPrev<=1.0 && g1d->Eval(cp)>=1.0){cpExcluded=cpPrev;break;}
              limitPrev=g1d->Eval(cp);  cpPrev=cp;
           }
           gMvsCp[type]->SetPoint(Mi, Masses[Mi], pow(cpExcluded,2));  //show the limit as a function of C'²
//         gMvsCp[type]->SetPoint(Mi, Masses[Mi], cpExcluded);         //show the limit as a function of C'

           delete g1d;  // not needed anymore
        }
     }

     //now can take care of the plotting
     if(true){ 
        c1 = new TCanvas("c", "c",600,600);
        framework = new TH1F("Graph","Graph",1,150,1050);
        framework->SetStats(false);
        framework->SetTitle("");
        framework->GetXaxis()->SetTitle("M_{H} [GeV/c^{2}]");
        framework->GetYaxis()->SetTitle("C'^{2}_{95%}");
        framework->GetYaxis()->SetTitleOffset(1.40);
        framework->GetYaxis()->SetRangeUser(0.0,1.0);
        framework->Draw();
        pave->Draw("same");

        //expected bands
        TCutG* TGExpLimit1S  = GetErrorBand("1S", gMvsCp[3]->GetN(), gMvsCp[3]->GetX(), gMvsCp[3]->GetY(), gMvsCp[4]->GetY());  
        TCutG* TGExpLimit2S  = GetErrorBand("2S", gMvsCp[3]->GetN(), gMvsCp[3]->GetX(), gMvsCp[5]->GetY(), gMvsCp[6]->GetY());  TGExpLimit2S->SetFillColor(5);
        TGExpLimit2S->Draw("f same");
        TGExpLimit1S->Draw("f same");

        //expected
        gMvsCp[1]->SetLineColor(1);
        gMvsCp[1]->SetLineWidth(1);
        gMvsCp[1]->SetLineStyle(2);
        gMvsCp[1]->Draw("same");

        //observed
        gMvsCp[2]->SetLineColor(1);
        gMvsCp[2]->SetLineWidth(2);
        gMvsCp[2]->SetLineStyle(1);
        gMvsCp[2]->Draw("same");

        LEG = new TLegend(0.50,0.20,0.95,0.44);
        LEG->SetFillStyle(0);
        LEG->SetBorderSize(0);
        LEG->SetHeader(NULL);
        LEG->AddEntry(gMvsCp[1],"Expected @95% CL", "L");
        LEG->AddEntry(gMvsCp[2],"Observed @95% CL", "L");
        LEG->Draw();

        c1->SaveAs((Dir+"/Stength_FinalPlot_Cprime.png").c_str());
        c1->SaveAs((Dir+"/Stength_FinalPlot_Cprime.pdf").c_str());
        c1->SaveAs((Dir+"/Stength_FinalPlot_Cprime.C"  ).c_str());
     }

      ///////////////////////////////////////////////
     //Mass Versus Cprime limits
     ///////////////////////////////////////////////

      for(int observed=0;observed<=1;observed++){
         for(int mode=0; mode<=3; mode++){
            //mode=0 --> CPrime versus BRNew  interpolated from 2D limit(CPrime, BRNew)
            //mode=1 --> CPrime versus BRNew  interpolated from 1D limit(CPrime, BRNew=0)
            //mode=2 --> width  versus BRNew  interpolated from 2D limit(CPrime, BRNew)
            //mode=3 --> width  versus BRNew  interpolated from 1D limit(CPrime, BRNew=0)
////
            double Masses[] = {200,400,600,800};
            for(int Mi=0;Mi<sizeof(Masses)/sizeof(double);Mi++){
               double Mass = Masses[Mi];

               int       t2dI = 0;
               TGraph2D* t2d  = new TGraph2D( (sizeof(CPs)/sizeof(double)) * (sizeof(BRs)/sizeof(double))  );
               TGraph*   g1d  = new TGraph  ( (sizeof(CPs)/sizeof(double))  );
               for(unsigned int CPi=0;CPi<(sizeof(CPs)/sizeof(double));CPi++){                    
                  g1d->SetPoint(CPi, CPs[CPi], (gCPBR[ mapIndex(CPs[CPi], 0.0) ])[1+observed]->Eval(Mass));

                  for(unsigned int BRi=0;BRi<(sizeof(BRs)/sizeof(double));BRi++){
                     if ( (CPs[CPi] / sqrt(1-BRs[BRi]))>1.0 )continue; //skip point leading to a width larger than SM
                     t2d->SetPoint(t2dI, CPs[CPi], BRs[BRi], (gCPBR[ mapIndex(CPs[CPi], BRs[BRi]) ])[1+observed]->Eval(Mass));        
                     t2dI++;
               }}t2d->Set(t2dI);

               int I=0;      
               double   cp2[] = {1.0, 0.975, 0.95, 0.925, 0.9, 0.875, 0.85, 0.825, 0.8, 0.775, 0.75, 0.725, 0.7, 0.675, 0.65, 0.625, 0.6, 0.575, 0.55, 0.525, 0.5, 0.475, 0.45, 0.425, 0.4, 0.375, 0.35, 0.325, 0.3, 0.275, 0.25, 0.225, 0.2, 0.175, 0.15, 0.125, 0.1, 0.075, 0.05, 0.025, 0.0};
               double  brn[] = {0.995, 0.975, 0.95, 0.925, 0.9, 0.875, 0.85, 0.825, 0.8, 0.775, 0.75, 0.725, 0.7, 0.675, 0.65, 0.625, 0.6, 0.575, 0.55, 0.525, 0.5, 0.475, 0.45, 0.425, 0.4, 0.375, 0.35, 0.325, 0.3, 0.275, 0.25, 0.225, 0.2, 0.175, 0.15, 0.125, 0.1, 0.075, 0.05, 0.025, 0.0};
               TGraph2D* g2d = new TGraph2D( (sizeof(cp2)/sizeof(double)) * (sizeof(brn)/sizeof(double))  );
               for(unsigned int C=0;C<(sizeof( cp2)/sizeof(double));C++){
               for(unsigned int B=0;B<(sizeof(brn)/sizeof(double));B++){
                  double BR = brn[B];

                  //move from br=0 --> br!=0
                  double cp2_true = cp2[C]*(1-BR);
                  double sl       = g1d->Eval(sqrt(cp2_true));    //THERE WAS A BUG in THE PAS... since we were using g1d->Eval(sqrt(cp2[C]));
                  double sl_true  = sl / (1-BR);

                  if(mode==0){        g2d->SetPoint(I, cp2[C], BR, t2d->Interpolate(sqrt(cp2[C]), BR) );
                  }else if(mode==1){  g2d->SetPoint(I, cp2_true, BR, sl_true);
                  }else if(mode==2){  g2d->SetPoint(I, sqrt(cp2[C]), BR, t2d->Interpolate(sqrt(cp2[C]), BR) );
                  }else if(mode==3){  g2d->SetPoint(I, sqrt(cp2_true), BR, sl_true);  //THERE WAS A BUG IN THE PAS.... we were using sqrt(cp2[C]), BR, sl_true);
                  }
                  I++;
               }}g2d->Set(I);

               c1 = new TCanvas("c", "c",600,600);
               c1->SetLogz(true);      
               c1->SetRightMargin(0.17);
//               framework2d = new TH2F("Graph","Graph",1,mode<=1?0.0:0.3,1, 1,0,mode%2==0?0.5:1.0);
               framework2d = new TH2F("Graph","Graph",1,mode<=1?0.0:0.3,1, 1,0,1.0);
               framework2d->SetStats(false);
               framework2d->SetTitle("");
               framework2d->GetXaxis()->SetTitle(mode<=1?"c'^{2}":"#Gamma/#Gamma_{SM}");
               framework2d->GetYaxis()->SetTitle("BR_{new}");
               framework2d->GetYaxis()->SetTitleOffset(1.40);
               framework2d->Draw("");
               pave->Draw("same");

               TH2D* h2d = g2d->GetHistogram();
               h2d->SetMaximum(10);
               h2d->SetMinimum(1E-1);
               h2d->GetZaxis()->SetTitle(observed==0?"Expected #sigma_{95%} / #sigma_{TH}":"Observed #sigma_{95%} / #sigma_{TH}");
               h2d->GetZaxis()->SetTitleOffset(1.33);
               h2d->Draw("COLZ same");

               TH1D* h2dLimit = new TH1D("ExcludedArea", "ExcludedArea", h2d->GetNbinsX(), 0.0, 1.0);  
               for(int x=0;x<=h2d->GetNbinsX();x++){
                  double limit = -1;
                  for(int y=0;y<=h2d->GetNbinsY()+1;y++){
                     double bin = h2d->GetBinContent(x,y);
//                     if(x==15)printf("Dir=%s Mass=%f: Mode=%i cprime2 = %f  brnew=%f  limit=%f\n", Dir.c_str(), Mass, mode, h2d->GetXaxis()->GetBinCenter(x), h2d->GetYaxis()->GetBinCenter(y), bin );
//                     if(mode<=1 && (bin>=1 || h2d->GetYaxis()->GetBinCenter(y)>=1-h2d->GetXaxis()->GetBinCenter(x))){limit = std::max(0.0, h2d->GetYaxis()->GetBinLowEdge(y) ); break;}
//                     if(mode>=2 && (bin<=1 && y<h2d->GetNbinsY())){limit = std::max(0.0, h2d->GetYaxis()->GetBinUpEdge(y) );}
                       if(bin>0.0 && bin<=1.0){limit = h2d->GetYaxis()->GetBinUpEdge(y);}
                  }
                  h2dLimit->SetBinContent(x,limit);
               } 
               h2dLimit->SetLineColor(1);
               h2dLimit->SetLineWidth(1);
               h2dLimit->SetFillStyle(3654);
               h2dLimit->SetFillColor(1);
               h2dLimit->Draw("HIST same");
               
               if(mode<=1){
                  TLine* W50 = new TLine(0.0, 1.0, 1.0, 0.8);
                  W50->SetLineWidth(1); W50->SetLineStyle(3); W50->SetLineColor(1);    W50->Draw("same C");

                  TLine* W10 = new TLine(0.0, 1.00, 1.0, 0.0);
                  W10->SetLineWidth(1); W10->SetLineStyle(3); W10->SetLineColor(1);    W10->Draw("same C");

                  TLine* W05 = new TLine(0.0, 1.0, 0.5, 0.0);
                  W05->SetLineWidth(1); W05->SetLineStyle(3); W05->SetLineColor(1);    W05->Draw("same C");

                  TLine* W01 = new TLine(0.0, 1.0, 0.1, 0.0);
                  W01->SetLineWidth(1); W01->SetLineStyle(3); W01->SetLineColor(1);    W01->Draw("same C");
               }

//               //add marker on the plots
//               for(unsigned int C=0;C<(sizeof( cp2)/sizeof(double));C++){
//               for(unsigned int B=0;B<(sizeof(brn)/sizeof(double));B++){
//                  double BR = brn[B];
//                  double sl = g1d->Eval(sqrt(cp2[C]));
//                  //move from br=0 --> br=BR
//                  double cp2_true = cp2[C]*(1-BR);
//                  double sl_true = sl/(1-BR);
//                  //TMarker* m = new TMarker(mode<=1?cp2_true:sqrt(cp2[C]), BR, 20); m->SetMarkerSize(1.0); m->SetMarkerColor(1);      m->Draw("same");
//               }}


               char massStr[512]; sprintf(massStr, "%04.0f", Mass);
               if(mode==0)sprintf(massStr, "%s_2DInterpol", massStr);
               if(mode==1)sprintf(massStr, "%s_1DInterpol", massStr);
               if(mode==2)sprintf(massStr, "%s_Width", massStr);
               if(mode==3)sprintf(massStr, "%s_WidthFrom1DInterpol", massStr);
               if(observed!=0)sprintf(massStr, "%s_Obs", massStr);
               c1->SaveAs((Dir+"/Stength_FinalPlot2D_"+massStr+".png").c_str());
               c1->SaveAs((Dir+"/Stength_FinalPlot2D_"+massStr+".pdf").c_str());
               c1->SaveAs((Dir+"/Stength_FinalPlot2D_"+massStr+".C").c_str());
            }//end of mass Loop
         }//end of mode loop
      }//end of observed loop

   }//end of directory loop
}






