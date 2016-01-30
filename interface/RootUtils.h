#ifndef _rootutils_h_
#define _rootutils_h_

#include "TH1.h"
#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TText.h"
#include "TF1.h"
#include "TCutG.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TPaveText.h"
#include "THStack.h"

#include <vector>
#include <list>
#include <iterator>
#include <algorithm>
#include <unordered_map>
#include <iostream>
#include <string>
#include <regex>

#include "UserCode/llvv_fwk/interface/JSONWrapper.h"

namespace utils
{
  namespace root
  {
    //
    void fixExtremities(TH1* h,bool addOverflow, bool addUnderflow){
      if(!h) return;
      
      if(addUnderflow){
	double fbin  = h->GetBinContent(0) + h->GetBinContent(1);
	double fbine = sqrt(h->GetBinError(0)*h->GetBinError(0) + h->GetBinError(1)*h->GetBinError(1));
	h->SetBinContent(1,fbin);
	h->SetBinError(1,fbine);
	h->SetBinContent(0,0);
	h->SetBinError(0,0);
      }
      
      if(addOverflow){  
	int nbins = h->GetNbinsX();
	double fbin  = h->GetBinContent(nbins) + h->GetBinContent(nbins+1);
	double fbine = sqrt(h->GetBinError(nbins)*h->GetBinError(nbins) + h->GetBinError(nbins+1)*h->GetBinError(nbins+1));

	h->SetBinContent(nbins,fbin);
	h->SetBinError(nbins,fbine);
        h->GetXaxis()->SetCanExtend(false);  //this is needed, otherwise it doubles the bins of histo with labels at the next lines
	h->SetBinContent(nbins+1,0);
	h->SetBinError(nbins+1,0);
      }
    }


    void checkSumw2(TH1 *h){
       if(h && h->GetDefaultSumw2()) h->Sumw2();  
    }

   bool getMatchingKeyword(JSONWrapper::Object& process, std::vector<std::string>& keywords, std::string& matching){
      matching = "";
      if(keywords.size()<=0)return true;
      if(!process.isTag("keys"))return false;
      std::vector<JSONWrapper::Object> dsetkeywords = process["keys"].daughters();
      for(size_t ikey=0; ikey<dsetkeywords.size(); ikey++){
         for(unsigned int i=0;i<keywords.size();i++){if(std::regex_match(dsetkeywords[ikey].toString(),std::regex(keywords[i]))){matching=dsetkeywords[ikey].toString(); return true;}}
      }
      return false;
   }

   bool matchKeyword(JSONWrapper::Object& process, std::vector<std::string>& keywords){
      std::string unused;
      return getMatchingKeyword(process, keywords, unused); 
   }


   std::string dropBadCharacters(std::string in){
      while(in.find("*")!=std::string::npos)in.replace(in.find("*"),1,"");
      while(in.find("#")!=std::string::npos)in.replace(in.find("#"),1,"");
      while(in.find("{")!=std::string::npos)in.replace(in.find("{"),1,"");
      while(in.find("}")!=std::string::npos)in.replace(in.find("}"),1,"");
      while(in.find("(")!=std::string::npos)in.replace(in.find("("),1,"");
      while(in.find(")")!=std::string::npos)in.replace(in.find(")"),1,"");
      while(in.find("^")!=std::string::npos)in.replace(in.find("^"),1,"");
      while(in.find("/")!=std::string::npos)in.replace(in.find("/"),1,"-");
      return in;
   }


   void setStyle(JSONWrapper::Object& SingleProcess, TH1* hist){
      hist->SetLineColor  (1);
      if(SingleProcess.isTag("color" ) )hist->SetMarkerColor((int)SingleProcess[ "color"].toDouble()); else hist->SetMarkerColor(1);
      if(SingleProcess.isTag("color" ) )hist->SetFillColor  ((int)SingleProcess[ "color"].toDouble()); else hist->SetFillColor  (0);
      if(SingleProcess.isTag("lcolor") )hist->SetLineColor  ((int)SingleProcess["lcolor"].toDouble());
      if(SingleProcess.isTag("mcolor") )hist->SetMarkerColor((int)SingleProcess["mcolor"].toDouble()); 
      if(SingleProcess.isTag("fcolor") )hist->SetFillColor  ((int)SingleProcess["fcolor"].toDouble()); 
      if(SingleProcess.isTag("lwidth") )hist->SetLineWidth  ((int)SingleProcess["lwidth"].toDouble());// else hist->SetLineWidth  (1);
      if(SingleProcess.isTag("lstyle") )hist->SetLineStyle  ((int)SingleProcess["lstyle"].toDouble());// else hist->SetLinStyle  (1);
      if(SingleProcess.isTag("fill"  ) )hist->SetFillColor  ((int)SingleProcess["fill"  ].toDouble());
      if(SingleProcess.isTag("marker") )hist->SetMarkerStyle((int)SingleProcess["marker"].toDouble());// else hist->SetMarkerStyle(1);
      if(SingleProcess.isTag("msize")  )hist->SetMarkerSize (      SingleProcess["msize"].toDouble());// else the Size is the defaoult one
   }

    void setStyleFromKeyword(std::string keyword, JSONWrapper::Object& SingleProcess, TH1* hist){
      hist->SetLineColor  (1);
      if(SingleProcess.isTagFromKeyword(keyword, "color" ) )hist->SetMarkerColor((int)SingleProcess.getDoubleFromKeyword(keyword, "color", 1)); else hist->SetMarkerColor(1);
      if(SingleProcess.isTagFromKeyword(keyword, "color" ) )hist->SetFillColor  ((int)SingleProcess.getDoubleFromKeyword(keyword, "color", 1)); else hist->SetFillColor  (0);
      if(SingleProcess.isTagFromKeyword(keyword, "lcolor") )hist->SetLineColor  ((int)SingleProcess.getDoubleFromKeyword(keyword, "lcolor", 1));
      if(SingleProcess.isTagFromKeyword(keyword, "mcolor") )hist->SetMarkerColor((int)SingleProcess.getDoubleFromKeyword(keyword, "mcolor", 1));
      if(SingleProcess.isTagFromKeyword(keyword, "fcolor") )hist->SetFillColor  ((int)SingleProcess.getDoubleFromKeyword(keyword, "fcolor", 1));
      if(SingleProcess.isTagFromKeyword(keyword, "lwidth") )hist->SetLineWidth  ((int)SingleProcess.getDoubleFromKeyword(keyword, "lwidth", 1));// else hist->SetLineWidth  (1);
      if(SingleProcess.isTagFromKeyword(keyword, "lstyle") )hist->SetLineStyle  ((int)SingleProcess.getDoubleFromKeyword(keyword, "lstyle", 1));// else hist->SetLinStyle  (1);
      if(SingleProcess.isTagFromKeyword(keyword, "fill"  ) )hist->SetFillColor  ((int)SingleProcess.getDoubleFromKeyword(keyword, "fill", 1));
      if(SingleProcess.isTagFromKeyword(keyword, "marker") )hist->SetMarkerStyle((int)SingleProcess.getDoubleFromKeyword(keyword, "marker", 1));// else hist->SetMarkerStyle(1);
      if(SingleProcess.isTagFromKeyword(keyword, "msize")  )hist->SetMarkerSize (     SingleProcess.getDoubleFromKeyword(keyword, "msize", 1));// else the Size is the defaoult one
   }
  

   double getXsecXbr(JSONWrapper::Object& SingleProcess){
      double toReturn = 1.0;
      if(!SingleProcess.isTag("data"))return toReturn;
      toReturn*= SingleProcess["data"][0].getDouble("xsec",1.0);
      if(SingleProcess["data"][0].isTag("br")){
         std::vector<JSONWrapper::Object> BRs = SingleProcess["data"][0]["br"].daughters();
         for(size_t ipbr=0; ipbr<BRs.size(); ipbr++){toReturn*=BRs[ipbr].toDouble();}
      }
       return toReturn;
   }



   // function that add the TPaveText on the current canvas with the "CMS Preliminary...." on top of the Histograms. For split Lumi
   void DrawPreliminary(double iLumi, double iEcm, double L, double B, double R, double T, bool isSim=false, bool preliminary=true){
        //TOP RIGHT OUT-FRAME
        char LumiText[1024];  sprintf(LumiText, "%.1f %s^{-1} (%.0f TeV)", iLumi>100?iLumi/1000:iLumi, iLumi>100?"fb":"pb", iEcm);       
        TPaveText* T1 = new TPaveText(1.0-R-0.50, 1.0-T, 1.02-R, 1.0-0.005, "NDC");
        T1->SetTextFont(43); T1->SetTextSize(23);   T1->SetTextAlign(31);
        T1->SetFillColor(0); T1->SetFillStyle(0);   T1->SetBorderSize(0);
        T1->AddText(LumiText);  T1->Draw();

        //TOP LEFT IN-FRAME
        TPaveText* T2 = new TPaveText(L+0.005, 1.0-T-0.05, L+0.20, 1.0-T-0.005, "NDC");
        T2->SetTextFont(63); T2->SetTextSize(30);   T2->SetTextAlign(11);
        T2->SetFillColor(0); T2->SetFillStyle(0);   T2->SetBorderSize(0);
        T2->AddText("CMS"); T2->Draw();

        if(preliminary){ //Bellow CMS
        TPaveText* T3 = new TPaveText(L+0.005, 1.0-T-0.085, L+0.20, 1.0-T-0.035, "NDC");
        T3->SetTextFont(53); T3->SetTextSize(23);   T3->SetTextAlign(11);
        T3->SetFillColor(0); T3->SetFillStyle(0);   T3->SetBorderSize(0);
        T3->AddText("Preliminary"); T3->Draw();
        }

        //if(Text!=""){  //TOP right IN-FRAME
        //TPaveText* T4 = new TPaveText(1.0-R-0.50, 1.0-T-0.06, 1.0-R, 1.0-T-0.01, "NDC");
        //T4->SetTextFont(43); T4->SetTextSize(23);   T4->SetTextAlign(32);
        //T4->SetFillColor(0); T4->SetFillStyle(0);   T4->SetBorderSize(0);
        //T4->AddText(Text.c_str());  T4->Draw();
        //}

        if(isSim){ //Right of CMS
        TPaveText* T5 = new TPaveText(L+0.085, 1.0-T-0.05, L+0.32, 1.0-T-0.005, "NDC");
        T5->SetTextFont(53); T5->SetTextSize(23);   T5->SetTextAlign(11);
        T5->SetFillColor(0); T5->SetFillStyle(0);   T5->SetBorderSize(0);
        T5->AddText("Simulation"); T5->Draw();
        }
   }

   void DrawPreliminary(double iLumi, double iEcm, TAttPad* pad=NULL, bool isSim=false,   bool preliminary=true){
      if     (pad )DrawPreliminary(iLumi, iEcm, pad->GetLeftMargin(), pad->GetBottomMargin(), pad->GetRightMargin(), pad->GetTopMargin() , isSim, preliminary);
      else if(gPad)DrawPreliminary(iLumi, iEcm, gPad->GetLeftMargin(),gPad->GetBottomMargin(),gPad->GetRightMargin(),gPad->GetTopMargin(), isSim, preliminary);
      else         DrawPreliminary(iLumi, iEcm, 0.15, 0.15, 0.15, 0.15                                                                   , isSim, preliminary);
   }

   TObject* GetObjectFromPath(TDirectory* File, std::string Path, bool GetACopy=false)
   {
      size_t pos = Path.find("/");
      if(pos < 256){
         std::string firstPart = Path.substr(0,pos);
         std::string endPart   = Path.substr(pos+1,Path.length());
         TDirectory* TMP = (TDirectory*)File->Get(firstPart.c_str());
         if(TMP!=NULL){
            TObject* TMP2 =  GetObjectFromPath(TMP,endPart,GetACopy);
            return TMP2;
         }
         return NULL;
      }else{
         TObject* TMP = File->Get(Path.c_str());
         if(GetACopy){	return TMP->Clone();
         }else{            return TMP;
         }
      }
   }

  }
}

#endif
