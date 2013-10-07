#include <iostream>
#include <boost/shared_ptr.hpp>
#include "Math/GenVector/Boost.h"

#include "UserCode/llvv_fwk/src/tdrstyle.C"
#include "UserCode/llvv_fwk/src/JSONWrapper.cc"
#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/RootUtils.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TString.h"
#include "TList.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TObjArray.h"
#include "THStack.h"
#include "TGraphErrors.h"

#include<iostream>
#include<fstream>
#include<map>
#include<algorithm>
#include<vector>
#include<set>

using namespace std;

TString outUrl("./");
TString suffixUrl(""); 
TString inFileUrl("");
TString systFileUrl("");
TString jsonFileUrl("");
TString histo("finalevtflow");
TString signalTag("");
std::set<TString> systVars;
std::vector<int> binsToProject;
std::vector<std::string> channels;
float lumiUnc(0.044);
float selEffUnc(0.02);
float iEcm(8);
bool doPowers = true;

//wrapper for a projected shape for a given set of cuts
struct Shape_t
{
  TH1F* data, *totalBckg, *totalSplusB, *signal;
  std::vector<TH1F*> bckg;
  //  std::vector<TH1F*> signalMassPoints; 
  
  //the first key corresponds to the proc name
  //the second key is the name of the variation: e.g. jesup, jesdown, etc.
  std::map<TString,std::map<TString, TH1F*> > bckgVars;
  std::set<TString> dataDrivenBckg;
  std::map<TString, TH1F*> signalVars;
  
  //cross section and uncertainties
  std::map<TString,std::pair<float,float> > crossSections;

  //signal rate uncertainties
  std::map<TString,std::pair<float,float> > rateUncs;
};


//
void printHelp();
Shape_t getShapeFromFile(TFile* inF, TString ch,JSONWrapper::Object &Root,TFile *systF=0);
void getYieldsFromShapes(const map<TString, Shape_t> &allShapes);
void convertShapesToDataCards(const map<TString, Shape_t> &allShapes);
void saveShapeForMeasurement(TH1F *h, TDirectory *oDir,TString syst="");
TString convertNameForDataCard(TString title);
TString convertNameForFileName(TString histoName);
TString convertMassPointNameForFileName(TString title);
float getIntegratedSystematics(TH1F *h,const std::map<TString, TH1F*> &hSysts, std::map<TString,float> &rateSysts);
std::map<TString,float> getDYUncertainties(TString ch);

//
std::map<TString,float> getDYUncertainties(TString ch)
{
  //assign the bin to use
  int dybin=-1;
  if(ch=="ee")   dybin=0;
  if(ch=="mumu") dybin=1;
  if(ch=="emu")  dybin=2;

  //build the uncertainty map (for an inclusive channel I'm averaging the unc - it doesn't really matter as it will be combined in the end)
  std::map<TString,float> dysfUncs; 
  float stat[]        = {0.064,      0.038,     0.079};     dysfUncs["stat"]         = (dybin==-1 ? (stat[0]+stat[1]+stat[2])/3             : stat[dybin]);
  float jes[]         = {0.021,      0.019,     0.041};     dysfUncs["jes"]          = (dybin==-1 ? (jes[0]+jes[1]+jes[2])/3                : jes[dybin]);
  float jer[]         = {0.014,      0.014,     0.013 };    dysfUncs["jer"]          = (dybin==-1 ? (jer[0]+jer[1]+jer[2])/3                : jer[dybin]);
  float pu[]          = {0.02,       0.01,      0.008 };    dysfUncs["pu"]           = (dybin==-1 ? (pu[0]+pu[1]+pu[2])/3                   : pu[dybin]);
  float mcsignal[]    = {0.087,      0.010,     0.053 };    dysfUncs["mcsignal"]     = (dybin==-1 ? (mcsignal[0]+mcsignal[1]+mcsignal[2])/3 : mcsignal[dybin]);
  float q2[]          = {0.014,      0.014,     0.001 };    dysfUncs["q2"]           = (dybin==-1 ? (q2[0]+q2[1]+q2[2])/3                   : q2[dybin]);
  float meps[]        = {0.011,      0.011,     0.001 };    dysfUncs["meps"]         = (dybin==-1 ? (meps[0]+meps[1]+meps[2])/3             : meps[dybin]);
  float dy_template[] = {0.08,       0.06,      0.01 };  
  if(dybin==-1)
    {
      dysfUncs["dy_ee_template"]   = dy_template[0];
      dysfUncs["dy_mumu_template"] = dy_template[1];
      dysfUncs["dy_emu_template"]  = dy_template[2];
    }
  else
    dysfUncs["dy_"+ch+"_template"] = dy_template[dybin];

  return dysfUncs;
}


//
TString convertNameForDataCard(TString title)
{
  if(title=="VV")                                             return "vv";
  if(title=="QCD")                                            return "qcd";
  if(title=="W#rightarrow l#nu")                              return "w"; 
  if(title=="W,multijets")                                    return "wjets";
  if(title=="other t#bar{t}")                                 return "otherttbar";
  if(title=="Z#rightarrow ll")                                return "dy"; 
  if(title=="Single top")                                     return "st";
  if(title=="t#bar{t}V")                                      return "ttv";
  if(title=="t#bar{t}")                                       return "ttbar";
  if(title=="#splitline{H^{+}#rightarrow tb}{[180 GeV]}")     return "HTB";   signalTag = "HTB" ;
  if(title=="#splitline{H^{+}#rightarrow#tau#nu}{[180 GeV]}") return "TBH";   signalTag = "TBH" ;
  if(title=="#splitline{H^{+}#rightarrow tb}{[200 GeV]}")     return "HTB";   signalTag = "HTB" ;
  if(title=="#splitline{H^{+}#rightarrow#tau#nu}{[200 GeV]}") return "TBH";   signalTag = "TBH" ;
  if(title=="#splitline{H^{+}#rightarrow tb}{[220 GeV]}")     return "HTB";   signalTag = "HTB" ;
  if(title=="#splitline{H^{+}#rightarrow#tau#nu}{[220 GeV]}") return "TBH";   signalTag = "TBH" ;
  if(title=="#splitline{H^{+}#rightarrow tb}{[240 GeV]}")     return "HTB";   signalTag = "HTB" ;
  if(title=="#splitline{H^{+}#rightarrow tb}{[250 GeV]}")     return "HTB";   signalTag = "HTB" ;
  if(title=="#splitline{H^{+}#rightarrow#tau#nu}{[250 GeV]}") return "TBH";   signalTag = "TBH" ;
  if(title=="#splitline{H^{+}$#rightarrow tb}{[260 GeV]}")    return "HTB";   signalTag = "HTB" ;
  if(title=="#splitline{H^{+}#rightarrow tb}{[280 GeV]}")     return "HTB";   signalTag = "HTB" ;
  if(title=="#splitline{H^{+}#rightarrow tb}{[300 GeV]}")     return "HTB";   signalTag = "HTB" ;
  if(title=="#splitline{H^{+}#rightarrow#tau#nu}{[300 GeV]}") return "TBH";   signalTag = "TBH" ;
  if(title=="#splitline{H^{+}#rightarrow tb}{[350 GeV]}")     return "HTB";   signalTag = "HTB" ;
  if(title=="#splitline{H^{+}#rightarrow tb}{[400 GeV]}")     return "HTB";   signalTag = "HTB" ;
  if(title=="#splitline{H^{+}#rightarrow tb}{[500 GeV]}")     return "HTB";   signalTag = "HTB" ;
  if(title=="#splitline{H^{+}#rightarrow tb}{[600 GeV]}")     return "HTB";   signalTag = "HTB" ;
  if(title=="#splitline{H^{+}#rightarrow tb}{[700 GeV]}")     return "HTB";   signalTag = "HTB" ; 
  return title;
}
// 
TString convertNameForFileName(TString histoName)
{
  if(histoName=="finalevtflow0") return "_0btag";
  if(histoName=="finalevtflow1") return "_1btag";
  if(histoName=="finalevtflow2") return "_2btags";
  if(histoName=="finalevtflow3") return "_3btags";
  if(histoName=="finalevtflow4") return "_4btags";
  if(histoName=="finalevtflow5") return "_geq5btags";
  return ""; 
}

// 
TString convertMassPointNameForFileName(TString title)
{

  if(title=="#splitline{H^{+}#rightarrow tb}{[180 GeV]}")     return "_m180_";
  if(title=="#splitline{H^{+}#rightarrow tb}{[200 GeV]}")     return "_m200_";
  if(title=="#splitline{H^{+}#rightarrow tb}{[220 GeV]}")     return "_m220_";
  if(title=="#splitline{H^{+}#rightarrow tb}{[240 GeV]}")     return "_m240_";
  if(title=="#splitline{H^{+}#rightarrow tb}{[250 GeV]}")     return "_m250_";
  if(title=="#splitline{H^{+}#rightarrow tb}{[260 GeV]}")     return "_m260_";
  if(title=="#splitline{H^{+}#rightarrow tb}{[280 GeV]}")     return "_m280_";
  if(title=="#splitline{H^{+}#rightarrow tb}{[300 GeV]}")     return "_m300_";
  if(title=="#splitline{H^{+}#rightarrow tb}{[300 GeV]}")     return "_m350_";
  if(title=="#splitline{H^{+}#rightarrow tb}{[300 GeV]}")     return "_m400_";
  if(title=="#splitline{H^{+}#rightarrow tb}{[300 GeV]}")     return "_m500_";
  if(title=="#splitline{H^{+}#rightarrow tb}{[300 GeV]}")     return "_m600_";
  if(title=="#splitline{H^{+}#rightarrow tb}{[300 GeV]}")     return "_m700_";



  return ""; 
}


//
float getIntegratedSystematics(TH1F *h,const std::map<TString, TH1F*> &hSysts, std::map<TString,float> &rateSysts, int bin)
{
  if(h==0) return -1;
  float rate=h->GetBinContent(bin);
  float varUp(0),varDown(0);
  for(std::map<TString, TH1F*>::const_iterator it=hSysts.begin(); it!=hSysts.end(); it++)
    {
      float var=it->second->GetBinContent(bin)/rate;
      if(isnan(var) || isinf(var)) continue;
      if(var>1) varUp += pow(var-1,2);
      else      varDown += pow(1-var,2);
    } 
  for(std::map<TString, float>::iterator it=rateSysts.begin(); it!=rateSysts.end(); it++)
    {
      varUp += pow(it->second,2);
      varDown += pow(it->second,2);
    }
  varUp=sqrt(varUp);
  varDown=sqrt(varDown);
  float var=0.5*(varUp+varDown);
  return var*rate;
}


//
void printHelp()
{
  printf("Options\n");
  printf("--out       --> output director\n");
  printf("--suffix    --> suffix to append to datacard filenames\n");
  printf("--in        --> input file from plotter\n");
  printf("--syst      --> input file with syst shapes\n");
  printf("--json      --> json file with the sample descriptor\n");
  printf("--histo     --> name of histogram to be used\n");
  printf("--noPowers --> Do not use powers of 10 for numbers in tables\n");
  printf("--bins      --> list of bins to be used (they must be comma separated without space)\n");
  printf("--ch        --> list of channels to be used (they must be comma separated without space)\n");
}

//
Shape_t getShapeFromFile(TFile* inF, TString ch, JSONWrapper::Object &Root, TFile *systF)
{
  Shape_t shape; 
  shape.totalBckg=NULL;
  shape.signal=NULL;
  shape.data=NULL;

  std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
  for(unsigned int i=0;i<Process.size();i++)
    {
      TString procCtr(""); procCtr+=i;
      TString proc=(Process[i])["tag"].toString();
      cout << "--------------------" << endl;
      cout << "Tag: " << proc << endl;  

      TDirectory *pdir = (TDirectory *)inF->Get(proc);         
      if(pdir==0) continue;
      
      bool isData(Process[i]["isdata"].toBool());
      bool isSignal(Process[i]["issignal"].toBool());
      int color(1);       if(Process[i].isTag("color" ) ) color  = (int)Process[i]["color" ].toInt();
      int lcolor(color);  if(Process[i].isTag("lcolor") ) lcolor = (int)Process[i]["lcolor"].toInt();
      int mcolor(color);  if(Process[i].isTag("mcolor") ) mcolor = (int)Process[i]["mcolor"].toInt();
      int fcolor(color);  if(Process[i].isTag("fcolor") ) fcolor = (int)Process[i]["fcolor"].toInt();
      int lwidth(1);      if(Process[i].isTag("lwidth") ) lwidth = (int)Process[i]["lwidth"].toInt();
      int lstyle(1);      if(Process[i].isTag("lstyle") ) lstyle = (int)Process[i]["lstyle"].toInt();
      int fill(1001);     if(Process[i].isTag("fill"  ) ) fill   = (int)Process[i]["fill"  ].toInt();
      int marker(20);     if(Process[i].isTag("marker") ) marker = (int)Process[i]["marker"].toInt();
  
      TH1F* syst = (TH1F*) pdir->Get("optim_systs");
      for(int ivar = 1; ivar<= (syst==0? 1 : syst->GetNbinsX());ivar++)
	{
	  TString varName = (syst==0? "" : syst->GetXaxis()->GetBinLabel(ivar));
	  if(!varName.IsNull()) systVars.insert(varName);
	  
	  TString histoName = ch;  if(!ch.IsNull()) histoName += "_"; histoName += histo+varName ;
	  TH1F* hshape = (TH1F*) pdir->Get( histoName );
	  if(hshape==0) continue;
	  cout << "Hshape nbins: " << hshape->GetXaxis()->GetNbins() << endl;
	  //project out required bins (set the others to 0)
	  if(binsToProject.size()) {
	    for(int ibin=1; ibin<=hshape->GetXaxis()->GetNbins(); ibin++) { 
	      if(find(binsToProject.begin(),binsToProject.end(),ibin) != binsToProject.end()) continue;
	      hshape->SetBinContent(ibin,0);
	      hshape->SetBinError(ibin,0);
	    }
	  }
	  else for(int ibin=1; ibin<=hshape->GetXaxis()->GetNbins(); ibin++) binsToProject.push_back(ibin); 

	  //format shape
	  utils::root::fixExtremities(hshape,true,true);
	  hshape->SetDirectory(0);  
	  hshape->SetTitle(proc);
	  hshape->SetFillColor(color); 
	  hshape->SetLineColor(lcolor); 
	  hshape->SetMarkerColor(mcolor);
	  hshape->SetFillStyle(fill);  
	  hshape->SetFillColor(fcolor);
	  hshape->SetLineWidth(lwidth); 
	  hshape->SetMarkerStyle(marker); 
	  hshape->SetLineStyle(lstyle);
	
	 //save in structure
	  if(isData){
	    if(varName=="")  shape.data=hshape;
	    else continue;
	  }else if(isSignal){
	    if(varName==""){
	      float thxsec = Process[i]["data"].daughters()[0]["xsec"].toDouble();
	      float thxsecunc=0;
	      if(Process[i]["data"].daughters()[0].isTag("xsecunc") )  thxsecunc = Process[i]["data"].daughters()[0]["xsecunc"].toDouble();
	      shape.crossSections[proc]=std::pair<float,float>(thxsec,thxsecunc);
	      shape.signal=hshape;
	      //	      shape.signalMassPoints.push_back(hshape);

	      //get alternative shapes for signal from systematics file
	      if(systF)
		{
		  TString signalVars[]={"powheg","q2up","q2down","mepsup","mepsdown"};
		  for(size_t isigvar=0; isigvar<sizeof(signalVars)/sizeof(TString); isigvar++)
		    {
		      TH1F *hmcsig=(TH1F *)systF->Get("t#bar{t}syst"+signalVars[isigvar]+"/"+histoName);
		      // TH1F *normHisto=hshape; 
		      //if(signalVars[isigvar].Contains("meps") || signalVars[isigvar].Contains("q2"))
		      //normHisto=(TH1F *)systF->Get("t#bar{t}systspincorrel/"+histoName);
		      if(hmcsig==0) { cout << "Skipping null variation: " << signalVars[isigvar] << endl;  continue; }
		      for(int ibin=1; ibin<=hshape->GetXaxis()->GetNbins(); ibin++) { 
			if(find(binsToProject.begin(),binsToProject.end(),ibin) != binsToProject.end()) continue;
			hmcsig->SetBinContent(ibin,0);
			hmcsig->SetBinError(ibin,0);
		      }
		      hmcsig->SetDirectory(0); 
		      //Double_t sf=hshape->Integral()/hmcsig->Integral();
		      hmcsig->Scale(hshape->Integral()/hmcsig->Integral()); 
		      hmcsig->SetName(hmcsig->GetName()+TString("mcsignalup"));
		      hmcsig->SetTitle(proc);

		      //check this rates from theortical point of view. in prep the variation is very high...
		      // if(signalVars[isigvar]=="q2up")           shape.rateUncs["q2_rate"]=std::pair<float,float>(sf,sf);
		      // else if(signalVars[isigvar]=="q2down")    shape.rateUncs["q2_rate"].second=sf;
		      // if(signalVars[isigvar]=="mepsup")         shape.rateUncs["meps_rate"]=std::pair<float,float>(sf,sf);
		      // else if(signalVars[isigvar]=="mepsdown")  shape.rateUncs["meps_rate"].second=sf;
		      
		      //if variation corresponds already to a signed variation save it directly
		      //otherwise create an artificial symmetric variation to build the +/- envelope of the nominal shape
		      if(signalVars[isigvar].EndsWith("up") || signalVars[isigvar].EndsWith("down"))  
			{
			  shape.signalVars[signalVars[isigvar]]=hmcsig;
			}
		      else
			{
			  shape.signalVars[signalVars[isigvar]+"up"]   = hmcsig;
			  TH1F *hmcsigdown=(TH1F *) hmcsig->Clone(hmcsig->GetName()+TString("mcsignaldown"));
			  for(int ibin=1; ibin<=hmcsigdown->GetXaxis()->GetNbins();ibin++)
			    {
			      float var=hmcsig->GetBinContent(ibin)-hshape->GetBinContent(ibin);
			      float newVal(hshape->GetBinContent(ibin)-var);
			      hmcsigdown->SetBinContent(ibin,newVal);
			    }
			  shape.signalVars[signalVars[isigvar]+"down"] = hmcsigdown;
			}
		    }
		}
            }
	    else{
	      shape.signalVars[varName]=hshape;
            }
	  }else if(proc.Contains("Z#rightarrow ll")){
	    if(varName==""){
	      shape.dataDrivenBckg.insert(proc);

	      //get the parametrized dy uncertainties
	      std::map<TString,float> dyUnc=getDYUncertainties(ch);
	      
	      //set the estimated stat uncertainty
	      for(int ibin=1; ibin<=hshape->GetXaxis()->GetNbins(); ibin++) hshape->SetBinError(ibin,hshape->GetBinContent(ibin)*dyUnc["stat"]); 
	      shape.bckg.push_back(hshape);

	      //set the other uncertainties with variations of the shape
	      for(std::map<TString,float>::iterator dyUncIt=dyUnc.begin(); dyUncIt!=dyUnc.end(); dyUncIt++)
		{
		  if(dyUncIt->first=="stat") continue; //this is set separately
		  if(dyUncIt->second<=0)     continue;
		  systVars.insert(dyUncIt->first+"up");  systVars.insert(dyUncIt->first+"down");

		  TString newName(hshape->GetName()); newName+=dyUncIt->first;
		  TH1F *dyUncUp   = (TH1F *) hshape->Clone(newName+"up");   dyUncUp->SetDirectory(0);
		  TH1F *dyUncDown = (TH1F *) hshape->Clone(newName+"down"); dyUncDown->SetDirectory(0);
		  for(int ibin=1; ibin<=hshape->GetXaxis()->GetNbins(); ibin++) 
		    {
		      dyUncUp  ->SetBinContent(ibin,std::min(2*hshape->GetBinContent(ibin), std::max(0.01*hshape->GetBinContent(ibin), hshape->GetBinContent(ibin)*(1+dyUncIt->second))));
		      dyUncDown->SetBinContent(ibin,std::min(2*hshape->GetBinContent(ibin), std::max(0.01*hshape->GetBinContent(ibin), hshape->GetBinContent(ibin)*(1-dyUncIt->second))));
		    }
		  shape.bckgVars[proc][dyUncIt->first+"up"]=dyUncUp;
		  shape.bckgVars[proc][dyUncIt->first+"down"]=dyUncDown;
		}
	    }
	  }
	  else{
	    if(varName==""){
	      shape.bckg.push_back(hshape);
	      float thxsec = Process[i]["data"].daughters()[0]["xsec"].toDouble();
	      float thxsecunc=0;
	      if(Process[i]["data"].daughters()[0].isTag("xsecunc") )  thxsecunc = Process[i]["data"].daughters()[0]["xsecunc"].toDouble();
	      shape.crossSections[proc]=std::pair<float,float>(thxsec,thxsecunc);
	    }
	    else{
	      shape.bckgVars[proc][varName]=hshape;
	    }
	  }
	}
    }

  //compute the total background
  for(size_t i=0; i<shape.bckg.size(); i++)
    {
      if(shape.totalBckg==0) { shape.totalBckg = (TH1F *)shape.bckg[i]->Clone(ch+"_"+histo+"_total"); shape.totalBckg->SetDirectory(0); shape.totalBckg->SetTitle("total bckg"); }
      else                   { shape.totalBckg->Add(shape.bckg[i]); }
    }

  //total prediction
  if(shape.totalBckg && shape.signal)
    {
      shape.totalSplusB = (TH1F *) shape.totalBckg->Clone(ch+"_"+histo+"_totalsplusb"); 
      shape.totalSplusB->SetTitle("total");
      shape.totalSplusB->SetDirectory(0);
      shape.totalSplusB->Add(shape.signal);
    }
  
  //all done
  return shape;
}


//
void getYieldsFromShapes(const map<TString, Shape_t> &allShapes)
{
  FILE* pFile = fopen(outUrl+"CrossSectionYields"+convertNameForFileName(histo)+"_"+suffixUrl+".tex","w");

  TH1F *dataTempl=allShapes.begin()->second.data;
  const std::vector<TH1F *> &bckgTempl=allShapes.begin()->second.bckg;
  for(std::vector<int>::iterator bIt = binsToProject.begin(); bIt != binsToProject.end(); bIt++)
    {

      cout << "pointer " << dataTempl << endl;
      cout << "histo" << dataTempl->GetTitle();
      cout << ", with " << dataTempl->GetNbinsX();
      cout << " bins, choosing bin " << (*bIt) << endl;
 	
      TString cat=dataTempl->GetXaxis()->GetBinLabel(*bIt);
        
      //table header
      fprintf(pFile,"\\begin{center}\n\\caption{Event yields expected for background and signal processes and observed in data for the %s category. The uncertainty associated to the limited statistics in the MC is quoted separately from the other systematic uncertainties.}\n\\label{tab:table}\n",cat.Data());
      TString Ccol   = "\\begin{tabular}{|c|";
      TString Cval   = "Channel ";
      for(std::map<TString,Shape_t>::const_iterator cIt=allShapes.begin(); cIt!=allShapes.end(); cIt++) {
	TString ch(cIt->first); if(ch.IsNull()) ch="inclusive"; 
	Ccol += "l";
	TString icol(ch/*+"-"+cat*/); utils::TLatexToTex(icol);
	Cval += " & "+icol+" ";      
      }
      Ccol += "}\\hline\\hline\n"; fprintf(pFile,"%s",Ccol.Data());
      Cval += "\\\\\\hline\n";      fprintf(pFile,"%s",Cval.Data());

      //event yields
      std::map<TString, TString > CByields;
      TString CSyields,CSpByields,CDyields;
      int ich(0);
      for(std::map<TString,Shape_t>::const_iterator cIt=allShapes.begin(); cIt!= allShapes.end(); cIt++,ich++) 
	{
	  const Shape_t &shape=cIt->second;
	  
	  //data
	  if(ich==0) CDyields = "data ";
	  CDyields += " & "; CDyields += (int) shape.data->GetBinContent(*bIt);
	  
	  float totalSyst(0);
	  
	  //signal
	  TString sigProc=shape.signal->GetTitle();
	  if(ich==0) { CSyields = sigProc; utils::TLatexToTex(CSyields); }
	  std::map<TString,float> sigRateSysts;
	  if(shape.crossSections.find(sigProc) != shape.crossSections.end())
	    {
	      std::pair<float,float> xsec=shape.crossSections.find(sigProc)->second;
	      sigRateSysts["xsec"]=xsec.second/xsec.first;
	  }
	  for(std::map<TString,std::pair<float,float> >::const_iterator rIt = shape.rateUncs.begin(); rIt!=shape.rateUncs.end(); rIt++)
	      sigRateSysts[rIt->first]=0.5*(fabs(rIt->second.first-1)+fabs(rIt->second.second-1));
	  sigRateSysts["lumi"]=lumiUnc;
	  sigRateSysts["seleff"]=selEffUnc;
	  float sSyst=getIntegratedSystematics(shape.signal,shape.signalVars,sigRateSysts,*bIt);
	  totalSyst += pow(sSyst,2); 
	  CSyields += " & ";
	  CSyields += utils::toLatexRounded( shape.signal->GetBinContent(*bIt), shape.signal->GetBinError(*bIt), sSyst, doPowers);
	  
	  //background
	  for(std::vector<TH1F *>::const_iterator bckgIt=bckgTempl.begin(); bckgIt!=bckgTempl.end(); bckgIt++)
	    {
	      TString proc=(*bckgIt)->GetTitle();
	      bool procFound(false);
	      for(std::vector<TH1F *>::const_iterator bckgItt=shape.bckg.begin(); bckgItt!=shape.bckg.end(); bckgItt++)
		{
		  if(proc!=TString((*bckgItt)->GetTitle())) continue;
		  float bSyst(-1);
		  if(true)
		    {
		      std::map<TString,float> rateSysts;
		      if(shape.dataDrivenBckg.find(proc)==shape.dataDrivenBckg.end())
			{
			  if(shape.crossSections.find(proc) != shape.crossSections.end())
			    {
			      std::pair<float,float> xsec=shape.crossSections.find(proc)->second;
			      rateSysts["xsec"]=xsec.second/xsec.first;
			    }
			  rateSysts["lumi"]=lumiUnc;
			}
		      std::map<TString, TH1F*> bckgVars;
		      if(shape.bckgVars.find(proc)!=shape.bckgVars.end()) bckgVars =shape.bckgVars.find(proc)->second;
		      bSyst=getIntegratedSystematics(*bckgItt,bckgVars,rateSysts,*bIt);
		      totalSyst += pow(bSyst,2); 
		    }
		  
		  if(ich==0) { CByields[proc]=proc; utils::TLatexToTex(CByields[proc]); }
		  CByields[proc] += " & "; CByields[proc] += utils::toLatexRounded( (*bckgItt)->GetBinContent(*bIt), (*bckgItt)->GetBinError(*bIt), bSyst, doPowers); 
		  procFound=true;
		  break;
		}
	      if(procFound) continue;
	      CByields[proc] += " & ";
	    }	  
	  
	  //signal + background
	  totalSyst=sqrt(totalSyst);
	  if(ich==0) { CSpByields = shape.totalSplusB->GetTitle(); utils::TLatexToTex(CSpByields); }
	  CSpByields += " & "; CSpByields += utils::toLatexRounded( shape.totalSplusB->GetBinContent(*bIt), shape.totalSplusB->GetBinError(*bIt),totalSyst, doPowers);
	}
      for(std::map<TString,TString>::iterator cbyIt=CByields.begin(); cbyIt!=CByields.end(); cbyIt++)
	{ 
	  cbyIt->second += "\\\\\n";  
	  fprintf(pFile,"%s",cbyIt->second.Data()); 
	}
      CSyields += "\\\\\\hline\n";       fprintf(pFile,"%s",CSyields.Data()); 
      CSpByields += "\\\\\\hline\n";     fprintf(pFile,"%s",CSpByields.Data());
      CDyields += "\\\\\\hline\n";       fprintf(pFile,"%s",CDyields.Data());      

      //close table
      fprintf(pFile,"\\hline\\end{tabular}\n\\end{center}\n");
      fprintf(pFile,"\n\n\n\n");
    }
  fclose(pFile);
}

//
void saveShapeForMeasurement(TH1F *h, TDirectory *oDir,TString syst)
{
  if(h==0 || oDir==0) return;
  oDir->cd();
  TString proc=convertNameForDataCard(h->GetTitle());
  if(syst.IsNull())
    {
      if(proc=="data") h->Write("data_obs");
      else {
	h->Write(proc);
	
	//build also the statistical uncertainty shapes
	//for each bin set content as val +/- statErr (beware however of negative and extremely large values)
	TString statSystName(proc+"_"); 
	//	if(proc=="signal") statSystName=signalTag+"_";//"ttbar_";
	statSystName+=oDir->GetTitle(); 
	statSystName+="_stat";
	TH1* statup   = (TH1 *)h->Clone(statSystName+"Up");
	TH1* statdown = (TH1 *)h->Clone(statSystName+"Down");
	for(int ibin=1; ibin<=statup->GetXaxis()->GetNbins(); ibin++){
	  statup  ->SetBinContent(ibin,std::min(2*h->GetBinContent(ibin), std::max(0.01*h->GetBinContent(ibin), statup  ->GetBinContent(ibin) + statup  ->GetBinError(ibin))));
	  statdown->SetBinContent(ibin,std::min(2*h->GetBinContent(ibin), std::max(0.01*h->GetBinContent(ibin), statdown->GetBinContent(ibin) - statdown->GetBinError(ibin))));
	}
	statup  ->Write(proc+"_"+statSystName+"Up");
	statdown->Write(proc+"_"+statSystName+"Down");
      }
    }
  else
    {
      TString systName(proc+"_");
      //systName+=oDir->GetTitle();
      //systName+="_";
      systName+=syst;
      systName.ReplaceAll("down","Down");
      systName.ReplaceAll("up","Up");
      h->Write(systName);
    }
}


//
void convertShapesToDataCards(const map<TString, Shape_t> &allShapes)
{
  
  //  for(size_t currentPoint=0; currentPoint<allShapes.signalMassPoints.size(); ++currentPoint){

    TFile *fOut = TFile::Open(outUrl+"CrossSectionShapes"+convertNameForFileName(histo)+"_"+suffixUrl+".root","RECREATE");
    for(std::map<TString, Shape_t>::const_iterator it=allShapes.begin(); it!=allShapes.end(); it++)
      {
	TString ch(it->first); if(ch.IsNull()) ch="inclusive";
	TDirectory *oDir=fOut->mkdir(ch);
	
	TString shapesFile("DataCard_"+ch+convertNameForFileName(histo)+"_"+suffixUrl+".dat");
	const Shape_t &shape=it->second;
	
	FILE* pFile = fopen(outUrl+shapesFile,"w");
	
	fprintf(pFile, "imax 1\n");
	fprintf(pFile, "jmax *\n");
	fprintf(pFile, "kmax *\n");
	fprintf(pFile, "-------------------------------\n");
	TString shapesFileName("CrossSectionShapes"+convertNameForFileName(histo)+"_"+suffixUrl+".root");
	fprintf(pFile, "shapes * * %s %s/$PROCESS %s/$PROCESS_$SYSTEMATIC\n",shapesFileName.Data(), ch.Data(), ch.Data());
	fprintf(pFile, "-------------------------------\n");
	
	//observations
	fprintf(pFile, "bin a\n");
	fprintf(pFile, "observation %f\n",shape.data->Integral());
	fprintf(pFile, "-------------------------------\n");
	saveShapeForMeasurement(shape.data,oDir);
	
	//process rows
	fprintf(pFile,"%30s ", "bin");
	fprintf(pFile,"%6s ","a");
	for(size_t j=0; j<shape.bckg.size(); j++) fprintf(pFile,"%6s ", "a");
	fprintf(pFile,"\n");
	
	fprintf(pFile,"%30s ", "process");
	fprintf(pFile,"%6s ", convertNameForDataCard(shape.signal->GetTitle()).Data());
	for(size_t j=0; j<shape.bckg.size(); j++) fprintf(pFile,"%6s ", convertNameForDataCard(shape.bckg[j]->GetTitle()).Data());
	fprintf(pFile,"\n");
	
	fprintf(pFile,"%30s ", "process");
	fprintf(pFile,"%6s ","0");
	for(size_t j=0; j<shape.bckg.size(); j++) fprintf(pFile,"%6d ", int(j+1));
	fprintf(pFile,"\n");
	
	fprintf(pFile,"%30s ", "rate");
	fprintf(pFile,"%6.3f ",shape.signal->Integral());
	saveShapeForMeasurement(shape.signal,oDir);
	for(size_t j=0; j<shape.bckg.size(); j++) { fprintf(pFile,"%6.3f ", shape.bckg[j]->Integral()); saveShapeForMeasurement(shape.bckg[j],oDir); }
	fprintf(pFile,"\n");
	fprintf(pFile, "-------------------------------\n");
	
	//systematics
	
	//lumi
	fprintf(pFile,"%30s_%dTeV %10s","lumi",int(iEcm),"lnN");
	fprintf(pFile,"%6.3f ",1+lumiUnc);
	for(size_t j=0; j<shape.bckg.size(); j++) {
	  if(shape.dataDrivenBckg.find(shape.bckg[j]->GetTitle()) != shape.dataDrivenBckg.end()) fprintf(pFile,"%6s ","-");
	  else                                                                                   fprintf(pFile,"%6.3f ",1+lumiUnc);
	}
	fprintf(pFile,"\n");
	
	//sel eff
	fprintf(pFile,"%30s_%dTeV %10s","seleff",int(iEcm),"lnN");
	fprintf(pFile,"%6.3f ",1+selEffUnc);
	for(size_t j=0; j<shape.bckg.size(); j++) {
	  if(shape.dataDrivenBckg.find(shape.bckg[j]->GetTitle()) != shape.dataDrivenBckg.end()) fprintf(pFile,"%6s ","-");
	  else                                                                                   fprintf(pFile,"%6.3f ",1+selEffUnc);
	}
	fprintf(pFile,"\n");
	
	//diepton BR
	fprintf(pFile,"%35s %10s","br","lnN");
	fprintf(pFile,"%6.3f ",1.017);
	for(size_t j=0; j<shape.bckg.size(); j++) {
	  if(convertNameForDataCard(shape.bckg[j]->GetTitle())!="ttbar") fprintf(pFile,"%6s ","-");
	  else                                                           fprintf(pFile,"%6.3f ",1.017);
	}
	fprintf(pFile,"\n");
	
	//rate systematics
	for(std::map<TString,std::pair<float,float> >::const_iterator rIt = shape.rateUncs.begin(); rIt!=shape.rateUncs.end(); rIt++)
	  {
	    fprintf(pFile,"%35s %10s",rIt->first.Data(),"lnN");
	    fprintf(pFile,"%6.3f ",1+0.5*(fabs(rIt->second.first-1)+fabs(rIt->second.second-1)));
	    for(size_t j=0; j<shape.bckg.size(); j++) {
	      fprintf(pFile,"%6s ","-");
	    }
	  }
	
	//th.xsec
	// fprintf(pFile,"%35s %10s ", "theoryUncXS_ttbar", "lnN");
	//       std::pair<float,float> ttbarXsec=shape.crossSections.find(shape.signal->GetTitle())->second;
	//       fprintf(pFile,"%6.3f ",1.0+ttbarXsec.second/ttbarXsec.first);
	//       for(size_t j=0; j<shape.bckg.size(); j++) {
	// 	if(convertNameForDataCard(shape.bckg[j]->GetTitle())!="ttbar") fprintf(pFile,"%6s ","-");
	// 	else                                                           fprintf(pFile,"%6.3f ",1.0+ttbarXsec.second/ttbarXsec.first);
	//       }
	//       fprintf(pFile,"\n");
	
	for(size_t j=0; j<shape.bckg.size(); j++)
	  {
	    TString proc(convertNameForDataCard(shape.bckg[j]->GetTitle()));
	    if(proc=="ttbar") continue;
	    std::pair<float,float> procXsec=shape.crossSections.find(shape.bckg[j]->GetTitle())->second;
	    if(procXsec.second<=0) continue;
	    if(shape.dataDrivenBckg.find(shape.bckg[j]->GetTitle()) != shape.dataDrivenBckg.end()) continue;
	    
	    TString uncName("theoryUncXS_"+proc);
	    fprintf(pFile,"%35s %10s ", uncName.Data(), "lnN");
	    fprintf(pFile,"%6s ","-");
	    for(size_t k=0; k<shape.bckg.size(); k++) {
	      if(k!=j) fprintf(pFile,"%6s ","-");
	      else if(shape.dataDrivenBckg.find(shape.bckg[k]->GetTitle()) != shape.dataDrivenBckg.end()) fprintf(pFile,"%6s ","-");
	      else                                                                                        fprintf(pFile,"%6.3f ",1.0+procXsec.second/procXsec.first);
	    }
	    fprintf(pFile,"\n");
	  }
	
	//fakes
	fprintf(pFile,"%35s %10s ", "fakes", "lnN");
	fprintf(pFile,"%6s ","-");
	for(size_t j=0; j<shape.bckg.size(); j++) {
	  TString name=convertNameForDataCard(shape.bckg[j]->GetTitle());
	  if(name!="qcd" && name!="ttbar" && name !="w")  fprintf(pFile,"%6s ","-");
	  else                                            fprintf(pFile,"%6s ","2.0");
	}
	fprintf(pFile,"\n");
	
	//systematics described by shapes
	for(std::set<TString>::iterator it=systVars.begin(); it!=systVars.end(); it++)
	  {
	    if(it->EndsWith("down")) continue;
	    TString systName(*it);
	    if(systName.EndsWith("up")) systName.Remove(systName.Length()-2,2);
	    
	    bool systIsValid(false);
	    TString systLine("");
	    char systBuf[500];
	    sprintf(systBuf,"%35s %10s ", systName.Data(), "shape");   systLine+=systBuf; memset( systBuf, 0, sizeof(systBuf) );
	    if(shape.signalVars.find(*it)==shape.signalVars.end())     { sprintf(systBuf,"%6s ","-"); systLine+=systBuf; memset( systBuf, 0, sizeof(systBuf) ); }
	    else if(shape.signalVars.find(*it)->second->Integral()==0) { sprintf(systBuf,"%6s ","-"); systLine+=systBuf; memset( systBuf, 0, sizeof(systBuf) ); }
	    else 
	      { 
		systIsValid=true;
		sprintf(systBuf,"%6s ","1"); 
		systLine+=systBuf; memset( systBuf, 0, sizeof(systBuf) );
	      saveShapeForMeasurement(shape.signalVars.find(systName+"up")->second,oDir,systName+"up"); 
	      saveShapeForMeasurement(shape.signalVars.find(systName+"down")->second,oDir,systName+"down"); 
	      }
	    for(size_t j=0; j<shape.bckg.size(); j++) {
	      if(shape.bckgVars.find(shape.bckg[j]->GetTitle())==shape.bckgVars.end())                                                                { sprintf(systBuf,"%6s ","-"); systLine+=systBuf; memset( systBuf, 0, sizeof(systBuf) ); }
	      //else if(shape.dataDrivenBckg.find(shape.bckg[j]->GetTitle()) != shape.dataDrivenBckg.end())                                           { sprintf(systBuf,"%6s ","-"); systLine+=systBuf; memset( systBuf, 0, sizeof(systBuf) ); }
	      else if(shape.bckgVars.find(shape.bckg[j]->GetTitle())->second.find(*it)==shape.bckgVars.find(shape.bckg[j]->GetTitle())->second.end()) { sprintf(systBuf,"%6s ","-"); systLine+=systBuf; memset( systBuf, 0, sizeof(systBuf) ); }
	      else if(shape.bckgVars.find(shape.bckg[j]->GetTitle())->second.find(*it)->second->Integral()==0)                                        { sprintf(systBuf,"%6s ","-"); systLine+=systBuf; memset( systBuf, 0, sizeof(systBuf) ); }
	      else                                                                                                                                  
		{ 
		  systIsValid=true;
		  sprintf(systBuf,"%6s ","1"); 
		  systLine+=systBuf; memset( systBuf, 0, sizeof(systBuf) );
		  saveShapeForMeasurement(shape.bckgVars.find(shape.bckg[j]->GetTitle())->second.find(systName+"up")->second,oDir,systName+"up"); 
		  saveShapeForMeasurement(shape.bckgVars.find(shape.bckg[j]->GetTitle())->second.find(systName+"down")->second,oDir,systName+"down"); 
		}
	    }
	    if(systIsValid) fprintf(pFile,"%s \n",systLine.Data());
	  }
	
	//MC statistics (is also systematic but written separately, it is saved at the same time as the nominal shape)
	TString myName=convertNameForDataCard(shape.signal->GetTitle()).Data();
	fprintf(pFile,"%35s %10s ", (myName/*signalTag*/+"_"+ch+"_stat").Data(), "shape");
	fprintf(pFile,"%6s ","1");
	for(size_t j=0; j<shape.bckg.size(); j++) {
	  //	if(convertNameForDataCard(shape.bckg[j]->GetTitle())!="ttbar") fprintf(pFile,"%6s ","-");
	  //	else                                                           fprintf(pFile,"%6s ","1");
	  fprintf(pFile,"%6s","-");
	}
	fprintf(pFile,"\n");
	
	for(size_t j=0; j<shape.bckg.size(); j++)
	  {
	    TString proc(convertNameForDataCard(shape.bckg[j]->GetTitle()));
	    //	  if(proc=="ttbar") continue;
	    //if(shape.dataDrivenBckg.find(shape.bckg[j]->GetTitle()) != shape.dataDrivenBckg.end()) continue;
	    
	    fprintf(pFile,"%35s %10s ", (proc+"_"+ch+"_stat").Data(), "shape");
	    fprintf(pFile,"%6s ","-");
	    for(size_t k=0; k<shape.bckg.size(); k++) {
	      if(k!=j) fprintf(pFile,"%6s ","-");
	      //else if(shape.dataDrivenBckg.find(shape.bckg[k]->GetTitle()) != shape.dataDrivenBckg.end()) fprintf(pFile,"%6s ","-");
	      else     fprintf(pFile,"%6s ","1");
	    }
	    fprintf(pFile,"\n");
	  }
	
	//all done
	fprintf(pFile,"\n");
	fclose(pFile);
      }      
    
    fOut->Close();
    
    // }
}      

//
int main(int argc, char* argv[])
{
  setTDRStyle();

  //get input arguments
  for(int i=1;i<argc;i++){
    string arg(argv[i]);
    if(arg.find("--help")        !=string::npos)              { printHelp(); return -1;} 
    else if(arg.find("--out")     !=string::npos && i+1<argc) { outUrl = argv[i+1]; outUrl+="/";  i++;  printf("out = %s\n", outUrl.Data());  }
    else if(arg.find("--suffix") !=string::npos && i+1<argc) { suffixUrl = argv[i+1];  i++;  printf("suffix = %s\n", suffixUrl.Data());  }
    else if(arg.find("--in")     !=string::npos && i+1<argc)  { inFileUrl = argv[i+1];  i++;  printf("in = %s\n", inFileUrl.Data());  }
    else if(arg.find("--syst")   !=string::npos && i+1<argc)  { systFileUrl = argv[i+1];  i++;  printf("syst = %s\n", systFileUrl.Data());  }
    else if(arg.find("--json")   !=string::npos && i+1<argc)  { jsonFileUrl  = argv[i+1];  i++;  printf("json = %s\n", jsonFileUrl.Data()); }
    else if(arg.find("--histo")  !=string::npos && i+1<argc)  { histo     = argv[i+1];  i++;  printf("histo = %s\n", histo.Data()); }
    else if(arg.find("--noPowers" )!=string::npos){ doPowers= false;    }
    else if(arg.find("--bins")   !=string::npos && i+1<argc)  { char* pch = strtok(argv[i+1],",");printf("bins to use are : ");while (pch!=NULL){int b; sscanf(pch,"%d",&b); printf(" %d ",b); binsToProject.push_back(b);  pch = strtok(NULL,",");}printf("\n"); i++; }
    else if(arg.find("--ch")     !=string::npos && i+1<argc)  { char* pch = strtok(argv[i+1],",");printf("ch to use are : ");  while (pch!=NULL){printf(" %s ",pch); channels.push_back(pch);  pch = strtok(NULL,",");}printf("\n"); i++; }
  }
  if(jsonFileUrl.IsNull() || inFileUrl.IsNull() || histo.IsNull()) { printHelp(); return -1; }
  if(channels.size()==0) { channels.push_back("ee"); channels.push_back("mumu"); channels.push_back("emu"); channels.push_back(""); }

 
  //get the shapes
  TFile *inF = TFile::Open(inFileUrl);
  JSONWrapper::Object jsonF(jsonFileUrl.Data(), true);
  TFile *systF=0;
  if(!systFileUrl.IsNull()) systF=TFile::Open(systFileUrl);
  std::map<TString, Shape_t> shapes;
  for(std::vector<string>::iterator cIt = channels.begin(); cIt != channels.end(); cIt++) shapes[*cIt]=getShapeFromFile(inF, *cIt, jsonF,systF);
  inF->Close();
  //  if(!systF) systF->Close();

  //print the tables/datacards
  getYieldsFromShapes(shapes);
  
  //convert shapes to datacards
  convertShapesToDataCards(shapes);

  //  //combine the datacards and run profile likelihood analysis
  //  TString combCardCmd("combineCards.py "),plrAnalysisCmd("runPLRanalysis --in "); 
  //  int icard(1);
  //  for(std::map<TString, Shape_t>::const_iterator it=shapes.begin(); it!=shapes.end(); it++)
  //    {
  //      TString ch(it->first); if(ch.IsNull() || ch=="inclusive") continue;
  //      combCardCmd += "Name"; combCardCmd += icard; combCardCmd +="="+outUrl+"DataCard_"+ch+".dat ";
  //      plrAnalysisCmd += outUrl+"DataCard_"+ch+".dat,";
  //      icard++;
  //    }
  //  combCardCmd += " > "+outUrl+"DataCard_combined.dat";
  //  plrAnalysisCmd += outUrl+"DataCard_combined.dat";
  //  gSystem->Exec("mv DataCard* " + outUrl);
  //  gSystem->Exec(combCardCmd.Data());
  //  gSystem->Exec(plrAnalysisCmd.Data());
  //  gSystem->Exec("mv PLR* " + outUrl);
}

