#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooConstVar.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooSimPdfBuilder.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "RooHistPdf.h"
#include "RooNLLVar.h"
#include "RooProdPdf.h"
#include "RooProfileLL.h"
#include "RooKeysPdf.h"

#include "TList.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1D.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TSystem.h"
#include "TString.h"
#include "TH1F.h"
#include "TMath.h"
#include "TFile.h"
#include "TObjArray.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "UserCode/llvv_fwk/interface/tdrstyle.h"

#include<vector>
#include<sstream>

using namespace RooFit;
using namespace std;

bool doSyst(false);
TString stdUrl(""),outUrl(""),dyReplacementUrl("");
int smoothOrder(1);

void printHelp();
void showDYFitResults(RooRealVar &x, RooDataHist &data, RooAbsPdf &model,RooAbsPdf &dyModel, RooRealVar &ndysf,RooNLLVar &nll,
		      TString tag, TString caption, TString xvar);

void fixTemplate(TH1 *h);


//
void fixTemplate(TH1 *h){
  if(h==0) return;
  for(int ibin=1; ibin<=h->GetXaxis()->GetNbins(); ibin++)	
    if(h->GetBinContent(ibin)<=0)
      h->SetBinContent(ibin,0.01); 
}

//
void printHelp()
{
  printf("--in      --> input file with standard control plots from showControlPlots\n");
  printf("--ttrep   --> input file with control plots from showControlPlots for the syst. json (optional)\n");
  printf("--out     --> output directory (optional)\n");
  printf("--syst    --> will run systematics also\n");
  printf("--smooth  --> smooth histograms before templating\n");
  printf("command line example: fitDYforTop --in std_plotter.root --ttrep syst_plotter.root\n");
}

//
int main(int argc,char *argv[])
{
  for(int i=1;i<argc;i++)
    {
      string arg(argv[i]);
      if(arg.find("--help")!=string::npos) { printHelp();    return 0; }
      if(arg.find("--in")!=string::npos)                  { stdUrl=argv[i+1];           gSystem->ExpandPathName(stdUrl);            i++;  printf("in      = %s\n", stdUrl.Data()); }
      if(arg.find("--out")!=string::npos)                 { outUrl=argv[i+1];           gSystem->ExpandPathName(outUrl);  gSystem->Exec("mkdir -p "+outUrl); i++;  printf("out     = %s\n", outUrl.Data()); }
      if(arg.find("--ttrep")!=string::npos)               { dyReplacementUrl=argv[i+1]; gSystem->ExpandPathName(dyReplacementUrl);  i++;  printf("ttRep   = %s\n", dyReplacementUrl.Data()); }
      if(arg.find("--syst")!=string::npos)                { doSyst=true;                                                            printf("Will run systematics\n"); }
      if(arg.find("--smooth")!=string::npos)              { smoothOrder=5;                                                          printf("Will smooth templates with %d-order interpol.\n",smoothOrder); }
    }
  if(stdUrl=="") { printHelp(); return 0; }

  //report to produce
  stringstream report;
  
  //global definitions
  setTDRStyle();
  RooMsgService::instance().setSilentMode(true);
  RooMsgService::instance().getStream(0).removeTopic(Minimization);
  RooMsgService::instance().getStream(1).removeTopic(Minimization);
  RooMsgService::instance().getStream(1).removeTopic(Fitting);
  RooMsgService::instance().getStream(1).removeTopic(Eval);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
  RooMsgService::instance().getStream(0).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);

  //systematics to consider for the templates
  std::vector<TString> systVars;
  systVars.push_back("");
  if(doSyst)
    {
      //instrumental variations
      systVars.push_back("jesup");
      systVars.push_back("jesdown");
      systVars.push_back("puup");
      systVars.push_back("pudown");
      systVars.push_back("jerup");
      systVars.push_back("jerdown");
      systVars.push_back("umetup");
      systVars.push_back("umetdown");
      systVars.push_back("topptup");
      systVars.push_back("topptdown");

      //theoretical variations
      //systVars.push_back("systmcatnlo");
      systVars.push_back("systq2up");
      systVars.push_back("systq2down");
      systVars.push_back("systmepsup");
      systVars.push_back("systmepsdown");
      systVars.push_back("systtunep11");
      systVars.push_back("systtunep11tev");
      systVars.push_back("systtunep11mpihi");
      systVars.push_back("systtunep11nocr");
      systVars.push_back("systpowhegpy");
      systVars.push_back("systpowheghw");
      systVars.push_back("mcstatsup");
      systVars.push_back("mcstatsdown");
    }
  const size_t nsystVars=systVars.size();
  
  TString procs[]={"VV", "Single top", "Z#rightarrow ll","t#bar{t}","other t#bar{t}", "W,multijets"};
  TString ch[]               ={"eeeq1jets",                    "ee",                    "mumueq1jets",                   "mumu",                    "emueq1jets",            "emu"};//,      "ee", "mumu"};
  size_t runNsystVars[]      ={nsystVars,                              nsystVars,               nsystVars,                                nsystVars,                 nsystVars,                       nsystVars};//, 1,1};
  TString signalRegionHisto[]={"ee_eq1jetsdilarccosine",       "ee_dilarccosine",       "mumu_eq1jetsdilarccosine",       "mumu_dilarccosine",       "emu_eq1jetsmtsum",      "emu_mtsum"};//, "ee_osbtagdilarccosine", "mumu_osbtagdilarccosine"};
  TString templateHisto[]    ={"ee_eq1jetslowmetdilarccosine", "ee_lowmetdilarccosine", "mumu_eq1jetslowmetdilarccosine", "mumu_lowmetdilarccosine", "emu_eq1jetsmtsum",      "emu_mtsum"};//, "ee_osbtagdilarccosine", "mumu_osbtagdilarccosine"};
  TString templateTitle[]    ={"Z#rightarrow ee (=1 jets)",    "Z#rightarrow ee",       "Z#rightarrow #mu#mu (=1 jets)",  "Z#rightarrow #mu#mu",     "Z#rightarrow #tau#tau (=1 jets)", "Z#rightarrow #tau#tau"};//, "Z#rightarrow ee", "Z#rightarrow #mu#mu"};
  TString templateName[]     ={"dytoee",                       "dytoee",                "dytomumu",                       "dytomumu",                "dytoemu",               "dytoemu"};// ,"dytoee", "dytomumu"};
  size_t nchs=sizeof(ch)/sizeof(TString);

  TH1F *dysfH = new TH1F("dysf",";Category;SF_{DY}",nchs,0,nchs);
  dysfH->SetDirectory(0);
  TH1F *dysfbiasH = new TH1F("dysfbias",";Category;SF_{DY} bias",nchs,0,nchs);
  dysfbiasH->SetDirectory(0);

  for(size_t ich=0; ich<nchs; ich++)
    {

      //
      // DY MODEL
      //
      TFile *dyReplacementFile=TFile::Open(stdUrl);
      TString dataDir("data");
      if(ch[ich].Contains("emu")) dataDir="Z#rightarrow ll";      
      
      //DY->tautau from mu->tau replacement
      //if(ch[ich].BeginsWith("emu") && dyReplacementUrl.IsNull()) continue;
      //TFile *dyReplacementFile=TFile::Open( (ch[ich].BeginsWith("emu")) ? dyReplacementUrl : stdUrl);
      //get the DY template (data-driven)
      //TString dataDir("data");
      //if(ch[ich].BeginsWith("emu") ) dataDir="Z#rightarrow#tau#tausystdata";
            
      TH1F *dyReplacementHisto=(TH1F *) dyReplacementFile->Get(dataDir+"/"+templateHisto[ich]);
      dyReplacementHisto->SetDirectory(0); 
      dyReplacementHisto->SetName(templateName[ich]+"model");
      dyReplacementHisto->SetTitle(templateTitle[ich]);
      dyReplacementFile->Close();

      //
      // DATA AND OTHER PROCESSES
      //
     
      //open standard file
      TFile *stdFile=TFile::Open(stdUrl);
      
      //get histogram for data
      TH1F *dataHisto=(TH1F *) stdFile->Get("data/"+signalRegionHisto[ich])->Clone("data");
      dataHisto->SetDirectory(0);
      dataHisto->SetTitle("data");
           
      //get histograms for MC
      TH1F *dyMCHisto=0;
      std::vector<TH1F *> otherProcsHisto;
      const size_t nprocs=sizeof(procs)/sizeof(TString);
      for(size_t iproc=0; iproc<nprocs; iproc++)
	{
	  
	  for(size_t ivar=0; ivar<runNsystVars[ich]; ivar++)
	    {

	      //baseline shape
	      TH1F *histo = (TH1F *) stdFile->Get(procs[iproc]+"/"+signalRegionHisto[ich]);
	      if(histo==0) continue;
	      histo->SetDirectory(0);

	      //variations (only worth considering for ttbar)
	      if(procs[iproc]=="t#bar{t}" && ivar>0)
		{
		  cout << "[WARNING] signal will be replaced from alternative simulation for " << systVars[ivar] << endl;

		  TH1F *repHisto=0;
		  if(systVars[ivar].Contains("syst"))
		    {
		      TFile *systReplacementFile=TFile::Open( dyReplacementUrl );
		      systReplacementFile->cd();
		      repHisto=(TH1F *) systReplacementFile->Get("t#bar{t}"+systVars[ivar]+"/"+signalRegionHisto[ich]);
		      if(repHisto) repHisto->SetDirectory(0);
		      systReplacementFile->Close();
		      stdFile->cd();
		    }
		  else
		    {
		      repHisto=(TH1F *) stdFile->Get("t#bar{t}/"+signalRegionHisto[ich]+systVars[ivar]);
		      if(repHisto) repHisto->SetDirectory(0);
		    }
		  
		  if(repHisto==0) cout << "[WARNING] using default: couldn't find replacement for " << systVars[ivar] << " using default" << endl;
		  else
		    {
		      repHisto->Scale(histo->Integral()/repHisto->Integral());
		      histo=repHisto;
		    }
		}
	      
	      //now save histogram in the appropriate category
	      if(procs[iproc].Contains("Z#rightarrow"))
		{
		  if(ivar==0)   
		    {
		      dyMCHisto = (TH1F *)histo->Clone(templateName[ich]+"mc");
		      dyMCHisto->SetTitle(templateTitle[ich]+" (MC)");
		      dyMCHisto->SetDirectory(0);		  
		    }
		}
	      else 
		{
		  if(otherProcsHisto.size()==ivar)
		    {
		      histo = (TH1F *) histo->Clone("otherprocs"+systVars[ivar]);
		      histo->SetDirectory(0);
		      histo->SetTitle("Other processes");
		      otherProcsHisto.push_back( histo );
		    }
		  else otherProcsHisto[ivar]->Add(histo);
		}
	    }
	}
      //close the file
      stdFile->Close();
      
      //
      // Now fit
      //
      Double_t totalData   ( dataHisto->Integral(1,dataHisto->GetXaxis()->GetNbins()) );
      Double_t dyExpected  ( dyMCHisto->Integral(1,dyMCHisto->GetXaxis()->GetNbins()) );
      Double_t uncOthers(0);
      Double_t totalOthers ( otherProcsHisto[0]->IntegralAndError(1,otherProcsHisto[0]->GetXaxis()->GetNbins(),uncOthers) );
      uncOthers = sqrt( pow(uncOthers,2)+pow(0.044*totalOthers,2)+pow(0.1*totalOthers,2) ); //MC stats + Lumi + 10% on th. xsec

      //define variable
      RooRealVar x("x",dataHisto->GetXaxis()->GetTitle(), dataHisto->GetXaxis()->GetXmin(), dataHisto->GetXaxis()->GetXmax());
      if(ch[ich].Contains("emu")) x.setRange(0,400.);
      
      //data
      RooDataHist* sumData = new RooDataHist("sumData", "sumData", RooArgList(x), dataHisto);
      
      //data (or MC) driven template for DY
      fixTemplate(dyReplacementHisto);
      RooDataHist* dataTemplate = new RooDataHist("dataTemplate", "dataTemplate", RooArgList(x), dyReplacementHisto, smoothOrder );
      RooHistPdf modelDataTemplate("modelDataTemplate", "modelDataTemplate", RooArgSet(x), *dataTemplate);
      RooRealVar ndysf("SF_{DY}","dyyieldssfactor",0.91,0.5,3.0);
      RooFormulaVar ndy("N_{DY}","@0*@1",RooArgSet(RooConst(dyExpected),ndysf));

      //MC based template for other processes
      fixTemplate(otherProcsHisto[0]);
      RooDataHist *otherTemplate = new RooDataHist("otherTemplate", "otherTemplate", RooArgList(x), otherProcsHisto[0], smoothOrder);
      RooHistPdf modelOtherTemplate("modelOtherTemplate", "modelOtherTemplate", RooArgSet(x), *otherTemplate);
      RooRealVar nother("N_{other}","otheryields",totalOthers,0,totalData);
      RooGaussian other_constraint("otherconstraintpdf","otherconstraintpdf", nother, RooConst(totalOthers), RooConst(uncOthers));
      
      //the model of the data (DY+other)xconstraints
      RooAddPdf shapeModelData("shapeModelData","signal+background", RooArgSet(modelDataTemplate,modelOtherTemplate),RooArgSet(ndy,nother));
      RooProdPdf constrShapeModelData("constrShapeModelData","shape*constrain",RooArgSet(shapeModelData,other_constraint));
      
      //fit data
      constrShapeModelData.fitTo(*sumData,Extended(kTRUE),Constrain(nother),Save(kTRUE));
      RooNLLVar *nll = (RooNLLVar*) constrShapeModelData.createNLL(*sumData,CloneData(kFALSE),Extended(kTRUE),Constrain(nother));
      showDYFitResults(x,*sumData,constrShapeModelData,modelDataTemplate,ndysf,*nll,
       		       signalRegionHisto[ich],"CMS preliminary",dyReplacementHisto->GetXaxis()->GetTitle());
      

      //report the result
      report << endl
	     << "FIT REPORT FOR " << signalRegionHisto[ich] << endl
	     << "---------------------------------------------------" << endl
	     << "[Observed] " << totalData << endl 
	     << "---------------------------------------------------" << endl
	     << "[DY expected] " <<  dyExpected << endl 
	     << "[Others expected]:" << totalOthers << " +/- " << uncOthers << endl
	     << "---------------------------------------------------" << endl;
      if(fabs(ndysf.getError()/ndysf.getVal())>0.5) 
	report << "Fit did not converge properly most probably..." << endl;
      report << "[SF-DY] "      << ndysf.getVal() << " +/- " << ndysf.getError()  << endl
	     << "[DY fit]"      << ndy.getVal() << " +/- "  << ndy.getVal()*ndysf.getError()/ndysf.getVal() << endl
	     << "[Others fit] " << nother.getVal() << " +/- " << nother.getError() << endl
	     << "---------------------------------------------------" << endl;
      
      //save to histogram
      dysfH->GetXaxis()->SetBinLabel(ich+1,ch[ich]);
      if(signalRegionHisto[ich].Contains("osbtag"))  dysfH->GetXaxis()->SetBinLabel(ich+1,ch[ich]+"osbtag");
      dysfH->SetBinContent(ich+1,ndysf.getVal());
      dysfH->SetBinError(ich+1,ndysf.getError());
	
      //
      // MC closure tests + systematic variations
      //
      report << "[#Delta SF-DY in closure test]" << endl; 


      TH1 *nominalMCHisto=(TH1 *) otherProcsHisto[0]->Clone("total"+ch[ich]+"mc");
      nominalMCHisto->Add(dyMCHisto);
      nominalMCHisto->Scale( totalData/nominalMCHisto->Integral(1,nominalMCHisto->GetXaxis()->GetNbins()) );
      RooDataHist* sumMC   = new RooDataHist("sumMC",   "sumMC",   RooArgList(x), nominalMCHisto);

      Double_t nominalMCSF(1.0);
      for(size_t ivar=0; ivar<runNsystVars[ich]; ivar++)
	{
	  //for MC statistics : modify bin contents according to variation
	  if(systVars[ivar].Contains("mcstat"))
	    {
	      int sign(systVars[ivar].EndsWith("up")?+1:-1);
	      for(int ibin=1; ibin<=otherProcsHisto[ivar]->GetXaxis()->GetNbins(); ibin++)
		{
		  float val=max(otherProcsHisto[ivar]->GetBinContent(ibin)+sign*otherProcsHisto[ivar]->GetBinError(ivar),0.);
		  otherProcsHisto[ivar]->SetBinContent(ibin,val);
		}	 
	    }
	  otherProcsHisto[ivar]->Scale( totalOthers / otherProcsHisto[ivar]->Integral(1,otherProcsHisto[ivar]->GetXaxis()->GetNbins()) );
	  fixTemplate(otherProcsHisto[ivar]);

	  //replace the "others" component in the model
	  RooDataHist *otherMCTemplate = new RooDataHist("otherMCTemplate", "otherMCTemplate", RooArgList(x), otherProcsHisto[ivar],smoothOrder);
	  RooHistPdf modelOtherMCTemplate("modelOtherMCTemplate", "modelOtherMCTemplate", RooArgSet(x), *otherMCTemplate);
	  RooAddPdf shapeModelMC("shapeModelMC","signal+background",RooArgSet(modelDataTemplate,modelOtherMCTemplate),RooArgSet(ndy,nother));
	  RooProdPdf constrShapeModelMC("constrShapeModelMC","shape*constrain",RooArgSet(shapeModelMC,other_constraint)); 

	  //fit
	  ndysf.setVal(1.0);
	  nother.setVal(totalOthers);
	  constrShapeModelMC.fitTo(*sumMC,Extended(kTRUE),SumW2Error(kTRUE), Constrain(nother),Save(kTRUE));

	  if(ivar==0) {
	    nominalMCSF=ndysf.getVal();
	    //save to bias histogram
	    dysfbiasH->GetXaxis()->SetBinLabel(ich+1,ch[ich]);
	    if(signalRegionHisto[ich].Contains("osbtag"))  dysfbiasH->GetXaxis()->SetBinLabel(ich+1,ch[ich]+"osbtag");
	    dysfbiasH->SetBinContent(ich+1,nominalMCSF-1.0);
	    dysfbiasH->SetBinError(ich+1,0);
	  }
	  Double_t relDifference(nominalMCSF);
	  relDifference=(ndysf.getVal()/nominalMCSF-1)*100;
	  Double_t absDifference(ndysf.getVal()-nominalMCSF);
	 
	  //report the result
	  char resBuf[1000];
	  //if(ivar==0) sprintf(resBuf," Nominal: %3.4f",      relDifference);
	  //else        sprintf(resBuf," %s variation %3.4f",systVars[ivar].Data(),relDifference);
	  if(ivar==0) sprintf(resBuf," Nominal: %3.4f",      absDifference);
	  else        sprintf(resBuf," %s variation %3.4f",systVars[ivar].Data(),absDifference);
	  report << resBuf << endl;

	  //save plot for nominal fit
	  RooNLLVar *nll = (RooNLLVar*) constrShapeModelMC.createNLL(*sumMC,CloneData(kFALSE),Extended(kTRUE),Constrain(nother));
	  showDYFitResults(x,*sumMC,constrShapeModelMC,modelDataTemplate,ndysf,*nll,
			   "mc"+signalRegionHisto[ich]+systVars[ivar],"CMS simulation",dyReplacementHisto->GetXaxis()->GetTitle());
	}  
    }
  
  //all done
  cout << report.str() << endl;

  //save in file
  TFile *fOut=TFile::Open("top_dysf.root","RECREATE");
  dysfH->Write();
  dysfbiasH->Write();
  fOut->Close();
}


//
void showDYFitResults(RooRealVar &x,
		      RooDataHist &data,
		      RooAbsPdf &model,
		      RooAbsPdf &dyModel,
		      RooRealVar &ndysf,
		      RooNLLVar &nll,
		      TString tag, 
		      TString caption, 
		      TString xvar)
{
  TCanvas *cnv = new TCanvas(tag+"c",tag+"c",600,600);
  cnv->cd();
  RooPlot *frame = x.frame();
  if(tag.Contains("mc")) data.plotOn(frame,Name("data"),DataError(RooAbsData::SumW2));
  else                   data.plotOn(frame,Name("data"));
  frame->GetXaxis()->SetTitle( xvar );     frame->GetXaxis()->SetTitleOffset(1.0);
  frame->GetYaxis()->SetTitle( "Events" ); frame->GetYaxis()->SetTitleOffset(1.2);

  model.plotOn(frame,Name("total"),MoveToBack());
  model.plotOn(frame,Components(dyModel),FillStyle(1001),FillColor(kGray),LineColor(kGray+1),DrawOption("lf"),Name("dytemplate"),MoveToBack());

  frame->Draw();
  Double_t ymax=frame->GetMaximum();
  frame->GetYaxis()->SetRangeUser(0.01,max(float(ymax*1.4),float(1.0)));

  TLegend *leg=new TLegend(0.2,0.82,0.55,0.85);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetNColumns(3);
  leg->AddEntry("data","data","p");
  leg->AddEntry("total","Total","l");
  leg->AddEntry("dytemplate","DY template","f");
  leg->Draw("same");

  TPaveText *pave = new TPaveText(0.2,0.85,0.5,0.92,"NDC");
  pave->SetBorderSize(0);
  pave->SetTextAlign(12);
  pave->SetFillStyle(0);
  pave->SetTextFont(42);
  pave->AddText(caption);
  char buf[200]; 
  if(tag.BeginsWith("mc")) sprintf(buf,"SF_{DY}=%3.2f",ndysf.getVal());
  else                     sprintf(buf,"SF_{DY}=%3.2f #pm %3.2f",ndysf.getVal(),ndysf.getError());
  pave->AddText(buf);
  pave->Draw("same");
	  
  //display likelihood as inset
  TPad *npad = new TPad("llpad","ll", 0.7, 0.68, 0.95, 0.93);
  npad->Draw();
  npad->cd();
  frame = ndysf.frame();
  nll.plotOn(frame,ShiftToZero(),Name("ll"));
  frame->GetXaxis()->SetTitle("SF_{DY}");  
  frame->GetXaxis()->SetTitleOffset(1.0);
  frame->GetXaxis()->SetTitleSize(0.1);
  frame->GetXaxis()->SetLabelSize(0.08);
  frame->GetYaxis()->SetTitle("-log(L/L_{max})");
  frame->GetYaxis()->SetTitleOffset(1.2);
  frame->GetYaxis()->SetTitleSize(0.1);
  frame->GetYaxis()->SetLabelSize(0.08);
  frame->GetYaxis()->SetRangeUser(0,5);
  frame->Draw();
	      
  //update and save
  cnv->Modified();
  cnv->Update();
  cnv->SaveAs(outUrl+"/"+tag+".png");
  cnv->SaveAs(outUrl+"/"+tag+".C");
  cnv->SaveAs(outUrl+"/"+tag+".pdf");
}



