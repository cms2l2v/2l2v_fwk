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

#include "UserCode/llvv_fwk/src/tdrstyle.C"

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
      /*
	systVars.push_back("puup");
	systVars.push_back("pudown");
	systVars.push_back("jesup");
	systVars.push_back("jesdown");
	systVars.push_back("jerup");
	systVars.push_back("jerdown");
	systVars.push_back("signaltunep11");
	systVars.push_back("signaltunep11tev");
	systVars.push_back("signalmcatnlo");
      */
      systVars.push_back("signalscaleup");
      systVars.push_back("signalscaledown");
      systVars.push_back("signalmatchingup");
      systVars.push_back("signalmatchingdown");
      systVars.push_back("mcstatsup");
      systVars.push_back("mcstatsdown");
    }
  const size_t nsystVars=systVars.size();
  
  TString procs[]={"VV", "Single top", "W,multijets","Z#rightarrow ll","t#bar{t}","other t#bar{t}"};

  TString ch[]               ={"eeeq1jets",                    "ee",                    "mumueq1jets",                   "mumu",                    "emueq1jets",            "emu"};
  size_t runNsystVars[]      ={1,                              nsystVars,               1,                                nsystVars,                 1,                       nsystVars};
  TString signalRegionHisto[]={"ee_eq1jetsdilarccosine",       "ee_dilarccosine",       "mumu_eq1jetsdilarccosine",       "mumu_dilarccosine",       "emu_eq1jetsmtsum",      "emu_mtsum"};
  TString templateHisto[]    ={"ee_eq1jetslowmetdilarccosine", "ee_lowmetdilarccosine", "mumu_eq1jetslowmetdilarccosine", "mumu_lowmetdilarccosine", "emu_eq1jetsmtsum",      "emu_mtsum"};
  TString templateTitle[]    ={"Z#rightarrow ee (=1 jets)",    "Z#rightarrow ee",       "Z#rightarrow #mu#mu (=1 jets)",  "Z#rightarrow #mu#mu",     "Z#rightarrow #tau#tau (=1 jets)", "Z#rightarrow #tau#tau"};
  TString templateName[]     ={"dytoee",                       "dytoee",                "dytomumu",                       "dytomumu",                "dytoemu",               "dytoemu", "dytoemu"};
  bool rebin[]               ={false,                           false,                   false,                            false,                    false,                   false};
  const size_t nchs=sizeof(ch)/sizeof(TString);

  TH1F *dysfH = new TH1F("dysf",";Category;SF_{DY}",nchs,0,nchs);
  dysfH->SetDirectory(0);

  for(size_t ich=0; ich<nchs; ich++)
    {

      //
      // DY MODEL
      //
      TFile *dyReplacementFile=TFile::Open(stdUrl);
      TString dataDir("data");
      if(ch[ich].BeginsWith("emu")) dataDir="Z#rightarrow ll";      
      
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
      if(rebin[ich]) dyReplacementHisto->Rebin();
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
      if(rebin[ich])   dataHisto->Rebin();
           
      //get histograms for MC
      std::vector<TH1F *> dyReplacementMCHisto, dyMCHisto, otherProcsHisto, mcSumHisto;
      const size_t nprocs=sizeof(procs)/sizeof(TString);
      for(size_t ivar=0; ivar<runNsystVars[ich]; ivar++)
	{
	  TString ivarName=systVars[ivar];
	  TString ivarOrigName=ivarName;
	  if(ivarName.Contains("signal") || ivarName.Contains("mcstats")) ivarName="";
	  for(size_t iproc=0; iproc<nprocs; iproc++)
	    {
	      TH1F *histo = (TH1F *) stdFile->Get(procs[iproc]+"/"+signalRegionHisto[ich]+ivarName);
	      cout << procs[iproc] << " " << histo << endl;
	      if(histo==0) continue;
	      if(ivarOrigName.Contains("signal") && procs[iproc]=="t#bar{t}")
		{
		  cout << "[WARNING] signal will be replaced from alternative simulation for " << ivarOrigName << endl;
		  TFile *systReplacementFile=TFile::Open( dyReplacementUrl );
		  systReplacementFile->cd();
		  TH1F *repHisto=0;
		  if(ivarOrigName.Contains("powheg"))       repHisto = (TH1F *) systReplacementFile->Get("t#bar{t}systpowheg/"+signalRegionHisto[ich]);
		  if(ivarOrigName.Contains("spincorrel"))   repHisto = (TH1F *) systReplacementFile->Get("t#bar{t}systspincorrel/"+signalRegionHisto[ich]);
		  if(ivarOrigName.Contains("spinmcatnlo"))  repHisto = (TH1F *) systReplacementFile->Get("t#bar{t}systmcatnlo/"+signalRegionHisto[ich]);
		  if(ivarOrigName.Contains("matchingup"))   repHisto = (TH1F *) systReplacementFile->Get("t#bar{t}systmepsup/"+signalRegionHisto[ich]);
		  if(ivarOrigName.Contains("matchingdown")) repHisto = (TH1F *) systReplacementFile->Get("t#bar{t}systmepsdown/"+signalRegionHisto[ich]);
		  if(ivarOrigName.Contains("scaleup"))      repHisto = (TH1F *) systReplacementFile->Get("t#bar{t}systq2up/"+signalRegionHisto[ich]);
		  if(ivarOrigName.Contains("scaledown"))    repHisto = (TH1F *) systReplacementFile->Get("t#bar{t}systq2down/"+signalRegionHisto[ich]);
		  if(repHisto==0)
		      cout << "[WARNING] using default: couldn't find replacement for " << ivarOrigName << endl;
		  else
		    histo=repHisto;
		  histo->SetDirectory(0);
		  systReplacementFile->Close();
		  stdFile->cd();
		}
	      if(rebin[ich]) histo->Rebin();
	      
	      //force the normalization of the shape
	      if(ivar>0 && !procs[iproc].Contains("Z#rightarrow"))
		{
		  TH1F *nominal=(TH1F *) stdFile->Get(procs[iproc]+"/"+signalRegionHisto[ich]);
		  if(nominal==0) continue;
		  float sf=nominal->Integral()/histo->Integral();
		  if(isnan(sf) || isinf(sf)) continue;
		  histo->Scale(sf);
		}

	      //now save histogram in the appropriate category
	      if(procs[iproc].Contains("Z#rightarrow"))   
		{
		  cout << "DY" << endl;
		  histo = (TH1F *)histo->Clone(templateName[ich]+ivarOrigName+"mc");
		  histo->SetTitle(templateTitle[ich]+" (MC)");
		  histo->SetDirectory(0);		  
		  dyMCHisto.push_back( histo );
		  
		  TH1F *repHistoMC = (TH1F *) stdFile->Get(procs[iproc]+"/"+templateHisto[ich]+ivarName);
		  repHistoMC       = (TH1F *) repHistoMC->Clone(templateName[ich]+ivarOrigName+"repmc");
		  if(rebin[ich]) 
		    repHistoMC->Rebin();
		  repHistoMC->SetDirectory(0);
		  dyReplacementMCHisto.push_back( repHistoMC );
		}
	      else if(otherProcsHisto.size()==ivar)
		{
		  histo = (TH1F *) histo->Clone("otherprocs"+ivarOrigName);
		  histo->SetDirectory(0);
		  histo->SetTitle("Other processes");
 		  otherProcsHisto.push_back( histo );
		}
	      else otherProcsHisto[ivar]->Add(histo);
	      
	      //total mc
	      if(histo==0) continue;
	      if(mcSumHisto.size()<=ivar)
		{
		  histo = (TH1F *)histo->Clone("mcsum"+templateName[ich]+ivarOrigName);
		  histo->SetDirectory(0);
		  histo->SetTitle("Total MC");
		  mcSumHisto.push_back( histo );
		}
	      else mcSumHisto[ivar]->Add(histo);
	    }

	  //for MC statistics : modify bin contents according to variation
	  if(ivarOrigName.Contains("mcstat"))
	    {
	      int sign(ivarOrigName.EndsWith("up")?+1:-1);
	      for(int ibin=1; ibin<=mcSumHisto[ivar]->GetXaxis()->GetNbins(); ibin++)
		{
		  float val=max(mcSumHisto[ivar]->GetBinContent(ibin)+sign*mcSumHisto[ivar]->GetBinError(ivar),0.);
		  mcSumHisto[ivar]->SetBinContent(ibin,val);
		}
	    }
	}

      //close the file
      stdFile->Close();
      
      //
      // Now fit
      //
      cout << dataHisto << " " << dyMCHisto.size() << " " << otherProcsHisto.size() << endl;
      float totalData=dataHisto->Integral(1,dataHisto->GetXaxis()->GetNbins());
      float dyExpected=dyMCHisto[0]->Integral(1,dyMCHisto[0]->GetXaxis()->GetNbins());
      float totalOthers=otherProcsHisto[0]->Integral(1,otherProcsHisto[0]->GetXaxis()->GetNbins());
      float uncOthers=float(sqrt(totalOthers)); //max(float(0.04*totalOthers),float(sqrt(totalOthers)));  //allow 30% unc. on backgrounds

      //define variable
      RooRealVar x("x","x", dataHisto->GetXaxis()->GetXmin(), dataHisto->GetXaxis()->GetXmax());
      if(signalRegionHisto[ich]=="emu_eq1jetsmtsum") x.setRange(0,350.);

      //data
      RooDataHist* sumData = new RooDataHist("sumData", "sumData", RooArgList(x), dataHisto);
      
      //data (or MC) driven template for DY
      //protect against bins with 0 entries...
      for(int ibin=1; ibin<=dyReplacementHisto->GetXaxis()->GetNbins(); ibin++)	
	if(dyReplacementHisto->GetBinContent(ibin)<=0)
	  dyReplacementHisto->SetBinContent(ibin,0.01); 
      RooDataHist* dataTemplate = new RooDataHist("dataTemplate", "dataTemplate", RooArgList(x), dyReplacementHisto, smoothOrder );
      RooHistPdf modelDataTemplate("modelDataTemplate", "modelDataTemplate", RooArgSet(x), *dataTemplate);
      RooRealVar ndysf("SF_{DY}","dyyieldssfactor",0.91,0.5,3.0);
      RooFormulaVar ndy("N_{DY}","@0*@1",RooArgSet(RooConst(dyExpected),ndysf));

      //MC based template for other processes
      RooDataHist *otherMCTemplate = new RooDataHist("otherMCTemplate", "otherMCTemplate", RooArgList(x), otherProcsHisto[0], smoothOrder);
      RooHistPdf modelOtherMCTemplate("modelOtherMCTemplate", "modelOtherMCTemplate", RooArgSet(x), *otherMCTemplate);
      RooRealVar nother("N_{other}","otheryields",totalOthers,0,totalData);
      RooGaussian other_constraint("otherconstraintpdf","otherconstraintpdf", nother, RooConst(totalOthers), RooConst(uncOthers));
      
      //the model of the data (DY+other)xconstraints
      RooAddPdf shapeModelData("shapeModelData","signal+background", RooArgSet(modelDataTemplate,modelOtherMCTemplate),RooArgSet(ndy,nother));
      RooProdPdf constrShapeModelData("constrShapeModelData","shape*constrain",RooArgSet(shapeModelData,other_constraint));
      
      report << endl
	     << "FIT REPORT FOR " << signalRegionHisto[ich] << endl
	     << "---------------------------------------------------" << endl
	     << "[Observed] " << totalData << endl 
	     << "---------------------------------------------------" << endl
	     << "[DY expected] " <<  dyExpected << endl 
	     << "[Others expected]:" << totalOthers << " +/- " << uncOthers << endl
	     << "---------------------------------------------------" << endl;

      //fit data
      constrShapeModelData.fitTo(*sumData,Extended(kTRUE),Constrain(nother),Save(kTRUE));
      RooNLLVar *nll = (RooNLLVar*) constrShapeModelData.createNLL(*sumData,CloneData(kFALSE),Extended(kTRUE),Constrain(nother));
      if(fabs(ndysf.getError()/ndysf.getVal())>0.5) 
	report << "Fit did not converge properly most probably..." << endl;
      report << "[SF-DY] "      << ndysf.getVal() << " +/- " << ndysf.getError()  << endl
	     << "[DY fit]"      << ndy.getVal() << " +/- "  << ndy.getVal()*ndysf.getError()/ndysf.getVal() << endl
	     << "[Others fit] " << nother.getVal() << " +/- " << nother.getError() << endl
	     << "---------------------------------------------------" << endl;

      showDYFitResults(x,*sumData,constrShapeModelData,modelDataTemplate,ndysf,*nll,
       		       signalRegionHisto[ich],"CMS preliminary",dyReplacementHisto->GetXaxis()->GetTitle());
            
      //save to histogram
      dysfH->GetXaxis()->SetBinLabel(ich+1,ch[ich]);
      dysfH->SetBinContent(ich+1,ndysf.getVal());
      dysfH->SetBinError(ich+1,ndysf.getError());

      //do the closure tests on MC for different variations 
      report << "[#Delta SF-DY in closure test]" << endl; 
      float nominalMCSF(0);
      for(size_t ivar=0; ivar<runNsystVars[ich]; ivar++)
	{
	  RooDataHist* sumMC   = new RooDataHist("sumMC",   "sumMC",   RooArgList(x), mcSumHisto[0]);

	  //float totalMC=mcSumHisto[0]->Integral();
	  float ndyExp=dyMCHisto[ivar]->Integral(); 
	  float totalOthers=otherProcsHisto[ivar]->Integral();
	  float uncOthers=max(float(sqrt(totalOthers)),float(0.04*totalOthers)); 
	  	  
	  RooDataHist* mcTemplate = new RooDataHist("mcTemplate", "mcTemplate", RooArgList(x), dyReplacementMCHisto[ivar],smoothOrder );
	  RooHistPdf modelMCTemplate("modelMCTemplate", "modelMCTemplate", RooArgSet(x), *mcTemplate);
	  RooRealVar ndysf("SF_{DY}","dyyieldssfactor",1.0,0.5,3.0);
	  RooFormulaVar ndy("N_{DY}","@0*@1",RooArgSet(RooConst(ndyExp),ndysf));

	  RooDataHist *otherMCTemplate = new RooDataHist("otherMCTemplate", "otherMCTemplate", RooArgList(x), otherProcsHisto[ivar],smoothOrder);
	  RooHistPdf modelOtherMCTemplate("modelOtherMCTemplate", "modelOtherMCTemplate", RooArgSet(x), *otherMCTemplate);
	  RooRealVar nother("N_{other}","otheryields",totalOthers,0,totalData);
	  RooGaussian other_constraint("otherconstraintpdf","otherconstraintpdf", nother, RooConst(totalOthers), RooConst(uncOthers));
       	  
	  RooAddPdf shapeModelMC("shapeModelMC","signal+background",RooArgSet(modelMCTemplate,modelOtherMCTemplate),RooArgSet(ndy,nother));
	  RooProdPdf constrShapeModelMC("constrShapeModelMC","shape*constrain",RooArgSet(shapeModelMC,other_constraint)); 

	  constrShapeModelMC.fitTo(*sumMC,Extended(kTRUE),SumW2Error(kTRUE), Constrain(nother),Save(kTRUE));

	  //report the result
	  char resBuf[1000];
	  if(ivar==0) sprintf(resBuf," Nominal: %3.4f",      ndy.getVal()/dyMCHisto[0]->Integral()-nominalMCSF);
	  else        sprintf(resBuf," %s variation %3.4f %s",systVars[ivar].Data(),(ndy.getVal()/dyMCHisto[0]->Integral()-nominalMCSF)*100,"%");
	  report << resBuf << endl;

	  //save plot for nominal fit
	  if(ivar==0) {
	    nominalMCSF=ndysf.getVal();
	    RooNLLVar *nll = (RooNLLVar*) constrShapeModelMC.createNLL(*sumMC,CloneData(kFALSE),Extended(kTRUE),Constrain(nother));
	    showDYFitResults(x,*sumMC,constrShapeModelMC,modelMCTemplate,ndysf,*nll,
			     "mc"+signalRegionHisto[ich],"CMS simulation",dyReplacementMCHisto[ivar]->GetXaxis()->GetTitle());
   
	  }    
	}  
    }
  
  //all done
  cout << report.str() << endl;

  //save in file
  TFile *fOut=TFile::Open("top_dysf.root","RECREATE");
  dysfH->Write();
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
  TPad *npad = new TPad("llpad","ll", 0.6, 0.65, 0.88, 0.93);
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



