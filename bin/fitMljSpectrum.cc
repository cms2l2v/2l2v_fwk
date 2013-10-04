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
#include "UserCode/llvv_fwk/src/JSONWrapper.cc"

#include<vector>
#include<sstream>


using namespace std;
using namespace RooFit;


stringstream report; 
int iEcm               = 8;
float dataLumi         = 19701;
float baseRelUnc       = 0.027;

class Shape_t
{
public:
  Shape_t() : data(0), ddBkg(0), mcBkg(0), ddMCbkg(0), signal(0), pseudoData(0) { }
  ~Shape_t() { };

  //void finalize();
  //void runFit(TH1 *data, TH1 *correct, TH1 *model, TH1 *misassigned, bool isData, TString tag);

  TH1F *data, *ddBkg, *mcBkg, *ddMCbkg, *signal, *pseudoData;
  std::map<TString, TH1F *> signalSysts;

  RooFitResult *fitRes;
  Double_t ncorrectMC,ncorrectMC_err,ncorrectData,ncorrectData_err;
  Double_t nwrongMC,nwrongMC_err,nwrongData,nwrongData_err;
  Double_t fcorrectMC,fcorrectMC_err,fcorrectData,fcorrectData_err;
  Double_t sfcorrect, sfcorrect_err;
  Double_t rho;
};

std::map<TString, Shape_t> allShapes;
std::vector<TString> allChannels;


//
void printHelp();
void fixTemplate(TH1 *h);
void fillMljData(TString url,JSONWrapper::Object &Root,bool isSyst);

//TCanvas *plot(TObjArray *mc, TH1 *data, TObjArray *spimpose=0, bool doLog=false, bool norm=false, bool doDiff=false);

//                                                                                                                                                                                          
void printHelp()
{
  printf("--help     --> print this\n");
  printf("--in       --> input directory with the summary trees\n");
  printf("--json     --> json file with the description of the samples\n");
  printf("--systJson --> json file with the description of the samples used for systematic uncertainties\n");
  printf("--iLumi    --> total integrated luminosity\n");
  printf("--use      --> use previously created report\n");
}

//
void fixTemplate(TH1 *h)
{
  if(h==0) return;
  //do not let the templates to have 0 in any bin as it will prevent the binned fit to converge properly
  for(int ibin=1; ibin<=h->GetXaxis()->GetNbins(); ibin++)
    {
      if(h->GetBinContent(ibin)>0) continue;
      h->SetBinContent(ibin,0.00001/h->GetBinWidth(ibin));
    }
}

//
void fillMljData(TString url,JSONWrapper::Object &Root,bool isSyst)
{
  TFile *inF=TFile::Open(url);
  if(inF==0 || inF->IsZombie()) return;
  
  std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
  for(unsigned int i=0;i<Process.size();i++)
    {
      TString proc=(Process[i])["tag"].toString();
      TDirectory *pdir = (TDirectory *)inF->Get(proc);
      bool isData(Process[i]["isdata"].toBool());
      
      //iterate over the channels
      for(std::vector<TString>::iterator it=allChannels.begin(); it!=allChannels.end(); it++)
	{
	  TH1F *nomH  =(TH1F *)pdir->Get(*it+"_mlj");
	  TH1F *sigH  =(TH1F *)pdir->Get(*it+"_correctmlj"); fixTemplate(sigH);
	  TH1F *mcBkgH=(TH1F *)pdir->Get(*it+"_wrongmlj");   fixTemplate(mcBkgH);
	  TH1F *ddBkgH=(TH1F *)pdir->Get(*it+"_rotmlj");     fixTemplate(ddBkgH);
	  
	  //create a new shape if non existing
	  if(allShapes.find(*it)==allShapes.end()) { 
	    Shape_t newShape; 
	    newShape.data    = (TH1F *) nomH->Clone(*it+"_data");         newShape.data->Reset("ICE");    newShape.data->SetDirectory(0);
	    newShape.signal  = (TH1F *) sigH->Clone(*it+"_signal");       newShape.signal->Reset("ICE");  newShape.signal->SetDirectory(0);
	    newShape.ddBkg   = (TH1F *) ddBkgH->Clone(*it+"_ddbkg");      newShape.ddBkg->Reset("ICE");   newShape.ddBkg->SetDirectory(0);
	    newShape.mcBkg   = (TH1F *) mcBkgH->Clone(*it+"_mcbkg");      newShape.mcBkg->Reset("ICE");   newShape.mcBkg->SetDirectory(0);
	    newShape.ddMCbkg = (TH1F *) ddBkgH->Clone(*it+"_ddmcbkg");    newShape.ddMCbkg->Reset("ICE"); newShape.ddMCbkg->SetDirectory(0);
	    allShapes[*it]=newShape; 
	  }

	  //update shape
	  Shape_t &shape=allShapes.find(*it)->second;
	  if(isData)      { shape.data->Add(nomH); shape.ddBkg->Add(ddBkgH); continue; }
	  else if(isSyst) { shape.signalSysts[proc]=(TH1F *)nomH->Clone(*it+"_"+proc); shape.signalSysts[proc]->SetDirectory(0); continue; }
	  else {
	    shape.signal->Add(sigH);
	    shape.mcBkg->Add( mcBkgH );
	    shape.ddMCbkg->Add( ddBkgH );
	    TH1F *systH = (TH1F*) pdir->Get("optim_systs");
	    if(systH==0) continue;
            for(int ivar=1; ivar<systH->GetNbinsX(); ivar++)
	      {
		TString varName( systH->GetXaxis()->GetBinLabel(ivar) );
		TH1F *sigVarH  =(TH1F *)pdir->Get(*it+"_correctmlj"+varName);
		if(sigVarH==0) continue; 
		sigVarH->SetDirectory(0);
		shape.signalSysts[varName] = sigVarH;
	      }
	  }
	}
    }
}

//
int main(int argc, char* argv[])
{
  // load framework libraries                                                                                          
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();
  
  //configure
  TString url(""), json(""), systUrl(""), systJson("");
  for(int i=1;i<argc;i++)
    {
      string arg(argv[i]);
      if(arg.find("--help")!=string::npos)                  { printHelp();    return 0; }
      if(arg.find("--nom")!=string::npos && i+1<argc)       { url=argv[i+1];     gSystem->ExpandPathName(url);     i++; printf("url       = %s\n", url.Data()); }
      if(arg.find("--json")!=string::npos && i+1<argc)      { json=argv[i+1];                                      i++; printf("json      = %s\n", json.Data()); }
      if(arg.find("--syst")!=string::npos && i+1<argc)      { systUrl=argv[i+1]; gSystem->ExpandPathName(systUrl); i++; printf("systUrl   = %s\n", systUrl.Data()); }
      if(arg.find("--systJson")!=string::npos && i+1<argc)  { systJson=argv[i+1];                                  i++; printf("systJson  = %s\n", systJson.Data()); }
      if(arg.find("--iLumi")!=string::npos && i+1<argc)     { sscanf(argv[i+1],"%f",&dataLumi);                    i++; printf("lumi      = %f\n", dataLumi); }
    }
  if(url=="" || json=="") { printHelp(); return 0;}
  
  //general plotting style
  setTDRStyle();
//   gStyle->SetPadTopMargin   (0.06);
//   gStyle->SetPadBottomMargin(0.12);
//   gStyle->SetPadRightMargin (0.16);
//   gStyle->SetPadLeftMargin  (0.14);
//   gStyle->SetTitleSize(0.04, "XYZ");
//   gStyle->SetTitleXOffset(1.1);
//   gStyle->SetTitleYOffset(1.45);
//   gStyle->SetPalette(1);
//   gStyle->SetNdivisions(505);
  
  //categories to be fit
  allChannels.push_back("ee");      
  allChannels.push_back("eeeq2jets");
  allChannels.push_back("eeeq3jets");
  allChannels.push_back("eeeq4jets");
  allChannels.push_back("mumu");      
  allChannels.push_back("mumueq2jets");
  allChannels.push_back("mumueq3jets");
  allChannels.push_back("mumueq4jets");
  allChannels.push_back("emu");      
  allChannels.push_back("emueq2jets");
  allChannels.push_back("emueq3jets");
  allChannels.push_back("emueq4jets");

  //readout the shapes from the files
  JSONWrapper::Object jsonF(json.Data(), true);
  fillMljData(url,jsonF,false);
  if(systJson!="") { JSONWrapper::Object systJsonF(systJson.Data(), true); fillMljData(systUrl,systJsonF,true); }

  //fit all categories
  for(size_t icat=0; icat<allChannels.size(); icat++)
    {
      Shape_t &ishape=allShapes[ allChannels[icat] ];
      cout << "Analyzing [" << allChannels[icat] << "]" << endl
	   << ishape.data  << " "<< ishape.ddBkg << " " << ishape.ddMCbkg << " " << ishape.signal << endl
	   << "with " << ishape.signalSysts.size() << " signal systematics" << endl;
    }
      

  //show report
}



		 
//
// MljFitResults_t runFit(TH1 *data, TH1 *correct, TH1 *model, TH1 *misassigned, bool isData, TString tag)
// {
//   MljFitResults_t r;

//   int njets(2);
//   if(tag.Contains("eq3jets")) njets=3;
//   if(tag.Contains("eq4jets")) njets=4;
  
//   //predictions for the fraction of correct assignments
//   r.ncorrectMC=correct->IntegralAndError(1,correct->GetXaxis()->GetNbins(),r.ncorrectMC_err);
//   r.nwrongMC=misassigned->IntegralAndError(1,misassigned->GetXaxis()->GetNbins(),r.nwrongMC_err);
//   r.fcorrectMC=r.ncorrectMC/(r.ncorrectMC+r.nwrongMC);
//   r.fcorrectMC_err=sqrt(pow(r.ncorrectMC_err*(2*r.ncorrectMC+r.nwrongMC),2)+pow(r.ncorrectMC*r.nwrongMC_err,2))/pow(r.ncorrectMC+r.nwrongMC,2);
  
//   RooRealVar mlj          (tag+"_mlj",         "mlj",         data->GetXaxis()->GetXmin(), data->GetXaxis()->GetXmax());
//   RooDataHist dataHist    (tag+"_datahist",    "datahist",    mlj, RooFit::Import(*data));
//   RooDataHist correctHist (tag+"_correcthist", "correcthist", mlj, RooFit::Import(*correct));
//   RooHistPdf correctPdf   (tag+"_correctpdf",  "correctpdf",  RooArgSet(mlj), correctHist);
//   RooDataHist modelHist   (tag+"_modelhist",   "modelhist",   mlj, RooFit::Import(*model));
//   RooHistPdf modelPdf     (tag+"_modelpdf",    "modelpdf",    RooArgSet(mlj), modelHist);
//   RooRealVar ntwq         (tag+"_ntwq",        "ntwq",       r.ncorrectMC/(2*njets), 0, dataHist.sumEntries() );
//   // RooRealVar ncorrect    (tag+"_ncorrect",    "ncorrect",    r.ncorrectMC,  r.ncorrectMC*0.5, r.ncorrectMC*2);
//   RooFormulaVar ncorrect  (tag+"_ncorrect",    "@0*@1", RooArgSet(ntwq,RooConst(2*njets)));
//   RooRealVar nother       (tag+"_nother",        "nother",       r.nwrongMC/(2*njets),0, dataHist.sumEntries() );
//   //  RooRealVar nwrong       (tag+"_nwrong",      "nwrong",      r.nwrongMC,    r.nwrongMC*0.5,   r.nwrongMC*2);
//   RooFormulaVar nwrong    (tag+"_nwrong",    "@0*@1", RooArgSet(nother,RooConst(2*njets)));
//   RooFormulaVar fcorrect  (tag+"_fcorrect","@0/(@0+@1)",RooArgSet(ncorrect,nwrong));
//   RooFormulaVar fwrong    (tag+"_fwrong","1-@0",RooArgSet(fcorrect));
//   RooAddPdf     shapeModel("shapemodel",   "shapemodel", RooArgSet(correctPdf,modelPdf), RooArgSet(ncorrect,nwrong));
//   RooFitResult* fitRes = shapeModel.fitTo(dataHist,RooFit::Save(kTRUE), Extended(kTRUE), RooFit::Optimize(kTRUE), RooFit::Minos(kTRUE), RooFit::InitialHesse(kTRUE), SumW2Error(kFALSE) );
  
//   r.fitRes           = fitRes;
//   r.ncorrectData     = ncorrect.getVal();
//   r.ncorrectData_err = (ntwq.getError()/ntwq.getVal())*ncorrect.getVal();
//   r.nwrongData       = nwrong.getVal();
//   r.nwrongData_err   = (nother.getError()/nother.getVal())*nwrong.getVal();
//   r.rho              = fitRes->correlation(ntwq,nother);
//   r.fcorrectData     = ntwq.getVal()/(ntwq.getVal()+nother.getVal());
//   r.fcorrectData_err = sqrt( 
// 			    pow(nother.getVal()*ntwq.getError(),2) 
// 			    + pow( ntwq.getVal()*nother.getError(),2) 
// 			    + 2* ntwq.getVal() * nother.getVal() * r.rho * ntwq.getError() * nother.getError() 
// 			    ) 
//     / pow( ntwq.getVal()+nother.getVal(),4);
//   r.sfcorrect        = r.fcorrectData/r.fcorrectMC;
//   r.sfcorrect_err    = sqrt(pow(r.fcorrectData_err*r.fcorrectMC,2)+pow(r.fcorrectData*r.fcorrectMC_err,2))/pow(r.fcorrectMC,2);

//   //save fit values
//   float ntwqVal(ntwq.getVal()),     ntwqVal_err(ntwq.getError());
//   float notherVal(nother.getVal()), notherVal_err(nother.getError());

//   if(isData)
//     {
//       TCanvas *c=new TCanvas("c","c",600,600);
//       c->cd();

//       RooNLLVar *nll = (RooNLLVar *) shapeModel.createNLL(dataHist);

//       RooMinuit minuit(*nll); 
//       minuit.setStrategy(2);
//       minuit.setPrintLevel(1);
//       minuit.setErrorLevel(0.5);
//       minuit.hesse();
//       minuit.migrad();

//       RooPlot *cplot = minuit.contour(ntwq,nother,1,2,3);
//       cplot->Draw();
//       cplot->GetXaxis()->SetTitle("N(t#rightarrow Wq)");
//       cplot->GetXaxis()->SetTitleOffset(1.1);
//       cplot->GetXaxis()->SetRangeUser(ntwqVal-4*ntwqVal_err,ntwqVal+4*ntwqVal_err);
//       cplot->GetYaxis()->SetTitle("N(background)");
//       cplot->GetYaxis()->SetTitleOffset(1.9);
//       cplot->GetYaxis()->SetRangeUser(notherVal-4*notherVal_err,notherVal+4*notherVal_err);

//       TLegend *leg=new TLegend(0.845,0.2,0.99,0.99,"", "brNDC");
//       leg->SetFillColor(0);
//       leg->SetFillStyle(0);
//       leg->SetTextAlign(12);
//       leg->SetLineColor(0);
//       leg->SetBorderSize(0);
//       leg->SetTextFont(42);
//        TIter next(c->GetListOfPrimitives());
//       TObject *obj=0;
//       while((obj = next()))
// 	{
// 	  TString name=obj->GetName();
// 	  TString title("");
// 	  int fill=0, color=1;
// 	  if(name.Contains("n1")) { title="1-#sigma"; fill=1001; color=kAzure+9; }
// 	  if(name.Contains("n2")) { title="2-#sigma"; fill=3001; color=kOrange-2; }
// 	  if(name.Contains("n3")) { title="3-#sigma"; fill=3002; color=kGreen+3; }
// 	  if(fill==0) continue;
// 	  TGraph *gr=(TGraph *)obj;
// 	  gr->SetFillStyle(fill);
// 	  gr->SetFillColor(color);
// 	  gr->SetLineColor(1);
// 	  gr->SetDrawOption("f");
// 	  leg->AddEntry(gr,title,"f");
// 	}

//       TString header(tag);
//       header.ReplaceAll("mu","#mu");
//       header.ReplaceAll("eq"," =");
//       header.ReplaceAll("jets"," jets");
//       leg->SetHeader(header);
//       leg->Draw();

//       TPaveText *T = new TPaveText(0.1,0.995,0.84,0.95, "NDC");
//       T->SetFillColor(0); T->SetFillStyle(0);  T->SetLineColor(0); T->SetBorderSize(0);   T->SetTextAlign(22);
//       char Buffer[1024];
//       sprintf(Buffer, "CMS preliminary, #sqrt{s}=%d TeV, #scale[0.5]{#int} L=%.1f fb^{-1}, N(t#rightarrow Wq)=", iEcm, dataLumi/1000);
//       TString buf=toLatexRounded(ntwqVal,ntwqVal_err);
//       buf.ReplaceAll("\\pm","#pm");
//       buf.ReplaceAll("$","");
//       T->AddText(Buffer+buf);
//       T->Draw("same");

//       c->SaveAs(tag+"_contour_mlj.png");
//       c->SaveAs(tag+"_contour_mlj.pdf");
//       c->SaveAs(tag+"_contour_mlj.C");
//       delete c;
//     }

//   return r;
// }

//
// void fitMljData(TString reportUrl)
// {
 
//   report << "\\begin{center}"
// 	 << "\\caption{Results of the fit for the fraction of correct $t\\rightarrow Wq$ assignments in the $M_{lj}$ spectrum}" << endl 
// 	 << "\\begin{tabular}{lccc} \\hline" << endl
// 	 << "Category & $f_{correct}^{MC}$ & $f_{correct}^{data}$ & $SF_{correct}$ \\\\\\hline\\hline" << endl;

//   TString ch[]={"ee","emu","mumu"};
//   TString cats[]={"eq2jets","eq3jets","eq4jets"};
//   std::map<TString,std::vector<std::pair<TString,Float_t> > > systMap;
//   std::map<TString, TH1F *> correctH,misassignedH, dataH, dataResH,unmatchedH,wrongH,mcmodelH;
//   for(size_t ich=0; ich<sizeof(ch)/sizeof(TString); ich++)
//     {
//       for(size_t icat=0; icat<sizeof(cats)/sizeof(TString); icat++)
// 	{
// 	  TString tag(ch[ich]+cats[icat]);
	  
// 	  //get histos from file
// 	  TFile *inF=TFile::Open(reportUrl);

// 	  TH1F *data        = (TH1F *)inF->Get("data"+tag+"_mlj");	  
// 	  formatPlot(data,        1,1,1,20,0,true,true,1,0,1);
// 	  fixTemplate(data);
// 	  if(data) { data->SetTitle("data"); data->Rebin(); }
	  
// 	  TH1F *correct     = (TH1F *)inF->Get(tag+"_correctmlj");        
// 	  formatPlot(correct,     614,1,1,1,1001,true,true,1,614,614); 
// 	  fixTemplate(correct);
// 	  if(correct)     { correct->Scale(dataLumi); correct->SetTitle("correct"); correct->Rebin(); }
	  
// 	  TH1F *unmatched   = (TH1F *)inF->Get(tag+"_unmatchedmlj");      
// 	  formatPlot(unmatched,   824,1,1,1,1001,true,true,1,824,824); 
// 	  fixTemplate(unmatched);
// 	  if(unmatched)   { unmatched->Scale(dataLumi); unmatched->SetTitle("unmatched"); unmatched->Rebin(); }
	  
// 	  TH1F *misassigned = (TH1F *)inF->Get(tag+"_misassignedmlj");    
// 	  formatPlot(misassigned, 822,1,1,1,1001,true,true,1,822,822); 
// 	  fixTemplate(misassigned);
// 	  if(misassigned) { misassigned->Scale(dataLumi); misassigned->SetTitle("misassigned"); misassigned->Rebin();}
	  
// 	  TH1F *dy          = (TH1F *)inF->Get("dy"+tag+"_mlj");          
// 	  formatPlot(dy,          1,1,1,20,0,true,true,1,1,1);    
// 	  fixTemplate(dy);
// 	  if(dy)          { dy->Scale(dataLumi);      dy->SetTitle("#splitline{Z#rightarrow ll}{(MC)}"); dy->Rebin(); }

// 	  TString ztag(ch[ich]+"z"+cats[icat]);
// 	  TH1F *dymodel     = (TH1F *)inF->Get("data"+ztag+"_mlj");          
// 	  formatPlot(dymodel,          831,1,1,20,1001,true,true,1,831,831);    
// 	  fixTemplate(dymodel);
// 	  if(dymodel)      {  dymodel->SetTitle("#splitline{Z#rightarrow ll}{(data)}"); dymodel->Rebin(); }

// 	  TH1F *dymcmodel     = (TH1F *)inF->Get("dy"+ztag+"_mlj");          
// 	  formatPlot(dymcmodel,          831,1,1,1,1001,true,true,1,831,831);    
// 	  fixTemplate(dymcmodel);
// 	  if(dymcmodel)      {  dy->Scale(dataLumi); dymcmodel->SetTitle("#splitline{Z#rightarrow ll}{(pseudo-data)}"); dymcmodel->Rebin();}

// 	  TH1F *wrong       = (TH1F *)inF->Get(tag+"_wrongmlj");          
// 	  formatPlot(wrong,       592,1,1,1,1001,true,true,1,592,592); 
// 	  fixTemplate(wrong);
// 	  if(wrong)       { wrong->Scale(dataLumi);   wrong->SetTitle("wrong"); wrong->Rebin(); }

// 	  TH1F *model       = (TH1F *)inF->Get("data"+tag+"_rotmlj");     
// 	  // TH1F *model       = (TH1F *)inF->Get("data"+tag+"_swapmlj");     
// 	  formatPlot(model,       1,1,1,20,0,true,true,1,0,1);
// 	  fixTemplate(model);
// 	  if(model)       { model->SetTitle("#splitline{wrong}{(data)}"); model->Rebin(); }

// 	  TH1F *mcmodel       = (TH1F *)inF->Get(tag+"_rotmlj");     
// 	  formatPlot(mcmodel,       1,1,1,20,0,true,true,1,0,1);
// 	  fixTemplate(mcmodel);
// 	  if(mcmodel)       { 
// 	    mcmodel->SetTitle("#splitline{wrong}{(pseudo-data)}"); 
// 	    mcmodel->Rebin();
// 	    if(model) mcmodel->Scale(model->Integral()/mcmodel->Integral());
// 	  }

// 	  //correct low mass bins with mc truth prediction
// 	  int binNorm = model->FindBin(200);
// 	  float sf=model->Integral(binNorm,model->GetXaxis()->GetNbins())/misassigned->Integral(binNorm,misassigned->GetXaxis()->GetNbins());
// 	  float mcsf=mcmodel->Integral(binNorm,mcmodel->GetXaxis()->GetNbins())/misassigned->Integral(binNorm,misassigned->GetXaxis()->GetNbins());
// 	  for(int ibin=1; ibin<=model->FindBin(50); ibin++)
// 	    {
// 	      model->SetBinContent(ibin,misassigned->GetBinContent(ibin)*sf);
// 	      model->SetBinError(ibin,misassigned->GetBinError(ibin)*sf);
// 	      mcmodel->SetBinContent(ibin,misassigned->GetBinContent(ibin)*mcsf);
// 	      mcmodel->SetBinError(ibin,misassigned->GetBinError(ibin)*mcsf);
// 	    }

// 	  //closure test
// 	  TH1F *totalMC=(TH1F *) correct->Clone("totalmc");
// 	  totalMC->Add(misassigned);
// 	  MljFitResults_t mcr=runFit(totalMC,correct,mcmodel,misassigned,false,tag);
// 	  delete totalMC;

// 	  //fit to data
// 	  MljFitResults_t r=runFit(data,correct,model,misassigned,true,tag);
// 	  report << ch[ich] << " " << cats[icat] << " & "
// 	    //<< toLatexRounded(mcr.sfcorrect,  mcr.sfcorrect_err) << " & "
// 		 << toLatexRounded(r.fcorrectMC,   r.fcorrectMC_err) << " & "
// 		 << toLatexRounded(r.fcorrectData, r.fcorrectData_err) << " & "
// 		 << toLatexRounded(r.sfcorrect,    r.sfcorrect_err) << " \\\\" << endl;

// 	  //evaluate the systematics
// 	  typedef std::pair<TString,TString> SystKey_t;
// 	  std::map<TString, SystKey_t> systList;
// 	  systList["JES"]    = SystKey_t("jesup","jesdown");
// 	  systList["JER"]    = SystKey_t("jerup","jerdown");
// 	  systList["Pileup"] = SystKey_t("puup","pudown");
// 	  systList["ME-PS"]  = SystKey_t("mepsup","mepsdown");
// 	  systList["Q^2"]    = SystKey_t("q2up","q2down");
// 	  systList["Signal"] = SystKey_t("powheg","mcatnlo");
// 	  systList["DY"]     = SystKey_t("dyup","dydown");
// 	  for(std::map<TString, SystKey_t>::iterator sIt = systList.begin(); sIt!=systList.end(); sIt++)
// 	    {
	      
// 	      TH1F *correctUp     = correct;
// 	      TH1F *correctDown   = correct;
// 	      TH1F *dataUp       = data;
// 	      TH1F *dataDown       = data;
// 	      if(sIt->first!="DY")
// 		{
// 		  correctUp = (TH1F *)inF->Get(tag+"_correctmlj_"+sIt->second.first);        
// 		  if(correctUp==0) continue;
// 		  correctUp->SetDirectory(0);
// 		  fixTemplate(correctUp);
// 		  correctUp->Scale(correct->Integral()/correctUp->Integral());
// 		  correctUp->Rebin(); 
		  
// 		  correctDown = (TH1F *)inF->Get(tag+"_correctmlj_"+sIt->second.second);
// 		  if(correctDown==0) continue;
// 		  correctDown->SetDirectory(0);
// 		  fixTemplate(correctUp);
// 		  correctDown->Scale(correct->Integral()/correctDown->Integral());
// 		  correctDown->Rebin(); 
// 		}
// 	      else
// 		{
// 		  dataUp=(TH1F *) data->Clone("dataup");     dataUp->Add(dy,0.10);
// 		  dataDown=(TH1F *) data->Clone("datadown"); dataDown->Add(dy,-0.10);
// 		}

// 	      MljFitResults_t rUp   = runFit(dataUp,correctUp,model,misassigned,false,tag);
// 	      MljFitResults_t rDown = runFit(dataDown,correctDown,model,misassigned,false,tag);
	      
// 	      float varUp   = 100*(rUp.sfcorrect/r.sfcorrect-1);
// 	      float varDown = 100*(rDown.sfcorrect/r.sfcorrect-1);
// 	      float var=0.5*(fabs(varUp)+fabs(varDown));
// 	      if(sIt->first=="Signal") var=fabs(rUp.sfcorrect/r.sfcorrect-1); //var=fabs(rUp.sfcorrect/rDown.sfcorrect-1);
// 	      systMap[tag].push_back(std::pair<TString,Float_t>(sIt->first,var));
// 	    }
// 	  inF->Close();

// 	  //show distributions

// 	  //before the fit
// 	  TObjArray stackList; stackList.Add(misassigned);  stackList.Add(correct); 
// 	  TCanvas *c=plot(&stackList,data,0,false);
// 	  c->SaveAs(tag+"_pre_mlj.png");
// 	  c->SaveAs(tag+"_pre_mlj.pdf");
// 	  c->SaveAs(tag+"_pre_mlj.C");
	  
// 	  //postfit
// 	  TH1F *correctFit=(TH1F *)correct->Clone(correct->GetName()+TString("post"));  correctFit->Scale(r.sfcorrect);
// 	  TH1F *misassignedFit=(TH1F *)model->Clone(model->GetName()+TString("post"));  misassignedFit->Scale((misassigned->Integral()/model->Integral())*( 1-r.fcorrectData)/(1-r.fcorrectMC));
// 	  stackList.Clear(); stackList.Add(misassignedFit); stackList.Add(correctFit);
// 	  c=plot(&stackList,data,0,false);
// 	  c->SaveAs(tag+"_post_mlj.png");
// 	  c->SaveAs(tag+"_post_mlj.pdf");
// 	  c->SaveAs(tag+"_post_mlj.C");

// 	  //postfit subtracted
// 	  TH1F *dataRes=(TH1F *)data->Clone("datares"); 
// 	  dataRes->SetTitle("#splitline{data}{(residuals)}");
// 	  dataRes->Add(misassignedFit,-1);
// 	  stackList.Clear(); stackList.Add(correctFit);
// 	  c=plot(&stackList,dataRes,0,false,false,true);
// 	  c->SaveAs(tag+"_postres_mlj.png");
// 	  c->SaveAs(tag+"_postres_mlj.pdf");
// 	  c->SaveAs(tag+"_postres_mlj.C");
	  
// 	  //misassignment composition in MC
// 	  stackList.Clear();    stackList.Add(unmatched);   stackList.Add(wrong);
// 	  c=plot(&stackList,mcmodel,0,false,true);
// 	  c->SaveAs(tag+"_wrong_mlj.png");
// 	  c->SaveAs(tag+"_wrong_mlj.pdf");
// 	  c->SaveAs(tag+"_wrong_mlj.C");

// 	  stackList.Clear();      stackList.Add(wrong);
// 	  c=plot(&stackList,model,0,false,true);
// 	  c->SaveAs(tag+"_wrongmodel_mlj.png");
// 	  c->SaveAs(tag+"_wrongmodel_mlj.pdf");
// 	  c->SaveAs(tag+"_wrongmodel_mlj.C");

// 	  c=plot(&stackList,dy,0,false,true);
// 	  c->SaveAs(tag+"_wrongvsdy_mlj.png");
// 	  c->SaveAs(tag+"_wrongvsdy_mlj.pdf");
// 	  c->SaveAs(tag+"_wrongvsdy_mlj.C");

// 	  //save inclusive histos
// 	  TString key("");
// 	  if(tag.Contains("emu")) key="of";
// 	  else                    key="sf";
// 	  if(correctH.find(key)==correctH.end())
// 	    {
// 	      correctH[key]=correctFit;
// 	      misassignedH[key]=misassignedFit;
// 	      dataH[key]=data;
// 	      dataResH[key]=dataRes;
// 	      unmatchedH[key]=unmatched;
// 	      wrongH[key]=wrong;
// 	      mcmodelH[key]=mcmodel;
// 	    }
// 	  else
// 	    {
// 	      correctH[key]->Add(correctFit);
// 	      misassignedH[key]->Add(misassignedFit);
// 	      dataH[key]->Add(data);
// 	      dataResH[key]->Add(dataRes);
// 	      unmatchedH[key]->Add(unmatched);
// 	      wrongH[key]->Add(wrong);
// 	      mcmodelH[key]->Add(mcmodel);
// 	    }
	  
// 	}
//     }

//   report << "\\hline\\end{tabular}\n\\end{center}" << endl;


//   //show the inclusive plots
//   for(std::map<TString,TH1F *>::iterator it = correctH.begin(); it != correctH.end(); it++)
//     {
//       TString key=it->first;
//       TObjArray stackList;
      
//       //postfit
//       stackList.Add(misassignedH[key]); stackList.Add(correctH[key]);
//       TCanvas *c=plot(&stackList,dataH[key],0,false);
//       c->SaveAs(key+"_post_mlj.png");
//       c->SaveAs(key+"_post_mlj.pdf");
//       c->SaveAs(key+"_post_mlj.C");
      
//       //postfit subtracted
//       stackList.Clear(); stackList.Add(correctH[key]);
//       c=plot(&stackList,dataResH[key],0,false,false,true);
//       c->SaveAs(key+"_postres_mlj.png");
//       c->SaveAs(key+"_postres_mlj.pdf");
//       c->SaveAs(key+"_postres_mlj.C");
	  
//       //misassignment composition in MC
//       stackList.Clear();    stackList.Add(unmatchedH[key]);   stackList.Add(wrongH[key]);
//       c=plot(&stackList,mcmodelH[key],0,false,true);
//       c->SaveAs(key+"_wrong_mlj.png");
//       c->SaveAs(key+"_wrong_mlj.pdf");
//       c->SaveAs(key+"_wrong_mlj.C");
//     }

//   //systs table
//   report << endl;
//   for(std::map<TString,std::vector<std::pair<TString,Float_t> > >::iterator cIt=systMap.begin(); cIt!= systMap.end(); cIt++)
//     {
//       if(cIt==systMap.begin())
// 	{
// 	  report << "Source";
// 	  for(std::vector<std::pair<TString,Float_t> >::iterator sIt=cIt->second.begin(); sIt!=cIt->second.end(); sIt++)
// 	    report << " & " << sIt->first;
// 	  report << " & Total \\\\\\hline" << endl;
// 	}

//       float total(0);
//       report << cIt->first;
//       for(std::vector<std::pair<TString,Float_t> >::iterator sIt=cIt->second.begin(); sIt!=cIt->second.end(); sIt++)
// 	{
// 	  report << " & " << sIt->second;
// 	  total += pow(sIt->second,2);
// 	}
//       report << " & " << sqrt(total) << " \\\\" << endl;
//     }
// }

// //
// TCanvas *plot(TObjArray *mc, TH1 *data, TObjArray *spimpose, bool doLog, bool norm, bool doDiff)
// {
//   TCanvas* c = new TCanvas("c","c",800,800);
//   if(data==0) return c;

//   //main plot
//   c->cd();
//   TPad* t1 = new TPad("t1","t1", 0.0, 0.20, 1.0, 1.0);
//   t1->Draw();  t1->cd();
//   if(doLog) t1->SetLogy(true);

//   TLegend* leg  = new TLegend(0.845,0.2,0.99,0.99,"", "brNDC"); 
//   leg->AddEntry(data,data->GetTitle(), "P");

//   //normalization factor
//   float normFactor(1.0);
//   float totalPred(0);
//   for(int i=0; i<mc->GetEntriesFast(); i++) totalPred += ((TH1 *) mc->At(i))->Integral();
//   if(norm && totalPred>0) normFactor=data->Integral()/totalPred;
  
//   //prediction
//   THStack* stack = new THStack("MC","MC");
//   TH1 *totalMC=0;
//   for(int i=0; i<mc->GetEntriesFast(); i++)
//     {
//       TH1 *hist=(TH1 *)mc->At(i);
//       if(hist==0) continue;
//       hist=(TH1 *) hist->Clone();
//       hist->Scale(normFactor);
//       stack->Add(hist, "HIST");               
//       leg->AddEntry(hist, hist->GetTitle(), "F");

//       if(totalMC==0) { totalMC=(TH1 *) hist->Clone("totalpred"); totalMC->SetDirectory(0); }
//       else           { totalMC->Add(hist); }
//     }

//   for(int ibin=1; ibin<=totalMC->GetXaxis()->GetNbins(); ibin++)
//     {
//       Double_t error=sqrt(pow(totalMC->GetBinError(ibin),2)+pow(totalMC->GetBinContent(ibin)*baseRelUnc,2));
//       totalMC->SetBinError(ibin,error);
//     }
//   totalMC->SetFillStyle(3427);
//   totalMC->SetFillColor(kGray+1);
//   totalMC->SetMarkerStyle(1);

//   stack->Draw("");
//   stack->GetXaxis()->SetTitle(totalMC->GetXaxis()->GetTitle());
//   stack->GetYaxis()->SetTitle(totalMC->GetYaxis()->GetTitle());
//   if(norm) stack->GetYaxis()->SetTitle(totalMC->GetYaxis()->GetTitle()+ TString(" (a.u.)"));
//   stack->SetMinimum(data->GetMinimum()<0? data->GetMinimum()*1.50: 1e-2);
//   stack->SetMaximum(max(totalMC->GetMaximum()*1.10,data->GetMaximum()*1.10));

//   //overlay total uncertainty band
//   totalMC->SetLineColor(kGray);
//   totalMC->SetFillStyle(3001);
//   totalMC->SetFillColor(kGray);
//   totalMC->SetMarkerColor(kGray);
//   totalMC->SetMarkerStyle(1);
//   totalMC->Draw("e2same");

//   //for comparison
//   if(spimpose)
//     {
//       for(int i=0; i<spimpose->GetEntriesFast(); i++)
// 	{
// 	  TH1 *hist=(TH1 *)spimpose->At(i);
// 	  if(hist==0) continue;
// 	  hist = (TH1 *)hist->Clone();
// 	  hist->Scale(normFactor);
// 	  hist->Draw("histsame");
// 	  leg->AddEntry(hist, hist->GetTitle(), "L");
// 	}
//     }

//   ///data
//   data->Draw("e1 same");

//   //caption
//   TString tag=data->GetName();
//   TString header("ee");
//   if(tag.Contains("mumu"))    header="#mu#mu";
//   if(tag.Contains("emu"))     header="e#mu";
//   if(tag.Contains("eq2jets")) header+=",=2 jets";
//   if(tag.Contains("eq3jets")) header+=",=3 jets";
//   if(tag.Contains("eq4jets")) header+=",=4 jets";
//   leg->SetHeader(header);
//   leg->SetFillColor(0); leg->SetFillStyle(0); leg->SetLineColor(0);  leg->SetTextFont(42); leg->Draw("same");

//   //chi-2,K-S tests
//   TPaveText *pave = new TPaveText(0.6,0.8,0.8,0.9,"NDC");
//   pave->SetBorderSize(0);
//   pave->SetFillStyle(0);
//   pave->SetTextAlign(32);
//   pave->SetTextFont(42);
//   char buf[100];
//   sprintf(buf,"#chi^{2}/ndof : %3.2f", data->Chi2Test(totalMC,"WWCHI2/NDF") );
//   pave->AddText(buf);
//   sprintf(buf,"K-S prob : %3.2f", data->KolmogorovTest(totalMC,"X") );
//   pave->AddText(buf);
//   pave->Draw();

//   //title
//   TPaveText* T = new TPaveText(0.1,0.995,0.84,0.95, "NDC");
//   T->SetFillColor(0); T->SetFillStyle(0);  T->SetLineColor(0);  T->SetTextAlign(22);
//   char Buffer[1024]; 
//   sprintf(Buffer, "CMS preliminary, #sqrt{s}=%d TeV, #scale[0.5]{#int} L=%.1f fb^{-1}", iEcm, dataLumi/1000);
//   T->AddText(Buffer);
//   T->Draw("same");
//   T->SetBorderSize(0);
  
//   //comparison 
//   c->cd();
//   TPad* t2 = new TPad("t2","t2", 0.0, 0.0, 1.0, 0.2);
//   t2->Draw();
//   t2->cd();
//   t2->SetGridy(true);
//   t2->SetPad(0,0.0,1.0,0.2);
//   t2->SetTopMargin(0);
//   t2->SetBottomMargin(0.5);

//   //mc stats
//   TH1 *denRelUncH=(TH1 *) totalMC->Clone("totalpredrelunc");
//   for(int xbin=1; xbin<=denRelUncH->GetXaxis()->GetNbins(); xbin++)
//     {
//       if(denRelUncH->GetBinContent(xbin)==0) continue;
//       Double_t val=1.0;
//       Double_t err=denRelUncH->GetBinError(xbin)/denRelUncH->GetBinContent(xbin);
//       if(doDiff) 
// 	{
// 	  val=0.0;
// 	  err=denRelUncH->GetBinError(xbin)*1.044;
// 	}
//       denRelUncH->SetBinContent(xbin,val);
//       denRelUncH->SetBinError(xbin,err);
//     }
//   TGraphErrors *denRelUnc=new TGraphErrors(denRelUncH);
//   denRelUnc->SetLineColor(1);
//   denRelUnc->SetFillStyle(3001);
//   denRelUnc->SetFillColor(kGray);
//   denRelUnc->SetMarkerColor(1);
//   denRelUnc->SetMarkerStyle(1);
//   denRelUncH->Reset("ICE");       
//   denRelUncH->Draw();
//   denRelUnc->Draw("3");

//   TH1 *ratio=(TH1 *)data->Clone("totalpredrel");
//   if(doDiff)  ratio->Add(totalMC,-1);
//   else        ratio->Divide(totalMC);
//   ratio->SetMarkerStyle(20);
//   ratio->Draw("e1 same");

//   float yscale = (1.0-0.2)/(0.18-0);       
//   if(doDiff)
//     {
//       denRelUncH->GetYaxis()->SetTitle("Data-Exp.");
//       denRelUncH->SetMinimum(-100);
//       denRelUncH->SetMaximum(100);
//     }
//   else
//     {
//       denRelUncH->GetYaxis()->SetTitle("Data/Exp.");
//       denRelUncH->SetMinimum(0.4);
//       denRelUncH->SetMaximum(1.6);
//     }
//   denRelUncH->GetXaxis()->SetTitle("");
//   denRelUncH->GetXaxis()->SetTitleOffset(1.3);
//   denRelUncH->GetXaxis()->SetLabelSize(0.033*yscale);
//   denRelUncH->GetXaxis()->SetTitleSize(0.036*yscale);
//   denRelUncH->GetXaxis()->SetTickLength(0.03*yscale);
//   denRelUncH->GetYaxis()->SetTitleOffset(0.3);
//   denRelUncH->GetYaxis()->SetTitleFont(42);
//   denRelUncH->GetYaxis()->SetNdivisions(5);
//   denRelUncH->GetYaxis()->SetLabelSize(0.033*yscale);
//   denRelUncH->GetYaxis()->SetTitleSize(0.036*yscale);
   
//   //all done here
//   c->Modified();
//   c->Update();
//   return c;
// }
 





//   //calibrate method
//   /*
//   typedef std::pair<Double_t,Int_t> PseudoExperimentKey_t;
//   std::vector<PseudoExperimentKey_t > lumiScenarios;
//   //   lumiScenarios.push_back(PseudoExperimentKey_t(36,100));
//   //   lumiScenarios.push_back(PseudoExperimentKey_t(50,100));
//   //   lumiScenarios.push_back(PseudoExperimentKey_t(100,100));
//   //   lumiScenarios.push_back(PseudoExperimentKey_t(200,100));
//   //   lumiScenarios.push_back(PseudoExperimentKey_t(400,50));
//   //   lumiScenarios.push_back(PseudoExperimentKey_t(600,25));
//   //   lumiScenarios.push_back(PseudoExperimentKey_t(1000,10));
//   lumiScenarios.push_back(PseudoExperimentKey_t(dataLumi,5));
//   for(size_t ijet=0; ijet<njbins; ijet++)
//     {
//       //init the misassignment measurement handler from the mlj spectrum
//       exclusiveMisMeasurements[jetBins[ijet]] = new MisassignmentMeasurement("${CMSSW_BASE}/src/CMGTools/HtoZZ2l2nu/data/GR_R_42_V23_AK5PFchs_Uncertainty.txt");
//       for(int icat=0; icat<ncats; icat++)  misMeasurement->setBiasCorrections(cats[icat],0.0);

//       for(size_t ipe=0; ipe<lumiScenarios.size(); ipe++)
// 	{
// 	  PseudoExperimentKey_t pe=lumiScenarios[ipe];
// 	  bool freezeResults=(pe.first==200);
// 	  runCalibration(url,pe.first,pe.second,jetBins[ijet],freezeResults );
// 	}
//     }
//   */

	  
// // 	  TString iname(cats[icat]+jetBins[ijet]);
// // 	  TGraphErrors *gr = new TGraphErrors;
// // 	  gr->SetName(iname);
// // 	  std::pair<int,TString> key(jetBins[ijet],cats[icat]);
// // 	  calibCurves[key]=gr;
// // 	}




// /*


//   typedef std::pair<Double_t,Int_t> PseudoExperimentKey_t;
//   std::vector<PseudoExperimentKey_t > lumiScenarios;
//   if(syst=="")// || syst=="metcut") 
//     {
//       lumiScenarios.push_back(PseudoExperimentKey_t(36,100));
//       lumiScenarios.push_back(PseudoExperimentKey_t(50,100));
//       lumiScenarios.push_back(PseudoExperimentKey_t(100,100));
//       lumiScenarios.push_back(PseudoExperimentKey_t(200,100));
//       lumiScenarios.push_back(PseudoExperimentKey_t(400,50));
//       lumiScenarios.push_back(PseudoExperimentKey_t(600,25));
//       lumiScenarios.push_back(PseudoExperimentKey_t(1000,10));
//       lumiScenarios.push_back(PseudoExperimentKey_t(dataLumi,5));
//     }
//   else
//     {
//       lumiScenarios.push_back(PseudoExperimentKey_t(200,50));
//     }



//   //display the results
//   setStyle();
//   gStyle->SetOptFit(0);
//   TCanvas *biasevolc = new TCanvas("biasevolc","biasevolc",true);
//   biasevolc->SetCanvasSize(1200,1200);
//   biasevolc->SetWindowSize(1200,1200);
//   biasevolc->Divide(2,2);

//   report << "Jet bin &" << flush;
//   for(int icat=0; icat<4; icat++){
//     report << catTits[icat] << flush;
//     if(icat<3) report << " & " << flush; 
//   };
//   report << "\\\\\\hline" << endl;

//   for(size_t ijet=0; ijet<3; ijet++)
//     {
//       report << "$" << jetBinTitles[ijet] << "$ &" << flush;
//       for(int icat=0; icat<4; icat++)
// 	{
// 	  TPad *p=(TPad *) biasevolc->cd(icat+1);
// 	  std::pair<int,TString> key(jetBins[ijet],cats[icat]);
// 	  TGraphErrors *gr = calibCurves[key];
// 	  gr->SetTitle(jetBinTitles[ijet]);
// 	  gr->SetMarkerStyle(20+ijet*2);
// 	  gr->SetFillStyle(0);
// 	  gr->SetFillColor(0);
// 	  gr->Draw(ijet==0? "ap" : "p");
// 	  if(ijet==0)
// 	    {
// 	      gr->GetXaxis()->SetTitle("Integrated luminosity (pb^{-1})");
// 	      gr->GetYaxis()->SetTitle("Average bias");
// 	      gr->GetYaxis()->SetRangeUser(-0.25,0.25);
// 	      p->SetLogx(true);

// 	      TPaveText *pave = new TPaveText(0.20,0.8,0.35,0.9,"NDC");
// 	      pave->SetBorderSize(0);
// 	      pave->SetFillStyle(0);
// 	      pave->SetTextFont(42);
// 	      pave->AddText(catTits[icat]);
// 	      pave->Draw();
// 	    }
// 	  if(icat==0 && ijet==2)
// 	    {
// 	      TLegend *leg = p->BuildLegend();
// 	      formatForCmsPublic(p,leg,"CMS simulation",3);
// 	    }

// 	  //fit the bias using high lumi pseudo-experiments
// 	  gr->Fit("pol0","FMQR0+","",100,dataLumi);
// 	  Float_t avgBias=gr->GetFunction("pol0")->GetParameter(0);
// 	  Float_t avgBiasError=gr->GetFunction("pol0")->GetParError(0);
// 	  char buf[100];
// 	  sprintf(buf,"%3.3f $\\pm$ %3.3f", avgBias, avgBiasError);
// 	  report << buf;
// 	  if(icat<3) report << " & " << flush; 
// 	}

//       report << "\\\\" << endl;
//     }
//   biasevolc->Modified();
//   biasevolc->Update();
//   biasevolc->SaveAs("MisassignmentBiasEvol.C");
//   biasevolc->SaveAs("MisassignmentBiasEvol.png");
//   biasevolc->SaveAs("MisassignmentBiasEvol.pdf");

//   report << endl;
//   for(size_t ijet=0; ijet<3; ijet++)
//     {
//       report << "$" << jetBinTitles[ijet] << "$ &" << flush;
//       for(int icat=0; icat<4; icat++)
// 	{
// 	  report << catTits[icat] << flush;
// 	  if(icat<3) report << " & " << flush; 
// 	}
//       report << "\\\\\\hline" << endl;


//       report << "$f_{correct}^{MC}(signal) $ &" << flush;
//       for(int icat=0; icat<4; icat++)
// 	{
// 	  std::pair<int,TString> key(jetBins[ijet],cats[icat]);
// 	  char buf[100];
// 	  sprintf(buf,"%3.3f $\\pm$ %3.3f", fCorrectSigOnlyTrue[key][0],fCorrectSigOnlyTrue[key][1]);
// 	  report << buf;
// 	  if(icat<3) report << " & " << flush; 	  
// 	}
//       report << "\\\\" << endl;

//       report << "$f_{correct}^{MC}$ &" << flush;
//       for(int icat=0; icat<4; icat++)
// 	{
// 	  std::pair<int,TString> key(jetBins[ijet],cats[icat]);
// 	  char buf[100];
// 	  sprintf(buf,"%3.3f $\\pm$ %3.3f", fCorrectTrue[key][0],fCorrectTrue[key][1]);
// 	  report << buf;
// 	  if(icat<3) report << " & " << flush; 	  
// 	}
//       report << "\\\\" << endl;

//       report << "$f_{correct}^{data}$ &" << flush;
//       for(int icat=0; icat<4; icat++)
// 	{
// 	  std::pair<int,TString> key(jetBins[ijet],cats[icat]);
// 	  char buf[100];
// 	  sprintf(buf,"%3.3f $\\pm$ %3.3f", fCorrectData[key][0],fCorrectData[key][1]);
// 	  report << buf;
// 	  if(icat<3) report << " & " << flush; 	  
// 	}
//       report << "\\\\" << endl;
//     }
//   report << endl;
//   cout << report.str() << endl;


//   //runPileupCheck(url);

// */


	  
// //
// void runCalibration(TString url, JSONWrapper::Object &jsonF, Double_t lumi, int maxPE, Int_t jetbin, bool freezeResults)
// {


//   /*

//   ZZ2l2nuSummaryHandler evHandler, ensembleHandler;
//   // MisassignmentMeasurement *misMeasurement = exclusiveMisMeasurements[jetbin];
//   // misMeasurement->resetHistograms(true);
//   TString peTag("Jets: "); if(jetbin==JETINCLUSIVE) peTag += "inclusive"; else peTag += jetbin; peTag += " Lumi:"; peTag += int(lumi);
//   for(int iPE=1; iPE<=maxPE; iPE++)
//     {
//       if(iPE%5==0) { printf("\r[ %s ]  completed %d/%d ",peTag.Data(),iPE,maxPE); cout << flush; }


	  

//   //
//   // map the samples per type
//   //
//   //signal
  
//   TString syst("");
//   if(syst=="matchingdown")                    procYields["TTJets_matchingdown"]=0;
//   else if(syst=="matchingup")                 procYields["TTJets_matchingup"]=0;
//   else if(syst=="scaledown")                  procYields["TTJets_scaledown"]=0;
//   else if(syst=="scaleup")                    procYields["TTJets_scaleup"]=0;
//   else if(syst=="flavor")                     procYields["TT2l_R0"]=0;
//   else if(syst!="singletop" && syst!="notop") procYields["TTJets_signal"]=0;
//   //other sources of t->Wb
//   if(syst!="notop")
//     {
//       procYields["SingleTbar_tW"]=0;
//       procYields["SingleT_tW"]=0;
//       procYields["SingleTbar_t"]=0;
//       procYields["SingleT_t"]=0;
//       procYields["SingleTbar_s"]=0;
//       procYields["SingleT_s"]=0;
//       procYields["TTJets"]=0;
//     }
//   //other backgrounds
//   procYields["WW"]=0;
//   procYields["ZZ"]=0;
//   procYields["WZ"]=0;
//   procYields["DYJetsToLL"]=0;

//   //
//   // map the trees with the events for the pseudo-experiments
//   //
//   //  cout << "[Filling trees]" << endl;
//   TFile *inF = TFile::Open(url);
//   for(std::map<TString,Int_t>::iterator it = procYields.begin(); it!= procYields.end(); it++)
//     {
//       TString tname=it->first + "/data";
//       TTree *t = (TTree *) inF->Get(tname);
//       if(t==0) continue;
//       Float_t xsecWeight;
//       t->GetBranch("xsecWeight")->SetAddress(&xsecWeight);
//       t->GetEntry(0);
//       procTrees[it->first]  = t;
//       double sf=1.0;
//       if(it->first.Contains("DY") ) sf=(0.5*1.0+0.25*1.73+0.25*1.68);
//       if(it->first.Contains("DY") && syst=="dyup")   sf*=1.15;
//       if(it->first.Contains("DY") && syst=="dydown") sf*=0.85;
//       procYields[it->first] = sf*lumi*xsecWeight*t->GetEntriesFast();
//       //      cout << "\t" << it->first <<  " " <<  procYields[it->first] << " evts" << endl; 
//     }
//       bool initEnsemble(true);
//       for(std::map<TString,TTree*>::iterator it = procTrees.begin(); it!= procTrees.end(); it++)
// 	{

// 	  bool hasSignal(it->first.Contains("TTJets") || it->first.Contains("SingleT") );
// 	  TTree *t=it->second;
// 	  bool result = evHandler.attachToTree(t);
// 	  if(!result) continue;
	  
// 	  //clone tree structure only and release to base directory
// 	  //NB the addresses of the branches are connected to the original tree
// 	  //   when filling the tree the values in the evHandler.evSummary_ structure will be copied
// 	  //   for details cf. http://root.cern.ch/root/html/TTree.html#TTree:CloneTree
// 	  if(initEnsemble)
// 	    {
// 	      ensembleHandler.initTree(t->CloneTree(0),false);
// 	      ensembleHandler.getTree()->SetDirectory(0);
// 	      initEnsemble=false;
// 	    }
	  
// 	  int nevts   = evHandler.getEntries();
// 	  if(nevts==0) continue;

// 	  //generate randomly according to expectations
// 	  float nevExpected = procYields[it->first];
// 	  Int_t nevToGenerate=gRandom->Poisson(nevExpected);
// 	  std::set<int> evtList;
// 	  for(int iev=0; iev<nevToGenerate; iev++)
// 	    {
// 	      int ientry=gRandom->Uniform(0,nevts);
// 	      if(evtList.find(ientry)!=evtList.end()) continue;
// 	      evtList.insert(ientry);
// 	      evHandler.getEntry(ientry);
// 	      evHandler.evSummary_.nmeasurements=1;
// 	      evHandler.evSummary_.evmeasurements[0]=hasSignal;
// 	      ensembleHandler.fillTree();
// 	    }	  
// 	}
      
//       //take control of the filled ensemble and run the misassignment measurement
//       ensembleHandler.attachToTree( ensembleHandler.getTree() );
//       misMeasurement->measureMisassignments( ensembleHandler, 180, 40, false, jetbin, syst );

//       //reset tree used ensemble
//       ensembleHandler.getTree()->Delete("all");
//     }
  

//   //
//   // apply a similar procedure to data
//   //
//   TString data[]={"DoubleMuMay10ReReco",  "DoubleElectronMay10ReReco",  "MuEGMay10ReReco",
// 		  "DoubleMuPromptRecov4", "DoubleElectronPromptRecov4", "MuEGPromptRecov4",
// 		  "DoubleMu05AugReReco",  "DoubleElectron05AugReReco",  "MuEG05AugReReco",
// 		  "DoubleMuPromptRecov6", "DoubleElectronPromptRecov6", "MuEGPromptRecov6"};
//   const size_t ndata=sizeof(data)/sizeof(TString);
//   bool initEnsemble(true);
//   for(size_t idata=0; idata<ndata; idata++)
//     {
      
//       TTree *t=(TTree *) inF->Get(data[idata]+"/data");
//       bool result = evHandler.attachToTree(t);
//       if(!result) continue;
//       if(initEnsemble)
// 	{
// 	  ensembleHandler.initTree(t->CloneTree(0),false);
// 	  ensembleHandler.getTree()->SetDirectory(0);
// 	  initEnsemble=false;
// 	}
//       int nevts   = evHandler.getEntries();
//       if(nevts==0) continue;
//       for(int iev=0; iev<nevts; iev++)
// 	{
// 	  evHandler.getEntry(iev);
// 	  ensembleHandler.fillTree();
// 	}	  
//     }
//   ensembleHandler.attachToTree( ensembleHandler.getTree() );
//   TString datamode(syst=="metcut"?syst:"");
//   misMeasurement->measureMisassignments( ensembleHandler, 180, 40, true, jetbin, datamode );
//   ensembleHandler.getTree()->Delete("all");
//   misMeasurement->finishMonitoringHistograms();

//   //save the bias (MC truth and data are retrieved if results are to be frozen)
//   for(int icat=0; icat<4; icat++)
//     {
//       std::pair<int,TString> key(jetbin,cats[icat]);
//       if(freezeResults)
// 	{
// 	  fCorrectTrue[key] = misMeasurement->getTrueCorrectPairsFraction(cats[icat]);
// 	  fCorrectData[key] = misMeasurement->getCorrectPairsFraction(cats[icat]);
// 	  fCorrectSigOnlyTrue[key] = misMeasurement->getTrueCorrectPairsFraction(cats[icat],true);
// 	}

//       //if the number of pseudo-experiments is large get the bias estimate from the gaussian fit
//       TH1 *biasH = misMeasurement->controlHistos.getHisto("bias",cats[icat]);
//       Float_t fCorrBias=biasH->GetMean();
//       Float_t fCorrBiasErr=biasH->GetMeanError();
//       if(maxPE>20)
// 	{
// 	  TF1 *biasFit=biasH->GetFunction("gaus");
// 	  fCorrBias=biasFit->GetParameter(1);
// 	  fCorrBiasErr=biasFit->GetParError(1);
// 	}
      
//       TGraphErrors *gr = calibCurves[key];
//       Int_t ipt=gr->GetN();
//       gr->SetPoint(ipt,lumi,fCorrBias);
//       gr->SetPointError(ipt,0,fCorrBiasErr);



//     }



//   if(!freezeResults) return;
      
//   //display the results
//   char titBuf[250];
//   sprintf(titBuf,"CMS preliminary, #sqrt{s}=7 TeV, #int L=%3.1f fb^{-1}",float(dataLumi/1000));
//   setStyle();
//   gStyle->SetOptFit(0);
//   TCanvas *mljc = new TCanvas("mljc","mljc",true);
//   mljc->SetCanvasSize(1200,1200);
//   mljc->SetWindowSize(1200,1200);
//   mljc->Divide(2,2);
//   TCanvas *misc = new TCanvas("misc","misc",true);
//   misc->SetCanvasSize(1200,1200);
//   misc->SetWindowSize(1200,1200);
//   misc->Divide(2,2);
//   TCanvas *misresc = new TCanvas("misresc","misc",true);
//   misresc->SetCanvasSize(1200,1200);
//   misresc->SetWindowSize(1200,1200);
//   misresc->Divide(2,2);
//   TCanvas *biasc = new TCanvas("biasc","biasc",true);
//   biasc->SetCanvasSize(1200,1200);
//   biasc->SetWindowSize(1200,1200);
//   biasc->Divide(2,2);
//   TCanvas *pullc = new TCanvas("pullc","pullc",true);
//   pullc->SetCanvasSize(1200,1200);
//   pullc->SetWindowSize(1200,1200);
//   pullc->Divide(2,2);

//   for(int icat=0; icat<4; icat++)
//     {
//       //get the histograms
//       TH1 *mcWrongH = misMeasurement->controlHistos.getHisto("avgwrongmlj",cats[icat]);
//       formatPlot(mcWrongH,809,1,1,0,1001,true,false,1,809,809);
//       mcWrongH->SetTitle("Wrong assignments");
//       mcWrongH->Scale(dataLumi/lumi);

//       TH1 *mcCorrectH= misMeasurement->controlHistos.getHisto("avgcorrectmlj",cats[icat]);
//       formatPlot(mcCorrectH,614,1,1,0,1001,true,false,1,614,614);
//       mcCorrectH->SetTitle("Correct assignments");
//       mcCorrectH->Scale(dataLumi/lumi);

//       TH1 *dataH=misMeasurement->controlHistos.getHisto("datainclusivemlj",cats[icat]);
//       dataH->SetTitle("data");

//       TH1 *mcmodelH=misMeasurement->controlHistos.getHisto("avgwrongmodelmlj",cats[icat]);
//       mcmodelH->SetTitle("MC model");
//       formatPlot(mcmodelH,8,9,2,1,0,true,false,8,8,8);
//       mcmodelH->Scale(dataLumi/lumi);
      
//       TH1 *datamodelH = misMeasurement->controlHistos.getHisto("datawrongmodelmlj",cats[icat]);
//       datamodelH->SetTitle("data model");
//       formatPlot(datamodelH,1,9,2,1,0,true,false,1,1,1);
      
//       TH1 *mcSubtractedH = (TH1 *)mcCorrectH->Clone("mcres"+cats[icat]);
//       mcSubtractedH->Add(mcWrongH);
//       mcSubtractedH->Add(mcmodelH,-1);
//       formatPlot(mcSubtractedH,8,9,2,1,0,true,false,8,8,8);
//       mcSubtractedH->SetTitle("MC subtracted");

//       TH1 *dataSubtractedH = (TH1 *)dataH->Clone("datares"+cats[icat]);
//       dataSubtractedH->Add(datamodelH,-1);
//       dataSubtractedH->SetTitle("data subtracted");

//       TH1 *biasH = misMeasurement->controlHistos.getHisto("bias",cats[icat]);
//       formatPlot(biasH,1,1,2,20,0,true,false,1,1,1);
      
//       TH1 *pullH = misMeasurement->controlHistos.getHisto("pull",cats[icat]);
//       formatPlot(pullH,1,1,2,20,0,true,false,1,1,1);
      
//       TList nullList;
//       TList dataList;   dataList.Add(dataH);
//       TList mcList;     mcList.Add(mcWrongH);       mcList.Add(mcCorrectH);    
//       TList modelList;  modelList.Add(mcmodelH);    modelList.Add(datamodelH);
      
//       Tpad *p=(TPad *) mljc->cd(icat+1);
//       TLegend *leg = showPlotsAndMCtoDataComparison(p,mcList,nullList,dataList);     
//       TPad *sp=(TPad *)p->cd(1);
//       if(icat==0)  formatForCmsPublic(sp,leg,titBuf,3);
//       else         leg->Delete();
//       sp->SetLogx(true);
//       sp->SetLogy(false);
//       TPaveText *pave = new TPaveText(0.20,0.8,0.35,0.9,"NDC");
//       pave->SetBorderSize(0);
//       pave->SetFillStyle(0);
//       pave->SetTextFont(42);
//       pave->AddText(catTits[icat]);
//       pave->Draw();
//       sp=(TPad *)p->cd(2);
//       sp->SetLogx(true);
      
//       p=(TPad *)misc->cd(icat+1);
//       leg = showPlotsAndMCtoDataComparison( p,mcList,modelList,dataList);     
//       sp=(TPad *)p->cd(1);
//       if(icat==0)  formatForCmsPublic(sp,leg,titBuf,4);
//       else         leg->Delete();
//       sp->SetLogx(true);
//       sp->SetLogy(false);
//       pave->DrawClone();
//       sp=(TPad *)p->cd(2);
//       sp->SetLogx(true);

//       p=(TPad *)misresc->cd(icat+1);
//       TList mcCorList;    mcCorList.Add(mcCorrectH); 
//       TList mcresList;    mcresList.Add(mcSubtractedH);
//       TList dataresList;  dataresList.Add(dataSubtractedH);
//       leg = showPlotsAndMCtoDataComparison(p,mcCorList,mcresList,dataresList);     
//       sp=(TPad *)p->cd(1);
//       if(icat==0)  formatForCmsPublic(sp,leg,titBuf,2);
//       else         leg->Delete();
//       sp->SetLogx(true);
//       sp->SetLogy(false);
//       pave->DrawClone();
//       sp=(TPad *)p->cd(2);
//       sp->SetLogx(true);

//       p=(TPad *)biasc->cd(icat+1);
//       biasH->Draw("e1");
//       leg=p->BuildLegend();
//       if(icat==0)  formatForCmsPublic(p,leg,"CMS simulation",2);
//       leg->Delete();
//       pave->DrawClone();

//       p=(TPad *)pullc->cd(icat+1);
//       pullH->Draw("e1");
//       leg=p->BuildLegend();
//       if(icat==0)  formatForCmsPublic(p,leg,"CMS simulation",2);
//       leg->Delete();
//       pave->DrawClone();
//     }

//   TString postfix("");
//   if(jetbin>0) { postfix+="_"; postfix += jetbin; }

//   mljc->Modified();
//   mljc->Update();
//   mljc->SaveAs("Mlj"+postfix+".C");
//   mljc->SaveAs("Mlj"+postfix+".png");
//   mljc->SaveAs("Mlj"+postfix+".pdf");

//   misc->Modified();
//   misc->Update();
//   misc->SaveAs("MljWithModel"+postfix+".C");
//   misc->SaveAs("MljWithModel"+postfix+".png");
//   misc->SaveAs("MljWithModel"+postfix+".pdf");

//   misresc->Modified();
//   misresc->Update();
//   misresc->SaveAs("MljSubtracted"+postfix+".C");
//   misresc->SaveAs("MljSubtracted"+postfix+".png");
//   misresc->SaveAs("MljSubtracted"+postfix+".pdf");

//   biasc->Modified();
//   biasc->Update();
//   biasc->SaveAs("MisassignmentBias"+postfix+".C");
//   biasc->SaveAs("MisassignmentBias"+postfix+".png");
//   biasc->SaveAs("MisassignmentBias"+postfix+".pdf");

//   pullc->Modified();
//   pullc->Update();
//   pullc->SaveAs("MisassignmentPull"+postfix+".C");
//   pullc->SaveAs("MisassignmentPull"+postfix+".png");
//   pullc->SaveAs("MisassignmentPull"+postfix+".pdf");

//   inF->Close();    
// 	      */
// }



// //
// void runPileupCheck(TString url)
// {
//   /*
//   SelectionMonitor controlHistos;
  
//   //generate fake data with different pileup average
//   double massAxis[]={10,20,30,40,50,60,70,80,90,100,115,130,145,160,185,200,250,300,400,500,1000,2000};
//   int nMassBins=sizeof(massAxis)/sizeof(double)-1;
//   for(int i=0; i<12; i++)
//     {
//       TString name("gen"); name+=i;
//       controlHistos.addHistogram( new TH1D(name+"pu",name+" PU",36,-0.5,35.5) );
//       TH1 *h=controlHistos.getHisto(name+"pu");
//       for(int ibin=1; ibin<=h->GetXaxis()->GetNbins(); ibin++)
// 	{
// 	  Float_t binCenter=h->GetXaxis()->GetBinCenter(ibin);
// 	  h->SetBinContent(ibin,TMath::Poisson(binCenter,ibin));
// 	}
//       controlHistos.addHistogram( new TH1D(name+"mlj","Lepton-jet spectrum;Invariant Mass(l,j) [GeV/c^{2}];Lepton-jet pairs",nMassBins,massAxis) );
//     }

//   //now build the mlj spectrum
//   TFile *inF = TFile::Open(url);
//   for(std::map<TString,Int_t>::iterator it = procYields.begin(); it!= procYields.end(); it++)
//     {
//       TString tname=it->first + "/data";
//       TTree *t = (TTree *) inF->Get(tname);
//       if(t==0) continue;

//       top::EventSummaryHandler evSummaryHandler;
//       if( !evSummaryHandler.attachToTree( t) ) continue;
//       const Int_t totalEntries=evSummaryHandler.getEntries();

//       //1. fill the PU in MC
//       controlHistos.addHistogram( new TH1D(it->first + "pumc",it->first , 36,-0.5,35.5) );
//       TH1 *MCH=controlHistos.getHisto(it->first+"pumc");
//       float procWeight=dataLumi*totalEntries;
//       for(int ientry=0; ientry<totalEntries; ientry++)
// 	{
// 	  evSummaryHandler.getEntry(ientry);
// 	  top::EventSummary_t &ev = evSummaryHandler.getEvent();
// 	  procWeight *= ev.xsecWeight;
// 	  MCH->Fill(ev.ngenpu);
// 	}
      
//       //2. instantiate 12 reweighters
//       edm::LumiReWeighting *LumiWeights[12];
//       for(int i=0; i<12; i++)
// 	{
// 	  std::vector< float > MC_distr,Lumi_distr;
// 	  TString name("gen"); name+=i;
// 	  TH1 *LumiH=controlHistos.getHisto(name+"pu");
// 	  for(int ibin=1; ibin<=MCH->GetXaxis()->GetNbins(); ibin++)
// 	    {
// 	      MC_distr.push_back( MCH->GetBinContent(i) );
// 	      Lumi_distr.push_back( LumiH->GetBinContent(i) );
// 	    }
// 	  LumiWeights[i] = new edm::LumiReWeighting(MC_distr,Lumi_distr);
// 	}


//       //3. construct the weighted spectra
//       std::vector<float> weight(12,0);
//       for(int ientry=0; ientry<totalEntries; ientry++)
// 	{
// 	  evSummaryHandler.getEntry(ientry);
// 	  top::EventSummary_t &ev = evSummaryHandler.getEvent();
// 	  top::PhysicsEvent_t phys = getPhysicsEventFrom(ev);
	  
// 	  //get the weights
// 	  for(int i=0; i<12; i++)  weight[i]=LumiWeights[i]->weight( ev.ngenpu )*procWeight;

// 	  //inclusive mlj spectrum
// 	  for(top::PhysicsObjectJetCollection::iterator it=phys.jets.begin(); it!=phys.jets.end(); it++)
// 	    {
//               if(it->pt()<=30 || fabs(it->eta())>=2.5) continue;
             
// 	      for(int i=0; i<12; i++) 
// 		{
// 		  TString name("gen"); name+=i;
// 		  for(int ilep=0; ilep<2; ilep++)
// 		    {
// 		      LorentzVector lj=phys.leptons[ilep]+*it;
// 		      float mlj=lj.mass();
// 		      controlHistos.fillHisto(name+"mlj","all",mlj,weight[i],true);
// 		    }
// 		}
// 	    }
	  
// 	}

//       //delete the lumi reweighters
//       for(int i=0; i<12; i++) delete LumiWeights[i];
//     }
  
//   //now display the results
//   TCanvas *mljc = new TCanvas("mljc","mljc",true);
//   mljc->SetCanvasSize(600,600);
//   mljc->SetWindowSize(600,600);
//   TGraphErrors *mljavg = new TGraphErrors; mljavg->SetMarkerStyle(20);
//   TGraphErrors *mljrms = new TGraphErrors; mljrms->SetMarkerStyle(20);
//   for(int i=0; i<12; i++) 
//     {
//       TString name("gen"); name+=i;
//       TH1 *h=controlHistos.getHisto(name+"mlj","all");
//       h->SetLineColor(kGreen-9+i);
//       h->SetMarkerColor(kGreen-9+i);
//       h->SetFillStyle(0);
//       h->SetMarkerStyle(1);
//       h->Draw(i==0?"hist":"histsame");

//       mljavg->SetPoint(i,i,h->GetMean());
//       mljavg->SetPointError(i,i,h->GetMeanError());

//       mljrms->SetPoint(i,i,h->GetRMS());
//       mljrms->SetPointError(i,i,h->GetRMSError());
//     }
//   TLegend *leg=mljc->BuildLegend();
//   formatForCmsPublic(mljc,leg,"CMS simulation",10);
//   mljc->Modified();
//   mljc->Update();
//   mljc->SaveAs("MljPU.C");
//   mljc->SaveAs("MljPU.png");
//   mljc->SaveAs("MljPU.pdf");

//   TCanvas *mljmomc = new TCanvas("mljmomc","mljmomc",true);
//   mljmomc->SetCanvasSize(600,600);
//   mljmomc->SetWindowSize(600,600);
//   mljmomc->Divide(1,2);
//   mljmomc->cd(1);
//   mljavg->Draw("ap");
//   mljavg->GetXaxis()->SetTitle("Average pileup");
//   mljavg->GetYaxis()->SetTitle("Average M_{lj}");
//   mljmomc->cd(2);
//   mljrms->Draw("ap");
//   mljrms->GetXaxis()->SetTitle("Average pileup");
//   mljrms->GetYaxis()->SetTitle("RMS M_{lj}");
//   mljmomc->Modified();
//   mljmomc->Update();
//   mljmomc->SaveAs("MljPUMomenta.C");
//   mljmomc->SaveAs("MljPUMomenta.png");
//   mljmomc->SaveAs("MljPUMomenta.pdf");

//   */
// }
