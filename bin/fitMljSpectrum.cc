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

#include "THStack.h"
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

#include "UserCode/llvv_fwk/interface/MacroUtils.h"
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

struct FitResult_t
{
  Double_t ncorrect,ncorrect_err,nwrong,nwrong_err,rho,fcorrect,fcorrect_err;
};

class MljAnalyzer
{
public:
  MljAnalyzer() : data(0), ddBkg(0), mcBkg(0), ddMCbkg(0), signal(0), pseudoData(0) { }
  ~MljAnalyzer() { };

  void prepare();
  void analyze(TString tag);

  TH1F *data, *ddBkg, *mcBkg, *ddMCbkg, *signal, *pseudoData;
  std::map<TString, TH1F *> signalSysts;
  
  Double_t ncorrectMC,ncorrectMC_err, nwrongMC,nwrongMC_err, fcorrectMC,fcorrectMC_err;
  FitResult_t dataFit;
  std::map<TString, Double_t> systUnc;

private :
  
  FitResult_t runFit( TH1F *data, TH1F *signal, Double_t signalNorm, TH1F *ddBkg, Double_t ddBkgNorm, Double_t maxComb,bool weightedEvents=kFALSE);

};

std::map<TString, MljAnalyzer> allShapes;
std::vector<TString> allChannels;


//
void printHelp();
void fixTemplate(TH1 *h);
void fillMljData(TString url,JSONWrapper::Object &Root,bool isSyst);


//                                                                                                                                                                                          
void printHelp()
{
  printf("--help     --> print this\n");
  printf("--in       --> input file with data and nominal templates\n");
  printf("--json     --> json file with the description of the samples\n");
  printf("--syst     --> input file with syst samples\n");
  printf("--systJson --> json file with the description of the samples used for systematic uncertainties\n");
  printf("--iLumi    --> total integrated luminosity\n");
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
      if(pdir==0) { cout << "Failed to find " << proc << endl; continue; }

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
	    MljAnalyzer newShape; 
	    newShape.data    = (TH1F *) nomH->Clone(*it+"_data");       newShape.data->Reset("ICE");    newShape.data->SetDirectory(0);
	    newShape.signal  = (TH1F *) nomH->Clone(*it+"_signal");     newShape.signal->Reset("ICE");  newShape.signal->SetDirectory(0);
	    newShape.ddBkg   = (TH1F *) nomH->Clone(*it+"_ddbkg");      newShape.ddBkg->Reset("ICE");   newShape.ddBkg->SetDirectory(0);
	    newShape.mcBkg   = (TH1F *) nomH->Clone(*it+"_mcbkg");      newShape.mcBkg->Reset("ICE");   newShape.mcBkg->SetDirectory(0);
	    newShape.ddMCbkg = (TH1F *) nomH->Clone(*it+"_ddmcbkg");    newShape.ddMCbkg->Reset("ICE"); newShape.ddMCbkg->SetDirectory(0);
	    allShapes[*it]=newShape; 
	  }

	  //update shape
	  MljAnalyzer &shape=allShapes.find(*it)->second;
	  if(isData)      { shape.data->Add(nomH); shape.ddBkg->Add(ddBkgH); continue; }
	  else if(isSyst && nomH) { 
	    TH1F *h=(TH1F *)sigH->Clone(*it+"_"+proc);
	    h->SetDirectory(0);
	    shape.signalSysts[proc]=h;
	    continue; 
	  }
	  else {
	    shape.signal->Add(sigH);
	    shape.mcBkg->Add( mcBkgH );
	    shape.ddMCbkg->Add( ddBkgH );
	    TH1F *systH = (TH1F*) pdir->Get("optim_systs");
	    if(systH==0) continue;
            for(int ivar=1; ivar<=systH->GetNbinsX(); ivar++)
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

  //close file
  inF->Close();
}

//
int main(int argc, char* argv[])
{
  // load framework libraries                                                                                          
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();
  
  //keep RooFit quiet                                                                                                                                                          
  RooMsgService::instance().setSilentMode(true);
  RooMsgService::instance().getStream(0).removeTopic(Minimization);
  RooMsgService::instance().getStream(1).removeTopic(Minimization);
  RooMsgService::instance().getStream(1).removeTopic(Fitting);
  RooMsgService::instance().getStream(0).removeTopic(Eval);
  RooMsgService::instance().getStream(1).removeTopic(Eval);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
  RooMsgService::instance().getStream(1).removeTopic(ObjectHandling);
  RooMsgService::instance().getStream(1).removeTopic(DataHandling);

  //configure
  TString url(""), json(""), systUrl(""), systJson("");
  for(int i=1;i<argc;i++)
    {
      string arg(argv[i]);
      if(arg.find("--help")!=string::npos)                  { printHelp();    return 0; }
      if(arg.find("--in")!=string::npos && i+1<argc)        { url=argv[i+1];      gSystem->ExpandPathName(url);      i++; printf("url       = %s\n", url.Data()); }
      if(arg.find("--json")!=string::npos && i+1<argc)      { json=argv[i+1];     gSystem->ExpandPathName(json);     i++; printf("json      = %s\n", json.Data()); }
      if(arg.find("--systJson")!=string::npos && i+1<argc)  { systJson=argv[i+1]; gSystem->ExpandPathName(systJson); i++; printf("systJson  = %s\n", systJson.Data()); }     
      else if(arg.find("--syst")!=string::npos && i+1<argc)      { systUrl=argv[i+1];  gSystem->ExpandPathName(systUrl);  i++; printf("systUrl   = %s\n", systUrl.Data()); }
      if(arg.find("--iLumi")!=string::npos && i+1<argc)     { sscanf(argv[i+1],"%f",&dataLumi);                      i++; printf("lumi      = %f\n", dataLumi); }
    }
  if(url=="" || json=="") { printHelp(); return 0;}
  
  //general plotting style
  setTDRStyle();
  gStyle->SetPadTickX(0);
  gStyle->SetPadTickY(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

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
      TString key(allChannels[icat]);
      MljAnalyzer &ishape=allShapes[ key ];
      ishape.prepare();
      ishape.analyze(key);
    }
    
  //show report
  cout << report.str() << endl;
}

//
void MljAnalyzer::prepare()
{
  if(signal==0 || mcBkg==0) return;
  TString name(signal->GetName()); name.ReplaceAll("signal","pseudodata");
  pseudoData=(TH1F *) signal->Clone(name);
  pseudoData->Add(mcBkg);
  pseudoData->SetDirectory(0);

  //shift the nominal template according to the statistics
  TH1F* statsUp=(TH1F *)signal->Clone(signal->GetName()+TString("statsUp"));     statsUp->SetDirectory(0);
  TH1F* statsDown=(TH1F *)signal->Clone(signal->GetName()+TString("statsDown")); statsDown->SetDirectory(0);
  for(int ibin=1; ibin<=statsUp->GetXaxis()->GetNbins(); ibin++){
    statsUp  ->SetBinContent(ibin,std::min(2*signal->GetBinContent(ibin), std::max(0.01*signal->GetBinContent(ibin), statsUp  ->GetBinContent(ibin) + statsUp  ->GetBinError(ibin))));
    statsDown->SetBinContent(ibin,std::min(2*signal->GetBinContent(ibin), std::max(0.01*signal->GetBinContent(ibin), statsDown->GetBinContent(ibin) - statsDown->GetBinError(ibin))));
  }          
  signalSysts["statsup"]=statsUp;
  signalSysts["statsdown"]=statsDown;
}  


//
void MljAnalyzer::analyze(TString tag)
{
  //identify number of jets
  int njets(2);
  if(tag.Contains("eq3jets")) njets=3;
  if(tag.Contains("eq4jets")) njets=4;

  //max. combinations per jet
  int maxCombs(2*njets);

  //expectations from simulation
  ncorrectMC     = signal->IntegralAndError(1,signal->GetXaxis()->GetNbins(),ncorrectMC_err);
  nwrongMC       = mcBkg ->IntegralAndError(1,mcBkg->GetXaxis()->GetNbins(), nwrongMC_err);
  fcorrectMC     = ncorrectMC/(ncorrectMC+nwrongMC);
  fcorrectMC_err = sqrt(pow(ncorrectMC_err*(2*ncorrectMC+nwrongMC),2)+pow(ncorrectMC*nwrongMC_err,2))/pow(ncorrectMC+nwrongMC,2);
  
  //perform the fit to the data
  dataFit = runFit( data, signal, ncorrectMC/maxCombs, ddBkg, nwrongMC/maxCombs, maxCombs ); 
  Double_t sfcorrect        = dataFit.fcorrect/fcorrectMC;
  Double_t sfcorrect_err    = sqrt(pow(dataFit.fcorrect_err*fcorrectMC,2)+pow(dataFit.fcorrect*fcorrectMC_err,2))/pow(fcorrectMC,2);
  
  //systematics (evaluated from the difference wrt to MC)
  for(std::map<TString, TH1F *>::iterator it=signalSysts.begin(); it!=signalSysts.end(); it++)
    {
      it->second->Scale(signal->Integral()/it->second->Integral());
      FitResult_t iSystFit = runFit( pseudoData, it->second, ncorrectMC/maxCombs, mcBkg, nwrongMC/maxCombs, maxCombs, true );
      systUnc[it->first]=iSystFit.fcorrect;
    }

  //build a report
  report << endl
	 << "[" << tag << "]" << endl
	 << "f_{correct}^{MC}:\t" << fcorrectMC << " \\pm " << fcorrectMC_err << endl
	 << "f_{correct}:\t" << dataFit.fcorrect << " \\pm " << dataFit.fcorrect_err << endl
	 << "#Kappa^{f_{correct}}_{data/MC}:\t" << sfcorrect << " \\pm " << sfcorrect_err << endl
	 << "Systematic uncertainties" << endl
	 << "------------------------" << endl;

  std::map<TString,std::vector<TString > > systTable;
  systTable["MC stats"].push_back("statsup");             systTable["MC stats"].push_back("statsdown");
  systTable["JES"].push_back("jesup");                    systTable["JES"].push_back("jesdown");
  systTable["JER"].push_back("jerup");                    systTable["JER"].push_back("jerdown");
  systTable["PU"].push_back("puup");                      systTable["PU"].push_back("pudown");
  systTable["uMET"].push_back("umetup");                  systTable["uMET"].push_back("umetdown");
  systTable["p_{T}(top)"].push_back("topptup");           systTable["p_{T}(top)"].push_back("topptdown");
  systTable["ME-PS"].push_back("t#bar{t}systmepsdown");   systTable["ME-PS"].push_back("t#bar{t}systmepsup");
  systTable["Q^{2}"].push_back("t#bar{t}systq2down");     systTable["Q^{2}"].push_back("t#bar{t}systq2up");
  systTable["Had"].push_back("t#bar{t}systpowhegpy");     systTable["Had"].push_back("t#bar{t}systpowheghw");
  systTable["Signal"].push_back("t#bar{t}systpowhegpy");  systTable["Signal"].push_back("t#bar{t}172.5");
  systTable["UE"].push_back("t#bar{t}systtunep11");       systTable["UE"].push_back("t#bar{t}systtunep11mpihi"); systTable["UE"].push_back("t#bar{t}systtunep11tev");
  systTable["CR"].push_back("t#bar{t}systtunep11");       systTable["CR"].push_back("t#bar{t}systtunep11nocr");
  systTable["m_{top}"].push_back("t#bar{t}166.5");        systTable["m_{top}"].push_back("t#bar{t}172.5");
  float totalSyst(0);
  for( std::map<TString,std::vector<TString > >::iterator it = systTable.begin();
       it!= systTable.end(); 
       it++)
    {
      int ncontribs(0);
      float totalUnc(0);

      int firstContrib(0);
      float nomRef(fcorrectMC);
      if(it->first=="Had" || it->first=="UE" || it->first=="CR" || it->first=="m_{top}") { firstContrib=1; nomRef= systUnc[ it->second[0] ]; }

      for(size_t icontrib=firstContrib; icontrib<it->second.size(); icontrib++)
	{
	  if(systUnc.find(it->second[icontrib])==systUnc.end()) continue;
	  totalUnc += fabs(nomRef-systUnc[ it->second[icontrib] ]);
	  ncontribs++;
	}
      if(ncontribs<1) continue;

      float weight(1.0);
      if(it->first=="m_{top}") weight=6.0;
      //totalUnc=(totalUnc/(weight*ncontribs));
      totalUnc=(totalUnc/(weight*ncontribs))*100/nomRef;
      totalSyst += pow(totalUnc,2);
      report << it->first << " " << totalUnc << "%" << endl;
    }
  report << "TOTAL: " << sqrt(totalSyst) << endl;

  //
  // SHOW FIT RESULT
  //
  TCanvas *c=new TCanvas("c","c",600,600);
  c->cd();
  TH1F *dataBinCorr   = (TH1F *) data->Clone("databincorr");     dataBinCorr->SetDirectory(0);
  TH1F *signalBinCorr = (TH1F *) signal->Clone("signalbincorr"); signalBinCorr->SetDirectory(0);
  signalBinCorr->SetFillColor(614);
  signalBinCorr->SetFillStyle(1001);
  signalBinCorr->SetLineColor(614);
  signalBinCorr->Scale(dataFit.ncorrect/signalBinCorr->Integral());
  TH1F *wrongBinCorr  = (TH1F *) ddBkg->Clone("ddbkgbincorr");   wrongBinCorr->SetDirectory(0);
  wrongBinCorr->SetFillColor(822);
  wrongBinCorr->SetFillStyle(1001);
  wrongBinCorr->SetLineColor(822);
  wrongBinCorr->Scale(dataFit.nwrong/wrongBinCorr->Integral());
  for(int ibin=1; ibin<=dataBinCorr->GetXaxis()->GetNbins(); ibin++)
    {
      float binWidth=dataBinCorr->GetXaxis()->GetBinWidth(ibin);
      dataBinCorr->SetBinContent(ibin, dataBinCorr->GetBinContent(ibin)/binWidth);
      dataBinCorr->SetBinError  (ibin, dataBinCorr->GetBinError(ibin)  /binWidth);
      signalBinCorr->SetBinContent(ibin, signalBinCorr->GetBinContent(ibin)/binWidth);
      signalBinCorr->SetBinError  (ibin, signalBinCorr->GetBinError(ibin)  /binWidth);
      wrongBinCorr->SetBinContent(ibin, wrongBinCorr->GetBinContent(ibin)/binWidth);
      wrongBinCorr->SetBinError  (ibin, wrongBinCorr->GetBinError(ibin)  /binWidth);
    }
  
  THStack *stack=new THStack;
  stack->Add(wrongBinCorr,"hist");
  stack->Add(signalBinCorr,"hist");
  stack->Draw();
  dataBinCorr->Draw("e1 same");
  stack->GetXaxis()->SetTitle("Lepton-jet invariant mass [GeV]");
  stack->GetXaxis()->SetNdivisions(208);
  TString ytitle("Events x"); ytitle += (2*njets); ytitle += " / bin  [GeV^{-1}]"; 
  stack->GetYaxis()->SetTitle(ytitle);
  stack->GetYaxis()->SetNdivisions(208);
  stack->GetYaxis()->SetTitleOffset(1.3);

  TLegend *leg=new TLegend(0.4,0.9,0.94,0.94);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->SetNColumns(3);
  leg->AddEntry(dataBinCorr,   "data",          "p");
  leg->AddEntry(wrongBinCorr,  "other pairs",   "f");
  leg->AddEntry(signalBinCorr, "correct pairs", "f");
  leg->Draw("same");
 
  TPaveText *pave = new TPaveText(0.12,0.995,0.84,0.95, "NDC");
  pave->SetFillColor(0); 
  pave->SetFillStyle(0);  
  pave->SetLineColor(0); 
  pave->SetBorderSize(0);
  pave->SetTextAlign(12);
  pave->SetTextSize(0.04);
  char Buffer[1024];
  sprintf(Buffer, "CMS preliminary, #sqrt{s}=%d TeV, #scale[0.5]{#int} L=%.1f fb^{-1}", iEcm, dataLumi/1000);
  pave->AddText(Buffer);
  pave->Draw("same");

  pave = new TPaveText(0.82,0.995,0.95,0.95, "NDC");
  pave->SetFillColor(0);
  pave->SetFillStyle(0);
  pave->SetLineColor(0);
  pave->SetBorderSize(0);
  pave->SetTextFont(42);
  pave->SetTextAlign(12);
  pave->SetTextSize(0.03);
  TString header(tag);
  header.ReplaceAll("mu","#mu");
  header.ReplaceAll("eq",",=");
  header.ReplaceAll("jets"," jets");
  pave->AddText(header);
  pave->Draw("same");

  //display the subtracted data
  c->cd();
  TPad *npad = new TPad("llpad","ll", 0.42, 0.4, 0.95, 0.9);
  npad->Draw();
  npad->cd();
  npad->SetFillStyle(0);
  TH1F *dataSub=(TH1F *) dataBinCorr->Clone(tag+"_dataSub");  dataSub->SetDirectory(0);
  dataSub->Add(wrongBinCorr,-1);
  signalBinCorr->Draw("hist");
  signalBinCorr->GetXaxis()->SetTitle("Lepton-jet invariant mass [GeV]");    signalBinCorr->GetXaxis()->SetTitleOffset(1.0);
  signalBinCorr->GetYaxis()->SetTitle("Lepton-jet pairs / bin [GeV^{-1}]");  signalBinCorr->GetYaxis()->SetTitleOffset(1.3);
  signalBinCorr->GetXaxis()->SetNdivisions(208);
  signalBinCorr->GetYaxis()->SetNdivisions(208);
  dataSub->Draw("e1 same");
  signalBinCorr->GetYaxis()->SetRangeUser(dataSub->GetMinimum()*1.15,dataSub->GetMaximum()*1.15);

  pave = new TPaveText(0.3,0.6,0.9,0.9,"NDC");
  pave->SetBorderSize(0);
  pave->SetTextAlign(12);
  pave->SetFillStyle(0);
  pave->SetTextFont(42);
  pave->SetTextSize(0.055);
  TString buf=utils::toLatexRounded(dataFit.ncorrect,dataFit.ncorrect_err);
  buf.ReplaceAll("\\pm","#pm");
  buf.ReplaceAll("$","");
  pave->AddText("N(t#rightarrow Wq)="+buf);
  buf=utils::toLatexRounded(sfcorrect,sfcorrect_err);
  buf.ReplaceAll("\\pm","#pm");
  buf.ReplaceAll("$","");
  pave->AddText("#Kappa_{data/MC}^{N(t#rightarrow Wq)}="+buf);
  float chi2(0),ndof(0);
  for(int ibin=1; ibin<=dataSub->GetXaxis()->GetNbins(); ibin++)
    {
      ndof++;
      chi2 += 
	TMath::Power(dataSub->GetBinContent(ibin)-signalBinCorr->GetBinContent(ibin),2) /
	(TMath::Power(dataSub->GetBinError(ibin),2)+TMath::Power(signalBinCorr->GetBinContent(ibin),2));
    }
  char chi2Buf[1000];
  sprintf(chi2Buf,"#chi^{2}/ndof : %3.1f /%3.0f", chi2, ndof);
  //sprintf(chi2Buf,"#chi^{2}/ndof : %3.1f /%3.0f  p-val: %3.2f", chi2,ndof, TMath::Prob(chi2,ndof) );
  pave->AddText(chi2Buf);
  pave->Draw("same");

  c->cd();
  c->Modified();
  c->Update();
  c->SaveAs(tag+"_fit_mlj.png");
  c->SaveAs(tag+"_fit_mlj.pdf");
  c->SaveAs(tag+"_fit_mlj.C");
}

//
FitResult_t MljAnalyzer::runFit( TH1F *data, TH1F *signal, Double_t signalNorm, TH1F *ddBkg, Double_t ddBkgNorm, Double_t maxComb, bool weightedEvents)
{
  //create the model
  RooRealVar mlj          ("mlj",         "mlj",         data->GetXaxis()->GetXmin(), data->GetXaxis()->GetXmax());
  RooDataHist dataHist    ("datahist",    "datahist",    mlj, RooFit::Import(*data,kFALSE)); 
  RooDataHist correctHist ("correcthist", "correcthist", mlj, RooFit::Import(*signal,kFALSE));
  RooHistPdf correctPdf   ("correctpdf",  "correctpdf",  RooArgSet(mlj), correctHist);
  RooDataHist modelHist   ("modelhist",   "modelhist",   mlj, RooFit::Import(*ddBkg,kFALSE));
  RooHistPdf modelPdf     ("modelpdf",    "modelpdf",    RooArgSet(mlj), modelHist);
  RooRealVar ntwq         ("ntwq",        "ntwq",        signalNorm, 0, dataHist.sumEntries() );
  RooFormulaVar ncorrect  ("ncorrect",    "@0*@1",       RooArgSet(ntwq,RooConst(maxComb)));
  RooRealVar nother       ("nother",      "nother",      ddBkgNorm,0, dataHist.sumEntries() );
  RooFormulaVar nwrong    ("nwrong",      "@0*@1",       RooArgSet(nother,RooConst(maxComb)));
  RooFormulaVar fcorrect  ("fcorrect",    "@0/(@0+@1)",  RooArgSet(ncorrect,nwrong));
  RooFormulaVar fwrong    ("fwrong",      "1-@0",RooArgSet(fcorrect));
  RooAddPdf     shapeModel("shapemodel",   "shapemodel", RooArgSet(correctPdf,modelPdf), RooArgSet(ncorrect,nwrong));
  
  //fit and save result to return
  FitResult_t r;
  RooFitResult *fitRes = weightedEvents ? 
    shapeModel.fitTo(dataHist,RooFit::Save(kTRUE), Extended(kTRUE), RooFit::Optimize(kTRUE), RooFit::InitialHesse(kTRUE),SumW2Error(kTRUE)):
    shapeModel.fitTo(dataHist,RooFit::Save(kTRUE), Extended(kTRUE), RooFit::Optimize(kTRUE), RooFit::InitialHesse(kTRUE));

  r.ncorrect           = ncorrect.getVal();
  r.ncorrect_err       = (ntwq.getError()/ntwq.getVal())*ncorrect.getVal();
  r.nwrong             = nwrong.getVal();
  r.nwrong_err         = (nother.getError()/nother.getVal())*nwrong.getVal();
  r.rho                = fitRes->correlation(ntwq,nother);
  r.fcorrect           = ntwq.getVal()/(ntwq.getVal()+nother.getVal());
  r.fcorrect_err       = sqrt( 
			      pow(nother.getVal()*ntwq.getError(),2) 
			      + pow( ntwq.getVal()*nother.getError(),2) 
			      + 2* ntwq.getVal() * nother.getVal() * r.rho * ntwq.getError() * nother.getError() 
			      ) 
    / pow( ntwq.getVal()+nother.getVal(),2);

  return r;
}
