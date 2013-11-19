#include "TChain.h"
#include "TSystem.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TLegend.h"

#include <iostream>

#include "UserCode/llvv_fwk/interface/tdrstyle.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooKeysPdf.h"
#include "RooNDKeysPdf.h"
#include "RooProdPdf.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooWorkspace.h"
#include "RooStats/SPlot.h"

using namespace RooFit ;
using namespace std;

class EWKZp2JSummary{
public:
  EWKZp2JSummary() { };
  ~EWKZp2JSummary() { };
  Float_t ch,weight,cnorm,mjj,pt1,pt2,eta1,eta2,qg1,qg2,spt,ystar,hardpt,ncjv15,htcjv15,pt3,ystar3,mva;
  void setAsProxyTo(TTree *c)
  {
    c->SetBranchAddress("ch",      &ch);
    c->SetBranchAddress("weight",  &weight);
    c->SetBranchAddress("cnorm",   &cnorm);
    c->SetBranchAddress("mjj",     &mjj);
    c->SetBranchAddress("pt1",     &pt1);
    c->SetBranchAddress("pt2",     &pt2);
    c->SetBranchAddress("eta1",    &eta1);
    c->SetBranchAddress("eta2",    &eta2);
    c->SetBranchAddress("qg1",     &qg1);
    c->SetBranchAddress("qg2",     &qg2);
    c->SetBranchAddress("spt",     &spt);
    c->SetBranchAddress("ystar",   &ystar);
    c->SetBranchAddress("hardpt",  &hardpt);
    c->SetBranchAddress("ncjv15",  &ncjv15);
    c->SetBranchAddress("htcjv15", &htcjv15);
    c->SetBranchAddress("pt3",     &pt3);
    c->SetBranchAddress("ystar3",  &ystar3);
    c->SetBranchAddress("mva",     &mva);
  }
};

Float_t lumi=19704;
TString llDir="~/work/ewkzp2j_5311/ll/";
TString gDir="~/work/ewkzp2j_5311/g/data/qt_tight";
std::map<TString,std::pair<float, std::vector<TString> > > procs;
TString wsUrl="";
EWKZp2JSummary ev;

void printHelp();
void getAvailableSamples();


//
void printHelp()
{
  printf("Options\n");
  printf("--ll       --> directory with the dilepton summary trees\n");
  printf("--g        --> directory with the photon summary trees\n");
  printf("--w        --> an already existing workspace\n");
}


//TString generateWorkspace(Float_t minMjj=750);
//void showTheUnfoldedResult(RooDataSet *dataw_sig, RooDataSet *mcSig, RooRealVar *var, int nbins,TString m_cuts,bool showLegend);
//void drawCMSheader();
//void drawCaption(TString tag);

//
/*
void drawCMSheader()
{
  TPaveText *pt =new TPaveText(0.1,0.95,0.9,0.99,"brNDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->SetTextAlign(12);
  char buf[500];
  sprintf(buf,"CMS preliminary, #sqrt{s}=8 TeV, #scale[0.5]{#int}L=%2.1f fb^{-1}",lumi/1000.);
  pt->AddText(buf);
  pt->Draw();
}
*/
//
 /*
void drawCaption(TString tag)
{
  TPaveText *pt =new TPaveText(0.6,0.9,0.9,0.93,"brNDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->SetTextAlign(12);
  pt->SetTextFont(42);
  pt->AddText(tag);
  pt->Draw();
}
 */


//
void getAvailableSamples()
{ 
  TString dirsToScan[]={llDir,gDir};
  const size_t nDirsToScan=sizeof(dirsToScan)/sizeof(TString);
  for(size_t i=0; i<nDirsToScan; i++)
    {
      TString dir=dirsToScan[i];
      
      //open the directory
      gSystem->ExpandPathName(dir);
      void *dirp = gSystem->OpenDirectory(dir);
      const char *dirEntry;
      
      while ((dirEntry = gSystem->GetDirEntry(dirp)))
	{ 
	  TString fName(dir +"/" + TString(dirEntry));
	  if(!fName.Contains(".root") || !fName.Contains("_summary")) continue;
	  bool isData(fName.Contains("Data8TeV"));
	  if(!isData && (fName.Contains("filt22") || fName.Contains("filt111")) ) continue;
	  if(!isData && fName.Contains("JetsToLL_50toInf")) continue;

	  //filter the process name
	  TString tag(dirEntry);
	  tag.ReplaceAll(".root","");
	  Ssiz_t i=tag.Index("summary")-2;
	  if(tag.Contains("filt")) i = tag.Index("filt")-2;
	  TString substr(tag(i));
	  if(!substr.IsAlpha() && TString(tag(i-1))!="v") tag.Remove(i,2);
	  
	  //check the initial number of events
	  TFile *inF=TFile::Open(fName);
	  TTree *t=(TTree *)inF->Get("ewkzp2j");
	  if(t->GetEntriesFast()==0) { inF->Close(); continue; }
	  ev.setAsProxyTo(t);
	  t->GetEntry(1);
	  float cnorm=ev.cnorm;
	  inF->Close();
	  
	  //add to the map
	  if( procs.find(tag) == procs.end() )
	    {
	      std::pair<float, std::vector<TString> > obj(cnorm, std::vector<TString>(1,fName));
	      procs[tag]=obj;
	    }
	  else
	    {
	      procs[tag].first = (isData ? 1.0 : procs[tag].first+cnorm);
	      procs[tag].second.push_back(fName);
	    }
	}
    }


  //report what has been found
  cout << "[getAvailableSamples] found the following valid samples for analysis" << endl;
  for(std::map<TString,std::pair<float, std::vector<TString> > >::iterator it = procs.begin();
      it!= procs.end(); 
      it++)
    cout << "\t" << it->first << " contains " << it->second.second.size() << " files and self normalization is " << 1./it->second.first << endl;
  cout << "\t...will use these ntuples to build the workspace" << endl;
}

//
void generateWorkspace(Float_t minMjj)
{
  //
  cout << "[generateWorkspace] Instantiating RooWorkspace::w" << endl;
  RooRealVar mjj("mjj","Dijet invariant mass",0,5000);
  RooRealVar pt1("pt1","Leading jet transverse momentum [GeV]",0,1000);
  RooRealVar pt2("pt2","Second jet transvesrse momentum [GeV]",0,1000);
  RooRealVar mva("mva","MVA discriminant",-1.5,1.5);
  //  RooRealVar hardpt("hardpt","Hard process transverse balance [GeV]",0,200);
  RooRealVar pt3("pt3","Third jet transverse momentum [GeV]",0,200);
  RooRealVar ystar3("ystar3","y*_{3}=y_{j3}-(y_{j1}+y_{j2})/2",-6,6);
  RooRealVar ncjv15("ncjv15","Jet multiplicity",0,10);
  RooRealVar htcjv15("htcjv15","H_{T}=#Sigma p_{T} [GeV]",0,250);
  RooRealVar w("w","w",0,1000);
  
  RooDataSet data("data","data",RooArgSet(mjj,mva,pt1,pt2,ystar3,htcjv15,pt3,ncjv15,w),"w");
  RooDataSet sig ("sig", "sig", RooArgSet(mjj,mva,pt1,pt2,ystar3,htcjv15,pt3,ncjv15,w),"w");
  RooDataSet bkg ("bkg", "bkg", RooArgSet(mjj,mva,pt1,pt2,ystar3,htcjv15,pt3,ncjv15,w),"w");
  
  //read the ntuples to dataset
  for(std::map<TString,std::pair<float, std::vector<TString> > >::iterator it=procs.begin(); it!=procs.end(); it++)
    {
      //check nature of the ntuple
      bool isMC(it->first.Contains("MC8TeV"));
      bool isSignal(false);
      bool isDataDriven(it->first.Contains("filt22"));
      bool isData(it->first.Contains("Data8TeV") && !isDataDriven);

      //assign dataset according to nature
      float cnorm=it->second.first;
      RooDataSet *d=&data;
      if(it->first.Contains("DYJJ")){	
	d=&sig;
	isSignal=true;
      }
      else if (isMC || isDataDriven) d=&bkg;

      //read chain of ntuples
      TChain procCh("ewkzp2j");
      for(std::vector<TString>::iterator fIt=it->second.second.begin(); fIt!=it->second.second.end(); fIt++) procCh.AddFile(*fIt);
      ev.setAsProxyTo(&procCh);
      for(int i=0; i<procCh.GetEntries(); i++)
	{
	  procCh.GetEntry(i);
	  if( fabs(ev.ch)!=22 && fabs(ev.ch)!=11*11 && fabs(ev.ch)!=13*13 ) continue;
	  if( ev.mjj<minMjj )        continue;

	  mjj.setVal(ev.mjj);
	  pt1.setVal(ev.pt1);
	  pt2.setVal(ev.pt2);
	  mva.setVal(ev.mva);
	  //hardpt.setVal(ev.hardpt);
	  ystar3.setVal(ev.ystar3);
	  htcjv15.setVal(ev.htcjv15);
	  pt3.setVal(ev.pt3);
	  ncjv15.setVal(ev.ncjv15);
	  
	  float wgt = 1.0;
	  if(isMC)         wgt=ev.weight*lumi/cnorm;	     
	  if(isDataDriven) wgt=ev.weight;
	  w.setVal(wgt);

	  d->add(RooArgSet(mjj,mva,pt1,pt2,/*hardpt,*/ystar3,htcjv15,pt3,ncjv15,w),w.getVal());
	}
      cout << "\t" << procCh.GetEntries() << " events from " << it->first 
	   << " added as data?" << isData << " dataDriven?" << isDataDriven << " MC?" << isMC << " signal?" << isSignal << endl;  
    }

  //import all to workspace
  RooWorkspace *ws=new RooWorkspace("w");  
  ws->import(sig);
  ws->import(bkg);
  ws->import(data);

    cout << "\t Generating signal PDF" << endl;
  RooKeysPdf sigPdf("sigpdf","sigpdf",mva,sig,RooKeysPdf::MirrorBoth,2) ;
  cout << "\t Generating background PDF" << endl;
  RooKeysPdf bkgPdf("bkgpdf","bkgpdf",mva,bkg,RooKeysPdf::MirrorBoth,2) ;
  
  cout << "Fitting model" << endl;
  RooRealVar sigYield("sigYield","Signal",     data.sumEntries()*0.1, 0,  data.sumEntries()*2);
  RooRealVar bkgYield("bkgYield","Background", data.sumEntries()*0.9, 0, data.sumEntries()*2);
  RooAddPdf addmodel("addmodel","Model", RooArgList(sigPdf, bkgPdf), RooArgList(sigYield,bkgYield));
  addmodel.fitTo(data,Extended());
  ws->import(addmodel);  

  //FIXME do this for the residual backgrounds
  // RooRealVar mu("musig", "average", sig.sumEntries());
  // RooRealVar sigma("sigma", "sigma",0.1*sig.sumEntries());
  //RooGaussian gauss("gauss","gaussian PDF", sigYield, mu, sigma);
  //RooProdPdf model("model","Model",RooArgSet(addmodel,gauss));
  //model.fitTo(data,Extended());
  //ws->import(model);  
  
  //save in a workspace
  wsUrl="EWKZp2J_workspace.root";
  ws->writeToFile(wsUrl);
  cout << "Result saved @ " << wsUrl << endl;
}


//
int main(int argc, char* argv[])
{
  setTDRStyle();

  //get input arguments
  for(int i=1;i<argc;i++){
    string arg(argv[i]);
    if(arg.find("--help")        !=string::npos)              { printHelp(); return -1;} 
    else if(arg.find("--ll")     !=string::npos && i+1<argc)  { llDir  = argv[i+1];        i++;  printf("llDir = %s\n", llDir.Data()); }
    else if(arg.find("--g")      !=string::npos && i+1<argc)  { gDir  = argv[i+1];         i++;  printf("gDir = %s\n",  gDir.Data()); }
    else if(arg.find("--w")      !=string::npos && i+1<argc)  { wsUrl = argv[i+1];         i++;  printf("ws = %s\n",    wsUrl.Data()); }
  }

  if(wsUrl==""){
    getAvailableSamples();
    generateWorkspace(750);
  }

  return 0;
}





//
/*
void showTheUnfoldedResult(RooDataSet *dataw_sig, RooDataSet *mcSig, RooRealVar *var, int nbins,TString m_cuts,bool showLegend)
{
  if(var==0) return;
  var->setBins(nbins);
  RooPlot *frame = var->frame();

  //show unfolded data
  RooDataSet *reducedData=(RooDataSet *)dataw_sig->reduce(RooArgSet(*var),m_cuts);
  reducedData->plotOn(frame,DataError(RooAbsData::SumW2),Name("unfolded"));

  //compare to model normalized to the unfolded data
  RooDataSet *reducedMC = (RooDataSet*) mcSig->reduce(RooArgSet(*var),m_cuts) ;
  RooDataHist mcSigBin(TString("binsig")+var->GetName(),"binsig",RooArgSet(*var),*reducedMC);
  RooHistPdf mcSigPdf(TString("pdfsig")+var->GetName(),"pdfsig",*var,mcSigBin) ;
  mcSigPdf.plotOn(frame,LineColor(kBlue),DrawOption("l"),Name("sim"),FillStyle(0));
  //,
  //		  Normalization(reducedData->sumEntries(),RooAbsReal::NumEvent));
  //		  Normalization(reducedData->sumEntries()/reducedMC->sumEntries(),RooAbsReal::RelativeExpeted));
  frame->Draw() ;
  frame->GetYaxis()->SetTitle("Events");
  frame->GetYaxis()->SetTitleOffset(1.2);

  if(showLegend){
    TLegend *leg=new TLegend(0.6,0.75,0.9,0.9,"","brNDC");
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextAlign(12);
    leg->AddEntry("unfolded","Unfolded","p");
    leg->AddEntry("sim","Expected","f");
    leg->Draw();
  }
}


*/

//
/*
void runEWKZp2Junfolding(TString wsUrl, TString baseCut)
{
  RooMsgService::instance().setSilentMode(true);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);

  //if(wsUrl=="") wsUrl=generateWorkspace(dir,minMjj);
  gSystem->ExpandPathName(wsUrl);
  cout << "Retrieving workspace from " << wsUrl << endl;
  TFile *wsF=TFile::Open(wsUrl);
  RooWorkspace *ws=(RooWorkspace *)wsF->Get("w");
  wsF->Close();

  //get objects from workspace
  RooDataSet *data     = (RooDataSet *)ws->data("data");
  RooDataSet *mcSig        = (RooDataSet *)ws->data("sig");
  RooDataSet *mcBkg        = (RooDataSet *)ws->data("bkg");
  RooAbsPdf *constraint    = ws->pdf("mvagauss");
  RooAbsPdf *model         = ws->pdf("mvamodel");
  RooRealVar *sigYield     = ws->var("mvasigYield");
  RooRealVar *bkgYield     = ws->var("bkgYield");
  model->fitTo(*data,Extended(),Constrain(*constraint),Cut(baseCut));
       
 
  TCanvas *cfit=new TCanvas(disc+"fit",disc+"fit",600,600);
  cfit->cd();
  RooPlot *frame=ws->var(disc)->frame(Bins(50),Range(-1.5,1.5));
  data->plotOn(frame,DrawOption("p"),MarkerStyle(20),Cut(baseCut),Name("data"));
  model->plotOn(frame, Components( *(ws->pdf(disc+"bkgpdf")) ), FillStyle(0), LineColor(1), LineWidth(2), DrawOption("LF"), FillStyle(1001), FillColor(kGray), MoveToBack(),Cut(baseCut),Name("bkg"));
  model->plotOn(frame,FillStyle(0),MoveToBack(),Cut(baseCut),Name("total"));
  
  frame->Draw();
  frame->GetYaxis()->SetTitle("Events");
  frame->GetYaxis()->SetTitleOffset(1.2);
  cfit->Modified();
  cfit->Update();
  drawCMSheader();

  TLegend *leg=new TLegend(0.7,0.8,0.9,0.95,"","brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->AddEntry("data","data","p");
  leg->AddEntry("bkg","#Sigma bkg","f");
  leg->AddEntry("total","EWK Z+2j","f");
  leg->Draw();

  frame->GetYaxis()->SetRangeUser(0.9,1e3);
  cfit->SetLogy();

  return;

  // Now we use the SPlot class to add SWeights to our data set
  // based on our model and our yield variables
  RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot from " + disc, 
					       *data , model, 
					       RooArgList(*sigYield,*bkgYield ) );
  ws->import(*data, Rename("dataWithSWeights"));
  data = (RooDataSet*) ws->data("dataWithSWeights");
  data->Print("v");
  RooDataSet *dataw_sig = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,disc+"sigYield_sw") ;
  RooDataSet *bkgw_sig = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,disc+"bkgYield_sw") ;

  //show all
  TString j3cut("ncjv15>0"); if(baseCut!="") j3cut += "&&"+baseCut;
  //TString j3cut(baseCut);

  TCanvas *c=new TCanvas("ystar3c","ystar3c",1200,600);
  c->Divide(2,1);
  c->cd(1);  showTheUnfoldedResult(dataw_sig,mcSig,ws->var("ystar3"),8,j3cut,true); drawCMSheader(); drawCaption("Signal-like");
  c->cd(2);  showTheUnfoldedResult(bkgw_sig,mcBkg,ws->var("ystar3"),8,j3cut,false); drawCaption("Background-like");

  c=new TCanvas("ncjvc15","ncjvc15",1200,600);
  c->Divide(2,1);
  c->cd(1);  showTheUnfoldedResult(dataw_sig,mcSig,ws->var("ncjv15"),5,baseCut,true); drawCMSheader();drawCaption("Signal-like");
  c->cd(2);  showTheUnfoldedResult(bkgw_sig,mcBkg,ws->var("ncjv15"),5,baseCut,false);drawCaption("Background-like");

  c=new TCanvas("htcjvc15","htcjvc15",1200,600);
  c->Divide(2,1);
  c->cd(1);  showTheUnfoldedResult(dataw_sig,mcSig,ws->var("htcjv15"),10,j3cut,true); drawCMSheader();drawCaption("Signal-like");
  c->cd(2);  showTheUnfoldedResult(bkgw_sig,mcBkg,ws->var("htcjv15"),10,j3cut,false);drawCaption("Background-like");

  c=new TCanvas("jetgapc","jetgapc",1200,600);
  c->Divide(2,1);
  c->cd(1);  showTheUnfoldedResult(dataw_sig,mcSig,ws->var("maxcjpt"),10,baseCut,true); drawCMSheader();drawCaption("Signal-like");
  c->cd(2);  showTheUnfoldedResult(bkgw_sig,mcBkg,ws->var("maxcjpt"),10,baseCut,false);drawCaption("Background-like");

}

*/
