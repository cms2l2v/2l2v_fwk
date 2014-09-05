#include "TList.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1D.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TSystem.h"
#include "TProfile.h"
#include "THStack.h"
#include "TLegend.h"
#include "TObjArray.h"
#include "TSpline.h"

#include "UserCode/llvv_fwk/interface/tdrstyle.h"
#include "UserCode/llvv_fwk/interface/dFtM.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include <vector>

using namespace std;
using namespace RooFit;

void printHelp();
void buildSFbHistos(TString wsurl);

/**
   @author P. Musella
*/
class GraphToTF1 {
public:
  GraphToTF1( TString name, TGraph * g ) { sp_ = new TSpline3(name,g); };
  double operator() (double *x, double *p) {
    return sp_->Eval( x[0] ) - p[0];
  };
  TSpline * sp_;
};



//
void printHelp()
{
  printf("--help       --> print this\n");
  printf("--flav       --> flavour configuration file\n");
  printf("--btag       --> tag efficiencies configuration file\n");
  printf("--in         --> input file with b-tag multiplicity distribution\n");
  printf("--ws         --> input file with a workspace\n");
  printf("--fit        --> fit type code (SF_b=0,SF_q=1,SF_b_vs_SF_q=2)\n");
  printf("Example: rundFtM --fit 0 --in plotter.root --flav flavParams_cfg.json --btag csvL_cfg.json\n");
  printf("         rundFtM --fit 0 --ws RooWorkspace.root\n"); 
}

//
int main(int argc, char* argv[])
{
  // load framework libraries                                                                                          
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  setTDRStyle();
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(0);
  
  TString url(""),wsurl("");
  TString flavParFile(""),btagParFile("");
  Int_t fitType(dFtM::FIT_SFb);
  for(int i=1;i<argc;i++)
    {
      string arg(argv[i]);
      if(arg.find("--help")!=string::npos)              { printHelp();	  return 0; }
      if(arg.find("--in")!=string::npos && i+1<argc)    { url=argv[i+1];         gSystem->ExpandPathName(url);         i++;  printf("in      = %s\n", url.Data());              }
      if(arg.find("--ws")!=string::npos && i+1<argc)    { wsurl=argv[i+1];       gSystem->ExpandPathName(wsurl);       i++;  printf("ws      = %s\n", wsurl.Data());            }
      if(arg.find("--flav")!=string::npos && i+1<argc)  { flavParFile=argv[i+1]; gSystem->ExpandPathName(flavParFile); i++;  printf("flavParFile = %s\n", flavParFile.Data());       }
      if(arg.find("--btag")!=string::npos && i+1<argc)  { btagParFile=argv[i+1]; gSystem->ExpandPathName(btagParFile); i++;  printf("btagFile  = %s\n", btagParFile.Data());    }
      if(arg.find("--fit")!=string::npos && i+1<argc)   { sscanf(argv[i+1],"%d",&fitType);                             i++;  printf("fitType  = %d\n", fitType);                }
    } 
  if(url=="" && wsurl=="")                            { printHelp(); return 0; }
  if(url!="" && (flavParFile=="" || btagParFile=="")) { printHelp(); return 0; }
  if(wsurl!="")                                       { buildSFbHistos(wsurl); return 0;}
  
  //output directory set to the same directory where root file is stored
  TString outDir(gSystem->DirName(url));

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
  
  //run fit
  dFtM fitter(fitType, flavParFile, btagParFile, url);
  //fitter.printConfiguration(cout);
  fitter.fit();
  fitter.save(outDir+("/"+fitter.getWP()+"_workspace.root").c_str());

  //show result
  std::map<TString,TH1F *> bmultHistos=fitter.getExtendedBtagMultiplicityHistograms();
  TCanvas *c=new TCanvas("c","c",1215,750);
  TString ch[]={"ee","mumu","emu"};
  std::vector<TPad *> pads;
  for(size_t ich=0; ich<3; ich++)
    {
      c->cd();
      TPad *p=new TPad("p"+ch[ich],"p"+ch[ich],0.0,1.0-float(ich+1)/3.,1.0,1-float(ich)/3.);
      pads.push_back(p);
      p->SetLeftMargin(0.05);
      p->SetRightMargin(0.05);
      p->SetBottomMargin(ich==2? 0.2 :0);
      p->SetTopMargin(ich==0? 0.1 : 0);
      p->SetTicks(0,0);
      p->Draw();
      p->cd();

      TH1F *h=bmultHistos[ ch[ich]+"_exp" ];
      h->Draw("hist");
      h->GetYaxis()->SetTitleOffset(0.4);
      h->GetYaxis()->SetTitleSize(0.06);
      h->GetYaxis()->SetLabelSize(0.06);
      h->GetXaxis()->SetTitleSize(0.07);
      h->GetXaxis()->SetLabelSize(0.07);
      Float_t ymin(1.0), ymax(1.2*TMath::Max(h->GetMaximum(),bmultHistos[ ch[ich] ]->GetMaximum()));
      h->GetYaxis()->SetRangeUser(ymin,ymax);

      //separate categories with a dashed line (first pad will have the jet kinematics categories as well)
      TLine *l=new TLine(0,ymin,0,ymax);
      l->SetLineColor(kGray);
      l->SetLineStyle(7);
      l->SetLineWidth(2);
      for(Int_t xbin=3; xbin<=h->GetXaxis()->GetNbins(); xbin+=3)
	{
	  Float_t xval(h->GetXaxis()->GetBinUpEdge(xbin));
	  if(xbin!=h->GetXaxis()->GetNbins()) l->DrawLine(xval,ymin,xval,ymax);

	  if(ich) continue;

	  xval=h->GetXaxis()->GetBinLowEdge(xbin-1);
	  TPaveText *pt=new TPaveText(xval,ymax*0.9,xval+1,ymax*0.95);
	  pt->SetBorderSize(0);
	  pt->SetFillStyle(0);
	  pt->SetTextAlign(12);
	  pt->SetTextSize(0.05);
	  pt->SetTextFont(42);
	  TString cat("k"); cat += xbin/3;
	  std::pair<Int_t,Int_t> jetKinematics=fitter.getJetKinematicsForCat(cat.Data());
	  char buf[25];
	  sprintf(buf,"(%d,%d)",jetKinematics.first,jetKinematics.second);
	  pt->AddText(buf);
	  pt->Draw();
	}

      //add background expectations
      bmultHistos[ ch[ich]+"_bckgexp" ]->Draw("histsame");

      //superimpose data
      bmultHistos[ ch[ich] ]->Draw("e1same");

      p->RedrawAxis();

      //label with category
      TPaveText *pt=ich==0 ? new TPaveText(0.85,0.7,0.95,0.8,"brNDC"):new TPaveText(0.85,0.8,0.95,0.9,"brNDC");
      pt->SetBorderSize(0);
      pt->SetFillStyle(0);
      pt->SetTextAlign(12);
      pt->SetTextColor(kBlue);
      pt->SetTextFont(42);
      TString label("["); label += h->GetTitle(); label += " events]";
      pt->AddText(label);
      pt->Draw();

      if(ich) continue;

      pt=new TPaveText(0.02,0.95,0.6,0.98,"brNDC");
      pt->SetBorderSize(0);
      pt->SetFillStyle(0);
      pt->SetTextAlign(12);
      pt->SetTextSize(0.07);
      pt->AddText("CMS preliminary, #sqrt{s}=8 TeV");
      pt->Draw();
  
      TLegend *leg=new TLegend(0.7,0.95,0.99,0.98);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.07);
      leg->AddEntry(bmultHistos[ ch[ich] ],"data","p");
      leg->AddEntry(h,"signal","f");
      leg->AddEntry(bmultHistos[ ch[ich]+"_bckgexp" ],"background","f");
      leg->SetNColumns(3);
      leg->Draw();
    }
  
  
  c->Modified();  c->Update();
  c->SaveAs(outDir+("/bmultextended_"+fitter.getWP()+".png").c_str());
  for(size_t ich=0; ich<3; ich++) pads[ich]->SetLogy();
  c->Modified();  c->Update();
  c->SaveAs(outDir+("/bmultextended_"+fitter.getWP()+"_log.png").c_str());
}  

//
void buildSFbHistos(TString wsurl)
{

  //FIXME: HARDCODED!
  Float_t ptCats[]={30,50,80,120,210,320,500};
  std::map<TString,Float_t> wps;
  wps["csv0"]=0.405;
  wps["csv1"]=0.594;
  wps["csv2"]=0.783;
  wps["csv3"]=0.819;
  wps["csv4"]=0.855;

  //overview efficiency measurements
  std::vector<TGraphErrors *> ptSFgr,ptDiscgr;
  std::vector<TF1 *> ptSFsplines;
  for(size_t i=0; i<sizeof(ptCats)/sizeof(Float_t); i++)
    {
      ptSFgr.push_back( new TGraphErrors );
      TString name("eff"); name+=ptCats[i];
      ptSFgr[i]->SetName(name);
      char buf[20];
      sprintf(buf,"%3.0f<p_{T}<%3.0f",ptCats[i],ptCats[i+1]);
      ptSFgr[i]->SetTitle(buf);
      ptSFgr[i]->SetMarkerStyle(20+i);
      ptSFgr[i]->SetLineWidth(2);
      ptSFgr[i]->SetLineColor(1);

      name="disc"; name+=ptCats[i];
      ptDiscgr.push_back ((TGraphErrors *)ptSFgr[i]->Clone(name));
      ptDiscgr[i]->SetLineWidth(1);
    }

  TObjArray *tkns=wsurl.Tokenize(",");
  for(Int_t iurl=0; iurl<tkns->GetEntriesFast(); iurl++)
    {
      TString url(tkns->At(iurl)->GetName());

      TFile *inF=TFile::Open(url);
      if(inF==0) continue;
      if(inF->IsZombie()) { inF->Close(); continue; }

      TString pf("_"); pf+= iurl;
      TString outDir(gSystem->DirName(url));
      TString tagger(gSystem->BaseName(outDir));

      RooWorkspace *w=(RooWorkspace *) inF->Get("w");
      RooRealVar *mcstat=w->var("mcstat");

      TGraphErrors *mcebGr=new TGraphErrors;
      mcebGr->SetName("mcebGr"+pf); mcebGr->SetMarkerStyle(24); mcebGr->SetLineWidth(2); mcebGr->SetLineColor(1); mcebGr->SetFillStyle(0);
      TGraphErrors *ebGr=(TGraphErrors *)mcebGr->Clone("eb"+pf);
      ebGr->SetMarkerStyle(20);
      TGraphErrors *sfbGr=(TGraphErrors *)mcebGr->Clone("sfb"+pf);
      sfbGr->SetMarkerStyle(20);
      
      const RooArgSet *poi=w->set("poi");
      RooFIter iter = poi->fwdIterator();
      RooAbsArg *fitPar;
      while((fitPar = iter.next())) 
	{
	  TString name(fitPar->GetName());
	  if(!name.Contains("SFb")) continue;

	  RooRealVar *iSFb=(RooRealVar *)w->var(name);
	  
	  name.ReplaceAll("SFb","");
	  name.ReplaceAll("k","k_");
	  RooAbsReal *ieb=w->function("absepsb_"+name);

	  w->loadSnapshot("prefit");
	  Float_t mceb( ieb->getVal() );
	  mcstat->setVal(1.0);
	  Float_t mcebUnc( mceb-ieb->getVal() );
	  mcstat->setVal(0.0);

	  w->loadSnapshot("postfit");
	  Float_t sfb( iSFb->getVal() );
	  Float_t sfbUnc( iSFb->getError() );
	  Float_t eb( mceb*sfb );
	  Float_t ebUnc( mceb*sfbUnc );

	  TString catStr(name); catStr.ReplaceAll("k_","");
	  Int_t icat(catStr.Atoi()-1);
	  Float_t pt=0.5*(ptCats[icat+1]+ptCats[icat]);
	  Float_t ptUnc=0.5*(ptCats[icat+1]-ptCats[icat]);

	  if(icat==0) continue;

	  //fill plots
	  Int_t np=mcebGr->GetN();
	  mcebGr->SetPoint(np,pt,mceb);
	  mcebGr->SetPointError(np,ptUnc,mcebUnc);
	  ebGr->SetPoint(np,pt,eb);
	  ebGr->SetPointError(np,ptUnc,ebUnc);
	  sfbGr->SetPoint(np,pt,sfb);
	  sfbGr->SetPointError(np,ptUnc,sfbUnc);
	  
	  //summary plots
	  np=ptSFgr[icat]->GetN();
	  ptSFgr[icat]->SetPoint(np,wps[tagger],eb);
	  ptSFgr[icat]->SetPointError(np,0,ebUnc);
	}    
  
      //absolute efficiency and scale factor
      TCanvas *effc= new TCanvas("effc"+pf,"effc"+pf,500,500);
      effc->cd();
      TPad* p = new TPad("effp0"+pf,"effp0"+pf,0,0.5,1.0,1.0); p->Draw();   p->cd();
      p->SetTopMargin(0.08); p->SetBottomMargin(0); p->SetLeftMargin(0.12); p->SetRightMargin(0.05);
      mcebGr->Draw("ap");  
      ebGr->Draw("p");
      mcebGr->GetXaxis()->SetTitle("Transverse momentum [GeV]");
      mcebGr->GetXaxis()->SetTitleSize(0.08);  
      mcebGr->GetXaxis()->SetLabelSize(0.07); 
      mcebGr->GetYaxis()->SetTitle("b-tagging efficiency");  
      mcebGr->GetYaxis()->SetTitleSize(0.08);  
      mcebGr->GetYaxis()->SetLabelSize(0.07);  
      mcebGr->GetYaxis()->SetTitleOffset(0.7);       
      Double_t maxVal=mcebGr->GetYaxis()->GetXmax();
      Double_t minVal=mcebGr->GetYaxis()->GetXmin();
      mcebGr->GetYaxis()->SetRangeUser(minVal*0.95,maxVal*1.15);
      
      TPaveText *pt=new TPaveText(0.12,0.95,0.6,0.98,"brNDC");
      pt->SetBorderSize(0);
      pt->SetFillStyle(0);
      pt->SetTextAlign(12);
      pt->SetTextSize(0.07);
      pt->AddText("CMS preliminary, #sqrt{s}=8 TeV");
      pt->Draw();

      TLegend *leg=new TLegend(0.6,0.95,0.95,0.98);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.07);
      leg->AddEntry(mcebGr,"MC","p");
      leg->AddEntry(ebGr,"data","p");
      leg->SetNColumns(2);
      leg->Draw();
      
      effc->cd();
      p = new TPad("effp1"+pf,"effp1"+pf,0,0.,1.0,0.5,0); p->Draw(); p->cd();
      p->SetTopMargin(0.0); p->SetBottomMargin(0.2); p->SetLeftMargin(0.12); p->SetRightMargin(0.05); 
      sfbGr->Draw("ap"); 
      sfbGr->GetXaxis()->SetTitle("Transverse momentum [GeV]");
      sfbGr->GetXaxis()->SetTitleSize(0.08);     
      sfbGr->GetXaxis()->SetLabelSize(0.07); 
      sfbGr->GetYaxis()->SetRangeUser(0.8,1.2);
      sfbGr->GetYaxis()->SetTitle("Data/MC scale");        
      sfbGr->GetYaxis()->SetTitleSize(0.08);    
      sfbGr->GetYaxis()->SetLabelSize(0.07);  
      sfbGr->GetYaxis()->SetTitleOffset(0.7);      
	 
      // gStyle->SetOptFit(1111);
      //TF1 *fitFunc=new TF1("ffunc"+pf,"[0]*(1.0+[1]*x)/(1.0+[2]*x)",sfbGr->GetXaxis()->GetXmin(),sfbGr->GetXaxis()->GetXmax());
      //fitFunc->SetParLimits(0,0.5,2.0);
      //sfbGr->Fit(fitFunc,"RQM+");

      effc->cd();
      effc->Modified();
      effc->Update();
      effc->SaveAs(outDir+"/effc_"+tagger+".png");
      
      //all done with this file
      inF->Close();
    }

  


  //absolute efficiency and scale factor
  TCanvas *overvieweffc= new TCanvas("overvieweffc","overviewc",500,500);
  overvieweffc->cd();
  //  overvieweffc->SetTopMargin(0.05);
  // overvieweffc->SetBottomMargin(0.1);
  // overvieweffc->SetRightMargin(0.05);
  //overvieweffc->SetLeftMargin(0.1);

  TLegend *leg=new TLegend(0.15,0.15,0.9,0.3);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  bool notFilled(true);
  for(size_t i=0; i<ptSFgr.size()-2; i++){
    if(ptSFgr[i]->GetN()==0) { ptSFsplines.push_back(0); continue; }

    ptSFgr[i]->Draw( (notFilled ? "ap" : "p") );
    leg->AddEntry(ptSFgr[i],ptSFgr[i]->GetTitle(),"p");
    ptSFgr[i]->GetXaxis()->SetTitle("CSV discriminator");
    ptSFgr[i]->GetYaxis()->SetTitle("b efficiency");
    ptSFgr[i]->GetYaxis()->SetRangeUser(0,1);
    notFilled=false;

    TString name("effspline"); name+=i;
    GraphToTF1 grToTF1=GraphToTF1( name, ptSFgr[i] );
    TF1 *func = new TF1("f"+name,grToTF1,ptSFgr[i]->GetX()[0],ptSFgr[i]->GetX()[ptSFgr[i]->GetN()-1],1,"GraphToTF1");
    func->SetParameter(0,0.);
    func->SetLineColor(1);
    func->SetLineStyle(7);
    func->Draw("same");
    ptSFsplines.push_back(func);
  }
  leg->SetNColumns(2);
  leg->Draw();
      
  TPaveText *pt=new TPaveText(0.14,0.96,0.6,0.98,"brNDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.04);
  pt->AddText("CMS preliminary, #sqrt{s}=8 TeV");
  pt->Draw();
  
  overvieweffc->Modified();
  overvieweffc->Update();
  overvieweffc->SaveAs("overvieweff.png");

  //absolute efficiency and scale factor
  TCanvas *overviewdiscc= new TCanvas("overviewdiscc","overviewc",500,500);
  overviewdiscc->cd();

  leg=new TLegend(0.15,0.15,0.9,0.3);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  notFilled=true;
  for(size_t i=0; i<ptSFgr.size()-2; i++){

    if(ptSFsplines[i]==0) continue;

    for(Float_t x=wps["csv0"]; x<wps["csv4"]; x+=0.02)
      {
	Double_t disc=ptSFsplines[i]->Derivative(x);
	Int_t np=ptDiscgr[i]->GetN();
	ptDiscgr[i]->SetPoint(np,x,disc);
      }
    ptDiscgr[i]->Draw( (notFilled ? "alp" : "lp") );
    ptDiscgr[i]->GetXaxis()->SetTitle("CSV discriminator");
    ptDiscgr[i]->GetYaxis()->SetTitle("d#varepsilon_{b}/dCSV");
    leg->AddEntry(ptDiscgr[i],ptDiscgr[i]->GetTitle(),"p");
    //ptDiscgr[i]->GetYaxis()->SetRangeUser(0,1);
    notFilled=false;
  }
  leg->SetNColumns(2);
  leg->Draw();
      
  pt=new TPaveText(0.14,0.96,0.6,0.98,"brNDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.04);
  pt->AddText("CMS preliminary, #sqrt{s}=8 TeV");
  pt->Draw();
  
  overviewdiscc->Modified();
  overviewdiscc->Update();
  overviewdiscc->SaveAs("overviewdisc.png");






}
