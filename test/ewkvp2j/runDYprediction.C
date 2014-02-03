#include "TFile.h"
#include "TH2F.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TMath.h"

#include <fstream>
#include<strstream>
#include <iostream> 
#include <string> 
#include <map>

using namespace std;

TString outDir("Img/");
strstream report;
TStyle *tdrStyle;
TGraph *mjjWgtGr;
std::map<std::pair<TString,TString>,float > hzzSFs;
std::map<TString,float> hzzEWKPhotonSFs;

struct Shape_t
{
  TH1* data;
  std::map<TString, TH1 *> mcbckg, bckg, signal;
  TH1 *ewkBkg;
};

void setTDRStyle();
void drawCMSHeader();
void checkShape(TH1* h,bool isData);
void showShape(const Shape_t &shape,TString SaveName,int anMode);
void drawEvolutionFor(TH1 *data, TH1 *totalPred, TH1 *mcPred, const std::map<TString,TH1 *>&sig, TString SaveName);
TGraphErrors *getEvolutionWithErrors(TH1 *h);
enum AnalysisMode {DISCOVERY,SEARCH};
void runDYprediction(TString llFile="~/work/ewkzp2j_5311/plotter.root",
                     TString gammaFile="~/work/ewkzp2j_5311/plotter_g_tight_qt.root",
                     TString looseGammaFile="~/work/ewkzp2j_5311/plotter_g_loose_qt.root",
                     TString closureFile="~/work/ewkzp2j_5311/closure_dy_closure_g_qt_tight.root",
		     int anMode=DISCOVERY);
TH1 *subtractEWKbkg(TH1 *gData, TH1 *ewkBkg,int anMode=DISCOVERY);


//
void drawCMSHeader()
{
  //cms header
  TPaveText *pave = new TPaveText(0.11,0.94,0.65,0.99, "brNDC");
  pave->SetFillColor(0);
  pave->SetFillStyle(0);  
  pave->SetLineColor(0);
  pave->SetTextAlign(12);
  pave->SetBorderSize(0);
  char buf[500];
  sprintf(buf, "CMS preliminary, #sqrt{s}=8 TeV, #scale[0.5]{#int} L=%.1f fb^{-1}", 19700./1000);
  pave->AddText(buf);
  pave->Draw("same");
}


//
void checkShape(TH1* h,bool isData)
{
  if(h==0) return;

  TString hname(h->GetName());
  bool isHighReg( hname.Contains("highmjj") || hname.Contains("092") || hname.Contains("100"));
  isHighReg |= (hname.EndsWith("1pt") || hname.EndsWith("2pt") || hname.EndsWith("deta") || hname.EndsWith("seta"));// || hname.EndsWith("MLP"));
  isHighReg |= (hname.Contains("qgmva"));
  if(hname.Contains("vbfcjv")) isHighReg=false;

  int rebinFactor(2);
  if( h->InheritsFrom("TH2") ){
    if(!isData)
      {	
	//do not let templates with 0 entries
	for(int xbin=1; xbin<=h->GetXaxis()->GetNbins(); xbin++)
	  for(int ybin=1; ybin<=h->GetYaxis()->GetNbins(); ybin++)
	    {
	      float z=h->GetBinContent(xbin,ybin);
	      if(z==0) h->SetBinContent(xbin,ybin,1e-5);
	    }
      }
  }
  else{
    if(isHighReg) ((TH2 *)h)->Rebin(rebinFactor); 
   
    if(!isData){
      //do not let templates with 0 entries
      for(int xbin=1; xbin<=h->GetXaxis()->GetNbins(); xbin++)
	{
	  float y=h->GetBinContent(xbin);
	  if(y==0) h->SetBinContent(xbin,1e-5);
	}
    }
    
    //add underflow
    double fbin  = h->GetBinContent(0) + h->GetBinContent(1);
    double fbine = TMath::Sqrt(h->GetBinError(0)*h->GetBinError(0)
			+ h->GetBinError(1)*h->GetBinError(1));
    h->SetBinContent(1,fbin);
    h->SetBinError(1,fbine);
    h->SetBinContent(0,0);
    h->SetBinError(0,0);
    
    //add overflow
    int nbins = h->GetNbinsX();
    fbin  = h->GetBinContent(nbins) + h->GetBinContent(nbins+1);
    fbine = sqrt(h->GetBinError(nbins)*h->GetBinError(nbins) 
		 + h->GetBinError(nbins+1)*h->GetBinError(nbins+1));
    h->SetBinContent(nbins,fbin);
    h->SetBinError(nbins,fbine);
    h->SetBinContent(nbins+1,0);
    h->SetBinError(nbins+1,0);
  }
}


//
void addToShape(Shape_t &a, Shape_t &b,int sign=+1)
{
  if(a.data && b.data)           a.data->Add(b.data,sign);
  for(std::map<TString, TH1 *>::iterator it = a.bckg.begin(); it != a.bckg.end(); it++)
    {
      if(b.bckg.find(it->first)==b.bckg.end()) continue;
      it->second->Add( b.bckg[it->first], sign );
    }
  for(std::map<TString, TH1 *>::iterator it = b.bckg.begin(); it != b.bckg.end(); it++)
    {
      if(a.bckg.find(it->first)!=a.bckg.end()) continue;
      TH1 *h=(TH1 *)it->second->Clone();
      h->SetDirectory(0);
      h->Scale(sign);
      a.bckg[it->first]=h;
    }
  for(std::map<TString, TH1 *>::iterator it = a.mcbckg.begin(); it != a.mcbckg.end(); it++)
    {
      if(b.mcbckg.find(it->first)==b.mcbckg.end()) continue;
      it->second->Add( b.mcbckg[it->first], sign );
    }
  for(std::map<TString, TH1 *>::iterator it = b.mcbckg.begin(); it != b.mcbckg.end(); it++)
    {
      if(a.mcbckg.find(it->first)!=a.mcbckg.end()) continue;
      TH1 *h=(TH1 *)it->second->Clone();
      h->SetDirectory(0);
      h->Scale(sign);
      a.mcbckg[it->first]=h;
    }

  for(std::map<TString, TH1 *>::iterator it = a.signal.begin(); it != a.signal.end(); it++)
    {
      if(b.signal.find(it->first)==b.signal.end()) continue;
      it->second->Add( b.signal[it->first], sign );
    }
  for(std::map<TString, TH1 *>::iterator it = b.signal.begin(); it != b.signal.end(); it++)
    {
      if(a.signal.find(it->first)!=b.signal.end()) continue;
      TH1 *h=(TH1 *)it->second->Clone();
      h->SetDirectory(0);
      h->Scale(sign);
      a.signal[it->first]=h;
    }
}

//
Shape_t cloneShape(Shape_t &orig,TString newName)
{
  Shape_t newShape;
  newShape.data      = (TH1 *) orig.data->Clone(newName+"data");           newShape.data->SetDirectory(0);
  for(std::map<TString, TH1 *>::iterator it = orig.bckg.begin(); it != orig.bckg.end(); it++)
    {
      TH1 *h=(TH1 *)it->second->Clone(newName+it->first); h->SetDirectory(0);
      newShape.bckg[it->first]=h;
    }
  for(std::map<TString, TH1 *>::iterator it = orig.mcbckg.begin(); it != orig.mcbckg.end(); it++)
    {
      TH1 *h=(TH1 *)it->second->Clone(newName+it->first); h->SetDirectory(0);
      newShape.mcbckg[it->first]=h;
    }

  for(std::map<TString, TH1 *>::iterator it = orig.signal.begin(); it != orig.signal.end(); it++)
    {
      TH1 *h=(TH1 *)it->second->Clone(newName+it->first); h->SetDirectory(0);
      newShape.signal[it->first]=h;
    }
  return newShape;
}


//
void runDYprediction(TString llFile,TString gammaFile,TString looseGammaFile,TString closureFile,int anMode)
{
  setTDRStyle();

  mjjWgtGr=0;
  mjjWgtGr=new TGraph;
  mjjWgtGr->SetPoint(0,150,0.957);
  mjjWgtGr->SetPoint(1,250,0.992);
  mjjWgtGr->SetPoint(2,350,1.015);
  mjjWgtGr->SetPoint(3,455,1.016);
  mjjWgtGr->SetPoint(4,600,1.011);
  mjjWgtGr->SetPoint(5,850,0.998);
  mjjWgtGr->SetPoint(6,1200,0.940);
  mjjWgtGr->SetPoint(7,1600,0.875);
  mjjWgtGr->SetPoint(8,2000,0.783);


  //open files
  TFile *llIn=TFile::Open(llFile);  
  TFile *gIn=TFile::Open(gammaFile);
  //TString gammaRawFile(gammaFile);
  // gammaRawFile.ReplaceAll("_qt","_raw");
  //  TFile *gRawIn=TFile::Open(gammaRawFile);
  TFile *loosegIn=0;
  if(looseGammaFile!="")
    {
      loosegIn=TFile::Open(looseGammaFile);
      if(loosegIn!=0 && loosegIn->IsZombie()) loosegIn=0;
    }
  TFile *closeF=0;
  if(closureFile!="") 
    {
      closeF=TFile::Open(closureFile);
      if(closeF!=0 && closeF->IsZombie()) closeF=0;
    }

  //
  // histograms to retrieve from the files
  //
  std::vector<std::string> histos;
  std::vector<std::string> cats;  
  std::vector<string> dilprocs;
  std::vector<int> dilColors;
  std::vector<std::string> dilSignal;
  if(anMode==DISCOVERY){
    //std control plots
    histos.push_back("vbfmjj");
    /*
    histos.push_back("vbfcandjet1eta");
    histos.push_back("vbfcandjet2eta");
    histos.push_back("vbfcandjetdeta");      
    histos.push_back("vbfcandjetseta");      
    histos.push_back("vbfcandjetetaprod");
    histos.push_back("vbfdphijj");
    histos.push_back("vbfcandjet1pt");
    histos.push_back("vbfcandjet2pt");
    histos.push_back("vbfqgmva1");
    histos.push_back("vbfqgmva2");
    histos.push_back("vbfspt");
    histos.push_back("vbfystar");
    histos.push_back("vbfhardpt");
    histos.push_back("met");
    histos.push_back("axialmet");
    histos.push_back("vbfystar3");
    histos.push_back("vbfcjv15");
    histos.push_back("vbfhtcjv15");
    histos.push_back("vbfmaxcjvjpt");
    */
    //histos.push_back("Fisher");
    histos.push_back("BDTD");
    histos.push_back("MLP");
    /*
    //soft hadronic activity control plots
    histos.push_back("softjetsvsdetajj");
    histos.push_back("softhtvsdetajj");
    histos.push_back("softinjetsvsdetajj");
    histos.push_back("softinhtvsdetajj");
    histos.push_back("softjetsvsmjj");
    histos.push_back("softhtvsmjj");
    histos.push_back("softinjetsvsmjj");
    histos.push_back("softinhtvsmjj");
    histos.push_back("softjetsvsdphijj");
    histos.push_back("softhtvsdphijj");
    histos.push_back("softinjetsvsdphijj");
    histos.push_back("softinhtvsdphijj");
    histos.push_back("softjetsvsnvtx");
    histos.push_back("softhtvsnvtx");
    histos.push_back("softinjetsvsnvtx");
    histos.push_back("softinhtvsnvtx");
    histos.push_back("vbfystarvsmjj");
    */
    //shape analysis distributions
    //histos.push_back("dijet_deta_shapes");
    //histos.push_back("Fisher_shapes");
    // histos.push_back("BDTD_shapes");
    // histos.push_back("MLP_shapes");
    

    cats.push_back("");
    cats.push_back("mjjq016");
    cats.push_back("mjjq033");
    cats.push_back("mjjq049");
    cats.push_back("mjjq066");
    cats.push_back("mjjq083");
    cats.push_back("mjjq092");
    cats.push_back("mjjq100");
    cats.push_back("mjjgt092");
    cats.push_back("highmjj");
    cats.push_back("highhardpt");
    cats.push_back("lowhardpt");

    //
    //processes to retrieve from the files 
    //
    dilprocs.push_back("VV");                dilColors.push_back(592);
    dilprocs.push_back("Top");               dilColors.push_back(824);
    dilprocs.push_back("W#rightarrow l#nu"); dilColors.push_back(809);
    dilprocs.push_back("Z#rightarrow ll");   dilColors.push_back(831);
    dilprocs.push_back("data");              dilColors.push_back(1);

    dilSignal.push_back("EWK ll2j");
  }
  else{
    histos.push_back("met");
    histos.push_back("balance");
    histos.push_back("mt");
    histos.push_back("mtNM1");
    histos.push_back("mindphijmet");
    histos.push_back("mindphijmetNM1");
    histos.push_back("balanceNM1");
    histos.push_back("axialmet");
    histos.push_back("axialmetNM1");
    histos.push_back("metcount");
    histos.push_back("met_shapes");
    histos.push_back("mt_shapes");

    cats.push_back("");
    cats.push_back("eq0jets");
    cats.push_back("geq1jets");
    cats.push_back("vbf");
    
    dilprocs.push_back("EWK ll#nu#nujj");          dilColors.push_back(804);
    dilprocs.push_back("EWK lljj");                dilColors.push_back(835);
    dilprocs.push_back("WW#rightarrow 2l2#nu");    dilColors.push_back(592);
    dilprocs.push_back("ZZ#rightarrow 2l2#nu");    dilColors.push_back(590);
    dilprocs.push_back("WZ#rightarrow 3l#nu");     dilColors.push_back(596);
    dilprocs.push_back("Top");                     dilColors.push_back(8);
    dilprocs.push_back("W,multijets");             dilColors.push_back(809);
    dilprocs.push_back("Z#rightarrow ll");         dilColors.push_back(831);
    dilprocs.push_back("data");                    dilColors.push_back(1);

    //dilSignal.push_back("ggH(600)#rightarrow ZZ");
    //dilSignal.push_back("qqH(600)#rightarrow ZZ");
  }

  //
  // channels and categories to analyze
  //
  std::vector<string> ch;
  ch.push_back("ee");
  ch.push_back("mumu");
    

  std::string gproc("data (#gamma)");
  float m_gclosurePDFunc(0.02); //2% unc. from PDF variations in closure test

  const size_t nhistos=histos.size();
  const size_t nchs=ch.size();
  const size_t ncats=cats.size();
  const size_t ndilprocs=dilprocs.size();
  const size_t nDilSignals=dilSignal.size();

  //
  // DILEPTON SHAPES AND DATA
  //
  std::vector<string> shapesKeys;
  std::map<string,Shape_t> shapesMap;
  for(size_t ich=0; ich<nchs; ich++)
    {
      for(size_t icat=0; icat<ncats; icat++)
	{
	  for(size_t ih=0; ih<nhistos; ih++)
	    {

	      //initiate a new shape
	      Shape_t m_shape;
	      for(size_t iproc=0; iproc<ndilprocs; iproc++)
		{
		  bool isData(dilprocs[iproc].find("data") != string::npos);		 

		  string hname=dilprocs[iproc]+"/"+ch[ich]+cats[icat]+"_"+histos[ih];
		  TH1 *h=(TH1 *)llIn->Get(hname.c_str());
		  if(h==0)  { cout << "Missing " << hname << endl; continue; }
		  h->SetDirectory(0);
		  h->SetTitle(dilprocs[iproc].c_str());
		  checkShape(h,isData);
		  
		  //save as data or background accordingly
		  if(isData)
		    {
		      m_shape.data=h;
		    }
		  else
		    {
		      h->SetFillColor(dilColors[iproc]);
		      if(dilprocs[iproc].find("Z#rightarrow ll")==string::npos) m_shape.bckg[dilprocs[iproc]]=h;
		      m_shape.mcbckg[dilprocs[iproc]]=h;
		    }

		  //retrieve the signal also
		  if(iproc>0) continue;
		  for(size_t isig=0; isig<nDilSignals; isig++)
		    {
		      string hsigname = dilSignal[isig]+"/"+ch[ich]+cats[icat]+"_"+histos[ih];
		      TH1 *hsig       = (TH1 *)llIn->Get(hsigname.c_str());
		      if(hsig==0) { cout << "Missing " << hsigname << endl; continue; }
		      checkShape(hsig,false);
		      TString title(dilSignal[isig].c_str());
		      if(title=="EWK ll2j") title="EWK Zjj";
		      hsig->SetTitle(title);
		      hsig->SetDirectory(0);
		      hsig->SetLineColor(1+isig/2);
		      hsig->SetLineStyle(1+isig%2);
		      hsig->SetLineWidth(2);
		      m_shape.signal[dilSignal[isig]]=hsig;
		    }
		}

	      //save dilepton shape to map
	      string key=ch[ich]+cats[icat]+"_"+histos[ih];
	      shapesMap[key]=m_shape;
	      shapesKeys.push_back(key);
	    }
	}
    }

  //
  // GAMMA SHAPES FROM PHOTON DATA
  //
  std::string gEwkBkgProc="EWK #gamma+2j";	      
  if(anMode==SEARCH) gEwkBkgProc="V#gamma";
  std::map<string,Shape_t>  gShapesMap;
  for(size_t ich=0; ich<nchs; ich++)
    {
      for(size_t icat=0; icat<ncats; icat++)
	{
	  for(size_t ih=0; ih<nhistos; ih++)
	    {
	      string hname=gproc+"/"+ch[ich]+cats[icat]+"_"+histos[ih];
	      TH1 *h=(TH1 *)gIn->Get(hname.c_str());
	      if(h==0) { cout << " Missing: " << hname <<endl; continue; }
	      bool isTH2( h->InheritsFrom("TH2") );
	      checkShape(h,false);
	      h->SetDirectory(0);
	      h->SetName( TString("g") + h->GetName() );
	      h->SetTitle("QCD Zjj(data)");
	      if(anMode==SEARCH) h->SetTitle("Z#rightarrow ll (data)");
	      // correct uncertainty of the prediction adding in quadrature to the stat
	      //  - | tight - loose | difference
	      //  - | tight - pure | difference in MC (if available)
	      //  - x% for the PDF variations on the closure test

	      TH1 *loose2tightH=0;
	      if(loosegIn){
		loose2tightH=(TH1 *)loosegIn->Get(hname.c_str());
		checkShape(loose2tightH,false);
		loose2tightH->SetDirectory(0);
		loose2tightH->Add(h,-1);
	      }

	      TGraph *pure2tightGr=0;
	      if(closeF) pure2tightGr=(TGraph *)closeF->Get(("ll"+cats[icat]+"_"+histos[ih]+"_ptunc").c_str());
	      
	      if(loose2tightH)
		{
		  if(isTH2){
		    for(int ibin=1; ibin<=loose2tightH->GetXaxis()->GetNbins(); ibin++){
		      for(int jbin=1; jbin<=loose2tightH->GetYaxis()->GetNbins(); jbin++){
			
			float y=loose2tightH->GetYaxis()->GetBinCenter(jbin);
			float z=h->GetBinContent(ibin,jbin);
			
			float statUnc(h->GetBinError(ibin,jbin));
			float loose2tightUnc(fabs(loose2tightH->GetBinContent(ibin,jbin)));
			float pure2tightUnc(pure2tightGr ? (1-pure2tightGr->Eval(y)) : 0.);   pure2tightUnc *= z;
			float closurePDFunc(m_gclosurePDFunc);                                closurePDFunc *= z;
			
			float newerr=sqrt(pow(statUnc,2)+pow(loose2tightUnc,2)+pow(pure2tightUnc,2)+pow(closurePDFunc,2));
			h->SetBinError(ibin,jbin,newerr);
		      }
		    }
		  }
		  else{
		    for(int ibin=1; ibin<=loose2tightH->GetXaxis()->GetNbins(); ibin++){
		      
		      float x=loose2tightH->GetXaxis()->GetBinCenter(ibin);
		      float y=h->GetBinContent(ibin);
		      
		      float statUnc(h->GetBinError(ibin));
		      float loose2tightUnc(fabs(loose2tightH->GetBinContent(ibin)));
		      float pure2tightUnc(pure2tightGr ? fabs(1-pure2tightGr->Eval(x)) : 0.);  pure2tightUnc *= y;
		      float closurePDFunc(m_gclosurePDFunc);                                   closurePDFunc *= y;
		      
		      float newerr=sqrt(pow(statUnc,2)+pow(loose2tightUnc,2)+pow(pure2tightUnc,2)+pow(closurePDFunc,2));
		      h->SetBinError(ibin,newerr);
		    }
		  }
		}
	      
	      Shape_t m_shape;
	      m_shape.data=h; 
	      m_shape.ewkBkg=(TH1 *)gIn->Get( (gEwkBkgProc+"/"+ch[ich]+cats[icat]+"_"+histos[ih]).c_str() );
	      if(m_shape.ewkBkg){
		m_shape.ewkBkg->SetDirectory(0);
		m_shape.ewkBkg->SetName(TString("mcewkg")+h->GetName());
		m_shape.ewkBkg->SetTitle("EWK #gammajj");
		if(anMode==SEARCH) m_shape.ewkBkg->SetTitle("EWK #gamma+E_{T}^{miss}");
	      }
	      string key=ch[ich]+cats[icat]+"_"+histos[ih];
	      gShapesMap[key]=m_shape;
	    }
	}
    }
    
  //all done with the input files
  llIn->Close();
  gIn->Close();
  //gRawIn->Close();
  if(loosegIn) loosegIn->Close();
  if(closeF) closeF->Close();

  //derive EWK photon scale factors from balance side band
  if(anMode==SEARCH)
    {
      for(size_t ich=0; ich<nchs; ich++)
	{
	  string key=ch[ich]+"_balance";
	  if(gShapesMap.find(key)==gShapesMap.end()) { cout << "Warning: balance histogram is not found, no data-driven scale factor will be derived for photons" << endl;  continue; }

	  Shape_t &gShape=gShapesMap[key];
	  TH1 *data=gShape.data;
	  TH1 *ewk=gShape.ewkBkg;
	  if(data==0 || ewk==0) { cout << "Warning data or ewk photons are null, no data-driven scale factor will be derived for photons" << endl; continue; }
	  Int_t ibin=data->GetXaxis()->FindBin(2);
	  Int_t ebin=data->GetXaxis()->GetNbins();
	  Float_t dataCts(data->Integral(ibin,ebin)), ewkCts(ewk->Integral(ibin,ebin));
	  if(ewkCts==0) { cout << "Warning EWK expected contribution is null, no data-driven scale factor will be derived for photons" << endl; continue; }
	  hzzEWKPhotonSFs[ ch[ich] ]=dataCts/ewkCts;
	  cout <<  "EWK g+jets scale factor for " << ch[ich] << " is " << hzzEWKPhotonSFs[ ch[ich] ] << endl; 
	}
    }

  //
  //PRODUCE FINAL HISTOGRAMS
  //
  TFile *gOut         = TFile::Open("gamma_out.root","RECREATE");
  TDirectory *gOutDir = gOut->mkdir("Instr. background (data)");
  for(size_t ishape=0; ishape<shapesKeys.size(); ishape++)
   {
     std::map<string,Shape_t>::iterator it=shapesMap.find(shapesKeys[ishape]); 
     if(it==shapesMap.end()){
       cout << shapesKeys[ishape] << " can't be found" << endl;
     }
     
     TH1 *corrGammaH=0;
     if(gShapesMap.find(it->first)==gShapesMap.end()) { cout << "[Warning] " << it->first << " is not found in for gammas" << endl; continue; }

     Shape_t &gShape=gShapesMap[it->first];
     corrGammaH=(TH1 *)gShape.data->Clone((it->first+"corrg").c_str());

     //subtract background if available
     TH1 *corrGammaSubUncH=subtractEWKbkg(corrGammaH,gShape.ewkBkg,anMode);

     corrGammaH->SetDirectory(0);
     corrGammaH->SetFillColor(831);     
     //it->second.bckg["Instr. background (data)"]=corrGammaH;
     it->second.bckg["Instr. background (data)"]=corrGammaSubUncH;

     showShape(it->second,"final_"+it->first,anMode);
     gOutDir->cd();

     TString keyToWrite(it->first.c_str());
     corrGammaH->Write(keyToWrite+"_nosubunc");
     if(corrGammaSubUncH) corrGammaSubUncH->Write( keyToWrite );
     
     //build the inclusive shape starting from one of the dilepton channels
     if(it->first.find("mumu")!= string::npos)
       {
	 TString keyToGet(it->first);
	 keyToGet=keyToGet.ReplaceAll("mumu","ee");
	 if(shapesMap.find(keyToGet.Data())!=shapesMap.end()){
	   Shape_t &eeShape=shapesMap[keyToGet.Data()];
	   keyToGet.ReplaceAll("ee","");
	   Shape_t totalShape=cloneShape(it->second,keyToGet);
	   addToShape(totalShape,eeShape);
	   showShape(totalShape,keyToGet.Data(),anMode);
	 }
       }
   }

  gOut->Close();
}

//
TH1 *subtractEWKbkg(TH1 *gDataOrig, TH1 *ewkBkg,int anMode){
  if(ewkBkg==0) return 0;

  TString hname(gDataOrig->GetName());
  TString ch("");
  if(hname.BeginsWith("ee")) ch="ee";
  if(hname.BeginsWith("mumu")) ch="mumu";
  Float_t ewkBkgSF(1.0);
  if(hzzEWKPhotonSFs.find(ch)!=hzzEWKPhotonSFs.end()) ewkBkgSF=hzzEWKPhotonSFs[ch]; 
  ewkBkg->Scale(ewkBkgSF);

  TH1 *gData=(TH1 *)gDataOrig->Clone(gDataOrig->GetName()+TString("_sub"));
  gData->SetDirectory(0);
  gData->SetFillStyle(0);

  //this will include the stat unc of the EWK bkg + 30% uncertainty on the nominal value
  TH1 *gDataSubUnc=(TH1 *)gData->Clone(gData->GetName()+TString("_subunc"));
  gDataSubUnc->Reset("ICE");
  gDataSubUnc->SetDirectory(0);

  bool is2D(gData->InheritsFrom("TH2"));  
  if(is2D){
    cout << gDataOrig->GetName() << " " << gDataOrig->GetTitle() << endl;
    for(int xbin=1; xbin<=gData->GetXaxis()->GetNbins(); xbin++)
      for(int ybin=1; ybin<=gData->GetYaxis()->GetNbins(); ybin++)
	{
	  Double_t gCounts   = gData->GetBinContent(xbin,ybin);
	  Double_t ewkCounts = ewkBkg->GetBinContent(xbin,ybin);
	  Double_t ewkSubErr = sqrt(pow(ewkBkg->GetBinError(xbin,ybin),2)+pow(ewkCounts*0.30,2));
	  Double_t finalGerr = sqrt(pow(ewkSubErr,2)+pow(gData->GetBinError(xbin,ybin),2));

	  gData->SetBinContent(xbin,ybin,TMath::Max(gCounts-ewkCounts,1.e-5));
	  gData->SetBinError(xbin,ybin,ewkSubErr);

	  gDataSubUnc->SetBinContent(xbin,ybin,TMath::Max(gCounts-ewkCounts,1.e-5));
	  gDataSubUnc->SetBinError(xbin,ybin,finalGerr);	  
	}
  }else{
    for(int xbin=1; xbin<=gData->GetXaxis()->GetNbins(); xbin++)
      {
	Double_t gCounts   = gData->GetBinContent(xbin);
	Double_t ewkCounts = ewkBkg->GetBinContent(xbin);
	Double_t ewkSubErr = sqrt(pow(ewkBkg->GetBinError(xbin),2)+pow(ewkCounts*0.30,2));
	Double_t finalGerr = sqrt(pow(ewkSubErr,2)+pow(gData->GetBinError(xbin),2));

	gData->SetBinContent(xbin,TMath::Max(gCounts-ewkCounts,1.e-5));
	gData->SetBinError(xbin,ewkSubErr);
	  
	gDataSubUnc->SetBinContent(xbin,TMath::Max(gCounts-ewkCounts,1.e-5));
	gDataSubUnc->SetBinError(xbin,finalGerr);
      }
  }
  if(anMode==DISCOVERY){
    gDataSubUnc->Scale(gDataOrig->Integral()/gDataSubUnc->Integral());
    gData->Scale(gDataOrig->Integral()/gData->Integral());
  }

  //illustrate subtraction performed
  if(!is2D)
    {
      TCanvas *c=new TCanvas("c","c",600,600);

      gDataOrig->Draw("hist"); 
      gDataOrig->SetFillColor(30);
      gDataOrig->GetYaxis()->SetRangeUser(1e-1,gDataOrig->GetMaximum()*1.25);

      gDataSubUnc->Draw("e2 same"); 
      gDataSubUnc->SetFillStyle(1001); 
      gDataSubUnc->SetFillColor(831);
      gDataSubUnc->SetMarkerStyle(1); 
      gDataSubUnc->SetMarkerColor(33);

      gData->Draw("e2same");
      gData->SetFillColor(30);
      gData->SetMarkerColor(gData->GetFillColor());
      gData->SetMarkerStyle(1);
      gData->SetMarkerColor(30);

      ewkBkg->Draw("hist same");  
      ewkBkg->SetFillStyle(1001); 
      ewkBkg->SetFillColor(804);

      //gDataOrig->SetLineColor(1); 
      //gDataOrig->SetLineStyle(9);
      //gDataOrig->SetFillStyle(0);

      TLegend* leg  = new TLegend(0.15,0.83,0.93,0.92);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->SetNColumns(4);
      leg->AddEntry(gDataOrig,"Raw #gammajj (data)","l");
      leg->AddEntry(gData,"Final #gammajj (data)","f");
      leg->AddEntry(gDataSubUnc,"EWK sub. uncertainty","f");
      leg->AddEntry(ewkBkg,"EWK #gamma jj","f");
      leg->Draw("same");
      drawCMSHeader();

      c->SetLogy();
      c->SaveAs(outDir+"/"+TString(ewkBkg->GetName())+TString("_ewksub.png"));
      delete c;
    }
  
  return gDataSubUnc;
}

//
void showShape(const Shape_t &shape,TString outName,int anMode)
{
  float mjjCen(0);
  TString pname(shape.data->GetName());
  TString mjjCat("ll events, ");
  if(outName.Contains("_ee"))           mjjCat="ee events, ";
  if(outName.Contains("_mumu"))         mjjCat="#mu#mu events, ";
  if(pname.Contains("mjjq016"))          { mjjCen=250*0.5;        mjjCat+="M_{jj}<250 GeV";     }
  else if(pname.Contains("mjjq033"))     { mjjCen=(350+250)*0.5;  mjjCat+="250<M_{jj}<350 GeV"; }
  else if(pname.Contains("mjjq049"))     { mjjCen=(350+450)*0.5;  mjjCat+="350<M_{jj}<450 GeV"; }
  else if(pname.Contains("mjjq066"))     { mjjCen=(450+550)*0.5;  mjjCat+="450<M_{jj}<550 GeV"; }
  else if(pname.Contains("mjjq083"))     { mjjCen=(750+550)*0.5;  mjjCat+="550<M_{jj}<750 GeV"; }
  else if(pname.Contains("mjjq092"))     { mjjCen=(1000+750)*0.5; mjjCat+="750<M_{jj}<1000 GeV"; }
  else if(pname.Contains("mjjgt092"))    { mjjCen=1200;           mjjCat+="M_{jj}>750 GeV"; }
  else if(pname.Contains("highmjj"))     { mjjCen=1700;           mjjCat+="M_{jj}>1250 GeV"; }
  else if(pname.Contains("lowhardpt"))   { mjjCen=0;              mjjCat+="#scale[0.5]{Hard p_{T}}<50 GeV"; }
  else if(pname.Contains("highhardpt"))  { mjjCen=0;              mjjCat+="#scale[0.5]{Hard p_{T}}>50 GeV"; }
  else if(pname.Contains("eq0jets"))     { mjjCen=0;              mjjCat="=0 jets"; }
  else if(pname.Contains("geq1jets"))    { mjjCen=0;              mjjCat="#geq1 jets"; }
  else if(pname.Contains("vbf"))         { mjjCen=0;              mjjCat="VBF"; }
  else if(anMode==DISCOVERY)             { mjjCen=0;              mjjCat+="inclusive M_{jj}"; }
  else                                   { mjjCen=0;              mjjCat+="inclusive";        }

  bool is2D(shape.data->InheritsFrom("TH2"));
  if(is2D)  cout << "2D will only be profiled: " << shape.data->GetName() << endl;
  if(outName.BeginsWith("_")) outName="inc"+outName;

  TCanvas* c1 = new TCanvas(outName,outName,600,600);
  c1->cd();

  TPad* t1 = new TPad("t1","t1", 0.0, 0.20, 1.0, 1.0);  t1->Draw();  t1->cd();
  t1->SetTopMargin(0.2); 
  if(outName.Contains("Likelihood") || outName.Contains("Fisher"))// || outName.Contains("BDTD") || outName.Contains("MLP") )
    if(!outName.Contains("092") && !outName.Contains("100") && !outName.Contains("highmjj"))
      t1->SetLogy(true);

  TLegend* leg  = new TLegend(0.15,0.83,0.93,0.92);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  int ncols(shape.bckg.size()>5 ? shape.bckg.size()/2 : shape.bckg.size()+2);
  leg->SetNColumns(ncols);
  THStack *stack = new THStack("stack","stack");
  
  //main signal is the first
  TH1 *hsignal=shape.signal.begin()->second;
  
  //start with the backgrounds
  TH1 *instrBckg=0;
  for(std::map<TString,TH1 *>::const_reverse_iterator it = shape.bckg.rbegin(); it != shape.bckg.rend(); it++)
    {
      TString itit(it->second->GetTitle());
      if(anMode==DISCOVERY && itit.Contains("QCD Zjj")) { instrBckg=it->second; continue; }
      else if (itit.Contains("Z#rightarrow ll (data)")) { instrBckg=it->second; continue; }
      stack->Add(it->second,"HIST");
      leg->AddEntry(it->second,it->second->GetTitle(),"F");
    }

  TH1 *totalPredMCBkg=(TH1 *)shape.data->Clone(TString("totalmcbkg_")+shape.data->GetName());
  totalPredMCBkg->Reset("ICE");
  totalPredMCBkg->SetDirectory(0);
  for(std::map<TString,TH1 *>::const_reverse_iterator it = shape.mcbckg.rbegin(); it != shape.mcbckg.rend(); it++)
    {
      TString itit(it->second->GetTitle());
      totalPredMCBkg->Add(it->second);
    }

  //add the instrumental background 
  //  apply correction factor to the prediction to match data counts assuming signal strength of 1
  //  data=others+k.dy <=> k=(data-others)/dy
  float totalData(0), totalInstr(0), totalOthers(0),mcdy(0);
  float sfdy(1.0),sfmcdy(1.0);
  if(anMode==DISCOVERY)
    {
      TH1 *totalPred = (TH1 *)stack->GetStack()->At( stack->GetStack()->GetEntriesFast()-1 );
      if(!is2D)
	{
	  totalData=shape.data->Integral("width");
	  totalInstr=instrBckg->Integral("width");
	  totalOthers=hsignal->Integral("width")+totalPred->Integral("width");
	  mcdy=shape.mcbckg.find("Z#rightarrow ll")->second->Integral("width");
	}
      else
	{
	  Int_t nybins=((TH2 *)shape.data)->GetYaxis()->GetNbins();
	  totalData=((TH2 *)shape.data)->Integral();//1,1,1,nybins);
	  totalInstr=((TH2 *)instrBckg)->Integral();//1,1,1,nybins);
	  totalOthers=((TH2 *)hsignal)->Integral()-totalPred->Integral();//1,1,1,nybins)+totalPred->Integral(1,1,1,nybints);
	  mcdy=((TH2 *)shape.mcbckg.find("Z#rightarrow ll")->second)->Integral();//1,10,1,nybins);
	}
      
      sfdy=( totalInstr > 0 ? (totalData-totalOthers)/totalInstr : 1.0);
      sfmcdy=(totalInstr>0 ? mcdy/totalInstr : 1.0);
      //      if(mjjWgtGr && mjjCen!=0) { sfmcdy *= mjjWgtGr->Eval(mjjCen);  }
    }
  else{

    //for discovery mode we normalize from a sideband
    std::pair<TString,TString> key("","");
    if(outName.Contains("mumu"))          key.first="mumu";
    else if(outName.Contains("ee"))       key.first="ee";
    if(outName.Contains("eq0jets"))       key.second="eq0jets";
    else if(outName.Contains("geq1jets")) key.second="geq1jets";
    else if(outName.Contains("vbf"))      key.second="vbf";
    
    if( outName.Contains("_met") && !outName.Contains("NM1") && !outName.Contains("shapes") && !outName.Contains("count")) {
      //if( outName.Contains("mindphijmet") && !outName.Contains("NM1") ) {
      //if( outName.Contains("balance") && !outName.Contains("NM1") ) {
      //Int_t maxBin=shape.data->GetXaxis()->FindBin(0.5);
      Int_t maxBin=shape.data->GetXaxis()->FindBin(50);
      if(key.second=="eq0jets") maxBin=shape.data->GetXaxis()->GetNbins();
      totalData=shape.data->Integral(1,maxBin,"width");
      totalInstr=instrBckg->Integral(1,maxBin,"width");
      sfdy=totalData/totalInstr;
      hzzSFs[key]=sfdy;
      cout << maxBin << " " << totalData << " " << totalInstr << " " << sfdy << endl;
    }
    sfdy=hzzSFs[key];
    sfmcdy=hzzSFs[key];
  }

  //if(outName.Contains("_shapes")) { sfdy=0.5*(1.0+sfdy); sfmcdy=0.5*(1.0+sfmcdy); }
  cout << outName << " " << " DD:" << sfdy << " MC:" << sfmcdy << endl;  
  cout << instrBckg << " " << stack << endl;
  instrBckg->Scale(sfdy);
  //instrBckg->Scale(sfmcdy);
  stack->Add(instrBckg,"HIST");
  leg->AddEntry(instrBckg,instrBckg->GetTitle(),"F");

  //add the signal
  if(anMode==DISCOVERY){
    hsignal->SetFillColor(804);
    hsignal->SetLineStyle(1);
    hsignal->SetLineColor(1);
    stack->Add(hsignal,"HIST");
    leg->AddEntry(hsignal,hsignal->GetTitle(),"F");
  }
  cout << stack->GetStack()->GetEntriesFast() << endl;
  //update total prediction
  TH1 *totalPredBkg=(TH1 *) (stack->GetStack()->At( stack->GetStack()->GetEntriesFast()-1 )->Clone(TString("totalbkg_")+shape.data->GetName()) );
  totalPredBkg->SetDirectory(0);

  TH1 *totalPred=(TH1 *) totalPredBkg->Clone( TString("totalbkgps_")+shape.data->GetName() );
  cout << hsignal << endl;
  // totalPred->Add(hsignal);

  cout << __LINE__ << endl;
  if(!is2D)
    {
      //draw the stack
      TH1D *frame=(TH1D *) shape.data->Clone("frame");
      frame->Reset("ICE");       
      frame->Draw();
      frame->GetXaxis()->SetTitle(totalPred->GetXaxis()->GetTitle());
      frame->GetYaxis()->SetTitle(totalPred->GetYaxis()->GetTitle());
      frame->SetMinimum(0.01);
      frame->SetMaximum(1.3*max(totalPred->GetMaximum(),shape.data->GetMaximum()));
      if(anMode==SEARCH)
	{
	  if(!outName.Contains("mindphijmet")) t1->SetLogy(true);
	  if((outName.Contains("_met") && !outName.Contains("count")) || outName.Contains("_mt") || outName.Contains("mjj")) {
	    t1->SetLogx(true);
	    frame->GetXaxis()->SetRangeUser(frame->GetXaxis()->GetBinCenter(2),frame->GetXaxis()->GetXmax());
	  }
	}
      
  cout << __LINE__ << endl;
      stack->Draw("histsame");
        cout << __LINE__ << endl;
      //draw the uncertainty band
      //      TGraphAsymmErrors *mcgr=new TGraphAsymmErrors(anMode==DISCOVERY ? totalPred : totalPredBkg);
      TGraphAsymmErrors *mcgr=new TGraphAsymmErrors(totalPredBkg);
      mcgr->SetFillStyle(3001);//3427);
      mcgr->SetFillColor(kGray+1);
      mcgr->SetMarkerStyle(1);
      mcgr->Draw("e2p");
        cout << __LINE__ << endl;
      //draw the data
      shape.data->Draw("e1same");
      leg->AddEntry(shape.data,shape.data->GetTitle(),"P");
        cout << __LINE__ << endl;
      //superimpose the signals 
      for(std::map<TString , TH1 *>::const_iterator it=shape.signal.begin(); it!= shape.signal.end(); it++)
	{
	  TH1 *h=(TH1 *)it->second->Clone();
	  h->SetDirectory(0);
	  h->SetFillStyle(0);
	  if(anMode==SEARCH) h->SetLineColor(1);
	  else h->SetLineColor(h->GetFillColor());
	  h->SetLineWidth(3);
	  h->Draw("histsame");
	  if(it!=shape.signal.begin() || anMode==SEARCH) leg->AddEntry(h,h->GetTitle(),"L");
	}
        cout << __LINE__ << endl;
      //caption
      drawCMSHeader();
      leg->Draw("same");  
      
      //compute the chi2 and K-S test for the total prediction
      TPaveText *pave = new TPaveText(0.65,0.955,0.975,0.97,"NDC");
      pave->SetBorderSize(0);
      pave->SetFillStyle(0);
      pave->SetTextFont(42);
      pave->SetTextAlign(32);
      pave->SetTextSize(0.035);
      char buf[1000];
      //sprintf(buf,"#splitline{%s}{#chi^{2}/ndof : %3.2f , K-S prob : %3.2f}",mjjCat.Data(), shape.data->Chi2Test(totalPred,"UWCHI2/NDF"),shape.data->KolmogorovTest(totalPred,"") );
      //pave->AddText(buf);
      pave->AddText(mjjCat);
      pave->Draw();
        cout << __LINE__ << endl;
      //ratio canvas
      c1->cd();
      TPad* t2 = new TPad("t2","t2", 0.0, 0.0, 1.0, 0.2);    
      t2->Draw();
      t2->cd();  
      t2->SetTopMargin(0); 
      t2->SetBottomMargin(0.2);
  cout << __LINE__ << endl;
      TH1 *ratio = (TH1*)shape.data->Clone("RatioHistogram");
      ratio->SetDirectory(0);
      if(anMode==DISCOVERY) ratio->Divide((TH1 *)stack->GetStack()->At( stack->GetStack()->GetEntriesFast()-1 ) );
      else                  ratio->Divide(totalPredBkg);
        cout << __LINE__ << endl;
      TGraphAsymmErrors *denRelUnc=new TGraphAsymmErrors;
      denRelUnc->SetLineColor(1);
      denRelUnc->SetFillStyle(3001);
      denRelUnc->SetFillColor(kGray);
      denRelUnc->SetMarkerColor(1);
      denRelUnc->SetMarkerStyle(1);
      for(int ip=0; ip<mcgr->GetN(); ip++)
	{
	  Double_t x,y;
	  mcgr->GetPoint(ip,x,y);
	  Double_t xLo=mcgr->GetErrorXlow(ip);
	  Double_t xHi=mcgr->GetErrorXhigh(ip);
	  Double_t yLo=y!=0?mcgr->GetErrorYlow(ip)/y:0.;
	  Double_t yHi=y!=0?mcgr->GetErrorYhigh(ip)/y:0.;
	  denRelUnc->SetPoint(ip,x,1.0);
	  denRelUnc->SetPointError(ip,xLo,xHi,yLo,yHi);
	}
  cout << __LINE__ << endl;
      TH1D *ratioFrame=(TH1D *) ratio->Clone("ratioframe");
      ratioFrame->Reset("ICE");       
      ratioFrame->Draw();
      if(anMode==SEARCH && ((outName.Contains("_met") && !outName.Contains("count")) || outName.Contains("_mt") || outName.Contains("mjj")) ) {
	t2->SetLogx(true);
	ratioFrame->GetXaxis()->SetRangeUser(frame->GetXaxis()->GetBinCenter(2),ratioFrame->GetXaxis()->GetXmax());
      }
        cout << __LINE__ << endl;
      float yscale = (1.0-0.2)/(0.18-0);       
      ratioFrame->GetYaxis()->SetTitle("Data/#Sigma Bckg");
      ratioFrame->SetMinimum(0.7);
      ratioFrame->SetMaximum(1.3);
      ratioFrame->GetXaxis()->SetTitle("");
      ratioFrame->GetXaxis()->SetTitleOffset(1.3);
      ratioFrame->GetXaxis()->SetLabelSize(0.033*yscale);
      ratioFrame->GetXaxis()->SetTitleSize(0.036*yscale);
      ratioFrame->GetXaxis()->SetTickLength(0.03*yscale);
      ratioFrame->GetYaxis()->SetTitleOffset(0.3);
      ratioFrame->GetYaxis()->SetNdivisions(5);
      ratioFrame->GetYaxis()->SetLabelSize(0.033*yscale);
      ratioFrame->GetYaxis()->SetTitleSize(0.036*yscale);
      denRelUnc->Draw("3");
      ratio->Draw("e1 same");
        cout << __LINE__ << endl;
      c1->cd();
      c1->Modified();
      c1->Update();

      c1->SaveAs(outDir+"/"+outName+".png");
      c1->SaveAs(outDir+"/"+outName+".pdf");
      c1->SaveAs(outDir+"/"+outName+".C");
  cout << __LINE__ << endl;
      //draw sub-pad with subtracted data
      t1->cd();
      TPad *t11=new TPad("t11","t11",0.67,0.48,0.94,0.83);
      t11->Draw();
      t11->cd();
      t11->SetFillStyle(0);
      t11->SetBottomMargin(0.25);
      t11->SetLeftMargin(0.25);
      
      TH1 *resData=(TH1 *)shape.data->Clone("residualdata");
      resData->SetDirectory(0);
      resData->Add(totalPredBkg,-1);

      //draw first the signal
      int isig(0);
      for(std::map<TString , TH1 *>::const_iterator it=shape.signal.begin(); it!= shape.signal.end(); it++,isig++)
	{
	  TH1 *h=it->second;
	  if(isig==0) h=(TH1 *)h->Clone();
	  h->Draw(isig==0 ? "hist" : "histsame");
	  
	  if(isig>0) continue;
	  //h->GetXaxis()->SetTitle("");
	  h->GetXaxis()->SetLabelSize(0.08);
	  h->GetYaxis()->SetTitle("Data-#Sigma Bckg");
	  h->GetYaxis()->SetTitleSize(0.08);
	  h->GetXaxis()->SetTitleSize(0.08);
	  h->GetYaxis()->SetLabelSize(0.08);
	  h->GetYaxis()->SetRangeUser(resData->GetMinimum()*1.1,resData->GetMaximum()*1.2);
	}
      
      //now the data
      resData->Draw("e1same");
      
      c1->cd();
      c1->Modified();
      c1->Update();

      c1->SaveAs(outDir+"/"+outName+"_sub.png");
      c1->SaveAs(outDir+"/"+outName+"_sub.pdf");
      c1->SaveAs(outDir+"/"+outName+"_sub.C");
    }
  
  //  if(drawEvol)
  drawEvolutionFor(shape.data,totalPredBkg,totalPredMCBkg,shape.signal,outDir+"/"+outName+"_evol");
}


//
void drawEvolutionFor(TH1 *data, TH1 *totalPred, TH1 *mcPred, const std::map<TString,TH1 *>&sig, TString SaveName)
{
  TCanvas *c=new TCanvas(SaveName,SaveName,600,600);
  TPad *p=new TPad("p1","p1",0,0.2,1,1);
  p->Draw();
  p->cd();
  
  bool isTH2evol(data->InheritsFrom("TH2") );

  TLegend *leg = isTH2evol ? new TLegend(0.15,0.7,0.5,0.92) : new TLegend(0.6,0.1,0.95,0.4);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  //  if(!isTH2evol) leg->SetNColumns(3);

  //plot the signal
  /*
  int i(0);
  std::vector<TGraphErrors *> sig_cdf;
  for(std::map<TString, TH1 *>::const_iterator it = sig.begin(); it!=sig.end(); it++,i++){
    TGraphErrors *cdf=getEvolutionWithErrors(it->second);
    cdf->SetTitle(it->first);
    cdf->SetFillStyle(1001);
    cdf->SetMarkerColor(804-i*2);
    cdf->SetFillColor(804-i*2);
    cdf->SetLineColor(804-1*2);
    sig_cdf.push_back(cdf);
    cdf->Draw(i==0?"al3":"l3");
    if(i>0) continue;
    if(!isTH2evol ) cdf->GetYaxis()->SetTitle("Gap fraction");
    else           cdf->GetYaxis()->SetTitle( "<" + TString( it->second->GetYaxis()->GetTitle() ) +">" );
    cdf->GetXaxis()->SetTitle(it->second->GetXaxis()->GetTitle());
    //cdf->GetXaxis()->SetTitleSize(0);
    //cdf->GetXaxis()->SetLabelSize(0);
    //cdf->GetYaxis()->SetTitleSize(0.06);
    //cdf->GetYaxis()->SetLabelSize(0.04);
    //cdf->GetYaxis()->SetTitleOffset(0.75);
    //cdf->GetYaxis()->SetRangeUser(0.04,1.04);
    leg->AddEntry(cdf,cdf->GetTitle(),"f");
  }
  */
  //plot the mc based prediction
  TGraphErrors *mcPred_cdf=0;
  if(mcPred){
    mcPred_cdf=getEvolutionWithErrors(mcPred);
    if(!isTH2evol ) mcPred_cdf->GetYaxis()->SetTitle("Gap fraction");
    else            mcPred_cdf->GetYaxis()->SetTitle( "<" + TString( mcPred->GetYaxis()->GetTitle() ) +">" );
    mcPred_cdf->GetXaxis()->SetTitle(mcPred->GetXaxis()->GetTitle());   
    mcPred_cdf->SetFillStyle(1001);
    mcPred_cdf->SetFillColor(33);
    mcPred_cdf->SetLineColor(33);
    mcPred_cdf->SetTitle("Total expected (MC)");
    mcPred_cdf->Draw("al3");
    leg->AddEntry(mcPred_cdf,mcPred_cdf->GetTitle(),"f");
  }

  //total prediction
  TGraphErrors *totalPred_cdf=getEvolutionWithErrors(totalPred);
  totalPred_cdf->SetMarkerStyle(24);
  totalPred_cdf->SetTitle("Total expectations (data)");
  totalPred_cdf->Draw("p");
  leg->AddEntry(totalPred_cdf,totalPred_cdf->GetTitle(),"p");

  TGraphErrors *data_cdf=getEvolutionWithErrors(data);
  data_cdf->SetMarkerStyle(20);
  data_cdf->SetTitle("data");
  data_cdf->Draw("p");
  leg->AddEntry(data_cdf,data_cdf->GetTitle(),"p");

  //caption
  drawCMSHeader();
  leg->Draw();


      //compute the chi2 and K-S test for the total prediction
  TPaveText *pave = new TPaveText(0.65,0.955,0.975,0.985,"NDC");
  pave->SetBorderSize(0);
  pave->SetFillStyle(0);
  pave->SetTextFont(42);
  pave->SetTextAlign(32);
  TString pname(data->GetName());
  TString mjjCat("ee/#mu#mu events, ");
  if(SaveName.BeginsWith("ee"))        mjjCat="ee events, ";
  if(SaveName.BeginsWith("mumu"))      mjjCat="#mu#mu events, ";
  if(pname.Contains("mjjq016"))       mjjCat+="M_{jj}<250 GeV";
  else if(pname.Contains("mjjq033"))  mjjCat+="250<M_{jj}<350 GeV";
  else if(pname.Contains("mjjq049"))  mjjCat+="350<M_{jj}<450 GeV"; 
  else if(pname.Contains("mjjq066"))  mjjCat+="450<M_{jj}<550 GeV"; 
  else if(pname.Contains("mjjq083"))  mjjCat+="550<M_{jj}<750 GeV"; 
  else if(pname.Contains("mjjq092"))  mjjCat+="750<M_{jj}<1000 GeV"; 
  else if(pname.Contains("mjjgt092")) mjjCat+="M_{jj}750 GeV"; 
  else if(pname.Contains("highmjj"))  mjjCat+="M_{jj}>1250 GeV"; 
  else                                mjjCat+="inclusive M_{jj}"; 
  pave->AddText(mjjCat);
  pave->Draw();


  //compare data to prediction(s)
  c->cd();
  p=new TPad("p2","p2",0,0,1,0.2);
  p->Draw();
  p->cd();
  p->SetTopMargin(0);
  p->SetBottomMargin(0.2);

  TGraphErrors *gr=new TGraphErrors;
  gr->SetName("data2pred");
  for(int ip=0; ip<data_cdf->GetN(); ip++)
    {
      Double_t x,yn;
      data_cdf->GetPoint(ip,x,yn); 
      Double_t xerr=data_cdf->GetErrorX(ip); 
      Double_t ynerr=data_cdf->GetErrorY(ip);
      Double_t yd;
      totalPred_cdf->GetPoint(ip,x,yd); 
      Double_t yderr=totalPred_cdf->GetErrorY(ip); 
      if(yd==0) continue;
      int np=gr->GetN();
      gr->SetPoint(np, x, yn/yd);
      gr->SetPointError(np,xerr,sqrt(pow(ynerr*yd,2)+pow(yderr*yn,2))/(yd*yd));
    }
  gr->SetMarkerStyle(24);
  gr->Draw("ap");
  float yscale = (1.0-0.2)/(0.18-0);
  gr->GetYaxis()->SetTitle("Data/#Sigma Bckg");
  gr->SetMinimum(0.7);
  gr->SetMaximum(1.3);
  gr->GetXaxis()->SetTitle("");
  gr->GetXaxis()->SetTitleOffset(1.3);
  gr->GetXaxis()->SetLabelSize(0.033*yscale);
  gr->GetXaxis()->SetTitleSize(0.036*yscale);
  gr->GetXaxis()->SetTickLength(0.03*yscale);
  gr->GetYaxis()->SetTitleOffset(0.3);
  gr->GetYaxis()->SetNdivisions(5);
  gr->GetYaxis()->SetLabelSize(0.033*yscale);
  gr->GetYaxis()->SetTitleSize(0.036*yscale);

  c->cd();
  c->SaveAs(SaveName+".png");
  c->SaveAs(SaveName+".pdf");
  c->SaveAs(SaveName+".C");
}

//
TGraphErrors *getEvolutionWithErrors(TH1 *h)
{
  TGraphErrors *gr=new TGraphErrors;
  gr->SetName(h->GetName()+TString("_evol"));
  
  bool isTH2( h->InheritsFrom("TH2") );

  Int_t nbins=h->GetXaxis()->GetNbins();
  for(Int_t ibin=1; ibin<=nbins; ibin++)
    {
      Double_t x=h->GetXaxis()->GetBinCenter(ibin);
      Double_t xerr=h->GetXaxis()->GetBinWidth(ibin)/2;
      Double_t val(0),val_err(0);
        
      if(!isTH2)
	{
	  //for TH1 we compute the gap fraction
	  Double_t total=h->Integral(1,nbins);
	  val=h->IntegralAndError(1,ibin,val_err)/total;
	  val_err/=total;
	}
      else
	{    
	  //for TH2 we project the average y as function of x
	  TH2* h2d=(TH2 *)h;
	  TH1 *py=h2d->ProjectionY(TString(h->GetName())+"_py",ibin,ibin);
	  if(py->GetEntries()<2) continue; 
	  val=py->GetMean();
	  val_err=py->GetMeanError();
	  delete py;
	}
      
      if(val==0) continue;
      Int_t ip=gr->GetN();
      gr->SetPoint(ip,x,val);
      gr->SetPointError(ip,xerr,val_err);
    }
  
  return gr;
}



void setTDRStyle() {
  tdrStyle = new TStyle("tdrStyle","Style for P-TDR");
  
  // For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

  // For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);
  
  // For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);
  
  // For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);
  
  tdrStyle->SetEndErrorSize(2);
  // tdrStyle->SetErrorMarker(20);
  tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);
  
  //For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  //For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);
  
  // For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);
  
  // Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);
  
  // For the Global title:
  
  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);
  
  // For the axis titles:
  
  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.25);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset
  
  // For the axis labels:
  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");
  
  // For the axis:
  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);
  
  // Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);
  
  // Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);
  
  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);
  
  tdrStyle->cd();
  
  gStyle->SetErrorX(0.5);
  gStyle->SetPadTopMargin   (0.06);
  gStyle->SetPadBottomMargin(0.12);
  //  gStyle->SetPadRightMargin (0.16);
  tdrStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin  (0.14);
  gStyle->SetTitleSize(0.04, "XYZ");
  gStyle->SetTitleXOffset(1.1);
  gStyle->SetTitleYOffset(1.45);
  gStyle->SetPalette(1);
  gStyle->SetNdivisions(505);
}
