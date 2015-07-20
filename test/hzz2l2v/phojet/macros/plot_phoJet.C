#include "setTDRStyle.C"

// Exponentially modified Gaussian
Double_t expgaus(Double_t *x, Double_t *par) {

      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      Double_t np = 50.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;
      bool x_initial = false;

      // MP shift correction
      //   mpc = par[1] - mpshift * par[0]; 

      // Range of convolution integral
      //      xlow = x[0] - sc * par[2];
      //xupp = x[0] + sc * par[2];
      xlow = par[4];
      xupp=x[0];
      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;

         //cout << xx << endl;
	 fland = (par[3]/2.)*TMath::Exp( (par[3]/2.)*(2.*par[1]+par[3]*par[2]*par[2]-2.*xx) ); 
         sum += fland * TMath::Erfc( (par[1]+par[3]*par[2]*par[2]-xx)/(sqrt(2.)*par[2]) );
	   //errorFun(xx,&par[5]); //TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;


         //cout << xx << endl;
         //cout << "=========" << endl;
	 fland = (par[3]/2.)*TMath::Exp( (par[3]/2.)*(2.*par[1]+par[3]*par[2]*par[2]-2.*xx) ); 
         sum += fland * TMath::Erfc( (par[1]+par[3]*par[2]*par[2]-xx)/(sqrt(2.)*par[2]) );
         //cout << "Par 3  " << par[3] << " x value " << *x << endl;
      }



      return (par[0] *step * sum ); //* invsq2pi / par[3]);
}


TH1D* readHist(TString nameHist,TString nameFile, int rebin)
{
  TFile* file = new TFile(nameFile,"read");
  file->ls();
 if (!(file->GetListOfKeys()->Contains(nameHist)) ) return 0;

 TH1D* hist = (TH1D*)file->Get(nameHist);
 // if (hist==0) return;
 hist->GetSumw2();
 // hist->SetLineWidth(2);
 if(rebin>0) hist->Rebin(rebin);
 hist->GetXaxis()->SetTitleSize(.055);
 hist->GetYaxis()->SetTitleSize(.055);
 hist->GetXaxis()->SetLabelSize(.05);
 hist->GetYaxis()->SetLabelSize(.05);
 hist->SetStats(kFALSE);
 return hist;
}

TH2D* readHist2D(TString nameHist,TString nameFile, int rebin)
{
  TFile* file = new TFile(nameFile,"read");

  if (!(file->GetListOfKeys()->Contains(nameHist)) ) return 0;

 TH2D* hist = (TH2D*)file->Get(nameHist);
 hist->GetSumw2();
 // hist->SetLineWidth(2);
 // if(rebin>0) hist->RebinX(rebin); hist->RebinY(rebin);
 hist->GetXaxis()->SetTitleSize(.055);
 hist->GetYaxis()->SetTitleSize(.055);
 hist->GetXaxis()->SetLabelSize(.05);
 hist->GetYaxis()->SetLabelSize(.05);
 hist->SetStats(kFALSE);

 return hist;
}

TH3D* readHist3D(TString nameHist,TString nameFile, int rebin)
{

 TFile* file = new TFile(nameFile);

 TH3D* hist = (TH3D*)file->Get(nameHist);
 hist->GetSumw2();
 // hist->SetLineWidth(2);
 // if(rebin>0) hist->RebinX(rebin); hist->RebinY(rebin);
 hist->GetXaxis()->SetTitleSize(.055);
 hist->GetYaxis()->SetTitleSize(.055);
 hist->GetXaxis()->SetLabelSize(.05);
 hist->GetYaxis()->SetLabelSize(.05);
 hist->SetStats(kFALSE);
 return hist;
}

TCanvas* getaCanvas(TString name)
{

  TCanvas* aCanvas = new TCanvas(name,"",292,55,500,700);//,"",181,237,1575,492);

  aCanvas->SetFillColor(0);
  aCanvas->SetBottomMargin(0.125);
  aCanvas->SetLeftMargin(0.125);
  aCanvas->SetFrameFillColor(0);
  aCanvas->SetFrameBorderMode(0);
  aCanvas->SetFrameLineWidth(2);
  return aCanvas;
}

TCanvas* getanotherCanvas(TString name)
{

  TCanvas* aCanvas = new TCanvas(name,"",292,55,700,465);//,"",181,237,1575,492);

  aCanvas->SetFillColor(0);
  aCanvas->SetBottomMargin(0.125);
  aCanvas->SetLeftMargin(0.125);
  aCanvas->SetFrameFillColor(0);
  aCanvas->SetFrameBorderMode(0);
  aCanvas->SetFrameLineWidth(2);
  return aCanvas;
}

TLegend *legend() {

 TLegend *leg2 = new TLegend(0.52,0.67,0.92,0.90);
 leg2->SetFillStyle(0);
 leg2->SetBorderSize(0);
 leg2->SetTextSize(0.03);
 leg2->SetTextFont(42);

 return leg2;

}

TPad *getaPad_up(TString name){

  TPad *pad1 = new TPad("name", "The pad with the function",0.05,0.4,0.95,0.1);
  //   pad1->Draw();

   //   pad1->Range(-112.6742,-73.17708,1143.438,551.3021);
   pad1->SetFillColor(0);
   pad1->SetBorderMode(0);
   pad1->SetBorderSize(2);
   pad1->SetGridx();
   pad1->SetGridy();
   pad1->SetLeftMargin(0.1271439);
   pad1->SetRightMargin(0.07307979);
   pad1->SetTopMargin(0.08215179);
   pad1->SetBottomMargin(0.117181);
   pad1->SetFrameBorderMode(0);
   pad1->SetFrameBorderMode(0);
  
   return pad1;

}

TPad *getaPad_dn(TString name){

  
 TPad *pad2 = new TPad("pad2", "The pad with the histogram",0.05,0.95,0.95,0.5);
 //   pad2->Draw();

   //   pad2->Range(-153.3652,-2.142584,1185.427,5.367464);
   pad2->SetFillColor(0);
   pad2->SetBorderMode(0);
   pad2->SetBorderSize(2);
   pad2->SetGridx();
   pad2->SetGridy();
   pad2->SetLeftMargin(0.1215511);
   pad2->SetRightMargin(0.07867263);
   pad2->SetTopMargin(0.04892967);
   pad2->SetBottomMargin(0.1521407);
   pad2->SetFrameBorderMode(0);
   pad2->SetFrameBorderMode(0);

   return pad2;

}

TGraphAsymmErrors* drawEff(bool fitcurve, TString hname1, TString hname2, TString filename, TString header, TString xtitle, int rebin, bool asymBin, float xlow, float xhigh, float xmin, float xmax, int icol, int imark, TString draw, double mup){//,double y, double c) {

  TH1D *ref = readHist(hname1, filename, rebin);
  TH1D *sel = readHist(hname2, filename, rebin);

  ref->SetName(hname1+header);
  sel->SetName(hname2+header);

  TGraphAsymmErrors *Eff = new TGraphAsymmErrors();
  if (asymBin) {
    Double_t xbins[16]={0.,22.,34.,42.,52.,65.,75.,88.,100.,150.,210.,280.,360.,500.,700.,1000.};

    sel->Rebin(16,"nsel",xbins);
    ref->Rebin(16,"nref",xbins);
    Eff->BayesDivide(nsel, nref);
  } else {
    Eff->BayesDivide(sel, ref);
    //Eff->Divide(sel,ref,"cl=0.683 b(1,1) mode");
    //Eff->Divide(sel,ref,"cl=0.683 b(1,1) mode v n");

  }

  Double_t Nbins=(double)ref->GetNbinsX();
  Double_t Nconv=(xhigh-xlow)/100.;

  if (Nbins>=Nconv) {
    cout << "Histo bins >= Nconvolution steps; OK" << endl;
  } else {
    cout << "Histo bins < Nconvolution steps!!! Not OK" << endl;
  }


  // gStyle->SetOptStat(1);
  gStyle->SetOptFit(0);

  TF1 *fermiFunction = new TF1("fermiFunction",expgaus,xlow,xhigh,5);

  Double_t params[5] = {0.6,mup,10.,1.,xlow};
  fermiFunction->SetParameters(params);
  fermiFunction->SetParNames("#epsilon","#mu","#sigma","#lambda","xlow");
  if (hname2 != "L1JetAnalysis/SumEt60"  ){ 
    fermiFunction->FixParameter(4,xlow);
    // fermiFunction->SetParLimits(1,100.,150.);
  }
  //  fermiFunction->SetParLimits(1,100.,150.);
  fermiFunction->SetParLimits(0,0.0,1.0);
  fermiFunction->SetLineColor(icol);
  fermiFunction->SetLineWidth(3);

  if (fitcurve) {
    for(int i=0; i != 1 ; i++){  
      Eff->Fit("fermiFunction","RS+");
    }
  }

  Eff->SetLineWidth(2);
  Eff->SetMarkerStyle(imark); Eff->SetMarkerSize(0.5);
  Eff->SetMarkerColor(icol);
  Eff->SetLineColor(icol);
  Eff->GetXaxis()->SetTitle(xtitle);
  Eff->GetYaxis()->SetTitle("Trigger Efficiency");
  // Eff->GetXaxis()->SetLabelSize(0.03);
  // Eff->GetYaxis()->SetLabelSize(0.03);
  Eff->GetXaxis()->SetTitleSize(0.05);
  Eff->GetYaxis()->SetTitleSize(0.05);
  Eff->GetXaxis()->SetTitleOffset(1.);

  Eff->GetYaxis()->SetRangeUser(0.,1.1);
  Eff->GetXaxis()->SetRangeUser((xmin+0.01),xmax);
  Eff->GetXaxis()->SetMoreLogLabels();

  Double_t xpoint;
  Double_t ypoint;

  TF1 *fit = Eff->GetFunction("fermiFunction");
  if (draw=="") {
    Eff->Draw("AP");
  } else {
    Eff->Draw("P");
  }
  
  return Eff;

}

// Scale Histogram errors according to prescale: err'~ sqrt(prescale) ?????
TH1D *hTemp(TH1D* h, Int_t prescale) { 

  TH1D *h_temp=h->Clone();

  Double_t err1, err2;

  Int_t nB=h->GetNbinsX();
  for (Int_t i=0; i<nB; i++) {

    err1=h->GetBinError(i);
    // if (err1!=0) {cout << "Bin error = " << err1 << endl;}
    err2=err1*sqrt(prescale);
    h_temp->SetBinError(i,err2);
    // if (err1!=0) {
    //   cout << "Bin error after prescaling of = " << prescale << ", is = " << err2 << endl;
    // }
  }

  return h_temp;
}

void plot_phojet(TString tag="eevbf", Bool_t prescale=false) {
  /* tag string can be:
     tag="sel" ; trigger before any offline cut (just 1 offline photon)
     tag="sel_met" ; 1photon + pfMET>100 offline
     tag="vbf"
  */
  TString ifile="plotter.root";

  const Int_t nT=16;

  TString HLT_trigName[nT]={"HLT_Photon160_PFMET40_v1",
			    "HLT_Photon150_PFMET40_v1",
			    "HLT_Photon135_PFMET40_v1",
			    "HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_PFMET40_v1",
			    "HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_PFMET40_v1",
			    "HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_PFMET40_v1",
			    "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PFMET40_v1",
			    "HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_PFMET40_v1",
			    "HLT_Photon160_VBF_v1",
			    "HLT_Photon150_VBF_v1",
			    "HLT_Photon135_VBF_v1",
			    "HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_VBF_v1",
			    "HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF_v1",
			    "HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF_v1",
			    "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF_v1",
			    "HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_VBF_v1"
};


  Int_t HLT_prescale[nT]={1,1,1,3,6,17,700,4000,1,1,1,2,3,6,700,4000};

  Double_t thr[nT]={160.,150.,135.,90.,75.,50.,36.,22.,160.,150.,135.,90.,75.,50.,36.,22.};

  setTDRStyle();
  
  // MET turn-on [PFMET40 ]
  TH1D *h_met_pho160=readHist("sel_pho160_met",ifile,0);
  TH1D *h_met_pho160_trig=readHist("sel_pho160_trig135_met",ifile,0);

   // Plot eff turn-on
  TCanvas *cm=getanotherCanvas("cm");
  //c0->Divide(2,1);
  //c0->cd(1);
  gPad->SetGridx(); gPad->SetGridy();

  Int_t ebin=0;
  
  TLegend *leg=legend();

  TGraphAsymmErrors *ieff=drawEff(0,"sel_photon_VBF_met","sel_photon_VBF_PFMET_met",ifile,"Photon trigger PFMET40 (VBF x pfMET)/VBF","pfMET (GeV)",ebin,1,20.,450.,0.1,1000.01,kBlue-1,20,"",40.);
  ieff->GetXaxis()->SetRangeUser(0.,1000.);
  TGraphAsymmErrors *ieff2=drawEff(0,"sel_pho160_met","sel_pho160_trig135_met",ifile,"Photon trigger Photon135_PFMET40","pfMET (GeV)",ebin,1,20.,500.,0.1,1000.01,kRed-2,20,"SAME",40.);

  leg->AddEntry(ieff,"Ref: VBF trigger; Sel: VBF OR PFMET40","LP");
  leg->AddEntry(ieff2,"Ref: Photon135; Sel: Photon135_PFMET40","LP");
  
  TGraphAsymmErrors *ieff3=drawEff(0,"sel_met","sel_met40_met",ifile,"PFMET40","pfMET (GeV)",ebin,1,30.,1000.,0.1,1000.01,kRed,20,"SAME",40.);
  leg->Draw("SAME");
  
  cm->Update();
  
  //return;
  // Photon triggre vs photon pt
  
  Int_t rbin=8;

  Int_t iS=0;
  // if (tag=="sel_Njet0") { iS=1; }

  // Photon after pre-selection (At least 1 good photon)
  TH1D *h_sel_phopt=readHist(tag+"_"+"phopt",ifile,rbin);
  TH1D *h_sel_phoeta=readHist(tag+"_"+"phoeta",ifile,0);

  TH1D *h_sel_allTrg_phopt=readHist(tag+"_allTrg_"+"phopt",ifile,rbin);
  TH1D *h_sel_allTrg_phoeta=readHist(tag+"_allTrg_"+"phoeta",ifile,0);

  // Photon after Trigger selection
  Bool_t exist[nT];
  
  TH1D *h_sel_phopt_trg[nT];
  for (Int_t i=iS; i<nT; i++) {

    exist[i]=true;

    h_sel_phopt_trg[i]=readHist(tag+"_"+HLT_trigName[i]+"_phopt",ifile,rbin);
    if (h_sel_phopt_trg[i]==0) {
      std::cout << "Historgam does not exist!!" << std::endl;
      exist[i]=false;

    }
  }

  // pfMET
  // TH1D *h_sel_met=readHist(tag+"_"+"met",ifile,rbin);
  // TH1D *h_sel_allTrg_met=readHist(tag+"_allTrg_"+"met",ifile,rbin);

  //  return;
  
  // If use prescaling  
  TH1D *h_sel_phopt_pretrg[nT];

  gStyle->SetOptStat(1);
  
  Int_t color[nT]={1,2,30,4,7,9,12,40,50,3,25,35,46,52,58,67}; //kRed,kMagenta,kViolet,kBlue,kCyan,kGreen,kSpring,kOrange,kBlack};

  //  if (tag=="sel") {
  TCanvas *c=getanotherCanvas("c");
  c->Divide(2,1);
  c->cd(1); // photon pt
  gPad->SetLogy(); 
  gPad->SetLogx();
  
  TLegend *leg=legend();
  leg->AddEntry(h_sel_phopt,"Offline p_{T}^{#gamma}","L");
  
  h_sel_phopt->Draw("EHIST");
  //h_sel_phopt->Fit("pol2","FR+","EHIST",60.,600.);
  // TF1 *f1=h_sel_phopt->GetFunction("pol2");
  //f1->Draw("SAME");
  h_sel_phopt->GetYaxis()->SetRangeUser(0.5,1000000.);
  h_sel_phopt->GetXaxis()->SetRangeUser(30.01,1000.01);
 
  leg->AddEntry(h_sel_allTrg_phopt,"Offline p_{T}^{#gamma} - after Trigger","L");
  
  h_sel_allTrg_phopt->Draw("EHISTSAME");
  h_sel_allTrg_phopt->SetLineColor(kBlue-1);
  
  h_sel_allTrg_phopt->GetYaxis()->SetRangeUser(0.5,1000000.);
  h_sel_allTrg_phopt->GetXaxis()->SetRangeUser(30.01,1000.01);
  
  leg->Draw("SAME");

  if (tag=="sel") {
    c->cd(2); //photon eta
  
    TLegend *leg=legend();
    leg->AddEntry(h_sel_phoeta,"Offline #eta^{#gamma}","L");
    
    h_sel_phoeta->Draw("EHIST");
    // h_sel_phoeta->GetYaxis()->SetRangeUser(0.5,100000.);
    h_sel_phoeta->GetXaxis()->SetRangeUser(0.,3.0);
    
    leg->AddEntry(h_sel_allTrg_phoeta,"Offline #eta^{#gamma} - after Trigger","L");
    
    h_sel_allTrg_phoeta->Draw("EHISTSAME");
    h_sel_allTrg_phoeta->SetLineColor(kBlue-1);
    
    leg->Draw("SAME");
  }
  
  c->Update();
  // c->SaveAs("canvas4.C");

  // return;

  // Plot eff turn-on
  TCanvas *c0=getanotherCanvas("c0");
  c0->Divide(2,1);
  c0->cd(1);
  gPad->SetGridx(); gPad->SetGridy();
  
  TGraphAsymmErrors *ieff=drawEff(0,tag+"_"+"phopt",tag+"_allTrg_"+"phopt",ifile,"All photon triggers","offline photon p_{T} (GeV)",5,0,30.,1000.,0.1,1000.01,kBlue-1,20,"",55.);

  if (tag=="sel") {
    c0->cd(2);
    gPad->SetGridx(); gPad->SetGridy();
    
    TGraphAsymmErrors *etaEff = new TGraphAsymmErrors();
    etaEff->BayesDivide(h_sel_allTrg_phoeta, h_sel_phoeta);
    etaEff->Draw("AP");
    etaEff->GetYaxis()->SetRangeUser(0.,1.1);
    etaEff->GetXaxis()->SetRangeUser(0.,3.);
  }
  c0->Update();

  //return;

  

  // Draw individual Photon Triggers
  TGraphAsymmErrors *eff[nT];

  TCanvas *c_trg=getanotherCanvas("ctrg");
  gPad->SetLogy(); 
  //  gPad->SetLogx();

  TLegend *leg=legend();

  for (Int_t i=iS; i<nT; i++) {

    //h_sel_phopt_trg[i]->Sumw2();
    if (!exist[i]) continue;
    
    if (prescale) {

      h_sel_phopt_pretrg[i] = hTemp(h_sel_phopt_trg[i],HLT_prescale[i]);
      h_sel_phopt_pretrg[i]->Draw("EHISTSAME");
      h_sel_phopt_pretrg[i]->SetLineColor(color[i]);
      
      leg->AddEntry(h_sel_phopt_pretrg[i],HLT_trigName[i],"L");

      if (iS==0) {
	h_sel_phopt_pretrg[0]->GetYaxis()->SetRangeUser(0.5,100000.);
	h_sel_phopt_pretrg[0]->GetXaxis()->SetRangeUser(30.01,1000.01);
      } else {
	h_sel_phopt_pretrg[4]->GetYaxis()->SetRangeUser(0.5,100000.);
	h_sel_phopt_pretrg[4]->GetXaxis()->SetRangeUser(30.01,1000.01);
      }
    } else {
      h_sel_phopt_trg[i]->Draw("EHISTSAME");
      h_sel_phopt_trg[i]->SetLineColor(color[i]);

      leg->AddEntry(h_sel_phopt_trg[i],HLT_trigName[i],"L");
      if (iS==0) {
	h_sel_phopt_trg[0]->GetYaxis()->SetRangeUser(0.5,100000.);
	h_sel_phopt_trg[0]->GetXaxis()->SetRangeUser(30.01,1000.01);
      } else {
	h_sel_phopt_trg[4]->GetYaxis()->SetRangeUser(0.5,100000.);
	h_sel_phopt_trg[4]->GetXaxis()->SetRangeUser(30.01,1000.01);
      }     
     
    }// no prescale
  
    leg->Draw("SAME");
  }
  c_trg->Update();

  // return;
  
  Double_t xmax=600.01;

  // Draw individual eff turn-ons
  TCanvas *ct=getanotherCanvas("ct");
  ct->Divide(2,1);
  
  ct->cd(2);
  gPad->SetGridx(); gPad->SetGridy();
  // gPad->SetLogx();
  
  TLegend *leg=legend();
  //rbin=4;
  for (Int_t i=0; i<3; i++) {

    if (!exist[i]) continue;
    
    if (i==0) {
      eff[i]=drawEff(1,tag+"_"+"phopt",tag+"_eff_"+HLT_trigName[i]+"_phopt",ifile,"photon triggers","offline photon p_{T} (GeV)",rbin,0,80.,xmax,0.1,xmax,color[i],20,"",thr[i]);
    } else {
      eff[i]=drawEff(1,tag+"_"+"phopt",tag+"_eff_"+HLT_trigName[i]+"_phopt",ifile,"photon triggers","offline photon p_{T} (GeV)",rbin,0,80.,xmax,0.1,xmax,color[i],20,"SAME",thr[i]);
    }
    leg->AddEntry(eff[i],HLT_trigName[i],"LP");
  }
  leg->Draw("SAME");

  ct->cd(1);
  //TCanvas *ct2=getanotherCanvas("ct2");
  // gPad->SetLogx();
  
  TLegend *leg=legend();
  // rbin=5;
  gPad->SetGridx(); gPad->SetGridy();
  for (Int_t i=3; i<8; i++) {

    if (!exist[i]) continue;
    
    if (i==3) {
      eff[i]=drawEff(0,tag+"_"+"phopt",tag+"_eff_"+HLT_trigName[i]+"_phopt",ifile,"photon triggers","offline photon p_{T} (GeV)",rbin,0,30.,xmax,0.1,xmax,color[i],20,"",thr[i]);
    } else {
      eff[i]=drawEff(0,tag+"_"+"phopt",tag+"_eff_"+HLT_trigName[i]+"_phopt",ifile,"photon triggers","offline photon p_{T} (GeV)",rbin,0,30.,xmax,0.1,xmax,color[i],20,"SAME",thr[i]);
    }
    leg->AddEntry(eff[i],HLT_trigName[i],"LP");
  }
  leg->Draw("SAME");
  
  ct->Update();


  // VBF
   TCanvas *ct2=getanotherCanvas("ct2");
  ct2->Divide(2,1);
  
  ct2->cd(2);
  gPad->SetGridx(); gPad->SetGridy();
  //gPad->SetLogx();
 
  TLegend *leg=legend();
  // rbin=10;
  
  for (Int_t i=8; i<11; i++) {

    if (!exist[i]) continue;
    
    if (i==8) {
      eff[i]=drawEff(1,tag+"_"+"phopt",tag+"_eff_"+HLT_trigName[i]+"_phopt",ifile,"photon triggers","offline photon p_{T} (GeV)",rbin,0,80.,xmax,0.1,xmax,color[i],20,"",thr[i]);
    } else {
      eff[i]=drawEff(1,tag+"_"+"phopt",tag+"_eff_"+HLT_trigName[i]+"_phopt",ifile,"photon triggers","offline photon p_{T} (GeV)",rbin,0,80.,xmax,0.1,xmax,color[i],20,"SAME",thr[i]);
    }
    leg->AddEntry(eff[i],HLT_trigName[i],"LP");
  }
  leg->Draw("SAME");

  ct2->cd(1);
  //TCanvas *ct2=getanotherCanvas("ct2");
  
  TLegend *leg=legend();
  //rbin=5;
  gPad->SetGridx(); gPad->SetGridy();
  //gPad->SetLogx();
   
  for (Int_t i=11; i<nT; i++) {

    if (!exist[i]) continue;
    
    if (i==11) {
      eff[i]=drawEff(0,tag+"_"+"phopt",tag+"_eff_"+HLT_trigName[i]+"_phopt",ifile,"photon triggers","offline photon p_{T} (GeV)",rbin,0,30.,xmax,0.1,xmax,color[i],20,"",thr[i]);
    } else {
      eff[i]=drawEff(0,tag+"_"+"phopt",tag+"_eff_"+HLT_trigName[i]+"_phopt",ifile,"photon triggers","offline photon p_{T} (GeV)",rbin,0,30.,xmax,0.1,xmax,color[i],20,"SAME",thr[i]);
    }
    leg->AddEntry(eff[i],HLT_trigName[i],"LP");
  }
  leg->Draw("SAME");
  
  ct2->Update();

  }
