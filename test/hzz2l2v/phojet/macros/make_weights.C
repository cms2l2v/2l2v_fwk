#include "setTDRStyle.C"


Double_t fermipow(Double_t *x, Double_t *par) {

  const Int_t n=6;
  Double_t thr[n]={22,36,50,75,90,300};

  Double_t fland;
  Double_t fpow, ffermi;

  Double_t alpha,beta,gamma,delta;
  
  Int_t j=1;
  
  for (Int_t i=0; i<n; i++) {

    alpha=par[j]; beta=par[j+1];
    gamma=par[j+2]; delta=par[j+3];
    
    fpow = TMath::Power( (x[0]/alpha), beta);
    ffermi = (TMath::Exp(-x[0]/alpha) / (1.+TMath::Exp( (gamma-x[0])/delta )));

    //fland += TMath::Floor(x[0]) * fpow * ffermi;
    if (x[0]>=thr[i]) {
      fland += fpow * ffermi;
    } else { fland += 0.;}
    
    j+=4;
  }
  //return fland;
  return par[0]*fland;
}

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


TH1F* readHist(TString nameHist,TString nameFile, TString nameDir, int rebin)
{
 TFile* file = new TFile(nameFile);
 //file->ls();
 TDirectory *dir=(TDirectory*)file->Get(nameDir);
 
 TH1F* hist = (TH1F*)dir->Get(nameHist);
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
 TFile* file = new TFile(nameFile);

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

  TCanvas* aCanvas = new TCanvas(name,"",292,55,500,405);//,"",181,237,1575,492);

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

  TCanvas* aCanvas = new TCanvas(name,"",292,55,1100,405);//,"",181,237,1575,492);

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
 leg2->SetTextSize(0.04);
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

TGraphAsymmErrors* drawEff(TString hname1, TString hname2, TString filename, TString header, TString xtitle, int rebin, bool asymBin, float xlow, float xhigh, float xmin, float xmax, int icol, int imark, TString draw, double mup){//,double y, double c) {

  TH1D *ref = readHist(hname1, filename, rebin);
  TH1D *sel = readHist(hname2, filename, rebin);

  ref->SetName(hname1+header);
  sel->SetName(hname2+header);

  TGraphAsymmErrors *Eff = new TGraphAsymmErrors();
  if (asymBin) {
    Double_t xbins[14]={0.,4.,6.,8.,10.,12.,14.,18.
		      22.,28.,34.,42.,52.,100.};

    sel->Rebin(14,"nsel",xbins);
    ref->Rebin(14,"nref",xbins);
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

  Double_t params[5] = {1.,mup,20.,1.,xlow};
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
  
  for(int i=0; i != 1 ; i++){  
    Eff->Fit("fermiFunction","RS+");
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

TH1D *hNorm(TH1D *h, double xsec, double nevtmc, double lum ) {
  // xsec is in pb , lum is in fb-1
  // h->Sumw2();

  double weight;
  weight=((xsec * lum * 1000.)/nevtmc);
  
  h->Scale(weight);
  h->Sumw2();
  
  return h;
  
}

TGraphErrors *makeGraph(TH1* htest) {

  TH1D *h=htest->Clone();
  
  Double_t err;
  
  Int_t nB=h->GetNbinsX();
  const Int_t n=nB;
  
  Double_t x[n],y[n];
  Double_t xe[n],ye[n];

  TGraphErrors *gr;
  
  for (Int_t i=0; i<nB; i++) {

    x[i]=h->GetBinCenter(i); xe[i]=h->GetBinWidth(i);
    cout << "bin center = " << x[i] << " and bin widht = " << xe[i] << endl;
    y[i]=h->GetBinContent(i); ye[i]=h->GetBinError(i);
    cout << " bin content " << y[i] << " and bin content error = " << ye[i] << endl; 
    
  }

  gr=new TGraphErrors(n,x,y,xe,ye);

  return gr;

}


void make_weights(TString tag="sel", TString var="_qt",  Bool_t asym=false) {

  Bool_t isData=true;
  
  TString ifile="plotter_2015_12_04.root";
  
  TString dydir, dydir2; // DY+jets
  if (isData) { dydir="emudata"; }
  else {
    dydir="Ztoll_M-50";
    dydir2="Ztoll_M-10to50";
  }
  
  TString gdir, gdir2, gdir3, gdir4, gdir5; // G+jets
  if (isData) {gdir="photondata";}
  else {
    gdir="GJets_HT-40to100";
    gdir2="GJets_HT-100to200";
    gdir3="GJets_HT-200to400";
    gdir4="GJets_HT-400to600";
    gdir5="GJets_HT-600toInf";
  }
  
  Bool_t bfit=false;
  
  double xmin,xmax;
  double ymin, ymax;
  ymin=0.5; ymax=100000000.;
  
  if (var=="_qt") {
    // bfit=true;
    xmin=45.01; xmax=1000.01;
  } else if ( (var=="_met") || var=="_mt") {
    xmin=1.; xmax=1000.0;
  } 
 
  setTDRStyle();
  
  // TString cat1,cat2,hcat[2];
  // cat1="eq0jets"; hcat[1]="geq1jets"; hcat[2]="vbf";

  TString hcat[3]={"eq0jets","geq1jets","vbf"};
  
  TString cat="ee";
  
  Int_t rbin=5;
  
  TH1F *h_zpt_ee_1=readHist(cat+hcat[0]+var,ifile,dydir,rbin); //h_zpt_ee->Scale(dyweight);
  TH1F *h_zpt_ee_2=readHist(cat+hcat[1]+var,ifile,dydir,rbin);
  TH1F *h_zpt_ee_3=readHist(cat+hcat[2]+var,ifile,dydir,rbin);
  h_zpt_ee_1->SetFillColor(0);
  h_zpt_ee_2->SetFillColor(0);
  h_zpt_ee_3->SetFillColor(0);

  // Adddy2
  // TH1F *h2_zpt_ee_1=readHist(cat+hcat[0]+var,ifile,dydir2,rbin); //h_zpt_ee->Scale(dyweight);
  // TH1F *h2_zpt_ee_2=readHist(cat+hcat[1]+var,ifile,dydir2,rbin);
  // TH1F *h2_zpt_ee_3=readHist(cat+hcat[2]+var,ifile,dydir2,rbin);
  // h_zpt_ee_1->Add(h2_zpt_ee_1);
  // h_zpt_ee_2->Add(h2_zpt_ee_2);
  // h_zpt_ee_3->Add(h2_zpt_ee_3);

  
    
  TCanvas *c=getanotherCanvas("c");
  c->Divide(3,1);

  c->cd(1);
  gPad->SetLogy(); gPad->SetLogx();
  TLegend *leg=legend();
  leg->SetHeader(cat+" : "+hcat[0]);

  h_zpt_ee_1->Draw("EHIST");
  if (bfit) h_zpt_ee_1->Fit("pol2","FR+","EHISTSAME",60.,500.);
  h_zpt_ee_1->GetXaxis()->SetRangeUser(xmin,xmax);
  h_zpt_ee_1->GetYaxis()->SetRangeUser(ymin,ymax);
  
  leg->Draw("SAME");
  
  c->cd(2); 
  gPad->SetLogy(); gPad->SetLogx();
  TLegend *leg=legend();
  leg->SetHeader(cat+" : "+hcat[1]);
  
  h_zpt_ee_2->Draw("EHIST");
  if (bfit) h_zpt_ee_2->Fit("pol2","FR+","EHIST",60.,500.);
  h_zpt_ee_2->GetXaxis()->SetRangeUser(xmin,xmax);
  h_zpt_ee_2->GetYaxis()->SetRangeUser(ymin,ymax);
  leg->Draw("SAME");

  c->cd(3); 
  gPad->SetLogy(); gPad->SetLogx();
  TLegend *leg=legend();
  leg->SetHeader(cat+" : "+hcat[2]);
  
  h_zpt_ee_3->Draw("EHIST");
  if (bfit) h_zpt_ee_3->Fit("pol2","FR+","",60.,500.);
  h_zpt_ee_3->GetXaxis()->SetRangeUser(xmin,xmax);
  h_zpt_ee_3->GetYaxis()->SetRangeUser(ymin,ymax);
  leg->Draw("SAME");

  // DIMUONS
  cat="mumu";

  std::cout << "=0jets: " << cat+hcat[0]+var << std::endl;
  std::cout << ">=2: " << cat+hcat[1]+var << std::endl;
  std::cout << "VBF: " << cat+hcat[2]+var << std::endl;
  
  TH1F *h_zpt_mm_1=readHist(cat+hcat[0]+var,ifile,dydir,rbin);//h_zpt_mm->Scale(dyweight);
  TH1F *h_zpt_mm_2=readHist(cat+hcat[1]+var,ifile,dydir,rbin);//h_zpt_mm_2->Scale(dyweight);
  TH1F *h_zpt_mm_3=readHist(cat+hcat[2]+var,ifile,dydir,rbin);//h_zpt_mm_3->Scale(dyweight);
  h_zpt_mm_1->SetFillColor(0);
  h_zpt_mm_2->SetFillColor(0);
  h_zpt_mm_3->SetFillColor(0);

 // Adddy2
  // TH1F *h2_zpt_mm_1=readHist(cat+hcat[0]+var,ifile,dydir2,rbin); //h_zpt_mm->Scale(dyweight);
  // TH1F *h2_zpt_mm_2=readHist(cat+hcat[1]+var,ifile,dydir2,rbin);
  // TH1F *h2_zpt_mm_3=readHist(cat+hcat[2]+var,ifile,dydir2,rbin);
  // h_zpt_mm_1->Add(h2_zpt_mm_1);
  // h_zpt_mm_2->Add(h2_zpt_mm_2);
  // h_zpt_mm_3->Add(h2_zpt_mm_3);

  
  TCanvas *cm=getanotherCanvas("cm");
  cm->Divide(3,1);

  cm->cd(1);
  gPad->SetLogy(); gPad->SetLogx();
   TLegend *leg=legend();
  leg->SetHeader(cat+" : "+hcat[0]);
  
  h_zpt_mm_1->Draw("EHIST");
  h_zpt_mm_1->GetXaxis()->SetRangeUser(xmin,xmax);
  h_zpt_mm_1->GetYaxis()->SetRangeUser(ymin,ymax);
  leg->Draw("SAME");
  
  cm->cd(2); 
  gPad->SetLogy(); gPad->SetLogx();
  TLegend *leg=legend();
  leg->SetHeader(cat+" : "+hcat[1]);
  
  h_zpt_mm_2->Draw("EHIST");
  h_zpt_mm_2->GetXaxis()->SetRangeUser(xmin,xmax);
  h_zpt_mm_2->GetYaxis()->SetRangeUser(ymin,ymax);
  leg->Draw("SAME");

  cm->cd(3); 
  gPad->SetLogy(); gPad->SetLogx();
  TLegend *leg=legend();
  leg->SetHeader(cat+" : "+hcat[2]);
  
  h_zpt_mm_3->Draw("EHIST");
  h_zpt_mm_3->GetXaxis()->SetRangeUser(xmin,xmax);
  h_zpt_mm_3->GetYaxis()->SetRangeUser(ymin,ymax);
  leg->Draw("SAME");

  // Photon+jet weight
  //double gweight= ((224012.*20.*1000.)/12365798.);

  cat="ee";
  //  var="_phopt";
  
  // Photon + Jet
  std::cout << "=0jets: " << cat+hcat[0]+var << std::endl;
  std::cout << ">=2: " << cat+hcat[1]+var << std::endl;
  std::cout << "VBF: " << cat+hcat[2]+var << std::endl;
  
  TH1F *h_zpt_gj_1=readHist(cat+hcat[0]+var,ifile,gdir,rbin);// h_zpt_gj->Scale(gweight);//h_zpt_gj->Sumw2();
  TH1F *h_zpt_gj_2=readHist(cat+hcat[1]+var,ifile,gdir,rbin); //h_zpt_gj_2->Scale(gweight);//h_zpt_gj_2->Sumw2();
  TH1F *h_zpt_gj_3=readHist(cat+hcat[2]+var,ifile,gdir,rbin);// h_zpt_gj_3->Scale(gweight);//h_zpt_gj_3->Sumw2();
  h_zpt_gj_1->SetFillColor(0);
  h_zpt_gj_2->SetFillColor(0);
  h_zpt_gj_3->SetFillColor(0);

  if (!isData) {
    
  TH1F *h2_zpt_gj_1=readHist(cat+hcat[0]+var,ifile,gdir2,rbin);
  TH1F *h2_zpt_gj_2=readHist(cat+hcat[1]+var,ifile,gdir2,rbin); 
  TH1F *h2_zpt_gj_3=readHist(cat+hcat[2]+var,ifile,gdir2,rbin);
  h2_zpt_gj_1->SetFillColor(0);
  h2_zpt_gj_2->SetFillColor(0);
  h2_zpt_gj_3->SetFillColor(0);
 
  TH1F *h3_zpt_gj_1=readHist(cat+hcat[0]+var,ifile,gdir3,rbin);
  TH1F *h3_zpt_gj_2=readHist(cat+hcat[1]+var,ifile,gdir3,rbin); 
  TH1F *h3_zpt_gj_3=readHist(cat+hcat[2]+var,ifile,gdir3,rbin);
  h3_zpt_gj_1->SetFillColor(0);
  h3_zpt_gj_2->SetFillColor(0);
  h3_zpt_gj_3->SetFillColor(0);

  TH1F *h4_zpt_gj_1=readHist(cat+hcat[0]+var,ifile,gdir4,rbin);
  TH1F *h4_zpt_gj_2=readHist(cat+hcat[1]+var,ifile,gdir4,rbin); 
  TH1F *h4_zpt_gj_3=readHist(cat+hcat[2]+var,ifile,gdir4,rbin);
  h4_zpt_gj_1->SetFillColor(0);
  h4_zpt_gj_2->SetFillColor(0);
  h4_zpt_gj_3->SetFillColor(0);

  TH1F *h5_zpt_gj_1=readHist(cat+hcat[0]+var,ifile,gdir5,rbin);
  TH1F *h5_zpt_gj_2=readHist(cat+hcat[1]+var,ifile,gdir5,rbin); 
  TH1F *h5_zpt_gj_3=readHist(cat+hcat[2]+var,ifile,gdir5,rbin);
  h5_zpt_gj_1->SetFillColor(0);
  h5_zpt_gj_2->SetFillColor(0);
  h5_zpt_gj_3->SetFillColor(0);
  
  h_zpt_gj_1->Add(h2_zpt_gj_1); h_zpt_gj_1->Add(h3_zpt_gj_1); h_zpt_gj_1->Add(h4_zpt_gj_1); h_zpt_gj_1->Add(h5_zpt_gj_1);
  h_zpt_gj_2->Add(h2_zpt_gj_2); h_zpt_gj_2->Add(h3_zpt_gj_2); h_zpt_gj_2->Add(h4_zpt_gj_2); h_zpt_gj_2->Add(h5_zpt_gj_2);
  h_zpt_gj_3->Add(h2_zpt_gj_3); h_zpt_gj_3->Add(h3_zpt_gj_3); h_zpt_gj_3->Add(h4_zpt_gj_3); h_zpt_gj_3->Add(h5_zpt_gj_3);

  } // if MC
  
  cat="g+jet";

  TCanvas *cg=getanotherCanvas("cg");
  cg->Divide(3,1);
   
  cg->cd(1);
  gPad->SetLogy(); gPad->SetLogx();
  TLegend *leg=legend();
  leg->SetHeader(cat+" : "+hcat[0]);

  // fermiFunction->SetParameters(params);
  
  h_zpt_gj_1->Draw("EHIST");
 
  h_zpt_gj_1->GetXaxis()->SetRangeUser(xmin,xmax);
  //h_zpt_gj->GetYaxis()->SetRangeUser(ymin,ymax);
  leg->Draw("SAME");
  
  cg->cd(2); 
  gPad->SetLogy(); gPad->SetLogx();
  TLegend *leg=legend();
  leg->SetHeader(cat+" : "+hcat[1]);

  //fermiFunction->SetParameters(params);
  
  h_zpt_gj_2->Draw("EHIST");
 
  h_zpt_gj_2->GetXaxis()->SetRangeUser(xmin,xmax);
  //h_zpt_gj_2->GetYaxis()->SetRangeUser(ymin,ymax);
  leg->Draw("SAME");

  cg->cd(3); 
  gPad->SetLogy(); gPad->SetLogx();
  TLegend *leg=legend();
  leg->SetHeader(cat+" : "+hcat[2]);

  // fermiFunction->SetParameters(params);
  
  h_zpt_gj_3->Draw("EHIST");

  h_zpt_gj_3->GetXaxis()->SetRangeUser(xmin,xmax);
  //h_zpt_gj_3->GetYaxis()->SetRangeUser(ymin,ymax);
  leg->Draw("SAME");

  // Superimpose all
  rbin=1;
  
 Double_t xbins[27] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100,
  			125, 150, 175, 200, 250, 300, 400, 500, 1000,1500};
  
  // hee1->Rebin(26,"nhee1",xbins);  hee2->Rebin(26,"nhee2",xbins);  hee3->Rebin(26,"nhee3",xbins); 
  // hmm1->Rebin(26,"nhmm1",xbins);  hmm2->Rebin(26,"nhmm2",xbins);  hmm3->Rebin(26,"nhmm3",xbins); 

  
  TCanvas *call=getanotherCanvas("call");
  call->Divide(3,1);
  
  // TH1F *h_ee1=h_zpt_ee_1->Clone(); h_ee1->Rebin(rbin);
  // TH1F *h_mm1=h_zpt_mm_1->Clone(); h_mm1->Rebin(rbin);
  // TH1F *h_gj1=h_zpt_gj_1->Clone(); h_gj1->Rebin(rbin);
  if (asym) {
    h_zpt_ee_1->Rebin(26,"h_ee1",xbins);
    h_zpt_mm_1->Rebin(26,"h_mm1",xbins);
    h_zpt_gj_1->Rebin(26,"h_gj1",xbins);
  } else {
      TH1F *h_ee1=h_zpt_ee_1->Clone(); h_ee1->Rebin(rbin);
      TH1F *h_mm1=h_zpt_mm_1->Clone(); h_mm1->Rebin(rbin);
      TH1F *h_gj1=h_zpt_gj_1->Clone(); h_gj1->Rebin(rbin);
    }
  // h_ee1->Scale(1./h_ee1->Integral());
  // h_mm1->Scale(1./h_mm1->Integral());
  // h_gj1->Scale(1./h_gj1->Integral());

  if (asym) {
    h_zpt_ee_2->Rebin(26,"h_ee2",xbins);
    h_zpt_mm_2->Rebin(26,"h_mm2",xbins);
    h_zpt_gj_2->Rebin(26,"h_gj2",xbins);
  } else {
    TH1F *h_ee2=h_zpt_ee_2->Clone(); h_ee2->Rebin(rbin);
    TH1F *h_mm2=h_zpt_mm_2->Clone(); h_mm2->Rebin(rbin);
    TH1F *h_gj2=h_zpt_gj_2->Clone(); h_gj2->Rebin(rbin);
  }
  // h_ee2->Scale(1./h_ee2->Integral());
  // h_mm2->Scale(1./h_mm2->Integral());
  // h_gj2->Scale(1./h_gj2->Integral());

 if (asym) {
    h_zpt_ee_3->Rebin(26,"h_ee3",xbins);
    h_zpt_mm_3->Rebin(26,"h_mm3",xbins);
    h_zpt_gj_3->Rebin(26,"h_gj3",xbins);
  } else {
   TH1F *h_ee3=h_zpt_ee_3->Clone(); h_ee3->Rebin(rbin);
   TH1F *h_mm3=h_zpt_mm_3->Clone(); h_mm3->Rebin(rbin);
   TH1F *h_gj3=h_zpt_gj_3->Clone(); h_gj3->Rebin(rbin);
 }
  // h_ee3->Scale(1./h_ee3->Integral());
  // h_mm3->Scale(1./h_mm3->Integral());
  // h_gj3->Scale(1./h_gj3->Integral());

  call->cd(1);
  gPad->SetLogy();

  TLegend *leg=legend();
  leg->SetHeader("=0 jets");
  leg->AddEntry(h_ee1,"ee","LP");
  leg->AddEntry(h_mm1,"#mu#mu","LP");
  leg->AddEntry(h_gj1,"#gamma+jet","LP");
  
  h_ee1->Draw("ehist");
  h_ee1->GetXaxis()->SetRangeUser(xmin,xmax);
  h_ee1->SetLineColor(kBlue);
  h_ee1->SetMarkerColor(kBlue); h_ee1->SetLineWidth(2); h_ee1->SetMarkerStyle(20); h_ee1->SetMarkerSize(0.3);
  h_mm1->Draw("ehistsame");
  h_mm1->SetLineColor(kBlue-6);  
  h_mm1->SetMarkerColor(kBlue-6); h_mm1->SetLineWidth(2); h_mm1->SetMarkerStyle(22); h_mm1->SetMarkerSize(0.3);
  h_gj1->Draw("ehistsame");
  h_gj1->SetLineColor(kRed+1);
  h_gj1->SetMarkerColor(kRed+1); h_gj1->SetLineWidth(2); h_gj1->SetMarkerStyle(29); h_gj1->SetMarkerSize(0.3);

  leg->Draw("SAME");

  call->cd(2);
  gPad->SetLogy();

  TLegend *leg=legend();
  leg->SetHeader(">= 1jets");
  leg->AddEntry(h_ee2,"ee","LP");
  leg->AddEntry(h_mm2,"#mu#mu","LP");
  leg->AddEntry(h_gj2,"#gamma+jet","LP");
  
  h_ee2->Draw("ehist");
  h_ee2->GetXaxis()->SetRangeUser(xmin,xmax);
  h_ee2->SetLineColor(kBlue);
  h_ee2->SetMarkerColor(kBlue); h_ee2->SetLineWidth(2); h_ee2->SetMarkerStyle(20); h_ee2->SetMarkerSize(0.3);
  h_mm2->Draw("ehistsame");
  h_mm2->SetLineColor(kBlue-6);  
  h_mm2->SetMarkerColor(kBlue-6); h_mm2->SetLineWidth(2); h_mm2->SetMarkerStyle(22); h_mm2->SetMarkerSize(0.3);
  h_gj2->Draw("ehistsame");
  h_gj2->SetLineColor(kRed+1);
  h_gj2->SetMarkerColor(kRed+1); h_gj2->SetLineWidth(2); h_gj2->SetMarkerStyle(29); h_gj2->SetMarkerSize(0.3);

  leg->Draw("SAME");

  call->cd(3);
  gPad->SetLogy();

  TLegend *leg=legend();
  leg->SetHeader("VBF");
  leg->AddEntry(h_ee3,"ee","LP");
  leg->AddEntry(h_mm3,"#mu#mu","LP");
  leg->AddEntry(h_gj3,"#gamma+jet","LP");
  
  h_ee3->Draw("ehist");
  h_ee3->GetXaxis()->SetRangeUser(xmin,10000000.);
  h_ee3->SetLineColor(kBlue);
  h_ee3->SetMarkerColor(kBlue); h_ee3->SetLineWidth(2); h_ee3->SetMarkerStyle(20); h_ee3->SetMarkerSize(0.3);
  h_mm3->Draw("ehistsame");
  h_mm3->SetLineColor(kBlue-6);  
  h_mm3->SetMarkerColor(kBlue-6); h_mm3->SetLineWidth(2); h_mm3->SetMarkerStyle(22); h_mm3->SetMarkerSize(0.3);
  h_gj3->Draw("ehistsame");
  h_gj3->SetLineColor(kRed+1);
  h_gj3->SetMarkerColor(kRed+1); h_gj3->SetLineWidth(2); h_gj3->SetMarkerStyle(29); h_gj3->SetMarkerSize(0.3);

  leg->Draw("SAME");
 
  // return;
 
  // REDRAW in one CANVAS
  TCanvas *call1 = new TCanvas("call1");
  call1->SetWindowSize(500,600);
  call1->cd();
  
  //distributions
  TPad *t1 = new TPad("p1","p1",0,0.3,1.0,1.0);
  t1->Draw();
  t1->cd();
  t1->SetLogy(); t1->SetLogx();
  t1->SetTopMargin(0.08);
  t1->SetBottomMargin(0);
  t1->SetRightMargin(0.05);

  h_ee1->SetTitle("ee");
  h_ee1->GetYaxis()->SetLabelSize(0.06);
  h_ee1->GetYaxis()->SetTitleSize(0.05);      
  h_ee1->GetYaxis()->SetTitleOffset(0.9);
  h_ee1->GetYaxis()->SetLabelSize(0.05);
  h_ee1->GetXaxis()->SetRangeUser(25.01,1000.01); 
  h_ee1->GetYaxis()->SetTitle("Events / 100 pb^{-1}");
  // h_ee1->SetLineColor(1);
  // h_ee1->SetMarkerColor(1);
  // h_ee1->SetMarkerStyle(1);
  // h_ee1->SetFillStyle(3001);
  // h_ee1->SetFillColor(5);
  h_ee1->Draw("ehistsame");

  h_mm1->SetTitle("#mu#mu");
  h_mm1->Draw("ehistsame");
  // h_mm1->SetMarkerColor(9);
  // h_mm1->SetMarkerStyle(20);
  // h_mm1->SetMarkerSize(0.8);
  h_gj1->SetTitle("#gamma + jets");
  h_gj1->Draw("ehistsame");
  // h_gj1->SetFillColor(kRed+1); h_gj1->SetFillStyle(3001);
  
  TPaveText *pave = new TPaveText(0.7,0.85,0.95,0.9,"brNDC");
  pave->SetBorderSize(0);
  pave->SetFillStyle(0);
  pave->SetTextAlign(32);
  pave->SetTextFont(42);
  pave->SetTextSize(0.05);
  pave->SetTextColor(kBlue);
  pave->AddText("=0 jets");
  pave->Draw();

  TLegend *leg=new TLegend(0.6,0.6,0.9,0.8);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextAlign(12);
  leg->SetTextSize(0.05);
  leg->SetNColumns(2);
  leg->AddEntry(h_ee1,"ee","L");
  leg->AddEntry(h_mm1,"#mu#mu","L");
  leg->AddEntry(h_gj1,"#gamma+jets","L");
  leg->Draw("same");

  //closure
  call1->cd();
  TPad *t2 = new TPad("p2","p2",0,0.0,1.0,0.3);
  t2->SetTopMargin(0);
  t2->SetBottomMargin(0.25);
  t2->SetRightMargin(0.05);
  t2->Draw();
  t2->cd();
   t2->SetLogx();
  t2->SetGridx(); t2->SetGridy();
  
  leg = new TLegend(0.6,0.6,0.9,0.8,"","brNDC");
  leg->SetBorderSize(0);
  leg->SetFillStyle(3001);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.11);
  leg->SetTextAlign(12);

  TH1D *hee1=h_ee1->Clone(); hee1->Divide(h_gj1);
  TH1D *hmm1=h_mm1->Clone(); hmm1->Divide(h_gj1);
   
  hee1->Draw("e"); hee1->SetLineColor(kBlue); hee1->SetMarkerColor(kBlue);hee1->SetLineWidth(2);
  hmm1->Draw("esame"); hmm1->SetLineColor(kBlue-6); hmm1->SetMarkerColor(kBlue-6);hmm1->SetLineWidth(2);

  hee1->GetYaxis()->SetRangeUser(-0.1,0.5);
   //  hee1->GetYaxis()->SetRangeUser(-0.2,1.74);
  hee1->GetXaxis()->SetTitle(hee1->GetXaxis()->GetTitle());
  hee1->GetXaxis()->SetLabelSize(0.12);
  hee1->GetXaxis()->SetTitleSize(0.14);
  hee1->GetXaxis()->SetTitleOffset(0.8);
  hee1->GetYaxis()->SetLabelSize(0.12);
  hee1->GetYaxis()->SetNdivisions(5);
  hee1->GetYaxis()->SetTitleSize(0.12);
  //gr->GetYaxis()->SetTitle("Pred. rel. bias");
  hee1->GetYaxis()->SetTitle("Weight");
  hee1->GetYaxis()->SetTitleOffset(0.3);
  hee1->GetXaxis()->SetRangeUser(25.01,1000.01); 

  leg->AddEntry(hee1,"ee","L");
  leg->AddEntry(hmm1,"#mu#mu","L");
  leg->Draw();
  leg->SetNColumns(4);

  

  // 2
  TCanvas *call2 = new TCanvas("call2");
  call2->SetWindowSize(500,600);
  call2->cd();
  
  //distributions
  TPad *t1 = new TPad("p1","p1",0,0.3,1.0,1.0);
  t1->Draw();
  t1->cd();
  t1->SetLogy(); t1->SetLogx();
  t1->SetTopMargin(0.08);
  t1->SetBottomMargin(0);
  t1->SetRightMargin(0.05);

  h_ee2->SetTitle("ee");
  h_ee2->GetYaxis()->SetLabelSize(0.06);
  h_ee2->GetYaxis()->SetTitleSize(0.05);      
  h_ee2->GetYaxis()->SetTitleOffset(0.9);
  h_ee2->GetYaxis()->SetLabelSize(0.05);
  h_ee2->GetXaxis()->SetRangeUser(25.01,1000.01); 
  h_ee2->GetYaxis()->SetTitle("Events / 100 pb^{-1}");
 
  h_ee2->Draw("ehistsame");

  h_mm2->SetTitle("#mu#mu");
  h_mm2->Draw("ehistsame");
  // h_mm2->SetMarkerColor(9);
  // h_mm2->SetMarkerStyle(20);
  // h_mm2->SetMarkerSize(0.8);
  h_gj2->SetTitle("#gamma + jets");
  h_gj2->Draw("ehistsame");

  TPaveText *pave = new TPaveText(0.7,0.85,0.95,0.9,"brNDC");
  pave->SetBorderSize(0);
  pave->SetFillStyle(0);
  pave->SetTextAlign(32);
  pave->SetTextFont(42);
  pave->SetTextSize(0.05);
  pave->SetTextColor(kBlue);
  pave->AddText(">=1 jet");
  pave->Draw();

  TLegend *leg=new TLegend(0.6,0.6,0.9,0.8);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextAlign(12);
  leg->SetTextSize(0.05);
  leg->SetNColumns(2);
  leg->AddEntry(h_ee2,"ee","L");
  leg->AddEntry(h_mm2,"#mu#mu","L");
  leg->AddEntry(h_gj1,"#gamma+jets","L");
  leg->Draw("same");

  //closure
  call2->cd();
  TPad *t2 = new TPad("p2","p2",0,0.0,1.0,0.3);
  t2->SetTopMargin(0);
  t2->SetBottomMargin(0.25);
  t2->SetRightMargin(0.05);
  t2->Draw();
  t2->cd();
  t2->SetLogx();
  t2->SetGridx(); t2->SetGridy();
  
  leg = new TLegend(0.6,0.6,0.9,0.8,"","brNDC");
  leg->SetBorderSize(0);
  leg->SetFillStyle(3001);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.11);
  leg->SetTextAlign(12);

  TH1D *hee2=h_ee2->Clone(); hee2->Divide(h_gj2);
  TH1D *hmm2=h_mm2->Clone(); hmm2->Divide(h_gj2);
 
  hee2->Draw("e"); hee2->SetLineColor(kBlue); hee2->SetMarkerColor(kBlue);hee2->SetLineWidth(2);
  hmm2->Draw("esame"); hmm2->SetLineColor(kBlue-6); hmm2->SetMarkerColor(kBlue-6);hmm2->SetLineWidth(2);

  hee2->GetYaxis()->SetRangeUser(-0.2,0.2);
   //  hee2->GetYaxis()->SetRangeUser(-0.2,1.74);
  hee2->GetXaxis()->SetTitle(hee2->GetXaxis()->GetTitle());
  hee2->GetXaxis()->SetLabelSize(0.12);
  hee2->GetXaxis()->SetTitleSize(0.14);
  hee2->GetXaxis()->SetTitleOffset(0.8);
  hee2->GetYaxis()->SetLabelSize(0.12);
  hee2->GetYaxis()->SetNdivisions(5);
  hee2->GetYaxis()->SetTitleSize(0.12);
  //gr->GetYaxis()->SetTitle("Pred. rel. bias");
  hee2->GetYaxis()->SetTitle("Weight");
  hee2->GetYaxis()->SetTitleOffset(0.3);
  hee2->GetXaxis()->SetRangeUser(25.01,1000.01); 

  leg->AddEntry(hee2,"ee","L");
  leg->AddEntry(hmm2,"#mu#mu","L");
  leg->Draw();
  leg->SetNColumns(4);

  // 3
  TCanvas *call3 = new TCanvas("call3");
  call3->SetWindowSize(500,600);
  call3->cd();
  
  //distributions
  TPad *t1 = new TPad("p1","p1",0,0.3,1.0,1.0);
  t1->Draw();
  t1->cd();
  t1->SetLogy(); t1->SetLogx();
  t1->SetTopMargin(0.08);
  t1->SetBottomMargin(0);
  t1->SetRightMargin(0.05);

  h_ee3->SetTitle("ee");
  h_ee3->GetYaxis()->SetLabelSize(0.06);
  h_ee3->GetYaxis()->SetTitleSize(0.05);      
  h_ee3->GetYaxis()->SetTitleOffset(0.9);
  h_ee3->GetYaxis()->SetLabelSize(0.05);
  h_ee3->GetXaxis()->SetRangeUser(25.01,1000.01); 
  h_ee3->GetYaxis()->SetTitle("Events / 100 pb^{-1}");
 
  h_ee3->Draw("ehistsame");

  h_mm3->SetTitle("#mu#mu");
  h_mm3->Draw("ehistsame");
  // h_mm3->SetMarkerColor(9);
  // h_mm3->SetMarkerStyle(20);
  // h_mm3->SetMarkerSize(0.8);
  h_gj3->SetTitle("#gamma + jets");
  h_gj3->Draw("ehistsame");

  TPaveText *pave = new TPaveText(0.7,0.85,0.95,0.9,"brNDC");
  pave->SetBorderSize(0);
  pave->SetFillStyle(0);
  pave->SetTextAlign(32);
  pave->SetTextFont(42);
  pave->SetTextSize(0.05);
  pave->SetTextColor(kBlue);
  pave->AddText("VBF");
  pave->Draw();

  TLegend *leg=new TLegend(0.6,0.6,0.9,0.8);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextAlign(12);
  leg->SetTextSize(0.05);
  leg->SetNColumns(2);
  leg->AddEntry(h_ee3,"ee","L");
  leg->AddEntry(h_mm3,"#mu#mu","L");
  leg->AddEntry(h_gj1,"#gamma+jets","L");
  leg->Draw("same");

  //closure
  call3->cd();
  TPad *t2 = new TPad("p2","p2",0,0.0,1.0,0.3);
  t2->SetTopMargin(0);
  t2->SetBottomMargin(0.25);
  t2->SetRightMargin(0.05);
  t2->Draw();
  t2->cd();
  t2->SetLogx();
  t2->SetGridx(); t2->SetGridy();
  
  leg = new TLegend(0.6,0.6,0.9,0.8,"","brNDC");
  leg->SetBorderSize(0);
  leg->SetFillStyle(3001);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.11);
  leg->SetTextAlign(12);


  TH1D *hee3=h_ee3->Clone(); hee3->Divide(h_gj3);
  TH1D *hmm3=h_mm3->Clone(); hmm3->Divide(h_gj3);
   
  hee3->Draw("e"); hee3->SetLineColor(kBlue); hee3->SetMarkerColor(kBlue);hee3->SetLineWidth(2);
  hmm3->Draw("esame"); hmm3->SetLineColor(kBlue-6); hmm3->SetMarkerColor(kBlue-6);hmm3->SetLineWidth(2);

  hee3->GetYaxis()->SetRangeUser(-0.4,0.2);
   //  hee3->GetYaxis()->SetRangeUser(-0.2,1.74);
  hee3->GetXaxis()->SetTitle(hee3->GetXaxis()->GetTitle());
  hee3->GetXaxis()->SetLabelSize(0.12);
  hee3->GetXaxis()->SetTitleSize(0.14);
  hee3->GetXaxis()->SetTitleOffset(0.8);
  hee3->GetYaxis()->SetLabelSize(0.12);
  hee3->GetYaxis()->SetNdivisions(5);
  hee3->GetYaxis()->SetTitleSize(0.12);
  //gr->GetYaxis()->SetTitle("Pred. rel. bias");
  hee3->GetYaxis()->SetTitle("Weight");
  hee3->GetYaxis()->SetTitleOffset(0.3);
  hee3->GetXaxis()->SetRangeUser(25.01,1000.01); 

  leg->AddEntry(hee3,"ee","L");
  leg->AddEntry(hmm3,"#mu#mu","L");
  leg->Draw();
  leg->SetNColumns(4);
   
  // return;
  
  // printWeight(hee1);
  // printWeight(hee1);

  // return;
  // TH1F *n_hee1=hee1->Clone();
  // gDirectory->ls(); return;
  
  // TFile *outfile=new TFile("photonWeights_Spring15_ght.root","recreate");
  //  TDirectory *dir=outfile->mkdir("photonWgt_graphs",""); dir->cd();
  
  // Graphs
  TCanvas *gff=getanotherCanvas("gff");
  gff->Divide(3,1);

  gff->cd(1);
  gPad->SetGridy(); gPad->SetGridx();
  
  // GRaph1
  TGraphErrors *gr_e1=new TGraphErrors(hee1); gr_e1->SetName("eeeq0jets_qt_datafitwgts");
  TGraphErrors *gr_m1=new TGraphErrors(hmm1); gr_m1->SetName("mumueq0jets_qt_datafitwgts");
  
  
  TLegend *leg=legend();
  leg->SetHeader("0jet category");
  leg->AddEntry(gr_e1,"ee","LP");
  leg->AddEntry(gr_m1,"#mu#mu","LP");

  gr_e1->Draw("AP");
  gr_e1->SetMarkerSize(0.1); gr_e1->SetMarkerColor(4);
  gr_m1->Draw("P");
  gr_m1->SetMarkerSize(0.1); gr_m1->SetMarkerColor(2);

  leg->Draw("SAME");

  gff->cd(2);
  gPad->SetGridy(); gPad->SetGridx();

  // GRaph1
  TGraphErrors *gr_e2=new TGraphErrors(hee2);gr_e2->SetName("eegeq1jets_qt_datafitwgts");
  TGraphErrors *gr_m2=new TGraphErrors(hmm2);gr_m2->SetName("mumugeq1jets_qt_datafitwgts");
  
  TLegend *leg=legend();
  leg->SetHeader(">=1jet category");
  leg->AddEntry(gr_e2,"ee","LP");
  leg->AddEntry(gr_m2,"#mu#mu","LP");

  gr_e2->Draw("AP");
  gr_e2->SetMarkerSize(0.2); gr_e2->SetMarkerColor(4);
  gr_m2->Draw("P");
  gr_m2->SetMarkerSize(0.2); gr_m2->SetMarkerColor(2);

  leg->Draw("SAME");


  gff->cd(3);
  gPad->SetGridy(); gPad->SetGridx();

// GRaph1
  TGraphErrors *gr_e3=new TGraphErrors(hee3);gr_e3->SetName("eevbf_qt_datafitwgts");
  TGraphErrors *gr_m3=new TGraphErrors(hmm3);gr_m3->SetName("mumuvbf_qt_datafitwgts");
  
  TLegend *leg=legend();
  leg->SetHeader("VBF category");
  leg->AddEntry(gr_e3,"ee","LP");
  leg->AddEntry(gr_m3,"#mu#mu","LP");

  gr_e3->Draw("AP");
  gr_e3->SetMarkerSize(0.2); gr_e3->SetMarkerColor(4);
  gr_m3->Draw("P");
  gr_m3->SetMarkerSize(0.2); gr_m3->SetMarkerColor(2);

  leg->Draw("SAME");

  TFile *outfile;

  if (isData) {
    outfile=new TFile("photonWeights_Run2015Data.root","recreate");
  } else {
    outfile=new TFile("photonWeights_Spring15MC.root","recreate");
  }
  
  gr_e1->Write(); gr_m1->Write();
  gr_e2->Write(); gr_m2->Write();
  gr_e3->Write(); gr_m3->Write();

  outfile->Write(); outfile->Close();
  
}






