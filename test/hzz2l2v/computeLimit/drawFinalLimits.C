

int theMassPoints[12] = {300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 2500, 3000};

TCutG* GetErrorBand(string name, TGraph* Low, TGraph* High)
{
   TCutG* cutg = new TCutG(name.c_str(),Low->GetN()+High->GetN()+2);
   cutg->SetFillColor(kGreen-7);
   cutg->SetLineStyle(0);
   cutg->SetLineColor(0);
   int I = 0;
   for(int i=0;i<Low->GetN();i++){  cutg->SetPoint(I,Low ->GetX()[i]               , Low ->GetY()[i]               );I++; }
                                    cutg->SetPoint(I,Low ->GetX()[Low ->GetN()-1]  , Low ->GetY()[Low ->GetN()-1]  );I++;
                                    cutg->SetPoint(I,High->GetX()[High->GetN()-1]  , High->GetY()[High->GetN()-1]  );I++;
   for(int i=0;i<High->GetN() ;i++){cutg->SetPoint(I,High->GetX()[High->GetN()-1-i], High->GetY()[High->GetN()-1-i]);I++;}
   return cutg;
}

TGraph *createTheLimitTGraph(TString inputFile,TString name, float limitType){
  cout << "will do " << name << endl;
  TChain *chain = new TChain("limit");
  chain->Add(inputFile);
  float x[12], y[12];
  int nbOfPoints=0;
  for (int i=0; i<12 ; i++){
    TH1F *expectedLimit = new TH1F("expected","",30000,0,10000);
    TString theCut = Form("mh==%i&&abs(quantileExpected-%f)<0.01",theMassPoints[i],limitType);
    chain->Draw("limit>>expected",theCut);
    float theLimit =  expectedLimit->GetMean()/1000;
    cout << "mass " <<  theMassPoints[i] << "=" << theLimit << endl;
    delete expectedLimit;
    if (theLimit>0) {
      x[nbOfPoints] = theMassPoints[i];
      y[nbOfPoints] = theLimit;
      nbOfPoints++;
    }
  }
  TGraph *theLimitGraph= new TGraph(nbOfPoints,x,y);
  theLimitGraph->SetName(name);
  return theLimitGraph;
}

void drawTheCanvas(TString type){
  TGraph *exp100 = createTheLimitTGraph("cards_SB13TeV_SM_cp100.00_brn0.00/LimitTree_"+type+"_blinded.root", "exp100", 0.5);
  TGraph *obs100 = createTheLimitTGraph("cards_SB13TeV_SM_cp100.00_brn0.00/LimitTree_"+type+".root", "obs100", -1);

  TGraph *exp100_1Sdown = createTheLimitTGraph("cards_SB13TeV_SM_cp100.00_brn0.00/LimitTree_"+type+"_blinded.root", "exp100_1Sdown", 0.84);
  TGraph *exp100_1Sup = createTheLimitTGraph("cards_SB13TeV_SM_cp100.00_brn0.00/LimitTree_"+type+"_blinded.root", "exp100_1Sup", 0.16);

  TGraph *exp100_2Sdown = createTheLimitTGraph("cards_SB13TeV_SM_cp100.00_brn0.00/LimitTree_"+type+"_blinded.root", "exp100_2Sdown", 0.975);
  TGraph *exp100_2Sup = createTheLimitTGraph("cards_SB13TeV_SM_cp100.00_brn0.00/LimitTree_"+type+"_blinded.root", "exp100_2Sup", 0.025);

  TCutG* TGExpLimit1S  = GetErrorBand("1S", exp100_1Sup, exp100_1Sdown);
  TCutG* TGExpLimit2S  = GetErrorBand("2S", exp100_2Sup, exp100_2Sdown);
  TGExpLimit2S->SetFillColor(5);

  TGraph *exp10 = createTheLimitTGraph("cards_SB13TeV_SM_cp10.00_brn0.00/LimitTree_"+type+"_blinded.root", "exp10", 0.5);
  TGraph *obs10 = createTheLimitTGraph("cards_SB13TeV_SM_cp10.00_brn0.00/LimitTree_"+type+".root", "obs10", -1);

  TGraph *exp5 = createTheLimitTGraph("cards_SB13TeV_SM_cp5.00_brn0.00/LimitTree_"+type+"_blinded.root", "exp5", 0.5);
  TGraph *obs5 = createTheLimitTGraph("cards_SB13TeV_SM_cp5.00_brn0.00/LimitTree_"+type+".root", "obs5", -1);

  TCanvas *c0 = new TCanvas("c0","coucou",800,800);
  c0->SetGridx();
  c0->SetGridy();
  c0->SetLogy();


  exp100->SetMinimum(0.001);
  exp100->SetLineWidth(2);
  exp100->SetLineColor(kBlack);
  exp100->SetLineStyle(2);
  exp100->Draw("AC");
 
 	//X/Y axis legend
 	string prod = "";
 	if(type == "ggH") prod = "gg";
 	if(type == "qqH") prod = "qq";
 	if(type == "ppH") prod = "pp";
 	exp100->GetYaxis()->SetTitle((string("#sigma_{95%} (") + prod +" #rightarrow H #rightarrow ZZ) (pb)").c_str());
 	exp100->GetYaxis()->SetTitleOffset(1.40);
 	exp100->GetXaxis()->SetTitle("M_{H} [GeV]");
 	exp100->GetXaxis()->SetTitleOffset(1.20);
  exp100->GetYaxis()->SetRangeUser(1E-3,6);

  /*exp100_1Sdown->SetLineWidth(2);
  exp100_1Sdown->SetLineColor(kGreen);
  exp100_1Sdown->Draw("same");

  exp100_1Sup->SetLineWidth(2);
  exp100_1Sup->SetLineColorkBlue);
  exp100_1Sup->Draw("same");*/
  TGExpLimit2S->Draw("fc same");
  TGExpLimit1S->Draw("fc same");

  exp100->Draw("same");
  obs100->SetLineWidth(2);
  obs100->SetLineColor(kBlack);
  obs100->Draw("same");

  exp10->SetLineWidth(2);
  exp10->SetLineColor(kRed);
  exp10->SetLineStyle(2);
  exp10->Draw("same");

  obs10->SetLineWidth(2);
  obs10->SetLineColor(kRed);
  obs10->Draw("same");

  exp5->SetLineWidth(2);
  exp5->SetLineColor(kBlue);
  exp5->SetLineStyle(2);
  exp5->Draw("same");

  obs5->SetLineWidth(2);
  obs5->SetLineColor(kBlue);
  obs5->Draw("same");

	TLegend* LEG = new TLegend(0.35,0.73,0.60,0.93);
  LEG->SetHeader("");
  LEG->SetFillColor(0);
  LEG->SetFillStyle(0);
  LEG->SetTextFont(42);
  LEG->SetBorderSize(0);
  LEG->AddEntry(exp100  , "median expected"  ,"L");
  LEG->AddEntry(TGExpLimit1S  , "expected #pm 1#sigma"  ,"F");
  LEG->AddEntry(TGExpLimit2S  , "expected #pm 2#sigma"  ,"F");
	LEG->AddEntry(obs100  , "observed"  ,"LP");
  LEG->Draw();

  TLegend* LEG2 = new TLegend(0.60,0.73,0.85,0.93);
  LEG2->SetHeader("");
  LEG2->SetFillColor(0);
  LEG2->SetFillStyle(0);
  LEG2->SetTextFont(42);
  LEG2->SetBorderSize(0);
  LEG2->AddEntry(obs5  , "#Gamma = 5 GeV"  ,"L");
  LEG2->AddEntry(obs10  , "#Gamma = 10 GeV"  ,"L");
  LEG2->AddEntry(obs100  , "#Gamma = 100 GeV"  ,"L");
  LEG2->Draw();



  c0->Print("limit_"+type+".png");
}

void drawFinalLimits(){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  drawTheCanvas("ggH");
  drawTheCanvas("qqH");
  drawTheCanvas("ppH");
}
