{
//Superimposes model with and without EW corrections ; TProfiles for EWK corrections as a function of mZZ ans ptZ.
cout<<"test1"<<endl;
TFile f1("./plotter.root");
f1.cd("ZZ#rightarrow 2l2#nu");
cout<<"test2"<<endl;

TH1F* LO_MZZ = (TH1F*) gDirectory->Get("LO_Nevent_vs_Mzz");
cout<<"test3"<<endl;
TH1F* NLO_MZZ = (TH1F*) gDirectory->Get("NLO_Nevent_vs_Mzz");
cout<<"test4"<<endl;
TH1F* LO_PTZ = (TH1F*) gDirectory->Get("ll_LO_Nevent_vs_ZpT");
cout<<"test5"<<endl;
TH1F* NLO_PTZ = (TH1F*) gDirectory->Get("ll_NLO_Nevent_vs_ZpT");
cout<<"test6"<<endl;

LO_MZZ->Rebin(10);
NLO_MZZ->Rebin(10);
LO_PTZ->Rebin(10);
NLO_PTZ->Rebin(10);

TCanvas* can = new TCanvas("can","EWK for mZZ",1600,950);
gStyle->SetOptStat(0);
can->Divide(2);
can->cd(1);
gPad->SetLogy();

for(int i =1 ; i < LO_MZZ->GetNbinsX()+1 ; i++){
if (NLO_MZZ->GetBinContent(i) > LO_MZZ->GetBinContent(i)) cout << "Problem at bin " << i << endl;
if(i==LO_MZZ->GetNbinsX())cout << "Total size : " << i << endl;
}

TGraphAsymmErrors *Corr_MZZ = new TGraphAsymmErrors(1);
Corr_MZZ->Divide(NLO_MZZ, LO_MZZ);

legend = new TLegend(0.6,0.7,0.89,0.89);
legend->SetFillColor(kWhite);
legend->SetLineColor(kWhite);

NLO_MZZ->SetLineColor(2);

legend->AddEntry(LO_MZZ,"EW LO","L");
legend->AddEntry(NLO_MZZ,"EW NLO","L");

LO_MZZ->Scale(1./LO_MZZ->Integral());
NLO_MZZ->Scale(1./NLO_MZZ->Integral());

LO_MZZ->SetTitle("M_{ZZ} for LO and NLO EWK");
LO_MZZ->GetYaxis()->SetTitle("N/N_{tot}/bin");

LO_MZZ->DrawCopy("hist");
NLO_MZZ->DrawCopy("histsame");
legend->Draw();

can->cd(2);

Corr_MZZ->SetTitle("EWK K-factor vs M_{ZZ}");
Corr_MZZ->GetXaxis()->SetTitle("M_{ZZ}");
Corr_MZZ->GetYaxis()->SetTitle("K");

Corr_MZZ->Draw();

can->SaveAs("EwkCorrectionsPlots/mZZ_EWK.root");
can->SaveAs("EwkCorrectionsPlots/mZZ_EWK.png");
can->SaveAs("EwkCorrectionsPlots/mZZ_EWK.pdf");


TCanvas* can2 = new TCanvas("can2","EWK for ptZ",1600,950);

gStyle->SetOptStat(0);
can2->Divide(2);
can2->cd(1);
gPad->SetLogy();


for(int i =1 ; i < LO_PTZ->GetNbinsX()+1 ; i++){
if (NLO_PTZ->GetBinContent(i) > LO_PTZ->GetBinContent(i)) cout << "Problem at bin " << i << endl;
if(i==LO_PTZ->GetNbinsX())cout << "Total size : " << i << endl;
}

TGraphAsymmErrors *Corr_PTZ = new TGraphAsymmErrors(1);
Corr_PTZ->Divide(NLO_PTZ, LO_PTZ);

NLO_PTZ->SetLineColor(2);

LO_PTZ->Scale(1./LO_PTZ->Integral());
NLO_PTZ->Scale(1./NLO_PTZ->Integral());

LO_PTZ->SetTitle("p_{T,Z} for LO and NLO EWK");
LO_PTZ->GetYaxis()->SetTitle("N/N_{tot}/bin");

LO_PTZ->DrawCopy("hist");
NLO_PTZ->DrawCopy("histsame");
legend->Draw();

can2->cd(2);

Corr_PTZ->SetTitle("EWK K-factor vs p_{T,Z}");
Corr_PTZ->GetXaxis()->SetTitle("P_{T,Z}");
Corr_PTZ->GetYaxis()->SetTitle("K");

Corr_PTZ->Draw();

can2->SaveAs("EwkCorrectionsPlots/pTZ_EWK.root");
can2->SaveAs("EwkCorrectionsPlots/pTZ_EWK.png");
can2->SaveAs("EwkCorrectionsPlots/pTZ_EWK.pdf");

}
