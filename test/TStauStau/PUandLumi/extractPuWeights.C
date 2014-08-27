#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TH1F.h>

void extractPuWeights(TString pufile){
  
  TFile* in = new TFile(pufile,"READ");
  TH1F*  pu = (TH1F*) in->Get("pileup");
  
  ofstream out;
  out.open("pileupVector.txt");
  out << "Non normalized array: " << endl;
  out << "{ ";
  for(int i=1; i<=pu->GetNbinsX(); ++i) out << pu->GetBinContent(i) << ", ";  
  out << " }" << endl;
  pu->Scale(1./pu->Integral());
  out << "Normalized array: " << endl;
  out << "{ ";
  int k(0);
  for(int i=1; i<=pu->GetNbinsX(); ++i){ k++; out << pu->GetBinContent(i) << ", "; }
  out << " }" << endl;
  out << "k is " << k << endl; 
  out.close();
} 
