#include "TCanvas.h"
#include "TPaveText.h"
#include "TString.h"
#include "TH1F.h"
#include "TFile.h"

#include <map>

using namespace std;

void showHeavyFlavorClosureTestResults(TString url="~/work/top_5311/hfc_closure_-1.0/HFCClosureTest.root")
{
  std::map<std::string, TH1F *> closureTests;
  closureTests["ee"]=0;
  closureTests["mumu"]=0;
  closureTests["emu"]=0;
  closureTests["inclusive"]=0;


  TFile *inF=TFile::Open(url);
  for(std::map<std::string, TH1F *>::iterator rit=closureTests.begin(); rit!=closureTests.end(); rit++)
    {
      rit->second=(TH1F *)inF->Get( ("bias_"+rit->first).c_str() );
      rit->second->SetDirectory(0);
    }
  inF->Close();

  
  
  //show results
  TCanvas *c = new TCanvas("c","c",1000,1000);
  c->Divide(2,2);
  int icat=0;
  for(std::map<std::string, TH1F *>::iterator rit=closureTests.begin(); rit!=closureTests.end(); rit++,icat++)
    {
      string tag=rit->first;
      c->cd(icat+1);
      rit->second->Draw("e1");
      rit->second->GetXaxis()->SetNdivisions(5);
      rit->second->GetYaxis()->SetTitleOffset(1.0);
      rit->second->Fit("gaus","L+","e1same");
      
      //category
      TPaveText *pt = new TPaveText(0.2,0.85,0.45,0.94,"brNDC");
      pt->SetBorderSize(0);
      pt->SetFillColor(0);
      pt->SetFillStyle(0);
      pt->SetTextFont(42);
      TString caption=tag.c_str();
      caption=caption.ReplaceAll("mu","#mu");
      pt->AddText(caption);
      pt->Draw();
      
      //overall header
      if(icat==0)
	{
	  pt = new TPaveText(0.1,0.96,0.9,1.0,"brNDC");
	  pt->SetBorderSize(0);
	  pt->SetFillColor(0);
	  pt->SetFillStyle(0);
	  pt->SetTextAlign(12);
	  pt->AddText("CMS simulation, #sqrt{s}=8 TeV");
	  pt->Draw();
	}
    }

  c->SaveAs("HFCClosureTest.png");
  c->SaveAs("HFCClosureTest.pdf");
  c->SaveAs("HFCClosureTest.C");
  
}
