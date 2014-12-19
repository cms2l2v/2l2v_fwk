#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TCanvas.h"

#include <iostream>
#include <sstream>
#include <string>


void makeCountour()
{
  gStyle->SetOptStat(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetMarkerStyle(15);
  gStyle->SetMarkerSize(0.25);
  gStyle->SetTextFont(42);
  gStyle->SetMarkerColor(37);


  int minStauM = 50;
  int maxStauM = 500;
  int stepStauM = 10;

  int minNeutM = 0;
  int maxNeutM = 480;
  int stepNeutM = 10;


  TCanvas c1("c1", "c1", 800, 600);
  TFile outFile("TStauStau_limit.root","RECREATE");

  TH2D histo2D("twodplot", "twodplot", (maxStauM-minStauM)/stepStauM + 1, minStauM - stepStauM/2, maxStauM + stepStauM/2, (maxNeutM-minNeutM)/stepNeutM + 1, minNeutM - stepNeutM/2, maxNeutM + stepNeutM/2);

  for(int stauM = minStauM; stauM <= maxStauM; stauM += stepStauM)
  {
    for(int neutM = minNeutM; neutM <= maxNeutM; neutM += stepNeutM)
    {
      std::stringstream temp;
      temp << "./Results/";
      temp << "higgsCombineS" << stauM << "-N" << neutM << ".Asymptotic.mH120.root";

      std::string fileName;
      temp >> fileName;
      std::cout << fileName << std::endl;

      TFile file(fileName.c_str(), "READ");

      if(file.IsOpen())
      {
        TTree *treelimit = (TTree*)file.Get("limit");

        TH1D* obs = new TH1D("obs","",100,0,100);
        treelimit->Draw("limit>>obs", "quantileExpected==-1");
        TH1D* expectedm2 = new TH1D("expectedm2","",100,0,100);
        treelimit->Draw("limit>>expectedm2", "quantileExpected>0.02 && quantileExpected<0.03");
        TH1D* expectedm1 = new TH1D("expectedm1","",100,0,100);
        treelimit->Draw("limit>>expectedm1", "quantileExpected>0.15 && quantileExpected<0.16");
        TH1D* expected = new TH1D("expected","",100,0,100);
        treelimit->Draw("limit>>expected", "quantileExpected==0.5");
        TH1D* expectedp1 = new TH1D("expectedp1","",100,0,100);
        treelimit->Draw("limit>>expectedp1", "quantileExpected>0.83 && quantileExpected<0.84");
        TH1D* expectedp2 = new TH1D("expectedp2","",100,0,100);
        treelimit->Draw("limit>>expectedp2", "quantileExpected>0.97 && quantileExpected<0.98");// */

        double limit = expected->GetMean();
        if(limit <= 0)
          limit = 2e8;
        else
        {
          limit /= 30;
        }

        std::cout << stauM << ", " << neutM << ": " << limit << std::endl;

        histo2D.Fill(stauM, neutM, limit);

        file.Close();
      }
      else
      {
//        histo2D.Fill(stauM, neutM, 500);
      }
    }
  }

  outFile.cd();

  double contours[2];
  contours[0] = 2;
  contours[1] = 10000;

  int colors[2] = {2,4}; //red, blue,black
  gStyle->SetPalette(2,colors);

  histo2D.SetContour(2);
  histo2D.SetContourLevel(0,1); //value for your first level
  histo2D.SetContourLevel(1,1e8); //non-existing high level
  histo2D.SetFillColor(2);
  histo2D.Draw("cont3");

  histo2D.Write("twodplot");

  TH2D *histo2D2=(TH2D*) histo2D.Clone("histo2Dclone");
  histo2D2->SetContour(2);
  histo2D2->SetContourLevel(0,-1e3); //non existing low level
  histo2D2->SetContourLevel(1,1); //value for your second level
  histo2D2->SetFillColor(4);
  histo2D2->Draw("colz");
  histo2D2->Draw("cont3 same");

  gStyle->SetTextFont(42);
  c1.SaveAs("contour_1overmu.png");

  outFile.Write();
  outFile.Close();
}
