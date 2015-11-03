#ifndef _readtree_hh_
#define _readtree_hh_

#include <TString.h>
#include <TGraph.h>
#include <TH1F.h>

#include <map>
#include <vector>

enum FlavourSplitting {NOFLAVOURSPLITTING=0, UDSGSPLITTING=1, CSPLITTING=4, BSPLITTING=5 };
enum GenWeightMode { NOGENWGT=0, GENWEIGHT=1 };

Int_t getSecVtxBinToFill(Double_t firstVtxMass, Double_t secondVtxMass,Int_t nJets);
Int_t getBtagCatBinToFill(Int_t nBtags, Int_t nJets);
void ReadTree(TString filename,
	      TString outname,
	      Int_t channelSelection=13,
	      Int_t chargeSelection=0,
	      TH1F* normH=0,
	      Bool_t isTTbar=false,
	      FlavourSplitting flavourSplitting=NOFLAVOURSPLITTING,
	      GenWeightMode genWgtMode=NOGENWGT,
	      TGraph* puWgtGr=NULL, TGraph* puUpWgtGr=NULL, TGraph* puDownWgtGr=NULL);
std::map<Int_t,Double_t> lumiPerRun();
std::vector<double> getJetResolutionScales(double pt, double eta, double genjpt);
std::vector<double> getLeptonSelectionScaleFactor(int l_id, double l_pt, double l_eta, bool isData);

#endif
