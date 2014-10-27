// <author>Cristóvão B da Cruz e Silva</author>
// <email>c.beirao@cern.ch</email>
// <date>2014-09-30</date>

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <cassert>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <fstream>

#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TTree.h"
#include "TCut.h"
#include "TEventList.h"
#include "TInterpreter.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include "UserCode/llvv_fwk/interface/JSONWrapper.h"
#include "UserCode/llvv_fwk/interface/llvvObjects.h"

#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h"

void printHelp();

//
// Output codes:
// 0 - Everything has run fine
// 1 - Invalid arguments
// 2 - Problem parsing the arguments
//
int main(int argc, char** argv)
{
  AutoLibraryLoader::enable();

  std::string jsonFile;
  std::string ttree = "Events";
  std::string outFile;
  std::string inFile;
//  std::vector<std::string> plotExt;
  bool verbose = false;

  // Parse the command line options
  for(int i = 1; i < argc; ++i)
  {
    std::string arg(argv[i]);

    if(arg.find("--help") != std::string::npos)
    {
      printHelp();
      return 0;
    }

    if(arg.find("--json") != std::string::npos)
    {
      jsonFile = argv[i+1];
      ++i;
    }

    if(arg.find("--ttree") != std::string::npos)
    {
      ttree = argv[i+1];
      ++i;
    }

    if(arg.find("--inFile") != std::string::npos)
    {
      inFile = argv[i+1];
      ++i;
    }

    if(arg.find("--outFile") != std::string::npos)
    {
      outFile = argv[i+1];
      ++i;
    }

    if(arg.find("--verbose") != std::string::npos)
    {
      verbose = true;
    }
  }
  if(inFile == "")
  {
    std::cout << "You should define the input file" << std::endl;
    return 2;
  }
  if(outFile == "")
  {
    std::cout << "You should define where to save the modified file" << std::endl;
    return 2;
  }

  llvvLeptonCollection *selLeptons = 0;
  llvvTauCollection *selTaus = 0;
  bool selected = false;
  int tauIndex = -1, leptonIndex = -1;
  double SVFitMass = -1, SVFitMass_old = -1;
  llvvMet* met = 0;

  TFile  in(inFile.c_str(), "READ");
  TTree* inTree = (TTree*)in.Get(ttree.c_str());
  Long64_t nentries = inTree->GetEntries();
  inTree->SetAutoSave(0);

  TFile out(outFile.c_str(), "RECREATE");
  TTree* outTree = inTree->CopyTree("");
  outTree->SetAutoSave(0);

//  TBranch* new_SVfitBranch = (TBranch*)outTree->GetBranch("SVFitMass")->Clone("SVFitMass_final");

  inTree->SetBranchAddress("tauIndex", &tauIndex);
  inTree->SetBranchAddress("leptonIndex", &leptonIndex);
  inTree->SetBranchAddress("selected", &selected);
  inTree->SetBranchAddress("selLeptons", &selLeptons);
  inTree->SetBranchAddress("selTaus", &selTaus);
  inTree->SetBranchAddress("met", &met);// */
  inTree->SetBranchAddress("SVFitMass", &SVFitMass_old);

/*  if(outTree->FindBranch("SVFitMass"))
  {
    std::cout << "Found a pre-existing SVFitMass Branch, deleting it!" << std::endl;
    TBranch* b = (TBranch*)outTree->GetBranch("SVFitMass");
    outTree->GetListOfBranches()->Remove(b);
  }// */

  TBranch* SVFitBranch = outTree->Branch("SVFitMass_final", &SVFitMass);

  for(Long64_t i = 0; i < nentries; ++i)
  {
    inTree->GetEntry(i);

    if(selected && selLeptons->size() > 0 && selTaus->size() > 0)
    {
      auto selLepton = (*selLeptons)[leptonIndex];
      auto selTau    = (*selTaus)[tauIndex];

      TMatrixD covMET(2, 2);
      std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
      svFitStandalone::Vector measuredMET(met->px(), met->py(), 0);
      covMET[0][0] = met->sigx2;
      covMET[0][1] = met->sigxy;
      covMET[1][0] = met->sigxy;
      covMET[1][1] = met->sigy2;

/*      if(verbose)
      {
        std::cout << "met:" << std::endl;
        std::cout << "  sigx2: " << met->sigx2 << std::endl;
        std::cout << "  sigy2: " << met->sigy2 << std::endl;
        std::cout << "  sigxy: " << met->sigxy << std::endl;
        std::cout << "lepton:" << std::endl;
        std::cout << "  id: " << selLepton.id << std::endl;
        std::cout << "  px: " << selLepton.px() << std::endl;
        std::cout << "  py: " << selLepton.py() << std::endl;
        std::cout << "  pz: " << selLepton.pz() << std::endl;
        std::cout << "  E : " << selLepton.E()  << std::endl;
        std::cout << "  M : " << selLepton.mass() << std::endl;
        std::cout << "tau:" << std::endl;
        std::cout << "  px: " << selTau.px() << std::endl;
        std::cout << "  py: " << selTau.py() << std::endl;
        std::cout << "  pz: " << selTau.pz() << std::endl;
        std::cout << "  E : " << selTau.E()  << std::endl;
        std::cout << "  M : " << selTau.mass() << std::endl;
      }// */

      if(abs(selLepton.id) == 11)
        measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToElecDecay, svFitStandalone::LorentzVector(selLepton.px(), selLepton.py(), selLepton.pz(), selLepton.E())));
      else
        measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToMuDecay, svFitStandalone::LorentzVector(selLepton.px(), selLepton.py(), selLepton.pz(), selLepton.E())));
      measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToHadDecay, svFitStandalone::LorentzVector(selTau.px(), selTau.py(), selTau.pz(), selTau.E())));

      SVfitStandaloneAlgorithm SVfit_algo(measuredTauLeptons, measuredMET, covMET, 0);
      SVfit_algo.addLogM(false);
      SVfit_algo.integrateVEGAS();
      if(SVfit_algo.isValidSolution())
        SVFitMass = SVfit_algo.mass();
      else
        SVFitMass = -1;
    }
    else
      SVFitMass = -1;

    SVFitBranch->Fill();
  }

  outTree->Write("", TObject::kOverwrite);

  out.Close();
  return 0;
}

//SVFitMass

void printHelp()
{
  std::cout << "addSVfit help - There are the following options:" << std::endl << std::endl;

  std::cout << "--help    -->  Print this help message" << std::endl;
//  std::cout << "--json    -->  Configuration file, should define which files to run on" << std::endl;
  std::cout << "--inFile  -->  Path to the file to modify" << std::endl;
  std::cout << "--outFile -->  Path to where to save the new file with SVfit added" << std::endl;
  std::cout << "--ttree   -->  Name of the ttree with the list of events (default: Events)" << std::endl;
//  std::cout << "--plotExt -->  Extension format with which to save the plots, repeat this command if multiple formats are desired" << std::endl;
//  std::cout << "--verbose -->  Set to verbose mode, the current step will be printed out to screen" << std::endl;

//  std::cout << std::endl << "Example command:" << std::endl << "\trunCutOptimizer --json optimization_options.json --outFile ./OUT/" << std::endl;
  return;
}
