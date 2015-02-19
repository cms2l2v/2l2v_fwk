// <author>Cristóvão B da Cruz e Silva</author>
// <email>c.beirao@cern.ch</email>
// <date>2014-09-30</date>
// <summary>This is a template tool that allows a user to add a variable to a summary tree/ntuple</summary>

// <description>
// In the current functional example code a collection of leptons and a collection of taus are retrieved.
// From these collections, the selected lepton and tau are chosen and the mass of their originating particle is computed with the SVfit algorithm.
// The mass is then saved in the new tree.
//
// To add your own variables to your own trees/ntuples, edit this file by adding the desired capability.
// In order to do this, search for the "EDIT HERE" string, which indicates the minimum number of places where code needs to be added for a functional executable.
// I recommend comenting out/deleting the current example code unless you need it since it is fairly slow and depends on the fact that you have branches in the input tree/ntuple with the same name. This mentioned code is also found near the "EDIT HERE" tags.
//
// To run on a full collection of trees/nTuples, use the script: [ToDo]
// </description>


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

  std::string ttree = "Events";
  std::string outFile;
  std::string inFile;

  // Parse the command line options
  for(int i = 1; i < argc; ++i)
  {
    std::string arg(argv[i]);

    if(arg.find("--help") != std::string::npos)
    {
      printHelp();
      return 0;
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

  // Define some control variables, variables that need to be read from the tree or the new variable(s) to be computed
  llvvLeptonCollection *selLeptons = 0;
  llvvTauCollection *selTaus = 0;
  bool selected = false;
  int tauIndex = -1, leptonIndex = -1;
  double SVFitMass = -1, SVFitMass_old = -1;
  llvvMet* met = 0;
  // EDIT HERE

  // The input tree is read from the input file
  TFile  in(inFile.c_str(), "READ");
  TTree* inTree = (TTree*)in.Get(ttree.c_str());
  Long64_t nentries = inTree->GetEntries();
  inTree->SetAutoSave(0);

  // The input tree is copied to the output file
  TFile out(outFile.c_str(), "RECREATE");
  TTree* outTree = inTree->CopyTree("");
  outTree->SetAutoSave(0);

  // Link the variables that are to be read to their branches in the input tree
  inTree->SetBranchAddress("tauIndex", &tauIndex);
  inTree->SetBranchAddress("leptonIndex", &leptonIndex);
  inTree->SetBranchAddress("selected", &selected);
  inTree->SetBranchAddress("selLeptons", &selLeptons);
  inTree->SetBranchAddress("selTaus", &selTaus);
  inTree->SetBranchAddress("met", &met);
  inTree->SetBranchAddress("SVFitMass", &SVFitMass_old);
  // EDIT HERE

  // Create the new branch in the new tree for the new variable
  TBranch* SVFitBranch = outTree->Branch("SVFitMass_final", &SVFitMass);
  // EDIT HERE

  for(Long64_t i = 0; i < nentries; ++i)
  {
    inTree->GetEntry(i);
    // WARNING: Never ever ever use a continue that will skip events on this outer loop! YOU HAVE BEEN WARNED

    // For each entry in the tree compute the new variable
    // EDIT HERE
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
      // Remember that when using a selection, as in this case, you should define a default value for the variable and always fill the new branch or else the new branch will not have the same number of entries as the previous branches, loosing the association of which event each entry in the new branch beongs to.


    // Fill the branch with the new variable
    SVFitBranch->Fill();
    // EDIT HERE
  }

  outTree->Write("", TObject::kOverwrite);

  out.Close();
  return 0;
}

void printHelp()
{
  std::cout << "addSVfit help - There are the following options:" << std::endl << std::endl;

  std::cout << "--help    -->  Print this help message" << std::endl;
  std::cout << "--inFile  -->  Path to the file to modify" << std::endl;
  std::cout << "--outFile -->  Path to where to save the new file with SVfit added" << std::endl;
  std::cout << "--ttree   -->  Name of the ttree with the list of events (default: Events)" << std::endl;

  return;
}
