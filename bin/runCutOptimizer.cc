// <author>Cristóvão B da Cruz e Silva</author>
// <email>c.beirao@cern.ch</email>
// <date>2014-06-02</date>
// <summary>Script that runs an iterative cut optimization procedure. Cut values are scanned and the value with the maximum FOM is used.</summary>

#include <iostream>
#include <string>

#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TTree.h"

#include "UserCode/llvv_fwk/interface/JSONWrapper.h"

// REMEMBER TO IGNORE ANY DATA DEFINED IN THE JSON, YOU SHOULD NEVER OPTIMIZE CUTS ON DATA
// IT MIGHT BE OK TO USE THE DATA IF USING A METHOD WITH A DATA DRIVEN BACKGROUND ESTIMATION (CONFIRM THIS BEFORE USING IT HERE)

void printHelp();

int main(int argc, char** argv)
{
  // Parse the command line options
  for(int i = 1; i < argc; ++i)
  {
    std::string arg(argv[i]);

    // First see if help info is being requested, if so, print help and then exit
    if(arg.find("--help") != std::string::npos)
    {
      printHelp();
      return 0;
    }
  }
}

void printHelp()
{
  std::cout << "runCutOptimizer help - There are the following options:" << std::endl << std::endl;

  std::cout << "--help    -->  Print this help message" << std::endl;
  std::cout << "--json    -->  File with list of processes to include, follows the same conventions as other jsons in this framework" << std::endl;
  std::cout << "--iLumi   -->  Integrated luminosity to be used for the MC rescale, given in pb-1" << std::endl;
  std::cout << "--inDir   -->  Path to the directory containing the summary trees" << std::endl;
  std::cout << "--outDir  -->  Path to the directory where to output plots and tables (will be created if it doesn't exist)" << std::endl;
  std::cout << "--cutVars -->  File with list of variables with which to perform cut optimization. Check ... for an example syntax" << std::endl; // TODO, place her my example syntax

  std::cout << std::endl << "Example command:" << std::endl << "\trunCutOptimizer --json analysis_samples.json --iLumi 8000 --cutVars variables.json --inDir /path/to/merged/nTuples --outDir ./OUT/" << std::endl;
  return;
}

