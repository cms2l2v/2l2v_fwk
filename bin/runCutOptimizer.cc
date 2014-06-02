// <author>Cristóvão B da Cruz e Silva</author>
// <email>c.beirao@cern.ch</email>
// <date>2014-06-02</date>
// <summary>Script that runs an iterative cut optimization procedure. Cut values are scanned and the value with the maximum FOM is used.</summary>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <cassert>
#include <unordered_map>

#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TTree.h"

#include "UserCode/llvv_fwk/interface/JSONWrapper.h"

// REMEMBER TO IGNORE ANY DATA DEFINED IN THE JSON, YOU SHOULD NEVER OPTIMIZE CUTS ON DATA
// IT MIGHT BE OK TO USE THE DATA IF USING A METHOD WITH A DATA DRIVEN BACKGROUND ESTIMATION (CONFIRM THIS BEFORE USING IT HERE)

void printHelp();
std::vector<std::pair<std::string,TChain*>> getChainsFromJSON(JSONWrapper::Object& json, std::string RootDir, std::string type="BG");


std::unordered_map<std::string,bool> FileExists;


//
// Output codes:
// 0 - Everything has run fine
// 1 - Invalid arguments
// 2 - Problem parsing the arguments
//

int main(int argc, char** argv)
{
  double iLumi = 0;
  std::string jsonFile;
  std::string inDir;
  std::string outDir = "./OUT/";
//  std::string cutVars;

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

    if(arg.find("--json") != std::string::npos)
    {
      jsonFile = argv[i+1];
      // TODO: Check if the file exists
      ++i;
    }

    if(arg.find("--iLumi") != std::string::npos)
    {
      std::stringstream converter(argv[i+1]);
      converter >> iLumi;
      if(converter.fail())
      {
        std::cout << "The integrated luminosity should be a number." << std::endl;
        std::cout << "For more information, consult the help (\"runCutOptimizer --help\")" << std::endl;
        return 2;
      }
      if(iLumi <= 0)
      {
        std::cout << "You should give a positive, non-zero, integrated luminosity." << std::endl;
        std::cout << "For more information, consult the help (\"runCutOptimizer --help\")" << std::endl;
        return 1;
      }
      ++i;
    }

    if(arg.find("--inDir") != std::string::npos)
    {
      inDir = argv[i+1];
      // TODO: Check if input directory exists
      ++i;
    }

    if(arg.find("--outDir") != std::string::npos)
    {
      outDir = argv[i+1];
      ++i;
    }

/*    if(arg.find("--cutVars") != std::string::npos)
    {
      cutVars = argv[i+1];
      // TODO: Check if file exists
      ++i;
    }*/
  }

  if(iLumi == 0 || jsonFile == "" || inDir == "")// || cutVars == "")
  {
    std::cout << "You should define at least the following arguments: iLumi, json, inDir" << std::endl;//, cutVars" << std::endl;
    std::cout << "For more information, consult the help (\"runCutOptimizer --help\")" << std::endl;
    return 2;
  }

  std::cout << "Running for an inntegrated luminosity of " << iLumi << ", using the samples defined in " << jsonFile << " from " << inDir << " directory" << std::endl;

  system(("mkdir -p " + outDir).c_str());

  std::cout << "Loading JSON" << std::endl;
  JSONWrapper::Object json(jsonFile, true);
  std::cout << "Done Loading" << std::endl;

  // Eventually change this to get the tree name from the config, like that it will not necessarily have to be named events
  std::vector<std::pair<std::string,TChain*>> BG_samples  = getChainsFromJSON(json, inDir, "BG");
  std::vector<std::pair<std::string,TChain*>> SIG_samples = getChainsFromJSON(json, inDir, "SIG");

  std::cout << "Found " << BG_samples.size()  << " background processes:" << std::endl;
  for(auto process = BG_samples.begin(); process != BG_samples.end(); ++process)
  {
    std::cout << "\t" << process->first << " with " << process->second->GetEntries() << " entries." << std::endl;
  }
  std::cout << "Found " << SIG_samples.size() << " signal processes:"     << std::endl;
  for(auto process = SIG_samples.begin(); process != SIG_samples.end(); ++process)
  {
    std::cout << "\t" << process->first << std::endl;
  }

  std::cout << "The list of ignored files, either missing or corrupt, can be found below:" << std::endl;
  for(auto key = FileExists.begin(); key != FileExists.end(); ++key)
  {
    if(!key->second)
      std::cout << "  " << key->first << std::endl;
  }
}

std::vector<std::pair<std::string,TChain*>> getChainsFromJSON(JSONWrapper::Object& json, std::string RootDir, std::string type)
{
  std::vector<std::pair<std::string,TChain*>> retVal;
  std::pair<std::string,TChain*> tempVal;
  std::vector<std::string> missingFiles;

  std::string treename = "Events"; // Change this to get the name from the json
  std::string customExtension = "_summary"; // Change this to get from the json

  if(type != "BG" && type != "SIG")
  {
    std::cout << "Unknown sample type requested." << std::endl;
    assert(type == "BG" || type == "SIG");
  }

  std::vector<JSONWrapper::Object> processes = json["proc"].daughters();
  for(auto process = processes.begin(); process != processes.end(); ++process)
  {
    bool isData = (*process)["isdata"].toBool();
    bool isSig  = !isData && (*process).isTag("spimpose") && (*process)["spimpose"].toBool();
    bool isMC   = !isData && !isSig;

    if(isData) // Here we are enforcing for the data samples to not even be present, might not make sense for a data-driven background estimation
      continue;

    if(type == "SIG" && isMC) // Here we make sure we are only processing the requested processes
      continue;
    if(type == "BG" && isSig)
      continue;

    tempVal.first = (*process).getString("tag", "Sample");
    tempVal.second = new TChain(treename.c_str(), tempVal.first.c_str());

    std::string filtExt;
    if((*process).isTag("mctruthmode"))
    {
      std::stringstream buf;
      buf << "_filt" << (*process)["mctruthmode"].toInt();
      buf >> filtExt;
    }

    std::vector<JSONWrapper::Object> samples = (*process)["data"].daughters();
    for(auto sample = samples.begin(); sample != samples.end(); ++sample)
    {
      int nFiles = (*sample).getInt("split", 1);
      for(int file = 0; file < nFiles; ++file)
      {
        std::string segmentExt;
        if(nFiles != 1)
        {
          std::stringstream buf;
          buf << "_" << file;
          buf >> segmentExt;
        }

        std::string fileName = RootDir + "/" + (*sample).getString("dtag", "") + (*sample).getString("suffix", "") + segmentExt + filtExt + customExtension + ".root";
        TFile* file = new TFile(fileName.c_str());
        bool& fileExists = FileExists[fileName];

        if(!file || file->IsZombie() || !file->IsOpen() || file->TestBit(TFile::kRecovered))
        {
          fileExists = false;
          file->Close();
          delete file;
          continue;
        }
        else
        {
          fileExists = true;
          file->Close();
          delete file;
        }

        // Chain the valid files together
        tempVal.second->Add(fileName.c_str());
      }
    }
    retVal.push_back(tempVal);
  }

  return retVal;
}

void printHelp()
{
  std::cout << "runCutOptimizer help - There are the following options:" << std::endl << std::endl;

  std::cout << "--help    -->  Print this help message" << std::endl;
  std::cout << "--json    -->  File with list of processes to include, follows the same conventions as other jsons in this framework" << std::endl;
  std::cout << "--iLumi   -->  Integrated luminosity to be used for the MC rescale, given in pb-1" << std::endl;
  std::cout << "--inDir   -->  Path to the directory containing the summary trees" << std::endl;
  std::cout << "--outDir  -->  Path to the directory where to output plots and tables (will be created if it doesn't exist)" << std::endl;
//  std::cout << "--cutVars -->  File with list of variables with which to perform cut optimization. Check ... for an example syntax" << std::endl; // TODO, place her my example syntax
// Moved cutVars to the json :D

  std::cout << std::endl << "Example command:" << std::endl << "\trunCutOptimizer --json analysis_samples.json --iLumi 8000 --cutVars variables.json --inDir /path/to/summary/trees --outDir ./OUT/" << std::endl;
  return;
}

