// <author>Cristóvão B da Cruz e Silva</author>
// <email>c.beirao@cern.ch</email>
// <date>2014-08-21</date>

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


struct SummaryInfo
{
  std::string sample;
  std::vector<std::string> OK;
  std::vector<std::string> Lost;
  std::vector<std::string> Zombie;
  std::vector<std::string> Recovered;
  std::vector<std::string> Empty;
};


void printHelp()
{
  std::cout << "checkNTuples help" << std::endl << std::endl;
  std::cout << "This tool checks a directory, using a json file that describes samples, to make sure all the nTuples have been correctly merged and are ready to be used." << std::endl << "Available options:" << std::endl;

  std::cout << "--help    -->  Print this help message" << std::endl;
  std::cout << "--json    -->  Configuration file, should define which files to run on, where they are located, the integrated luminosity and the variables to optimize the cuts on" << std::endl;
  std::cout << "--dir     -->  Path to the directory where to output plots and tables (will be created if it doesn't exist)" << std::endl;
  std::cout << "--log     -->  If this option is present, a file with the given name will be created and a detailed log output will be created in it" << std::endl;

  std::cout << std::endl << "Example command:" << std::endl << "\tcheckNTuples --json all_samples.json --dir ./mergedFiles/" << std::endl;
  return;
}

int main(int argc, char** argv)
{
  AutoLibraryLoader::enable();

  std::string jsonFile;
  std::string dir;
  std::string log;

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

    if(arg.find("--dir") != std::string::npos)
    {
      dir = argv[i+1];
      ++i;
    }

    if(arg.find("--log") != std::string::npos)
    {
      log = argv[i+1];
      ++i;
    }
  }

  if(jsonFile == "")
  {
    std::cout << "You should define at least the following arguments: json" << std::endl;
    std::cout << "For more information, consult the help (\"checkNTuples --help\")" << std::endl;
    return 2;
  }

  if(dir == "")
  {
    std::cout << "You should define at least the following arguments: dir" << std::endl;
    std::cout << "For more information, consult the help (\"checkNTuples --help\")" << std::endl;
    return 2;
  }

  JSONWrapper::Object json(jsonFile, true);
  std::vector<JSONWrapper::Object> processes = json["proc"].daughters();
  std::vector<SummaryInfo> summary;

  for(auto process = processes.begin(); process != processes.end(); ++process)
  {
    std::vector<JSONWrapper::Object> samples = (*process)["data"].daughters();
    for(auto sample = samples.begin(); sample != samples.end(); ++sample)
    {
      int split = (*sample).getInt("split", 1);
      SummaryInfo mySummary;
      mySummary.sample = (*sample).getString("dtag", "");

      for(int file_n = 0; file_n < split; ++file_n)
      {
        std::string segmentExt;
        if(split>1)
        {
          std::stringstream buf;
          buf << "_" << file_n;
          buf >> segmentExt;
        }

        std::string fileName = dir + "/" + (*sample).getString("dtag", "") + (*sample).getString("suffix", "") + segmentExt + ".root";

        TFile* file = new TFile(fileName.c_str());

        //bool fileDoesntExist = false;
        //bool fileIsZombie = false;
        //bool fileIsRecovered = false;
        //bool fileIsEmpty = false;

        if(!file || !file->IsOpen())
        {
//          fileDoesntExist = true;
          mySummary.Lost.push_back(fileName);
        }
        else
        {
          if(file->IsZombie())
          {
//            fileIsZombie = true;
            mySummary.Zombie.push_back(fileName);
          }
          else
          {
            if(file->TestBit(TFile::kRecovered))
            {
//              fileIsRecovered = true;
              mySummary.Recovered.push_back(fileName);
            }
            else
            {
              TTree* events = (TTree*)file->Get("Events");
              if(events == NULL)
              {
//                fileIsEmpty = true;
                mySummary.Empty.push_back(fileName);
              }
              else
              {
                mySummary.OK.push_back(fileName);
              }
            }
          }
        }
        file->Close();

/*        std::cout << fileName << " ";
        if(fileDoesntExist)
          std::cout << "not found.";
        else
        {
          if(fileIsZombie)
            std::cout << "is zombie.";
          else
          {
            if(fileIsRecovered)
              std::cout << "is recovered";
            else
            {
              if(fileIsEmpty)
                std::cout << "is empty";
              else
                std::cout << "seems to be fine";
            }
          }
        }
        std::cout << std::endl;*/
      }
      summary.push_back(mySummary);
    }
  }


  if(log != "")
  {
    std::ofstream out(log, std::ofstream::out);
    out << "Logfile for the command: " << argv << std::endl;
    out.close();
  }

  std::cout << "These were the results:" << std::endl;
  for(auto sum = summary.begin(); sum != summary.end(); ++sum)
  {
    std::cout << sum->sample << ": " << sum->OK.size() << " OK; " << sum->Lost.size()+sum->Zombie.size()+sum->Recovered.size()+sum->Empty.size() << " failed" << std::endl;
    if(log != "")
    {
      std::ofstream out(log, std::ofstream::out | std::ofstream::app);
      out << "Sample: " << sum->sample << std::endl;
      out << "---- Total files: " << sum->OK.size()+sum->Lost.size()+sum->Zombie.size()+sum->Recovered.size()+sum->Empty.size() << std::endl;
      out << "---- Files that seem to be OK: " << sum->OK.size() << std::endl;
      out << "---- Files that were lost (not produced or deleted): " << sum->Lost.size() << std::endl;
      out << "---- Zombie files: " << sum->Zombie.size() << std::endl;
      out << "---- Recovered files: " << sum->Recovered.size() << std::endl;
      out << "---- Empty files: " << sum->Empty.size() << std::endl;
      out << std::endl;
      out.close();
    }
  }

  if(log != "")
  {
    std::ofstream out(log, std::ofstream::out | std::ofstream::app);
    out << std::endl << "Summary of files:" << std::endl;

    for(auto sum = summary.begin(); sum != summary.end(); ++sum)
    {
      out << "-> Sample: " << sum->sample << std::endl;

      if(sum->OK.size() != 0)
      {
        out << "---- Files without any found issues:" << std::endl;
        for(auto file = sum->OK.begin(); file != sum->OK.end(); ++file)
        {
          out << "------ " << *file << std::endl;
        }
      }

      if(sum->Lost.size() != 0)
      {
        out << "---- Lost files:" << std::endl;
        for(auto file = sum->Lost.begin(); file != sum->Lost.end(); ++file)
        {
          out << "------ " << *file << std::endl;
        }
      }

      if(sum->Zombie.size() != 0)
      {
        out << "---- Zombie files:" << std::endl;
        for(auto file = sum->Zombie.begin(); file != sum->Zombie.end(); ++file)
        {
          out << "------ " << *file << std::endl;
        }
      }

      if(sum->Recovered.size() != 0)
      {
        out << "---- Recovered files:" << std::endl;
        for(auto file = sum->Recovered.begin(); file != sum->Recovered.end(); ++file)
        {
          out << "------ " << *file << std::endl;
        }
      }

      if(sum->Empty.size() != 0)
      {
        out << "---- Empty files:" << std::endl;
        for(auto file = sum->Empty.begin(); file != sum->Empty.end(); ++file)
        {
          out << "------ " << *file << std::endl;
        }
      }
    }

    out.close();
  }
}

