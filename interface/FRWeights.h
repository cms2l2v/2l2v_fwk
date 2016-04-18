#ifndef FRWEIGHTS_H
#define FRWEIGHTS_H

#include <vector>
#include <list>
#include <iterator>
#include <algorithm>
#include <map>
#include <iostream>
#include <string>
#include <regex>
#include <string>
#include "TFile.h"
#include "TGraphErrors.h"

class FRWeights
{
    public:
        FRWeights();
        ~FRWeights();
   
        bool init(const std::string& WeightsFileName);
        double getWeight(const std::string& cat ,const std::string& bin, const std::string& var,const std::string& wrt ,const double& pT);

    private:
        TFile* WeightsFile;
        std::map<std::string, TGraphErrors*> FRWeightGraphs; 
};

#endif

