#ifndef FRWEIGHTS_H
#define FRWEIGHTS_H

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
  
};

#endif

