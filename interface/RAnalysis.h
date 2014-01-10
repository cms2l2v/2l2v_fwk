#ifndef _ranalysis_
#define _ranalysis_

#include "UserCode/llvv_fwk/interface/BtagUncertaintyComputer.h"
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/llvv_fwk/interface/DataEventSummaryHandler.h"

#include "TGraphErrors.h"
#include "TRandom2.h"
#include "TString.h"

#include <map>

class RAnalysis
{
  
 public:
  RAnalysis(SmartSelectionMonitor &mon,std::vector<TString> &vars);
  void setTaggingCorrectors(std::map<std::pair<TString,TString>, std::pair<TGraphErrors *,TGraphErrors *> > &btagEffCorr) { btagEffCorr_=&btagEffCorr; }
  void prepareAnalysis(data::PhysicsObjectCollection_t &leptons, data::PhysicsObjectCollection_t &jets);
  void analyze(data::PhysicsObjectCollection_t &leptons, data::PhysicsObjectCollection_t &jets, float weight, TString var, bool isTopMC,TString ctrl="");

  ~RAnalysis(){};
    
 private:
  
  inline Int_t getKinCategory(Float_t j1pt, Float_t j2pt=-9999)
    {
      Int_t defaultCat(-1);
      for(size_t i=0; i<btvKinCatsHi_.size(); i++)
	{
	  if(j1pt<btvKinCatsLo_[i].first || j1pt>=btvKinCatsHi_[i].first) continue;
	  if(j2pt<0) return i;
	  if(j2pt>=btvKinCatsLo_[i].second && j2pt<btvKinCatsHi_[i].second) return i;
	}
      return defaultCat;
    }

  
  SmartSelectionMonitor *mon_;
  data::PhysicsObjectCollection_t rotLeptons_;
  TRandom2 rndGen_;
  BTagSFUtil btsfutil_;
  std::map<std::pair<TString,TString>, std::pair<TGraphErrors *,TGraphErrors *> > *btagEffCorr_;

  std::vector<std::pair<Float_t,Float_t > > btvKinCatsHi_,btvKinCatsLo_;
  std::vector<std::pair<TString,Float_t> > diffTaggers_;
};

#endif
