#ifndef _toppolarizationanalysis_h_
#define _toppolarizationanalysis_h_

#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/llvv_fwk/interface/DataEventSummaryHandler.h"

#include <map>

/**
   @class LeptonPlusJetsEvent
   @short run HitFit for a l+jets event
   @author P. Silva
 */
class LeptonPlusJetsEvent
{
 public:
  LeptonPlusJetsEvent(data::PhysicsObjectCollection_t &selJets, data::PhysicsObjectCollection_t &selLeptons, data::PhysicsObject_t &met);
  ~LeptonPlusJetsEvent() { };
  
  void fit();
  bool pass(TString step)  
  { 
    if(status_.find(step)==status_.end()) return true;
    return status_[step];
  }
  void addStepStatus(TString step,bool status)
  {
    status_[step]=status;
  }
  
  data::PhysicsObject_t recb1_,recb2_,recq1_,recq2_,recl_,recnu_;
  
 private:

  std::map<TString,bool> status_;
  LorentzVector fit_top[2], fit_antitop[2];
  Bool_t fitStatus[2];
  Double_t fitChi2[2];
};

/**
   @class TopPolarizationAnalysis
   @short performs the polariztion analysis 
   @author P.Silva
 */

class TopPolarizationAnalysis
{
 public:
  TopPolarizationAnalysis(SmartSelectionMonitor &mon,std::vector<TString> & vars);
  void analyze(data::PhysicsObjectCollection_t &selJets, data::PhysicsObjectCollection_t &selLeptons, data::PhysicsObject_t &met, data::PhysicsObjectCollection_t &gen);
  ~TopPolarizationAnalysis() {};
 private:
  SmartSelectionMonitor *mon_;
  std::vector<TString> selSteps_,vars_;
};

#endif
