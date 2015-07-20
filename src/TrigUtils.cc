#include <iostream>
#include <iomanip>
#include <string>
#include <random>

#include "UserCode/llvv_fwk/interface/TrigUtils.h"

namespace trigUtils
{


bool applyPrescale(float prop) {

  // When runnin on MC only
  
  bool hasPhotonTrigger(false);
  
  std::random_device rd;
  std::mt19937 gen(rd());
    
  std::bernoulli_distribution bern(prop);
  hasPhotonTrigger=bern(gen);

  return hasPhotonTrigger;
}

 bool passPhotonTriggerMC(fwlite::ChainEvent& ev, float &triggerThreshold,
			float &triggerPrescale, bool &prescale ){
   
   edm::TriggerResultsByName tr = ev.triggerResultsByName("HLT");
   if( !tr.isValid() ) return false;

   //prescale=false;
   
    bool hasPhotonTrigger(false);
    triggerPrescale = 1.0; 
    triggerThreshold = 0.0;

    float prop=(1./triggerPrescale);
      
    std::string successfulPath="";
    if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon300_NoHE_*")){
      hasPhotonTrigger=true;
      triggerThreshold=300;
      triggerPrescale=1;

      if (prescale) {
	prop=(1./triggerPrescale);
	hasPhotonTrigger=applyPrescale(prop);
      }
    }
    else if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon250_NoHE_*")){
      hasPhotonTrigger=true;
      triggerThreshold=250;
      triggerPrescale=1;

      if (prescale) {
	prop=(1./triggerPrescale);
	hasPhotonTrigger=applyPrescale(prop);
      }
    }
    else if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon160_*")){
      hasPhotonTrigger=true;
      triggerThreshold=160;
      triggerPrescale=1;
   
      if (prescale) {
	prop=(1./triggerPrescale);
	hasPhotonTrigger=applyPrescale(prop);
      }

    }
    else if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon150_*")){
      hasPhotonTrigger=true;
      triggerThreshold=150;
      triggerPrescale=1;

      if (prescale) {
	prop=(1./triggerPrescale);
	hasPhotonTrigger=applyPrescale(prop);
      }
    }
    else if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon135_*")){
      hasPhotonTrigger=true;
      triggerThreshold=135;
      triggerPrescale=1;

      if (prescale) {
	prop=(1./triggerPrescale);
	hasPhotonTrigger=applyPrescale(prop);
      }
    }
    else if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_*")){
      hasPhotonTrigger=true;
      triggerThreshold=120;
      triggerPrescale=1;

      if (prescale) {
	prop=(1./triggerPrescale);
	hasPhotonTrigger=applyPrescale(prop);
      }
    }
    else if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_*")){
      hasPhotonTrigger=true;
      triggerThreshold=90;
      triggerPrescale=3;

      if (prescale) {
	prop=(1./triggerPrescale);
	hasPhotonTrigger=applyPrescale(prop);
      }
    }
    else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_*")){
      hasPhotonTrigger=true;
      triggerThreshold=75;
      triggerPrescale=6;

      if (prescale) {
	prop=(1./triggerPrescale);
	hasPhotonTrigger=applyPrescale(prop);
      }
    }
    else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_*")){
      hasPhotonTrigger=true;
      triggerThreshold=50;
      triggerPrescale=17;

      if (prescale) {
	prop=(1./triggerPrescale);
	hasPhotonTrigger=applyPrescale(prop);
      }
    }
    else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_*")){
      hasPhotonTrigger=true;
      triggerThreshold=36;
      triggerPrescale=700;

      if (prescale) {
	prop=(1./triggerPrescale);
	hasPhotonTrigger=applyPrescale(prop);
      }
    }
    else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_*")){
      hasPhotonTrigger=true;
      triggerThreshold=22;
      triggerPrescale=4000;

      if (prescale) {
	prop=(1./triggerPrescale);
	hasPhotonTrigger=applyPrescale(prop);
      }
    }

    return hasPhotonTrigger; 
 }

// Plot selPhoton pT spectra for each Photon trigger of the Control region
  void photonControlSample(fwlite::ChainEvent& iEvent, pat::Photon& photon,
			    SmartSelectionMonitor& mon, TString tag
			    ) {
  
  // pat::PhotonCollection selPhotons;

  bool hasPhotonTrigger(false);
  float triggerThreshold(0);
  //float triggerPrescale(1.0); 

  double weight = 1.0;

  double pt=photon.pt();
  double eta=photon.superCluster()->eta();

  edm::TriggerResultsByName tr = iEvent.triggerResultsByName("HLT");
  if( !tr.isValid() ) return;

  std::string successfulPath="";
  if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon300_NoHE_*")){
    hasPhotonTrigger=true;
    triggerThreshold=300;
    // triggerPrescale=1;
  }
  else if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon250_NoHE_*")){
    hasPhotonTrigger=true;
    triggerThreshold=250;
    // triggerPrescale=1;
  }
  
  else if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon160_*")){
    hasPhotonTrigger=true;
    triggerThreshold=160;
  }
  else if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon150_*")){
    hasPhotonTrigger=true;
    triggerThreshold=150;
  }
  else if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon135_*")){
    hasPhotonTrigger=true;
    triggerThreshold=135;
  }
  
  else if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_*")){
    hasPhotonTrigger=true;
    triggerThreshold=120;
    // triggerPrescale=1;
  }
  else if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_*")){
    hasPhotonTrigger=true;
    triggerThreshold=92;
    //  triggerPrescale=2;
  }
  else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_*")){
    hasPhotonTrigger=true;
    triggerThreshold=77;
    // triggerPrescale=3;
  }
  else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_*")){
    hasPhotonTrigger=true;
    triggerThreshold=50;
    // triggerPrescale=6;
  }
  else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_*")){
    hasPhotonTrigger=true;
    triggerThreshold=36;
    // triggerPrescale=700;
  }
  else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_*")){
    hasPhotonTrigger=true;
    triggerThreshold=22;
    // triggerPrescale=4000;
  }
  
  
  if (hasPhotonTrigger && successfulPath!="") {
    //    weight=(1./triggerPrescale);
    // std::cout << "Trigger path is: " << tag+"_"+successfulPath << std::endl;
    mon.fillHisto("phopt", tag+"_"+successfulPath, pt, weight);
    mon.fillHisto("phoeta", tag+"_"+successfulPath, eta, weight);
    if (pt>triggerThreshold) {
      mon.fillHisto("phopt", tag+"_allTrg", pt, weight);
      mon.fillHisto("phoeta", tag+"_allTrg", eta, weight);
    }
  }
}

  void photonControlEff(fwlite::ChainEvent& iEvent, pat::Photon& photon, 
			    SmartSelectionMonitor& mon, TString tag
			    ) {
  

  double weight = 1.0;

  double pt=photon.pt();

  edm::TriggerResultsByName tr = iEvent.triggerResultsByName("HLT");
  if( !tr.isValid() ) return;

  std::string successfulPath="";
  if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon300_NoHE_*")){
    mon.fillHisto("phopt", tag+"_eff_"+successfulPath, pt, weight);
  }
  if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon250_NoHE_*")){
    mon.fillHisto("phopt", tag+"_eff_"+successfulPath, pt, weight);
  }
  if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon160_*")){
    mon.fillHisto("phopt", tag+"_eff_"+successfulPath, pt, weight);
  }
  if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon150_*")){
    mon.fillHisto("phopt", tag+"_eff_"+successfulPath, pt, weight);
   }
  if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon135_*")){
    mon.fillHisto("phopt", tag+"_eff_"+successfulPath, pt, weight);
  }
  if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_*")){
    mon.fillHisto("phopt", tag+"_eff_"+successfulPath, pt, weight);
  }
  if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_*")){
    mon.fillHisto("phopt", tag+"_eff_"+successfulPath, pt, weight);
  }
  if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_*")){
    mon.fillHisto("phopt", tag+"_eff_"+successfulPath, pt, weight);
  }
  if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_*")){
    mon.fillHisto("phopt", tag+"_eff_"+successfulPath, pt, weight);
  }
  if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_*")){
    mon.fillHisto("phopt", tag+"_eff_"+successfulPath, pt, weight);
  }
  if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_*")){
    mon.fillHisto("phopt", tag+"_eff_"+successfulPath, pt, weight);
  }
 

  }



}
