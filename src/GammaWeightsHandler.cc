#include "UserCode/llvv_fwk/interface/GammaWeightsHandler.h"

using namespace std;

//
GammaWeightsHandler::GammaWeightsHandler(const edm::ParameterSet &runProcess,bool forceAllToData)
{
  //cfg
  bool isMC = runProcess.getParameter<bool>("isMC");
  if(forceAllToData) isMC=false;
  std::vector<std::string> gammaPtWeightsFiles =  runProcess.getParameter<std::vector<std::string> >("weightsFile");  
  if(gammaPtWeightsFiles.size()==0) return;
  TString wgtName("qt");
  TString wgtType( isMC ? "mcfitwgts" : "datafitwgts");
    
  //categories to consider, add more if needed but keep these ones 
  TString cats[]   =  {"eq0jets","eq1jets","eq2jets","geq3jets","vbf","geq1jets","novbf","mjjq100","mjjq092","mjjq083","mjjq066","mjjq049","mjjq033","mjjq016"};
  dilCats_.push_back("ee"); dilCats_.push_back("mumu");
  
  //retrieve from file
  TString gammaPtWeightsFile(gammaPtWeightsFiles[0].c_str());
  gSystem->ExpandPathName(gammaPtWeightsFile);
  TFile *fwgt=TFile::Open(gammaPtWeightsFile);
  if(fwgt)
    {
      cout << "[GammaWeightsHandler] retrieving weights from: " << gammaPtWeightsFile << endl;
      
      std::map<TString, TGraph*> iWgtsH;
      for(size_t ic=0; ic<sizeof(cats)/sizeof(TString); ic++)
	{
	  for(size_t id=0; id<dilCats_.size(); id++)
	    {
	      TString key = dilCats_[id] + cats[ic];

	      //weights
	      TString hname= key + "_" + wgtName + "_" + wgtType;
	      TGraph *h = (TGraph *) fwgt->Get(hname);
	      if(h!=0) wgtsH_[key] = h;
	      
	      //mass shape
	      hname= key+"_zmass"; 
	      TH1 *massh = (TH1 *) fwgt->Get(hname);
	      if(massh!=0) { massh->SetDirectory(0); zmassH_[key]= massh; }
	    }
	}
      fwgt->Close();
    }
  
  if(wgtsH_.size()==0) return; 
  std::cout << "[GammaWeightsHandler] gamma spectrum will be reweighted using distributions found in "  
	    << gammaPtWeightsFiles.size() 
	    << " files" 
	    << std::endl;
}

//
LorentzVector GammaWeightsHandler::getMassiveP4(LorentzVector &gamma,TString evCategoryLabel)
{
  //generate a mass from the line shape (0 if not available)
  float mass(0);
  if(zmassH_.find(evCategoryLabel)!=zmassH_.end())
    {
      if(zmassH_[evCategoryLabel]->Integral())
	while(fabs(mass-91)>15) 
	  mass = zmassH_[evCategoryLabel]->GetRandom();
    }
  return LorentzVector(gamma.px(),gamma.py(),gamma.pz(),sqrt(pow(mass,2)+pow(gamma.energy(),2)));
}

LorentzVector GammaWeightsHandler::getMassiveP4(LorentzVectorF &gamma,TString evCategoryLabel)
{
  //generate a mass from the line shape (0 if not available)
  float mass(0);
  if(zmassH_.find(evCategoryLabel)!=zmassH_.end())
    {
      if(zmassH_[evCategoryLabel]->Integral())
        while(fabs(mass-91)>15)
          mass = zmassH_[evCategoryLabel]->GetRandom();
    }
  return LorentzVector(gamma.px(),gamma.py(),gamma.pz(),sqrt(pow(mass,2)+pow(gamma.energy(),2)));
}

float GammaWeightsHandler::getWeightFor(LorentzVector &gamma, TString evCategoryLabel)
{
  //get the weight (1.0 if not available)
  float weight(1.0);
  if(wgtsH_.find(evCategoryLabel) != wgtsH_.end())
    {
      TGraph *h = wgtsH_[evCategoryLabel];
      weight=h->Eval(gamma.pt());
      if(weight<0) weight=0;
    }
    
  return weight;
}

float GammaWeightsHandler::getWeightFor(LorentzVectorF &gamma, TString evCategoryLabel)
{
  //get the weight (1.0 if not available)
  float weight(1.0);
  if(wgtsH_.find(evCategoryLabel) != wgtsH_.end())
    {
      TGraph *h = wgtsH_[evCategoryLabel];
      weight=h->Eval(gamma.pt());
      if(weight<0) weight=0;
    }
    
  return weight;
}


//
GammaWeightsHandler::~GammaWeightsHandler()
{
}
