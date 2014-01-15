#include "UserCode/llvv_fwk/interface/TopPolarizationAnalysis.h"

using namespace std;

//
LeptonPlusJetsEvent::LeptonPlusJetsEvent( data::PhysicsObjectCollection_t &selJets, data::PhysicsObjectCollection_t &selLeptons, data::PhysicsObject_t &met)
{
  bool passPreSelection(selJets.size()>=4 && selLeptons.size()==1);
  addStepStatus("presel",passPreSelection);
  if(!passPreSelection) return;

  recb1_=selJets[0];
  recb2_=selJets[1];
  recq1_=selJets[2];
  recq2_=selJets[3];
  recl_=selLeptons[0];
  recnu_=met;
  bool passBtagSelection(recb1_.getVal("csv")>0.405 && recb1_.getVal("csv")>0.405);
  addStepStatus("fit",passBtagSelection);
}

//
void LeptonPlusJetsEvent::fit()
{
  if(!pass("presel")) return;
  addStepStatus("fit",false);
}

//
TopPolarizationAnalysis::TopPolarizationAnalysis(SmartSelectionMonitor &mon,std::vector<TString> & vars)
  : mon_(&mon), vars_(vars)
{
  selSteps_.push_back("gen");selSteps_.push_back("fid");selSteps_.push_back("trig");selSteps_.push_back("presel");selSteps_.push_back("btag");selSteps_.push_back("fit");
  for(size_t istep=0; istep<selSteps_.size(); istep++) mon_->addHistogram(new TH2F(selSteps_[istep]+"acc",";Top pseudo-rapidity; Top azimuthal angle [rad]",5,0,5,5,0,3.2) );
}


//
void TopPolarizationAnalysis::analyze(data::PhysicsObjectCollection_t &selJets, data::PhysicsObjectCollection_t &selLeptons, data::PhysicsObject_t &met, data::PhysicsObjectCollection_t &gen)
{
  //fit event kinematics
  LeptonPlusJetsEvent ljEvent(selJets,selLeptons,met);
  ljEvent.fit();

  //
  // consider top had and top lep per event
  // check if objects from top had or top lep are in fiducial range
  // check if objects from top had or top lep are reconstructed and selected
  // classify as:
  //                    | Top Had   | Top Lep
  //   in fiducial      |   0x    1 | 0x    1 0000 0000
  //   reco             |   0x   10 | 0x   10 0000 0000
  //   presel           |   0x  100 | 0x  100 0000 0000
  //   b-tag            |   0x 1000 | 0x 1000 0000 0000
  //   matching fit leg |   0x10000 | 0x10000 0000 0000
  //   2*3*2=20 categories to understand acceptance steps
  //   last category is the one more interesting in the end: 4 event categories
  //
  
  //
  //GENERATOR LEVEL TRUTH	
  //
  //data::PhysicsObject_t &genb1=ljEvent.recb1_.get("genMatch");
  //data::PhysicsObject_t &genb2=ljEvent.recb2_.get("genMatch");
  //data::PhysicsObject_t &genq1=ljEvent.recq1_.get("genMatch");
  // data::PhysicsObject_t &genq2=ljEvent.recq2_.get("genMatch");
  //LorentzVector genHadCand1(genb1+genq1+genq2),  genHadCand2(genb2+genq1+genq2);
  //double        genHadMass1(genHadCand1.mass()), genHadMass2(genHadCand2.mass());
  LorentzVector gentop(0,0,0,0), genantitop(0,0,0,0);
  for(size_t igen=0; igen<gen.size(); igen++)
    {
      int pdgId=gen[igen].get("id");
      if(abs(pdgId)!=6) continue;
      if(pdgId==6) gentop=gen[igen];
      else         genantitop=gen[igen];
    }
  LorentzVector ttbar=gentop+genantitop;
  TString ch("l"); 
  std::vector<TString> genMttbarCat(1,ch);
  if(ttbar.mass()<400)       genMttbarCat.push_back(ch+"0to400");
  else if(ttbar.mass()<1000) genMttbarCat.push_back(ch+"400to1000");	
  else genMttbarCat.push_back(ch+"gt1000");
	
  //acceptance matrices
  for(size_t istep=0; istep<selSteps_.size(); istep++)
    {
      if(!ljEvent.pass(selSteps_[istep])) continue;
      mon_->fillHisto(selSteps_[istep]+"acc",genMttbarCat,fabs(gentop.eta()),fabs(gentop.phi()));
      mon_->fillHisto(selSteps_[istep]+"acc",genMttbarCat,fabs(genantitop.eta()),fabs(genantitop.phi()));
    }
}
