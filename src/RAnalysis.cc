#include "UserCode/llvv_fwk/interface/RAnalysis.h"
#include <Math/VectorUtil.h>

using namespace std;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >::BetaVector BetaVector;

//
RAnalysis::RAnalysis(SmartSelectionMonitor &mon,std::vector<TString> & vars)
  : mon_(&mon)
{

  //mlj spectrum
  Float_t mljAxis[]={0,   10, 20, 30, 40, 50, 60, 70,  80,  90,
		     100,110,120,130,140,150,160, 170, 180, 190,
		     200,210,220,230,240,250,260, 270, 280, 290,
		     300,310,320,330,340,350,360, 370, 380, 390,
		     400,450,500,550,600,700,800, 1000};
  const size_t nMlj=sizeof(mljAxis)/sizeof(Float_t)-1;
  for(size_t i=0; i<vars.size(); i++)
    {
      TH1F *h=new TH1F("mlj"+vars[i],";Lepton-jet invariant Mass [GeV];Lepton-jet pairs",nMlj,mljAxis);
      mon_->addHistogram(h);
    }
}

//
void RAnalysis::prepareAnalysis(data::PhysicsObjectCollection_t &leptons, data::PhysicsObjectCollection_t &jets)
{
  rotLeptons_.clear();
  for(size_t ilep=0; ilep<2; ilep++)
    {
      int itry=0;
      do
	{
	  
	  data::PhysicsObject_t rotLepton(leptons[ilep]);
	  itry++;
	  if(itry>1000) { cout << "Failed to rotate lepton:" << itry << endl;  rotLeptons_.push_back( rotLepton ); break;  }
	  
	  //rotate lepton
	  double en    = rotLepton.E();
	  double pabs  = rotLepton.P();
	  double phi   = rndGen_.Uniform(0,2*TMath::Pi());
	  double theta = TMath::ACos( rndGen_.Uniform(-1,1) );
	  rotLepton.SetPxPyPzE(pabs*TMath::Cos(phi)*TMath::Sin(theta),pabs*TMath::Sin(phi)*TMath::Sin(theta),pabs*TMath::Cos(theta),en);
	  
	  //require selectable kinematics
	  if( TMath::Abs(rotLepton.Eta())>2.4 || rotLepton.Pt()<20 ) continue;	  
	  
	  //require object separation wrt to jets
	  double minDR(1000);
	  for(data::PhysicsObjectCollection_t::iterator jit = jets.begin(); jit != jets.end(); jit++)
	    {
	      double dR = deltaR(jit->eta(),jit->phi(),rotLepton.eta(),rotLepton.phi());
	      if(dR>minDR) continue;
	      minDR=dR;
	    }
	  if(minDR<0.4) continue;
	  
	  //compatible with object selection
	  rotLeptons_.push_back(rotLepton);
	  break;
	} while( 1 );
    }
}


//
void RAnalysis::analyze(data::PhysicsObjectCollection_t &leptons, data::PhysicsObjectCollection_t &jets, float weight, TString var, bool isTopMC)
{

  //check channel
  TString ch("");
  int lid1(leptons[0].get("id")), lid2(leptons[1].get("id"));
  if     (abs(lid1)*abs(lid2)==11*11)  ch="ee";
  else if(abs(lid1)*abs(lid2)==11*13)  ch="emu";
  else if(abs(lid1)*abs(lid2)==13*13)  ch="mumu";
  if(ch=="") return;

  //build categories
  std::vector<TString>    cats(1,ch);
  if(jets.size()==2)      cats.push_back(ch+"eq2jets");
  else if(jets.size()==3) cats.push_back(ch+"eq3jets");
  else if(jets.size()==4) cats.push_back(ch+"eq4jets");
  else return;

  for(size_t ijet=0; ijet<jets.size(); ijet++){

    //parton match
    const data::PhysicsObject_t &genParton=jets[ijet].getObject("gen");
    int genPartonId=genParton.info.find("id")->second;
    
    //flavor match
    const data::PhysicsObject_t &genJet=jets[ijet].getObject("genJet");
    int flavId=genJet.info.find("id")->second;

    for(size_t ilep=0; ilep<2; ilep++){

      //lepton match
      const data::PhysicsObject_t &genLep=leptons[ilep].getObject("gen");
      int genLeptonId=genLep.info.find("id")->second;

      //check if assignment is correct
      int assignCode=(genLeptonId*genPartonId); 
      bool isCorrect(assignCode<0 && isTopMC && fabs(flavId)==5 );

      //fill the histograms
      LorentzVector lj=jets[ijet]+leptons[ilep];
      mon_->fillHisto("mlj"+var,                                  cats,lj.mass(),weight,true);
      mon_->fillHisto((isCorrect ? "correctmlj" : "wrongmlj")+var,cats,lj.mass(),weight,true);

      lj=jets[ijet]+rotLeptons_[ilep];
      mon_->fillHisto("rotmlj"+var,                               cats,lj.mass(),weight,true);
    }
  }

}




