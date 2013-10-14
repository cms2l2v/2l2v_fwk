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
      mon_->addHistogram( new TH1F("mlj"+vars[i],       ";Lepton-jet invariant Mass [GeV];Lepton-jet pairs",nMlj,mljAxis) );
      mon_->addHistogram( new TH1F("correctmlj"+vars[i],";Lepton-jet invariant Mass [GeV];Lepton-jet pairs",nMlj,mljAxis) );
    }
  mon_->addHistogram( new TH1F("rotmlj",  ";Lepton-jet invariant Mass [GeV];Lepton-jet pairs",nMlj,mljAxis) );
  mon_->addHistogram( new TH1F("wrongmlj",";Lepton-jet invariant Mass [GeV];Lepton-jet pairs",nMlj,mljAxis) );

  //flavor analysis
  mon_->addHistogram( new TH2F ("csv1vscsv2", "; Jet #1 CSV; Jet #2 CSV; Jets", 50, -0.2,1.2, 50, -0.2, 1.2) );
  TH2F *h2=(TH2F *)mon_->addHistogram( new TH2F("extrajetflavor",";Flavor;Event type", 4,0,4,9,0,9) );
  h2->GetXaxis()->SetBinLabel(1,"udsg");
  h2->GetXaxis()->SetBinLabel(2,"c");
  h2->GetXaxis()->SetBinLabel(3,"b");
  h2->GetXaxis()->SetBinLabel(4,"Total");
  mon_->addHistogram( new TH2F("extrajetpt",";p_{T} [GeV];Event type", 50,0,500,9,0,9) );
  mon_->addHistogram( new TH2F("matchedjetpt",";p_{T} [GeV];Event type", 50,0,500,9,0,9) );
  for(size_t i=0; i<5; i++){
    TString ijet("");
    if(i) ijet+=i;
    mon_->addHistogram( new TH1F("csv"+ijet,";CSV discriminator;Jets", 50, -0.2, 1.2) );
  }
  TH1 *hl=mon_->addHistogram( new TH1F("csvLbtagsextended",";b-tag multiplicity;Events", 3*3*5, 0.,3*3*5.) );
  TH1 *hm=mon_->addHistogram( new TH1F("csvMbtagsextended",";b-tag multiplicity;Events", 3*3*5, 0.,3*3*5.) );
  TH1 *ht=mon_->addHistogram( new TH1F("csvTbtagsextended",";b-tag multiplicity;Events", 3*3*5, 0.,3*3*5.) );
  for(int ibin=1; ibin<=hl->GetXaxis()->GetNbins(); ibin++)
    {
      TString label(""); label += (ibin-1)%5;
      hl->GetXaxis()->SetBinLabel(ibin,label);
      hm->GetXaxis()->SetBinLabel(ibin,label);
      ht->GetXaxis()->SetBinLabel(ibin,label);
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
  cats.push_back("all");

  int ncorrectAssignments(0);
  std::vector<int> extraJetFlavors;
  std::vector<float> extraJetPt,matchedJetPt;
  std::vector<bool> extraJetHasCSVL, extraJetHasCSVM,extraJetHasCSVT;
  std::vector<bool> matchedJetHasCSVL, matchedJetHasCSVM,matchedJetHasCSVT;
  int nCSVL(0), nCSVM(0), nCSVT(0);
  for(size_t ijet=0; ijet<jets.size(); ijet++){

    //parton match
    const data::PhysicsObject_t &genParton=jets[ijet].getObject("gen");
    int genPartonId=genParton.info.find("id")->second;
    
    //flavor match
    const data::PhysicsObject_t &genJet=jets[ijet].getObject("genJet");
    int flavId=genJet.info.find("id")->second;

    bool hasCSVL(jets[ijet].getVal("csv")>0.405);
    bool hasCSVM(jets[ijet].getVal("csv")>0.783);
    bool hasCSVT(jets[ijet].getVal("csv")>0.920);
    bool correctAssignmentFound(false);
    for(size_t ilep=0; ilep<2; ilep++){

      //lepton match
      const data::PhysicsObject_t &genLep=leptons[ilep].getObject("gen");
      int genLeptonId=genLep.info.find("id")->second;

      //check if assignment is correct
      int assignCode=(genLeptonId*genPartonId); 
      bool isCorrect(assignCode<0 && isTopMC && fabs(flavId)==5 );
      correctAssignmentFound |= isCorrect;

      //fill the histograms (don't weight by bin width to fit with unweighted events in RooFit)
      LorentzVector lj=jets[ijet]+leptons[ilep];
      mon_->fillHisto("mlj"+var,                                  cats,lj.mass(),weight);
      mon_->fillHisto((isCorrect ? "correctmlj" : "wrongmlj")+var,cats,lj.mass(),weight);

      lj=jets[ijet]+rotLeptons_[ilep];
      mon_->fillHisto("rotmlj"+var,                               cats,lj.mass(),weight);
    }

    if(ncorrectAssignments==2) correctAssignmentFound=false; //this should never happen
    ncorrectAssignments += correctAssignmentFound;

    if(!correctAssignmentFound){
      int jetFlavBin(0);
      if(abs(flavId)==4) jetFlavBin=1;
      if(abs(flavId)==5) jetFlavBin=2;
      extraJetFlavors.push_back(jetFlavBin);
      extraJetPt.push_back( jets[ijet].pt() );
      extraJetHasCSVL.push_back(hasCSVL);
      extraJetHasCSVM.push_back(hasCSVM);
      extraJetHasCSVT.push_back(hasCSVT);
    }else{
      matchedJetPt.push_back( jets[ijet].pt() );
      matchedJetHasCSVL.push_back(hasCSVL);
      matchedJetHasCSVM.push_back(hasCSVM);
      matchedJetHasCSVT.push_back(hasCSVT);
    }
    
    nCSVL += hasCSVL;
    nCSVM += hasCSVM;
    nCSVT += hasCSVT;
  }
  
  if(var!="") return;
  
  //b-tag correlation (leading b-tagged jets)
  mon_->fillHisto("csv1vscsv2",cats,jets[0].getVal("csv"),jets[1].getVal("csv"),weight);
  for(size_t ijet=0; ijet<jets.size(); ijet++){
    TString pf(""); pf += (ijet+1);
    mon_->fillHisto("csv",cats,jets[ijet].getVal("csv"),weight);
    mon_->fillHisto("csv"+pf,cats,jets[ijet].getVal("csv"),weight);
  }

  //pt and flavor of matched and extra jets in the event
  int addBin(0);
  if(jets.size()==3) addBin=3;
  if(jets.size()==4) addBin=6;
  mon_->fillHisto("extrajetflavor",     cats, 3,                        ncorrectAssignments+addBin,weight);
  for(size_t iej=0; iej<extraJetFlavors.size(); iej++)
    {
      mon_->fillHisto("extrajetflavor", cats, extraJetFlavors[iej], ncorrectAssignments+addBin,weight);		      
      mon_->fillHisto("extrajetpt",     cats, extraJetPt[iej],      ncorrectAssignments+addBin,weight);		      
      if(extraJetHasCSVL[iej]) mon_->fillHisto("extrajetcsvLpt",     cats, extraJetPt[iej],      ncorrectAssignments+addBin,weight);		      
      if(extraJetHasCSVM[iej]) mon_->fillHisto("extrajetcsvMpt",     cats, extraJetPt[iej],      ncorrectAssignments+addBin,weight);		      
      if(extraJetHasCSVT[iej]) mon_->fillHisto("extrajetcsvTpt",     cats, extraJetPt[iej],      ncorrectAssignments+addBin,weight);		      
    }
  for(size_t imj=0; imj<matchedJetPt.size(); imj++)
    {
      mon_->fillHisto("matchedjetpt", cats, matchedJetPt[imj], ncorrectAssignments+addBin,weight);		      
      if(matchedJetHasCSVL[imj]) mon_->fillHisto("matchedjetpt",     cats, matchedJetPt[imj],      ncorrectAssignments+addBin,weight);		      
      if(matchedJetHasCSVM[imj]) mon_->fillHisto("matchedjetpt",     cats, matchedJetPt[imj],      ncorrectAssignments+addBin,weight);		      
      if(matchedJetHasCSVT[imj]) mon_->fillHisto("matchedjetpt",     cats, matchedJetPt[imj],      ncorrectAssignments+addBin,weight);		      
    }

  //b-tag counting: used for extracting R
  addBin=0;
  if(jets.size()==3) addBin += 5;
  if(jets.size()==4) addBin += 10;
  if(ch=="mumu")   addBin += 15;
  if(ch=="emu")    addBin += 2*15;
  mon_->fillHisto("csvLbtagsextended",cats,nCSVL+addBin,weight);
  mon_->fillHisto("csvMbtagsextended",cats,nCSVM+addBin,weight);
  mon_->fillHisto("csvTbtagsextended",cats,nCSVT+addBin,weight);
}




