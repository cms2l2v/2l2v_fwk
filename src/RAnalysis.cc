#include "UserCode/llvv_fwk/interface/RAnalysis.h"
#include <Math/VectorUtil.h>
#include "UserCode/llvv_fwk/interface/BtagUncertaintyComputer.h"

using namespace std;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >::BetaVector BetaVector;

//
RAnalysis::RAnalysis(SmartSelectionMonitor &mon,std::vector<TString> & vars)
  : mon_(&mon),
    btagEffCorr_(0)
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

  TString ptcats[]={"","csvL","csvM","csvT"};
  for(size_t i=0; i<4; i++){
    if(i==0) mon_->addHistogram( new TH2F("extrajetpt",";p_{T} [GeV];Event type", 50,0,500,9,0,9) );
    else     mon_->addHistogram( (TH2F *) h2->Clone("extrajetflavor"+ptcats[i]) );
    mon_->addHistogram( new TH2F("matchedjetpt"+ptcats[i],";p_{T} [GeV];Event type", 50,0,500,9,0,9) );
    mon_->addHistogram( new TH2F("flavorjetpt"+ptcats[i],";p_{T} [GeV];Event type", 50,0,500,3,0,3) );
  }
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
  mon_->addHistogram( (TH1F *) hl->Clone("csvLbtagsextendedcorr") );
  mon_->addHistogram( (TH1F *) hm->Clone("csvMbtagsextendedcorr") );
  mon_->addHistogram( (TH1F *) ht->Clone("csvTbtagsextendedcorr") );

  //differential b-tag efficiency measurement
  Float_t ptCats[]={30,50,80,120,210,320};
  const size_t nPtCats(sizeof(ptCats)/sizeof(Float_t));
  Int_t totalCats(Int_t(TMath::Binomial(nPtCats,2))+nPtCats);
  mon_->addHistogram( new TH1F("btvkincats1",";Kinematics categories;",totalCats,0,totalCats) );
  mon_->addHistogram( new TH1F("btvkincats2",";Kinematics categories;",totalCats,0,totalCats) );
  diffTaggers_.push_back(std::pair<TString,Float_t>("csv",0.405));
  diffTaggers_.push_back(std::pair<TString,Float_t>("csv",0.594));
  diffTaggers_.push_back(std::pair<TString,Float_t>("csv",0.783));
  diffTaggers_.push_back(std::pair<TString,Float_t>("csv",0.819));
  diffTaggers_.push_back(std::pair<TString,Float_t>("csv",0.855));
  size_t nTaggers(diffTaggers_.size());
  Int_t catCtr(0);
  for(size_t icat=0; icat<nPtCats; icat++)
    {
      for(size_t jcat=0; jcat<=icat; jcat++, catCtr++)
	{
	  mon_->fillHisto("btvkincats1","all",catCtr,ptCats[icat]);
	  mon_->fillHisto("btvkincats2","all",catCtr,ptCats[jcat]);

	  //flavor counting for efficiency versus pT can be taken inclusively for each pT bin
	  if(jcat==0)
	    {
	      TString singleTag("kin"); singleTag += catCtr;
	      TH1F *hflavCts=(TH1F *)mon_->addHistogram( new TH1F("btvflavcounts"+singleTag,";Jet flavour;Jets",4,0.,4.) );
	      hflavCts->GetXaxis()->SetBinLabel(1,"l");
	      hflavCts->GetXaxis()->SetBinLabel(2,"c");
	      hflavCts->GetXaxis()->SetBinLabel(3,"b");
	      hflavCts->GetXaxis()->SetBinLabel(4,"Total");
	      for(size_t itagger=0; itagger<nTaggers; itagger++)
		{
		  TString tagToFill(diffTaggers_[itagger].first); tagToFill += itagger;
		  mon_->addHistogram( (TH1F *)hflavCts->Clone("btvflavcounts"+singleTag+tagToFill) );
		}
	    }

	  //define the kinematics category
	  TString tag("kin"); tag+=catCtr;
	  Float_t j1ptLo(ptCats[icat]);
	  Float_t j1ptHi(icat==nPtCats-1 ?  99999. : ptCats[icat+1]);
	  Float_t j2ptLo(ptCats[jcat]);
	  Float_t j2ptHi(jcat==nPtCats-1 ?  99999. : ptCats[jcat+1]);
	  btvKinCatsLo_.push_back( std::pair<Float_t,Float_t>(j1ptLo,j2ptLo) );
	  btvKinCatsHi_.push_back( std::pair<Float_t,Float_t>(j1ptHi,j2ptHi) );
	  
	  //flavour categories
	  TH1F *hfij=(TH1F *)mon_->addHistogram( new TH1F("btvfij"+tag,";Flavour category;Events",3*3,0,3*3) );
	  hfij->GetXaxis()->SetBinLabel(1,"ll");
	  hfij->GetXaxis()->SetBinLabel(2,"lc");
	  hfij->GetXaxis()->SetBinLabel(3,"lb");
	  hfij->GetXaxis()->SetBinLabel(4,"cl");
	  hfij->GetXaxis()->SetBinLabel(5,"cc");
	  hfij->GetXaxis()->SetBinLabel(6,"cb");
	  hfij->GetXaxis()->SetBinLabel(7,"bl");
	  hfij->GetXaxis()->SetBinLabel(8,"bc");
	  hfij->GetXaxis()->SetBinLabel(9,"bb");
	  
	  //tag counting
	  for(size_t itagger=0; itagger<nTaggers; itagger++)
	    {
	      TString tagToFill(diffTaggers_[itagger].first); tagToFill += itagger;
	      TH1F *htagCts=(TH1F *)mon_->addHistogram( new TH1F("btv"+tag+tagToFill,";b-tag multiplicity;Events", 3*3, 0.,3*3.) );
	      for(int ibin=1; ibin<=htagCts->GetXaxis()->GetNbins(); ibin++)
		{
		  TString label(""); label += (ibin-1)%3;
		  htagCts->GetXaxis()->SetBinLabel(ibin,label);
		}
	    }
	}
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
void RAnalysis::analyze(data::PhysicsObjectCollection_t &leptons, data::PhysicsObjectCollection_t &jets, float weight, TString var, bool isTopMC,TString ctrl)
{
  //check channel
  TString ch("");
  int lid1(leptons[0].get("id")), lid2(leptons[1].get("id"));
  if     (abs(lid1)*abs(lid2)==11*11)  ch="ee";
  else if(abs(lid1)*abs(lid2)==11*13)  ch="emu";
  else if(abs(lid1)*abs(lid2)==13*13)  ch="mumu";
  if(ch=="") return;
  if(ctrl!="") ch =ch+ctrl; 

  //build categories
  std::vector<TString>    cats(1,ch);
  if(jets.size()==2)      cats.push_back(ch+"eq2jets");
  else if(jets.size()==3) cats.push_back(ch+"eq3jets");
  else if(jets.size()==4) cats.push_back(ch+"eq4jets");
  else return;
  cats.push_back("all"+ctrl);
  
  int ncorrectAssignments(0);
  std::vector<int> extraJetFlavors;
  std::vector<float> extraJetPt,matchedJetPt;
  std::vector<bool> extraJetHasCSVL, extraJetHasCSVM,extraJetHasCSVT;
  std::vector<bool> matchedJetHasCSVL, matchedJetHasCSVM,matchedJetHasCSVT;
  int nCSVL(0), nCSVM(0), nCSVT(0);
  int nCSVLcorr(0), nCSVMcorr(0), nCSVTcorr(0);
  float matchedCSV1(-100), matchedCSV2(-100);
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
      if(rotLeptons_.size()>ilep){
	lj=jets[ijet]+rotLeptons_[ilep];
	mon_->fillHisto("rotmlj"+var,                               cats,lj.mass(),weight);
      }
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
      if(matchedJetPt.size()==1) matchedCSV1 = jets[ijet].get("csv");
      else                       matchedCSV2 = jets[ijet].get("csv");
    }
    
    nCSVL += hasCSVL;
    nCSVM += hasCSVM;
    nCSVT += hasCSVT;

    //apply corrections (if available)
    if(btagEffCorr_==0) {
      nCSVLcorr += hasCSVL;
      nCSVMcorr += hasCSVM;
      nCSVTcorr += hasCSVT;
      continue;
    }

    TString flavKey("udsg");
    if(abs(flavId)==4) flavKey="c";
    if(abs(flavId)==5) flavKey="b";
    float jetpt=min(jets[ijet].pt(),400.0);
    for(int itagger=0; itagger<3; itagger++)
      {
	std::pair<TString,TString> key("csvL",flavKey);
	bool hasBtagCorr(hasCSVL);
	if(itagger==1) { key.first="csvM"; hasBtagCorr=hasCSVM; }
	if(itagger==2) { key.first="csvT"; hasBtagCorr=hasCSVT; }
	if(btagEffCorr_->find(key)!=btagEffCorr_->end())
	  {
	    TGraphErrors *mceffGr=(*btagEffCorr_)[key].first;
	    TGraphErrors *sfGr=(*btagEffCorr_)[key].second;
	    if(mceffGr && sfGr){
	      float eff=mceffGr->Eval(jetpt);
	      float sf=sfGr->Eval(jetpt);
	      btsfutil_.modifyBTagsWithSF(hasBtagCorr,sf,eff);
	    }
	  }
	if(itagger==0)     nCSVLcorr += hasBtagCorr;
	if(itagger==1)     nCSVMcorr += hasBtagCorr;
	if(itagger==2)     nCSVTcorr += hasBtagCorr;
      }
  }
  

  if(var!="") return;

  for(size_t ijet=0; ijet<jets.size(); ijet++){

    //expected b-tag efficiency
    const data::PhysicsObject_t &genJet=jets[ijet].getObject("genJet");
    int flavId=genJet.info.find("id")->second;
    int flavBin(0);
    if(abs(flavId)==4) flavBin=1;
    if(abs(flavId)==5) flavBin=2;
    mon_->fillHisto("flavorjetpt",     cats, jets[ijet].pt(),      flavBin,weight);		      
    if(jets[ijet].getVal("csv")>0.405) mon_->fillHisto("flavorjetptcsvL", cats, jets[ijet].pt(), flavBin,weight);		      
    if(jets[ijet].getVal("csv")>0.783) mon_->fillHisto("flavorjetptcsvM", cats, jets[ijet].pt(), flavBin,weight);		      
    if(jets[ijet].getVal("csv")>0.920) mon_->fillHisto("flavorjetptcsvT", cats, jets[ijet].pt(), flavBin,weight);		      
    
    if(ijet>2) continue; 
    //b-tag correlation (leading b-tagged jets)
    mon_->fillHisto("csv",cats,jets[ijet].getVal("csv"),weight);
    TString pf(""); pf += (ijet+1);
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
      if(extraJetHasCSVL[iej]) mon_->fillHisto("extrajetflavorcsvL",     cats, extraJetFlavors[iej],      ncorrectAssignments+addBin,weight);		      
      if(extraJetHasCSVM[iej]) mon_->fillHisto("extrajetflavorcsvM",     cats, extraJetFlavors[iej],      ncorrectAssignments+addBin,weight);		      
      if(extraJetHasCSVT[iej]) mon_->fillHisto("extrajetflavorcsvT",     cats, extraJetFlavors[iej],      ncorrectAssignments+addBin,weight);		      
      mon_->fillHisto("extrajetpt",     cats, extraJetPt[iej],      ncorrectAssignments+addBin,weight);		      
    } 
  if(matchedJetPt.size()==2) mon_->fillHisto("csv1vscsv2",cats,matchedCSV1, matchedCSV2, weight);
  for(size_t imj=0; imj<matchedJetPt.size(); imj++)
    {
      mon_->fillHisto("matchedjetpt", cats, matchedJetPt[imj], ncorrectAssignments+addBin,weight);		      
      if(matchedJetHasCSVL[imj]) mon_->fillHisto("matchedjetptcsvL",     cats, matchedJetPt[imj],      ncorrectAssignments+addBin,weight);		      
      if(matchedJetHasCSVM[imj]) mon_->fillHisto("matchedjetptcsvM",     cats, matchedJetPt[imj],      ncorrectAssignments+addBin,weight);		      
      if(matchedJetHasCSVT[imj]) mon_->fillHisto("matchedjetptcsvT",     cats, matchedJetPt[imj],      ncorrectAssignments+addBin,weight);		      
    }

  //b-tag counting: used for extracting R
  addBin=0;
  if(jets.size()==3) addBin += 5;
  if(jets.size()==4) addBin += 10;
  if(ch.Contains("mumu"))   addBin += 15;
  if(ch.Contains("emu"))    addBin += 2*15;

  //nominal prediction
  mon_->fillHisto("csvLbtagsextended",cats,nCSVL+addBin,weight);
  mon_->fillHisto("csvMbtagsextended",cats,nCSVM+addBin,weight);
  mon_->fillHisto("csvTbtagsextended",cats,nCSVT+addBin,weight);

  //this is just MC with corrected efficiencies
  mon_->fillHisto("csvLbtagsextendedcorr",cats,nCSVLcorr+addBin,weight);
  mon_->fillHisto("csvMbtagsextendedcorr",cats,nCSVMcorr+addBin,weight);
  mon_->fillHisto("csvTbtagsextendedcorr",cats,nCSVTcorr+addBin,weight);


  //differential measurements
  std::vector<Int_t> ptIdx(2,-1);
  std::vector<Float_t> ptLead(2,0);
  for(size_t ijet=0; ijet<jets.size(); ijet++)
    {
      if( jets[ijet].pt()>ptLead[0] )
	{
	  ptIdx[1]=ptIdx[0]; ptLead[1]=ptLead[0];
	  ptIdx[0]=ijet;     ptLead[0]=jets[ijet].pt();
	}
      else if( jets[ijet].pt()>ptLead[1] )
	{
	  ptIdx[1]=ijet;     ptLead[1]=jets[ijet].pt();
	}
    }

  //differential flavours
  std::vector<Int_t> flavIds;
  std::map<TString,Int_t> nDiffTaggersForBTV;
  for(size_t ijet=0; ijet<2; ijet++)
    {
      Int_t jCat( getKinCategory(jets[ ptIdx[ijet] ].pt()) );

      const data::PhysicsObject_t &genJet=jets[ ptIdx[ijet] ].getObject("genJet");
      int flavId=genJet.info.find("id")->second;
      const data::PhysicsObject_t &genParton=jets[ ptIdx[ijet] ].getObject("gen");
      int genPartonId=genParton.info.find("id")->second;
      flavIds.push_back(flavId);

      std::vector<TString> taggersToFill(1,"");
      for(size_t itag=0; itag<diffTaggers_.size(); itag++)
	{
	  TString algo(diffTaggers_[itag].first);
	  Float_t wp(diffTaggers_[itag].second);
	  bool hasTag(jets[ ptIdx[ijet] ].getVal(algo)>wp);
	  TString tagToFill(algo);  tagToFill+= itag;	 
	  if(hasTag) taggersToFill.push_back( tagToFill );
	  if(nDiffTaggersForBTV.find( tagToFill )==nDiffTaggersForBTV.end() ) nDiffTaggersForBTV[ tagToFill ]=0;
	  nDiffTaggersForBTV[ tagToFill ] += hasTag;
	}
      
      TString baseName("btvflavcountskin"); baseName += jCat;
      for(size_t itagger=0; itagger<taggersToFill.size(); itagger++)
	{
	  TString name=baseName+taggersToFill[itagger];
	  mon_->fillHisto(name,cats,3,weight);
	  if(abs(flavId)==5)      mon_->fillHisto(name,cats,2,weight);
	  else if(abs(flavId)==4) mon_->fillHisto(name,cats,1,weight);
	  else                    mon_->fillHisto(name,cats,0,weight);
	}
    }
  
  //dijet category
  Int_t jjCat( getKinCategory(jets[ ptIdx[0] ].pt(),jets[ ptIdx[1] ].pt()) );
  TString jjCatStr("kin"); jjCatStr += jjCat;

  //determine the flavour category
  Int_t fijBin=1;
  if(abs(flavIds[0])==5)
    {
      if(abs(flavIds[1])==5)      fijBin=9;
      else if(abs(flavIds[1])==4) fijBin=8;
      else                        fijBin=7;
    }
  else if(abs(flavIds[0])==4)
    {
      if(abs(flavIds[1])==5)      fijBin=6;
      else if(abs(flavIds[1])==4) fijBin=5;
      else                        fijBin=4;
    }
  else
    {
      if(abs(flavIds[1])==5)      fijBin=3;
      else if(abs(flavIds[1])==4) fijBin=2;
      else                        fijBin=1;
    }
  mon_->fillHisto("btvfij"+jjCatStr,cats,fijBin-1,weight);

  //tag multplicity
  addBin=0;
  if(ch.Contains("mumu")) addBin +=3;
  if(ch.Contains("emu")) addBin +=2*3;
  for(std::map<TString,Int_t>::iterator it=nDiffTaggersForBTV.begin(); it!=nDiffTaggersForBTV.end(); it++)
    mon_->fillHisto("btv"+jjCatStr+it->first,cats,it->second+addBin,weight);
}





