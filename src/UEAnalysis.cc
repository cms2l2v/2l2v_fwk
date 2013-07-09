#include "UserCode/llvv_fwk/interface/UEAnalysis.h"

using namespace std;

//
UEAnalysis::UEAnalysis(SmartSelectionMonitor &mon)
  : mon_(&mon)
{
  mon_->addHistogram( new TH2F("thrustphiresponse",";Generated #phi(t#bar{t}) [rad]; #Delta #phi(t#bar{t}) [rad];Events",50,0,3.4,50,0,3.4) );
  mon_->addHistogram( new TH2F("thrustptresponse",";Generated p_{T}(t#bar{t}) [GeV]; #Delta p_{T}(t#bar{t})/p_{T} (t#bar{t}); Events",60,0,300,50,-50,50) );
  
  mon_->addHistogram( new TH1F("ptttbar",";t#bar{t} transverse momentum [GeV];Events",25,0,500));

  ueReg_.push_back("");
  ueReg_.push_back("away");
  ueReg_.push_back("toward");
  ueReg_.push_back("transverse");
  TString distStr("d^{2}N/d(#Delta#eta d#Delta#phi)");
  for(size_t ireg=0; ireg<ueReg_.size(); ireg++)
    {
      mon_->addHistogram( new TH1F("nch"+ueReg_[ireg],";Charged particles;"+distStr,200,0,200));
      mon_->addHistogram( new TH1F("ngench"+ueReg_[ireg],";Generated charged particles;"+distStr,200,0,200));
      mon_->addHistogram( new TH2F("nchvsngench"+ueReg_[ireg],";Charged particles;Generated charged particles;"+distStr,200,0,200,200,0,200));
      
      mon_->addHistogram( new TH1F("ptflux"+ueReg_[ireg],";Charged p_{T} flux [GeV];"+distStr,50,0,500));
      mon_->addHistogram( new TH1F("avgptflux"+ueReg_[ireg],";Average p_{T} flux [GeV];"+distStr,25,0,25));
      
      mon_->addHistogram( new TH2F("nchprofpt"+ueReg_[ireg],";t#bar{t} transverse momentum [GeV];Charged particles;"+distStr,20,0,250, 200,0,200));
      mon_->addHistogram( new TH2F("ptfluxprofpt"+ueReg_[ireg],";t#bar{t} transverse momentum [GeV];Charged p_{T} flux [GeV];"+distStr,20,0,250,50,0,500));
      mon_->addHistogram( new TH2F("avgptfluxprofpt"+ueReg_[ireg],";t#bar{t} transverse momentum [GeV];Average p_{T} flux [GeV];"+distStr,20,0,250,25,0,25));
      
      mon_->addHistogram( new TH2F("nchprofmt"+ueReg_[ireg],";t#bar{t} transverse mass [GeV];Charged particles;"+distStr,20,0,1000,200,0,200));
      mon_->addHistogram( new TH2F("ptfluxprofmt"+ueReg_[ireg],";t#bar{t} transverse mass [GeV];Charged p_{T} flux [GeV];"+distStr,20,0,1000,50,0,500));
      mon_->addHistogram( new TH2F("avgptfluxprofmt"+ueReg_[ireg],";t#bar{t} transverse mass [GeV];Average p_{T} flux [GeV];"+distStr,20,0,1000,25,0,25));
      
      if(ireg==0) {
	mon_->addHistogram( new TH2F("ptfluxprofphi",";#Delta#phi[^{0}];Charged p_{T} flux [GeV];"+distStr,20,0,180,100,0,100));
	mon_->addHistogram( new TH2F("ptfluxprofnvtx",";Vertices;Charged p_{T} flux [GeV];"+distStr,50,0,50,100,0,100));
      }
    }
}


//
void UEAnalysis::analyze(data::PhysicsObjectCollection_t &leptons, 
			 data::PhysicsObjectCollection_t &jets,
			 LorentzVector &met,
			 data::PhysicsObjectCollection_t &pf,
			 data::PhysicsObjectCollection_t &gen,
			 int nvtx,
			 float weight)
{
  float minPFpt(0.5), maxPFeta(2.1);
  float acceptance( 2*TMath::Pi() * maxPFeta );

  //check the event category
  int lid1(leptons[0].get("id")), lid2(leptons[1].get("id"));
  std::vector<TString> ch;
  if     (abs(lid1)*abs(lid2)==11*11 || abs(lid1)*abs(lid2)==13*13) ch.push_back("ll"); 
  else if(abs(lid1)*abs(lid2)==11*13) ch.push_back("emu");
  else return;
  ch.push_back("");

  //select the tag jets
  if(jets[0].pt()<30 || jets[1].pt()<30) return; 
  if(fabs(jets[0].eta())>2.5 || fabs(jets[1].eta())>2.5) return; 
  if(jets[0].getVal("supercsv")<0.531 || jets[1].getVal("supercsv")<0.531) return;

  //ttbar system reconstructed
  LorentzVector htlep=leptons[0]+leptons[1]+jets[0]+jets[1];
  LorentzVector rec_ttbar=htlep+met;

  //add category depending on the number of extra jets
  int nExtraJets(0);
  for(size_t ijet=2; ijet<jets.size(); ijet++){
    if(jets[ijet].pt()<20 || fabs(jets[ijet].eta())>2.5 ) continue;
    nExtraJets++;
  }
  if(nExtraJets==0) ch.push_back( ch[0]+"eq0j" );
  if(nExtraJets==1) ch.push_back( ch[0]+"eq1j" );
  if(nExtraJets>1)  ch.push_back( ch[0]+"geq2j" );

  //
  //GENERATOR LEVEL ANALYSIS
  //
  LorentzVector top,antitop;
  LorentzVector chLepton,antiChLepton;
  LorentzVector bquark,antibquark;
  for(size_t igen=0; igen<gen.size(); igen++)
    {
      if(gen[igen].get("status")!=3) continue;
      if(gen[igen].get("id")==6)  top=gen[igen];
      if(gen[igen].get("id")==-6) antitop=gen[igen];
      if(gen[igen].get("id")==5)  bquark=gen[igen];
      if(gen[igen].get("id")==-5) antibquark=gen[igen];
      if(gen[igen].get("id")==11||gen[igen].get("id")==13)   chLepton=gen[igen];
      if(gen[igen].get("id")==-11||gen[igen].get("id")==-13) antiChLepton=gen[igen];
    }
  LorentzVector gen_ttbar=top+antitop;
  mon_->fillHisto("ptttbar",  ch, rec_ttbar.pt(), weight);
  if(gen_ttbar.pt()>0)
    {
      mon_->fillHisto("thrustphiresponse",  ch, fabs(gen_ttbar.phi()), fabs(deltaPhi(rec_ttbar.phi(),gen_ttbar.phi())), weight);
      mon_->fillHisto("thrustptresponse",   ch, gen_ttbar.pt(),        rec_ttbar.pt()-gen_ttbar.pt(),                   weight);
    }

  std::vector<int>   genChCount(4,0);
  std::vector<float> genChFlux(4,0);
  for(size_t igen=0; igen<gen.size(); igen++)
    {
      if(gen[igen].get("status") !=1 || gen[igen].get("charge")==0) continue;

      //do not consider if matching the leptons or the bquarks (within R=0.5 as they'll fragment to jets)
      float pDphi(gen[igen].phi());
      if( fabs(deltaPhi(pDphi,bquark.phi()))<0.5)        continue;     
      if( fabs(deltaPhi(pDphi,antibquark.phi()))<0.5)    continue; 
      if( fabs(deltaPhi(pDphi,chLepton.phi()))<0.05)     continue; 
      if( fabs(deltaPhi(pDphi,antiChLepton.phi()))<0.05) continue; 

      //check if in acceptance
      if(gen[igen].pt()<minPFpt || fabs(gen[igen].eta())>maxPFeta) continue;
 
      //count this particle
      float  dphi=deltaPhi(pDphi,gen_ttbar.phi())*180/TMath::Pi();
      size_t regIdx=3;
      if(dphi>120) regIdx=1;
      if(dphi<60)  regIdx=2;
      genChCount[0]++;                  genChCount[regIdx]++;
      genChFlux[0] += gen[igen].pt();   genChFlux[regIdx] += gen[igen].pt();
    }
  
  
  //
  //RECONSTRUCTED LEVEL ANALYSIS: study UE with charged PF
  //
  std::vector<int>   chCount(4,0);
  std::vector<float> chFlux(4,0);
  for(size_t ipfn=0; ipfn<pf.size(); ipfn++)
    {
      if(pf[ipfn].get("charge")==0) continue;
      
      //remove if it belongs to a tag jet
      bool belongsToTagJet(false);
      for(size_t ijet=0; ijet<2; ijet++)
	{
	  size_t pfstart=jets[ijet].get("pfstart");
	  size_t pfend=jets[ijet].get("pfend");
	  if(ipfn>=pfstart && ipfn<=pfend) belongsToTagJet=true;
	}
      if(belongsToTagJet) continue;
      
      //remove if matching a selected lepton
      double minDRpfl(9999.);
      for(size_t ilep=0; ilep<2; ilep++)
	minDRpfl = TMath::Min( minDRpfl, deltaR(pf[ipfn],leptons[ilep]) );
      if(minDRpfl<0.05) continue;
      
      if(pf[ipfn].pt()<minPFpt || fabs(pf[ipfn].eta())>maxPFeta) continue;

      mon_->fillHisto("ptfluxprofnvtx",   ch, nvtx, pf[ipfn].pt(),  weight/acceptance);

      //do the counting respectively to the ttbar estimate
      float dphi=deltaPhi(pf[ipfn].phi(),rec_ttbar.phi())*180/TMath::Pi();
      size_t regIdx=3;
      if(dphi>120) regIdx=1;
      if(dphi<60)  regIdx=2;
      chCount[0]++;                  chCount[regIdx]++;
      chFlux[0] += pf[ipfn].pt();    chFlux[regIdx] += pf[ipfn].pt();
      mon_->fillHisto("ptfluxprofphi",    ch, dphi, pf[ipfn].pt(),  weight/acceptance);
    }
  
  //fill profiles
  for(size_t ireg=0; ireg<4; ireg++)
    {
      float cts( chCount[ireg] ), gencts( genChCount[ireg] );
      float flux( chFlux[ireg] ), genflux( genChFlux[ireg] );
      float normFlux(cts>0?flux/cts:0), genNormFlux( gencts>0 ? genflux/gencts : 0);
      
      mon_->fillHisto("nch"+ueReg_[ireg],              ch, cts,                      weight/acceptance);
      mon_->fillHisto("ngench"+ueReg_[ireg],           ch, gencts,                   weight/acceptance);
      mon_->fillHisto("nchvsngench"+ueReg_[ireg],      ch, cts,    gencts,           weight/acceptance);

      mon_->fillHisto("ptflux"+ueReg_[ireg],          ch, flux ,                    weight/acceptance);
      if(cts>0) mon_->fillHisto("avgptflux"+ueReg_[ireg],       ch, normFlux,       weight/acceptance);
      
      mon_->fillHisto("nchprofpt"+ueReg_[ireg],       ch, rec_ttbar.pt(), cts,      weight/acceptance);
      mon_->fillHisto("ptfluxprofpt"+ueReg_[ireg],    ch, rec_ttbar.pt(), flux,     weight/acceptance);

      if(cts>0) mon_->fillHisto("avgptfluxprofpt"+ueReg_[ireg], ch, rec_ttbar.pt(), normFlux, weight/acceptance);
    }
}
