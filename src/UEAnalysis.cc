#include "UserCode/llvv_fwk/interface/UEAnalysis.h"

using namespace std;

//
UEAnalysis::UEAnalysis(SmartSelectionMonitor &mon)
  : mon_(&mon)
{
  mon_->addHistogram( new TH2F("thrustphiresponse",";Generated #phi(t#bar{t}) [rad]; #Delta #phi(t#bar{t}) [rad];Events",50,0,3.4,50,0,3.4) );
  mon_->addHistogram( new TH2F("thrustptresponse",";Generated p_{T}(t#bar{t}) [GeV]; #Delta p_{T}(t#bar{t})/p_{T} (t#bar{t}); Events",60,0,300,50,-50,50) );
  
  
  mon_->addHistogram( new TH1F("ptttbar",";t#bar{t} transverse momentum [GeV];Events",25,0,250));
  mon_->addHistogram( new TH1F("mtttbar",";t#bar{t} transverse mass [GeV];Events",25,0,1000));

  //color flow
  mon_->addHistogram( new TH1F("allj1pull",";#Delta#theta_{t} [rad];Events",50,0,3.2) );
  mon_->addHistogram( new TH1F("chj1pull",";#Delta#theta_{t} [rad];Events",50,0,3.2) );
  mon_->addHistogram( new TH1F("allj2pull",";#Delta#theta_{t} [rad];Events",50,02,3.2) );
  mon_->addHistogram( new TH1F("chj2pull",";#Delta#theta_{t} [rad];Events",50,0,3.2) );


  mon_->addHistogram( new TH2F("alljpull",";#Delta#theta_{t} [rad];|#vec{t}|;Jets",50,0,3.2,100,0,0.1) );
  mon_->addHistogram( new TH2F("chjpull",";#Delta#theta_{t} [rad];|#vec{t}|;Jets",50,0,3.2,100,0,0.1) );

  ueReg_.push_back("");
  ueReg_.push_back("away");
  ueReg_.push_back("toward");
  ueReg_.push_back("transverse");
  TString distStr("d^{2}N/d(#Delta#eta d#Delta#phi)");
  for(size_t dir=0; dir<2; dir++)
    {
      TString dirPF(dir==0 ? "" : "bb");
      for(size_t ireg=0; ireg<ueReg_.size(); ireg++)
	{
	  mon_->addHistogram( new TH1F("nch"+dirPF+ueReg_[ireg],";Charged particles;"+distStr,100,0,100));
	  mon_->addHistogram( new TH1F("ptflux"+dirPF+ueReg_[ireg],";Charged p_{T} flux [GeV];"+distStr,50,0,100));
	  mon_->addHistogram( new TH1F("avgptflux"+dirPF+ueReg_[ireg],";Average p_{T} flux [GeV];"+distStr,25,0,5));
	  
	  mon_->addHistogram( new TH2F("nchprofpt"+dirPF+ueReg_[ireg],";t#bar{t} transverse momentum [GeV];Charged particles;"+distStr,25,0,250, 100,0,100));
	  mon_->addHistogram( new TH2F("ptfluxprofpt"+dirPF+ueReg_[ireg],";t#bar{t} transverse momentum [GeV];Charged p_{T} flux [GeV];"+distStr,25,0,250,50,0,100));
	  mon_->addHistogram( new TH2F("avgptfluxprofpt"+dirPF+ueReg_[ireg],";t#bar{t} transverse momentum [GeV];Average p_{T} flux [GeV];"+distStr,25,0,250,20,0,5));
	  
	  mon_->addHistogram( new TH2F("nchprofmt"+dirPF+ueReg_[ireg],";t#bar{t} transverse mass [GeV];Charged particles;"+distStr,25,0,1000,100,0,100));
	  mon_->addHistogram( new TH2F("ptfluxprofmt"+dirPF+ueReg_[ireg],";t#bar{t} transverse mass [GeV];Charged p_{T} flux [GeV];"+distStr,25,0,1000,50,0,100));
	  mon_->addHistogram( new TH2F("avgptfluxprofmt"+dirPF+ueReg_[ireg],";t#bar{t} transverse mass [GeV];Average p_{T} flux [GeV];"+distStr,25,0,1000,25,0,5));
	  
	  if(ireg==0) {
	    mon_->addHistogram( new TH2F("ptfluxprofphi"+dirPF,";#Delta#phi[^{0}];Charged p_{T} flux [GeV];"+distStr,25,0,180,50,0,100));
	    if(dir==0) mon_->addHistogram( new TH2F("ptfluxprofvtx",";Vertices;Charged p_{T} flux [GeV];"+distStr,50,0,50,50,0,100));
	  }
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

  //check the category
  int lid1(leptons[0].get("id")), lid2(leptons[1].get("id"));
  std::vector<TString> ch;
  if     (abs(lid1)*abs(lid2)==11*11 || abs(lid1)*abs(lid2)==13*13) ch.push_back("ll");
  else if(abs(lid1)*abs(lid2)==11*13) ch.push_back("emu");
  else return;

  //select the tag jets
  if(jets[0].pt()<30 || jets[1].pt()<30) return; 
  if(fabs(jets[0].eta())>2.5 || fabs(jets[1].eta())>2.5) return; 
  if(jets[0].getVal("supercsv")<0.531 || jets[1].getVal("supercsv")<0.531) return;

  //add category depending on the number of extra jets
  int nExtraJets(0);
  for(size_t ijet=2; ijet<jets.size(); ijet++){
    if(jets[ijet].pt()<20 || fabs(jets[ijet].eta())>2.5 ) continue;
    nExtraJets++;
  }
  if(nExtraJets==0) ch.push_back( ch[0]+"eq0j" );
  if(nExtraJets==1) ch.push_back( ch[0]+"eq1j" );
  if(nExtraJets>1)  ch.push_back( ch[0]+"geq2j" );

  //control the resolution
  LorentzVector top,antitop;
  for(size_t igen=0; igen<gen.size(); igen++)
    {
      if(gen[igen].get("id")==6)  top=gen[igen];
      if(gen[igen].get("id")==-6) antitop=gen[igen];
    }
  LorentzVector htlep=leptons[0]+leptons[1]+jets[0]+jets[1];
  LorentzVector rec_ttbar=htlep+met;
  LorentzVector gen_ttbar=top+antitop;
  mon_->fillHisto("ptttbar",  ch, rec_ttbar.pt(), weight);
  mon_->fillHisto("mtttbar",  ch, rec_ttbar.Mt(), weight);
  if(gen_ttbar.pt()>0)
    {
      mon_->fillHisto("thrustphiresponse",  ch, fabs(gen_ttbar.phi()), fabs(deltaPhi(rec_ttbar.phi(),gen_ttbar.phi())), weight);
      mon_->fillHisto("thrustptresponse",   ch, gen_ttbar.pt(),        rec_ttbar.pt()-gen_ttbar.pt(),                   weight);
    }

  //alternative direction defined by the bb system
  LorentzVector rec_bb=jets[0]+jets[1];
 
  //
  //study UE with charged PF
  //
  float acceptance( 2*TMath::Pi() * maxPFeta );
  std::vector<int>   chCount(4,0), chCountBB(4,0);
  std::vector<float> chFlux(4,0),  chFluxBB(4,0);
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

      //do the counting respectively to the bbbar estimate
      dphi=deltaPhi(pf[ipfn].phi(),rec_bb.phi())*180/TMath::Pi();
      regIdx=3;
      if(dphi>120) regIdx=1;
      if(dphi<60)  regIdx=2;
      chCountBB[0]++;                  chCountBB[regIdx]++;
      chFluxBB[0] += pf[ipfn].pt();    chFluxBB[regIdx] += pf[ipfn].pt();
      mon_->fillHisto("ptfluxprofphibb",    ch, dphi, pf[ipfn].pt(),  weight/acceptance);      
    }

  //fill profiles
  for(size_t dir=0; dir<2; dir++)
    {
      TString dirPF(dir==0 ? "" : "bb");
      for(size_t ireg=0; ireg<4; ireg++)
	{
	  float cts( dir==0 ? chCount[ireg] : chCountBB[ireg] );
	  float flux( dir==0 ? chFlux[ireg] :  chFluxBB[ireg] );
	  float normFlux(cts>0?flux/cts:0);
	  
	  mon_->fillHisto("nch"+dirPF+ueReg_[ireg],             ch, cts,                      weight/acceptance);
	  mon_->fillHisto("ptflux"+dirPF+ueReg_[ireg],          ch, flux ,                    weight/acceptance);
	  mon_->fillHisto("avgptflux"+dirPF+ueReg_[ireg],       ch, normFlux,                 weight/acceptance);
	  
	  mon_->fillHisto("nchprofpt"+dirPF+ueReg_[ireg],       ch, rec_ttbar.pt(), cts,      weight/acceptance);
	  mon_->fillHisto("ptfluxprofpt"+dirPF+ueReg_[ireg],    ch, rec_ttbar.pt(), flux,     weight/acceptance);
	  mon_->fillHisto("avgptfluxprofpt"+dirPF+ueReg_[ireg], ch, rec_ttbar.pt(), normFlux, weight/acceptance);
	  mon_->fillHisto("nchprofmt"+dirPF+ueReg_[ireg],       ch, rec_ttbar.Mt(), cts,      weight/acceptance);
	  mon_->fillHisto("ptfluxprofmt"+dirPF+ueReg_[ireg],    ch, rec_ttbar.Mt(), flux,     weight/acceptance);
	  mon_->fillHisto("avgptfluxprofmt"+dirPF+ueReg_[ireg], ch, rec_ttbar.Mt(), normFlux, weight/acceptance);
	}
    }

  //
  //color flow studies
  //
  std::vector<float> pullSummary=utils::cmssw::pullDeltaTheta(jets[0].pt()>jets[1].pt()?jets[0]:jets[1],
							      jets[0].pt()>jets[1].pt()?jets[1]:jets[0],
							      pf);
  if(pullSummary.size()==6)
    {
      mon_->fillHisto("allj1pull", ch, fabs(pullSummary[utils::cmssw::toJ1_ALL]), weight);
      mon_->fillHisto("allj2pull",  ch, fabs(pullSummary[utils::cmssw::toJ2_ALL]), weight);
      mon_->fillHisto("chj1pull", ch, fabs(pullSummary[utils::cmssw::toJ1_CH]), weight);
      mon_->fillHisto("chj2pull",  ch, fabs(pullSummary[utils::cmssw::toJ2_CH]), weight);
    }
  for(size_t ijet=0; ijet<2; ijet++)
    {
      const data::PhysicsObject_t &t_all=jets[ijet].getObject("t_all");
      const data::PhysicsObject_t &t_ch =jets[ijet].getObject("t_ch");
      mon_->fillHisto("alljpull", ch, fabs(t_all.py()), t_all.energy(), weight);
      mon_->fillHisto("chjpull",  ch, fabs(t_all.py()), t_ch.energy(),  weight);
    }
}
