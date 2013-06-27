#include "UserCode/llvv_fwk/interface/UEAnalysis.h"

using namespace std;

//
UEAnalysis::UEAnalysis(SmartSelectionMonitor &mon)
  : mon_(&mon)
{
  TString recType[]={"raw","t1","t12"};
  for(size_t i=0; i<3; i++){
    mon_->addHistogram( new TH2F("thrustphiresponse_"+recType[i],";Generated #phi(t#bar{t}) [rad]; #Delta #phi(t#bar{t}) [rad];Events",50,0,3.4,50,0,3.4) );
    mon_->addHistogram( new TH2F("thrustptresponse_"+recType[i],";Generated p_{T}(t#bar{t}) [GeV]; #Delta p_{T}(t#bar{t})/p_{T} (t#bar{t}); Events",60,0,300,50,-50,50) );
  }
  
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
  for(size_t ireg=0; ireg<ueReg_.size(); ireg++)
    {
      mon_->addHistogram( new TH1F("nch"+ueReg_[ireg],";Charged particles;Events",100,0,100));
      mon_->addHistogram( new TH1F("ptflux"+ueReg_[ireg],";Charged p_{T} flux [GeV];Events",50,0,100));
      mon_->addHistogram( new TH1F("avgptflux"+ueReg_[ireg],";Average p_{T} flux [GeV];Events",25,0,5));

      mon_->addHistogram( new TH2F("nchprofpt"+ueReg_[ireg],";t#bar{t} transverse momentum [GeV];Charged particles;Events",25,0,250, 100,0,100));
      mon_->addHistogram( new TH2F("ptfluxprofpt"+ueReg_[ireg],";t#bar{t} transverse momentum [GeV];Charged p_{T} flux [GeV];Events",25,0,250,50,0,100));
      mon_->addHistogram( new TH2F("avgptfluxprofpt"+ueReg_[ireg],";t#bar{t} transverse momentum [GeV];Average p_{T} flux [GeV];Events",25,0,250,20,0,5));

      mon_->addHistogram( new TH2F("nchprofmt"+ueReg_[ireg],";t#bar{t} transverse mass [GeV];Charged particles;Events",25,0,1000,100,0,100));
      mon_->addHistogram( new TH2F("ptfluxprofmt"+ueReg_[ireg],";t#bar{t} transverse mass [GeV];Charged p_{T} flux [GeV];Events",25,0,1000,50,0,100));
      mon_->addHistogram( new TH2F("avgptfluxprofmt"+ueReg_[ireg],";t#bar{t} transverse mass [GeV];Average p_{T} flux [GeV];Events",25,0,1000,25,0,5));

      if(ireg==0)
	mon_->addHistogram( new TH2F("ptfluxprofphi",";#Delta#phi[^{0}];Charged p_{T} flux [GeV];Events",25,0,180,50,0,100));
    }
}

//
void UEAnalysis::analyze(data::PhysicsObjectCollection_t &leptons, 
			 data::PhysicsObjectCollection_t &jets,
			 data::PhysicsObjectCollection_t &met, 
			 data::PhysicsObjectCollection_t &pf,
			 data::PhysicsObjectCollection_t &gen,
			 float weight)
{
  float minPFpt(0.5), maxPFeta(2.1); 

  //check the category
  int lid1(leptons[0].get("id")), lid2(leptons[1].get("id"));
  std::vector<TString> ch;
  if     (abs(lid1)*abs(lid2)==11*11 || abs(lid1)*abs(lid2)==13*13) ch.push_back("ll");
  else if(abs(lid1)*abs(lid2)==11*13) ch.push_back("emu");
  else return;

  //control the resolution
  LorentzVector top,antitop;
  for(size_t igen=0; igen<gen.size(); igen++)
    {
      if(gen[igen].get("id")==6)  top=gen[igen];
      if(gen[igen].get("id")==-6) antitop=gen[igen];
    }
  LorentzVector htlep=leptons[0]+leptons[1]+jets[0]+jets[1];
  LorentzVector gen_ttbar=top+antitop;
  if(gen_ttbar.pt()>0)
    {
      mon_->fillHisto("thrustphiresponse_raw",  ch, fabs(gen_ttbar.phi()), fabs(deltaPhi((htlep+met[0]).phi(),gen_ttbar.phi())), weight);
      mon_->fillHisto("thrustphiresponse_t1",   ch, fabs(gen_ttbar.phi()), fabs(deltaPhi((htlep+met[2]).phi(),gen_ttbar.phi())), weight);
      mon_->fillHisto("thrustphiresponse_t12",  ch, fabs(gen_ttbar.phi()), fabs(deltaPhi((htlep+met[3]).phi(),gen_ttbar.phi())), weight);
      mon_->fillHisto("thrustptresponse_raw",   ch, gen_ttbar.pt(),        (htlep+met[0]).pt()-gen_ttbar.pt(),                   weight);
      mon_->fillHisto("thrustptresponse_t1",    ch, gen_ttbar.pt(),        (htlep+met[2]).pt()-gen_ttbar.pt(),                   weight);
      mon_->fillHisto("thrustptresponse_t12",   ch, gen_ttbar.pt(),        (htlep+met[3]).pt()-gen_ttbar.pt(),                   weight);
    }

  //using raw met as the estimator
  LorentzVector rec_ttbar=htlep+met[0];

  //color flow
  std::vector<float> pullSummary=utils::cmssw::pullDeltaTheta(jets[0].pt()>jets[1].pt()?jets[0]:jets[1],
							      jets[0].pt()>jets[1].pt()?jets[1]:jets[0],
							      pf);
  
  //study UE with charged PF
  std::vector<int> chCount(4,0);
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
      
      float dphi=deltaPhi(pf[ipfn].phi(),rec_ttbar.phi())*180/TMath::Pi();
      size_t regIdx=3;
      if(dphi>120) regIdx=1;
      if(dphi<60)  regIdx=2;
      chCount[0]++;                  chCount[regIdx]++;
      chFlux[0] += pf[ipfn].pt();    chFlux[regIdx] += pf[ipfn].pt();
    
      mon_->fillHisto("ptfluxprofphi",    ch, dphi, pf[ipfn].pt(),  weight);
    }
  mon_->fillHisto("ptttbar",  ch, rec_ttbar.pt(), weight);
  mon_->fillHisto("mtttbar",  ch, rec_ttbar.Mt(), weight);

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

  for(size_t ireg=0; ireg<4; ireg++)
    {
      float cts(chCount[ireg]);
      float flux(chFlux[ireg]);
      float normFlux(cts>0?flux/cts:0);
      
      mon_->fillHisto("nch"+ueReg_[ireg],             ch, cts,                      weight);
      mon_->fillHisto("ptflux"+ueReg_[ireg],          ch, flux ,                    weight);
      mon_->fillHisto("avgptflux"+ueReg_[ireg],       ch, normFlux,                 weight);

      mon_->fillHisto("nchprofpt"+ueReg_[ireg],       ch, rec_ttbar.pt(), cts,      weight);
      mon_->fillHisto("ptfluxprofpt"+ueReg_[ireg],    ch, rec_ttbar.pt(), flux,     weight);
      mon_->fillHisto("avgptfluxprofpt"+ueReg_[ireg], ch, rec_ttbar.pt(), normFlux, weight);
      mon_->fillHisto("nchprofmt"+ueReg_[ireg],       ch, rec_ttbar.Mt(), cts,      weight);
      mon_->fillHisto("ptfluxprofmt"+ueReg_[ireg],    ch, rec_ttbar.Mt(), flux,     weight);
      mon_->fillHisto("avgptfluxprofmt"+ueReg_[ireg], ch, rec_ttbar.Mt(), normFlux, weight);
    }
}
