#include "UserCode/llvv_fwk/interface/UEAnalysis.h"

using namespace std;

//
UEAnalysis::UEAnalysis(SmartSelectionMonitor &mon)
  : mon_(&mon)
{
  Double_t ptttaxis[]={0,10.,20.,25.,30.,40.,55.,70.,90.,135.,500.};
  Int_t nptttBins=sizeof(ptttaxis)/sizeof(Double_t)-1;

  mon_->addHistogram( new TH1F("ptttbar",";t#bar{t} transverse momentum [GeV];Events",nptttBins,ptttaxis));
  mon_->addHistogram( new TH2F("phiresponse",";Generated #phi(t#bar{t}) [rad]; #Delta #phi(t#bar{t}) [rad];Events",50,0,3.4,50,0,3.4) );
  mon_->addHistogram( new TH2F("phivspt",";Generated t#bar{t} transverse momentum [GeV]; #Delta #phi(t#bar{t}) [rad];Events",nptttBins,ptttaxis,50,0.,3.4) );
  mon_->addHistogram( new TH2F("ptresponse",";Generated p_{T}(t#bar{t}) [GeV]; #Delta p_{T}(t#bar{t}) [GeV]; Events",nptttBins,ptttaxis,50,0.,50.) );

  mon_->addHistogram( new TH1F("metresponse",";E_{T}^{miss} response;Events",50,0,2));
  mon_->addHistogram( new TH2F("metvsnvtx",";Vertices;E_{T}^{miss} [GeV];Events",20,0,100,50,0,500));


  TString pfmatch[]={"","matched"};
  for(size_t i=0; i<2; i++){
    mon_->addHistogram( new TH1F(pfmatch[i]+"pfdz",";|#Delta z| [cm];Charged candidate",50,0,25) );
    mon_->addHistogram( new TH1F(pfmatch[i]+"pfsigdz",";|#Delta z| significance;Charged candidate",50,0,25) );
    mon_->addHistogram( new TH1F(pfmatch[i]+"pfd0",";|#Delta d_{0}| [cm];Charged candidate",50,0,2) );
    mon_->addHistogram( new TH1F(pfmatch[i]+"pfsigd0",";|#Delta d_{0}| significance;Charged candidate",50,0,25) );
  }

  ueReg_.push_back("");
  ueReg_.push_back("away");
  ueReg_.push_back("toward");
  ueReg_.push_back("transverse");
  TString incdistStr("Charged particles");
  TString distStr("d^{3}N/(d#Delta#eta d#Delta#phi d#Delta p_{t#bar{t}})");
  for(size_t ireg=0; ireg<ueReg_.size(); ireg++)
    {
      mon_->addHistogram( new TH1F("nch"+ueReg_[ireg],";Particles;Events",200,0,200));
      mon_->addHistogram( new TH1F("ngench"+ueReg_[ireg],";Particles;Events",200,0,200));
      mon_->addHistogram( new TH2F("nchvsngench"+ueReg_[ireg],";Particles;Generated particles;Events",200,0,200,200,0,200));
      
      mon_->addHistogram( new TH1F("ptflux"+ueReg_[ireg],";p_{T} flux [GeV];Events",50,0,500));
      mon_->addHistogram( new TH1F("avgptflux"+ueReg_[ireg],";Average p_{T} flux [GeV];Events",25,0,25));
      
      mon_->addHistogram( new TH2F("nchprofpt"+ueReg_[ireg],";t#bar{t} transverse momentum [GeV];Particles;"+distStr,nptttBins,ptttaxis, 200,0.,200.));
      mon_->addHistogram( new TH2F("nchprofavgptflux"+ueReg_[ireg],";Average p_{T} flux [GeV];Particles;"+distStr,25,0.,25., 200,0.,200.));
      mon_->addHistogram( new TH2F("ptfluxprofpt"+ueReg_[ireg],";t#bar{t} transverse momentum [GeV];p_{T} flux [GeV];"+distStr,nptttBins,ptttaxis,50,0.,500.));
      mon_->addHistogram( new TH2F("avgptfluxprofpt"+ueReg_[ireg],";t#bar{t} transverse momentum [GeV];Average p_{T} flux [GeV];"+distStr,nptttBins,ptttaxis,25,0.,25.));
      
      if(ireg==0) {
	mon_->addHistogram( new TH1F("nchlt15",";Particles;Events",200,0,200));
	mon_->addHistogram( new TH1F("nchgt15",";Particles;Events",200,0,200));
	//nvertex profiles
	mon_->addHistogram( new TH2F("rawnchprofnvtx",";Vertices;p_{T} flux [GeV];"+distStr,50,0,50,200,0,200));
	mon_->addHistogram( new TH2F("nchprofnvtx",";Vertices;p_{T} flux [GeV];"+distStr,50,0,50,200,0,200));
	mon_->addHistogram( new TH2F("matchednchprofnvtx",";Vertices;Charged candidate;Events",50,0,50,200,0,200) );
      }

      //phi profiles
      TString phiPF("");
      if(ireg==1) phiPF="0to25";
      if(ireg==2) phiPF="25to70";
      if(ireg==3) phiPF="gt70";
      mon_->addHistogram( new TH2F("nchprofphi"+phiPF,       ";#Delta#phi[^{0}];Particles;"+distStr,       40,-180,180,200,0.,200.));
      mon_->addHistogram( new TH2F("avgptfluxprofphi"+phiPF, ";#Delta#phi[^{0}];Average p_{T} flux [GeV];"+distStr,40,-180,180,25,0.,25.));
      mon_->addHistogram( new TH2F("ptfluxprofphi"+phiPF,";#Delta#phi[^{0}];p_{T} flux [GeV];"+distStr,40,-180,180,100,0,100));
    }

  //soft hadronic activity
  mon_->addHistogram( new TH1F("softleadpt",";Leading extra jet p_{T} [GeV];Events",10,0,250) );
  mon_->addHistogram( new TH1F("softht",";Extra jet H_{T} [GeV];Events",10,0,100) );
  TH1 *hsoft_inc=  mon_->addHistogram( new TH1F("nsoftjetsinc",";Extra jet multiplicity;Events",6,0,6));
  TH1 *hsoft_inc10to20=  mon_->addHistogram( new TH1F("nsoftjetsinc10to20",";Extra jet multiplicity;Events",6,0,6));
  TH1 *hsoft_inc20to30=  mon_->addHistogram( new TH1F("nsoftjetsinc20to30",";Extra jet multiplicity;Events",6,0,6));
  TH1 *hsoft_incgt30=  mon_->addHistogram( new TH1F("nsoftjetsincgt30",";Extra jet multiplicity;Events",6,0,6));
  TH1 *hsoft_out=  mon_->addHistogram( new TH1F("nsoftjetsout",";Extra jet multiplicity;Events",6,0,6));
  TH1 *hsoft_bb =  mon_->addHistogram( new TH1F("nsoftjets",   ";Extra jet multiplicity;Events",6,0,6));
  TH1 *hsoft_ll =  mon_->addHistogram( new TH1F("nsoftjetsll", ";Extra jet multiplicity;Events",6,0,6));
  TH2 *hsoft_inc_prof=(TH2 *)mon_->addHistogram( new TH2F("nsoftjetsincvsdetabb",";#Delta #eta (b,b');Extra jet multiplicity;Events",16,0,8,6,0,6));
  TH2 *hsoft_out_prof=(TH2 *)mon_->addHistogram( new TH2F("nsoftjetsoutvsdetabb",";#Delta #eta (b,b');Extra jet multiplicity;Events",16,0,8,6,0,6));
  TH2 *hsoft_bb_prof =(TH2 *)mon_->addHistogram( new TH2F("nsoftjetsvsdetabb",";#Delta #eta (b,b');Extra jet multiplicity;Events",16,0,8,6,0,6));
  TH2 *hsoft_ll_prof =(TH2 *)mon_->addHistogram( new TH2F("nsoftjetsvsdetall",";#Delta #eta (l,l');Extra jet multiplicity;Events",16,0,8,6,0,6));
  for(int ibin=1; ibin<=hsoft_inc->GetXaxis()->GetNbins(); ibin++){
    TString label("="); if(ibin==hsoft_inc->GetXaxis()->GetNbins()) label="#geq"; label += (ibin-1);
    hsoft_inc->GetXaxis()->SetBinLabel(ibin,label);
    hsoft_inc10to20->GetXaxis()->SetBinLabel(ibin,label);
    hsoft_inc20to30->GetXaxis()->SetBinLabel(ibin,label);
    hsoft_incgt30->GetXaxis()->SetBinLabel(ibin,label);
    hsoft_out->GetXaxis()->SetBinLabel(ibin,label);
    hsoft_bb->GetXaxis()->SetBinLabel(ibin,label);
    hsoft_ll->GetXaxis()->SetBinLabel(ibin,label);
    hsoft_inc_prof->GetYaxis()->SetBinLabel(ibin,label);
    hsoft_out_prof->GetYaxis()->SetBinLabel(ibin,label);
    hsoft_bb_prof->GetYaxis()->SetBinLabel(ibin,label);
    hsoft_ll_prof->GetYaxis()->SetBinLabel(ibin,label);
    TString pf(""); pf+= ibin;
    mon_->addHistogram( new TH1F("softpt"+pf,";Extra jet #"+pf+" p_{T} [GeV];Events",10,0,250) );
  }
  for(size_t i=0; i<=2; i++)
    {
      TString nbtags(""); nbtags+=i;
      mon_->addHistogram( new TH1F("softdeltaptrel"+nbtags+"t", ";#Delta_{rel}(j,j');Events",     10,0,1) );
      mon_->addHistogram( new TH1F("softdphijj"+nbtags+"t", ";#Delta#phi_{jj};Events",            10,0,3.5) );
    }

  //soft lepton activity
  hsoft_inc=mon_->addHistogram( new TH1F("nsoftleptons",";Soft leptons multiplicity;Events",6,0,6));
  for(int ibin=1; ibin<=hsoft_inc->GetXaxis()->GetNbins(); ibin++){
    TString label("="); if(ibin==hsoft_inc->GetXaxis()->GetNbins()) label="#geq"; label += (ibin-1);
    hsoft_inc->GetXaxis()->SetBinLabel(ibin,label);
  }
  mon_->addHistogram( new TH1F("softleptonsmll",";Soft dileptons mass [GeV]",20,0,10));


  //summary ntuple
  TString summaryTupleVarNames("ch:weight:normWeight");
  summaryTupleVarNames += ":gen_ptttbar:gen_phittbar:rec_ptttbar:rec_phittbar:rec_nsoftjets";
  summaryTupleVarNames += ":gen_nch_away:rec_nch_away:gen_nch_tow:rec_nch_tow:gen_nch_tran:rec_nch_tran";
  summaryTupleVarNames += ":gen_ptflux_away:rec_ptflux_away:gen_ptflux_tow:rec_ptflux_tow:gen_ptflux_tran:rec_ptflux_tran";
  summaryTupleVarNames += ":gen_avgptflux_away:rec_avgptflux_away:gen_avgptflux_tow:rec_avgptflux_tow:gen_avgptflux_tran:rec_avgptflux_tran";
  summaryTupleVarNames += ":nvtx:njets";
  summaryTupleVarNames += ":leadpt:trailerpt:st:sumpt";
  summaryTuple_ = new TNtuple("ue","ue",summaryTupleVarNames);
  summaryTuple_->SetDirectory(0);
  summaryTupleVars_ = new Float_t[summaryTupleVarNames.Tokenize(":")->GetEntriesFast()];
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
  int chIdx(abs(lid1)*abs(lid2));
  if     (chIdx==11*11 || chIdx==13*13) ch.push_back("ll"); 
  else if(chIdx==11*13) ch.push_back("emu");
  else return;
  ch.push_back("");

  //ttbar system reconstructed
  LorentzVector htlep=leptons[0]+leptons[1]+jets[0]+jets[1];
  LorentzVector rec_ttbar=htlep+met;
  {
    float val=rec_ttbar.pt();
    if(isnan(val))
      {
	cout << "l1:" << leptons[0] << endl;
	cout << "l2:" << leptons[1] << endl;
	cout << "j1:" << jets[0] << endl;
	cout << "j2:" << jets[1] << endl;
	cout << "met:" << met << endl;
	cout << "---->ttbar:" << rec_ttbar << endl;
      }
  }

  float detaBB=fabs(jets[0].eta()-jets[1].eta());
  float detaLL=fabs(leptons[0].eta()-leptons[1].eta());

  //add category depending on the number of extra jets
  float softHt(0);
  int nExtraJets(0), nExtraJetsInBB(0), nExtraJetsInLL(0);
  int nExtraJets10to20(0),nExtraJets20to30(0),nExtraJetsgt30(0);
  data::PhysicsObject_t *softj1=0, *softj2=0;
  for(size_t ijet=2; ijet<jets.size(); ijet++){

    if(fabs(jets[ijet].eta())>2.5 ) continue;

    if(jets[ijet].pt()>15 && jets[ijet].pt()<20) nExtraJets10to20++;
    if(jets[ijet].pt()>=20 && jets[ijet].pt()<30) nExtraJets20to30++;  
    if(jets[ijet].pt()>=30) nExtraJetsgt30++;

    if(jets[ijet].pt()<20) continue;
    nExtraJets++;
    softHt+=jets[ijet].pt();
    if(softj1==0)      softj1=&(jets[ijet]);
    else{
      if(jets[ijet].pt()>softj1->pt())      { softj2=softj1;        softj1=&(jets[ijet]); }
      else if(softj2==0)                    { softj2=&(jets[ijet]);                       }
      else if(jets[ijet].pt()>softj2->pt()) { softj2=&(jets[ijet]);                       }
    }
    if(jets[ijet].eta() > min(jets[0].eta(),jets[1].eta())       && jets[ijet].eta() < max(jets[0].eta(),jets[1].eta()) ) nExtraJetsInBB++;
    if(jets[ijet].eta() > min(leptons[0].eta(),leptons[1].eta()) && jets[ijet].eta() < max(leptons[0].eta(),leptons[1].eta()) ) nExtraJetsInLL++;
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
  LorentzVector genMet(0,0,0,0);
  for(size_t igen=0; igen<gen.size(); igen++)
    {
      if(gen[igen].get("status")!=3) continue;
      if(gen[igen].get("id")==6)  top=gen[igen];
      if(gen[igen].get("id")==-6) antitop=gen[igen];
      if(gen[igen].get("id")==5)  bquark=gen[igen];
      if(gen[igen].get("id")==-5) antibquark=gen[igen];
      if(gen[igen].get("id")==11||gen[igen].get("id")==13)   chLepton=gen[igen];
      if(gen[igen].get("id")==-11||gen[igen].get("id")==-13) antiChLepton=gen[igen];
      if(fabs(gen[igen].get("id"))==12||fabs(gen[igen].get("id"))==14||fabs(gen[igen].get("id"))==16) genMet+=gen[igen];
    }
  LorentzVector gen_ttbar=top+antitop;
  
  float const_rec_ttbar_pt( rec_ttbar.pt()>500 ? 500. : rec_ttbar.pt() );
  Int_t ptbin(mon_->getHisto("ptttbar","emu")->GetXaxis()->FindBin(const_rec_ttbar_pt) );
  float ptbinWidth(mon_->getHisto("ptttbar","emu")->GetXaxis()->GetBinWidth(ptbin));
  mon_->fillHisto("ptttbar",  ch, const_rec_ttbar_pt, weight/ptbinWidth);
  if(genMet.pt()>0) mon_->fillHisto("metresponse", ch, met.pt()/genMet.pt(), weight);
  mon_->fillHisto("metvsnvtx",   ch, nvtx, met.pt(), weight);

  if(gen_ttbar.pt()>0)
    {
      mon_->fillHisto("phiresponse",  ch, fabs(gen_ttbar.phi()), fabs(deltaPhi(rec_ttbar.phi(),gen_ttbar.phi())), weight);
      mon_->fillHisto("phivspt",      ch, gen_ttbar.pt(),        fabs(deltaPhi(rec_ttbar.phi(),gen_ttbar.phi())), weight/ptbinWidth);
      mon_->fillHisto("ptresponse",   ch, gen_ttbar.pt(),        fabs(rec_ttbar.pt()-gen_ttbar.pt()),             weight/ptbinWidth);
    }

  std::vector<int>   genChCount(4,0);
  std::vector<float> genChFlux(4,0);
  for(size_t igen=0; igen<gen.size(); igen++)
    {
      if(gen[igen].get("status") !=1 || gen[igen].get("charge")==0) continue;

      //do not consider if matching the leptons or the charged pf candidates associated to the b-tagged jets
      if( fabs(deltaR(gen[igen],chLepton))<0.1 )     continue; 
      if( fabs(deltaPhi(gen[igen],antiChLepton))<0.1 ) continue; 
      bool belongsToTagJet(false);
      for(size_t ijet=0; (ijet<2 && !belongsToTagJet); ijet++)
	{
	  size_t pfstart=jets[ijet].get("pfstart");
	  size_t pfend=jets[ijet].get("pfend");
	  for(size_t ipfn=pfstart; ipfn<=pfend; ipfn++)
	    {
	      if( fabs(deltaR(gen[igen],pf[ipfn]))>0.1 ) continue;
	      belongsToTagJet=true;
	      break;
	    }
	}
      if(belongsToTagJet) continue;
      
      //check if in acceptance
      if(gen[igen].pt()<minPFpt || fabs(gen[igen].eta())>maxPFeta) continue;
 
      //count this particle
      float dphi=fabs(deltaPhi(gen[igen].phi(),gen_ttbar.phi())*180/TMath::Pi());
      size_t regIdx=3;
      if(dphi>120) regIdx=1;
      if(dphi<60)  regIdx=2;
      genChCount[0]++;                  genChCount[regIdx]++;
      genChFlux[0] += gen[igen].pt();   genChFlux[regIdx] += gen[igen].pt();
    }
  
  
  //
  //RECONSTRUCTED LEVEL ANALYSIS: study UE with charged PF
  //
  std::vector<int> softMuonsIdx;
  int rawChCount(0),matchedChCount(0);
  float rawChFlux(0),matchedChFlux(0);
  std::vector<int>   chCount(4,0);
  std::vector<float> chFlux(4,0);
  const TH1 *nchprofphiH=mon_->getHisto("nchprofphi",ch[0]);
  int nphibins=nchprofphiH->GetXaxis()->GetNbins();
  std::vector<int> chCountPhi(nphibins,0);
  std::vector<float> chFluxPhi(nphibins,0);
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

      //fiducial cut
      if(pf[ipfn].pt()<minPFpt || fabs(pf[ipfn].eta())>maxPFeta) continue;
      
      //remove if matching a selected lepton
      double minDRpfl(9999.);
      for(size_t ilep=0; ilep<2; ilep++)
	minDRpfl = TMath::Min( minDRpfl, deltaR(pf[ipfn],leptons[ilep]) );
      if(minDRpfl<0.05) continue;
     
      float dz(pf[ipfn].getVal("dz")), sdz(pf[ipfn].getVal("sdz"));
      float d0(pf[ipfn].getVal("d0")), sd0(pf[ipfn].getVal("sd0"));
      std::vector<TString> matchCats(1,"");
      for(size_t igen=0; igen<gen.size(); igen++)
	{
	  if(gen[igen].get("status") !=1 || gen[igen].get("charge")==0) continue;
	  float dR=deltaR(gen[igen],pf[ipfn]);
	  if(dR>0.1) continue;
	  matchCats.push_back("matched");
	  pf[ipfn].set("gmatch",igen);
	  matchedChCount++;
	  matchedChFlux+=pf[ipfn].pt();
	  break;
	}
      for(size_t imatch=0; imatch<matchCats.size(); imatch++)
	{
	  mon_->fillHisto(matchCats[imatch]+"pfdz",    ch, fabs(dz),  weight);
	  mon_->fillHisto(matchCats[imatch]+"pfsigdz", ch, fabs(sdz), weight);
	  mon_->fillHisto(matchCats[imatch]+"pfd0",    ch, fabs(d0),  weight);
	  mon_->fillHisto(matchCats[imatch]+"pfsigd0", ch, fabs(sd0), weight);
	}

      rawChCount++;
      rawChFlux += pf[ipfn].pt(); 
      //if(fabs(dz)>1 || fabs(sdz)>3 || fabs(d0)>1 || fabs(sd0)>10) continue;
      if(fabs(sdz)>10 || fabs(sd0)>10) continue;
      
      if(fabs(pf[ipfn].get("id"))==13 && pf[ipfn].pt()>3) softMuonsIdx.push_back( ipfn );

      //do the counting respectively to the ttbar estimate
      float sigdphi=deltaPhi(pf[ipfn].phi(),rec_ttbar.phi())*180/TMath::Pi();
      float dphi=fabs(sigdphi);
      size_t regIdx=3;
      if(dphi>120) regIdx=1;
      if(dphi<60)  regIdx=2;
      chCount[0]++;                  chCount[regIdx]++;                
      chFlux[0] += pf[ipfn].pt();    chFlux[regIdx] += pf[ipfn].pt();  
      int iphibin=nchprofphiH->GetXaxis()->FindBin(sigdphi)-1;
      if(iphibin>=0 && iphibin<nphibins) { chCountPhi[iphibin]++; chFluxPhi[iphibin]+= pf[ipfn].pt(); }
    }

  int totalaway(0); 
  for(int ibin=nchprofphiH->GetXaxis()->FindBin(120); ibin<nphibins; ibin++) totalaway += chCountPhi[ibin];

  //n vertex profiles
  mon_->fillHisto("rawnchprofnvtx",    ch, nvtx, rawChCount,  weight/acceptance);	
  mon_->fillHisto("nchprofnvtx",       ch, nvtx, chCount[0],  weight/acceptance);	
  mon_->fillHisto("matchednchprofnvtx",ch, nvtx, matchedChCount,  weight/acceptance);	
  
  //phi profiles
  float phiProfAcceptance=acceptance*2*TMath::Pi()*nchprofphiH->GetXaxis()->GetBinWidth(1);
  std::vector<TString> phiPFcats(1,"");
  if(const_rec_ttbar_pt<25) phiPFcats.push_back("0to25");
  else if(const_rec_ttbar_pt<70) phiPFcats.push_back("25to70");
  else phiPFcats.push_back("gt70");
  for(size_t iphibin=0; iphibin<chCountPhi.size(); iphibin++)
    {
      float dphi=nchprofphiH->GetXaxis()->GetBinCenter(iphibin);
      for(size_t iphipf=0; iphipf<phiPFcats.size(); iphipf++){
	mon_->fillHisto("nchprofphi"+phiPFcats[iphipf], ch, dphi, chCountPhi[iphibin],  weight/phiProfAcceptance);
	if(chCountPhi[iphibin]>0){
	  mon_->fillHisto("ptfluxprofphi"+phiPFcats[iphipf],    ch, dphi, chFluxPhi[iphibin],  weight/phiProfAcceptance);
	  mon_->fillHisto("avgptfluxprofphi"+phiPFcats[iphipf], ch, dphi, chFluxPhi[iphibin]/chCountPhi[iphibin],  weight/phiProfAcceptance);
	}
      }
    }
  
  //fill profiles
  for(size_t ireg=0; ireg<4; ireg++)
    {
      float cts( chCount[ireg] ), gencts( genChCount[ireg] );
      float flux( chFlux[ireg] ), genflux( genChFlux[ireg] );
      float normFlux(cts>0?flux/cts:0), genNormFlux( gencts>0 ? genflux/gencts : 0);
      
      if(ireg==0){
	TString pf("lt15");
	if(nvtx>15) pf="gt15";
	mon_->fillHisto("nch"+pf,              ch, cts,                      weight/acceptance);
      }
      mon_->fillHisto("nch"+ueReg_[ireg],              ch, cts,                      weight/acceptance);
      mon_->fillHisto("ngench"+ueReg_[ireg],           ch, gencts,                   weight/acceptance);
      mon_->fillHisto("nchvsngench"+ueReg_[ireg],      ch, cts,    gencts,           weight/acceptance);

      mon_->fillHisto("ptflux"+ueReg_[ireg],          ch, flux ,                    weight/acceptance);
      if(cts>0) mon_->fillHisto("avgptflux"+ueReg_[ireg],       ch, normFlux,       weight/acceptance);
      
      mon_->fillHisto("nchprofpt"+ueReg_[ireg],       ch, const_rec_ttbar_pt, cts,      weight/(acceptance*ptbinWidth));
      mon_->fillHisto("ptfluxprofpt"+ueReg_[ireg],    ch, const_rec_ttbar_pt, flux,     weight/(acceptance*ptbinWidth));

      if(cts>0) {
	mon_->fillHisto("avgptfluxprofpt"+ueReg_[ireg], ch, const_rec_ttbar_pt, normFlux, weight/(acceptance*ptbinWidth));
	mon_->fillHisto("nchprofavgptflux"+ueReg_[ireg], ch, normFlux, cts,      weight/(acceptance*ptbinWidth));
      }
    }

  //
  // SOFT HADRONIC ACTIVITY
  //
  mon_->fillHisto("nsoftjetsinc",      ch,nExtraJets,        weight);
  mon_->fillHisto("nsoftjetsinc10to20",      ch,nExtraJets10to20,        weight);
  mon_->fillHisto("nsoftjetsinc20to30",      ch,nExtraJets20to30,        weight);
  mon_->fillHisto("nsoftjetsincgt30",      ch,nExtraJetsgt30,        weight);
  mon_->fillHisto("nsoftjetsout",      ch,nExtraJets-nExtraJetsInBB,        weight);
  mon_->fillHisto("nsoftjets",         ch,nExtraJetsInBB,    weight);
  mon_->fillHisto("nsoftjetsll",       ch,nExtraJetsInLL,    weight);
  mon_->fillHisto("nsoftjetsincvsdetabb", ch,detaBB, nExtraJets,        weight);
  mon_->fillHisto("nsoftjetsoutvsdetabb", ch,detaBB, nExtraJets-nExtraJetsInBB,        weight);
  mon_->fillHisto("nsoftjetsvsdetabb", ch,detaBB, nExtraJetsInBB, weight);
  mon_->fillHisto("nsoftjetsvsdetall", ch,detaLL, nExtraJetsInLL, weight);

  for(size_t ijet=2; ijet<jets.size(); ijet++){
    if(jets[ijet].pt()<15 || fabs(jets[ijet].eta())>2.5 ) continue;
    TString pf(""); pf+= nExtraJets;
    mon_->fillHisto("softpt"+pf,ch,jets[ijet].pt(),weight);
  }

  if(softj1)
    {
      mon_->fillHisto("softleadpt", ch, softj1->pt(), weight);
      mon_->fillHisto("softht",  ch, softHt, weight);

      if(softj1 && softj2)
	{
	  LorentzVector softjj=(*softj1) + (*softj2);
	  float dphijj=deltaPhi(softj1->phi(),softj2->phi());
	  float deltarel=softjj.pt()/(softj1->pt()+softj2->pt());
	  int softBtags( (softj1->getVal("csv")>0.405)+(softj2->getVal("csv")>0.405) );
	  if(softBtags>2) softBtags=2;
	  TString nbtags(""); nbtags += softBtags;
	  
	  mon_->fillHisto("softdeltaptrel"+nbtags+"t", ch, deltarel, weight);
	  mon_->fillHisto("softdphijj"+nbtags+"t", ch, fabs(dphijj), weight);
	}
    }

  //SOFT LEPTONIC ACTIVITY
  mon_->fillHisto("nsoftleptons",      ch,softMuonsIdx.size(),        weight);
  for(size_t isoft=0; isoft<softMuonsIdx.size(); isoft++)
    for(size_t jsoft=isoft+1; jsoft<softMuonsIdx.size(); jsoft++)
      {
	int llCharge( pf[isoft].get("charge") *  pf[jsoft].get("charge") );
	if(llCharge>=0) continue;
	LorentzVector ll=pf[isoft]+pf[jsoft];
	mon_->fillHisto("softleptonsmll", ch, ll.mass(), weight);
      }

  //summary
  summaryTupleVars_[ 0]=chIdx;
  summaryTupleVars_[ 1]=weight;
  summaryTupleVars_[ 2]=1;              //x-sec normalization will be set by fillSummaryTuple
  summaryTupleVars_[ 3]=gen_ttbar.pt();
  summaryTupleVars_[ 4]=gen_ttbar.phi();
  summaryTupleVars_[ 5]=rec_ttbar.pt();
  summaryTupleVars_[ 6]=rec_ttbar.phi();
  summaryTupleVars_[ 7]=nExtraJets;
  summaryTupleVars_[ 8]=genChCount[1];
  summaryTupleVars_[ 9]=chCount[1];
  summaryTupleVars_[10]=genChCount[2];
  summaryTupleVars_[11]=chCount[2];
  summaryTupleVars_[12]=genChCount[3];
  summaryTupleVars_[13]=chCount[3];
  summaryTupleVars_[14]=genChFlux[1];
  summaryTupleVars_[15]=chFlux[1];
  summaryTupleVars_[16]=genChFlux[2];
  summaryTupleVars_[17]=chFlux[2];
  summaryTupleVars_[18]=genChFlux[3];
  summaryTupleVars_[19]=chFlux[3];
  summaryTupleVars_[20]=genChCount[1]>0 ?  genChCount[1]/genChFlux[1] : 0;
  summaryTupleVars_[21]=chCount[1]>0    ?  chCount[1]/chFlux[1]       : 0;
  summaryTupleVars_[22]=genChCount[2]>0 ?  genChCount[2]/genChFlux[2] : 0;
  summaryTupleVars_[23]=chCount[2]>0    ?  chCount[2]/chFlux[2]       : 0;
  summaryTupleVars_[24]=genChCount[3]>0 ?  genChCount[3]/genChFlux[3] : 0;
  summaryTupleVars_[25]=chCount[3]>0    ?  chCount[3]/chFlux[3]       : 0;
  summaryTupleVars_[26]=nvtx;
  summaryTupleVars_[27]=nExtraJets;
  summaryTupleVars_[28]=leptons[0].pt();
  summaryTupleVars_[29]=leptons[1].pt();
  summaryTupleVars_[30]=leptons[0].pt()+leptons[1].pt();
  LorentzVector dil=leptons[0]+leptons[1];
  summaryTupleVars_[31]=dil.pt();
}
