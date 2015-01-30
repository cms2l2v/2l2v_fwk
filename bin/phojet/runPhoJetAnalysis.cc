// -*- C++ -*-
//
// Xin Shi <Xin.Shi@cern.ch>
// Fri Dec 5 11:24:56 CET 2014
// 
// Analysis for the photon + jet
// 

#include <iostream>

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/llvv_fwk/interface/PatUtils.h"
#include "UserCode/llvv_fwk/interface/HiggsUtils.h"
#include "UserCode/llvv_fwk/interface/BtagUncertaintyComputer.h"

#include "RecoJets/JetProducers/interface/PileupJetIdAlgo.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

// constants

const int MUON_PDGID = 13;
const int ELECTRON_PDGID = 11;


SmartSelectionMonitor
initHistograms(){
  SmartSelectionMonitor mon;
  // event selection
  TH1F *h=(TH1F*) mon.addHistogram(new TH1F ("eventflow", ";;Events", 9,0,9) );
  h->GetXaxis()->SetBinLabel(1,"raw");
  h->GetXaxis()->SetBinLabel(2,"#geq 2 iso leptons");
  h->GetXaxis()->SetBinLabel(3,"|M-91|<15");
  h->GetXaxis()->SetBinLabel(4,"p_{T}>55");
  h->GetXaxis()->SetBinLabel(5,"3^{rd}-lepton veto");
  h->GetXaxis()->SetBinLabel(6,"b-veto"); 
  h->GetXaxis()->SetBinLabel(7,"#Delta #phi(jet,E_{T}^{miss})>0.5");
  h->GetXaxis()->SetBinLabel(8,"E_{T}^{miss}>80");
  
  // pile up 
  mon.addHistogram(new TH1F("nvtx", ";Vertices;Events", 50, 0, 50) ); 

  // photon 
  mon.addHistogram(new TH1F("rho", ";Average energy density (#rho);Events", 100, 0, 100) ); 
  mon.addHistogram(new TH1F("npho", ";Number of Photons;Events", 20, 0, 20) ); 
  mon.addHistogram(new TH1F("phopt", ";Photon pT [GeV];Events", 100, 0, 1000) ); 
  mon.addHistogram(new TH1F("phoeta", ";Photon pseudo-rapidity;Events", 50, 0, 5) );
  
  // jet 
  mon.addHistogram(new TH1F("njet", ";Number of Jets;Events", 100, 0, 100) );
  mon.addHistogram(new TH1F("jetpt", ";Jet pT [GeV];Events", 100, 0, 1000) ); 
  mon.addHistogram(new TH1F("jeteta", ";Jet pseudo-rapidity;Events", 50, 0, 5) );
  mon.addHistogram(new TH1F("jetrawen", ";Jet raw energy;Events", 100, 0, 1000) );
  mon.addHistogram(new TH1F("jetnhf", ";Jet neutral hadron energy fraction;Events", 100, 0, 1) );
  mon.addHistogram(new TH1F("jetnef", ";Jet electromagnetic energy fraction;Events", 100, 0, 1) );
  mon.addHistogram(new TH1F("jetcef", ";Jet charged electromagnetic energy fraction;Events", 100, 0, 1) );
  mon.addHistogram(new TH1F("jetchf", ";Jet charged hadron energy fraction;Events", 100, 0, 1) );
  mon.addHistogram(new TH1F("jetnch", ";Jet charged multiplicity;Events", 100, 0, 100) );
  mon.addHistogram(new TH1F("jetnconst", ";Jet number of constitutes;Events", 100, 0, 100) );
  mon.addHistogram(new TH1F("jetpudsct", ";Jet pileup ID discriminant;Events", 100, -1, 1) );

  // met
  mon.addHistogram(new TH1F("met", ";Missing ET [GeV];Events", 100, 0, 1000) );

  // lepton
  // mon.addHistogram(new TH1F("leadpt", ";Transverse momentum [GeV];Events", 50,0,500) );
  mon.addHistogram(new TH1F("mindrlg", ";Min #Delta R(lepton, #gamma);Events", 100, 0, 10) );
  mon.addHistogram(new TH1F("nlep", ";Number of leptons;Events", 10, 0, 10) );
  mon.addHistogram(new TH1F("nexlep", ";Number of extra leptons;Events", 10, 0, 10) );
  mon.addHistogram(new TH1F("zmass", ";Mass [GeV];Events", 100, 40, 250) );
  mon.addHistogram(new TH1F("zy", ";Rapidity;Events", 50, 0, 3) );
  mon.addHistogram(new TH1F("zpt", ";Transverse momentum [GeV];Events", 100, 0, 1500));
  mon.addHistogram(new TH1F("qt", ";Transverse momentum [GeV];Events / (1 GeV)",1500,0,1500));
  mon.addHistogram(new TH1F("qtraw", ";Transverse momentum [GeV];Events / (1 GeV)",1500,0,1500));

  //extra leptons in the event
  mon.addHistogram(new TH1F("nextraleptons", ";Extra leptons;Events",4,0,4) );
  mon.addHistogram(new TH1F("thirdleptonpt", ";Transverse momentum;Events", 50,0,500) );
  mon.addHistogram(new TH1F("thirdleptoneta", ";Pseudo-rapidity;Events", 50,0,2.6) );
  mon.addHistogram(new TH1F("thirdleptonmt", ";Transverse mass(3^{rd} lepton,E_{T}^{miss}) [GeV];Events", 50,0,500) );

  mon.addHistogram( new TH1F("csv",      ";Combined Secondary Vertex;Jets",50,0.,1.) );
  mon.addHistogram( new TH1F("csvb",     ";Combined Secondary Vertex;Jets",50,0.,1.) );
  mon.addHistogram( new TH1F("csvc",     ";Combined Secondary Vertex;Jets",50,0.,1.) );
  mon.addHistogram( new TH1F("csvothers",";Combined Secondary Vertex;Jets",50,0.,1.) );
  
  mon.addHistogram(new TH1F("nbtags", ";b-tag multiplicity;Events",5,0,5) );
  

  return mon; 
}


pat::PhotonCollection
passPhotonSelection(SmartSelectionMonitor mon,
		    pat::PhotonCollection photons,
		    float triggerThreshold,
		    double rho){

  pat::PhotonCollection selPhotons;
  TString tag = "all";  
  double weight = 1.0;

  for(size_t ipho=0; ipho<photons.size(); ipho++) {
    pat::Photon photon = photons[ipho]; 
    double pt=photon.pt();
    double eta=photon.superCluster()->eta();
    mon.fillHisto("phopt", tag, pt, weight);
    mon.fillHisto("phoeta", tag, eta, weight);

    if( pt < triggerThreshold || fabs(eta)>1.4442 ) continue;

    bool passId = patUtils::passId(photon, rho, patUtils::llvvPhotonId::Tight);
    if(!passId) continue; 
    selPhotons.push_back(photon);
  }

  return selPhotons; 
}

bool passPFJetID(SmartSelectionMonitor mon,
		 std::string label,
		 pat::Jet jet){
  
  TString tag = "all";  
  double weight = 1.0;
  bool passID(false); 
  
  float rawJetEn(jet.correctedJet("Uncorrected").energy() );

  double eta=jet.eta();
 
  float nhf( (jet.neutralHadronEnergy() + jet.HFHadronEnergy())/rawJetEn );
  float nef( jet.neutralEmEnergy()/rawJetEn );
  float cef( jet.chargedEmEnergy()/rawJetEn );
  float chf( jet.chargedHadronEnergy()/rawJetEn );
  float nch    = jet.chargedMultiplicity();
  float nconst = jet.numberOfDaughters();
  
  mon.fillHisto("jetrawen", tag, rawJetEn, weight);
  mon.fillHisto("jetnhf", tag, nhf, weight);
  mon.fillHisto("jetnef", tag, nef, weight);
  mon.fillHisto("jetcef", tag, cef, weight);
  mon.fillHisto("jetchf", tag, chf, weight);
  mon.fillHisto("jetnch", tag, nch, weight);
  mon.fillHisto("jetnconst", tag, nconst, weight);

   // use the original for now
   if (label == "Loose") 
    passID = (nhf<0.99  && nef<0.99 && nconst>1); 

  if(fabs(eta)<2.4) {
    passID &= (chf>0 && nch>0 && cef<0.99);
  }

  return passID; 
  
}


bool passCutBasedPUJetID(SmartSelectionMonitor mon,
			 pat::Jet jet,
			 edm::ParameterSet pujetidparas,
			 reco::VertexCollection vtx){
  // Currently not for miniAOD yet
  bool passID(false); 

  std::unique_ptr<PileupJetIdAlgo> cutBasedPuJetIdAlgo;
  cutBasedPuJetIdAlgo.reset(new PileupJetIdAlgo(pujetidparas)); 

  // float jec=1./ev.jn_torawsf[ev.jn];
  // PileupJetIdentifier cutBasedPuIdentifier = cutBasedPuJetIdAlgo_->computeIdVariables(dynamic_cast<const reco::Jet*>(jet->originalObject()), jec, primVtx.get(), *vtxH.product(), true);

  // https://github.com/cms-analysis/flashgg/blob/master/MicroAODProducers/plugins/JetProducer.cc#L86
  // PileupJetIdentifier lPUJetId = pileupJetIdAlgo_->computeIdVariables(pjet.get(),vtx,*vertexCandidateMap,true);


  

  return passID; 
}

pat::JetCollection
passJetSelection(SmartSelectionMonitor mon,
		 pat::JetCollection jets,
		 edm::ParameterSet pujetidparas,
		 reco::VertexCollection vtx,
		 std::vector<patUtils::GenericLepton> selLeptons,
		 pat::PhotonCollection selPhotons){
  
  pat::JetCollection selJets;
  TString tag = "all";  
  double weight = 1.0;

  for(size_t ijet=0; ijet<jets.size(); ijet++){
    pat::Jet jet = jets[ijet]; 
    double pt=jet.pt();
    double eta=jet.eta();
    mon.fillHisto("jetpt", tag, pt, weight);
    mon.fillHisto("jeteta", tag, eta, weight);

    //mc truth for this jet
    const reco::GenJet* genJet=jet.genJet();
    TString jetType( genJet && genJet->pt()>0 ? "truejetsid" : "pujetsid" );

    //cross-clean with selected leptons and photons
    float minDRlj(9999.), minDRlg(9999.);
    for(size_t ilep=0; ilep<selLeptons.size(); ilep++)
      minDRlj = TMath::Min( minDRlj, deltaR(jets[ijet],selLeptons[ilep]) );
    for(size_t ipho=0; ipho<selPhotons.size(); ipho++)
      minDRlg = TMath::Min( minDRlg, deltaR(jets[ijet],selPhotons[ipho]) );

    if(minDRlj<0.4 || minDRlg<0.4) continue;
    
    if(pt<15 || fabs(eta)>4.7 ) continue;
    bool passPFloose = passPFJetID(mon, "Loose", jet); 
    if (!passPFloose) continue;  

    // bool passPuId = passCutBasedPUJetID(mon, jet, pujetidparas, vtx);
    // if (!passPuId) continue;
    float jetpudsct = jet.userFloat("pileupJetId:fullDiscriminant");
    mon.fillHisto("jetpudsct", tag, jetpudsct, weight);
    
    selJets.push_back(jet);
    std::sort(selJets.begin(), selJets.end(), utils::sort_CandidatesByPt);

  }
  return selJets; 
}

void passLeptonSelection(SmartSelectionMonitor mon,
			 std::vector<patUtils::GenericLepton> leptons,
			 pat::PhotonCollection selPhotons,
			 reco::VertexCollection vtx, 
			 std::vector<patUtils::GenericLepton> & selLeptons,
			 std::vector<patUtils::GenericLepton> & extraLeptons) {
  
  LorentzVector muDiff(0,0,0,0);
  for(size_t ilep=0; ilep<leptons.size(); ilep++) {
    bool passKin(true),passId(true),passIso(true);
    bool passLooseLepton(true), passSoftMuon(true);
    bool passSoftElectron(true), passVetoElectron(true);

    int lid=leptons[ilep].pdgId();
    
    //no muon corrections yet, which need the charge info

    //no need for charge info in the following process 
    lid=abs(lid);
    
    //veto nearby photon (loose electrons are many times photons...)
    double minDRlg(9999.);
    for(size_t ipho=0; ipho<selPhotons.size(); ipho++)
      minDRlg=TMath::Min(minDRlg,deltaR(leptons[ilep].p4(),selPhotons[ipho].p4()));

    mon.fillHisto("mindrlg", "all", minDRlg, 1.0);
    if(minDRlg<0.1) continue;

    //kinematics
    float leta = fabs(lid==ELECTRON_PDGID ?  leptons[ilep].el.superCluster()->eta() : leptons[ilep].eta());
    if(leta > (lid==ELECTRON_PDGID ? 2.5 : 2.4) ) passKin=false;
    if(lid==ELECTRON_PDGID && (leta > 1.4442 && leta < 1.5660)) passKin=false;

    passLooseLepton &= passKin;
    passSoftMuon    &= passKin;
    if(lid == MUON_PDGID){
      if(leptons[ilep].pt()< 10) passLooseLepton=false;
      if(leptons[ilep].pt()< 3)  passSoftMuon=false;
    } else if(lid==ELECTRON_PDGID){
      if(leptons[ilep].pt()<10) passLooseLepton=false;
    }
    if(leptons[ilep].pt()<20) passKin=false;

    //Cut based identification 
    passId = lid==ELECTRON_PDGID ? patUtils::passId(leptons[ilep].el, vtx[0], patUtils::llvvElecId::Tight) : patUtils::passId(leptons[ilep].mu, vtx[0], patUtils::llvvMuonId::Tight);
    
    //isolation
    passIso = lid==11?patUtils::passIso(leptons[ilep].el,  patUtils::llvvElecIso::Tight) : patUtils::passIso(leptons[ilep].mu,  patUtils::llvvMuonIso::Tight);

    if(passId && passIso && passKin)          selLeptons.push_back(leptons[ilep]);
    else if(passLooseLepton || passSoftMuon)  extraLeptons.push_back(leptons[ilep]);
    
  } // end lepton loop
} 

int main(int argc, char* argv[])
{
  // check arguments
  if(argc<2){ std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl; exit(0); }
  
  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();
  
  // configure the process
  const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");
  bool debug = runProcess.getParameter<bool>("debug");
  bool isMC = runProcess.getParameter<bool>("isMC");  
  double xsec = runProcess.getParameter<double>("xsec");
  int mctruthmode=runProcess.getParameter<int>("mctruthmode");
  edm::ParameterSet pujetidparas = runProcess.getParameter<edm::ParameterSet>("pujetidparas"); 
  std::vector<std::string> urls=runProcess.getUntrackedParameter<std::vector<std::string> >("input");
  // TString url = TString(argv[1]);
  // TString outFileUrl(gSystem->BaseName(url));
  // outFileUrl.ReplaceAll("_cfg.py","");
  // if(mctruthmode!=0) { outFileUrl += "_filt"; outFileUrl += mctruthmode; }
  // TString outdir=runProcess.getParameter<std::string>("outdir");
  // TString outUrl( outdir );
  // gSystem->Exec("mkdir -p " + outUrl);

  TString output=runProcess.getParameter<std::string>("output");
  
  // initiating histograms
  SmartSelectionMonitor mon = initHistograms();
  
  // get ready for the event loop
  fwlite::ChainEvent ev(urls);
  const size_t totalEntries= ev.size();

  //MC normalization (to 1/pb)
  double xsecWeight = xsec/totalEntries;

  if(!isMC) xsecWeight=1.0;
  if (debug) {
    printf("DEBUG: xsec= %f\n", xsec);
    printf("DEBUG: xsecWeight = %f\n", xsecWeight);
  }

  //b-tagging: beff and leff must be derived from the MC sample using the discriminator vs flavor
  //the scale factors are taken as average numbers from the pT dependent curves see:
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2012_Data_and_MC_EPS13_prescript
  BTagSFUtil btsfutil;
  float beff(0.68), sfb(0.99), sfbunc(0.015);
  float leff(0.13), sfl(1.05), sflunc(0.12);

  // make sure that histogram internally produced in 
  // lumireweighting are not destroyed when closing the file
  gROOT->cd();  

  higgs::utils::EventCategory eventCategoryInst(higgs::utils::EventCategory::EXCLUSIVE2JETSVBF); //jet(0,>=1)+vbf binning

  // ----------------------------------------------------------------------------------------
  // event loop  
  // loop on all the events
  printf("Progressing Bar     :0%%       20%%       40%%       60%%       80%%       100%%\n");
  printf("Scanning the ntuple :");
  int treeStep(totalEntries/50);
  
  TString tag("all");
  double weight(1.0); // for testing now. 

  for( size_t iev=0; iev<totalEntries; iev++){
    if(iev%treeStep==0){printf(".");fflush(stdout);}
    // load the event content from the EDM file
    ev.to(iev);
    
    //apply trigger and require compatibilitiy of the event with the PD
    float triggerPrescale(1.0); 
    float triggerThreshold(0.0); 
    bool hasPhotonTrigger = patUtils::passPhotonTrigger(ev, triggerThreshold, triggerPrescale);
    
    // only run on the events that pass our triggers
    if( !hasPhotonTrigger ) continue; 
    
    //load all the objects we will need to access
    reco::VertexCollection vtx;
    fwlite::Handle< reco::VertexCollection > vtxHandle; 
    vtxHandle.getByLabel(ev, "offlineSlimmedPrimaryVertices");
    if (vtxHandle.isValid() ) { vtx = *vtxHandle; }
    mon.fillHisto("nvtx", "all", vtx.size(), weight);

    double rho = 0;
    fwlite::Handle< double > rhoHandle;
    rhoHandle.getByLabel(ev, "fixedGridRhoFastjetAll");
    if(rhoHandle.isValid()){ rho = *rhoHandle;}
    mon.fillHisto("rho", "all", rho, weight);

    pat::PhotonCollection photons;
    fwlite::Handle< pat::PhotonCollection > photonsHandle;
    photonsHandle.getByLabel(ev, "slimmedPhotons");
    if(photonsHandle.isValid()){ photons = *photonsHandle;}
    mon.fillHisto("npho", "all", photons.size(), weight);

    pat::JetCollection jets;
    fwlite::Handle< pat::JetCollection > jetsHandle;
    jetsHandle.getByLabel(ev, "slimmedJets");
    if(jetsHandle.isValid()){ jets = *jetsHandle;}
    mon.fillHisto("njet", "all", jets.size(), weight);

    pat::METCollection mets;
    fwlite::Handle< pat::METCollection > metsHandle;
    metsHandle.getByLabel(ev, "slimmedMETs");
    if(metsHandle.isValid()){ mets = *metsHandle;}
    LorentzVector met = mets[0].p4(); 
    
    pat::ElectronCollection electrons;
    fwlite::Handle< pat::ElectronCollection > electronsHandle;
    electronsHandle.getByLabel(ev, "slimmedElectrons");
    if(electronsHandle.isValid()){ electrons = *electronsHandle;}

    pat::MuonCollection muons;
    fwlite::Handle< pat::MuonCollection > muonsHandle;
    muonsHandle.getByLabel(ev, "slimmedMuons");
    if(muonsHandle.isValid()){ muons = *muonsHandle;}
    
    // below follows the analysis of the main selection with n-1 plots
    tag = "sel";
    
    // select photons 
    pat::PhotonCollection selPhotons = passPhotonSelection(mon, photons, triggerThreshold, rho);
    if ( selPhotons.size() == 0) continue;  
    mon.fillHisto("npho", tag, selPhotons.size(), weight);
    mon.fillHisto("nvtx", tag, vtx.size(), weight);

    for(size_t ipho=0; ipho<selPhotons.size(); ipho++) {
      pat::Photon photon = selPhotons[ipho]; 
      mon.fillHisto("phopt", tag, photon.pt(), weight);
      mon.fillHisto("phoeta", tag, photon.superCluster()->eta(), weight);
      // mon.fillHisto("phor9", tag, photon.r9(), weight);
      // mon.fillHisto("phoiso", tag, photon.photonIso(), weight);
      mon.fillHisto("phohoe", tag, photon.hadTowOverEm(), weight);
      mon.fillHisto("elevto", tag, photon.hasPixelSeed(), weight);
      // mon.fillHisto("sigietaieta", tag, photon.sigmaIetaIeta(), weight);
      mon.fillHisto("sigietaieta", tag, photon.userFloat("sigmaIetaIeta_NoZS"), weight);
    }


    // merge electrons and muons
    std::vector<patUtils::GenericLepton> leptons;
    for(size_t l=0;l<electrons.size();l++){leptons.push_back(patUtils::GenericLepton(electrons[l]));}      
    for(size_t l=0;l<muons    .size();l++){leptons.push_back(patUtils::GenericLepton(muons    [l]));}      
    std::sort(leptons.begin(),   leptons.end(), utils::sort_CandidatesByPt);

    // select leptons
    std::vector<patUtils::GenericLepton> selLeptons, extraLeptons;
    passLeptonSelection(mon, leptons, selPhotons, vtx, selLeptons, extraLeptons); 

    std::sort(selLeptons.begin(),   selLeptons.end(), utils::sort_CandidatesByPt);
    std::sort(extraLeptons.begin(), extraLeptons.end(), utils::sort_CandidatesByPt);

    mon.fillHisto("nlep", tag, selLeptons.size(), weight);
    mon.fillHisto("nexlep", tag, extraLeptons.size(), weight);
  
    // select jets
    pat::JetCollection selJets = passJetSelection(mon, jets, pujetidparas,
						  vtx, selLeptons, selPhotons); 
    if ( selJets.size() == 0) continue;  
    int njets(0), nbtags(0); 
    float mindphijmet(9999.);

    for(size_t ijet=0; ijet<selJets.size(); ijet++) {
      pat::Jet jet = selJets[ijet]; 
      mon.fillHisto("jetpt", tag, jet.pt(), weight);
      mon.fillHisto("jeteta", tag, jet.eta(), weight);

      if(jet.pt()>30) {
	njets++;
	float dphijmet=fabs(deltaPhi(met.phi(), jet.phi()));
	if(dphijmet<mindphijmet) mindphijmet=dphijmet;
	if(fabs(jet.eta())<2.5){
	  bool hasCSVtag(jet.bDiscriminator("combinedSecondaryVertexBJetTags")>0.405);
	  //update according to the SF measured by BTV
	  if(isMC){
	    int flavId=jet.partonFlavour();
	    if(abs(flavId)==5)        btsfutil.modifyBTagsWithSF(hasCSVtag,sfb,beff);
	    else if(abs(flavId)==4)   btsfutil.modifyBTagsWithSF(hasCSVtag,sfb/5,beff);
	    else	              btsfutil.modifyBTagsWithSF(hasCSVtag,sfl,leff);
	  }
	  nbtags   += hasCSVtag;
	}
      }
    } 
    
    // met
    mon.fillHisto("met", tag, met.pt(), weight);

    // assign channel based on lepton flavors
    std::vector<TString> chTags;
    int dilId(1);
    LorentzVector boson(0,0,0,0);
    // only consider the photon selection case  
    if(hasPhotonTrigger && selPhotons.size()) {
      dilId=22;
      chTags.push_back("ee");
      chTags.push_back("mumu");
      boson = selPhotons[0].p4();
      weight *= triggerPrescale;
    }

    TString evCat=eventCategoryInst.GetCategory(selJets, boson);
    std::vector<TString> tags(1,"all"); // first element init 
    for(size_t ich=0; ich<chTags.size(); ich++){
      tags.push_back( chTags[ich] );
      tags.push_back( chTags[ich]+evCat );
    }

    mon.fillHisto("eventflow", tags, 0, weight);
    if(chTags.size()==0) continue;

    // baseline selection
    bool passMass(fabs(boson.mass()-91)<15);
    bool passQt(boson.pt()>55);
    bool passThirdLeptonVeto( selLeptons.size()==2 && extraLeptons.size()==0 );
    bool passBtags(nbtags==0);
    bool passMinDphijmet( njets==0 || mindphijmet>0.5);
    
    passMass=hasPhotonTrigger;
    passThirdLeptonVeto=(selLeptons.size()==0 && extraLeptons.size()==0);
      
    mon.fillHisto("eventflow", tags, 1, weight);
    // mon.fillHisto("nvtxraw", tags, vtx.size(), weight/puWeight);
    mon.fillHisto("nvtx", tags, vtx.size(), weight);
    mon.fillHisto("rho", tags, rho, weight);
    mon.fillHisto("zmass", tags,boson.mass(),weight); 
    mon.fillHisto("zy",    tags,fabs(boson.Rapidity()),weight); 
  
    // "|M-91|<15"  
    if(!passMass) continue; 
    mon.fillHisto("eventflow",tags, 2,weight);
    mon.fillHisto("zpt",      tags, boson.pt(),weight);

    //these two are used to reweight photon -> Z, the 3rd is a control
    mon.fillHisto("qt",       tags, boson.pt(),weight,true); 
    mon.fillHisto("qtraw",    tags, boson.pt(),weight/triggerPrescale,true); 

    // "p_{T}>55" 
    if(!passQt) continue; 
    mon.fillHisto("eventflow",tags,3,weight);
    int nExtraLeptons((selLeptons.size()-2)+extraLeptons.size());
    mon.fillHisto("nextraleptons",tags,nExtraLeptons,weight);
    if(nExtraLeptons>0){
      LorentzVector thirdLepton(selLeptons.size()>2 ?  selLeptons[1].p4() : extraLeptons[0].p4());
      double dphi=fabs(deltaPhi(thirdLepton.phi(),met.phi()));
      double mt=TMath::Sqrt(2*thirdLepton.pt()*met.pt()*(1-TMath::Cos(dphi)));
      mon.fillHisto("thirdleptonpt",tags,thirdLepton.pt(),weight);
      mon.fillHisto("thirdleptoneta",tags,fabs(thirdLepton.eta()),weight);
      mon.fillHisto("thirdleptonmt",tags,mt,weight);
    }

    // "3^{rd}-lepton veto" 
    if(!passThirdLeptonVeto) continue; 
    mon.fillHisto("eventflow",tags,4,weight);
    for(size_t ijet=0; ijet<selJets.size(); ijet++){
      if(selJets[ijet].pt()<30 || fabs(selJets[ijet].eta())>2.5) continue;
      
      float csv(selJets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"));
      mon.fillHisto( "csv",tags,csv,weight);
      if(!isMC) continue;
      int flavId=selJets[ijet].partonFlavour();
      TString jetFlav("others");
      if(abs(flavId)==5)      jetFlav="b";
      else if(abs(flavId)==4) jetFlav="c";
      mon.fillHisto( "csv"+jetFlav,tags,csv,weight);
    }
    
    mon.fillHisto( "nbtags",tags,nbtags,weight);
    
    
  } // end event loop 
  printf(" done.\n"); 
  
  //save control plots to file
  // outUrl += "/";
  // outUrl += outFileUrl + ".root";
  // printf("Results saved in %s\n", outUrl.Data());
  // TFile *ofile=TFile::Open(outUrl, "recreate");

  printf("Results saved in %s\n", output.Data());
  TFile *ofile=TFile::Open(output, "recreate");

  mon.Write();
  ofile->Close();
}  







