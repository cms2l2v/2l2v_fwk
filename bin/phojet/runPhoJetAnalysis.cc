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

  // make sure that histogram internally produced in 
  // lumireweighting are not destroyed when closing the file
  gROOT->cd();  

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
    for(size_t ijet=0; ijet<selJets.size(); ijet++) {
      pat::Jet jet = selJets[ijet]; 
      mon.fillHisto("jetpt", tag, jet.pt(), weight);
      mon.fillHisto("jeteta", tag, jet.eta(), weight);
    }
    
    // met
    mon.fillHisto("met", tag, met.pt(), weight);
       
    
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







