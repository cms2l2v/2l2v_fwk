// -*- C++ -*-
//
// Package:    TopMassDifferenceAnalyzer
// Class:      TopMassDifferenceAnalyzer
// 
/**\class TopMassDifferenceAnalyzer TopMassDifferenceAnalyzer.cc TopMassDifference/TopMassDifferenceAnalyzer/src/TopMassDifferenceAnalyzer.cc
 
 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Juerg Eugster,40 1-B06,+41227671504,
//         Created:  Thu Jun 20 20:44:25 CEST 2013
// $Id: TopMassDifferenceAnalyzer.cc,v 1.1 2013/06/21 16:13:42 jueugste Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//
// class declaration
//

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "UserCode/llvv_fwk/interface/MacroUtils.h"

#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>


class TopMassDifferenceAnalyzer : public edm::EDAnalyzer {
public:
  explicit TopMassDifferenceAnalyzer(const edm::ParameterSet&);
  ~TopMassDifferenceAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  
  // ----------member data ---------------------------
  TH1F *h_cutflow;
  TH1F *h_ptTop, *h_ptTopFS;
  TH1F *h_mt_gen,           *h_mt_decay,           *h_mt_decay_fs;
  TH1F *h_deltaMt_gen,      *h_deltaMt_decay,      *h_deltaMt_decay_fs;
  TH2F *h_deltaMt_gen_ttpt, *h_deltaMt_decay_ttpt, *h_deltaMt_decay_fs_ttpt;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TopMassDifferenceAnalyzer::TopMassDifferenceAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


TopMassDifferenceAnalyzer::~TopMassDifferenceAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TopMassDifferenceAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace std;
   using namespace edm;

  //gen particles
  edm::Handle< std::vector<reco::GenParticle> > genParticles;
  iEvent.getByLabel("genParticles", genParticles);
  std::vector<reco::GenParticle>::const_iterator genParticle;
  if(!genParticles.isValid())     cerr << "  WARNING: genParticles is not valid! " << endl;

  //gen jets
  edm::Handle< std::vector<reco::GenJet> > genJets;
  iEvent.getByLabel("ak5GenJets", genJets);
  std::vector<reco::GenJet>::const_iterator genJet;
  if(!genJets.isValid())     cerr << "  WARNING: genJets is not valid! " << endl;

  //analyze the generator level jets
  std::map<std::string, TLorentzVector> selParticles;
  for(genParticle = genParticles->begin(); genParticle != genParticles->end(); ++genParticle) {

    if(genParticle->status()!=3) continue;
    int pid=genParticle->pdgId();

    //analyze the top/anti-top
    if(abs(pid)==6){
      std::string name( pid==6 ? "t" : "tbar" );
      selParticles[ name ] = TLorentzVector(genParticle->px(), genParticle->py(), genParticle->pz(), genParticle->energy());

      const reco::Candidate *fsTop = utils::cmssw::getGeneratorFinalStateFor( &(*genParticle),false);
      h_ptTop->Fill(genParticle->pt());
      h_ptTopFS->Fill(fsTop->pt());

      //analyze the decay products of the top
      for(size_t idaughter=0; idaughter<genParticle->numberOfDaughters(); idaughter++){

	int daughter_pid=genParticle->daughter(idaughter)->pdgId();

	//W bosons
	if(abs(daughter_pid)==24){
	  std::string daughter_name( daughter_pid==24 ? "wplus" : "wminus" );
	  selParticles[ daughter_name ] = TLorentzVector( genParticle->daughter(idaughter)->px(), genParticle->daughter(idaughter)->py(), genParticle->daughter(idaughter)->pz(), genParticle->daughter(idaughter)->energy() );
	}

	//b-quarks
	else if(abs(daughter_pid)==5){
	  
	  std::string daughter_name( daughter_pid==5 ? "b" : "bbar" );
	  selParticles[ daughter_name ] = TLorentzVector( genParticle->daughter(idaughter)->px(), genParticle->daughter(idaughter)->py(), genParticle->daughter(idaughter)->pz(), genParticle->daughter(idaughter)->energy() );
	  
	  const reco::Candidate *fsDaughter = utils::cmssw::getGeneratorFinalStateFor(genParticle->daughter(idaughter),false);
	  daughter_name += "_fs";
	  selParticles[ daughter_name ] = TLorentzVector( fsDaughter->px(), fsDaughter->py(), fsDaughter->pz(), fsDaughter->energy() );
	}
      }
    }
  }

  //compute the kinematics of interest
  TLorentzVector t_gen        = selParticles["t"];
  TLorentzVector tbar_gen     = selParticles["tbar"];
  float deltaMt_gen           = t_gen.M()   - tbar_gen.M();

  TLorentzVector t_decay       = selParticles["b"]       + selParticles["wplus"];
  TLorentzVector tbar_decay    = selParticles["bbar"]    + selParticles["wminus"];
  float deltaMt_decay          = t_decay.M()             - tbar_decay.M();

  TLorentzVector t_decay_fs    = selParticles["b_fs"]    + selParticles["wplus"];
  TLorentzVector tbar_decay_fs = selParticles["bbar_fs"] + selParticles["wminus"];
  float deltaMt_decay_fs       = t_decay_fs.M()          - tbar_decay_fs.M();

  TLorentzVector ttbar         = selParticles["t"]+selParticles["tbar"];
  float ttpt                   = ttbar.Pt();

  //fill the histograms
  h_cutflow->Fill(0);

  h_mt_gen->Fill( t_gen.M() );               h_mt_gen->Fill( tbar_gen.M() );
  h_mt_decay->Fill( t_decay.M() );           h_mt_decay->Fill( tbar_decay.M() );
  h_mt_decay_fs->Fill( t_decay_fs.M() );     h_mt_decay_fs->Fill( tbar_decay_fs.M() );

  h_deltaMt_gen->Fill(deltaMt_gen);
  h_deltaMt_decay->Fill(deltaMt_decay);
  h_deltaMt_decay_fs->Fill(deltaMt_decay_fs);

  h_deltaMt_gen_ttpt->Fill(ttpt,deltaMt_gen);
  h_deltaMt_decay_ttpt->Fill(ttpt,deltaMt_decay);
  h_deltaMt_decay_fs_ttpt->Fill(ttpt,deltaMt_decay_fs);

}


// ------------ method called once each job just before starting event loop  ------------
void 
TopMassDifferenceAnalyzer::beginJob()
{
  //book the histograms
  edm::Service<TFileService> fs;
  h_cutflow          = fs->make<TH1F>("cutflow", "cutflow"    ,5,0,5);
  h_ptTop            = fs->make<TH1F>("pttop",        ";Generated top p_{T} [GeV];Events x2",100,100.,250.);
  h_ptTopFS           = fs->make<TH1F>("ptfstop",        ";Generated top p_{T} [GeV];Events x2",100,100.,250.);

  h_mt_gen           = fs->make<TH1F>("mt",        ";Generated top mass [GeV];Events x2",100,100.,250.);
  h_mt_decay         = fs->make<TH1F>("mtdecay",   ";Generated top mass [GeV];Events x2",100,100.,250.);
  h_mt_decay_fs      = fs->make<TH1F>("mtdecayFS", ";Generated top mass [GeV];Events x2",100,100.,250.);
  h_deltaMt_gen      = fs->make<TH1F>("deltaMt",        ";Generated top mass difference [GeV];Events",100,-10.,10);
  h_deltaMt_decay    = fs->make<TH1F>("deltaMtdecay",   ";Generated top mass difference [GeV];Events",100,-10.,10);
  h_deltaMt_decay_fs = fs->make<TH1F>("deltaMtdecayFS", ";Generated top mass difference [GeV];Events",100,-10.,10);
  h_deltaMt_gen_ttpt      = fs->make<TH2F>("deltaMt_ttpt",          ";t#bar{t} transverse momentum [GeV];Generated top mass difference [GeV];Events",25,0,250,100,-10,10);
  h_deltaMt_decay_ttpt    = fs->make<TH2F>("deltaMt_decay_ttpt",    ";t#bar{t} transverse momentum [GeV];Generated top mass difference [GeV];Events",25,0,250,100,-10,10);
  h_deltaMt_decay_fs_ttpt = fs->make<TH2F>("deltaMt_decay_fs_ttpt", ";t#bar{t} transverse momentum [GeV];Generated top mass difference [GeV];Events",25,0,250,100,-10,10);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TopMassDifferenceAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
TopMassDifferenceAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
TopMassDifferenceAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
TopMassDifferenceAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TopMassDifferenceAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TopMassDifferenceAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TopMassDifferenceAnalyzer);
