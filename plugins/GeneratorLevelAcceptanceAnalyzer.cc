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


class GeneratorLevelAcceptanceAnalyzer : public edm::EDAnalyzer {
public:
  explicit GeneratorLevelAcceptanceAnalyzer(const edm::ParameterSet&);
  ~GeneratorLevelAcceptanceAnalyzer();
  
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
  TH1F *h_cutflow, *h_cutflow1b, *h_ptTop;
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
GeneratorLevelAcceptanceAnalyzer::GeneratorLevelAcceptanceAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


GeneratorLevelAcceptanceAnalyzer::~GeneratorLevelAcceptanceAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
GeneratorLevelAcceptanceAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace std;
   using namespace edm;

   h_cutflow->Fill(0);
   h_cutflow1b->Fill(0);

   //
   // gen particles
   //
   edm::Handle< std::vector<reco::GenParticle> > genParticles;
   iEvent.getByLabel("genParticles", genParticles);
   if(!genParticles.isValid())     cerr << "  WARNING: genParticles is not valid! " << endl;
   TLorentzVector neutFlux(0,0,0,0);
   std::vector<int> chLeptonsId; 
   std::vector<TLorentzVector> chLeptons;
   std::vector<TLorentzVector> bQuarks;
   for(std::vector<reco::GenParticle>::const_iterator genParticle = genParticles->begin(); genParticle != genParticles->end(); ++genParticle) 
     {
       int status=genParticle->status();
       int pid=genParticle->pdgId();
       if(status==3 && abs(pid)==6) h_ptTop->Fill(genParticle->pt());
       if(status==3 && abs(pid)==5) {
	 TLorentzVector p4( genParticle->px(), genParticle->py(), genParticle->pz(), genParticle->energy() );
	 bQuarks.push_back(p4);
       }
       if(status==1)
	 {
	   TLorentzVector p4( genParticle->px(), genParticle->py(), genParticle->pz(), genParticle->energy() );
	   if(abs(pid)==11 || abs(pid)==13) {
	     if(p4.Pt()<20 || fabs(p4.Eta()>2.5) ) continue;
	     chLeptonsId.push_back( pid );
	     chLeptons.push_back(p4);
	   }
	   if(abs(pid)==12 || abs(pid)==14 || abs(pid)==16) neutFlux += p4;
	 }
     }

   //determine channel and leading leptons
   int leadId(0),trailerId(0);
   float leadPt(0),trailerPt(0);
   for(size_t i=0; i<chLeptons.size(); i++)
     {
       if(chLeptons[i].Pt()>leadPt)         { trailerId=leadId;         trailerPt=leadPt;            leadId=chLeptonsId[i]; leadPt=chLeptons[i].Pt(); }
       else if(chLeptons[i].Pt()>trailerPt) { trailerId=chLeptonsId[i]; trailerPt=chLeptons[i].Pt();                                                  }
     }
   int ch(leadId*trailerId);
   
   //
   // gen jets
   //
   edm::Handle< std::vector<reco::GenJet> > genJets;
   iEvent.getByLabel("ak5GenJets", genJets);
   if(!genJets.isValid())     cerr << "  WARNING: genJets is not valid! " << endl;
   std::vector<TLorentzVector> jets;
   std::vector<TLorentzVector> bJets;
   for(std::vector<reco::GenJet>::const_iterator genJet=genJets->begin(); genJet!=genJets->end(); genJet++)
     {
       TLorentzVector p4( genJet->px(), genJet->py(), genJet->pz(), genJet->energy() );
       if(p4.Pt()<25 || fabs(p4.Eta())>2.5) continue;

       bool matchesLepton(false);
       for(size_t i=0; i<chLeptons.size(); i++)
	 {
	   float dR=p4.DeltaR(chLeptons[i]);
	   if(dR>0.4) continue;
	   matchesLepton=true;
	   break;
	 }
       if(matchesLepton) continue;

       bool matchesBquark(false);
       for(size_t i=0; i<bQuarks.size(); i++)
	 {
	   float dR=p4.DeltaR(bQuarks[i]);
	   if(dR>0.4) continue;
	   matchesBquark=true;
	 }

       jets.push_back(p4);
       if(matchesBquark) bJets.push_back(p4);
     }

   //compute acceptance
   bool passLeptons(abs(ch)==11*11 || abs(ch)==11*13 || abs(ch)==13*13);
   bool passOS(ch<0);
   bool passMET(true);
   if(abs(ch)==11*11 || abs(ch)==13*13) passMET = (neutFlux.Pt()>30);

   if(passLeptons                                       ) h_cutflow->Fill(1);
   if(passLeptons                   && passMET          ) h_cutflow->Fill(2);
   if(passLeptons                   && passMET && passOS) h_cutflow->Fill(3);
   if(passLeptons && jets.size()>=1                     ) h_cutflow->Fill(4);
   if(passLeptons && jets.size()>=1 && passMET          ) h_cutflow->Fill(5);
   if(passLeptons && jets.size()>=1 && passMET && passOS) h_cutflow->Fill(6);
   if(passLeptons && jets.size()>=2                     ) h_cutflow->Fill(7);
   if(passLeptons && jets.size()>=2 && passMET          ) h_cutflow->Fill(8);
   if(passLeptons && jets.size()>=2 && passMET && passOS) h_cutflow->Fill(9);
   if(passLeptons && jets.size()>=3                     ) h_cutflow->Fill(10);
   if(passLeptons && jets.size()>=3 && passMET          ) h_cutflow->Fill(11);
   if(passLeptons && jets.size()>=3 && passMET && passOS) h_cutflow->Fill(12);
   if(passLeptons && jets.size()>=4                     ) h_cutflow->Fill(13);
   if(passLeptons && jets.size()>=4 && passMET          ) h_cutflow->Fill(14);
   if(passLeptons && jets.size()>=4 && passMET && passOS) h_cutflow->Fill(15);
   if(passLeptons && jets.size()>=5                     ) h_cutflow->Fill(16);
   if(passLeptons && jets.size()>=5 && passMET          ) h_cutflow->Fill(17);
   if(passLeptons && jets.size()>=5 && passMET && passOS) h_cutflow->Fill(18);

   if(passLeptons &&                   bJets.size()>=1                      ) h_cutflow1b->Fill(1);
   if(passLeptons &&                   bJets.size()>=1 && passMET           ) h_cutflow1b->Fill(2);
   if(passLeptons &&                   bJets.size()>=1 && passMET && passOS ) h_cutflow1b->Fill(3);
   if(passLeptons && jets.size()>=2 && bJets.size()>=1                      ) h_cutflow1b->Fill(4);
   if(passLeptons && jets.size()>=2 && bJets.size()>=1 && passMET           ) h_cutflow1b->Fill(5);
   if(passLeptons && jets.size()>=2 && bJets.size()>=1 && passMET && passOS ) h_cutflow1b->Fill(6);
   if(passLeptons && jets.size()>=3 && bJets.size()>=1                      ) h_cutflow1b->Fill(7);
   if(passLeptons && jets.size()>=3 && bJets.size()>=1 && passMET           ) h_cutflow1b->Fill(8);
   if(passLeptons && jets.size()>=3 && bJets.size()>=1 && passMET && passOS ) h_cutflow1b->Fill(9);
   if(passLeptons && jets.size()>=4 && bJets.size()>=1                      ) h_cutflow1b->Fill(10);
   if(passLeptons && jets.size()>=4 && bJets.size()>=1 && passMET           ) h_cutflow1b->Fill(11);
   if(passLeptons && jets.size()>=4 && bJets.size()>=1 && passMET && passOS ) h_cutflow1b->Fill(12);
   if(passLeptons && jets.size()>=5 && bJets.size()>=1                      ) h_cutflow1b->Fill(13);
   if(passLeptons && jets.size()>=5 && bJets.size()>=1 && passMET           ) h_cutflow1b->Fill(14);
   if(passLeptons && jets.size()>=5 && bJets.size()>=1 && passMET && passOS ) h_cutflow1b->Fill(15);
}


// ------------ method called once each job just before starting event loop  ------------
void 
GeneratorLevelAcceptanceAnalyzer::beginJob()
{
  //book the histograms
  edm::Service<TFileService> fs;
 
  h_ptTop            = fs->make<TH1F>("pttop",        ";Generated top transverse momentum [GeV];Events x2",100,0.,250.);

  h_cutflow          = fs->make<TH1F>("cutflow", "cutflow"    ,19,0,19);
  h_cutflow->GetXaxis()->SetBinLabel(1,"Generated");
  h_cutflow->GetXaxis()->SetBinLabel(2,"2l       ");
  h_cutflow->GetXaxis()->SetBinLabel(3,"2l  MET  ");
  h_cutflow->GetXaxis()->SetBinLabel(4,"2l  METOS");
  h_cutflow->GetXaxis()->SetBinLabel(5,"2l1j     ");
  h_cutflow->GetXaxis()->SetBinLabel(6,"2l1jMETOS");
  h_cutflow->GetXaxis()->SetBinLabel(7,"2l1j   OS");
  h_cutflow->GetXaxis()->SetBinLabel(8,"2l2j     ");
  h_cutflow->GetXaxis()->SetBinLabel(9,"2l2jMET  ");
  h_cutflow->GetXaxis()->SetBinLabel(10,"2l2jMETOS");
  h_cutflow->GetXaxis()->SetBinLabel(11,"2l3j     ");
  h_cutflow->GetXaxis()->SetBinLabel(12,"2l3jMET  ");
  h_cutflow->GetXaxis()->SetBinLabel(13,"2l3jMETOS");
  h_cutflow->GetXaxis()->SetBinLabel(14,"2l4j     ");
  h_cutflow->GetXaxis()->SetBinLabel(15,"2l4jMET  ");
  h_cutflow->GetXaxis()->SetBinLabel(16,"2l4jMETOS");
  h_cutflow->GetXaxis()->SetBinLabel(17,"2l5j     ");
  h_cutflow->GetXaxis()->SetBinLabel(18,"2l5jMET  ");
  h_cutflow->GetXaxis()->SetBinLabel(19,"2l5jMETOS");

  h_cutflow1b = fs->make<TH1F>("cutflow1b", "cutflow1b"    ,16,0,16);
  h_cutflow1b->GetXaxis()->SetBinLabel(1, "Generated  ");
  h_cutflow1b->GetXaxis()->SetBinLabel(2, "2l1j1b     ");
  h_cutflow1b->GetXaxis()->SetBinLabel(3, "2l1j1bMET  ");
  h_cutflow1b->GetXaxis()->SetBinLabel(4, "2l1j1bMETOS");
  h_cutflow1b->GetXaxis()->SetBinLabel(5, "2l2j1b     ");
  h_cutflow1b->GetXaxis()->SetBinLabel(6, "2l2j1bMET  ");
  h_cutflow1b->GetXaxis()->SetBinLabel(7, "2l2j1bMETOS");
  h_cutflow1b->GetXaxis()->SetBinLabel(8, "2l3j1b     ");
  h_cutflow1b->GetXaxis()->SetBinLabel(9, "2l3j1bMET  ");
  h_cutflow1b->GetXaxis()->SetBinLabel(10,"2l3j1bMETOS");
  h_cutflow1b->GetXaxis()->SetBinLabel(11,"2l4j1b     ");
  h_cutflow1b->GetXaxis()->SetBinLabel(12,"2l4j1bMET  ");
  h_cutflow1b->GetXaxis()->SetBinLabel(13,"2l4j1bMETOS");
  h_cutflow1b->GetXaxis()->SetBinLabel(14,"2l5j1b     ");
  h_cutflow1b->GetXaxis()->SetBinLabel(15,"2l5j1bMET  ");
  h_cutflow1b->GetXaxis()->SetBinLabel(16,"2l5j1bMETOS");
}


// ------------ method called once each job just after ending the event loop  ------------
void 
GeneratorLevelAcceptanceAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
GeneratorLevelAcceptanceAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
GeneratorLevelAcceptanceAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
GeneratorLevelAcceptanceAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
GeneratorLevelAcceptanceAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GeneratorLevelAcceptanceAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GeneratorLevelAcceptanceAnalyzer);
