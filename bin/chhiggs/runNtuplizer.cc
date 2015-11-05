//
// Pietro Vischia, <pietro.vischia@gmail.com>
//
// ttbar and charged Higgs analyses
//

#include <iostream>
#include <boost/shared_ptr.hpp>

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

//Load here all the dataformat that we will need
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "GeneratorInterface/LHEInterface/interface/LHEEvent.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
//#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h" //for svfit

#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/HiggsUtils.h"
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/llvv_fwk/interface/TMVAUtils.h"
#include "UserCode/llvv_fwk/interface/LeptonEfficiencySF.h"
#include "UserCode/llvv_fwk/interface/PDFInfo.h"
#include "UserCode/llvv_fwk/interface/MuScleFitCorrector.h"
#include "UserCode/llvv_fwk/interface/BtagUncertaintyComputer.h"
#include "UserCode/llvv_fwk/interface/GammaWeightsHandler.h"

#include "UserCode/llvv_fwk/interface/PatUtils.h"

#include "UserCode/llvv_fwk/interface/MiniEvent.h"


#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TEventList.h"
#include "TROOT.h"
#include "TNtuple.h"
#include <Math/VectorUtil.h>

using namespace std;

namespace utils
{
  namespace cmssw
  {
    std::vector<double> smearJER(double pt, double eta, double genPt)
    {
      std::vector<double> toReturn(3,pt);
      if(genPt<=0) return toReturn;
      
      // FIXME: These are the 8 TeV values.
      //
      eta=fabs(eta);
      double ptSF(1.0), ptSF_err(0.06);
      if(eta<0.5)                  { ptSF=1.052; ptSF_err=sqrt(pow(0.012,2)+pow(0.5*(0.062+0.061),2)); }
      else if(eta>=0.5 && eta<1.1) { ptSF=1.057; ptSF_err=sqrt(pow(0.012,2)+pow(0.5*(0.056+0.055),2)); }
      else if(eta>=1.1 && eta<1.7) { ptSF=1.096; ptSF_err=sqrt(pow(0.017,2)+pow(0.5*(0.063+0.062),2)); }
      else if(eta>=1.7 && eta<2.3) { ptSF=1.134; ptSF_err=sqrt(pow(0.035,2)+pow(0.5*(0.087+0.085),2)); }
      else if(eta>=2.3 && eta<5.0) { ptSF=1.288; ptSF_err=sqrt(pow(0.127,2)+pow(0.5*(0.155+0.153),2)); }
      
      toReturn[0]=TMath::Max(0.,(genPt+ptSF*(pt-genPt)));
      toReturn[1]=TMath::Max(0.,(genPt+(ptSF+ptSF_err)*(pt-genPt)));
      toReturn[2]=TMath::Max(0.,(genPt+(ptSF-ptSF_err)*(pt-genPt)));
      return toReturn;
    }

    //
    std::vector<double> smearJES(double pt, double eta, JetCorrectionUncertainty *jecUnc)
    {
      jecUnc->setJetEta(eta);
      jecUnc->setJetPt(pt);
      double relShift=fabs(jecUnc->getUncertainty(true));
      std::vector<double> toRet;
      toRet.push_back((1.0+relShift)*pt);
      toRet.push_back((1.0-relShift)*pt);
      return toRet;
    }
    
    void updateJEC(pat::JetCollection &jets, FactorizedJetCorrector *jesCor, JetCorrectionUncertainty *totalJESUnc, float rho, int nvtx,bool isMC)
    {
      for(size_t ijet=0; ijet<jets.size(); ijet++)
        {
          pat::Jet jet = jets[ijet];
          //correct JES
          LorentzVector rawJet = jet.correctedP4("Uncorrected");
          //double toRawSF=jet.correctedJet("Uncorrected").pt()/jet.pt();
          //LorentzVector rawJet(jet*toRawSF);
          jesCor->setJetEta(rawJet.eta());
          jesCor->setJetPt(rawJet.pt());
          jesCor->setJetA(jet.jetArea());
          jesCor->setRho(rho);
          jesCor->setNPV(nvtx);
          double newJECSF=jesCor->getCorrection();
          rawJet *= newJECSF;
          //jet.SetPxPyPzE(rawJet.px(),rawJet.py(),rawJet.pz(),rawJet.energy());
          jet.setP4(rawJet);

          //smear JER
          double newJERSF(1.0);
          if(isMC)
            {
              const reco::GenJet* genJet=jet.genJet();
              double genjetpt( genJet ? genJet->pt(): 0.);
              std::vector<double> smearJER=utils::cmssw::smearJER(jet.pt(),jet.eta(),genjetpt);
              newJERSF=smearJER[0]/jet.pt();
              rawJet *= newJERSF;
              //jet.SetPxPyPzE(rawJet.px(),rawJet.py(),rawJet.pz(),rawJet.energy());
              jet.setP4(rawJet);
              // FIXME: change the way this is stored (to not storing it)
              // //set the JER up/down alternatives 
              // jets[ijet].setVal("jerup",   smearJER[1] );
              // jets[ijet].setVal("jerdown", smearJER[2] );
            }
      
          // FIXME: change the way this is stored (to not storing it)
          ////set the JES up/down pT alternatives
          //std::vector<float> ptUnc=utils::cmssw::smearJES(jet.pt(),jet.eta(), totalJESUnc);
          //jets[ijet].setVal("jesup",    ptUnc[0] );
          //jets[ijet].setVal("jesdown",  ptUnc[1] );
      
          // FIXME: this is not to be re-set. Check that this is a desired non-feature.
          // i.e. check that the uncorrectedJet remains the same even when the corrected momentum is changed by this routine. 
          //to get the raw jet again
          //jets[ijet].setVal("torawsf",1./(newJECSF*newJERSF));  
        }
    }

    enum METvariations { NOMINAL, JERUP, JERDOWN, JESUP, JESDOWN, UMETUP, UMETDOWN, LESUP, LESDOWN };    
    //
    std::vector<LorentzVector> getMETvariations(LorentzVector &rawMETP4, pat::JetCollection &jets, std::vector<patUtils::GenericLepton> &leptons,bool isMC)
    {
      std::vector<LorentzVector> newMetsP4(9,rawMETP4);
      if(!isMC) return newMetsP4;
      
      LorentzVector nullP4(0,0,0,0);
      //recompute the clustered and unclustered fluxes with energy variations
      for(size_t ivar=1; ivar<=8; ivar++)
        {
          
          //leptonic flux
          LorentzVector leptonFlux(nullP4), lepDiff(nullP4);
          for(size_t ilep=0; ilep<leptons.size(); ilep++) {
            LorentzVector lepton = leptons[ilep].p4();
            double varSign( (ivar==LESUP ? 1.0 : (ivar==LESDOWN ? -1.0 : 0.0) ) );
            int id( abs(leptons[ilep].pdgId()) );
            double sf(1.0);
            if(id==13) sf=(1.0+varSign*0.01);
            if(id==11) {
              if(fabs(leptons[ilep].eta())<1.442) sf=(1.0+varSign*0.02);
              else                                sf=(1.0-varSign*0.05);
            }
            leptonFlux += lepton;
            lepDiff += (sf-1)*lepton;
          }
      
          //clustered flux
          LorentzVector jetDiff(nullP4), clusteredFlux(nullP4);
          for(size_t ijet=0; ijet<jets.size(); ijet++)
            {
              if(jets[ijet].pt()==0) continue;
              double jetsf(1.0);
              // FIXME: change the way this is stored (to not storing it)              
              /// if(ivar==JERUP)   jetsf=jets[ijet].getVal("jerup")/jets[ijet].pt();
              /// if(ivar==JERDOWN) jetsf=jets[ijet].getVal("jerdown")/jets[ijet].pt();
              /// if(ivar==JESUP)   jetsf=jets[ijet].getVal("jesup")/jets[ijet].pt();
              /// if(ivar==JESDOWN) jetsf=jets[ijet].getVal("jesdown")/jets[ijet].pt();
              //LorentzVector newJet( jets[ijet] ); newJet *= jetsf;
              LorentzVector newJet = jets[ijet].p4(); newJet *= jetsf;
              jetDiff       += (newJet-jets[ijet].p4());
              clusteredFlux += jets[ijet].p4();
            }
          LorentzVector iMet=rawMETP4-jetDiff-lepDiff;

          //unclustered flux
          if(ivar==UMETUP || ivar==UMETDOWN)
            {
              LorentzVector unclusteredFlux=-(iMet+clusteredFlux+leptonFlux);
              unclusteredFlux *= (ivar==UMETUP ? 1.1 : 0.9); 
              iMet = -clusteredFlux -leptonFlux - unclusteredFlux;
            }
      
          //save new met
          newMetsP4[ivar]=iMet;
        }
  
      //all done here
      return newMetsP4;
    }

}
}
bool passPFJetID(std::string label,
                 pat::Jet jet){
  
  bool passID(false); 
  
  float rawJetEn(jet.correctedJet("Uncorrected").energy() );

  double eta=jet.eta();
 
  float nhf( (jet.neutralHadronEnergy() + jet.HFHadronEnergy())/rawJetEn );
  float nef( jet.neutralEmEnergy()/rawJetEn );
  float cef( jet.chargedEmEnergy()/rawJetEn );
  float chf( jet.chargedHadronEnergy()/rawJetEn );
  float nch    = jet.chargedMultiplicity();
  float nconst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
  float muf(jet.muonEnergy()/rawJetEn); 

  // Set of cuts from the POG group: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data
  if(label=="Loose")
    passID = ( ((nhf<0.99 && nef<0.99 && nconst>1) && ((abs(eta)<=2.4 && chf>0 && nch>0 && cef<0.99) || abs(eta)>2.4)) && abs(eta)<=3.0 );
  if(label=="Tight")
    passID = ( ((nhf<0.90 && nef<0.90 && nconst>1) && ((abs(eta)<=2.4 && chf>0 && nch>0 && cef<0.90) || abs(eta)>2.4)) && abs(eta) <=3.0);
  
  // Should be added the abs(eta)>3.0 part, but we never consider such jets, so... Meh!
  
  return passID; 
  
}



bool hasLeptonAsDaughter(const reco::GenParticle p)
{
  bool foundL(false);
  if(p.numberOfDaughters()==0) return foundL;

  // cout << "Particle " << p.pdgId() << " with status " << p.status() << " and " << p.numberOfDaughters() << endl;
  const reco::Candidate *part = &p;
  // loop on the daughter particles to check if it has an e/mu as daughter
  while ((part->numberOfDaughters()>0)) {
    const reco::Candidate* DaughterPart = part->daughter(0);
    // cout << "\t\t Daughter: " << DaughterPart->pdgId() << " with status " << DaughterPart->status() << endl;
    if (fabs(DaughterPart->pdgId()) == 11 || fabs(DaughterPart->pdgId() == 13)){
      foundL = true;
      break;
    }
    part=DaughterPart;
  }
  return foundL;
}


bool hasWasMother(const reco::GenParticle  p)
{
  bool foundW(false);
  if(p.numberOfMothers()==0) return foundW;
  const reco::Candidate* part =&p; // (p.mother());
  // loop on the mother particles to check if it has a W as mother
  while ((part->numberOfMothers()>0)) {
    const reco::Candidate* MomPart =part->mother();
    if (fabs(MomPart->pdgId())==24){
      foundW = true;
      break;
    }
    part = MomPart;
  }
  return foundW;
}

bool hasTauAsMother(const reco::GenParticle  p)
{
  bool foundTau(false);
  if(p.numberOfMothers()==0) return foundTau;
  const reco::Candidate* part = &p; //(p.mother());
  // loop on the mother particles to check if it has a tau as mother
  while ((part->numberOfMothers()>0)) {
    const reco::Candidate* MomPart =part->mother();
    if (fabs(MomPart->pdgId())==15)// && MomPart->status() == 2) // Not sure the status check is needed.
      {
        foundTau = true;
        break;
      }
    part = MomPart;
  }
  return foundTau;
}

int main (int argc, char *argv[])
{
  //##############################################
  //########    GLOBAL INITIALIZATION     ########
  //##############################################

  // check arguments
  if (argc < 2)
    {
      std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl;
      exit (0);
    }
  
  // load framework libraries
  gSystem->Load ("libFWCoreFWLite");
  AutoLibraryLoader::enable ();
  
  // configure the process
  const edm::ParameterSet & runProcess = edm::readPSetsFrom (argv[1])->getParameter < edm::ParameterSet > ("runProcess");

  bool debug           = runProcess.getParameter<bool>  ("debug");
  bool runSystematics  = runProcess.getParameter<bool>  ("runSystematics");
  bool saveSummaryTree = runProcess.getParameter<bool>  ("saveSummaryTree");
  bool isMC            = runProcess.getParameter<bool>  ("isMC");
  double xsec          = runProcess.getParameter<double>("xsec");
  int mctruthmode      = runProcess.getParameter<int>   ("mctruthmode");
  TString dtag         = runProcess.getParameter<std::string>("dtag");

  const edm::ParameterSet& myVidElectronIdConf = runProcess.getParameterSet("electronidparas");
  const edm::ParameterSet& myVidElectronMainIdWPConf = myVidElectronIdConf.getParameterSet("tight");
  const edm::ParameterSet& myVidElectronVetoIdWPConf = myVidElectronIdConf.getParameterSet("loose");
  
  VersionedPatElectronSelector electronVidMainId(myVidElectronMainIdWPConf);
  VersionedPatElectronSelector electronVidVetoId(myVidElectronVetoIdWPConf);
  
  TString suffix = runProcess.getParameter < std::string > ("suffix");
  std::vector < std::string > urls = runProcess.getUntrackedParameter < std::vector < std::string > >("input");

  TString outUrl = runProcess.getParameter<std::string>("outfile");
 
  // Good lumi mask
  lumiUtils::GoodLumiFilter goodLumiFilter(runProcess.getUntrackedParameter<std::vector<edm::LuminosityBlockRange> >("lumisToProcess", std::vector<edm::LuminosityBlockRange>()));

  bool
    filterOnlySINGLEE  (false),
    filterOnlySINGLEMU (false);
  if (!isMC)
    {
      if (dtag.Contains ("SingleMuon")) filterOnlySINGLEMU = true;
      if (dtag.Contains ("SingleEle"))  filterOnlySINGLEE  = true;
    }

  // Reactivate when you will use exclusive W+jets and DY+Jets samples, in order to correctly merge the exclusive ones with the inclusive one 
  //bool isV0JetsMC   (isMC && (dtag.Contains ("DYJetsToLL_50toInf") || dtag.Contains ("WJets")));
  // Reactivate for diboson shapes  
  // bool isMC_ZZ      (isMC && (string (dtag.Data ()).find ("TeV_ZZ") != string::npos));
  // bool isMC_WZ      (isMC && (string (dtag.Data ()).find ("TeV_WZ") != string::npos));
  
  bool isTTbarMC    (isMC && (dtag.Contains("TTJets") || dtag.Contains("_TT_") )); // Is this still useful?
  bool isPromptReco (!isMC && dtag.Contains("Run2015B-PromptReco")); //"False" picks up correctly the new prompt reco (2015C) and MC
  bool isRun2015B   (!isMC && dtag.Contains("Run2015B"));
  bool isNLOMC      (isMC && (dtag.Contains("amcatnlo") || dtag.Contains("powheg")) );
  
  TString outTxtUrl = outUrl + ".txt";
  FILE *outTxtFile = NULL;
  if (!isMC) outTxtFile = fopen (outTxtUrl.Data (), "w");
  printf ("TextFile URL = %s\n", outTxtUrl.Data ());
  
  //tree info
  TString dirname = runProcess.getParameter < std::string > ("dirName");
  
  //systematics
  std::vector<TString> systVars(1,"");
  if(runSystematics && isMC)
    {
      systVars.push_back("jerup" );     systVars.push_back("jerdown"   );
      systVars.push_back("jesup" );     systVars.push_back("jesdown"   );
      //systVars.push_back("lesup" );   systVars.push_back("lesdown"   );
      systVars.push_back("leffup");     systVars.push_back("leffdown"  );
      systVars.push_back("puup"  );     systVars.push_back("pudown"   );
      systVars.push_back("umetup");     systVars.push_back("umetdown" );
      systVars.push_back("btagup");     systVars.push_back("btagdown" );
      systVars.push_back("unbtagup");   systVars.push_back("unbtagdown" );
      if(isTTbarMC) {systVars.push_back("topptuncup"); systVars.push_back("topptuncdown"); }
      //      systVars.push_back(); systVars.push_back();

      if(isTTbarMC) { systVars.push_back("pdfup"); systVars.push_back("pdfdown"); }
      cout << "Systematics will be computed for this analysis - this will take a bit" << endl;
    }

  size_t nSystVars(systVars.size());
  
  std::vector < std::string > allWeightsURL = runProcess.getParameter < std::vector < std::string > >("weightsFile");
  std::string weightsDir (allWeightsURL.size ()? allWeightsURL[0] : "");
  
  //##############################################
  //########    INITIATING HISTOGRAMS     ########
  //##############################################
  SmartSelectionMonitor mon;

  //generator level control : add an underflow entry to make sure the histo is kept
  //((TH1F*)mon.addHistogram( new TH1D( "higgsMass_raw",     ";Higgs Mass [GeV];Events", 500,0,1500) ))->Fill(-1.0,0.0001);
  
  // ensure proper normalization
  TH1D* normhist = (TH1D*) mon.addHistogram(new TH1D("initNorm", ";;Events", 5, 0., 5.));
  normhist->GetXaxis()->SetBinLabel (1, "Gen. Events");
  normhist->GetXaxis()->SetBinLabel (2, "Events");
  normhist->GetXaxis()->SetBinLabel (3, "PU central");
  normhist->GetXaxis()->SetBinLabel (4, "PU up");
  normhist->GetXaxis()->SetBinLabel (5, "PU down");
  
  for(Int_t igenjet=0; igenjet<5; igenjet++)
    {
      TString fidtag("fidcounter"); fidtag+=igenjet;
      mon.addHistogram(new TH1D(fidtag, ";Variation;Events", 200, 0., 200.)); 
    }
  


  //
  // STATISTICAL ANALYSIS
  //
  TH1D* Hoptim_systs = (TH1D*) mon.addHistogram (new TH1D ("optim_systs", ";syst;", nSystVars, 0, nSystVars));
  for (size_t ivar=0; ivar<nSystVars; ++ivar) Hoptim_systs->GetXaxis()->SetBinLabel(ivar+1, systVars[ivar]);
  
  //##############################################
  //######## GET READY FOR THE EVENT LOOP ########
  //##############################################
  size_t totalEntries(0);

  TFile *ofile = TFile::Open (outUrl, "recreate");
  
  //TFile* summaryFile = NULL;
  TTree* summaryTree = NULL; //ev->;
  MiniEvent_t miniEvent;

  if(saveSummaryTree)
    {
      TDirectory* cwd = gDirectory;
      //TString outSummaryUrl = runProcess.getParameter<std::string>("summaryfile");
      //if(outSummaryUrl=="") outSummaryUrl=outUrl;
      //std::string summaryFileName(outUrl); 
      //summaryFileName.replace(summaryFileName.find(".root", 0), 5, "_summary.root");
      
      //summaryFile = new TFile(summaryFileName.c_str(), "recreate");
      summaryTree = new TTree("minievents", "minievents");

      createMiniEventTree(summaryTree, miniEvent);
      
      //summaryTree->SetDirectory(summaryFile);  // This line is probably not needed
        
      cwd->cd();
    }
  
  
  //jet energy scale and uncertainties 
  TString jecDir = runProcess.getParameter < std::string > ("jecDir");
  gSystem->ExpandPathName (jecDir);
  FactorizedJetCorrector *jesCor = utils::cmssw::getJetCorrector (jecDir, isMC);
  JetCorrectionUncertainty *totalJESUnc = new JetCorrectionUncertainty ((jecDir + "/MC_Uncertainty_AK4PFchs.txt").Data ());
  
  //muon energy scale and uncertainties
  MuScleFitCorrector *muCor = NULL; // FIXME: MuScle fit corrections for 13 TeV not available yet (more Zs are needed) getMuonCorrector (jecDir, dtag);

  //lepton efficiencies
  LeptonEfficiencySF lepEff;
  
  // b-tagging 
  // Prescriptions taken from: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation74X

  // b-tagging working points for 50ns 
  //   (pfC|c)ombinedInclusiveSecondaryVertexV2BJetTags
  //      v2CSVv2L 0.605
  //      v2CSVv2M 0.890
  //      v2CSVv2T 0.970
  double
    btagLoose(0.605),
    btagMedium(0.890),
    btagTight(0.970);

  //b-tagging: scale factors
  //beff and leff must be derived from the MC sample using the discriminator vs flavor
  //the scale factors are taken as average numbers from the pT dependent curves see:
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2012_Data_and_MC_EPS13_prescript
  BTagSFUtil btsfutil;
  double beff(0.68), sfb(0.99), sfbunc(0.015);
  double leff(0.13), sfl(1.05), sflunc(0.12);


  // Btag SF and eff from https://indico.cern.ch/event/437675/#preview:1629681
  sfb = 0.861;
  // sbbunc =;
  beff = 0.559;


  TString
    electronIdMainTag("cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
    electronIdVetoTag("cutBasedElectronID-Spring15-25ns-V1-standalone-tight");

  //pileup weighting
  edm::LumiReWeighting * LumiWeights = NULL;
  utils::cmssw::PuShifter_t PuShifters;
  double PUNorm[] = { 1, 1, 1 };
  if (isMC)
    {
      std::vector<double> dataPileupDistributionDouble = runProcess.getParameter < std::vector < double >>("datapileup");
      std::vector<float> dataPileupDistribution;
      for (unsigned int i = 0; i < dataPileupDistributionDouble.size (); i++)
        {
          dataPileupDistribution.push_back (dataPileupDistributionDouble[i]);
        }
      std::vector<float> mcPileupDistribution;
      utils::getMCPileupDistributionFromMiniAODtemp(urls, dataPileupDistribution.size (), mcPileupDistribution);
      while(mcPileupDistribution.size() < dataPileupDistribution.size()) mcPileupDistribution.push_back(0.0);
      while(mcPileupDistribution.size() > dataPileupDistribution.size()) dataPileupDistribution.push_back(0.0);
      gROOT->cd ();             //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE
      LumiWeights = new edm::LumiReWeighting(mcPileupDistribution, dataPileupDistribution);
      PuShifters = utils::cmssw::getPUshifters(dataPileupDistribution, 0.05);
      utils::getPileupNormalization(mcPileupDistribution, PUNorm, LumiWeights, PuShifters);
    }
  
  gROOT->cd ();                 //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE
  
  //higgs::utils::EventCategory eventCategoryInst(higgs::utils::EventCategory::EXCLUSIVE2JETSVBF); //jet(0,>=1)+vbf binning
  
  //##############################################
  //########           EVENT LOOP         ########
  //##############################################
  //loop on all the events
  printf ("Progressing Bar     :0%%       20%%       40%%       60%%       80%%       100%%\n");
  
  int nMultiChannel(0);
  
  for(size_t f=0; f<urls.size();++f){
    TFile* file = TFile::Open(urls[f].c_str());
    fwlite::Event ev(file);
    printf ("Scanning the ntuple %2lu/%2lu :",f+1, urls.size());
    int iev(0);
    int treeStep (ev.size()/50);
    //DuplicatesChecker duplicatesChecker;
    //int nDuplicates(0);
    for (ev.toBegin(); !ev.atEnd(); ++ev)
      {
        iev++;
        totalEntries++;
        if (iev % treeStep == 0)
          {
            printf (".");
            if(!debug) fflush (stdout); // Otherwise debug messages are flushed
          }

        if(saveSummaryTree){
          // Save event header (in case I need to prepare an event list)
          miniEvent.run    = ev.id().run();
          miniEvent.lumi   = ev.luminosityBlock();
          miniEvent.event  = ev.id().event(); 
          miniEvent.isData = ev.isRealData();
        }


        edm::EventBase const & myEvent = ev;
        // Take into account the negative weights from some NLO generators (otherwise some phase space will be double counted)
        double weightGen(1.);

        std::vector < TString > tags (1, "all"); // Inclusive inclusiveness

        //GENERATOR LEVEL INFO
        miniEvent.isFiducial = true;  
        miniEvent.ttbar_nw=0;
        miniEvent.me_np=0;
        miniEvent.ngenj=0;
        miniEvent.ttbar_genId=0;

        
        if(isNLOMC)
          {
            int ngenJets(-1);
            //double weightGen(0.);
            //double weightLhe(0.);
            if(saveSummaryTree)
              {
                fwlite::Handle<int> genTtbarIdHandle;
                genTtbarIdHandle.getByLabel(ev, "categorizeGenTtbar", "genTtbarId");
                if(genTtbarIdHandle.isValid()) miniEvent.ttbar_genId=*genTtbarIdHandle;


                fwlite::Handle<reco::GenParticleCollection> genParticlesHandle;
                genParticlesHandle.getByLabel(ev, "prunedGenParticles");
                
                //require only one lepton (can be from tau, if tau not from hadron)
                int nLeptons(0);
                float lphi(0), leta(0);
                for (size_t ipart = 0; ipart < genParticlesHandle->size(); ++ipart) {
                  const reco::GenParticle & genIt = (*genParticlesHandle)[ipart];
                  if(!genIt.isPromptFinalState() && !genIt.isDirectPromptTauDecayProductFinalState()) continue;
                  int ID = abs(genIt.pdgId());
                  if(ID!=11 && ID!=13) continue;
                  if(genIt.pt()<20 || fabs(genIt.eta())>2.5) continue;
                  nLeptons++;
                  lphi=genIt.phi();
                  leta=genIt.eta();
                }
                if(nLeptons!=1) ngenJets=0;
                
                //require 1 jets not overlapping with lepton
                fwlite::Handle<std::vector<reco::GenJet> > genJetsHandle;
                genJetsHandle.getByLabel(ev, "ak4GenJetsCustom");
                for(std::vector<reco::GenJet>::const_iterator genJet=genJetsHandle->begin(); genJet!=genJetsHandle->end(); genJet++)
                  {
                    if(genJet->pt()<20 || fabs(genJet->eta())>2.5) continue;
                    float dR=deltaR(genJet->eta(),genJet->phi(),leta,lphi);
                    if(dR<0.4) continue;
                    miniEvent.genj_pt  [miniEvent.ngenj]=genJet->pt();
                    miniEvent.genj_eta [miniEvent.ngenj]=genJet->eta();
                    miniEvent.genj_phi [miniEvent.ngenj]=genJet->phi();
                    miniEvent.genj_mass[miniEvent.ngenj]=genJet->mass();
                    miniEvent.ngenj++;
                  }
                ngenJets = miniEvent.ngenj;
                
              }
            
            fwlite::Handle<GenEventInfoProduct> evt;
            evt.getByLabel(ev, "generator");
            if(evt.isValid())
              {
                weightGen = (evt->weight() > 0 ) ? 1. : -1. ;

                if(saveSummaryTree)
                  {
                    miniEvent.ttbar_allmepartons   = evt->nMEPartons();
                    miniEvent.ttbar_matchmepartons = evt->nMEPartonsFiltered();
                    miniEvent.ttbar_w[0]           = evt->weight();
                    miniEvent.ttbar_nw++;
                  }
              }
            
            if(saveSummaryTree)
              {
                fwlite::Handle<LHEEventProduct> lheEvtProd;
                lheEvtProd.getByLabel(ev, "externalLHEProducer");
                if(lheEvtProd.isValid())
                  {
                    double weightLhe=lheEvtProd->originalXWGTUP();

                    for(unsigned int i=0; i<lheEvtProd->weights().size();++i)
                      {
                        double asdde=lheEvtProd->weights()[i].wgt;
                        miniEvent.ttbar_w[miniEvent.ttbar_nw]=miniEvent.ttbar_w[0]*asdde/weightLhe;
                        miniEvent.ttbar_nw++;
                      }
                    
                    const lhef::HEPEUP &hepeup=lheEvtProd->hepeup();
                    miniEvent.me_id=hepeup.IDPRUP;
                    for(int ip=0;ip<hepeup.NUP; ++ip)
                      {
                        miniEvent.me_pid [miniEvent.me_np]=hepeup.IDUP[ip];
                        miniEvent.me_px  [miniEvent.me_np]=hepeup.PUP[ip][0];
                        miniEvent.me_py  [miniEvent.me_np]=hepeup.PUP[ip][1];
                        miniEvent.me_pz  [miniEvent.me_np]=hepeup.PUP[ip][2];
                        miniEvent.me_mass[miniEvent.me_np]=hepeup.PUP[ip][4];
                        miniEvent.me_np++;
                      }
                  }
                
                for(Int_t igenjet=0; igenjet<5; igenjet++)
                  {
                    TString fidtag("fidcounter"); fidtag+=igenjet;
                    mon.fillHisto(fidtag, tags, 0., miniEvent.ttbar_w[0]);
                    if(igenjet<=ngenJets)
                      {
                        for(Int_t iw=1; iw<miniEvent.ttbar_nw; iw++)
                          mon.fillHisto(fidtag, tags, double(iw), miniEvent.ttbar_w[iw]);
                      }
                  }
              }
            
          }
        
        
        //
        // DERIVE WEIGHTS TO APPLY TO SAMPLE
        //
        
        //pileup weight
        double weight           (1.0);
        double rawWeight        (1.0);
        double TotalWeight_plus (1.0);
        double TotalWeight_minus(1.0);
        double puWeight         (1.0);
        
        if(isNLOMC){
          weight *= weightGen;
          rawWeight *=weightGen;
        }

        reco::VertexCollection vtx;
        reco::Vertex goodPV;
        unsigned int nGoodPV(0);
        fwlite::Handle<reco::VertexCollection> vtxHandle;
        vtxHandle.getByLabel(ev, "offlineSlimmedPrimaryVertices");
        if(vtxHandle.isValid() ) vtx = *vtxHandle;
        // Clean up vertex collection
        for(size_t ivtx=0; ivtx<vtx.size(); ++ivtx)
          {
            if(utils::isGoodVertex(vtx[ivtx]))
              {
                if(nGoodPV==0) goodPV=vtx[ivtx];
                nGoodPV++;
              }
          }

        if(saveSummaryTree) miniEvent.nvtx = nGoodPV;
        if(nGoodPV==0) continue; // Do not store/analyze events without any primary vertex


        // Apply pileup reweighting
        if(isMC)
          {
            int ngenITpu(0), ngenITtruepu(0);
            fwlite::Handle < std::vector < PileupSummaryInfo > >puInfoH;
            puInfoH.getByLabel (ev, "slimmedAddPileupInfo");
            for (std::vector < PileupSummaryInfo >::const_iterator it = puInfoH->begin (); it != puInfoH->end (); it++)
              {
                if (it->getBunchCrossing () == 0){
                  ngenITpu += it->getPU_NumInteractions();
                  ngenITtruepu += it->getTrueNumInteractions();
                }
              }
            
            if(saveSummaryTree)
              {
                miniEvent.pu=ngenITpu;
                miniEvent.putrue=ngenITtruepu;
              }

            //ngenITpu = nGoodPV; // based on nvtx
            puWeight = LumiWeights->weight (ngenITpu) * PUNorm[0];
            weight *= puWeight;//Weight; //* puWeight;
            TotalWeight_plus =  PuShifters[utils::cmssw::PUUP]  ->Eval (ngenITpu) * (PUNorm[2]/PUNorm[0]);
            TotalWeight_minus = PuShifters[utils::cmssw::PUDOWN]->Eval (ngenITpu) * (PUNorm[1]/PUNorm[0]);
          }
        
        
        mon.fillHisto("initNorm", tags, 0., weightGen); // Should be all 1, but for NNLO samples there are events weighting -1
        mon.fillHisto("initNorm", tags, 1., weightGen); // Should be all 1, but for NNLO samples there are events weighting -1
        mon.fillHisto("initNorm", tags, 2., puWeight);
        mon.fillHisto("initNorm", tags, 3., TotalWeight_plus);
        mon.fillHisto("initNorm", tags, 4., TotalWeight_minus);

        if(saveSummaryTree) miniEvent.weight=weightGen;
        
        //##############################################   EVENT LOOP STARTS   ##############################################

	// Not needed anymore if you run on 08Oct2015
        ///// Orthogonalize Run2015B PromptReco+17Jul15 mix
        ///if(isRun2015B)
        ///  {
        ///    if(!patUtils::exclusiveDataEventFilter(ev.eventAuxiliary().run(), isMC, isPromptReco ) ) continue;
        ///  }
        
        // Skip bad lumi
        if(!goodLumiFilter.isGoodLumi(ev.eventAuxiliary().run(),ev.eventAuxiliary().luminosityBlock())) continue; 
        
        //apply trigger and require compatibilitiy of the event with the PD
        edm::TriggerResultsByName tr = ev.triggerResultsByName ("HLT");
        if (!tr.isValid ()){
          cout << "Trigger is not valid" << endl;
          return false;
        }
        
        if(debug){
          cout << "Printing trigger list" << endl;
          for(edm::TriggerNames::Strings::const_iterator trnames = tr.triggerNames().begin(); trnames!=tr.triggerNames().end(); ++trnames)
            cout << *trnames << endl;
          cout << "----------- End of trigger list ----------" << endl;
          return 0;
        }

        
        // Need either to simulate the HLT (https://twiki.cern.ch/twiki/bin/view/CMS/TopTrigger#How_to_easily_emulate_HLT_paths) to match triggers.
        bool elTrigger   (
                          isMC ? 
                          utils::passTriggerPatterns (tr, "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v*")
                          :
                          utils::passTriggerPatterns (tr, "HLT_Ele23_WPLoose_Gsf_v*")
                          );
        bool muTrigger   (
                          isMC ? 
                          utils::passTriggerPatterns (tr, "HLT_IsoMu17_eta2p1_v*")
                          :
                          utils::passTriggerPatterns (tr, "HLT_IsoMu18_v*")
                          );

        if(!isMC && muTrigger) mon.fillHisto("nvtx_singlemu_pileup", tags, nGoodPV, 1.);
        if(!isMC && elTrigger) mon.fillHisto("nvtx_singlee_pileup",  tags, nGoodPV, 1.);
        
        if(filterOnlySINGLEMU) {                    elTrigger = false; }
        if(filterOnlySINGLEE)  { muTrigger = false;                    }


        if(saveSummaryTree)
          {
            miniEvent.elTrigger=elTrigger;
            miniEvent.muTrigger=muTrigger;
          }
        
        if (!(elTrigger || muTrigger)) continue;         //ONLY RUN ON THE EVENTS THAT PASS OUR TRIGGERS

        // ------------ Apply MET filters ------------
        if(!patUtils::passMetFilters(ev, isPromptReco)) continue;
        
        
        //load all the objects we will need to access

        double rho = 0;
        fwlite::Handle<double> rhoHandle;
        rhoHandle.getByLabel(ev, "fixedGridRhoFastjetAll");
        if(rhoHandle.isValid() ) rho = *rhoHandle;
        
        if(saveSummaryTree) miniEvent.rho = rho;

        reco::GenParticleCollection gen;
        fwlite::Handle<reco::GenParticleCollection> genHandle;
        genHandle.getByLabel(ev, "prunedGenParticles");
        if(genHandle.isValid() ) gen = *genHandle;
        
        // Save time and don't load the rest of the objects when selecting by mctruthmode :)
        bool hasTop(false);
        int
          ngenLeptonsStatus3(0),
          ngenLeptonsNonTauSonsStatus3(0),
          ngenTausStatus3(0),
          ngenQuarksStatus3(0);
        //double tPt(0.), tbarPt(0.); // top pt reweighting - dummy value results in weight equal to 1 if not set in loop
        //float wgtTopPt(1.0), wgtTopPtUp(1.0), wgtTopPtDown(1.0);
        if(isMC)
          {
            // FIXME: Considering add support for different generators (based on PYTHIA6) for comparison.
            for(size_t igen=0; igen<gen.size(); igen++){
              // FIXME: Should pass to the new status scheme from: https://github.com/cms-sw/cmssw/pull/7791
              //            ////// if(!gen[igen].isHardProcess() && !gen[igen].isPromptFinalState()) continue;
              
              if(gen[igen].status() != 1 &&  gen[igen].status() !=2 && gen[igen].status() !=62 ) continue;
              int absid=abs(gen[igen].pdgId());
              // OK, so taus should be checked as status 2, and quarks as 71 or 23. More testing needed
              //if( absid==15 && hasWasMother(gen[igen]) ) cout << "Event " << iev << ", Particle " << igen << " has " << gen[igen].numberOfDaughters() << " daughters, pdgId " << gen[igen].pdgId() << " and status " << gen[igen].status() << ", mothers " << gen[igen].numberOfMothers() << ", pt " << gen[igen].pt() << ", eta " << gen[igen].eta() << ", phi " << gen[igen].phi() << ". isHardProcess is " << gen[igen].isHardProcess() << ", and isPromptFinalState is " << gen[igen].isPromptFinalState() << endl;
              
              
              //////            if(absid==6 && gen[igen].isHardProcess()){ // particles of the hardest subprocess 22 : intermediate (intended to have preserved mass)
              if(absid==6 && gen[igen].status()==62){ // particles of the hardest subprocess 22 : intermediate (intended to have preserved mass). Josh says 62 (last in chain)
                hasTop=true;

                // FIXME: Top pT reweighting. 13 TeV values not propagated yet, so not using.
                //if(isTTbarMC){
                //  if(gen[igen].get("id") > 0) tPt=gen[igen].pt();
                //  else                        tbarPt=gen[igen].pt();
                //}
              } 
              
              
              //if(!gen[igen].isPromptFinalState() ) continue;
              if( (gen[igen].status() != 1 && gen[igen].status()!= 2 ) || !hasWasMother(gen[igen])) continue;
              
              if((absid==11 || absid==13) && hasLeptonAsDaughter(gen[igen])) cout << "Electron or muon " << igen << " has " << gen[igen].numberOfDaughters() << " daughter which is a lepton." << endl;
              
              if((absid==11 || absid==13) && gen[igen].status()==1)
                {
                  ngenLeptonsStatus3++;
                  
                  if(!hasTauAsMother(gen[igen]))
                    ngenLeptonsNonTauSonsStatus3++;
                }
              if(absid==15 && gen[igen].status()==2 )
                {
                  ngenTausStatus3++; // This should be summed to ngenLeptonsStatus3 for the dilepton final states, not summed for the single lepton final states.
                  //    if(hasLeptonAsDaughter(gen[igen])) cout << "Tau " << igen << " has " << gen[igen].numberOfDaughters() << " daughter which is a lepton." << endl;
                }
              if(absid<=5              ) ngenQuarksStatus3++;
            }
            
            if(debug && (ngenTausStatus3==1 && ngenLeptonsStatus3==1 )  ) cout << "Event: " << iev << ". Leptons: " << ngenLeptonsStatus3 << ". Leptons notaus: " << ngenLeptonsNonTauSonsStatus3 << ". Taus: " << ngenTausStatus3 << ". Quarks: " << ngenQuarksStatus3 << endl;
            
            // Dileptons:
            //    ttbar dileptons --> 1
            //    ttbar other     --> 2
            if(mctruthmode==1 && (ngenLeptonsStatus3+ngenTausStatus3!=2 || !hasTop )) continue;
            if(mctruthmode==2 && (ngenLeptonsStatus3+ngenTausStatus3==2 || !hasTop )) continue;
            // FIXME: port tt+bb splitting from 8 TeV (check the reference to the matched genjet)
            //if(mcTruthMode==1 && (ngenLeptonsStatus3!=2 || !hasTop || ngenBQuarksStatus23>=4)) continue;
            //if(mcTruthMode==2 && (ngenLeptonsStatus3==2 || !hasTop || ngenBQuarksStatus23>=4)) continue;
            //if(mcTruthMode==3 && (ngenBQuarksStatus23<4 || !hasTop))                           continue;
            
            // lepton-tau:
            //    ttbar ltau      --> 3
            //    ttbar dileptons --> 4
            //    ttbar ljets     --> 5
            //    ttbar hadrons   --> 6
            if(mctruthmode==3 && (ngenLeptonsNonTauSonsStatus3!=1 || ngenTausStatus3!=1  || !hasTop )) continue; // This is bugged, as it is obvious
            if(mctruthmode==4 && (ngenLeptonsNonTauSonsStatus3!=2                        || !hasTop )) continue;
            if(mctruthmode==5 && (ngenLeptonsNonTauSonsStatus3+ngenTausStatus3!=1        || !hasTop )) continue;
            
            bool isHad(false);
            if( 
               (ngenLeptonsNonTauSonsStatus3!=1 || ngenTausStatus3!=1 ) &&
               (ngenLeptonsNonTauSonsStatus3!=2                      ) &&
               (ngenLeptonsNonTauSonsStatus3+ngenTausStatus3!=1      ) 
                )
              isHad=true;
            
            //if(mctruthmode==6 && (ngenLeptonsNonTauSonsStatus3!=0 || ngenTausStatus3!=0  || !hasTop )) continue;
            if(mctruthmode==6 && (!isHad || !hasTop )) continue;
            
          }
        if(debug) cout << "DEBUG: Event was not stopped by the ttbar sample categorization (either success, or it was not ttbar)" << endl;      
        
        // FIXME: Top pT reweighting to be reactivated as soon as corrections are released
        //      if(tPt>0 && tbarPt>0 && topPtWgt)
        //        {
        //          topPtWgt->computeWeight(tPt,tbarPt);
        //          topPtWgt->getEventWeight(wgtTopPt, wgtTopPtUp, wgtTopPtDown);
        //          wgtTopPtUp /= wgtTopPt;
        //          wgtTopPtDown /= wgtTopPt;
        //        }
        
        
        
        pat::MuonCollection muons;
        fwlite::Handle<pat::MuonCollection> muonsHandle;
        muonsHandle.getByLabel(ev, "slimmedMuons");
        if(muonsHandle.isValid() ) muons = *muonsHandle;
        
        pat::ElectronCollection electrons;
        fwlite::Handle<pat::ElectronCollection> electronsHandle;
        electronsHandle.getByLabel(ev, "slimmedElectrons");
        if(electronsHandle.isValid() ) electrons = *electronsHandle;
        
        pat::JetCollection jets;
        fwlite::Handle<pat::JetCollection>jetsHandle;
        jetsHandle.getByLabel(ev, "slimmedJets");
        if(jetsHandle.isValid() ) jets = *jetsHandle;

        pat::PhotonCollection photons;
        fwlite::Handle<pat::PhotonCollection> photonsHandle;
        photonsHandle.getByLabel(ev, "slimmedPhotons");
        if(photonsHandle.isValid() ) photons = *photonsHandle;
        
        pat::METCollection mets;
        fwlite::Handle<pat::METCollection> metsHandle;
        metsHandle.getByLabel(ev, "slimmedMETs");
        if(metsHandle.isValid() ) mets = *metsHandle;
        LorentzVector met = mets[0].p4 ();
        
        if(debug){
          // MET try:
          double mypt = mets[0].shiftedPt(pat::MET::METUncertainty::JetEnUp);
          cout << "MET = " << mets[0].pt() << ", JetEnUp: " << mypt << endl;
          LorentzVector myshiftedMet = mets[0].shiftedP4(pat::MET::METUncertainty::JetEnUp);
          cout << "MET = " << mets[0].pt() << ", JetEnUp: " << myshiftedMet.pt() << endl;
        }
        
        pat::TauCollection taus;
        fwlite::Handle<pat::TauCollection> tausHandle;
        tausHandle.getByLabel(ev, "slimmedTaus");
        if(tausHandle.isValid() ) taus = *tausHandle;
        
        // Reactivate when you will use exclusive W+jets and DY+Jets samples, in order to correctly merge the exclusive ones with the inclusive one
        // if (isV0JetsMC)
        //   {
        //     fwlite::Handle < LHEEventProduct > lheEPHandle;
        //     lheEPHandle.getByLabel (ev, "externalLHEProducer");
        //     mon.fillHisto ("nup", "", lheEPHandle->hepeup ().NUP, 1);
        //     if (lheEPHandle->hepeup ().NUP > 5)  continue;
        //     mon.fillHisto ("nupfilt", "", lheEPHandle->hepeup ().NUP, 1);
        //   }
        
        
        //
        //
        // BELOW FOLLOWS THE ANALYSIS OF THE MAIN SELECTION WITH N-1 PLOTS. Whatever that means
        //
        //
        
        
        
        //
        // LEPTON ANALYSIS
        //
        
        //start by merging electrons and muons
        std::vector<patUtils::GenericLepton> leptons;
        for(size_t l=0; l<electrons.size(); ++l) leptons.push_back(patUtils::GenericLepton (electrons[l] ));
        for(size_t l=0; l<muons.size(); ++l)     leptons.push_back(patUtils::GenericLepton (muons[l]     ));
        std::sort(leptons.begin(), leptons.end(), utils::sort_CandidatesByPt);

        LorentzVector muDiff(0., 0., 0., 0.);
        std::vector<patUtils::GenericLepton> selLeptons;
        for(size_t ilep=0; ilep<leptons.size (); ++ilep)
          {
            patUtils::GenericLepton& lepton = leptons[ilep];

            bool 
              passKin(true),     passId(true),     passIso(true);
            
            int lid(lepton.pdgId());
            
            //apply muon corrections
            if(abs(lid) == 13)
            {
              if(muCor)
                {
                  TLorentzVector p4(lepton.px(), lepton.py(), lepton.pz(), lepton.energy());
                  muCor->applyPtCorrection(p4, lid < 0 ? -1 : 1);
                  if(isMC) muCor->applyPtSmearing(p4, lid < 0 ? -1 : 1, false);
                  muDiff -= lepton.p4();
                  lepton.setP4(LorentzVector(p4.Px(), p4.Py(), p4.Pz(), p4.E()));
                  muDiff += lepton.p4();
                }
            }
            
          //no need for charge info any longer
          lid = abs(lid);
          TString lepStr(lid == 13 ? "mu" : "e");
          
          // no need to mess with photon ID // //veto nearby photon (loose electrons are many times photons...)
          // no need to mess with photon ID // double minDRlg(9999.);
          // no need to mess with photon ID // for(size_t ipho=0; ipho<selPhotons.size(); ipho++)
          // no need to mess with photon ID //   minDRlg=TMath::Min(minDRlg,deltaR(leptons[ilep].p4(),selPhotons[ipho].p4()));
          // no need to mess with photon ID // if(minDRlg<0.1) continue;
          
          //kinematics
          double leta(fabs(lid==11 ? lepton.el.superCluster()->eta() : lepton.eta()));
          
          // Main leptons kin
          if(lepton.pt() < 20.)                      passKin = false;
          if(leta > 2.1)                                    passKin = false;
          if(lid == 11 && (leta > 1.4442 && leta < 1.5660)) passKin = false; // Crack veto
          
          
          //Cut based identification 
          
          passId = lid == 11 ? patUtils::passId(electronVidVetoId, myEvent, lepton.el) : patUtils::passId(lepton.mu, goodPV, patUtils::llvvMuonId::StdLoose);

          //isolation
          passIso = lid == 11 ? true : patUtils::passIso(lepton.mu, patUtils::llvvMuonIso::Loose); // Electron iso is included within the ID

          if     (passKin     && passId     && passIso)     selLeptons.push_back(lepton);
          
        }
      std::sort(selLeptons.begin(),   selLeptons.end(),   utils::sort_CandidatesByPt);

      if(saveSummaryTree)
        {
          miniEvent.nl=0;

          for(size_t ilep=0; ilep<selLeptons.size(); ++ilep)
            {
              patUtils::GenericLepton& lepton = selLeptons[ilep];
              int lid(fabs(lepton.pdgId()));
              
              const reco::GenParticle* gen= lid==11 ? lepton.el.genLepton() : lepton.mu.genLepton();
              miniEvent.isPromptFinalState[miniEvent.nl] = gen ? int(gen->isPromptFinalState()) : 0;
              miniEvent.isDirectPromptTauDecayProductFinalState[miniEvent.nl] = gen ? int(gen->isDirectPromptTauDecayProductFinalState()) : 0;
              miniEvent.l_id[miniEvent.nl]=lepton.pdgId();
              miniEvent.l_tightId[miniEvent.nl]= lid == 11 ? int(patUtils::passId(electronVidMainId, myEvent, lepton.el)) : int(patUtils::passId(lepton.mu, goodPV, patUtils::llvvMuonId::StdTight));
              miniEvent.l_tightIso[miniEvent.nl] = lid == 11 ? 1 : int(patUtils::passIso(lepton.mu, patUtils::llvvMuonIso::Tight)); // Electron iso is included within the ID
              miniEvent.l_charge[miniEvent.nl]= lid==11 ? lepton.el.charge() : lepton.mu.charge();
              miniEvent.l_pt[miniEvent.nl]=lepton.pt();
              miniEvent.l_eta[miniEvent.nl]= lid==11 ? lepton.el.superCluster()->eta() : lepton.eta();
              miniEvent.l_phi[miniEvent.nl]=lepton.phi();
              miniEvent.l_mass[miniEvent.nl]=lepton.mass();
              miniEvent.nl++;
            }
        }


      LorentzVector recoMET = met;// FIXME REACTIVATE IT - muDiff;

      if(saveSummaryTree)
        {
          miniEvent.met_pt=met.pt();
          miniEvent.met_phi=met.phi();
        }
      
      //select the taus
      pat::TauCollection selTaus;
      int ntaus (0);
      for (size_t itau = 0; itau < taus.size(); ++itau)
        {
          pat::Tau& tau = taus[itau];
          if (tau.pt() < 20. || fabs (tau.eta()) > 2.3) continue;
          
          bool overlapWithLepton(false);
          for(int l=0; l<(int)selLeptons.size();++l){
            if(reco::deltaR(tau, selLeptons[l])<0.4){overlapWithLepton=true; break;}
          }
          if(overlapWithLepton) continue;
          
          //      if(!tau.isPFTau()) continue; // Only PFTaus // It should be false for slimmedTaus
          //      if(tau.emFraction() >=2.) continue;
          
          // Discriminators from https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV
          // "The tau passes the discriminator if pat::Tau::tauID("name") returns a value of 0.5 or greater"
          if(tau.tauID("decayModeFindingNewDMs")<0.5) continue; // High pt tau. Otherwise, OldDMs
          // Anyways, the collection of taus from miniAOD should be already afer decayModeFinding cut (the tag - Old or New - is unspecified in the twiki, though).
          // Consequently, there might be a small bias due to events that are cut by the OldDM and would not be cut by the NewDM
          if (tau.tauID ("byMediumCombinedIsolationDeltaBetaCorr3Hits")<0.5) continue; // See whether to us the new byMediumPileupWeightedIsolation3Hits that is available only for dynamic strip reconstruction (default in CMSSW_7_4_14)
          if (tau.tauID ("againstMuonTight3")                          <0.5) continue; // Medium working point not usable. Available values: Loose, Tight
          if (tau.tauID ("againstElectronMediumMVA5")                  <0.5) continue; // Tight working point not usable. Avaiable values: VLoose, Loose, Medium
          
          // Pixel hits cut (will be available out of the box in new MINIAOD production)
          {
            int nChHadPixelHits = 0;
            reco::CandidatePtrVector chCands = tau.signalChargedHadrCands();
            for(reco::CandidatePtrVector::const_iterator iter = chCands.begin(); iter != chCands.end(); iter++){
              pat::PackedCandidate const* packedCand = dynamic_cast<pat::PackedCandidate const*>(iter->get());
              int pixelHits = packedCand->numberOfPixelHits();
              if(pixelHits > nChHadPixelHits) nChHadPixelHits = pixelHits;
            }
            if(nChHadPixelHits==0) continue;
          }
          /////
          
          selTaus.push_back(tau);
          ntaus++;
        }
      std::sort (selTaus.begin(), selTaus.end(), utils::sort_CandidatesByPt);
      
      
      if(saveSummaryTree)
        {
          miniEvent.nt=0;
          for(size_t itau=0; itau<selTaus.size(); ++itau)
            {
              pat::Tau& tau = selTaus[itau];
              const reco::GenJet* gen = tau.genJet();
              miniEvent.t_charge[miniEvent.nt]=tau.charge();
              miniEvent.t_pt  [miniEvent.nt]=tau.pt();
              miniEvent.t_eta [miniEvent.nt]=tau.eta();
              miniEvent.t_phi [miniEvent.nt]=tau.phi();
              miniEvent.t_mass[miniEvent.nt]=tau.mass(); // No svfit yet
              reco::CandidatePtr leadChargedHadron = tau.leadChargedHadrCand();
              miniEvent.t_leadChHadP [miniEvent.nt]=leadChargedHadron->p();
              miniEvent.t_leadChHadPt[miniEvent.nt]=leadChargedHadron->pt();
              miniEvent.t_energy     [miniEvent.nt]=tau.energy();
              miniEvent.t_et         [miniEvent.nt]=tau.et();
              miniEvent.nt++;
            }
        }


      //
      //JET/MET ANALYSIS
      //
      if(debug) cout << "Now update Jet Energy Corrections" << endl;
      //add scale/resolution uncertainties and propagate to the MET      
      utils::cmssw::updateJEC(jets,jesCor,totalJESUnc,rho,nGoodPV,isMC);
      if(debug) cout << "Update also MET" << endl;
      std::vector<LorentzVector> newMet=utils::cmssw::getMETvariations(met/*recoMet*/,jets,selLeptons,isMC); // FIXME: Must choose a lepton collection. Perhaps loose leptons?
      met=newMet[utils::cmssw::METvariations::NOMINAL];
      if(debug) cout << "Jet Energy Corrections updated" << endl;
      
      // Select the jets. I need different collections because of tau cleaning, but this is needed only for the single lepton channels, so the tau cleaning is performed later.
      pat::JetCollection
        selJets, selBJets;
      double mindphijmet (9999.);
      for (size_t ijet = 0; ijet < jets.size(); ++ijet)
        {
          pat::Jet& jet = jets[ijet];
          
          if (jet.pt() < 15 || fabs (jet.eta()) > 3.0) continue; // Was 4.7 in eta. Tightened for computing time. 3.0 ensures that we don't cut associations with leptons (0.4 from 2.4)
          
          //mc truth for this jet
          const reco::GenJet * genJet = jet.genJet();
          TString jetType (genJet && genJet->pt() > 0 ? "truejetsid" : "pujetsid");

          //cross-clean with selected leptons and photons
          double minDRlj (9999.), minDRlg (9999.), minDRljSingleLep(9999.);

          for (size_t ilep = 0; ilep < selLeptons.size(); ilep++)
            minDRlj = TMath::Min(minDRlj, reco::deltaR (jet, selLeptons[ilep]));
          // don't want to mess with photon ID // for(size_t ipho=0; ipho<selPhotons.size(); ipho++)
          // don't want to mess with photon ID //   minDRlg = TMath::Min( minDRlg, deltaR(jets[ijet],selPhotons[ipho]) );
          //          if (minDRlj < 0.4 /*|| minDRlg<0.4 */ ) continue;
          
          for (size_t ilep = 0; ilep < selLeptons.size(); ilep++)
            minDRljSingleLep = TMath::Min(minDRljSingleLep, reco::deltaR (jet, selLeptons[ilep]));
          
          //jet id
          bool passPFloose = passPFJetID("Loose", jet); 
          // FIXME: check when pileup ID will come out
          //if (jets[ijet].pt() > 30)
          //  {
          //    mon.fillHisto (jetType, "", fabs (jets[ijet].eta()), 0);
          //    if (passPFloose)                        mon.fillHisto (jetType, "", fabs (jets[ijet].eta()), 1);
          //    if (passLooseSimplePuId)                mon.fillHisto (jetType, "", fabs (jets[ijet].eta()), 2);
          //    if (passPFloose && passLooseSimplePuId) mon.fillHisto (jetType, "", fabs (jets[ijet].eta()), 3);
          //  }
          if (!passPFloose || jet.pt() <25. || fabs(jet.eta()) > 2.5) continue;
          if (minDRlj < 0.4) continue;
          
          selJets.push_back(jet);
          
          double dphijmet = fabs (deltaPhi (met.phi(), jet.phi()));
          if (dphijmet < mindphijmet) mindphijmet = dphijmet;
          bool hasCSVtag(jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > btagMedium);
          
	  // Leave btags status update for the offline analysis
          /*
	    if (isMC)
            {
              int flavId = jet.partonFlavour();
              if      (abs(flavId)==5) btsfutil.modifyBTagsWithSF(hasCSVtag, sfb,   beff);
              else if (abs(flavId)==4) btsfutil.modifyBTagsWithSF(hasCSVtag, sfb/5, beff);
              else                     btsfutil.modifyBTagsWithSF(hasCSVtag, sfl,   leff);
            }
	  */
          if(!hasCSVtag) continue;
          selBJets.push_back(jet);
        }

      std::sort (selJets.begin(),  selJets.end(),  utils::sort_CandidatesByPt);
      std::sort (selBJets.begin(), selBJets.end(), utils::sort_CandidatesByPt);

      if(saveSummaryTree)
        {

          miniEvent.nj=0;
          for(size_t ijet=0; ijet<selJets.size(); ++ijet)
            {
              pat::Jet jet = selJets[ijet];
              //save jet
              const reco::Candidate* genParton = jet.genParton();
              const reco::GenJet* genJet=jet.genJet(); 
              miniEvent.j_area[miniEvent.nj]=jet.jetArea();
              miniEvent.j_pt  [miniEvent.nj]=jet.pt();
              miniEvent.j_mass[miniEvent.nj]=jet.mass();
              miniEvent.j_eta [miniEvent.nj]=jet.eta();
              miniEvent.j_phi [miniEvent.nj]=jet.phi();
              miniEvent.genj_pt[miniEvent.nj]=genJet ? genJet->pt() : 0;
              miniEvent.genj_mass[miniEvent.nj]=genJet ? genJet->mass() : 0;
              miniEvent.genj_eta[miniEvent.nj]=genJet ? genJet->eta() : 0;
              miniEvent.genj_phi[miniEvent.nj]=genJet ?  genJet->phi() : 0;
              bool hasCSVtag(jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > btagMedium);
              // Leave b-tag status update to the offline offline analysis
	      /*
if (isMC)
                {
                  int flavId = jet.partonFlavour();
                  if      (abs(flavId)==5) btsfutil.modifyBTagsWithSF(hasCSVtag, sfb,   beff);
                  else if (abs(flavId)==4) btsfutil.modifyBTagsWithSF(hasCSVtag, sfb/5, beff);
                  else                     btsfutil.modifyBTagsWithSF(hasCSVtag, sfl,   leff);
                }
	      */
              miniEvent.j_csv[miniEvent.nj]=jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
              miniEvent.j_isbtag    [miniEvent.nj]= int(hasCSVtag);
              miniEvent.j_flav      [miniEvent.nj]=jet.partonFlavour();
              miniEvent.j_hadflav   [miniEvent.nj]=jet.hadronFlavour();
              miniEvent.j_pid       [miniEvent.nj]=genParton ? genParton->pdgId() : 0;
              miniEvent.nj++;
            }
	}

      if(saveSummaryTree && selLeptons.size() > 0 && selJets.size()>2)
        {
          TDirectory* cwd = gDirectory;
          ofile->cd();
          summaryTree->Fill();
          cwd->cd();
        }
      
      } // End single file event loop 
    printf("\n");
    delete file;
  } // End loop on files
  

  if(nMultiChannel>0) cout << "Warning! There were " << nMultiChannel << " multi-channel events out of " << totalEntries << " events!" << endl;
  printf ("\n");

  //##############################################
  //########     SAVING HISTO TO FILE     ########
  //##############################################
  //save control plots to file
  printf ("Results save in %s\n", outUrl.Data());



  //save all to the file
  TDirectory* cwd = gDirectory;
  ofile->cd();
  mon.Write();
  if(saveSummaryTree)  summaryTree->Write();
  ofile->Close();
  delete ofile;
  cwd->cd();
  
  if (outTxtFile)
    fclose (outTxtFile);

  // Now that everything is done, dump the list of lumiBlock that we processed in this job
  if(!isMC){
    goodLumiFilter.FindLumiInFiles(urls);
    goodLumiFilter.DumpToJson(((outUrl.ReplaceAll(".root",""))+".json").Data());
  }
  
}
