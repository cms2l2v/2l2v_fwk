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

#include "EgammaAnalysis/ElectronTools/interface/ElectronEnergyCalibratorRun2.h"  
#include "EgammaAnalysis/ElectronTools/interface/PhotonEnergyCalibratorRun2.h" 


#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/HiggsUtils.h"
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/llvv_fwk/interface/TMVAUtils.h"
#include "UserCode/llvv_fwk/interface/rochcor2015.h"
#include "UserCode/llvv_fwk/interface/muresolution_run2.h"
#include "UserCode/llvv_fwk/interface/LeptonEfficiencySF.h"
#include "UserCode/llvv_fwk/interface/PDFInfo.h"
#include "UserCode/llvv_fwk/interface/MuScleFitCorrector.h"
#include "UserCode/llvv_fwk/interface/BTagCalibrationStandalone.h"
#include "UserCode/llvv_fwk/interface/BtagUncertaintyComputer.h"
#include "UserCode/llvv_fwk/interface/GammaWeightsHandler.h"

#include "UserCode/llvv_fwk/interface/PatUtils.h"


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
  //TString baseDir = runProcess.getParameter < std::string > ("dirName");
  //  if (mctruthmode != 0) //FIXME
  //    {
  //      outFileUrl += "_filt";
  //      outFileUrl += mctruthmode;
  //    }
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

  bool isV0JetsMC   (isMC && (dtag.Contains ("DYJetsToLL") || dtag.Contains ("WJets")));
  // Reactivate for diboson shapes  
  // bool isMC_ZZ      (isMC && (string (dtag.Data ()).find ("TeV_ZZ") != string::npos));
  // bool isMC_WZ      (isMC && (string (dtag.Data ()).find ("TeV_WZ") != string::npos));
  
  bool isTTbarMC    (isMC && (dtag.Contains("TTJets") || dtag.Contains("_TT_") )); // Is this still useful?
  bool isPromptReco (!isMC && dtag.Contains("PromptReco")); //"False" picks up correctly the new prompt reco (2015C) and MC
  bool isRun2015B   (!isMC && dtag.Contains("Run2015B"));
  bool isNLOMC      (isMC && (dtag.Contains("amcatnlo") || dtag.Contains("powheg")) );
  
  TString outTxtUrl = outUrl + ".txt";
  FILE *outTxtFile = NULL;
  if (!isMC) outTxtFile = fopen (outTxtUrl.Data (), "w");
  printf ("TextFile URL = %s\n", outTxtUrl.Data ());
  
  //tree info
  TString dirname = runProcess.getParameter < std::string > ("dirName");
  
  //systematics

  std::vector<string> jetVarNames = {""};;//, "_scale_jup","_scale_jdown", "_res_jup", "_res_jdown"}; // Whaddafuck
  
  std::vector<TString> systVars(1,"");
  if(runSystematics && isMC)
    {
      systVars.push_back("jerup" );     systVars.push_back("jerdown"    );
      systVars.push_back("jesup" );     systVars.push_back("jesdown"    );
      systVars.push_back("mesup" );     systVars.push_back("mesdown"    ); // Muon energy scale
      systVars.push_back("eesup" );     systVars.push_back("eesdown"    ); // Electron energy scale
      systVars.push_back("leffup");     systVars.push_back("leffdown"   );
      systVars.push_back("puup"  );     systVars.push_back("pudown"     );
      systVars.push_back("umetup");     systVars.push_back("umetdown"   );
      systVars.push_back("btagup");     systVars.push_back("btagdown"   );
      systVars.push_back("unbtagup");   systVars.push_back("unbtagdown" );
      systVars.push_back("thfactup");   systVars.push_back("thfactdown" );
      systVars.push_back("pdfup"); systVars.push_back("pdfdown"); 
      if(isTTbarMC) {systVars.push_back("topptuncup"); systVars.push_back("topptuncdown"); }

      cout << "Systematics will be computed for this analysis - this will take a bit" << endl;
    }

  size_t nSystVars(systVars.size());
  
  std::vector < std::string > allWeightsURL = runProcess.getParameter < std::vector < std::string > >("weightsFile");
  std::string weightsDir (allWeightsURL.size ()? allWeightsURL[0] : "");
  
  //  //shape uncertainties for dibosons
  //  std::vector<TGraph *> vvShapeUnc;
  //  if(isMC_ZZ || isMC_WZ)
  //    {
  //      TString weightsFile=weightsDir+"/zzQ2unc.root";
  //      TString dist("zzpt");
  //      if(isMC_WZ) { weightsFile.ReplaceAll("zzQ2","wzQ2"); dist.ReplaceAll("zzpt","wzpt"); }
  //      gSystem->ExpandPathName(weightsFile);
  //      TFile *q2UncF=TFile::Open(weightsFile);
  //      if(q2UncF){
  //    vvShapeUnc.push_back( new TGraph( (TH1 *)q2UncF->Get(dist+"_up") ) );
  //    vvShapeUnc.push_back( new TGraph( (TH1 *)q2UncF->Get(dist+"_down") ) );
  //    q2UncF->Close();
  //      }
  //    }


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

  mon.addHistogram( new TH1F( "metFilter_eventflow",     ";metEventflow",20,0,20) );

  //event selection - charged Higgs
  TH1D* h = (TH1D*) mon.addHistogram (new TH1D ("chhiggseventflowdilep", ";;Events", 6, 0., 6.));
  h->GetXaxis()->SetBinLabel (1, "#geq 2 iso leptons");
  h->GetXaxis()->SetBinLabel (2, "M_{ll} veto");
  h->GetXaxis()->SetBinLabel (3, "#geq 2 jets");
  h->GetXaxis()->SetBinLabel (4, "E_{T}^{miss}");
  h->GetXaxis()->SetBinLabel (5, "op. sign");
  h->GetXaxis()->SetBinLabel (6, "#geq 2 b-tags");
  h = (TH1D*) mon.addHistogram(new TH1D ("chhiggsalteventflowdilep", ";;Events", 8, 0., 8.));
  h->GetXaxis()->SetBinLabel (1, "#geq 2 iso leptons");
  h->GetXaxis()->SetBinLabel (2, "M_{ll} veto");
  h->GetXaxis()->SetBinLabel (3, "#geq 2 jets");
  h->GetXaxis()->SetBinLabel (4, "E_{T}^{miss}");
  h->GetXaxis()->SetBinLabel (5, "op. sign");
  h->GetXaxis()->SetBinLabel (6, "= 0 b-tags");
  h->GetXaxis()->SetBinLabel (7, "= 1 b-tags");
  h->GetXaxis()->SetBinLabel (8, "#geq 2 b-tags");
  h = (TH1D*) mon.addHistogram(new TH1D ("chhiggseventflowslep", ";;Events", 6, 0., 6.));
  h->GetXaxis()->SetBinLabel (1, "1 iso lepton");
  h->GetXaxis()->SetBinLabel (2, "#geq 2 jets");
  h->GetXaxis()->SetBinLabel (3, "E_{T}^{miss}");
  h->GetXaxis()->SetBinLabel (4, "#geq 1 b-tag");
  h->GetXaxis()->SetBinLabel (5, "1 #tau");
  h->GetXaxis()->SetBinLabel (6, "op. sign");
  h = (TH1D*) mon.addHistogram(new TH1D("chhiggsalteventflowslep", ";;Events", 7, 0., 7.));
  h->GetXaxis()->SetBinLabel (1, "1 iso lepton");
  h->GetXaxis()->SetBinLabel (2, "E_{T}^{miss}");
  h->GetXaxis()->SetBinLabel (3, "1 #tau");
  h->GetXaxis()->SetBinLabel (4, "op. sign");
  h->GetXaxis()->SetBinLabel (5, "= 0 b-tags");
  h->GetXaxis()->SetBinLabel (6, "= 1 b-tags");
  h->GetXaxis()->SetBinLabel (7, "#geq 2 b-tags");

  // event selection - cross section
  h = (TH1D*) mon.addHistogram (new TH1D ("xseceventflowdilep", ";;Events", 6, 0., 6.));
  h->GetXaxis()->SetBinLabel (1, "#geq 2 iso leptons");
  h->GetXaxis()->SetBinLabel (2, "M_{ll} veto");
  h->GetXaxis()->SetBinLabel (3, "#geq 2 jets");
  h->GetXaxis()->SetBinLabel (4, "E_{T}^{miss}");
  h->GetXaxis()->SetBinLabel (5, "op. sign");
  h->GetXaxis()->SetBinLabel (6, "#geq 2 b-tags");
  h = (TH1D*) mon.addHistogram(new TH1D ("xsecalteventflowdilep", ";;Events", 8, 0., 8.));
  h->GetXaxis()->SetBinLabel (1, "#geq 2 iso leptons");
  h->GetXaxis()->SetBinLabel (2, "M_{ll} veto");
  h->GetXaxis()->SetBinLabel (3, "#geq 2 jets");
  h->GetXaxis()->SetBinLabel (4, "E_{T}^{miss}");
  h->GetXaxis()->SetBinLabel (5, "op. sign");
  h->GetXaxis()->SetBinLabel (6, "= 0 b-tags");
  h->GetXaxis()->SetBinLabel (7, "= 1 b-tags");
  h->GetXaxis()->SetBinLabel (8, "#geq 2 b-tags");
  h = (TH1D*) mon.addHistogram(new TH1D ("xseceventflowslep", ";;Events", 6, 0., 6.));
  h->GetXaxis()->SetBinLabel (1, "1 iso lepton");
  h->GetXaxis()->SetBinLabel (2, "#geq 2 jets");
  h->GetXaxis()->SetBinLabel (3, "E_{T}^{miss}");
  h->GetXaxis()->SetBinLabel (4, "#geq 1 b-tag");
  h->GetXaxis()->SetBinLabel (5, "1 #tau");
  h->GetXaxis()->SetBinLabel (6, "op. sign");
  h = (TH1D*) mon.addHistogram(new TH1D("xsecalteventflowslep", ";;Events", 7, 0., 7.));
  h->GetXaxis()->SetBinLabel (1, "1 iso lepton");
  h->GetXaxis()->SetBinLabel (2, "E_{T}^{miss}");
  h->GetXaxis()->SetBinLabel (3, "1 #tau");
  h->GetXaxis()->SetBinLabel (4, "op. sign");
  h->GetXaxis()->SetBinLabel (5, "= 0 b-tags");
  h->GetXaxis()->SetBinLabel (6, "= 1 b-tags");
  h->GetXaxis()->SetBinLabel (7, "#geq 2 b-tags");

  h = (TH1D*) mon.addHistogram( new TH1D("nvtx_pileup"         , ";;Events", 100, 0., 100.));
  h = (TH1D*) mon.addHistogram( new TH1D("nvtx_singlemu_pileup", ";;Events", 100, 0., 100.));
  h = (TH1D*) mon.addHistogram( new TH1D("nvtx_singlee_pileup" , ";;Events", 100, 0., 100.));

    
  // Setting up control categories to be analyzed
  std::vector < TString > controlCats;
  controlCats.clear ();
  controlCats.push_back("step1");
  controlCats.push_back("step2");
  controlCats.push_back("step3");
  controlCats.push_back("step4");
  controlCats.push_back("step5");
  controlCats.push_back("step6");
  
  controlCats.push_back("altstep1");
  controlCats.push_back("altstep2");
  controlCats.push_back("altstep3");
  controlCats.push_back("altstep4");
  controlCats.push_back("altstep5");
  controlCats.push_back("altstep6");
  controlCats.push_back("altstep7");
  controlCats.push_back("altstep8");

  for (size_t k = 0; k < controlCats.size (); ++k)
    {
      TString icat (controlCats[k]);

      //pu control to be completed
      mon.addHistogram (new TH1D (icat+"nvtx",    ";Vertices;Events", 100, 0., 100.));
      mon.addHistogram (new TH1D (icat+"nvtxraw", ";Vertices;Events", 100, 0., 100.));
      mon.addHistogram (new TH1D (icat+"rho", ";#rho;Events"        ,  50, 0.,  25.));
      

      //tau control to be completed
      TH1* htaus = mon.addHistogram (new TH1D (icat + "ntaus", ";Tau multiplicity;Events", 5, 0., 5.));
      for (int ibin = 1; ibin <= htaus->GetXaxis ()->GetNbins (); ibin++)
        {
          TString label ("");
          if (ibin == h->GetXaxis ()->GetNbins ())
            label += "#geq";
          else
            label += "=";
          label += (ibin - 1);
          htaus->GetXaxis ()->SetBinLabel (ibin, label);
        }
      mon.addHistogram( new TH1D(icat+"tauleadpt",              ";p_{T}^{#tau};Events", 30,  0.,  300.  ));
      mon.addHistogram( new TH1D(icat+"tauleadeta",             ";#eta^{#tau};Events",  50, -2.6,   2.6 ));
      mon.addHistogram( new TH1D(icat+"tauinclusivept",         ";p_{T}^{#tau};Events", 30,  0.,  300.  ));
      mon.addHistogram( new TH1D(icat+"tauinclusiveeta",        ";#eta^{#tau};Events",  50, -2.6,   2.6 ));
      mon.addHistogram( new TH1D(icat+"tauinclusivecharge",     ";p_{T}^{#tau};Events",  5, -2.,    2.  ));
      mon.addHistogram( new TH1D(icat+"tauinclusivedz",         ";dz^{#tau};Events",    50,  0.,   10.  ));
      mon.addHistogram( new TH1D(icat+"tauinclusivevz",         ";vz^{#tau};Events",    50,  0.,   10.  ));
      mon.addHistogram( new TH1D(icat+"tauinclusiveemfraction", ";emf^{#tau};Events",   50,  0.,    5.  ));
      mon.addHistogram( new TH1D(icat+"tauinclusivedizeta"    , ";dZ^{#tau};Events",    50,  0.,   10.  ));

      //lepton control
      mon.addHistogram( new TH1D(icat+"inclusivept",      ";Transverse momentum [GeV];Events",             50, 0.,  500.  ));
      mon.addHistogram( new TH1D(icat+"leadleptonpt",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  )); 
      mon.addHistogram( new TH1D(icat+"leadleptoneta",    ";Pseudo-rapidity;Events",                       50, 0.,    2.5 ));
      mon.addHistogram( new TH1D(icat+"trailerleptonpt",  ";Transverse momentum [GeV];Events",             50, 0.,  500.  ));
      mon.addHistogram( new TH1D(icat+"trailerleptoneta", ";Pseudo-rapidity;Events",                       50, 0.,    2.5 ));
      mon.addHistogram( new TH1D(icat+"pte",              ";Electron transverse momentum [GeV];Events",    50, 0.,  500.  ));
      mon.addHistogram( new TH1D(icat+"ptmu",             ";Muon transverse momentum [GeV];Events",        50, 0.,  500.  ));
      mon.addHistogram( new TH1D(icat+"qt",               ";Transverse momentum [GeV];Events / (1 GeV)", 1500, 0., 1500.  ));
      mon.addHistogram( new TH1D(icat+"emva", "; e-id MVA; Electrons", 50, 0.95,1.0) );
      
      // Dilepton control
      mon.addHistogram( new TH1D(icat+"sumptll",      ";Sum of lepton transverse momenta [GeV];Events",                     75, 0.,  750.  ));
      mon.addHistogram( new TH1D(icat+"mll",          ";Dilepton invariant mass [GeV];Events",                              50, 0.,  500.  ));
      mon.addHistogram( new TH1D(icat+"ptll",         ";Dilepton transverse momentum [GeV];Events",                         75, 0.,  750.  ));
      mon.addHistogram( new TH1D(icat+"yll",          ";Rapidity;Events",                                                   50, 0.,    3.  ));
      mon.addHistogram( new TH1D(icat+"dilarccosine", ";#theta(l,l') [rad];Events",                                         64, 0.,    3.2 ));

      mon.addHistogram( new TH1D(icat+"mtsum",        ";M_{T}(l^{1},E_{T}^{miss})+M_{T}(l^{2},E_{T}^{miss}) [GeV];Events", 100, 0.,   10.  ));
      mon.addHistogram( new TH1D(icat+"ht",           ";H_{T} [GeV];Events",                                                75, 0., 1500.  ));
      mon.addHistogram( new TH1D(icat+"htb",          ";H_{T} (bjets) [GeV];Events",                                        75, 0., 1500.  ));
      mon.addHistogram( new TH1D(icat+"htnol",        "; H_[T] (no leptons) [GeV];Events",                                  75, 0., 1500.  ));
      mon.addHistogram( new TH1D(icat+"htbnol",       "; H_[T] (bjets, no leptons) [GeV];Events",                           75, 0., 1500.  ));

      // add plots for extra leptons in the event
      // add plots for third lepton pt etc

      mon.addHistogram (new TH1D(icat+"csv",           ";Combined Secondary Vertex;Jets",    50, 0.,    1. ));
      mon.addHistogram (new TH1D(icat+"csvb",          ";Combined Secondary Vertex;Jets",    50, 0.,    1. ));
      mon.addHistogram (new TH1D(icat+"csvc",          ";Combined Secondary Vertex;Jets",    50, 0.,    1. ));
      mon.addHistogram (new TH1D(icat+"csvothers",     ";Combined Secondary Vertex;Jets",    50, 0.,    1. ));

      mon.addHistogram (new TH1D(icat+"leadjetpt",     ";Transverse momentum [GeV];Events", 100, 0., 1000. ));
      mon.addHistogram (new TH1D(icat+"trailerjetpt",  ";Transverse momentum [GeV];Events", 100, 0., 1000. ));
      mon.addHistogram (new TH1D(icat+"fwdjeteta",     ";Pseudo-rapidity;Events",            60, 0.,    3. ));
      mon.addHistogram (new TH1D(icat+"leadjeteta",    ";Pseudo-rapidity;Events",            60, 0.,    3. ));
      mon.addHistogram (new TH1D(icat+"trailerjeteta", ";Pseudo-rapidity;Events",            60, 0.,    3. ));
      mon.addHistogram (new TH1D(icat+"cenjeteta",     ";Pseudo-rapidity;Events",            60, 0.,    3. ));
      
      TH1* hbtags   = mon.addHistogram(new TH1D(icat+"nbtags",   ";b-tag multiplicity;Events", 5, 0., 5. ));
      TH1* hjets    = mon.addHistogram(new TH1D(icat+"njets",    ";Jet multiplicity;Events",   6, 0., 6. ));
      for (int ibin = 1; ibin <= hjets->GetXaxis ()->GetNbins (); ibin++)
        {
          TString label ("");
          if (ibin == h->GetXaxis()->GetNbins() || (TString(h->GetName()).Contains("btags") && ibin == h->GetXaxis()->GetNbins()-1 ) )
            label += "#geq";
          else
            label += "=";
          label += (ibin - 1);
          hjets   ->GetXaxis()->SetBinLabel(ibin, label);
          hbtags  ->GetXaxis()->SetBinLabel(ibin, label);
        }
     
      mon.addHistogram (new TH1D (icat + "mindphijmet",    ";min #Delta#phi(jet,E_{T}^{miss});Events",     40,    0.,    4.  ));
      mon.addHistogram (new TH1D (icat + "mindphijmetNM1", ";min #Delta#phi(jet,E_{T}^{miss});Events",     40,    0.,    4.  ));
      mon.addHistogram (new TH1D (icat + "balance",        ";E_{T}^{miss}/q_{T};Events",                   25,    0.,    2.5 ));
      mon.addHistogram (new TH1D (icat + "balanceNM1",     ";E_{T}^{miss}/q_{T};Events",                   25,    0.,    2.5 ));
      mon.addHistogram (new TH1D (icat + "axialmet",       ";Axial missing transvere energy [GeV];Events", 50, -100.,  400.  ));
      mon.addHistogram (new TH1D (icat + "axialmetNM1",    ";Axial missing transvere energy [GeV];Events", 50, -100.,  400.  ));
      mon.addHistogram (new TH1D (icat + "met",            ";Missing transverse energy [GeV];Events",      50,    0., 1000.  ));
      mon.addHistogram (new TH1D (icat + "mt",             ";Transverse mass;Events",                      50,    0.,  500.  ));
      mon.addHistogram (new TH1D (icat + "mtresponse",     ";Transverse mass response;Events",            100,    0.,    2.  ));
      mon.addHistogram (new TH1D (icat + "mtcheckpoint",   ";Transverse mass [GeV];Events",               160,  150., 1750.  ));
      mon.addHistogram (new TH1D (icat + "metcheckpoint",  ";Missing transverse energy [GeV];Events",     100,    0.,  500.  ));

    } // End of loop on controlCats


  //
  // STATISTICAL ANALYSIS
  //
  TH1D* Hoptim_systs = (TH1D*) mon.addHistogram (new TH1D ("optim_systs", ";syst;", nSystVars, 0, nSystVars));
  for (size_t ivar=0; ivar<nSystVars; ++ivar) Hoptim_systs->GetXaxis()->SetBinLabel(ivar+1, systVars[ivar]);

  // Final distributions to compute systematics on
  for(size_t ivar=0; ivar<nSystVars; ++ivar) 
    {
      TString var=systVars[ivar];
      
      // dilepton or both
      mon.addHistogram(new TH1D("finalnbjets"         +var, ";b-jet multiplicity;Events",     6, 0.,   6. ));        
      mon.addHistogram(new TH1D("finalmt"             +var, ";Transverse mass [GeV];Events", 50, 0., 500. ));
      
      // lepton-tau
      mon.addHistogram(new TH1D("finaltaur"           +var, ";R^{#tau};Events",                    10,  0.,   1.   ));
      mon.addHistogram(new TH1D("finaltaupolarization"+var, ";Y^{#tau};Events",                    40, -1.,   3.   ));
      mon.addHistogram(new TH1D("finaldphilepmet"     +var, ";#Delta#phi(#tau_{h}-#it{l});Events", 60,  0.  , 3.15 ));
      mon.addHistogram(new TH1D("finaldphitaumet"     +var, ";#Delta#phi(#tau_{h}-MET);Events",    60,  0.,   3.15 ));
      mon.addHistogram(new TH1D("finaldphileptau"     +var, ";#Delta#phi(#it{l}-#tau_{h});Events", 60,  0.,   3.15 ));
      mon.addHistogram(new TH1D("finaltaupt"          +var, ";p_{T}^{#tau} [GeV];Events",          50,  0., 500.   ));
      mon.addHistogram(new TH1D("finalmutaumass"      +var, ";M_{#mu#tau_{h}} [GeV];Events",       50,  0., 500.   ));

    }
  
  //##############################################
  //######## GET READY FOR THE EVENT LOOP ########
  //##############################################
  size_t totalEntries(0);

  TFile* summaryFile = NULL;
  TTree* summaryTree = NULL; //ev->;
  //  
  //  if(saveSummaryTree)
  //    {
  //      TDirectory* cwd = gDirectory;
  //      std::string summaryFileName(outUrl); 
  //      summaryFileName.replace(summaryFileName.find(".root", 0), 5, "_summary.root");
  //      
  //      summaryFile = new TFile(summaryFileName.c_str() "recreate");
  //      
  //      summaryTree = new TTree("Events", "Events");
  //    KEY: TTreeMetaData;1
  //                         KEY: TTreeParameterSets;1
  //                                                   KEY: TTreeParentage;1
  //                                                                         KEY: TTreeEvents;1
  //                                                                                            KEY: TTreeLuminosityBlocks;1
  //                                                                                                                         KEY: TTreeRuns;
  //      summaryTree->SetDirectory(summaryFile);  // This line is probably not needed
  //      
  //      summmaryTree->Branch(
  //
  //      cwd->cd();
  //    }
  //
  
  
  //MC normalization (to 1/pb)
  if(debug) cout << "DEBUG: xsec: " << xsec << endl;



  // MET Correction level
  pat::MET::METCorrectionLevel metcor = pat::MET::METCorrectionLevel::Type1XY;

  //jet energy scale and uncertainties 
  TString jecDir = runProcess.getParameter < std::string > ("jecDir");
  gSystem->ExpandPathName (jecDir);
  FactorizedJetCorrector *jesCor = utils::cmssw::getJetCorrector(jecDir, isMC);
  TString pf(isMC ? "MC" : "DATA");
  JetCorrectionUncertainty *totalJESUnc = new JetCorrectionUncertainty ((jecDir + "/"+pf+"_Uncertainty_AK4PFchs.txt").Data ());
  
  //muon energy scale and uncertainties
  TString muscleDir = runProcess.getParameter<std::string>("muscleDir");
  gSystem->ExpandPathName(muscleDir);
  rochcor2015* muCor = new rochcor2015(); // This replaces the old MusScleFitCorrector that was used at RunI

  // Electron energy scale, based on https://twiki.cern.ch/twiki/bin/viewauth/CMS/EGMSmearer and adapted to this framework
  string EGammaEnergyCorrectionFile = "EgammaAnalysis/ElectronTools/data/76X_16DecRereco_2015";
  EpCombinationTool theEpCombinationTool;
  theEpCombinationTool.init((string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/weights/GBRForest_data_25ns.root").c_str(), "gedelectron_p4combination_25ns");  //got confirmation from Matteo Sani that this works for both data and MC 
  ElectronEnergyCalibratorRun2 ElectronEnCorrector(theEpCombinationTool, isMC, false, EGammaEnergyCorrectionFile);
  ElectronEnCorrector.initPrivateRng(new TRandom(1234));
  
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
  //sfb is not actually used as it's taken from btagCal
  // beff = 0.747; sfb = 0.899; //for Loose WP
  sfb = 0.861;
  // sbbunc =;
  beff = 0.559;
  
  // Setup calibration readers
  BTagCalibration btagCalib("CSVv2", string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/weights/btagSF_CSVv2.csv");
  BTagCalibrationReader btagCal   (&btagCalib, BTagEntry::OP_LOOSE, "mujets", "central");  // calibration instance, operating point, measurement type, systematics type
  BTagCalibrationReader btagCalUp (&btagCalib, BTagEntry::OP_LOOSE, "mujets", "up"     );  // sys up
  BTagCalibrationReader btagCalDn (&btagCalib, BTagEntry::OP_LOOSE, "mujets", "down"   );  // sys down
  BTagCalibrationReader btagCalL  (&btagCalib, BTagEntry::OP_LOOSE, "comb", "central");  // calibration instance, operating point, measurement type, systematics type
  BTagCalibrationReader btagCalLUp(&btagCalib, BTagEntry::OP_LOOSE, "comb", "up"     );  // sys up
  BTagCalibrationReader btagCalLDn(&btagCalib, BTagEntry::OP_LOOSE, "comb", "down"   );  // sys down




  TString
    electronIdMainTag("cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
    electronIdVetoTag("cutBasedElectronID-Spring15-25ns-V1-standalone-tight");

  //pileup weighting
  edm::LumiReWeighting * LumiWeights = NULL;
  utils::cmssw::PuShifter_t PuShifters;
  double PUNorm[] = { 1, 1, 1 };
  double xsecWeight(1.0);
  if (isMC)
    {
      std::vector<double> dataPileupDistributionDouble = runProcess.getParameter < std::vector < double >>("datapileup");
      std::vector<float> dataPileupDistribution;
      for (unsigned int i = 0; i < dataPileupDistributionDouble.size (); i++)
        {
          dataPileupDistribution.push_back (dataPileupDistributionDouble[i]);
        }
      std::vector<float> mcPileupDistribution;
      double totalNumEvent = utils::getMCPileupDistributionAndTotalEventFromMiniAOD(urls,dataPileupDistribution.size(), mcPileupDistribution);
      xsecWeight=xsec/totalNumEvent;
      //utils::getMCPileupDistributionFromMiniAOD(urls, dataPileupDistribution.size (), mcPileupDistribution);
      while(mcPileupDistribution.size() < dataPileupDistribution.size()) mcPileupDistribution.push_back(0.0);
      while(mcPileupDistribution.size() > dataPileupDistribution.size()) dataPileupDistribution.push_back(0.0);
      gROOT->cd ();             //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE
      LumiWeights = new edm::LumiReWeighting(mcPileupDistribution, dataPileupDistribution);
      PuShifters = utils::cmssw::getPUshifters(dataPileupDistribution, 0.05);
      utils::getPileupNormalization(mcPileupDistribution, PUNorm, LumiWeights, PuShifters);
    }
  
  gROOT->cd ();                 //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE
  
  //higgs::utils::EventCategory eventCategoryInst(higgs::utils::EventCategory::EXCLUSIVE2JETSVBF); //jet(0,>=1)+vbf binning
 

  patUtils::MetFilter metFiler;
  if(!isMC) { 
        metFiler.FillBadEvents(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/MetFilter/DoubleEG_RunD/DoubleEG_csc2015.txt");
        metFiler.FillBadEvents(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/MetFilter/DoubleEG_RunD/DoubleEG_ecalscn1043093.txt");
        metFiler.FillBadEvents(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/MetFilter/DoubleMuon_RunD/DoubleMuon_csc2015.txt");
        metFiler.FillBadEvents(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/MetFilter/DoubleMuon_RunD/DoubleMuon_ecalscn1043093.txt");
        metFiler.FillBadEvents(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/MetFilter/MuonEG_RunD/MuonEG_csc2015.txt");
        metFiler.FillBadEvents(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/MetFilter/MuonEG_RunD/MuonEG_ecalscn1043093.txt");
  }


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

        edm::EventBase const & myEvent = ev;
        // Take into account the negative weights from some NLO generators (otherwise some phase space will be double counted)
        double weightGen(1.);
        if(isNLOMC)
          {
            //double weightGen(0.);
            //double weightLhe(0.);
            
            fwlite::Handle<GenEventInfoProduct> evt;
            evt.getByLabel(ev, "generator");
            if(evt.isValid())
              {
                weightGen = (evt->weight() > 0 ) ? 1. : -1. ;
              }
            
            // FIXME: this is for PDF uncertainties, must reactivate it at some point.
            //fwlite::Handle<LHEEventProduct> lheEvtProd;
            //lheEvtProd.getByLabel(ev, "externalLHEProducer");
            //if(lheEvtProd.isValid())
            //  {
            //    weightLhe=lheEvtProd->originalXWGTUP();
            //    
            //   //for(unsigned int i=0; i<evet->weights().size();i++){
            //   //  double asdde=evet->weights()[i].wgt;
            //   //  EventInfo.ttbar_w[EventInfo.ttbar_nw]=EventInfo.ttbar_w[0]*asdde/asdd;
            //   //  EventInfo.ttbar_nw++;
            //   //}
            //  }
            //cout << "Event " << iev << " has genweight: " << weightGen << " and LHE weight " << weightLhe << endl;
            
          }
        
        
        std::vector < TString > tags (1, "all"); // Inclusive inclusiveness
        
        //
        // DERIVE WEIGHTS TO APPLY TO SAMPLE
        //
        
        //pileup weight
        double weight           (xsecWeight);
        double rawWeight        (1.0);
        double TotalWeight_plus (1.0);
        double TotalWeight_minus(1.0);
        double puWeight         (1.0);
        
        
        // This must remain deactivated if you use HT-binned samples (it was for pthat-binned samples)
        // if (isV0JetsMC)
        //   {
        //     fwlite::Handle < LHEEventProduct > lheEPHandle;
        //     lheEPHandle.getByLabel (ev, "externalLHEProducer");
        //     mon.fillHisto ("nup", "", lheEPHandle->hepeup ().NUP, 1);
        //     if (lheEPHandle->hepeup ().NUP > 5)  continue;
        //     mon.fillHisto ("nupfilt", "", lheEPHandle->hepeup ().NUP, 1);
        //   }
        
        // HT-binned samples stitching: https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2015#MC_and_data_samples
        if(isV0JetsMC)
          {
            // access generator level HT               
            fwlite::Handle<LHEEventProduct> lheEventProduct;
            lheEventProduct.getByLabel(ev, "externalLHEProducer");
            //edm::Handle<LHEEventProduct> lheEventProduct;
            //ev.getByLabel( 'externalLHEProducer', lheEventProduct);
            const lhef::HEPEUP& lheEvent = lheEventProduct->hepeup(); 
            std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;
            double lheHt = 0.;
            size_t numParticles = lheParticles.size();
            for ( size_t idxParticle = 0; idxParticle < numParticles; ++idxParticle ) {
              int absPdgId = TMath::Abs(lheEvent.IDUP[idxParticle]);
              int status = lheEvent.ISTUP[idxParticle];
              if ( status == 1 && ((absPdgId >= 1 && absPdgId <= 6) || absPdgId == 21) ) { // quarks and gluons
                lheHt += TMath::Sqrt(TMath::Power(lheParticles[idxParticle][0], 2.) + TMath::Power(lheParticles[idxParticle][1], 2.)); // first entry is px, second py
              }                                        
            }
            if(debug) cout << "Sample: " << dtag << ", lheHt: " << lheHt << ", scale factor from spreadsheet: " << patUtils::getHTScaleFactor(dtag, lheHt) << endl;
            weightGen *=   patUtils::getHTScaleFactor(dtag, lheHt);
          }         
        
        weight *= weightGen;
        rawWeight *=weightGen;
        
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

        // Apply pileup reweighting
        if(isMC)
          {
            int ngenITpu = 0;
            fwlite::Handle < std::vector < PileupSummaryInfo > >puInfoH;
            puInfoH.getByLabel (ev, "slimmedAddPileupInfo");
            for (std::vector < PileupSummaryInfo >::const_iterator it = puInfoH->begin (); it != puInfoH->end (); it++)
              {
                if (it->getBunchCrossing () == 0) ngenITpu += it->getTrueNumInteractions(); //it->getPU_NumInteractions ();
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
        
        //##############################################   EVENT LOOP STARTS   ##############################################

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
        bool eTrigger    (
                          isMC ? 
                          utils::passTriggerPatterns (tr, "HLT_Ele27_eta2p1_WP75_Gsf_v*")
                          :
                          utils::passTriggerPatterns (tr, "HLT_Ele27_eta2p1_WPLoose_Gsf_v*")
                          );
        bool muTrigger   (
                          utils::passTriggerPatterns (tr, "HLT_IsoMu20_v*", "HLT_IsoTkMu20_v*")
                          );

        if(!isMC && muTrigger) mon.fillHisto("nvtx_singlemu_pileup", tags, nGoodPV, 1.);
        if(!isMC && eTrigger)  mon.fillHisto("nvtx_singlee_pileup",  tags, nGoodPV, 1.);
        
        if(filterOnlySINGLEMU) {                    eTrigger = false; }
        if(filterOnlySINGLEE)  { muTrigger = false;                   }
        
        if (!(eTrigger || muTrigger)) continue;         //ONLY RUN ON THE EVENTS THAT PASS OUR TRIGGERS

        // ------------ Apply MET filters ------------
        int metFilterValue(metFiler.passMetFilterInt(ev));
        mon.fillHisto("metFilter_eventflow", "", metFilterValue, weight);
        if(metFilterValue!=0) continue; // MET Filters must be applied to both data and MC

        //load all the objects we will need to access

        double rho = 0;
        fwlite::Handle<double> rhoHandle;
        rhoHandle.getByLabel(ev, "fixedGridRhoFastjetAll");
        if(rhoHandle.isValid() ) rho = *rhoHandle;
        
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
        pat::MET met = mets[0];
        
        pat::TauCollection taus;
        fwlite::Handle<pat::TauCollection> tausHandle;
        tausHandle.getByLabel(ev, "slimmedTaus");
        if(tausHandle.isValid() ) taus = *tausHandle;
        

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

        LorentzVector
          muDiff(0., 0., 0., 0.),
          elDiff(0., 0., 0., 0.);
        std::vector<patUtils::GenericLepton> selLeptons;
        unsigned int nVetoE(0), nVetoMu(0);
        for(size_t ilep=0; ilep<leptons.size (); ++ilep)
          {
            patUtils::GenericLepton& lepton = leptons[ilep];

            bool 
              passKin(true),     passId(true),     passIso(true),
              passVetoKin(true), passVetoId(true), passVetoIso(true);
            
            int lid(lepton.pdgId());
            
            // Apply muon corrections
            if(abs(lid) == 13)
            {
              if(muCor)
                {
                  float qter;
                  TLorentzVector p4(leptons[ilep].px(),leptons[ilep].py(),leptons[ilep].pz(),leptons[ilep].energy());
                  if(isMC)
                    {
                      muCor->momcor_mc  (p4, lid<0 ? -1 :1, 0, qter);
                    }
                  else
                    {
                      muCor->momcor_data(p4, lid<0 ? -1 :1, 0, qter); 
                    }
                  muDiff -= leptons[ilep].p4();
                  leptons[ilep].setP4(LorentzVector(p4.Px(),p4.Py(),p4.Pz(),p4.E() ) );
                  muDiff += leptons[ilep].p4();
                }
            }
            
            // Apply electron corrections
            if(abs(lid)==11)
              {
                elDiff -= leptons[ilep].p4();                   
                ElectronEnCorrector.calibrate(leptons[ilep].el, ev.eventAuxiliary().run(), edm::StreamID::invalidStreamID()); 
                leptons[ilep] = patUtils::GenericLepton(leptons[ilep].el); //recreate the generic lepton to be sure that the p4 is ok
                elDiff += leptons[ilep].p4();                 
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
          if(lepton.pt() < 25.)                      passKin = false;
          if(leta > 2.1)                                    passKin = false;
          if(lid == 11 && (leta > 1.4442 && leta < 1.5660)) passKin = false; // Crack veto
          
          // Veto leptons kin
          if (lepton.pt () < 20)                      passVetoKin = false;
          if (leta > 2.1)                                    passVetoKin = false;
          if (lid == 11 && (leta > 1.4442 && leta < 1.5660)) passVetoKin = false; // Crack veto
          
          //Cut based identification 
          
          //std::vector<pat::Electron> dummyShit; dummyShit.push_back(leptons[ilep].el);
          
          passId      = lid==11 ? patUtils::passId(leptons[ilep].el, vtx[0], patUtils::llvvElecId::Tight, patUtils::CutVersion::ICHEP16Cut) : patUtils::passId(leptons[ilep].mu, vtx[0], patUtils::llvvMuonId::Tight, patUtils::CutVersion::ICHEP16Cut);
          passVetoId = lid==11 ? patUtils::passId(leptons[ilep].el, vtx[0], patUtils::llvvElecId::Loose, patUtils::CutVersion::ICHEP16Cut) : patUtils::passId(leptons[ilep].mu, vtx[0], patUtils::llvvMuonId::Loose, patUtils::CutVersion::ICHEP16Cut);

          //passId     = lid == 11 ? patUtils::passId(electronVidMainId, myEvent, lepton.el) : patUtils::passId(lepton.mu, goodPV, patUtils::llvvMuonId::StdTight);
          //passVetoId = lid == 11 ? patUtils::passId(electronVidVetoId, myEvent, lepton.el) : patUtils::passId(lepton.mu, goodPV, patUtils::llvvMuonId::StdLoose);

          //isolation
          passIso     = lid == 11 ? patUtils::passIso(leptons[ilep].el, patUtils::llvvElecIso::Tight, patUtils::CutVersion::ICHEP16Cut, 0.) : patUtils::passIso(lepton.mu, patUtils::llvvMuonIso::Tight, patUtils::CutVersion::ICHEP16Cut);
          passVetoIso = lid == 11 ? patUtils::passIso(leptons[ilep].el, patUtils::llvvElecIso::Loose, patUtils::CutVersion::ICHEP16Cut, 0.) : patUtils::passIso(lepton.mu, patUtils::llvvMuonIso::Loose, patUtils::CutVersion::ICHEP16Cut);

          if     (passKin     && passId     && passIso)     selLeptons.push_back(lepton);
          else if(passVetoKin && passVetoId && passVetoIso) lid==11 ? nVetoE++ : nVetoMu++;
          
        }
      std::sort(selLeptons.begin(),   selLeptons.end(),   utils::sort_CandidatesByPt);


      // Propagate lepton energy scale to MET
      met.setP4(met.p4() - muDiff - elDiff); //note this also propagates to all MET uncertainties
      met.setUncShift(met.px() - muDiff.px()*0.01, met.py() - muDiff.py()*0.01, met.sumEt() - muDiff.pt()*0.01, pat::MET::METUncertainty::MuonEnUp);   //assume 1% uncertainty on muon rochester
      met.setUncShift(met.px() + muDiff.px()*0.01, met.py() + muDiff.py()*0.01, met.sumEt() + muDiff.pt()*0.01, pat::MET::METUncertainty::MuonEnDown); //assume 1% uncertainty on muon rochester
      met.setUncShift(met.px() - elDiff.px()*0.01, met.py() - elDiff.py()*0.01, met.sumEt() - elDiff.pt()*0.01, pat::MET::METUncertainty::ElectronEnUp);   //assume 1% uncertainty on electron scale correction
      met.setUncShift(met.px() + elDiff.px()*0.01, met.py() + elDiff.py()*0.01, met.sumEt() + elDiff.pt()*0.01, pat::MET::METUncertainty::ElectronEnDown); //assume 1% uncertainty on electron scale correction
      
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
      
      //
      //JET/MET ANALYSIS
      //
      if(debug) cout << "Now update Jet Energy Corrections" << endl;
      //add scale/resolution uncertainties and propagate to the MET      
      utils::cmssw::updateJEC(jets,jesCor,totalJESUnc,rho,nGoodPV,isMC);

      /// if(debug) cout << "Update also MET" << endl;
      /// std::vector<LorentzVector> newMet=utils::cmssw::getMETvariations(met/*recoMet*/,jets,selLeptons,isMC); // FIXME: Must choose a lepton collection. Perhaps loose leptons?
      /// met=newMet[utils::cmssw::METvariations::NOMINAL];
      if(debug) cout << "Jet Energy Corrections updated" << endl;
      
      std::map<string, pat::JetCollection> selJetsVar;
      std::map<string, int   > njetsVar;
      std::map<string, int   > nbtagsVar;
      std::map<string, double> mindphijmetVar;
      for(unsigned int ivar=0;ivar<jetVarNames.size();ivar++){njetsVar[jetVarNames[ivar]] = 0;}  //initialize
      for(unsigned int ivar=0;ivar<jetVarNames.size();ivar++){mindphijmetVar[jetVarNames[ivar]] = 9999.0;}  //initialize
      nbtagsVar[""] = 0; nbtagsVar["_eff_bup"] = 0; nbtagsVar["_eff_bdown"] = 0;  //initialize
      
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
          
          for (size_t ilep = 0; ilep < selLeptons.size(); ilep++)
            minDRljSingleLep = TMath::Min(minDRljSingleLep, reco::deltaR (jet, selLeptons[ilep]));
          
          //jet id
          bool passPFloose = passPFJetID("Loose", jet); 
          if (!passPFloose || jet.pt() <30. || fabs(jet.eta()) > 2.5) continue;
          if (minDRlj < 0.4) continue;
          
          selJets.push_back(jet);
          
          double dphijmet = fabs (deltaPhi (met.phi(), jet.phi()));
          if (dphijmet < mindphijmet) mindphijmet = dphijmet;
          bool hasCSVtag(jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > btagMedium);
          bool hasCSVtagUp(hasCSVtag);
          bool hasCSVtagDown(hasCSVtag);
          
          //update according to the SF measured by BTV
          if(isMC){
            int flavId=jet.partonFlavour();  double eta=jet.eta();
            if      (abs(flavId)==5){  btsfutil.modifyBTagsWithSF(hasCSVtag    , btagCal   .eval(BTagEntry::FLAV_B   , eta, jet.pt()), beff);
              btsfutil.modifyBTagsWithSF(hasCSVtagUp  , btagCalUp .eval(BTagEntry::FLAV_B   , eta, jet.pt()), beff);
              btsfutil.modifyBTagsWithSF(hasCSVtagDown, btagCalDn .eval(BTagEntry::FLAV_B   , eta, jet.pt()), beff);
            }else if(abs(flavId)==4){  btsfutil.modifyBTagsWithSF(hasCSVtag    , btagCal   .eval(BTagEntry::FLAV_C   , eta, jet.pt()), beff);
              btsfutil.modifyBTagsWithSF(hasCSVtagUp  , btagCalUp .eval(BTagEntry::FLAV_C   , eta, jet.pt()), beff);
              btsfutil.modifyBTagsWithSF(hasCSVtagDown, btagCalDn .eval(BTagEntry::FLAV_C   , eta, jet.pt()), beff);
            }else{                     btsfutil.modifyBTagsWithSF(hasCSVtag    , btagCalL  .eval(BTagEntry::FLAV_UDSG, eta, jet.pt()), leff);
              btsfutil.modifyBTagsWithSF(hasCSVtagUp  , btagCalLUp.eval(BTagEntry::FLAV_UDSG, eta, jet.pt()), leff);
              btsfutil.modifyBTagsWithSF(hasCSVtagDown, btagCalLDn.eval(BTagEntry::FLAV_UDSG, eta, jet.pt()), leff);
            }
          }
          
          if(hasCSVtag    )nbtagsVar[""          ]++;
          if(hasCSVtagUp  )nbtagsVar["_eff_bup"  ]++;
          if(hasCSVtagDown)nbtagsVar["_eff_bdown"]++;
          
          for(unsigned int ivar=0;ivar<jetVarNames.size();ivar++){
            pat::Jet varJet = jet;
            if(ivar!=0) varJet.setP4(jet.p4() * jet.userFloat(jetVarNames[ivar]));
            selJetsVar[jetVarNames[ivar]].push_back(varJet);
            
            if(varJet.pt()>30){
              njetsVar[jetVarNames[ivar]]++;
              
              float dphijmet=fabs(deltaPhi(met.corP4(metcor).phi(), varJet.phi()));
              if(dphijmet<mindphijmetVar[jetVarNames[ivar]]) mindphijmetVar[jetVarNames[ivar]]=dphijmet;
            }
          }
          
          if(!hasCSVtag) continue;
          selBJets.push_back(jet);
        }
      
      //sort all jet collection by pT
      for(auto jetCollIt = selJetsVar.begin(); jetCollIt!=selJetsVar.end(); jetCollIt++){
        std::sort(jetCollIt->second.begin(), jetCollIt->second.end(), utils::sort_CandidatesByPt);
      }
      
      std::sort (selJets.begin(),  selJets.end(),  utils::sort_CandidatesByPt);
      std::sort (selBJets.begin(), selBJets.end(), utils::sort_CandidatesByPt);
      
      //
      // ASSIGN CHANNEL
      //
      std::vector < TString > chTags; chTags.clear();
      int 
        dilId (1),
        slepId(0);
      LorentzVector dileptonSystem (0, 0, 0, 0);
      if(selLeptons.size()>=2)
        {
          for (size_t ilep = 0; ilep < 2; ilep++)
            {
              dilId *= selLeptons[ilep].pdgId();
              int id(abs (selLeptons[ilep].pdgId()));
              dileptonSystem += selLeptons[ilep].p4();
            }
        }
      
      if(selLeptons.size()>0)
        slepId=selLeptons[0].pdgId();
      
      // Event classification. Single lepton triggers are used for offline selection of dilepton events. The "else if"s guarantee orthogonality
      bool 
        isSingleMu(false),
        isSingleE(false),
        isDoubleMu(false),
        isDoubleE(false),
        isEMu(false);
      int multiChannel(0);
      if      (abs(slepId)==13 && muTrigger && nVetoE==0 && nVetoMu==0){ isSingleMu = true; multiChannel++; chTags.push_back("singlemu");}
      else if (abs(slepId)==11 && eTrigger  && nVetoE==0 && nVetoMu==0){ isSingleE  = true; multiChannel++; chTags.push_back("singlee");}
      else if (abs(dilId)==121 && eTrigger                            ){ isDoubleE  = true; multiChannel++; chTags.push_back("ee");}
      else if (abs(dilId)==169 && muTrigger                           ){ isDoubleMu = true; multiChannel++; chTags.push_back("mumu");}
      else if (abs(dilId)==143 && muTrigger                           ){ isEMu      = true; multiChannel++; chTags.push_back("emu");}
      else if (abs(dilId)==143 && eTrigger  && !muTrigger             ){ isEMu      = true; multiChannel++; chTags.push_back("emu");} // Pick up the largest number of emu events possible, maintaining orthogonality
      
      // keep in mind the eventCategory thingy for more refined categorization // TString evCat=eventCategoryInst.GetCategory(selJets,dileptonSystem);
      //std::vector < TString > tags (1, "all");
      for (size_t ich = 0; ich < chTags.size(); ich++)
        {
          tags.push_back (chTags[ich]);
          //tags.push_back( chTags[ich]+evCat );
        }
      if(multiChannel>1) nMultiChannel++;

      // Dilepton full analysis
      if( isDoubleE || isEMu || isDoubleMu){
        
        mon.fillHisto("nvtx_pileup", tags, nGoodPV, weight);
        
        if(selLeptons.size()<2 || nGoodPV == 0) continue; // Save time
        // Apply lepton efficiencies
        //for(size_t ilep=0; ilep<2; ++ilep){
        //  int id (abs (selLeptons[ilep].pdgId()));
        //  weight *= isMC ? lepEff.getLeptonEfficiency(selLeptons[ilep].pt(), selLeptons[ilep].eta(), id, id == 11 ? "loose" : "loose").first : 1.0;
        //}
        
        // Event selection booleans
        bool passMllVeto(isEMu ? dileptonSystem.mass()>12. : (fabs(dileptonSystem.mass()-91.)>15 && dileptonSystem.mass()>12. ) );
        bool passJetSelection(selJets.size()>1);
        bool passMetSelection(met.pt()>40.);
        bool passOS(selLeptons[0].pdgId() * selLeptons[1].pdgId() < 0 );
        bool passBtagsSelection(selBJets.size()>1); // FIXME: differentiate chhiggs selection from cross-section selection
       
        // Setting up control categories and fill up event flow histo
        std::vector < TString > ctrlCats; ctrlCats.clear ();
                                                                                                 { ctrlCats.push_back("step1"); mon.fillHisto("xseceventflowdilep", tags, 0, weight); mon.fillHisto("chhiggseventflowdilep", tags, 0, weight); }
        if(passMllVeto   )                                                                       { ctrlCats.push_back("step2"); mon.fillHisto("xseceventflowdilep", tags, 1, weight); mon.fillHisto("chhiggseventflowdilep", tags, 1, weight); }
        if(passMllVeto && passJetSelection )                                                     { ctrlCats.push_back("step3"); mon.fillHisto("xseceventflowdilep", tags, 2, weight); mon.fillHisto("chhiggseventflowdilep", tags, 2, weight); }
        if(passMllVeto && passJetSelection && passMetSelection )                                 { ctrlCats.push_back("step4"); mon.fillHisto("xseceventflowdilep", tags, 3, weight); mon.fillHisto("chhiggseventflowdilep", tags, 3, weight); }
        if(passMllVeto && passJetSelection && passMetSelection && passOS )                       { ctrlCats.push_back("step5"); mon.fillHisto("xseceventflowdilep", tags, 4, weight); mon.fillHisto("chhiggseventflowdilep", tags, 4, weight); }
        if(passMllVeto && passJetSelection && passMetSelection && passOS && passBtagsSelection ) { ctrlCats.push_back("step6"); mon.fillHisto("xseceventflowdilep", tags, 5, weight); mon.fillHisto("chhiggseventflowdilep", tags, 5, weight); }
        

        bool passBtagsSelection_0(selBJets.size()==0);
        bool passBtagsSelection_1(selBJets.size()==1);
        bool passBtagsSelection_2(selBJets.size()>1);

                                                                                                   { ctrlCats.push_back("altstep1"); mon.fillHisto("xsecalteventflowdilep", tags, 0, weight); mon.fillHisto("chhiggsalteventflowdilep", tags, 0, weight);}
        if(passMllVeto   )                                                                         { ctrlCats.push_back("altstep2"); mon.fillHisto("xsecalteventflowdilep", tags, 1, weight); mon.fillHisto("chhiggsalteventflowdilep", tags, 1, weight);}
        if(passMllVeto && passJetSelection )                                                       { ctrlCats.push_back("altstep3"); mon.fillHisto("xsecalteventflowdilep", tags, 2, weight); mon.fillHisto("chhiggsalteventflowdilep", tags, 2, weight);}
        if(passMllVeto && passJetSelection && passMetSelection )                                   { ctrlCats.push_back("altstep4"); mon.fillHisto("xsecalteventflowdilep", tags, 3, weight); mon.fillHisto("chhiggsalteventflowdilep", tags, 3, weight);}
        if(passMllVeto && passJetSelection && passMetSelection && passOS )                         { ctrlCats.push_back("altstep5"); mon.fillHisto("xsecalteventflowdilep", tags, 4, weight); mon.fillHisto("chhiggsalteventflowdilep", tags, 4, weight);}
        if(passMllVeto && passJetSelection && passMetSelection && passOS && passBtagsSelection_0 ) { ctrlCats.push_back("altstep6"); mon.fillHisto("xsecalteventflowdilep", tags, 5, weight); mon.fillHisto("chhiggsalteventflowdilep", tags, 5, weight);}
        if(passMllVeto && passJetSelection && passMetSelection && passOS && passBtagsSelection_1 ) { ctrlCats.push_back("altstep6"); mon.fillHisto("xsecalteventflowdilep", tags, 5, weight); mon.fillHisto("chhiggsalteventflowdilep", tags, 5, weight);}
        if(passMllVeto && passJetSelection && passMetSelection && passOS && passBtagsSelection_2 ) { ctrlCats.push_back("altstep6"); mon.fillHisto("xsecalteventflowdilep", tags, 5, weight); mon.fillHisto("chhiggsalteventflowdilep", tags, 5, weight);}
        

        // Fill the control plots
        for(size_t k=0; k<ctrlCats.size(); ++k){
          TString icat(ctrlCats[k]);
          
          mon.fillHisto(icat+"nvtxraw", tags, nGoodPV, rawWeight);
          mon.fillHisto(icat+"nvtx",    tags, nGoodPV, weight   );
          mon.fillHisto(icat+"rho",     tags, rho,     weight   );

          // Lepton and dilepton control
          mon.fillHisto(icat+"leadpt",      tags, selLeptons[0].pt(),        weight);
          mon.fillHisto(icat+"trailerpt",   tags, selLeptons[1].pt(),        weight);
          mon.fillHisto(icat+"leadeta",     tags, fabs(selLeptons[0].eta()), weight);
          mon.fillHisto(icat+"trailereta",  tags, fabs(selLeptons[1].eta()), weight);

          double thetall(utils::cmssw::getArcCos<patUtils::GenericLepton>(selLeptons[0],selLeptons[1]));
          double sumpt(selLeptons[0].pt()+selLeptons[1].pt());
          // double mtsum(utils::cmssw::getMT<patUtils::GenericLepton,LorentzVector>(selLeptons[0],met)+utils::cmssw::getMT<patUtils::GenericLepton,LorentzVector>(selLeptons[1],met));
          double mtsum(higgs::utils::transverseMass(selLeptons[0].p4(),met.p4(),false)+higgs::utils::transverseMass(selLeptons[1].p4(),met.p4(),false));

          mon.fillHisto(icat+"yll",          tags, fabs(dileptonSystem.Rapidity()), weight);
          mon.fillHisto(icat+"mll",          tags, dileptonSystem.mass(),           weight);
          mon.fillHisto(icat+"ptll",         tags, dileptonSystem.pt(),             weight);
          mon.fillHisto(icat+"met",          tags, met.pt(),                        weight);
          mon.fillHisto(icat+"dilarccosine", tags, thetall,                         weight);
          mon.fillHisto(icat+"sumpt",        tags, sumpt,                           weight);
          mon.fillHisto(icat+"mtsum",        tags, mtsum,                           weight);
          mon.fillHisto(icat+"qt",           tags, dileptonSystem.pt(),             weight, true);
          // mon.fillHisto("qtraw",    tags, dileptonSystem.pt(),weight/triggerPrescale,true);                                                                                      
          if(selJets.size()>0){
          mon.fillHisto(icat+"leadjetpt",      tags, selJets[0].pt(),         weight);
          //mon.fillHisto(icat+"trailerpt",   tags, selLeptons[1].pt(),         weight);
          mon.fillHisto(icat+"leadjeteta",     tags, fabs (selJets[0].eta()), weight);
          //mon.fillHisto(icat+"trailereta",  tags, fabs (selLeptons[1].eta()), weight);
          }

          double
            ht    (sumpt),
            htb   (sumpt),
            htnol (0),
            htbnol(0);

          for(size_t ilep=0; ilep<2; ++ilep){
            double lepid(fabs(selLeptons[ilep].pdgId()));
            double leppt(selLeptons[ilep].pt());
            mon.fillHisto(icat+"inclusivept",       tags, leppt, weight);
            if(lepid==11) mon.fillHisto(icat+"pte", tags, leppt, weight);
            if(lepid==13) mon.fillHisto(icat+"ptmu",tags, leppt, weight);
          }


          // Tau control ??? 
          mon.fillHisto (icat+"ntaus",      tags, ntaus,                      weight);
          if(ntaus > 0){
            mon.fillHisto (icat+"tauleadpt", tags, selTaus[0].pt(), weight);
            mon.fillHisto (icat+"tauleadeta", tags, selTaus[0].eta(), weight);
          }

          for(size_t ijet=0; ijet<selJets.size(); ++ijet)
            {
              pat::Jet selJet= selJets[ijet];
              ht+=selJet.pt();
              htnol+=selJet.pt();
              if(ijet<selBJets.size()) // Hack in order to have only one loop. Exploits the fact that selBjets are a subset of selJets.
                {
                  htb+=selBJets[ijet].pt();
                  htbnol+=selBJets[ijet].pt();
                }
              
              double csv (selJet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
              mon.fillHisto ("csv", tags, csv, weight);
              if (!isMC) continue;
              int flavId = selJet.partonFlavour();
              TString jetFlav ("others");
              if (abs (flavId) == 5) jetFlav = "b";
              else if (abs (flavId) == 4) jetFlav = "c";
              mon.fillHisto ("csv" + jetFlav, tags, csv, weight);
            }
          
          ht+=met.pt();
          htb+=met.pt();
          htnol+=met.pt();
          htbnol+=met.pt();

          mon.fillHisto(icat+"ht",     tags, ht,              weight);
          mon.fillHisto(icat+"htb",    tags, htb,             weight);
          mon.fillHisto(icat+"htnol",  tags, htnol,           weight);
          mon.fillHisto(icat+"htbnol", tags, htbnol,          weight);
          
          mon.fillHisto(icat+"nbtags", tags, selBJets.size(), weight);
          mon.fillHisto(icat+"njets",  tags, selJets.size(), weight);
          
        }
        
        //
        // HISTOS FOR STATISTICAL ANALYSIS (include systematic variations)
        //
        //Fill histogram for posterior optimization, or for control regions
        
        if(passMllVeto && passOS)
          {

            for (size_t ivar = 0; ivar < nSystVars; ivar++)
              {
                TString var(systVars[ivar]);
                
                double iweight = weight;       //nominal
                
                //energy scale/resolution
                bool varyJesUp    (systVars[ivar] == "_jesup"   );
                bool varyJesDown  (systVars[ivar] == "_jesdown" );
                bool varyJerUp    (systVars[ivar] == "_jerup"   );
                bool varyJerDown  (systVars[ivar] == "_jerdown" );
                bool varyUmetUp   (systVars[ivar] == "_umetup"  );
                bool varyUmetDown (systVars[ivar] == "_umetdown");
                bool varyMesUp    (systVars[ivar] == "_mesup"   );
                bool varyMesDown  (systVars[ivar] == "_mesdown" );
                bool varyEesUp    (systVars[ivar] == "_eesup"   );
                bool varyEesDown  (systVars[ivar] == "_eesdown" );
                
                //pileup variations
                if (systVars[ivar] == "_puup")   iweight *= TotalWeight_plus;
                if (systVars[ivar] == "_pudown") iweight *= TotalWeight_minus;
                
                //btag
                bool varyBtagUp   (systVars[ivar]=="_btagup");
                bool varyBtagDown (systVars[ivar]=="_btagdown");
                
                //Here were the Q^2 variations on VV pT spectum
                
                
                //recompute MET/MT if JES/JER was varied
                LorentzVector newMET = mets[0].p4();
                
                if(varyJesUp)    newMET = mets[0].shiftedP4(pat::MET::METUncertainty::JetEnUp);
                if(varyJesDown)  newMET = mets[0].shiftedP4(pat::MET::METUncertainty::JetEnDown);
                if(varyJerUp)    newMET = mets[0].shiftedP4(pat::MET::METUncertainty::JetResUp);
                if(varyJerDown)  newMET = mets[0].shiftedP4(pat::MET::METUncertainty::JetResDown);
                if(varyUmetUp)   newMET = mets[0].shiftedP4(pat::MET::METUncertainty::UnclusteredEnUp);
                if(varyUmetDown) newMET = mets[0].shiftedP4(pat::MET::METUncertainty::UnclusteredEnDown);
                if(varyMesUp)    newMET = mets[0].shiftedP4(pat::MET::METUncertainty::MuonEnUp);
                if(varyMesDown)  newMET = mets[0].shiftedP4(pat::MET::METUncertainty::MuonEnDown);
                if(varyEesUp)    newMET = mets[0].shiftedP4(pat::MET::METUncertainty::ElectronEnUp);
                if(varyEesDown)  newMET = mets[0].shiftedP4(pat::MET::METUncertainty::ElectronEnDown);
                
                //if(varyLesUp)    newMET = met[utils::cmssw::LESUP]; //FIXME  must vary all leptons separately: MuonEnUp/MuonEnDown/ElectronEnUp/ElectronEnDown/TauEnUp/TauEnDown
                //if(varyLesDown)  newMET = met[utils::cmssw::LESDOWN];


                auto& selJets      = selJetsVar[""];        if(selJetsVar    .find(systVars[ivar].Data())!=selJetsVar    .end())selJets     = selJetsVar    [systVars[ivar].Data()];            
                auto& njets        = njetsVar [""];         if(njetsVar      .find(systVars[ivar].Data())!=njetsVar      .end())njets       = njetsVar      [systVars[ivar].Data()];
                auto& nbtags       = nbtagsVar[""];         if(nbtagsVar     .find(systVars[ivar].Data())!=nbtagsVar     .end())nbtags      = nbtagsVar     [systVars[ivar].Data()];
                auto& mindphijmet  = mindphijmetVar[""];    if(mindphijmetVar.find(systVars[ivar].Data())!=mindphijmetVar.end())mindphijmet = mindphijmetVar[systVars[ivar].Data()];

           
                
                pat::JetCollection finalSelJets;
                pat::JetCollection finalSelBJets;
                bool passLocalBveto (true);///passBtags);
                for (size_t ijet = 0; ijet < jets.size(); ijet++)
                  {
                    pat::Jet jet = jets[ijet];
                    double eta = jet.eta();
                    double pt = jet.pt();
                    if(isMC)
                      {
                        std::vector<float> varPt = utils::cmssw::smearJES(pt, eta, totalJESUnc);
                        if(varyJesUp)   pt = varPt[0] * jet.pt();
                        if(varyJesDown) pt = varPt[1] * jet.pt();
                        //  smearJER(float pt, float eta, float genPt)
                        //  float newJERSF(1.0);
                        //if(isMC)
                        //  {
                        //    const data::PhysicsObject_t &genJet=jets[ijet].getObject("genJet");
                        //    std::vector<float> smearJER=utils::cmssw::smearJER(jets[ijet].pt(),jets[ijet].eta(),genJet.pt());
                        //    newJERSF=smearJER[0];
                        //    rawJet *= newJERSF;
                        // if(varyJerUp)    pt=jets[ijet].getVal("jerup");
                        // if(varyJerDown)  pt=jets[ijet].getVal("jerdown");
                      }
                    
                    if (pt < 30 || fabs(eta) > 2.5) continue;
                    bool passPFloose = passPFJetID("Loose", jet); 
                    if (!passPFloose) continue;
                    
                    //cross-clean with selected leptons and photons
                    double minDRlj (9999.), minDRlg (9999.);
                    for (size_t ilep = 0; ilep < selLeptons.size(); ilep++)
                      minDRlj = TMath::Min (minDRlj, reco::deltaR (jet.p4(), selLeptons[ilep].p4()));
                    // don't want to mess with photon ID // for(size_t ipho=0; ipho<selPhotons.size(); ipho++)
                    // don't want to mess with photon ID //   minDRlg = TMath::Min( minDRlg, deltaR(jets[ijet].p4(),selPhotons[ipho].p4()) );
                    if (minDRlj < 0.4 /*|| minDRlg<0.4 */ ) continue;
                    
                    finalSelJets.push_back(jet);

                    int flavId = jet.partonFlavour();
                    bool hasCSVtag (jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > btagMedium);
                    if (varyBtagUp)
                      {
                        if (abs (flavId) == 5)      btsfutil.modifyBTagsWithSF(hasCSVtag, sfb + sfbunc,     beff);
                        else if (abs (flavId) == 4) btsfutil.modifyBTagsWithSF(hasCSVtag, sfb/5 + 2*sfbunc, beff);
                        else                        btsfutil.modifyBTagsWithSF(hasCSVtag, sfl + sflunc,     leff);
                      }
                    else if (varyBtagDown)
                      {
                        if (abs (flavId) == 5)      btsfutil.modifyBTagsWithSF(hasCSVtag, sfb - sfbunc,     beff);
                        else if (abs (flavId) == 4) btsfutil.modifyBTagsWithSF(hasCSVtag, sfb/5 - 2*sfbunc, beff);
                        else                        btsfutil.modifyBTagsWithSF(hasCSVtag, sfl - sflunc,     leff);
                      }
                    if(hasCSVtag)
                      finalSelBJets.push_back(jet);
                  }
                std::sort(finalSelJets.begin(),  finalSelJets.end(),  utils::sort_CandidatesByPt);
                std::sort(finalSelBJets.begin(), finalSelBJets.end(), utils::sort_CandidatesByPt);
                
                bool passFinalJetSelection(finalSelJets.size()>1);
                bool passFinalMetSelection(newMET.pt()>40.);
                bool passFinalBtagsSelection(finalSelBJets.size()>1);
                
                if(!passFinalJetSelection || !passFinalMetSelection || !passFinalBtagsSelection) continue;
                // Here fill stat plots
                mon.fillHisto("finalnbjets"+var, tags, finalSelBJets.size(), iweight);
                LorentzVector ttbarSystem(selLeptons[0].p4() + selLeptons[1].p4() + finalSelBJets[0].p4() + finalSelBJets[1].p4());
                double mt = higgs::utils::transverseMass(ttbarSystem,newMET,false);
                mon.fillHisto("finalmt"+var, tags, mt, iweight);
              }
          } // End stat analysis

      } // End dilepton full analysis
      
      
      // Single lepton full analysis
      //if(tags[1] == "singlemu" || tags[1] == "singlee"){
      if(isSingleMu || isSingleE){
        
        // Clean jet collection from selected taus
        pat::JetCollection
          selSingleLepJets, selSingleLepBJets;
        for (size_t ijet = 0; ijet < selJets.size(); ++ijet)
          {
            pat::Jet jet = selJets[ijet];
            
            double minDRtj(9999.);
            for(size_t itau=0; itau<selTaus.size(); ++itau)
              {
                minDRtj = TMath::Min(minDRtj, reco::deltaR(jet, selTaus[itau]));
              }
            if(minDRtj>0.4) selSingleLepJets.push_back(jet);
            
            bool hasCSVtag = (jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > btagMedium);
            if (isMC)
              {
                int flavId = jets[ijet].partonFlavour();
                if      (abs (flavId) == 5) btsfutil.modifyBTagsWithSF(hasCSVtag, sfb,   beff);
                else if (abs (flavId) == 4) btsfutil.modifyBTagsWithSF(hasCSVtag, sfb/5, beff);
                else                        btsfutil.modifyBTagsWithSF(hasCSVtag, sfl,   leff);
              }
            
            if(!hasCSVtag) continue;
            if(minDRtj>0.4) selSingleLepBJets.push_back(jets[ijet]);
          }
        
        std::sort(selSingleLepJets.begin(),  selSingleLepJets.end(),  utils::sort_CandidatesByPt);
        std::sort(selSingleLepBJets.begin(), selSingleLepBJets.end(), utils::sort_CandidatesByPt);
        
        
        mon.fillHisto("nvtx_pileup", tags, nGoodPV, weight);
        
        if(selLeptons.size()!=1 || nGoodPV==0) continue; // Veto requirement alredy applied during the event categoriziation
        //int id (abs (selLeptons[0].pdgId()));
        //weight *= isMC ? lepEff.getLeptonEfficiency(selLeptons[0].pt(), selLeptons[0].eta(), id, id == 11 ? "loose" : "tight").first : 1.0;        
        
        // Event selection booleans
        bool passJetSelection(selSingleLepJets.size()>1);
        bool passMetSelection(met.pt()>40.);
        bool passBtagsSelection(selSingleLepBJets.size()>0);
        bool passTauSelection(selTaus.size()==1);
        bool passOS(selTaus.size()>0 ? selLeptons[0].pdgId() * selTaus[0].pdgId() < 0 : 0);
        
        // Setting up control categories and fill up event flow histo
        std::vector < TString > ctrlCats;
        ctrlCats.clear ();
                                                                                                      { ctrlCats.push_back ("step1"); mon.fillHisto("xseceventflowslep", tags, 0, weight); mon.fillHisto("chhiggseventflowslep", tags, 0, weight); }
        if(passJetSelection   )                                                                       { ctrlCats.push_back ("step2"); mon.fillHisto("xseceventflowslep", tags, 1, weight); mon.fillHisto("chhiggseventflowslep", tags, 1, weight); }
        if(passJetSelection && passMetSelection )                                                     { ctrlCats.push_back ("step3"); mon.fillHisto("xseceventflowslep", tags, 2, weight); mon.fillHisto("chhiggseventflowslep", tags, 2, weight); }
        if(passJetSelection && passMetSelection && passBtagsSelection )                               { ctrlCats.push_back ("step4"); mon.fillHisto("xseceventflowslep", tags, 3, weight); mon.fillHisto("chhiggseventflowslep", tags, 3, weight); }
        if(passJetSelection && passMetSelection && passBtagsSelection && passTauSelection )           { ctrlCats.push_back ("step5"); mon.fillHisto("xseceventflowslep", tags, 4, weight); mon.fillHisto("chhiggseventflowslep", tags, 4, weight); }
        if(passJetSelection && passMetSelection && passBtagsSelection && passTauSelection && passOS ) { ctrlCats.push_back ("step6"); mon.fillHisto("xseceventflowslep", tags, 5, weight); mon.fillHisto("chhiggseventflowslep", tags, 5, weight); }
        

        bool passBtagsSelection_0(selSingleLepBJets.size()==0);
        bool passBtagsSelection_1(selSingleLepBJets.size()==1);
        bool passBtagsSelection_2(selSingleLepBJets.size()>1);

                                                                                    { ctrlCats.push_back("altstep1"); mon.fillHisto("xsecalteventflowslep", tags, 0, weight); mon.fillHisto("chhiggsalteventflowslep", tags, 0, weight); }
        if(passMetSelection)                                                        { ctrlCats.push_back("altstep2"); mon.fillHisto("xsecalteventflowslep", tags, 1, weight); mon.fillHisto("chhiggsalteventflowslep", tags, 1, weight); }
        if(passMetSelection && passTauSelection)                                    { ctrlCats.push_back("altstep3"); mon.fillHisto("xsecalteventflowslep", tags, 2, weight); mon.fillHisto("chhiggsalteventflowslep", tags, 2, weight); }
        if(passMetSelection && passTauSelection && passOS)                          { ctrlCats.push_back("altstep4"); mon.fillHisto("xsecalteventflowslep", tags, 3, weight); mon.fillHisto("chhiggsalteventflowslep", tags, 3, weight); }
        if(passMetSelection && passTauSelection && passOS && passBtagsSelection_0)  { ctrlCats.push_back("altstep5"); mon.fillHisto("xsecalteventflowslep", tags, 4, weight); mon.fillHisto("chhiggsalteventflowslep", tags, 4, weight); }
        if(passMetSelection && passTauSelection && passOS && passBtagsSelection_1)  { ctrlCats.push_back("altstep6"); mon.fillHisto("xsecalteventflowslep", tags, 5, weight); mon.fillHisto("chhiggsalteventflowslep", tags, 5, weight); }
        if(passMetSelection && passTauSelection && passOS && passBtagsSelection_2)  { ctrlCats.push_back("altstep7"); mon.fillHisto("xsecalteventflowslep", tags, 6, weight); mon.fillHisto("chhiggsalteventflowslep", tags, 6, weight); }

        // Fill the control plots
        for(size_t k=0; k<ctrlCats.size(); ++k){
          
          TString icat(ctrlCats[k]);
          mon.fillHisto(icat+"nvtxraw",    tags, nGoodPV,                            rawWeight);
          mon.fillHisto(icat+"nvtx",       tags, nGoodPV,                            weight   );
          mon.fillHisto(icat+"rho",        tags, rho,                                weight   );
          mon.fillHisto(icat+"leadpt",     tags, selLeptons[0].pt(),                 weight   );
          mon.fillHisto(icat+"trailerpt",  tags, selLeptons[1].pt(),                 weight   );
          mon.fillHisto(icat+"leadeta",    tags, fabs(selLeptons[0].eta()),          weight   );
          mon.fillHisto(icat+"trailereta", tags, fabs(selLeptons[1].eta()),          weight   );
          mon.fillHisto(icat+"ntaus",      tags, ntaus,                              weight   );
          mon.fillHisto(icat+"met",        tags, met.pt(),                           weight   );
          if(selSingleLepJets.size()>0){
            mon.fillHisto(icat+"leadjetpt",      tags, selSingleLepJets[0].pt(),         weight);
            //mon.fillHisto(icat+"trailerpt",   tags, selLeptons[1].pt(),         weight);
            mon.fillHisto(icat+"leadjeteta",     tags, fabs (selSingleLepJets[0].eta()), weight);
            //mon.fillHisto(icat+"trailereta",  tags, fabs (selLeptons[1].eta()), weight);
          }
          if(ntaus > 0){
            mon.fillHisto (icat+"tauleadpt", tags, selTaus[0].pt(),             weight);
            mon.fillHisto (icat+"tauleadeta", tags, selTaus[0].eta(),             weight);
          }
         
          
          mon.fillHisto(icat+"nbtags", tags, selSingleLepBJets.size(), weight);
          mon.fillHisto(icat+"njets",  tags, selSingleLepJets.size(), weight);
          // dilepton only           mon.fillHisto (icat+"zmass", tags, dileptonSystem.mass(),           weight);
          // dilepton only           mon.fillHisto (icat+"zy",    tags, fabs(dileptonSystem.Rapidity()), weight);
          // dilepton only           mon.fillHisto (icat+"zpt",   tags, dileptonSystem.pt(),             weight);
          // dilepton only           //these two are used to reweight photon -> Z, the 3rd is a control
          // dilepton only           mon.fillHisto (icat+"qt",    tags, dileptonSystem.pt(),             weight, true);
          // dilepton only           ///     mon.fillHisto("qtraw",    tags, dileptonSystem.pt(),weight/triggerPrescale,true); 
          
          for (size_t ijet = 0; ijet < selSingleLepJets.size(); ijet++)
            {
              if (selSingleLepJets[ijet].pt() < 30 || fabs (selSingleLepJets[ijet].eta()) > 2.5) continue;
              
              double csv (selSingleLepJets[ijet].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
              mon.fillHisto ("csv", tags, csv, weight);
              if (!isMC) continue;
              int flavId = selSingleLepJets[ijet].partonFlavour();
              TString jetFlav ("others");
              if (abs (flavId) == 5) jetFlav = "b";
              else if (abs (flavId) == 4) jetFlav = "c";
              mon.fillHisto ("csv" + jetFlav, tags, csv, weight);
            }
        }
        //
        // HISTOS FOR STATISTICAL ANALYSIS (include systematic variations)
        //
        //Fill histogram for posterior optimization, or for control regions
        
        if(passTauSelection && passOS)
          {
            for (size_t ivar = 0; ivar < nSystVars; ivar++)
              {
                TString var(systVars[ivar]);
                
                double iweight = weight;       //nominal
                
                //energy scale/resolution
                bool varyJesUp    (systVars[ivar] == "_jesup"   );
                bool varyJesDown  (systVars[ivar] == "_jesdown" );
                bool varyJerUp    (systVars[ivar] == "_jerup"   );
                bool varyJerDown  (systVars[ivar] == "_jerdown" );
                bool varyUmetUp   (systVars[ivar] == "_umetup"  );
                bool varyUmetDown (systVars[ivar] == "_umetdown");
                bool varyLesUp    (systVars[ivar] == "_lesup"   );
                bool varyLesDown  (systVars[ivar] == "_lesdown" );
                
                //pileup variations
                if (systVars[ivar] == "_puup")   iweight *= TotalWeight_plus;
                if (systVars[ivar] == "_pudown") iweight *= TotalWeight_minus;

                //btag
                bool varyBtagUp (systVars[ivar] == "_btagup");
                bool varyBtagDown (systVars[ivar] == "_btagdown");
                
                //Here were the Q^2 variations on VV pT spectum
                
                //recompute MET/MT if JES/JER was varied
                LorentzVector newMET = mets[0].p4();
                
                if(varyJesUp)    newMET = mets[0].shiftedP4(pat::MET::METUncertainty::JetEnUp);
                if(varyJesDown)  newMET = mets[0].shiftedP4(pat::MET::METUncertainty::JetEnDown);
                if(varyJerUp)    newMET = mets[0].shiftedP4(pat::MET::METUncertainty::JetResUp);
                if(varyJerDown)  newMET = mets[0].shiftedP4(pat::MET::METUncertainty::JetResDown);
                if(varyUmetUp)   newMET = mets[0].shiftedP4(pat::MET::METUncertainty::UnclusteredEnUp);
                if(varyUmetDown) newMET = mets[0].shiftedP4(pat::MET::METUncertainty::UnclusteredEnDown);
                //if(varyLesUp)    newMET = met[utils::cmssw::LESUP]; //FIXME  must vary all leptons separately: MuonEnUp/MuonEnDown/ElectronEnUp/ElectronEnDown/TauEnUp/TauEnDown
                //if(varyLesDown)  newMET = met[utils::cmssw::LESDOWN];
                
                pat::JetCollection finalSelSingleLepJets;
                pat::JetCollection finalSelSingleLepBJets;
                bool passLocalBveto (true);///passBtags);
                for (size_t ijet = 0; ijet < jets.size(); ijet++)
                  {
                    pat::Jet jet = jets[ijet];
                    double eta = jet.eta();
                    double pt = jet.pt();
                    if(isMC)
                      {
                        std::vector<float> varPt = utils::cmssw::smearJES(pt, eta, totalJESUnc);
                        if(varyJesUp)   pt = varPt[0] * jet.pt();
                        if(varyJesDown) pt = varPt[1] * jet.pt();
                        //  smearJER(float pt, float eta, float genPt)
                        //  float newJERSF(1.0);
                        //if(isMC)
                        //  {
                        //    const data::PhysicsObject_t &genJet=jets[ijet].getObject("genJet");
                        //    std::vector<float> smearJER=utils::cmssw::smearJER(jets[ijet].pt(),jets[ijet].eta(),genJet.pt());
                        //    newJERSF=smearJER[0];
                        //    rawJet *= newJERSF;
                        // if(varyJerUp)    pt=jets[ijet].getVal("jerup");
                        // if(varyJerDown)  pt=jets[ijet].getVal("jerdown");
                      }
                    
                    if (pt < 30 || fabs(eta) > 2.5) continue;
                    bool passPFloose = passPFJetID("Loose", jet); 
                    if (!passPFloose) continue;
                    
                    //cross-clean with selected leptons and photons
                    double minDRlj (9999.), minDRlg (9999.);
                    for (size_t ilep = 0; ilep < selLeptons.size(); ilep++)
                      minDRlj = TMath::Min (minDRlj, reco::deltaR (jet.p4(), selLeptons[ilep].p4()));
                    // don't want to mess with photon ID // for(size_t ipho=0; ipho<selPhotons.size(); ipho++)
                    // don't want to mess with photon ID //   minDRlg = TMath::Min( minDRlg, deltaR(jets[ijet].p4(),selPhotons[ipho].p4()) );
                    double minDRtj(9999.);
                    for(size_t itau=0; itau<selTaus.size(); ++itau)
                      {
                        minDRtj = TMath::Min(minDRtj, reco::deltaR(jet, selTaus[itau]));
                      }
                    if (minDRlj < 0.4 || minDRtj<0.4 ) continue;
                    
                    finalSelSingleLepJets.push_back(jet);
                    
                    int flavId(jet.partonFlavour());
                    bool hasCSVtag(jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > btagMedium);

                    if (varyBtagUp)
                      {
                        if (abs (flavId) == 5)      btsfutil.modifyBTagsWithSF(hasCSVtag, sfb + sfbunc,     beff);
                        else if (abs (flavId) == 4) btsfutil.modifyBTagsWithSF(hasCSVtag, sfb/5 + 2*sfbunc, beff);
                        else                        btsfutil.modifyBTagsWithSF(hasCSVtag, sfl + sflunc,     leff);
                      }
                    else if (varyBtagDown)
                      {
                        if (abs (flavId) == 5)      btsfutil.modifyBTagsWithSF(hasCSVtag, sfb - sfbunc,     beff);
                        else if (abs (flavId) == 4) btsfutil.modifyBTagsWithSF(hasCSVtag, sfb/5 - 2*sfbunc, beff);
                        else                        btsfutil.modifyBTagsWithSF(hasCSVtag, sfl - sflunc,     leff);
                      }
                    if(hasCSVtag)
                      finalSelSingleLepBJets.push_back(jet);
                  }
                std::sort(finalSelSingleLepJets.begin(),  finalSelSingleLepJets.end(),  utils::sort_CandidatesByPt);
                std::sort(finalSelSingleLepBJets.begin(), finalSelSingleLepBJets.end(), utils::sort_CandidatesByPt);

                bool passFinalJetSelection(finalSelSingleLepJets.size()>1);
                bool passFinalMetSelection(newMET.pt()>40.);
                bool passFinalBtagsSelection(finalSelSingleLepBJets.size()>0);
                
                if(!passFinalJetSelection || !passFinalMetSelection || !passFinalBtagsSelection) continue;
                // Here fill stat plots
                pat::Tau & tau = selTaus[0];
                reco::CandidatePtr leadChargedHadron = tau.leadChargedHadrCand();
                double tauR(leadChargedHadron->p() / tau.energy());  // Sic. It is momentum, not transverse momentum
                double tauY(2*leadChargedHadron->pt()/tau.et() - 1);
                //                Y= [ p_T^{trk} - (E_T - p_T^{trk} ]/E_T  = 2p_T^{trk}/E_T  - 1 . Which is practically Y = 2R' -1

                
                LorentzVector mutauSystem (0, 0, 0, 0);
                mutauSystem += selLeptons[0].p4();
                mutauSystem += tau.p4();
                
                mon.fillHisto("finalnbjets"         +var, tags, finalSelSingleLepBJets.size(), iweight);
                mon.fillHisto("finaltaur"           +var, tags, tauR, iweight);
                mon.fillHisto("finaltaupolarization"+var, tags, tauY, iweight);
                mon.fillHisto("finaldphilepmet"     +var, tags, fabs(deltaPhi(newMET.phi(), selLeptons[0].phi())), iweight);
                mon.fillHisto("finaldphitaumet"     +var, tags, fabs(deltaPhi(newMET.phi(), selTaus[0].phi())), iweight);
                mon.fillHisto("finaldphileptau"     +var, tags, fabs(deltaPhi(selLeptons[0].phi(), selTaus[0].phi())), iweight);
                mon.fillHisto("finaltaupt"          +var, tags, selTaus[0].pt(), iweight);
                mon.fillHisto("finalmutaumass"      +var, tags, mutauSystem.mass(), iweight);

                if(saveSummaryTree)
                  {
                    TDirectory* cwd = gDirectory;
                    summaryFile->cd();
                    summaryTree->Fill();
                    cwd->cd();
                  }
              }
          } // End stat analysis
        
      } // End single lepton full analysis

      } // End single file event loop 
    printf("\n");
    delete file;
  } // End loop on files
  
  if(saveSummaryTree)
    {
      TDirectory* cwd = gDirectory;
      summaryFile->cd();
      summaryTree->Write();
      summaryFile->Close();
      delete summaryFile;
      cwd->cd();
    }
  

  if(nMultiChannel>0) cout << "Warning! There were " << nMultiChannel << " multi-channel events out of " << totalEntries << " events!" << endl;
  printf ("\n");

  //##############################################
  //########     SAVING HISTO TO FILE     ########
  //##############################################
  //save control plots to file
  printf ("Results save in %s\n", outUrl.Data());

  //save all to the file
  TFile *ofile = TFile::Open (outUrl, "recreate");
  mon.Write();
  ofile->Close();

  if (outTxtFile)
    fclose (outTxtFile);

  // Now that everything is done, dump the list of lumiBlock that we processed in this job
  if(!isMC){
    goodLumiFilter.FindLumiInFiles(urls);
    goodLumiFilter.DumpToJson(((outUrl.ReplaceAll(".root",""))+".json").Data());
  }
  
}
