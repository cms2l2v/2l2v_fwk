#include <iostream>
#include <boost/shared_ptr.hpp>

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
//#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

//Load here all the dataformat that we will need
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
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
#include "UserCode/llvv_fwk/interface/TrigUtils.h"

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
#include "TLorentzVector.h"
#include <Math/VectorUtil.h>

using namespace std;

// Additional functions
LorentzVector getZCand(std::vector<patUtils::GenericLepton> selLeptons,LorentzVector& leadingLep,LorentzVector& trailerLep, bool& passBestZmass, double& BestMass,int& dilLep1,int& dilLep2, int& dilId,float rho){
   LorentzVector zlltmp;
   dilLep1=-1; dilLep2=-1; dilId=-1;
   BestMass=0;
   passBestZmass=false;
   LorentzVector zll(0.,0.,0.,0.);

   for(unsigned int l1=0   ;l1<selLeptons.size();l1++){
      int lid1=selLeptons[l1].pdgId();
      if(abs(lid1)==15)continue;

      bool passIso1 = abs(lid1)==11?patUtils::passIso(selLeptons[l1].el,  patUtils::llvvElecIso::Tight) : patUtils::passIso(selLeptons[l1].mu,  patUtils::llvvMuonIso::Tight);
      bool passLooseLepton1 = abs(lid1)==11?patUtils::passIso(selLeptons[l1].el,  patUtils::llvvElecIso::Loose) : patUtils::passIso(selLeptons[l1].mu,  patUtils::llvvMuonIso::Loose);

      // float relIso1 = utils::cmssw::relIso(selLeptons[l1].lep, rho);
      // if( relIso1>0.30 ) continue;

      if(!passIso1) continue;

      for(unsigned int l2=l1+1;l2<selLeptons.size();l2++){
         // float relIso2 = utils::cmssw::relIso(selLeptons[l2].lep, rho);
         // if( relIso2>0.30 ) continue;								 //ISO SUBLEADING

         int lid2 = selLeptons[l2].pdgId();

         if(abs(lid2)==15)continue;
         bool passIso2 = abs(lid2)==11?patUtils::passIso(selLeptons[l2].el,  patUtils::llvvElecIso::Tight) : patUtils::passIso(selLeptons[l2].mu,  patUtils::llvvMuonIso::Tight);
         bool passLooseLepton2 = abs(lid2)==11?patUtils::passIso(selLeptons[l2].el,  patUtils::llvvElecIso::Loose) : patUtils::passIso(selLeptons[l2].mu,  patUtils::llvvMuonIso::Loose);

         if(!passIso2) continue;

         if( selLeptons[l1].pt()<20 ) continue; 
         if( selLeptons[l2].pt()<20 ) continue; 
         //if( !((selLeptons[l1].pt()>=20 && selLeptons[l2].pt()>=10) || 
         //			(selLeptons[l1].pt()>=10 && selLeptons[l2].pt()>=20))) continue; //CUT ON PT
         if(abs(lid1)!=abs(lid2)) continue; 				 //SAME FLAVOUR PAIR
         if(lid1*lid2>=0) continue;					 //OPPOSITE SIGN
         if(deltaR(selLeptons[l1].p4(), selLeptons[l2].p4())<0.1) continue;
         zlltmp = (selLeptons[l1].p4()+selLeptons[l2].p4());
         if( fabs(zlltmp.mass() - 91.2) < fabs(BestMass-91.2) ){    //BEST MASS [76.2,106.2]
            dilLep1 = l1; 
            dilLep2 = l2;
            zll=selLeptons[l1].p4()+selLeptons[l2].p4();
            leadingLep=selLeptons[l1].p4();
            trailerLep=selLeptons[l2].p4();
            dilId = lid1 * lid2;
            BestMass=zll.mass();
            passBestZmass=true;
         }
      }
   }
   return zll;
}

LorentzVector getHiggsCand(std::vector<patUtils::GenericLepton> selLeptons, int dilLep1, int dilLep2, int& higgsCandL1, int& higgsCandL2, int& higgsCandId, int& HiggsShortId, vector<TString>& chTagsMain){
  higgsCandL1=-1;
  higgsCandL2=-1;
  LorentzVector higgsCand(0.,0.,0.,0.);
  HiggsShortId=-1;
  higgsCandId=0;
  
  for(int l=0   ;l<(int)selLeptons.size();l++){
    if(l==dilLep1 || l==dilLep2)continue;
    if(deltaR(selLeptons[l].p4(),  selLeptons[dilLep1].p4())<0.1)continue;
    if(deltaR(selLeptons[l].p4(),  selLeptons[dilLep2].p4())<0.1)continue;
    if(higgsCandL1<0){higgsCandL1=l;continue;}
    if(higgsCandL2<0 && deltaR(selLeptons[l].p4(),  selLeptons[higgsCandL1].p4())>=0.1){higgsCandL2=l;break;}//ordered in pT, so all done
  }
  string ChannelName = "none";   string signName = "";
  if(higgsCandL1>=0 && higgsCandL2>=0){
    higgsCandId=selLeptons[higgsCandL1].pdgId()*selLeptons[higgsCandL2].pdgId();
    higgsCand = LorentzVector(selLeptons[higgsCandL1].p4()+selLeptons[higgsCandL2].p4());
    if(higgsCandId<0){signName="_OS";}else{signName="_SS";}
    if(higgsCandId<0){HiggsShortId = 0;}else{HiggsShortId = 8;}
    if(abs(selLeptons[dilLep1].pdgId())==11){HiggsShortId += 0;}else{HiggsShortId += 4;}
    switch(abs(higgsCandId)){
    case 11*13:  ChannelName  = "elmu";  HiggsShortId+= 0; break;
    case 11*15:  ChannelName  = "elha";  HiggsShortId+= 1; break;
    case 13*15:  ChannelName  = "muha";  HiggsShortId+= 2; break;
    case 15*15:  ChannelName  = "haha";  HiggsShortId+= 3; break;
    default:     ChannelName  = "none";  HiggsShortId =-1; break;
    }
  }               
  chTagsMain.push_back(chTagsMain[chTagsMain.size()-1] + signName + ChannelName); 
  return higgsCand;
}


/*
LorentzVector getSVFit(pat::MET met, std::vector<patUtils::GenericLepton> selLeptons, int higgsCandL1, int higgsCandL2){       
	if(higgsCandL1<0 || higgsCandL2<0) return LorentzVector(0,0,0,0);

	TMatrixD covMET(2, 2); // PFMET significance matrix
	covMET[0][0] = met.sigx2;
	covMET[0][1] = met.sigxy;
	covMET[1][0] = met.sigxy;
	covMET[1][1] = met.sigy2;

	std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
	measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(abs(selLeptons[higgsCandL1].pdgId())==15?svFitStandalone::kTauToHadDecay:abs(selLeptons[higgsCandL1].pdgId())==11?svFitStandalone::kTauToElecDecay:svFitStandalone::kTauToMuDecay, svFitStandalone::LorentzVector(selLeptons[higgsCandL1].px(), selLeptons[higgsCandL1].py(), selLeptons[higgsCandL1].pz(), selLeptons[higgsCandL1].E()) ));
	measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(abs(selLeptons[higgsCandL2].pdgId())==15?svFitStandalone::kTauToHadDecay:abs(selLeptons[higgsCandL1].pdgId())==11?svFitStandalone::kTauToElecDecay:svFitStandalone::kTauToMuDecay, svFitStandalone::LorentzVector(selLeptons[higgsCandL2].px(), selLeptons[higgsCandL2].py(), selLeptons[higgsCandL2].pz(), selLeptons[higgsCandL2].E()) ));
	SVfitStandaloneAlgorithm algo(measuredTauLeptons, svFitStandalone::Vector(met.px(), met.py(), 0) , covMET, 0);
	algo.addLogM(false);
	algo.fit();
	if(algo.isValidSolution()){
		return algo.fittedDiTauSystem();
	}
	return LorentzVector(selLeptons[higgsCandL1]+selLeptons[higgsCandL2]);
};
*/


bool passHiggsCuts(std::vector<patUtils::GenericLepton> selLeptons, 
		   float rho, int higgsCandId, int higgsCandL1, int higgsCandL2, 
		   vector<float> isoElCut, vector<float> isoMuCut, const char* isoHaCut, vector<float> sumPtCut, 
		   bool requireId,
		   reco::VertexCollection vtx){
  
  if(higgsCandL1<0 || higgsCandL2<0)return false;
  
        // e - mu final state	
	if(abs(higgsCandId) == 11*13){
		patUtils::GenericLepton& lep1 = selLeptons[higgsCandL1];
		patUtils::GenericLepton& lep2 = selLeptons[higgsCandL2];

		LorentzVector lep1T, lep2T;
		lep1T = selLeptons[higgsCandL1].p4();
		lep2T = selLeptons[higgsCandL2].p4();

		cout << "ELE CUT ISO PASSED: " << isoElCut.at(0) << "MUO CUT ISO PASSED: " << isoMuCut.at(0) << endl;
		cout << "In the EMU final state: iso( lep1 )/iso( lep2 )--> "<< patUtils::relIso(lep1, rho) <<"/"<<patUtils::relIso(lep2, rho)<<endl;

		if(requireId){

		  int lid1=lep1.pdgId();
		  int lid2=lep2.pdgId();

		  bool passLooseLepton1 = lid1==11?patUtils::passId(lep1.el, vtx[0], patUtils::llvvElecId::Loose) : patUtils::passId(lep1.mu, vtx[0], patUtils::llvvMuonId::Loose);
		  bool passLooseLepton2 = lid1==11?patUtils::passId(lep2.el, vtx[0], patUtils::llvvElecId::Loose) : patUtils::passId(lep2.mu, vtx[0], patUtils::llvvMuonId::Loose);

		  if( ((lep1.pt()+lep2.pt()) >= sumPtCut.at(0))
		      &&  ( patUtils::relIso(lep1, rho) <= ( (abs(lep1.pdgId())==11)?isoElCut.at(0):isoMuCut.at(0) ) ) 
		      &&  ( patUtils::relIso(lep2, rho) <= ( (abs(lep2.pdgId())==11)?isoElCut.at(0):isoMuCut.at(0) ) ) 
		      && passLooseLepton1  && passLooseLepton2 ){
		    cout << "IMMMMMMMM IN iso( lep1 )/iso( lep2 )--> "<< patUtils::relIso(lep1, rho) <<"/"<<patUtils::relIso(lep2, rho)<<endl;
		    return true; 
		  } 
		}else{
		  if(((lep1.pt()+lep2.pt()) >= sumPtCut.at(0))
		     &&  (patUtils::relIso(lep1, rho) <= ( (abs(lep1.pdgId())==11)?isoElCut.at(0):isoMuCut.at(0) ) )
		     &&  (patUtils::relIso(lep2, rho) <= ( (abs(lep2.pdgId())==11)?isoElCut.at(0):isoMuCut.at(0) ) ) ){
		    return true;
		  }
		} 
	}

	// e - tau final state
	else if(abs(higgsCandId) == 11*15){
	  pat::Tau& tau = (abs(selLeptons[higgsCandL1].pdgId())==15)?selLeptons[higgsCandL1].tau:selLeptons[higgsCandL2].tau;
	  patUtils::GenericLepton& lep = (abs(selLeptons[higgsCandL1].pdgId())==15)?selLeptons[higgsCandL2]:selLeptons[higgsCandL1];
	  int tauIdx = (abs(selLeptons[higgsCandL1].pdgId())==15)?higgsCandL1:higgsCandL2;
	  LorentzVector tauT;
	  tauT = selLeptons[tauIdx].p4();
	  
	  float relIso1 = patUtils::relIso(lep, rho);
	  if(requireId){
	    
	    bool passId = patUtils::passId(lep.el, vtx[0], patUtils::llvvElecId::Tight);
	    
	    if((relIso1<=isoElCut.at(1)) 
	       && passId
	       && tau.tauID("againstElectronTightMVA5") 
	       && tau.tauID("againstMuonLoose3")
	       && tau.tauID(isoHaCut) 
	       && ((tau.pt()+lep.pt()) >= sumPtCut.at(1)) ){
	      return true;
	    } 
	  }else{
	    if((relIso1<=isoElCut.at(1)) 
	       && tau.tauID("againstElectronTightMVA5") 
	       && tau.tauID("againstMuonLoose3")
	       && tau.tauID(isoHaCut) 
	       && ((tau.pt()+lep.pt()) >= sumPtCut.at(1)) ){
	      return true;
	    }
	  } 
	}
	// mu - tau final state
	else if(abs(higgsCandId) == 13*15){
	  pat::Tau& tau = abs(selLeptons[higgsCandL1].pdgId())==15?selLeptons[higgsCandL1].tau:selLeptons[higgsCandL2].tau;
	  patUtils::GenericLepton& lep = abs(selLeptons[higgsCandL1].pdgId())==15?selLeptons[higgsCandL2]:selLeptons[higgsCandL1];
	  int tauIdx = abs(selLeptons[higgsCandL1].pdgId())==15?higgsCandL1:higgsCandL2;
	  LorentzVector tauT;
	  tauT = selLeptons[tauIdx].p4();
	  //cout << "In the MUTAU final state: passLooseId(lep)/tau.tauID("againstElectronLooseMVA5")/tau.tauID("againstMuonTight3") "<<
	  //passLooseId( lep )<<"/"<<tau.tauID("againstElectronLooseMVA5")<<"/"<<tau.tauID("againstMuonTight3")<<endl;
	  
	  float relIso1 = patUtils::relIso(lep, rho);
	  if(requireId){

	    bool passId = patUtils::passId(lep.mu, vtx[0], patUtils::llvvMuonId::Loose);

	    if((relIso1<=isoMuCut.at(1))
	       && passId
	       && tau.tauID("againstElectronLooseMVA5")
	       && tau.tauID("againstMuonTight3")
	       && tau.tauID(isoHaCut)
	       && ((tau.pt()+lep.pt()) >= sumPtCut.at(2)) ){
	      return true;
	    }    
		}else{
	    if((relIso1<=isoMuCut.at(1))
	       && tau.tauID("againstElectronLooseMVA5")
	       && tau.tauID("againstMuonTight3")
	       && tau.tauID(isoHaCut)
	       && ((tau.pt()+lep.pt()) >= sumPtCut.at(2)) ){
	      return true;
	    }    
	  } 
	}
	// tau - tau final state
	else if(abs(higgsCandId) == 15*15){
		pat::Tau& tau1 = selLeptons[higgsCandL1].tau;
		pat::Tau& tau2 = selLeptons[higgsCandL2].tau;
		if(((tau1.pt()+tau2.pt()) >= sumPtCut.at(3))
		   && tau1.tauID("againstElectronLooseMVA5") 
		   && tau1.tauID("againstMuonLoose3") 
		   && tau2.tauID("againstElectronLooseMVA5") 
		   && tau2.tauID("againstMuonLoose3") 
		   && tau1.tauID(isoHaCut) 
		   && tau2.tauID(isoHaCut)  
		   ){
		  return true;
		}
	}               
	return false;
} 

bool passHiggsCuts(std::vector<patUtils::GenericLepton> selLeptons, float rho, int higgsCandId, int higgsCandL1, int higgsCandL2, float isoElCut, float isoMuCut, const char* isoHaCut, float sumPtCut, bool requireId, reco::VertexCollection vtx){
	vector<float> isoElCutV; 
	vector<float> isoMuCutV; 
	vector<float> sumPtCutV; 
	isoElCutV.push_back(isoElCut);
	isoElCutV.push_back(isoElCut);
	isoMuCutV.push_back(isoMuCut);
	isoMuCutV.push_back(isoMuCut);
	sumPtCutV.push_back(sumPtCut);
	sumPtCutV.push_back(sumPtCut);
	sumPtCutV.push_back(sumPtCut);
	sumPtCutV.push_back(sumPtCut);
	bool passHiggs = passHiggsCuts(selLeptons, rho, higgsCandId, higgsCandL1, higgsCandL2, isoElCutV, isoMuCutV, isoHaCut, sumPtCutV, requireId, vtx);
	return passHiggs;
}


int main(int argc, char* argv[])
{
  //##############################################
  //########    GLOBAL INITIALIZATION     ########
  //##############################################

  // check arguments
  if(argc<2){ std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl; exit(0); }
  
  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();
  
  // configure the process
  const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");

  bool isMC = runProcess.getParameter<bool>("isMC");
  double xsec = runProcess.getParameter<double>("xsec");
  int mctruthmode=runProcess.getParameter<int>("mctruthmode");
  bool photonTriggerStudy = runProcess.getParameter<bool>("triggerstudy");
  TString dtag=runProcess.getParameter<std::string>("dtag");

  TString suffix=runProcess.getParameter<std::string>("suffix");
  std::vector<std::string> urls=runProcess.getUntrackedParameter<std::vector<std::string> >("input");

  TString outUrl = runProcess.getParameter<std::string>("outfile");

  //good lumi MASK
  lumiUtils::GoodLumiFilter goodLumiFilter(runProcess.getUntrackedParameter<std::vector<edm::LuminosityBlockRange> >("lumisToProcess", std::vector<edm::LuminosityBlockRange>()));

  bool filterOnlyEE(false), filterOnlyMUMU(false), filterOnlyEMU(false);
  if( isMC ) std::cout << "Is MC" << std::endl;
  if(!isMC)
    {
      if(dtag.Contains("DoubleEle")) filterOnlyEE=true;
      if(dtag.Contains("DoubleMu"))  filterOnlyMUMU=true;
      if(dtag.Contains("MuEG"))      filterOnlyEMU=true;
    }
  bool isSingleMuPD(!isMC && dtag.Contains("SingleMu"));  
  bool isV0JetsMC(false);//isMC && (dtag.Contains("DYJetsToLL_50toInf") || dtag.Contains("_WJets")));  #FIXME should be reactivated as soon as we have exclusive jet samples
  bool isWGmc(isMC && dtag.Contains("WG"));
  bool isZGmc(isMC && dtag.Contains("ZG"));
  bool isMC_GG  = isMC && ( string(dtag.Data()).find("GG" )  != string::npos);
  bool isMC_VBF = isMC && ( string(dtag.Data()).find("VBF")  != string::npos);
  bool isMC_VBF1000(isMC && dtag.Contains("VBFtoH1000toZZto2L2Nu"));
  bool isMC_125OnShell = isMC && (mctruthmode==521);
  if(isMC_125OnShell) mctruthmode=125;
  bool isMC_ZZ  = isMC && ( string(dtag.Data()).find("TeV_ZZ")  != string::npos);
  bool isMC_WZ  = isMC && ( string(dtag.Data()).find("TeV_WZ")  != string::npos);
  bool isMC_QCD = (isMC && dtag.Contains("QCD"));
  bool isMC_GJet = (isMC && dtag.Contains("GJet"));
 
  //Tag for Met Filter
  bool isPromptReco (!isMC && dtag.Contains("PromptReco")); //"False" picks up correctly the new prompt reco (2015C) and MC

  TString outTxtUrl= outUrl + ".txt";    
  FILE* outTxtFile = NULL;
  if(!isMC)outTxtFile = fopen(outTxtUrl.Data(), "w");
  printf("TextFile URL = %s\n",outTxtUrl.Data());

  //Blind the sensible region
  bool blindData=true;

  //tree info
  TString dirname = runProcess.getParameter<std::string>("dirName");

  //systematics
  bool runSystematics                        = runProcess.getParameter<bool>("runSystematics");
  std::vector<TString> varNames(1,"");
  if(runSystematics){
    varNames.push_back("_jerup");    varNames.push_back("_jerdown");
    varNames.push_back("_jesup");    varNames.push_back("_jesdown");  
    varNames.push_back("_umetup");   varNames.push_back("_umetdown");  
    varNames.push_back("_lesup");    varNames.push_back("_lesdown");  
    varNames.push_back("_puup");     varNames.push_back("_pudown");  
    varNames.push_back("_btagup");   varNames.push_back("_btagdown");
    if(isMC_ZZ)             { varNames.push_back("_zzptup");   varNames.push_back("_zzptdown");     }
    if(isMC_WZ)             { varNames.push_back("_wzptup");   varNames.push_back("_wzptdown");     }
    if(isMC_GG || isMC_VBF) { varNames.push_back("_lshapeup"); varNames.push_back("_lshapedown"); }
  }
  size_t nvarsToInclude=varNames.size();
  
  std::vector<std::string> allWeightsURL=runProcess.getParameter<std::vector<std::string> >("weightsFile");
  std::string weightsDir( allWeightsURL.size() ? allWeightsURL[0] : "");

  //shape uncertainties for dibosons
  std::vector<TGraph *> vvShapeUnc;
  if(isMC_ZZ || isMC_WZ)
    {
      TString weightsFile=weightsDir+"/zzQ2unc.root";
      TString dist("zzpt");
      if(isMC_WZ) { weightsFile.ReplaceAll("zzQ2","wzQ2"); dist.ReplaceAll("zzpt","wzpt"); }
      gSystem->ExpandPathName(weightsFile);
      TFile *q2UncF=TFile::Open(weightsFile);
      if(q2UncF){
	vvShapeUnc.push_back( new TGraph( (TH1 *)q2UncF->Get(dist+"_up") ) );
	vvShapeUnc.push_back( new TGraph( (TH1 *)q2UncF->Get(dist+"_down") ) );
	q2UncF->Close();
      }
    }

  
  //STANDARD MODEL
  double HiggsMass=0; string VBFString = ""; string GGString("");
  TF1 *decayProbPdf=new TF1("relbw","(2*sqrt(2)*[0]*[1]*sqrt(pow([0],2)*(pow([0],2)+pow([1],2)))/(TMath::Pi()*sqrt(pow([0],2)+sqrt(pow([0],2)*(pow([0],2)+pow([1],2))))))/(pow(pow(x,2)-pow([0],2),2)+pow([0]*[1],2))",0,2000);
  if(isMC_GG){  
    size_t GGStringpos =  string(dtag.Data()).find("GG");
    string StringMass = string(dtag.Data()).substr(GGStringpos+5,4);  sscanf(StringMass.c_str(),"%lf",&HiggsMass);
    GGString = string(dtag.Data()).substr(GGStringpos);  
  }else if(isMC_VBF){
    size_t VBFStringpos =  string(dtag.Data()).find("VBF");
    string StringMass = string(dtag.Data()).substr(VBFStringpos+6,4);  sscanf(StringMass.c_str(),"%lf",&HiggsMass);
    VBFString = string(dtag.Data()).substr(VBFStringpos);
  }
  if(mctruthmode==125) HiggsMass=124;

  //##############################################
  //########    INITIATING HISTOGRAMS     ########
  //##############################################
  SmartSelectionMonitor mon;
  printf("Definition of plots");

  //generator level control : add an underflow entry to make sure the histo is kept
  ((TH1F*)mon.addHistogram( new TH1F( "higgsMass_raw",     ";Higgs Mass [GeV];Events", 500,0,1500) ))->Fill(-1.0,0.0001);
  ((TH1F*)mon.addHistogram( new TH1F( "higgsMass_cpspint", ";Higgs Mass [GeV];Events", 500,0,1500) ))->Fill(-1.0,0.0001);

  mon.addHistogram( new TH1F( "wdecays",     ";W decay channel",5,0,5) );
  mon.addHistogram( new TH1F( "zdecays",     ";Z decay channel",6,0,6) );

  //event selection
  TH1F *h=(TH1F*) mon.addHistogram( new TH1F ("eventflow", ";;Events", 9,0,9) );
  h->GetXaxis()->SetBinLabel(1,"raw");
  h->GetXaxis()->SetBinLabel(2,"#geq 2 iso leptons");
  h->GetXaxis()->SetBinLabel(3,"|M-91|<15");
  h->GetXaxis()->SetBinLabel(4,"p_{T}>55");
  h->GetXaxis()->SetBinLabel(5,"3^{rd}-lepton veto");
  h->GetXaxis()->SetBinLabel(6,"b-veto"); 
  h->GetXaxis()->SetBinLabel(7,"#Delta #phi(jet,E_{T}^{miss})>0.5");
  h->GetXaxis()->SetBinLabel(8,"E_{T}^{miss}>80");

  // new event selection
  TH1 *h1=mon.addHistogram( new TH1F ("eventflow2", ";;Events", 10,0,10) );
  h1->GetXaxis()->SetBinLabel(1,"InitialEv");
  h1->GetXaxis()->SetBinLabel(2,"Nlep#geq2");
  h1->GetXaxis()->SetBinLabel(3,"Zmass");
  h1->GetXaxis()->SetBinLabel(4,"Zkin");
  h1->GetXaxis()->SetBinLabel(5,"Nlep+Ntau#geq4"); 
  h1->GetXaxis()->SetBinLabel(6,"Lep Veto");
  h1->GetXaxis()->SetBinLabel(7,"Btag Veto");
  h1->GetXaxis()->SetBinLabel(8,"#Delta #phi Z-MET");
  h1->GetXaxis()->SetBinLabel(9,"di-#tau Cand");

  TH1 *h2=mon.addHistogram( new TH1F ("yields", ";;Events", 17,0,17) );
  h2->GetXaxis()->SetBinLabel(1,"OS eee#mu");
  h2->GetXaxis()->SetBinLabel(2,"OS eee#tau");
  h2->GetXaxis()->SetBinLabel(3,"OS ee#mu#tau");
  h2->GetXaxis()->SetBinLabel(4,"OS ee#tau#tau");
  h2->GetXaxis()->SetBinLabel(5,"OS #mu#mue#mu");
  h2->GetXaxis()->SetBinLabel(6,"OS #mu#mue#tau");
  h2->GetXaxis()->SetBinLabel(7,"OS #mu#mu#mu#tau");
  h2->GetXaxis()->SetBinLabel(8,"OS #mu#mu#tau#tau");
  h2->GetXaxis()->SetBinLabel(9,"SS eee#mu");
  h2->GetXaxis()->SetBinLabel(10,"SS eee#tau");
  h2->GetXaxis()->SetBinLabel(11,"SS ee#mu#tau");
  h2->GetXaxis()->SetBinLabel(12,"SS ee#tau#tau");
  h2->GetXaxis()->SetBinLabel(13,"SS #mu#mue#mu");
  h2->GetXaxis()->SetBinLabel(14,"SS #mu#mue#tau");
  h2->GetXaxis()->SetBinLabel(15,"SS #mu#mu#mu#tau");
  h2->GetXaxis()->SetBinLabel(16,"SS #mu#mu#tau#tau");
  
  TH1 *h3=mon.addHistogram( new TH1F ("yieldsOS", ";;Events", 8,0,8) );
  h3->GetXaxis()->SetBinLabel(1,"OS eee#mu");
  h3->GetXaxis()->SetBinLabel(2,"OS eee#tau");
  h3->GetXaxis()->SetBinLabel(3,"OS ee#mu#tau");
  h3->GetXaxis()->SetBinLabel(4,"OS ee#tau#tau");
  h3->GetXaxis()->SetBinLabel(5,"OS #mu#mue#mu");
  h3->GetXaxis()->SetBinLabel(6,"OS #mu#mue#tau");
  h3->GetXaxis()->SetBinLabel(7,"OS #mu#mu#mu#tau");
  h3->GetXaxis()->SetBinLabel(8,"OS #mu#mu#tau#tau");
  
  //bjets control
  mon.addHistogram( new TH1F( "nbjets",  ";Number of #b-jets;Events", 6,0,6) );
  mon.addHistogram( new TH1F( "njets2",   ";Number of #jets;Events", 6,0,6) );
  mon.addHistogram( new TH1F( "bjetpt",  ";p_{T}^{#b-jet} (GeV);Events/10 GeV", 50,0,500) );
  mon.addHistogram( new TH1F( "bjetcsv", ";CSV b-tagged jet;Events", 50,0, 1) );
  
  //boson control
  mon.addHistogram( new TH1F( "zy",      		";y_{ll};Events", 50,-6,6) );
  mon.addHistogram( new TH1F( "zeta",    		";#eta_{ll};Events", 50,-10,10) );
  mon.addHistogram( new TH1F( "zpt",     		";p_{T}^{ll} (GeV) ;Events/10 GeV", 50,0,500) );
  mon.addHistogram( new TH1F( "zmass",   		";M_{ll} (GeV);Events/2 GeV", 80,20,180) );
  mon.addHistogram( new TH1F( "dPhi_ZMet",              ";#Delta#phi(ll,#slash{E}_{T});Events",50,-3,3));
  
  // zll boson control
  mon.addHistogram( new TH1F( "zlly",      		";y_{ll};Events", 50,-6,6) );
  mon.addHistogram( new TH1F( "zlleta",    		";#eta_{ll};Events", 50,-10,10) );
  mon.addHistogram( new TH1F( "zllpt",     		";p_{T}^{ll} (GeV) ;Events/10 GeV", 50,0,500) );
  mon.addHistogram( new TH1F( "zllmass",   		";M_{ll} (GeV);Events/2 GeV", 80,20,180) );
  mon.addHistogram( new TH1F( "dPhi_ZllMet",        ";#Delta#phi(ll,#slash{E}_{T});Events",50,-3,3));

  mon.addHistogram( new TH1F( "sumpt",            ";L_{T} (GeV);Events/10 GeV", 50,0,500) );
  mon.addHistogram( new TH1F( "dPhi_AZ",          ";#DeltaPhi(#tau#tau,ll);Events",50,-3,3));
  mon.addHistogram( new TH1F( "dPhi_AMet",        ";#Delta#phi(#tau#tau,#slash{E}_{T});Events",50,-3,3));
  mon.addHistogram( new TH1F( "met2",             ";#slash{E}_{T} (GeV);Events/10 GeV",50,0,500));
  
  mon.addHistogram( new TH1F( "Amet",             ";#slash{E}_{T} (GeV);Events/10 GeV",50,0,500));
  mon.addHistogram( new TH1F( "Anjets",           ";Number of Jets;Events",10,-0.5,9.5));
  mon.addHistogram( new TH1F( "Apt",              ";p_{T}^{#tau#tau} (GeV);Events/10 GeV",50,0,500));
  mon.addHistogram( new TH1F( "Hpt",              ";p_{T}^{ll#tau#tau} (GeV);Events/10 GeV",50,0,500));
  
  //double bins[]={0.,32.4,70.2,110.2,170.,306.,550.8,1785.};      //NEW BINNING 1
  //double bins[]={0.,32.4,70.2,110.2,189.,306.,550.8,1785.};        //NEW BINNING 2
  //double bins[]={0.1,32.4,70.2,110.2,170.,189.,306.,550.8,1785.}; //OLD BINNING
  double bins[]={5.0,32.4,70.2,110.2,170.,189.,306.,550.8,1785.}; //OLD BINNING
  
  mon.addHistogram( new TH1F( "Amass",            ";M_{#tau#tau} (GeV);Events",8,bins));
  mon.addHistogram( new TH1F( "Hmass",            ";M_{ll#tau#tau} (GeV);Events",8,bins));
  mon.addHistogram( new TH1F( "Amasssvfit",       ";SVFit M_{#tau#tau} (GeV);Events",8,bins));
  mon.addHistogram( new TH1F( "Hmasssvfit",       ";SVFit M_{ll#tau#tau} (GeV);Events",8,bins));
  
  //pu control
  mon.addHistogram( new TH1F( "nvtxA",";Vertices;Events",50,0,50) ); 
  mon.addHistogram( new TH1F( "nvtxB",";Vertices;Events",50,0,50) ); 
  mon.addHistogram( new TH1F( "nvtxC",";Vertices;Events",50,0,50) ); 
  mon.addHistogram( new TH1F( "nvtxD",";Vertices;Events",50,0,50) ); 

  mon.addHistogram( new TH1F( "nvtx",";Vertices;Events",50,0,50) ); 
  mon.addHistogram( new TH1F( "nvtxraw",";Vertices;Events",50,0,50) ); 
  mon.addHistogram( new TH1F( "rho",";#rho;Events",50,0,25) ); 

  mon.addHistogram( new TH1F( "nvtx2",";Vertices;Events",50,0,50) ); 
  mon.addHistogram( new TH1F( "nvtxraw2",";Vertices;Events",50,0,50) ); 
  mon.addHistogram( new TH1F( "rho2",";#rho;Events",50,0,25) );

  //tau control
  mon.addHistogram( new TH1F( "leadtaupt",     ";Transverse momentum [GeV];Events", 50,0,500) );
  TH1 *htaus=mon.addHistogram( new TH1F("ntaus",  ";Tau multiplicity;Events",5,0,5) );
  for(int ibin=1; ibin<=htaus->GetXaxis()->GetNbins(); ibin++)
    {
      TString label("");
      if(ibin==h->GetXaxis()->GetNbins()) label +="#geq";
      else                                label +="=";
      label += (ibin-1);
      htaus->GetXaxis()->SetBinLabel(ibin,label);
    } 

  // photon control
  mon.addHistogram(new TH1F("npho", ";Number of Photons;Events", 20, 0, 20) ); 
  mon.addHistogram(new TH1F("phopt", ";Photon pT [GeV];Events", 500, 0, 1000) ); 
  mon.addHistogram(new TH1F("phoeta", ";Photon pseudo-rapidity;Events", 50, 0, 5) );

  mon.addHistogram(new TH1F("bosonphi", ";Photon #phi;Events", 40, 0, 4) );
  mon.addHistogram(new TH1F("metphi", ";MET #phi;Events", 40, 0, 4) );
  mon.addHistogram(new TH1F("dphi_boson_met", ";#Delta #phi(#gamma,MET);Events", 40, 0, 4) );
  
  //lepton control hztautau
  mon.addHistogram( new TH1F( "nlep2", 	    ";Number of Leptons (e/#mu);Events", 10,0,10) );
  mon.addHistogram( new TH1F( "leadpt2",    ";p_{T}^{e/#mu} (GeV);Events/10 GeV", 50,0,500) );
  mon.addHistogram( new TH1F( "leadeta2",   ";#eta_{e/mu};Events", 50,-2.6,2.6) );
  mon.addHistogram( new TH1F( "trailerpt2", ";p_{T}^{e/#mu} (GeV) ;Events/10 GeV", 50,0,500) );
  mon.addHistogram( new TH1F( "trailereta2",";#eta_{e/#mu};Events", 50,-2.6,2.6) );
  mon.addHistogram( new TH1F( "leppt2",     ";p_{T}^{e/#mu} (GeV);Events/10 GeV", 50,0,500) );
  mon.addHistogram( new TH1F( "lepeta2",    ";#eta_{e/mu};Events", 50,-2.6,2.6) );




  // lepton control
  mon.addHistogram( new TH1F( "leadpt",     ";Transverse momentum [GeV];Events", 50,0,500) );
  mon.addHistogram( new TH1F( "leadeta",    ";Pseudo-rapidity;Events", 50,0,2.6) );
  mon.addHistogram( new TH1F( "trailerpt",  ";Transverse momentum [GeV];Events", 50,0,500) );
  mon.addHistogram( new TH1F( "trailereta", ";Pseudo-rapidity;Events", 50,0,2.6) );
  mon.addHistogram( new TH1F( "zy",         ";Rapidity;Events", 50,0,3) );
  mon.addHistogram( new TH1F( "zmass",      ";Mass [GeV];Events / 2 GeV", 100,40,240) );
  mon.addHistogram( new TH1F( "zpt",        ";Transverse momentum [GeV];Events",100,0,1500));
  Double_t zptaxis[]= {0,15,30,45,60,75,90,105,120,135,150,165,180,195,210,225,240,255,270,285,300,315,330,345,360,375,390,405,435,465,495,525,555,585,615,675,735,795,855,975};
  Int_t nzptAxis=sizeof(zptaxis)/sizeof(Double_t);
  mon.addHistogram( new TH1F( "zpt_rebin",  ";Transverse momentum [GeV];Events / GeV", nzptAxis-1,zptaxis));
  mon.addHistogram( new TH1F( "qmass",      ";Mass [GeV];Events / (1 GeV)",100,76,106));
  mon.addHistogram( new TH1F( "qt",         ";Transverse momentum [GeV];Events / (1 GeV)",1500,0,1500));
  mon.addHistogram( new TH1F( "qtraw",      ";Transverse momentum [GeV];Events / (1 GeV)",1500,0,1500));

  //extra leptons in the event
  mon.addHistogram( new TH1F( "nextraleptons", ";Extra leptons;Events",4,0,4) );
  mon.addHistogram( new TH1F( "thirdleptonpt", ";Transverse momentum;Events", 50,0,500) );
  mon.addHistogram( new TH1F( "thirdleptoneta", ";Pseudo-rapidity;Events", 50,0,2.6) );
  mon.addHistogram( new TH1F( "thirdleptonmt", ";Transverse mass(3^{rd} lepton,E_{T}^{miss}) [GeV];Events", 50,0,500) );


  mon.addHistogram( new TH1F("csv",      ";Combined Secondary Vertex;Jets",50,0.,1.) );
  mon.addHistogram( new TH1F("csvb",     ";Combined Secondary Vertex;Jets",50,0.,1.) );
  mon.addHistogram( new TH1F("csvc",     ";Combined Secondary Vertex;Jets",50,0.,1.) );
  mon.addHistogram( new TH1F("csvothers",";Combined Secondary Vertex;Jets",50,0.,1.) );
  TH1 *hbtags=mon.addHistogram( new TH1F("nbtags",   ";b-tag multiplicity;Events",5,0,5) );
  TH1 *hbtagsJP=mon.addHistogram( new TH1F("nbtagsJP",   ";b-tag multiplicity;Events",5,0,5) );
  mon.addHistogram( new TH1F("leadjetpt",    ";Transverse momentum [GeV];Events",50,0,1000) );
  mon.addHistogram( new TH1F("trailerjetpt", ";Transverse momentum [GeV];Events",50,0,1000) );
  mon.addHistogram( new TH1F("fwdjeteta",    ";Pseudo-rapidity;Events",25,0,5) );
  mon.addHistogram( new TH1F("cenjeteta",       ";Pseudo-rapidity;Events",25,0,5) );
  Double_t mjjaxis[32];
  mjjaxis[0]=0.01;
  for(size_t i=1; i<20; i++)  mjjaxis[i]   =50*i;        //0-1000
  for(size_t i=0; i<5; i++)   mjjaxis[20+i]=1000+100*i; //1000-1500
  for(size_t i=0; i<=5; i++)   mjjaxis[25+i]=1500+300*i; //1500-5000  
  mjjaxis[31]=5000;
  mon.addHistogram( new TH1F("vbfmjj"       , ";Dijet invariant mass [GeV];Events",31,mjjaxis) );
  mon.addHistogram( new TH1F("vbfdphijj"    , ";Azimuthal angle difference;Events",20,0,3.5) );
  mon.addHistogram( new TH1F("vbfdetajj"    , ";Pseudo-rapidity span;Events",20,0,10) );
  TH1 *hjets=mon.addHistogram( new TH1F("njets",  ";Jet multiplicity;Events",5,0,5) );
  for(int ibin=1; ibin<=hjets->GetXaxis()->GetNbins(); ibin++)
    {
      TString label("");
      if(ibin==h->GetXaxis()->GetNbins()) label +="#geq";
      else                                label +="=";
      label += (ibin-1);
      hjets->GetXaxis()->SetBinLabel(ibin,label);
      hbtags->GetXaxis()->SetBinLabel(ibin,label);
      hbtagsJP->GetXaxis()->SetBinLabel(ibin,label);
    } 

  mon.addHistogram( new TH1F( "mindphijmet",  ";min #Delta#phi(jet,E_{T}^{miss});Events",40,0,4) );
  mon.addHistogram( new TH1F( "mindphijmetNM1",  ";min #Delta#phi(jet,E_{T}^{miss});Events",40,0,4) );
  mon.addHistogram( new TH1D( "balance",      ";E_{T}^{miss}/q_{T};Events", 25,0,2.5) );
  mon.addHistogram( new TH1D( "balanceNM1",   ";E_{T}^{miss}/q_{T};Events", 25,0,2.5) );
  mon.addHistogram( new TH1F( "axialmet",     ";Axial missing transvere energy [GeV];Events", 50,-100,400) );
  mon.addHistogram( new TH1F( "axialmetNM1",   ";Axial missing transvere energy [GeV];Events", 50,-100,400) );
  Double_t metaxis[]={0,5,10,15,20,25,30,35,40,45,50,55,60,70,80,90,100,125,150,175,200,250,300,400,500};
  Int_t nmetAxis=sizeof(metaxis)/sizeof(Double_t);
  mon.addHistogram( new TH1F( "metpuppi",          ";Missing transverse energy [GeV];Events / GeV",nmetAxis-1,metaxis) ); //50,0,1000) );
  mon.addHistogram( new TH1F( "met",          ";Missing transverse energy [GeV];Events / GeV",nmetAxis-1,metaxis) ); //50,0,1000) );
  mon.addHistogram( new TH1F( "metNM1",        ";Missing transverse energy [GeV];Events / GeV",nmetAxis-1,metaxis) ); //50,0,1000) );
  Double_t mtaxis[]={100,120,140,160,180,200,220,240,260,280,300,325,350,375,400,450,500,600,700,800,900,1000,1500};
  Int_t nmtAxis=sizeof(mtaxis)/sizeof(Double_t);
  mon.addHistogram( new TH1F( "mt"  ,         ";Transverse mass [GeV];Events / GeV",nmtAxis-1,mtaxis) );
  mon.addHistogram( new TH1F( "mtNM1"  ,       ";Transverse mass [GeV];Events / GeV",nmtAxis-1,mtaxis) );
  mon.addHistogram( new TH1F( "mtresponse",   ";Transverse mass response [GeV];Events / GeV", 100,0,2) );
  mon.addHistogram( new TH1F( "mtcheckpoint"  ,         ";Transverse mass [GeV];Events",160,150,1750) );
  mon.addHistogram( new TH1F( "metcheckpoint" ,         ";Missing transverse energy [GeV];Events",100,0,500) );

  //Debug Plots Alessio
  mon.addHistogram( new TH1F(   "numbereeTrigger",    "Number of event passing the ee Trigger",  2, 0, 2) );
  mon.addHistogram( new TH1F( "numbermumuTrigger",  "Number of event passing the mumu Trigger",  2, 0, 2) );
  mon.addHistogram( new TH1F(  "numberemuTrigger",   "Number of event passing the emu Trigger",  2, 0, 2) );

  //
  // HISTOGRAMS FOR OPTIMIZATION and STATISTICAL ANALYSIS
  //
  //

  //NEED FOR ALL OPTIMIZATION
  TH1F* Hoptim_systs     =  (TH1F*) mon.addHistogram( new TH1F ("optim_systs"    , ";syst;", nvarsToInclude,0,nvarsToInclude) ) ;
  for(size_t ivar=0; ivar<nvarsToInclude; ivar++){
      Hoptim_systs->GetXaxis()->SetBinLabel(ivar+1, varNames[ivar]);
  }


  std::vector<double> optim_Cuts1_met;
  for(double met=50;met<140;met+=5) {  optim_Cuts1_met    .push_back(met);  }
  TH2F* Hoptim_cuts  =(TH2F*)mon.addHistogram(new TProfile2D("optim_cut",      ";cut index;variable",       optim_Cuts1_met.size(),0,optim_Cuts1_met.size(), 1, 0, 1)) ;
  Hoptim_cuts->GetYaxis()->SetBinLabel(1, "met>");
  for(unsigned int index=0;index<optim_Cuts1_met.size();index++){ Hoptim_cuts    ->Fill(index, 0.0, optim_Cuts1_met[index]);  }
     
  //##############################################
  //######## GET READY FOR THE EVENT LOOP ########
  //##############################################
  //MC normalization (to 1/pb)
  double xsecWeight = 1.0;
  if(isMC) xsecWeight=xsec/utils::getTotalNumberOfEvents(urls, false);//need to use the slow method in order to take NLO negative events into account

  //jet energy scale and uncertainties 
  TString jecDir = runProcess.getParameter<std::string>("jecDir");
  gSystem->ExpandPathName(jecDir);
  FactorizedJetCorrector *jesCor        = utils::cmssw::getJetCorrector(jecDir,isMC);
  JetCorrectionUncertainty *totalJESUnc = new JetCorrectionUncertainty((jecDir+"/MC_Uncertainty_AK5PFchs.txt").Data());
  
  //muon energy scale and uncertainties
  TString muscleDir = runProcess.getParameter<std::string>("muscleDir");
  gSystem->ExpandPathName(muscleDir);
  MuScleFitCorrector* muCor=NULL;//getMuonCorrector(muscleDir,dtag); //FIXME Not yet updated for run2

  //lepton efficiencies
  LeptonEfficiencySF lepEff;

  //b-tagging: beff and leff must be derived from the MC sample using the discriminator vs flavor
  //the scale factors are taken as average numbers from the pT dependent curves see:
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2012_Data_and_MC_EPS13_prescript
  BTagSFUtil btsfutil;
  float beff(0.68), sfb(0.99), sfbunc(0.015);
  float leff(0.13), sfl(1.05), sflunc(0.12);

  //pileup weighting
  edm::LumiReWeighting* LumiWeights = NULL;
  utils::cmssw::PuShifter_t PuShifters;
  double PUNorm[] = {1,1,1};
  if(isMC){
          std::vector<double> dataPileupDistributionDouble = runProcess.getParameter< std::vector<double> >("datapileup");
          std::vector<float> dataPileupDistribution; for(unsigned int i=0;i<dataPileupDistributionDouble.size();i++){dataPileupDistribution.push_back(dataPileupDistributionDouble[i]);}
          std::vector<float> mcPileupDistribution;

	  utils::getMCPileupDistributionFromMiniAOD(urls,dataPileupDistribution.size(), mcPileupDistribution);
          while(mcPileupDistribution.size()<dataPileupDistribution.size())  mcPileupDistribution.push_back(0.0);
          while(mcPileupDistribution.size()>dataPileupDistribution.size())dataPileupDistribution.push_back(0.0);
          gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE
          LumiWeights = new edm::LumiReWeighting(mcPileupDistribution,dataPileupDistribution);
          PuShifters=utils::cmssw::getPUshifters(dataPileupDistribution,0.05);
          utils::getPileupNormalization(mcPileupDistribution, PUNorm, LumiWeights, PuShifters);
  }

 
  gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE

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
  //DuplicatesChecker duplicatesChecker;
  //int nDuplicates(0)
  
  printf("Progressing Bar           :0%%       20%%       40%%       60%%       80%%       100%%\n");
  for(unsigned int f=0;f<urls.size();f++){
     TFile* file = TFile::Open(urls[f].c_str() );
     fwlite::Event ev(file);
     printf("Scanning the ntuple %2i/%2i :", (int)f+1, (int)urls.size());
     int iev=0;
     int treeStep(ev.size()/50);
     for(ev.toBegin(); !ev.atEnd(); ++ev){ iev++;
         if(iev%treeStep==0){printf(".");fflush(stdout);}
         float weight = xsecWeight;
         float shapeWeight = 1.0;
         double TotalWeight_plus = 1.0;
         double TotalWeight_minus = 1.0;
         float puWeight(1.0);

          //##############################################   EVENT LOOP STARTS   ##############################################
          //if(!isMC && duplicatesChecker.isDuplicate( ev.run, ev.lumi, ev.event) ) { nDuplicates++; continue; }

          //Skip bad lumi
          if(!goodLumiFilter.isGoodLumi(ev.eventAuxiliary().run(),ev.eventAuxiliary().luminosityBlock()))continue;

         //WEIGHT for Pileup
         if(isMC){          
	     int ngenITpu = 0;
	     fwlite::Handle< std::vector<PileupSummaryInfo> > puInfoH;
             puInfoH.getByLabel(ev, "slimmedAddPileupInfo");
             for(std::vector<PileupSummaryInfo>::const_iterator it = puInfoH->begin(); it != puInfoH->end(); it++){
                if(it->getBunchCrossing()==0)      { ngenITpu += it->getTrueNumInteractions(); } //getPU_NumInteractions(); }
             }
             puWeight          = LumiWeights->weight(ngenITpu) * PUNorm[0];
             TotalWeight_plus  = PuShifters[utils::cmssw::PUUP  ]->Eval(ngenITpu) * (PUNorm[2]/PUNorm[0]);
             TotalWeight_minus = PuShifters[utils::cmssw::PUDOWN]->Eval(ngenITpu) * (PUNorm[1]/PUNorm[0]);
             weight *= puWeight;
         }



          //apply trigger and require compatibilitiy of the event with the PD
          edm::TriggerResultsByName tr = ev.triggerResultsByName("HLT");
          if(!tr.isValid())return false;

         bool eeTrigger          = utils::passTriggerPatterns(tr, "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
         bool muTrigger          = utils::passTriggerPatterns(tr, "HLT_Mu34_TrkIsoVVL_v*");
         bool mumuTrigger        = utils::passTriggerPatterns(tr, "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*"); 
         bool emuTrigger         = utils::passTriggerPatterns(tr, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*", "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*");
         if(filterOnlyEE)   { mumuTrigger=false; emuTrigger=false;  }
         if(filterOnlyMUMU) { eeTrigger=false;   emuTrigger=false;  }
         if(isSingleMuPD)   { eeTrigger=false;   emuTrigger=false;  if( muTrigger && !mumuTrigger) mumuTrigger=true; else mumuTrigger=false; }
         if(filterOnlyEMU)  { eeTrigger=false;   mumuTrigger=false; }

          if(!(eeTrigger || mumuTrigger || emuTrigger))continue;  //ONLY RUN ON THE EVENTS THAT PASS OUR TRIGGERS
   
          //##############################################   EVENT PASSED THE TRIGGER   #######################################
          if( !isMC ){
             if( !metFiler.passMetFilter( ev, isPromptReco )) continue;
	  }
          //##############################################   EVENT PASSED MET FILTER   ####################################### 

          //load all the objects we will need to access
          reco::VertexCollection vtx;
          fwlite::Handle< reco::VertexCollection > vtxHandle; 
          vtxHandle.getByLabel(ev, "offlineSlimmedPrimaryVertices");
          if(vtxHandle.isValid()){ vtx = *vtxHandle;}

          double rho = 0;
          fwlite::Handle< double > rhoHandle;
          rhoHandle.getByLabel(ev, "fixedGridRhoFastjetAll");
          if(rhoHandle.isValid()){ rho = *rhoHandle;}

          reco::GenParticleCollection gen;
          fwlite::Handle< reco::GenParticleCollection > genHandle;
          genHandle.getByLabel(ev, "prunedGenParticles");
          if(genHandle.isValid()){ gen = *genHandle;}

          pat::MuonCollection muons;
          fwlite::Handle< pat::MuonCollection > muonsHandle;
          muonsHandle.getByLabel(ev, "slimmedMuons");
          if(muonsHandle.isValid()){ muons = *muonsHandle;}

          pat::ElectronCollection electrons;
          fwlite::Handle< pat::ElectronCollection > electronsHandle;
          electronsHandle.getByLabel(ev, "slimmedElectrons");
          if(electronsHandle.isValid()){ electrons = *electronsHandle;}

          pat::JetCollection jets;
          fwlite::Handle< pat::JetCollection > jetsHandle;
          jetsHandle.getByLabel(ev, "slimmedJets");
          if(jetsHandle.isValid()){ jets = *jetsHandle;}

          pat::PhotonCollection photons;
          fwlite::Handle< pat::PhotonCollection > photonsHandle;
          photonsHandle.getByLabel(ev, "slimmedPhotons");
          if(photonsHandle.isValid()){ photons = *photonsHandle;}
          
          pat::METCollection mets;
          fwlite::Handle< pat::METCollection > metsHandle;
          metsHandle.getByLabel(ev, "slimmedMETs");
          if(metsHandle.isValid()){ mets = *metsHandle;}
          LorentzVector met = mets[0].p4(); 

          pat::METCollection puppimets;
          fwlite::Handle< pat::METCollection > puppimetsHandle;
          puppimetsHandle.getByLabel(ev, "slimmedMETsPuppi");
          if(puppimetsHandle.isValid()){ puppimets = *puppimetsHandle;}
          LorentzVector puppimet = puppimets[0].p4(); 

          pat::TauCollection taus;
          fwlite::Handle< pat::TauCollection > tausHandle;
          tausHandle.getByLabel(ev, "slimmedTaus");
          if(tausHandle.isValid()){ taus = *tausHandle;}

         if(isV0JetsMC){
            fwlite::Handle< LHEEventProduct > lheEPHandle;
            lheEPHandle.getByLabel(ev, "externalLHEProducer");
            if(lheEPHandle.isValid()){
               mon.fillHisto("nup","",lheEPHandle->hepeup().NUP,1);
               if(lheEPHandle->hepeup().NUP>5) continue;
               mon.fillHisto("nupfilt","",lheEPHandle->hepeup().NUP,1);
            }else{
               printf("Handle to externalLHEProducer is invalid --> Can not ignore V0+Jet events from inclusive samples\n");
            }
         }


         //
         // DERIVE WEIGHTS TO APPLY TO SAMPLE
         //
         //
    
         if(isMC){
             //NLO weight:  This is needed because NLO generator might produce events with negative weights FIXME: need to verify that the total cross-section is properly computed
             fwlite::Handle< GenEventInfoProduct > genEventInfoHandle;
             genEventInfoHandle.getByLabel(ev, "generator");
             if(genEventInfoHandle.isValid()){ if(genEventInfoHandle->weight()<0){shapeWeight*=-1;}  }

     
           //final event weight
           weight *= shapeWeight;
         }

         //
         //
         // BELOW FOLLOWS THE ANALYSIS OF THE MAIN SELECTION WITH N-1 PLOTS
         //
         //

         //
         // PHOTON ANALYSIS
         //
         pat::PhotonCollection selPhotons;
         mon.fillHisto("npho", "trg", photons.size(), weight);
         for(size_t ipho=0; ipho<photons.size(); ipho++){
            pat::Photon photon = photons[ipho];
            mon.fillHisto("phopt", "trg", photon.pt(), weight);
            mon.fillHisto("phoeta", "trg", photon.eta(), weight);

            if(photon.pt()<20)continue;
            if(fabs(photon.superCluster()->eta())>1.4442 ) continue;
            if(!patUtils::passId(photon, rho, patUtils::llvvPhotonId::Tight)) continue;
            selPhotons.push_back(photon);
         }


         //
         // LEPTON ANALYSIS
         //
         
         //start by merging electrons and muons
         std::vector<patUtils::GenericLepton> leptons;
         for(size_t l=0;l<electrons.size();l++){leptons.push_back(patUtils::GenericLepton(electrons[l]));}      
         for(size_t l=0;l<muons    .size();l++){leptons.push_back(patUtils::GenericLepton(muons    [l]));}      
         std::sort(leptons.begin(),   leptons.end(), utils::sort_CandidatesByPt);

         LorentzVector muDiff(0,0,0,0); 
         std::vector<patUtils::GenericLepton> selLeptons, extraLeptons;
         for(size_t ilep=0; ilep<leptons.size(); ilep++)
           {
             bool passKin(true),passId(true),passIso(true);
             bool passLooseLepton(true), passSoftMuon(true), passSoftElectron(true), passVetoElectron(true);

             int lid=leptons[ilep].pdgId();

             //apply muon corrections
             if(abs(lid)==13)
               {
                 passSoftMuon=false;
                 if(muCor){
                   TLorentzVector p4(leptons[ilep].px(),leptons[ilep].py(),leptons[ilep].pz(),leptons[ilep].energy());
                   muCor->applyPtCorrection(p4 , lid<0 ? -1 :1 );
                   if(isMC) muCor->applyPtSmearing(p4, lid<0 ? -1 : 1, false);
                   muDiff -= leptons[ilep].p4();
                   leptons[ilep].setP4(LorentzVector(p4.Px(),p4.Py(),p4.Pz(),p4.E() ) );
                   muDiff += leptons[ilep].p4();
                 }
               }

             //no need for charge info any longer
             lid=abs(lid);
             TString lepStr( lid==13 ? "mu" : "e");

             //veto nearby photon (loose electrons are many times photons...)
             double minDRlg(9999.);
             for(size_t ipho=0; ipho<selPhotons.size(); ipho++)
               minDRlg=TMath::Min(minDRlg,deltaR(leptons[ilep].p4(),selPhotons[ipho].p4()));
             if(minDRlg<0.1) continue;

             //kinematics
             float leta = fabs(lid==11 ?  leptons[ilep].el.superCluster()->eta() : leptons[ilep].eta());
             if(leta> (lid==11 ? 2.5 : 2.4) )            passKin=false;
             if(lid==11 && (leta>1.4442 && leta<1.5660)) passKin=false;
             passLooseLepton &= passKin;
             passSoftMuon    &= passKin;
             if(lid==13){
               if(leptons[ilep].pt()<10) passLooseLepton=false;
               if(leptons[ilep].pt()<3)  passSoftMuon=false;
             }
             else if(lid==11){
               if(leptons[ilep].pt()<10) passLooseLepton=false;
             }
             if(leptons[ilep].pt()<20) passKin=false;

             //Cut based identification
         passId = lid==11?patUtils::passId(leptons[ilep].el, vtx[0], patUtils::llvvElecId::Tight) : patUtils::passId(leptons[ilep].mu, vtx[0], patUtils::llvvMuonId::Tight);
             passLooseLepton &= lid==11?patUtils::passId(leptons[ilep].el, vtx[0], patUtils::llvvElecId::Loose) : patUtils::passId(leptons[ilep].mu, vtx[0], patUtils::llvvMuonId::Loose);
         passSoftMuon &= lid==11? false : patUtils::passId(leptons[ilep].mu, vtx[0], patUtils::llvvMuonId::Soft);

             //isolation
         passIso = lid==11?patUtils::passIso(leptons[ilep].el,  patUtils::llvvElecIso::Tight) : patUtils::passIso(leptons[ilep].mu,  patUtils::llvvMuonIso::Tight);
         passLooseLepton &= lid==11?patUtils::passIso(leptons[ilep].el,  patUtils::llvvElecIso::Loose) : patUtils::passIso(leptons[ilep].mu,  patUtils::llvvMuonIso::Loose);

             if(passId && passIso && passKin)          selLeptons.push_back(leptons[ilep]); 
             else if(passLooseLepton || passSoftMuon)  extraLeptons.push_back(leptons[ilep]);

           }
           std::sort(selLeptons.begin(),   selLeptons.end(), utils::sort_CandidatesByPt);
           std::sort(extraLeptons.begin(), extraLeptons.end(), utils::sort_CandidatesByPt);
           LorentzVector recoMET = met - muDiff;
         
         //
         //JET/MET ANALYSIS
         //
         //add scale/resolution uncertainties and propagate to the MET      
         //utils::cmssw::updateJEC(jets,jesCor,totalJESUnc,rho,vtx.size(),isMC);  //FIXME if still needed
         //std::vector<LorentzVector> met=utils::cmssw::getMETvariations(recoMet,jets,selLeptons,isMC); //FIXME if still needed

         //select the jets
         pat::JetCollection selJets;
	 pat::JetCollection selBJets;
         int njets(0),nbjets(0),nbtags(0),nbtagsJP(0);
         float mindphijmet(9999.);
         for(size_t ijet=0; ijet<jets.size(); ijet++){
             if(jets[ijet].pt()<15 || fabs(jets[ijet].eta())>4.7 ) continue;

             //mc truth for this jet
             const reco::GenJet* genJet=jets[ijet].genJet();
             TString jetType( genJet && genJet->pt()>0 ? "truejetsid" : "pujetsid" );
             
             //cross-clean with selected leptons and photons
             double minDRlj(9999.),minDRlg(9999.);
             for(size_t ilep=0; ilep<selLeptons.size(); ilep++)
               minDRlj = TMath::Min( minDRlj, deltaR(jets[ijet],selLeptons[ilep]) );
             for(size_t ipho=0; ipho<selPhotons.size(); ipho++)
               minDRlg = TMath::Min( minDRlg, deltaR(jets[ijet],selPhotons[ipho]) );
             if(minDRlj<0.4 || minDRlg<0.4) continue;
             
             //jet id
             bool passPFloose = patUtils::passPFJetID("Loose", jets[ijet]);
             float PUDiscriminant = jets[ijet].userFloat("pileupJetId:fullDiscriminant");
             bool passLooseSimplePuId = patUtils::passPUJetID(jets[ijet]); //FIXME Broken in miniAOD V2 : waiting for JetMET fix. (Hugo)
             if(jets[ijet].pt()>30){
                 mon.fillHisto(jetType,"",fabs(jets[ijet].eta()),0);
                 if(passPFloose)                        mon.fillHisto(jetType,"",fabs(jets[ijet].eta()),1);
                 if(passLooseSimplePuId)                mon.fillHisto(jetType,"",fabs(jets[ijet].eta()),2);
                 if(passPFloose && passLooseSimplePuId) mon.fillHisto(jetType,"",fabs(jets[ijet].eta()),3);
             }
             if(!passPFloose /*|| !passLooseSimplePuId*/) continue; //FIXME PUJetID is broken in miniAOD V2 : waiting for JetMET fix (Hugo)
             selJets.push_back(jets[ijet]);

	     if(jets[ijet].pt()>20 && fabs(jets[ijet].eta())<2.4 && jets[ijet].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")>0.679){
	       selBJets.push_back(jets[ijet]);  
	       nbjets++;
	     }

             if(jets[ijet].pt()>30) {
               njets++;
               float dphijmet=fabs(deltaPhi(met.phi(), jets[ijet].phi()));
               if(dphijmet<mindphijmet) mindphijmet=dphijmet;
               if(fabs(jets[ijet].eta())<2.5){
                 bool hasCSVtag(jets[ijet].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")>0.423);
                 //update according to the SF measured by BTV
                 if(isMC){
                     int flavId=jets[ijet].partonFlavour();
                     if(abs(flavId)==5)        btsfutil.modifyBTagsWithSF(hasCSVtag,sfb,beff);
                     else if(abs(flavId)==4)   btsfutil.modifyBTagsWithSF(hasCSVtag,sfb/5,beff);
                     else		            btsfutil.modifyBTagsWithSF(hasCSVtag,sfl,leff);
                 }
                 
                 if( hasCSVtag ) nbtags++;
                 //nbtags   += hasCSVtag; 
               }
             }
           }
         std::sort(selJets.begin(), selJets.end(), utils::sort_CandidatesByPt);

         //select the taus
         pat::TauCollection selTaus;
         int ntaus(0);
         for(size_t itau=0; itau<taus.size(); ++itau){
           pat::Tau& tau = taus[itau];
           if(tau.pt()<20. || fabs(tau.eta()) >2.3) continue;
           
           //	bool overlapWithLepton(false);
           //	for(int l1=0; l1<(int)selLeptons.size();++l1){
           //	  if(deltaR(tau, selLeptons[l1])<0.1){overlapWithLepton=true; break;}
           //	}
           //	if(overlapWithLepton) continue;
           
           //	if(!tau.isPFTau()) continue; // Only PFTaus
           //	if(tau.emFraction() >=2.) continue;
           
           if(!tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"))continue;
           if(!tau.tauID("againstMuonTight3"))continue; 
           if(!tau.tauID("againstElectronMediumMVA5"))continue;
           
           selTaus.push_back(tau);
	   selLeptons.push_back(tau);
           ntaus++;
         }
         std::sort(selTaus.begin(), selTaus.end(), utils::sort_CandidatesByPt);
         

         //
         // ASSIGN CHANNEL
         //
         std::vector<TString> chTags;

	 LorentzVector leadingLep, trailerLep, zll;
         int dilId(1);
	 int dilLep1, dilLep2;
	 double BestMass;
	 bool passBestZmass;
	 zll = getZCand(selLeptons,leadingLep,trailerLep,passBestZmass,BestMass,dilLep1,dilLep2,dilId,rho);


         LorentzVector boson = zll;
         bool isDileptonCandidate = false;
         if(dilId!=-1){
             //check the channel
             if( abs(dilId)==121 && eeTrigger){   chTags.push_back("ll");   chTags.push_back("ee"); isDileptonCandidate=true; }
             if( abs(dilId)==169 && mumuTrigger){ chTags.push_back("ll"); chTags.push_back("mumu"); isDileptonCandidate=true; }
         }

         std::vector<TString> tags(1,"all");
         for(size_t ich=0; ich<chTags.size(); ich++){
           tags.push_back( chTags[ich] );
         }

         //////////////////////////
         //                      //
         //  BASELINE SELECTION  //
         //                      //
         //////////////////////////

//         bool passMass(fabs(boson.mass()-91)<15);
	 
	 //
	 // DILEPTON ANALYSIS
	 //
	 

	 /*
	 //apply data/mc correction factors
	 if(dilLep1>=0)weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[dilLep1].pt(), selLeptons[dilLep1].eta(), abs(selLeptons[dilLep1].id),  abs(selLeptons[dilLep1].id) ==11 ? "loose" : "loose" ).first : 1.0;
	 if(dilLep2>=0)weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[dilLep2].pt(), selLeptons[dilLep2].eta(), abs(selLeptons[dilLep2].id),  abs(selLeptons[dilLep2].id) ==11 ? "loose" : "loose" ).first : 1.0;
	 */	 

	 
	 if(!isDileptonCandidate) continue;
	 //if(examineThisEvent) cout << "dilepton candidate found" << endl; 
	 
	 bool passZmass = (passBestZmass && (fabs(zll.mass()-91.2)<15));
	 bool passZpt   = (zll.pt()>20);
         bool passMass = passZmass;

	 // if(examineThisEvent) cout << "zll.mass()-zll.pt(): " << zll.mass() << "-" << zll.pt() << endl; 
	 //if(examineThisEvent) cout << "passZbestMass-passZmass-passZpt: " << passBestZmass << "-" << passZmass << "-" << passZpt << endl; 
	 
         bool passQt(boson.pt()>55);

	 int higgsCandL1, higgsCandL2;
	 LorentzVector higgsCand;
	 int HiggsShortId, higgsCandId;
	 
	 higgsCand = getHiggsCand(selLeptons,dilLep1,dilLep2,higgsCandL1,higgsCandL2,higgsCandId,HiggsShortId,chTags);
	 
	 //check if the pair pass Lepton Veto
	 bool passLepVetoMain = true;
	 for(int l=0;l<(int)selTaus.size() && passLepVetoMain;l++){
	   if(l==dilLep1 || l==dilLep2 || l==higgsCandL1 || l==higgsCandL2) continue; //lepton already used in the dilepton pair or higgs candidate
	   if(abs(selLeptons[l].pdgId())==15){
	     pat::Tau&  tau = selTaus[l];
	     if(!tau.tauID("againstElectronLooseMVA5") ||
	 	!tau.tauID("againstMuonLoose3")    ||
	 	!tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits")  ) continue;                    
	     passLepVetoMain = false; break;
	   }else{
	     passLepVetoMain = false; break;
	   }
	 }
	 
	 //check if the pair pass b-jet veto
	 bool passBJetVetoMain = true;
	 for(int j1=0;j1<(int)selBJets.size();j1++){
	   if(deltaR(selBJets[j1], selLeptons[1])>0.4){passBJetVetoMain=false; break;}
	   if(deltaR(selBJets[j1], selLeptons[2])>0.4){passBJetVetoMain=false; break;}
	   if(higgsCandL1!=-1 && deltaR(selBJets[j1],selLeptons[higgsCandL1])>0.4){passBJetVetoMain=false; break;}
	   if(higgsCandL2!=-1 && deltaR(selBJets[j1],selLeptons[higgsCandL2])>0.4){passBJetVetoMain=false; break;}
	 }

	 //check how many additional light jets are present
	 int NCleanedJetMain = 0;
	 for(int j1=0;j1<(int)selJets.size();j1++){
	   if(dilLep1    !=-1 && deltaR(selJets[j1]   , selLeptons[dilLep1 ])<0.4) continue;
	   if(dilLep2    !=-1 && deltaR(selJets[j1]   , selLeptons[dilLep2 ])<0.4) continue;
	   if(higgsCandL1         !=-1 && deltaR(selJets[j1]   , selLeptons[higgsCandL1      ])<0.4) continue;
	   if(higgsCandL2         !=-1 && deltaR(selJets[j1]   , selLeptons[higgsCandL2      ])<0.4) continue;
	   NCleanedJetMain++;
	 }
	 
	 bool passDPhiCut    =  (fabs(deltaPhi(zll.phi(), met.phi()))>1.5);
	 bool passHiggsLoose = passHiggsCuts(selLeptons, rho, higgsCandId, higgsCandL1, higgsCandL2, 0.5, 0.5, "decayModeFinding", 0., false, vtx); 
	 bool passHiggsMain  = passHiggsCuts(selLeptons, rho, higgsCandId, higgsCandL1, higgsCandL2, 0.3, 0.3, "byLooseCombinedIsolationDeltaBetaCorr3Hits", 20., true, vtx);
	 
	 //SVFIT MASS
	 LorentzVector higgsCand_SVFit = higgsCand;
	 /*
	   commented for the time being
	 if(passZmass && passZpt && passDPhiCut && passHiggsLoose && passLepVetoMain && passBJetVetoMain){
	   higgsCand_SVFit = getSVFit(mets[0], selLeptons, higgsCandL1, higgsCandL2);  //compute svfit mass in a smart way
	   }
	 */
	 
	 //build the higgs candH
	 LorentzVector higgsCandH       = zll + higgsCand;
	 //	 LorentzVector higgsCandH_SVFit = zll + higgsCand_SVFit;
	 
         bool passThirdLeptonVeto( selLeptons.size()==2 && extraLeptons.size()==0 );
         bool passBtags(nbtags==0); 
         bool passMinDphijmet( njets==0 || mindphijmet>0.5);
         bool removeDump(false);


         //VBF Control plots to understand VBF Tail in Z mass shape
         std::vector<reco::GenParticle> VisLep;
         if(isMC_VBF1000){
		double filterEff = 0.16275;
                for( unsigned int k=0; k<gen.size(); ++k){	
			if( gen[k].isHardProcess() && ( abs( gen[k].pdgId() ) == 11 || abs( gen[k].pdgId() ) == 13 ) ) VisLep.push_back( gen[k] ); 
      		  }
		if( VisLep.size() == 2 ){
                        TLorentzVector Lep1( VisLep[0].px(), VisLep[0].py(), VisLep[0].pz(), VisLep[0].p() );
                        TLorentzVector Lep2( VisLep[1].px(), VisLep[1].py(), VisLep[1].pz(), VisLep[1].p() );	
			double MassVisZ = ( Lep1 + Lep2 ).M();
			if( MassVisZ > 150 ) removeDump=true; 
		}
                if(removeDump ) continue;
                weight /= filterEff;
	 }

	 
	 //
	 // NOW FOR THE CONTROL PLOTS
	 //
	 
	 // patUtils::GenericLepton leadingLep = selLeptons[0];
	 // patUtils::GenericLepton trailerLep = selLeptons[1];
	   
	 mon.fillHisto("eventflow2"       ,   tags,                 0, weight);
	 if(selLeptons.size()>=2){
	   mon.fillHisto("nlep"           ,   chTags, selLeptons.size(), weight);
	   mon.fillHisto("eventflow2"     ,   tags,                 1, weight);
	   mon.fillHisto("zllmass"          ,   tags, boson.mass(),    weight);
	   if(passZmass){
	     mon.fillHisto("eventflow2"   ,   tags,                 2, weight);
	     //pu control
	     mon.fillHisto("nvtx2"        ,   tags, vtx.size(),      weight);
	     mon.fillHisto("nvtxraw2"     ,   tags, vtx.size(),      weight/puWeight);
	     mon.fillHisto("rho2"         ,   tags, rho,       weight);
	     
	     //Z kinematics control
	     mon.fillHisto("leadpt2"      ,   tags, leadingLep.pt(), weight);      
	     mon.fillHisto("leadeta2"     ,   tags, leadingLep.eta(), weight);      
	     mon.fillHisto("trailerpt2"   ,   tags, trailerLep.pt(), weight);      
	     mon.fillHisto("trailereta2"  ,   tags, trailerLep.eta(), weight);      
	     mon.fillHisto("leppt2"       ,   tags, leadingLep.pt(), weight);      
	     mon.fillHisto("leppt2"       ,   tags, trailerLep.pt(), weight);      
	     mon.fillHisto("lepeta2"      ,   tags, leadingLep.eta(), weight);      
	     mon.fillHisto("lepeta2"      ,   tags, trailerLep.eta(), weight);      
	     
	     //analyze dilepton kinematics
	     mon.fillHisto("zllpt"         ,   tags, zll.pt(),      weight);      
	     mon.fillHisto("zlleta"        ,   tags, zll.eta(),     weight);
	     mon.fillHisto("zlly"          ,   tags, zll.Rapidity(),weight);
	     
	     if(passZpt){
	       mon.fillHisto("eventflow2",   tags,                 3, weight);
	       
	       mon.fillHisto("ntaus"           ,  chTags, selTaus.size(), weight);
	       mon.fillHisto("tauleadpt"       ,  tags,   selTaus.size()>0?selTaus[0].pt():-1,  weight);
	       mon.fillHisto("tauleadeta"      ,  tags,   selTaus.size()>0?selTaus[0].eta():-10, weight);
	       mon.fillHisto("tautrailerpt"    ,  tags,   selTaus.size()>1?selTaus[1].pt():-1,  weight);
	       mon.fillHisto("tautrailereta"   ,  tags,   selTaus.size()>1?selTaus[1].eta():-10, weight);
	       mon.fillHisto("taupt"           ,  chTags, selTaus.size()>0?selTaus[0].pt():-1, weight);
	       mon.fillHisto("taupt"           ,  chTags, selTaus.size()>0?selTaus[1].pt():-1, weight);
	       mon.fillHisto("taueta"          ,  tags,   selTaus.size()>0?selTaus[0].eta():-10, weight);
	       mon.fillHisto("taueta"          ,  tags,   selTaus.size()>0?selTaus[0].eta():-10, weight);
	       
	       if(selLeptons.size()>=4){
		 mon.fillHisto("eventflow2",   tags,                 4, weight);
		 if(passLepVetoMain){
		   mon.fillHisto("eventflow2", tags,                 5, weight);
		   mon.fillHisto("nbjets"    , chTags, nbjets,  weight);
		   mon.fillHisto("njets2"     , chTags, njets,   weight);
		   
		   if(passBJetVetoMain){
		     mon.fillHisto("eventflow2"	,   tags,                 6, weight);
		     
		     mon.fillHisto("dPhi_AZ"    , tags, deltaPhi(higgsCand.phi(), boson.phi()),    weight);
		     mon.fillHisto("dPhi_AMet"  , tags, deltaPhi(higgsCand.phi(), met.phi()),    weight);
		     mon.fillHisto("dPhi_ZMet"  , tags, deltaPhi(boson.phi(), met.phi()),    weight);
		     mon.fillHisto("met2"      	, tags, met.pt()         , weight);
		     
		     if(passDPhiCut){
		       mon.fillHisto("eventflow2",   tags,                 7, weight);
		       if(passHiggsLoose){
			 mon.fillHisto("sumpt",   tags, selLeptons[higgsCandL1].pt()+selLeptons[higgsCandL2].pt(), weight);
			 if(passHiggsMain){
			   mon.fillHisto("eventflow2"   ,tags,                 8, weight);
			   mon.fillHisto("yields"	,tags,                HiggsShortId, weight);
			   mon.fillHisto("yieldsOS"     ,tags,                HiggsShortId, weight);
			   
			   mon.fillHisto("Apt"       	, tags, higgsCand.pt(),    weight);
			   mon.fillHisto("Amass"           , tags, higgsCand.mass(),  weight);
			   mon.fillHisto("Amasssvfit"      , tags, higgsCand_SVFit.mass(),  weight);
			   mon.fillHisto("Hmass"           , tags, higgsCandH.mass(),  weight);
			   mon.fillHisto("Hpt"             , tags, higgsCandH.pt(),  weight);
			   // mon.fillHisto("Hmasssvfit"   , tags, higgsCandH_SVFit.mass(),  weight);
			   
			   mon.fillHisto("Anjets"    	, tags, NCleanedJetMain      , weight); 
			   mon.fillHisto("Amet"      	, tags, met.pt()         , weight);
			 } 
		       }
		     }
		   }
		 }
	       }
	     }
	   }  
	 }  
	 
         mon.fillHisto("eventflow",  tags,0,weight);
         mon.fillHisto("nvtxA",  tags,vtx.size(),1);
         mon.fillHisto("nvtxB",  tags,vtx.size(),xsecWeight);
         mon.fillHisto("nvtxC",  tags,vtx.size(),xsecWeight * puWeight);
         mon.fillHisto("nvtxD",  tags,vtx.size(),xsecWeight * puWeight * shapeWeight);

         mon.fillHisto("nvtxraw",  tags,vtx.size(),xsecWeight * shapeWeight);
         mon.fillHisto("nvtx",  tags,vtx.size(),weight);
         mon.fillHisto("rho",  tags,rho,weight);

         if(chTags.size()==0) continue;
         mon.fillHisto("eventflow",  tags,1,weight);
         mon.fillHisto("leadpt",      tags,selLeptons[0].pt(),weight); 
         mon.fillHisto("trailerpt",   tags,selLeptons[1].pt(),weight); 
         mon.fillHisto("leadeta",     tags,fabs(selLeptons[0].eta()),weight); 
         mon.fillHisto("trailereta",  tags,fabs(selLeptons[1].eta()),weight); 
         mon.fillHisto("ntaus", tags, ntaus,weight);
         if(ntaus>0) mon.fillHisto("leadtaupt", tags, selTaus[0].pt(),weight);
  
         mon.fillHisto("zmass", tags,boson.mass(),weight); 
         mon.fillHisto("zy",    tags,fabs(boson.Rapidity()),weight); 

         if(passMass){

           mon.fillHisto("eventflow",tags, 2,weight);
           mon.fillHisto("zpt",      tags, boson.pt(),weight);
           mon.fillHisto("zpt_rebin",tags, boson.pt(),weight,true);

           //these two are used to reweight photon -> Z, the 3rd is a control
           mon.fillHisto("qt",       tags, boson.pt(),weight,true); 
           mon.fillHisto("qtraw",    tags, boson.pt(),weight,true); 

           if(passQt){
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
             if(passThirdLeptonVeto){
               
               mon.fillHisto("eventflow",tags,4,weight);
               for(size_t ijet=0; ijet<selJets.size(); ijet++){
                 if(selJets[ijet].pt()<30 || fabs(selJets[ijet].eta())>2.5) continue;

                 float csv(selJets[ijet].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
                 mon.fillHisto( "csv",tags,csv,weight);
                 if(!isMC) continue;
                 int flavId=selJets[ijet].partonFlavour();
                 TString jetFlav("others");
                 if(abs(flavId)==5)      jetFlav="b";
                 else if(abs(flavId)==4) jetFlav="c";
                 mon.fillHisto( "csv"+jetFlav,tags,csv,weight);
               }
               mon.fillHisto( "nbtags",tags,nbtags,weight);
               mon.fillHisto( "nbtagsJP",tags,nbtagsJP,weight);
               
               if(passBtags){
                 mon.fillHisto("eventflow",tags,5,weight);

                 //include photon prediction from this point forward
                 //requires looping tag by tag as weights are category-specific
                 //the following relies on the hypothesis that the tags are ordered as follows: all, ch, ch+subtag, ch, ch+subtag, etc...
                 //so that the ch will be assigned the weight of its subtag and all will be the summ of all ch+subtag weights
                 std::vector<LorentzVector> massiveBoson(tags.size(),boson);
                 std::vector<float> photonWeights(tags.size(),1.0);

                 for(size_t itag=0; itag<tags.size(); itag++){		
                   //update the weight
                   TString icat=tags[itag];
                   float iweight(weight*photonWeights[itag]);
                   
                   LorentzVector iboson=massiveBoson[itag];

                   mon.fillHisto( "mindphijmet",icat,mindphijmet,iweight);
                   if(met.pt()>80) mon.fillHisto( "mindphijmetNM1",icat,mindphijmet,iweight);
                   if(passMinDphijmet){
                     mon.fillHisto("eventflow",icat,6,iweight);
                     
                     //this one is used to sample the boson mass: cuts may shape Z lineshape
                     mon.fillHisto("qmass",       tags, boson.mass(),weight); 
                     mon.fillHisto( "njets",icat,njets,iweight);

                     double b_dphi=fabs(deltaPhi(iboson.phi(),met.phi()));
                     mon.fillHisto( "metphi",icat,met.phi(),iweight,true);                                                                    
                     mon.fillHisto( "bosonphi",icat,iboson.phi(),iweight,true);                                                               
                     mon.fillHisto( "dphi_boson_met",icat,b_dphi,iweight,true);

                     if( isMC ) { 
			mon.fillHisto( "met",icat,met.pt(),iweight,true);
                     	mon.fillHisto( "metpuppi",icat,puppimet.pt(),iweight,true);
                    	mon.fillHisto( "balance",icat,met.pt()/iboson.pt(),iweight);
                     } else if ( !isMC && blindData ){
                        if( met.pt() < 100 ){ 
			  mon.fillHisto( "met",icat,met.pt(),iweight,true);
			  mon.fillHisto( "balance",icat,met.pt()/iboson.pt(),iweight);
                        } 
                        if( puppimet.pt() < 100 ){ mon.fillHisto( "metpuppi",icat,puppimet.pt(),iweight,true);}
                     }

                     TVector2 met2(met.px(),met.py());
                     TVector2 boson2(iboson.px(), iboson.py());
                     double axialMet(boson2*met2); axialMet/=-iboson.pt();
                     mon.fillHisto( "axialmet",icat,axialMet,iweight);
                     double mt=higgs::utils::transverseMass(iboson,met,true);

                     if( isMC ){
			mon.fillHisto( "mt",icat,mt,iweight,true);
		     } else if( !isMC && blindData && mt<325 ){
                        mon.fillHisto( "mt",icat,mt,iweight,true);               
                     }

                     if(met.pt()>optim_Cuts1_met[0]) 
                       {
                         mon.fillHisto( "mtcheckpoint",  icat, mt,       iweight, true);
                         mon.fillHisto( "metcheckpoint", icat, met.pt(), iweight, true);
                       }

                     if(met.pt()>80){
                       mon.fillHisto("eventflow",icat,7,iweight);
                       mon.fillHisto( "mtNM1",icat,mt,iweight,true);
                       mon.fillHisto( "balanceNM1",icat,met.pt()/iboson.pt(),iweight);
                       mon.fillHisto( "axialmetNM1",icat,axialMet,iweight);
                     }
                     if(mt>500){
                       mon.fillHisto( "metNM1",icat,met.pt(),iweight,true);
                     }

                     //pre-VBF control
                     if(njets>=2){
                       LorentzVector dijet=selJets[0].p4()+selJets[1].p4();
                       float deta=fabs(selJets[0].eta()-selJets[1].eta());
                       float dphi=fabs(deltaPhi(selJets[0].phi(),selJets[1].phi()));
                       float pt1(selJets[0].pt()),pt2(selJets[1].pt());
                       mon.fillHisto( "leadjetpt",icat,pt1,iweight);
                       mon.fillHisto( "trailerjetpt",icat,pt2,iweight);
                       if(pt1>30 && pt2>30){
                         float eta1(selJets[0].eta()),eta2(selJets[1].eta());
                         float fwdEta( fabs(eta1)>fabs(eta2) ? eta1 : eta2);
                         float cenEta( fabs(eta1)>fabs(eta2) ? eta2 : eta1);
                         mon.fillHisto("fwdjeteta",icat,fabs(fwdEta),  iweight);
                         mon.fillHisto("cenjeteta",icat,fabs(cenEta),  iweight);
                         mon.fillHisto("vbfdetajj",icat,deta,        iweight);
                         if(deta>4.0){
                           mon.fillHisto("vbfmjj",   icat,dijet.mass(),iweight,true);
                           if(dijet.mass()>500){
                             mon.fillHisto("vbfdphijj",icat,dphi,        iweight);
                           }
                         }
                       }
                     }
                   }
                 }
               }
             }
           }        
         }


         //
         // HISTOS FOR STATISTICAL ANALYSIS (include systematic variations)
         //
         //Fill histogram for posterior optimization, or for control regions
         for(size_t ivar=0; ivar<nvarsToInclude; ivar++){
           float iweight = weight;                            //nominal
           
           //energy scale/resolution
           bool varyJesUp( varNames[ivar]=="_jesup" );
           bool varyJesDown( varNames[ivar]=="_jesdown" );
           bool varyJerUp( varNames[ivar]=="_jerup" );
           bool varyJerDown( varNames[ivar]=="_jerdown" );
           bool varyUmetUp( varNames[ivar]=="_umetup" );
           bool varyUmetDown( varNames[ivar]=="_umetdown" );
           bool varyLesUp( varNames[ivar]=="_lesup" );
           bool varyLesDown( varNames[ivar]=="_lesdown" );
                   
           //pileup variations
           if(varNames[ivar]=="_puup") iweight *=TotalWeight_plus;
           if(varNames[ivar]=="_pudown") iweight *=TotalWeight_minus;
           
           //btag
           bool varyBtagUp( varNames[ivar]=="_btagup" );
           bool varyBtagDown( varNames[ivar]=="_btagdown" );
           
           //Q^2 variations on VV pT spectum
           if( ( (isMC_ZZ && (varNames[ivar]=="_zzptup" || varNames[ivar]=="_zzptdown")) || (isMC_WZ && (varNames[ivar]=="_wzptup" || varNames[ivar]=="_wzptdown") ) ) && vvShapeUnc.size()==2 )
             {
               size_t idx( varNames[ivar].EndsWith("up") ? 0 : 1 );
               TGraph *varGr=vvShapeUnc[idx];
               if(varGr==0) continue;
               std::vector<LorentzVector> vs;
               for(size_t ipart=0; ipart<gen.size(); ipart++)
                 {
                   int status=gen[ipart].status();
                   if(status!=3) continue;
                   int pid=gen[ipart].pdgId();
                   if(abs(pid)!=23 && abs(pid)!=24) continue;
                   vs.push_back( gen[ipart].p4() );
                 }
               if(vs.size()==2)
                 {
                   LorentzVector vv=vs[0]+vs[1];
                   iweight *= varGr->Eval(vv.pt());
                 }
             }
           

           //recompute MET/MT if JES/JER was varied
           LorentzVector    zvv = mets[0].p4();

           //FIXME
           if(varyJesUp)    zvv = mets[0].shiftedP4(pat::MET::METUncertainty::JetEnUp);
           if(varyJesDown)  zvv = mets[0].shiftedP4(pat::MET::METUncertainty::JetEnDown);
           if(varyJerUp)    zvv = mets[0].shiftedP4(pat::MET::METUncertainty::JetResUp);
           if(varyJerDown)  zvv = mets[0].shiftedP4(pat::MET::METUncertainty::JetResDown);
           if(varyUmetUp)   zvv = mets[0].shiftedP4(pat::MET::METUncertainty::UnclusteredEnUp);
           if(varyUmetDown) zvv = mets[0].shiftedP4(pat::MET::METUncertainty::UnclusteredEnDown);
           //if(varyLesUp)    zvv = met[utils::cmssw::LESUP]; //FIXME  must vary all leptons separately: MuonEnUp/MuonEnDown/ElectronEnUp/ElectronEnDown/TauEnUp/TauEnDown
           //if(varyLesDown)  zvv = met[utils::cmssw::LESDOWN];

           pat::JetCollection tightVarJets;
           int nbtagsjets = 0; 
           for(size_t ijet=0; ijet<jets.size(); ijet++){

             float eta=jets[ijet].eta();
             if( fabs(eta)>4.7 ) continue;
             float pt=jets[ijet].pt();

             //FIXME
             /*if(varyJesUp)    pt=jets[ijet].getVal("jesup");
             if(varyJesDown)  pt=jets[ijet].getVal("jesdown");
             if(varyJerUp)    pt=jets[ijet].getVal("jerup");
             if(varyJerDown)  pt=jets[ijet].getVal("jerdown");*/

             if( pt < 30 ) continue;
    
             //cross-clean with selected leptons and photons
             double minDRlj(9999.),minDRlg(9999.);
             for(size_t ilep=0; ilep<selLeptons.size(); ilep++)
               minDRlj = TMath::Min( minDRlj, deltaR(jets[ijet].p4(),selLeptons[ilep].p4()) );
             for(size_t ipho=0; ipho<selPhotons.size(); ipho++)
               minDRlg = TMath::Min( minDRlg, deltaR(jets[ijet].p4(),selPhotons[ipho].p4()) );
             if(minDRlj<0.4 || minDRlg<0.4) continue;
             
             //jet id
             bool         passPFloose = patUtils::passPFJetID("Loose", jets[ijet]);
             float     PUDiscriminant = jets[ijet].userFloat("pileupJetId:fullDiscriminant");
             bool passLooseSimplePuId = patUtils::passPUJetID(jets[ijet]); //Uses recommended value of HZZ, will update this as soon my analysis is done. (Hugo)
             if(!passPFloose || !passLooseSimplePuId) continue;
            
             //jet is selected
             tightVarJets.push_back(jets[ijet]);

             //check b-tag
             if( pt < 30 || fabs(eta) > 2.5 ) continue;
             if(!isMC && ivar>0 ) continue;
             //if(!varyBtagUp && !varyBtagDown) continue;
             int flavId=jets[ijet].partonFlavour();
             bool hasCSVtag (jets[ijet].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")>0.423);

             if(abs(flavId)==5)        btsfutil.modifyBTagsWithSF(hasCSVtag,sfb,beff);
             else if(abs(flavId)==4)   btsfutil.modifyBTagsWithSF(hasCSVtag,sfb/5,beff);
             else		            btsfutil.modifyBTagsWithSF(hasCSVtag,sfl,leff);

             /*if(varyBtagUp) {
               if(abs(flavId)==5)        btsfutil.modifyBTagsWithSF(hasCSVtag,sfb+sfbunc,beff);
               else if(abs(flavId)==4)   btsfutil.modifyBTagsWithSF(hasCSVtag,sfb/5+2*sfbunc,beff);
               else		      btsfutil.modifyBTagsWithSF(hasCSVtag,sfl+sflunc,leff);
             }
             else if(varyBtagDown) {
               if(abs(flavId)==5)        btsfutil.modifyBTagsWithSF(hasCSVtag,sfb-sfbunc,beff);
               else if(abs(flavId)==4)   btsfutil.modifyBTagsWithSF(hasCSVtag,sfb/5-2*sfbunc,beff);
               else		      btsfutil.modifyBTagsWithSF(hasCSVtag,sfl-sflunc,leff);
             }*/

             if(hasCSVtag) nbtagsjets++;
           }

           bool passLocalBveto( nbtagsjets == 0 );	
           bool isZsideBand    ( (boson.mass()>40  && boson.mass()<70) || (boson.mass()>110 && boson.mass()<200) );
           bool isZsideBandPlus( (boson.mass()>110 && boson.mass()<200) );
           bool passPreselection                 (passMass && passQt && passThirdLeptonVeto && passMinDphijmet && passLocalBveto);
           bool passPreselectionMbvetoMzmass     (            passQt && passThirdLeptonVeto && passMinDphijmet                  );         
         
           for(size_t ich=0; ich<chTags.size(); ich++){

             TString tags_full=chTags[ich]; 
             float chWeight(iweight);

             //update weight and mass for photons
             LorentzVector iboson(boson);
             
             //updet the transverse mass
             float mt =higgs::utils::transverseMass(iboson,zvv,true);

             //scan the MET cut and fill the shapes
             for(unsigned int index=0;index<optim_Cuts1_met.size();index++){             
               
               if(zvv.pt()>optim_Cuts1_met[index]){              
                   
                   float nrweight=chWeight;
          
                   if(passPreselection && ivar==0                                              )   mon.fillHisto("metcount", tags_full, index, nrweight);
                   if(passPreselection                                                         )   mon.fillHisto(TString("mt_shapes")+varNames[ivar],tags_full,index, mt,nrweight);
                   if(passPreselection                                                         )   mon.fillHisto(TString("met_shapes")+varNames[ivar],tags_full,index, zvv.pt(),nrweight);                    
                   if(passPreselectionMbvetoMzmass && passMass          && passLocalBveto      )   mon.fillHisto("mt_shapes_NRBctrl"+varNames[ivar],tags_full,index,0,nrweight);
                   if(passPreselectionMbvetoMzmass && isZsideBand       && passLocalBveto      )   mon.fillHisto("mt_shapes_NRBctrl"+varNames[ivar],tags_full,index,1,nrweight);
                   if(passPreselectionMbvetoMzmass && isZsideBandPlus   && passLocalBveto      )   mon.fillHisto("mt_shapes_NRBctrl"+varNames[ivar],tags_full,index,2,nrweight);
                   if(passPreselectionMbvetoMzmass && passMass          && !passLocalBveto     )   mon.fillHisto("mt_shapes_NRBctrl"+varNames[ivar],tags_full,index,3,nrweight);
                   if(passPreselectionMbvetoMzmass && isZsideBand       && !passLocalBveto     )   mon.fillHisto("mt_shapes_NRBctrl"+varNames[ivar],tags_full,index,4,nrweight);
                   if(passPreselectionMbvetoMzmass && isZsideBandPlus   && !passLocalBveto     )   mon.fillHisto("mt_shapes_NRBctrl"+varNames[ivar],tags_full,index,5,nrweight);

                 
               }
             }
           }
         }
     }
     printf("\n"); 
     delete file;
  } 
  //##############################################
  //########     SAVING HISTO TO FILE     ########
  //##############################################
  //save control plots to file
  printf("Results save in %s\n", outUrl.Data());
  
  //save all to the file
  TFile *ofile=TFile::Open(outUrl, "recreate");
  mon.Write();
  ofile->Close();

  if(outTxtFile)fclose(outTxtFile);

  //Now that everything is done, dump the list of lumiBlock that we processed in this job
  if(!isMC){
     goodLumiFilter.FindLumiInFiles(urls);
     goodLumiFilter.DumpToJson(((outUrl.ReplaceAll(".root",""))+".json").Data());
  }
}  

