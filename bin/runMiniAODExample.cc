#include <iostream>
#include <boost/shared_ptr.hpp>

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

//Load here all the dataformat that we will need
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
//#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h" //for svfit

#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/llvv_fwk/interface/llvvObjects.h"
#include "UserCode/llvv_fwk/interface/TMVAUtils.h"
#include "UserCode/llvv_fwk/interface/LeptonEfficiencySF.h"
#include "UserCode/llvv_fwk/interface/PDFInfo.h"
#include "UserCode/llvv_fwk/interface/MuScleFitCorrector.h"
#include "UserCode/llvv_fwk/interface/GammaWeightsHandler.h"
//#include "UserCode/llvv_fwk/interface/i2HDMUtils.h"

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

	std::vector<std::string> urls=runProcess.getParameter<std::vector<std::string> >("input");
	TString baseDir    = runProcess.getParameter<std::string>("dirName");
	TString url = TString(argv[1]);
	TString outFileUrl(gSystem->BaseName(url));
	outFileUrl.ReplaceAll(".py","");
	if(mctruthmode!=0) { outFileUrl += "_filt"; outFileUrl += mctruthmode; }
	TString outdir=runProcess.getParameter<std::string>("outdir");
	TString outUrl( outdir );
	gSystem->Exec("mkdir -p " + outUrl);
//	bool filterOnlyEE(false), filterOnlyMUMU(false);
//	if(!isMC)
//	{
//		if(url.Contains("DoubleEle")) filterOnlyEE=true;
//		if(url.Contains("DoubleMu"))  filterOnlyMUMU=true;
//	}
//	bool isSingleMuPD(!isMC && url.Contains("SingleMu"));  
//	bool isV0JetsMC(isMC && (url.Contains("DYJetsToLL_50toInf") || url.Contains("WJets")));
	bool isSignal(isMC && (url.Contains("VBFNLO") || url.Contains("lljj")) );

	TString outTxtUrl= outUrl + "/" + outFileUrl + ".txt";
	FILE* outTxtFile = NULL;
	if(!isMC)outTxtFile = fopen(outTxtUrl.Data(), "w");
	printf("TextFile URL = %s\n",outTxtUrl.Data());

	//lepton efficienciesAOA
	LeptonEfficiencySF lepEff;

	//systematics
	bool runSystematics                        = runProcess.getParameter<bool>("runSystematics");
	std::vector<TString> varNames(1,"");
	size_t nvarsToInclude(1);
	if(runSystematics && isMC)
	{
		//      varNames.push_back("_jerup"); varNames.push_back("_jerdown");
		//      varNames.push_back("_jesup"); varNames.push_back("_jesdown");
		//      varNames.push_back("_puup"); varNames.push_back("_pudown");
		if(isSignal)
		{
			//	  varNames.push_back("_q2up"); varNames.push_back("_q2down");
			//	  varNames.push_back("_pdfup"); varNames.push_back("_pdfdown");
			//	  varNames.push_back("_balanceup"); varNames.push_back("_balancedown");
		}
		nvarsToInclude=varNames.size();
		cout << nvarsToInclude << " systematics will be computed for this analysis" << endl;
	}

	//##############################################
	//########    INITIATING HISTOGRAMS     ########
	//##############################################
	SmartSelectionMonitor mon;

	//pileup control
	mon.addHistogram( new TH1F( "nvtx",";Vertices;Events",50,-0.5,49.5) ); 
	mon.addHistogram( new TH1F( "nvtxraw",";Vertices;Events",50,-0.5,49.5) ); 

	//##############################################
	//######## GET READY FOR THE EVENT LOOP ########
	//##############################################

	fwlite::ChainEvent ev(urls);
	const unsigned int totalEntries= ev.size();

	//MC normalization (to 1/pb)
	double xsecWeight = xsec/totalEntries;
	if(!isMC) xsecWeight=1.0;

	//pileup weighting: based on vtx for now...
	edm::LumiReWeighting* LumiWeights = NULL;
	utils::cmssw::PuShifter_t PuShifters;
	double PUNorm[] = {1,1,1};
	if(isMC){
		std::vector<double> dataPileupDistributionDouble = runProcess.getParameter< std::vector<double> >("datapileup");
		std::vector<float> dataPileupDistribution; for(unsigned int i=0;i<dataPileupDistributionDouble.size();i++){dataPileupDistribution.push_back(dataPileupDistributionDouble[i]);}
		std::vector<float> mcPileupDistribution;
		utils::getMCPileupDistributionFromMiniAOD(ev,dataPileupDistribution.size(), mcPileupDistribution);
		while(mcPileupDistribution.size()<dataPileupDistribution.size())  mcPileupDistribution.push_back(0.0);
		while(mcPileupDistribution.size()>dataPileupDistribution.size())dataPileupDistribution.push_back(0.0);
		LumiWeights= new edm::LumiReWeighting(mcPileupDistribution,dataPileupDistribution);
		PuShifters=utils::cmssw::getPUshifters(dataPileupDistribution,0.05);
		utils::getPileupNormalization(mcPileupDistribution, PUNorm, LumiWeights, PuShifters);
	}
	gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE

//	FILE* outTxtEvents = NULL;
//	outTxtEvents = fopen(outTxtUrl.Data(), "w");

	//##############################################
	//########           EVENT LOOP         ########
	//##############################################
	//loop on all the events
	printf("Progressing Bar     :0%%       20%%       40%%       60%%       80%%       100%%\n");
	printf("Scanning the ntuple :");
	unsigned int step(totalEntries/50);
	for( unsigned int iev=0; iev<totalEntries; iev++ ){
		if(iev%step==0){printf(".");fflush(stdout);}

		//##############################################   EVENT LOOP STARTS   ##############################################
		ev.to(iev); //load the event content from the EDM file

                reco::VertexCollection vtx;
		fwlite::Handle< reco::VertexCollection > vtxHandle;
		vtxHandle.getByLabel(ev, "offlineSlimmedPrimaryVertices");
		if(vtxHandle.isValid()){ vtx = *vtxHandle;}
                int nvtx = vtx.size();


                fwlite::Handle< std::vector<PileupSummaryInfo> > puInfoH;
                puInfoH.getByLabel(ev, "addPileupInfo");
                if(!puInfoH.isValid()){printf("collection PileupSummaryInfos with name addPileupInfo does not exist\n"); exit(0);}
                int ngenITpu = 0;
                for(std::vector<PileupSummaryInfo>::const_iterator it = puInfoH->begin(); it != puInfoH->end(); it++){
                   if(it->getBunchCrossing()==0)      { ngenITpu += it->getPU_NumInteractions(); }
                }



/*
		//get the collection of generated Particles
		fwlite::Handle< llvvGenEvent > genEventHandle;
		genEventHandle.getByLabel(ev, "llvvObjectProducersUsed");
		if(!genEventHandle.isValid()){cout<<"llvvGenEvent Object NotFound\n"; continue;}
		llvvGenEvent genEv = *genEventHandle;

		if(isV0JetsMC){ //drop V+1,2,3,4Jets part of V+0Jets in order to avoid double counting
			mon.fillHisto("nup","",genEv.nup,1);
			if(genEv.nup>5) continue;
			mon.fillHisto("nupfilt","",genEv.nup,1);
		}

		fwlite::Handle< llvvGenParticleCollection > genPartCollHandle;
		genPartCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
		if(!genPartCollHandle.isValid()){printf("llvvGenParticleCollection Object NotFound\n");  continue;}
		llvvGenParticleCollection gen = *genPartCollHandle;

		fwlite::Handle< llvvLeptonCollection > leptonCollHandle;
		leptonCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
		if(!leptonCollHandle.isValid()){printf("llvvLeptonCollection Object NotFound\n"); continue;}
		llvvLeptonCollection leptons = *leptonCollHandle;

		fwlite::Handle< llvvElectronInfoCollection > electronInfoCollHandle;
		electronInfoCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
		if(!electronInfoCollHandle.isValid()){printf("llvvElectronInfoCollection Object NotFound\n");  continue;}

		fwlite::Handle< llvvMuonInfoCollection > muonInfoCollHandle;
		muonInfoCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
		if(!muonInfoCollHandle.isValid()){printf("llvvMuonInfoCollection Object NotFound\n");  continue;}

		fwlite::Handle< llvvTauCollection > tauCollHandle;
		tauCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
		if(!tauCollHandle.isValid()){printf("llvvTauCollection Object NotFound\n");  continue;}
		llvvTauCollection taus = *tauCollHandle;

		fwlite::Handle< llvvTauCollection > boostedtauCollHandle;
		boostedtauCollHandle.getByLabel(ev, "llvvObjectProducersUsed", "boosted");
		if(!boostedtauCollHandle.isValid()){printf("llvvTauCollection Boosted Object NotFound\n");  continue;}
		llvvTauCollection boostedtaus = *boostedtauCollHandle;
                //merged the two tau collections, start by the boosted taus
                for(unsigned int i=0;i<taus.size();i++){boostedtaus.push_back(taus[i]);}
                taus = boostedtaus;

		fwlite::Handle< llvvJetCollection > jetCollHandle;
		jetCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
		if(!jetCollHandle.isValid()){printf("llvvJetCollection Object NotFound\n");  continue;}
		llvvJetExtCollection jets;
		for(unsigned int i=0;i<jetCollHandle->size();i++){jets.push_back(llvvJetExt((*jetCollHandle)[i]));}

		fwlite::Handle< llvvMet > metHandle;
		metHandle.getByLabel(ev, "llvvObjectProducersUsed", "pfMETPFlow"); 
		if(!metHandle.isValid()){printf("llvvMet Object NotFound\n");  continue;}
		llvvMet met = *metHandle;

		fwlite::Handle< std::vector<bool> > triggerBitsHandle;
		triggerBitsHandle.getByLabel(ev, "llvvObjectProducersUsed", "triggerBits");
		if(!triggerBitsHandle.isValid()){printf("triggerBits Object NotFound\n");  continue;}
		std::vector<bool> triggerBits = *triggerBitsHandle;

		fwlite::Handle< std::vector<int> > triggerPrescalesHandle;
		triggerPrescalesHandle.getByLabel(ev, "llvvObjectProducersUsed", "triggerPrescales");
		if(!triggerPrescalesHandle.isValid()){printf("triggerPrescales Object NotFound\n");  continue;}
		std::vector<int> triggerPrescales = *triggerPrescalesHandle;

		fwlite::Handle< double > rhoHandle;
		rhoHandle.getByLabel(ev, "kt6PFJets", "rho");
		if(!rhoHandle.isValid()){printf("rho Object NotFound\n");  continue;}
		double rho = *rhoHandle;

		fwlite::Handle< double > rho25Handle;
		rho25Handle.getByLabel(ev, "kt6PFJetsCentral", "rho");
		if(!rho25Handle.isValid()){printf("rho25 Object NotFound\n");  continue;}
		double rho25 = *rho25Handle;


		//require compatibilitiy of the event with the PD
		bool eeTrigger          = triggerBits[0];
		bool muTrigger          = triggerBits[6];
		bool mumuTrigger        = triggerBits[2] || triggerBits[3] || muTrigger;
		if(filterOnlyEE)   { mumuTrigger=false; }
		if(filterOnlyMUMU) { eeTrigger=false;   }
		if(isSingleMuPD)   { eeTrigger=false; if( mumuTrigger || !muTrigger ) mumuTrigger= false;  }
*/


		//pileup weight
		float weight = 1.0;
		//double TotalWeight_plus = 1.0;
		//double TotalWeight_minus = 1.0;
		float puWeight(1.0);

		if(isMC){
			puWeight          = LumiWeights->weight(ngenITpu) * PUNorm[0];
			weight            = xsecWeight*puWeight;
			//TotalWeight_plus  = PuShifters[utils::cmssw::PUUP  ]->Eval(genEv.ngenITpu) * (PUNorm[2]/PUNorm[0]);
			//TotalWeight_minus = PuShifters[utils::cmssw::PUDOWN]->Eval(genEv.ngenITpu) * (PUNorm[1]/PUNorm[0]);
		}


                std::vector<TString> chTags;
                chTags.push_back("all");

   	        mon.fillHisto("nvtx"        ,   chTags, nvtx,      weight);
		mon.fillHisto("nvtxraw"     ,   chTags, nvtx,      weight/puWeight);


/*

		//
		// BELOW FOLLOWS THE ANALYSIS OF THE MAIN SELECTION WITH N-1 PLOTS
		//

		//
		// LEPTON ANALYSIS
		//
		llvvTauLeptonCollection selLeptons;
		if(examineThisEvent) cout << "***LEPTON ANALYSIS STARTS***" << endl;
		if(examineThisEvent) cout << "Lepton size is: " << leptons.size() << endl;
		for(size_t ilep=0; ilep<leptons.size(); ilep++){
			bool passKin(true),passId(true),passIso(true);
			int lid=leptons[ilep].id;

			//apply muon corrections
			if(examineThisEvent) cout << "Lepton ID is: " << lid << endl;

			if(lid==13 && muCor){
				TLorentzVector p4(leptons[ilep].px(), leptons[ilep].py(), leptons[ilep].pz(), leptons[ilep].energy());
				muCor->applyPtCorrection(p4 , lid<0 ? -1 :1 );
				if(isMC) muCor->applyPtSmearing(p4, lid<0 ? -1 : 1, false);
				leptons[ilep].SetPxPyPzE(p4.Px(), p4.Py(), p4.Pz(), p4.Energy());
			}

			//no need for charge info any longer
			lid=abs(lid);
			TString lepStr( lid==13 ? "mu" : "e");

			//kinematics
			float leta = lid==11 ? leptons[ilep].electronInfoRef->sceta : leptons[ilep].eta();
			if(examineThisEvent) cout << "Lepton pt/eta is: " << leptons[ilep].pt() << "/" << leta << endl;
			if(leptons[ilep].pt()<10)                   passKin=false;
			if(leta> (lid==11 ? 2.5 : 2.4) )            passKin=false;
			if(lid==11 && (leta>1.4442 && leta<1.5660)) passKin=false;

			//id
			Int_t idbits = leptons[ilep].idbits;
			if(lid==11){
				if(leptons[ilep].electronInfoRef->isConv)              passId=false;
				bool isLoose = leptons[ilep].electronInfoRef->mvanontrigv0; 
				if(examineThisEvent) cout << "Lepton ID loose: " << isLoose << endl;
				if(!isLoose)                                   passId=false;
				else                                           passId=true;
			}
			else{
				bool isLoose    = ((idbits >> 8) & 0x1);
				if(!isLoose)                                   passId=false;
			}

			//isolation
			float relIso = utils::cmssw::relIso(leptons[ilep], rho);
			if(examineThisEvent) cout << "Lepton relIso: " << relIso << endl;
			if(examineThisEvent) cout << "passId/passIso/passKin: " << passId << "/" << passIso << "/" << passKin << endl;
			if( (lid==11 && relIso>0.40) || (lid!=11 && relIso>0.40) ) passIso=false;

			if(!passId || !passIso || !passKin) continue;
			selLeptons.push_back(llvvTauLepton(leptons[ilep]));
		}
		std::sort(selLeptons.begin(), selLeptons.end(), sort_llvvObjectByPt);

		//at this point check if it's worth continuing
		if(examineThisEvent) cout << "***LEPTON ANALYSIS IS FINISHED***" << endl;
		if(examineThisEvent) cout << "Lepton size is: " << selLeptons.size() << endl;

		//
		// DILEPTON ANALYSIS
		//
		LorentzVector leadingLep, trailerLep, zll, zlltmp;
		int dilLep1=-1, dilLep2=-1, dilId=-1;
		double BestMass=0;
		bool passZmass=false;
		if(examineThisEvent) cout << "***DILEPTON ANALYSIS STARTS***" << endl;
		//identify the best lepton pair
		for(unsigned int l1=0   ;l1<selLeptons.size();l1++){
			float relIso1 = utils::cmssw::relIso(selLeptons[l1].lep, rho);
			if(examineThisEvent) cout << "ID/Pt/ISO lep1: " << selLeptons[l1].id <<"/"<<selLeptons[l1].pt()<<"/"<<relIso1<< endl;
			if( relIso1>0.30 ) continue;
			for(unsigned int l2=l1+1;l2<selLeptons.size();l2++){
				float relIso2 = utils::cmssw::relIso(selLeptons[l2].lep, rho);
				if(examineThisEvent) cout << "ID/Pt/ISO lep2: " << selLeptons[l2].id <<"/"<<selLeptons[l2].pt()<<"/"<<relIso2<< endl;
				if(abs(selLeptons[l1].id)!=abs(selLeptons[l2].id)) continue; //only consider same flavor lepton pairs
				if(selLeptons[l1].id*selLeptons[l2].id>=0) continue; //only consider opposite charge lepton pairs
				if( !((selLeptons[l1].pt()>=20 && selLeptons[l2].pt()>=10) || (selLeptons[l1].pt()>=10 && selLeptons[l2].pt()>=20))) continue;
				if( relIso2>0.30 ) continue;
				if(deltaR(selLeptons[l1], selLeptons[l2])<0.1) continue;
				zlltmp = (selLeptons[l1]+selLeptons[l2]);
				if(examineThisEvent) cout << "reco Z mass: " << zlltmp.mass() << endl;
				if(fabs(zlltmp.mass() - 91.2) < fabs(BestMass-91.2) && zlltmp.mass()>60 && zlltmp.mass()<120){
					dilLep1 = l1; 
					dilLep2 = l2;
					zll=selLeptons[l1]+selLeptons[l2];
					//mon.fillHisto("isomu"      ,   "check",relIso1, weight);
					//mon.fillHisto("isoele"      ,  "check",relIso2, weight);
					leadingLep=selLeptons[l1];
					trailerLep=selLeptons[l2];
					dilId = selLeptons[l1].id * selLeptons[l2].id;
					BestMass=zll.mass();
					passZmass=true;
					if(examineThisEvent) cout << "final Z mass " << zll.mass()  << " leading/trailer pt " << leadingLep.pt() << "/" << trailerLep.pt() << endl;
				}
			}
		}
		if(examineThisEvent) cout << "***DILEPTON ANALYSIS IS FINISHED***" << endl;

                // **********************************************************************************************************************************
                //THIS MUST BE UNCOMMENTED FOR THE NORMAL ANALYSIS !!!
		//apply data/mc correction factors
		if(dilLep1>=0)weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[dilLep1].pt(), selLeptons[dilLep1].eta(), abs(selLeptons[dilLep1].id),  abs(selLeptons[dilLep1].id) ==11 ? "loose" : "loose" ).first : 1.0;
		if(dilLep2>=0)weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[dilLep2].pt(), selLeptons[dilLep2].eta(), abs(selLeptons[dilLep2].id),  abs(selLeptons[dilLep2].id) ==11 ? "loose" : "loose" ).first : 1.0;
                // **********************************************************************************************************************************

		std::vector<TString> chTags;
		bool isDileptonCandidate = false;
		chTags.push_back("all");
		if( abs(dilId)==121 && eeTrigger  ){ chTags.push_back("ee"); isDileptonCandidate=true; }
		if( abs(dilId)==169 && mumuTrigger){ chTags.push_back("mumu"); isDileptonCandidate=true; }
		if( !isDileptonCandidate           ) chTags.push_back("ct");
		
		bool passZpt = (zll.pt()>50);
		bool passZeta = true;//(fabs(zll.eta())<1.4442);

		//
		//JET/MET ANALYSIS
		//
		llvvJetExtCollection selJets, selJetsNoId, selBJets;
		int njets(0), nbjets(0);
		for(size_t ijet=0; ijet<jets.size(); ijet++){
			//correct jet
			float toRawSF=jets[ijet].torawsf;
			LorentzVector rawJet(jets[ijet]*toRawSF);
			jesCor->setJetEta(rawJet.eta());
			jesCor->setJetPt(rawJet.pt());
			jesCor->setJetA(jets[ijet].area);
			jesCor->setRho(rho);
			float newJECSF=jesCor->getCorrection();
			jets[ijet].SetPxPyPzE(rawJet.px(),rawJet.py(),rawJet.pz(),rawJet.energy());
			jets[ijet] *= newJECSF;
			jets[ijet].torawsf = 1./newJECSF;
			if(jets[ijet].pt()<15 || fabs(jets[ijet].eta())>4.7 ) continue;

			if(examineThisEvent) cout << "Good jet with pt>15 and eta<4.7" << endl;
			//bjets
			mon.fillHisto("bjetpt"    ,  chTags, jets[ijet].pt(),  weight);
			mon.fillHisto("bjetcsv"   ,  chTags, jets[ijet].origcsv,  weight);
			if(jets[ijet].pt()>20 && fabs(jets[ijet].eta())<2.4 && jets[ijet].origcsv>0.679){
				selBJets.push_back(jets[ijet]);  
				nbjets++;
			}

			//cross-clean with selected leptons and photons
			double minDRlj(9999.),minDRlg(9999.);
			for(size_t ilep=0; ilep<selLeptons.size(); ilep++)
				minDRlj = TMath::Min( (double) minDRlj, (double) deltaR(jets[ijet],selLeptons[ilep]) );
			if(examineThisEvent) cout << "dR jets-leptons (should be > 0.4)" << minDRlj << endl;
			if(minDRlj<0.4 || minDRlg<0.4) 
				continue;

			//jet id
			// float pumva=jets[ijet].puMVA;
			Int_t idbits=jets[ijet].idbits;
			bool passPFloose( ((idbits>>0) & 0x1));
			int puId( ( idbits >>3 ) & 0xf );
			bool passLoosePuId( ( puId >> 2) & 0x1);
			int simplePuId( ( idbits >>7 ) & 0xf );
			bool passLooseSimplePuId(  ( simplePuId >> 2) & 0x1);
			TString jetType( jets[ijet].genj.pt()>0 ? "truejetsid" : "pujetsid" );
			if(jets[ijet].pt()>30)
			{
				mon.fillHisto(jetType,chTags,fabs(jets[ijet].eta()),0);
				if(passPFloose)                        mon.fillHisto(jetType,chTags,fabs(jets[ijet].eta()),1);
				if(passPFloose && passLoosePuId)       mon.fillHisto(jetType,chTags,fabs(jets[ijet].eta()),2);
				if(passPFloose && passLooseSimplePuId) mon.fillHisto(jetType,chTags,fabs(jets[ijet].eta()),3);
				if(passLoosePuId)                      mon.fillHisto(jetType,chTags,fabs(jets[ijet].eta()),4);
				if(passLooseSimplePuId)                mon.fillHisto(jetType,chTags,fabs(jets[ijet].eta()),5);
			}

			//add scale/resolution uncertainties
			std::vector<float> smearPt=utils::cmssw::smearJER(jets[ijet].pt(),jets[ijet].eta(),jets[ijet].genj.pt());
			jets[ijet].jer     = isMC ? smearPt[0] : jets[ijet].pt();
			jets[ijet].jerup   = isMC ? smearPt[1] : jets[ijet].pt();
			jets[ijet].jerdown = isMC ? smearPt[2] : jets[ijet].pt();
			smearPt=utils::cmssw::smearJES(jets[ijet].pt(),jets[ijet].eta(), totalJESUnc);
			jets[ijet].jesup   = isMC ? smearPt[0] : jets[ijet].pt();
			jets[ijet].jesdown = isMC ? smearPt[1] : jets[ijet].pt();

			selJetsNoId.push_back(jets[ijet]);
			if(passPFloose && passLooseSimplePuId){
				selJets.push_back(jets[ijet]);
				if(jets[ijet].pt()>minJetPtToApply) njets++;
			}
		}
		std::sort(selJets.begin(), selJets.end(), sort_llvvObjectByPt);
		std::sort(selBJets.begin(), selBJets.end(), sort_llvvObjectByPt);

		if(examineThisEvent) cout << "number of bjets" << nbjets << endl;

		if(examineThisEvent) cout << "***TAU ANALYSIS STARTS***" << endl;
	
		//
		// TAU ANALYSIS
		//
		if(examineThisEvent) cout << "tau size " << taus.size() << endl;
                llvvTauCollection selTaus;
		for(size_t itau=0; itau<taus.size(); itau++){
			llvvTau& tau = taus[itau];
			if(examineThisEvent) cout << "tau n: " << itau << " pt/eta" << tau.pt() << "/" << fabs(tau.eta()) << endl;
			if(examineThisEvent) cout << "muonLoose/DM " << tau.passId(llvvTAUID::againstMuonLoose2)  << "/" << tau.passId(llvvTAUID::decayModeFinding) << endl;
			mon.fillHisto("taupt"           ,   chTags, tau.pt(), weight);  //check the tau pt. YO!
			if(tau.pt()<15.0 || fabs(tau.eta()) > 2.3) continue; 

			bool overalWithLepton=false;
			for(int l1=0   ;l1<(int)selLeptons.size();l1++){
				//mon.fillHisto("drtaulep"           ,   chTags, deltaR(tau, selLeptons[l1]), weight); //check the dR between taus and leptons. YO!
				if(examineThisEvent) cout << "dR(tau,lep) " << deltaR(tau, selLeptons[l1]) << endl;
				if(deltaR(tau, selLeptons[l1])<0.1){overalWithLepton=true; break;}
			}
			if(overalWithLepton) continue;
			if(examineThisEvent) cout << "No overlap tau-leptons...GOOD! " << overalWithLepton << endl;

			//mon.fillHisto("antimuon2","checkTau", tau.passId(llvvTAUID::againstMuonLoose2),weight);
			//mon.fillHisto("decaymode","checkTau", tau.passId(llvvTAUID::decayModeFinding), weight);
			if(!tau.passId(llvvTAUID::againstMuonLoose2)) continue; 
			if(!tau.passId(llvvTAUID::decayModeFinding)) continue;
			//if(!tau.passId(llvvTAUID::decayModeFindingNewDMs)) continue;

                        selTaus   .push_back(tau);
                        selLeptons.push_back(llvvTauLepton(tau)); //Dirty Trick to add selected taus to the selected leptons vector --> A cast to llvvTauLepton will be needed to get all the info
		}
                std::sort(selLeptons.begin(), selLeptons.end(), sort_llvvObjectByPt); //resort the vector since taus have been added at the back
                //Since we resort the vector we screwed up the index for dilLep1 and dilLep2, we need to reasign these index!
                if(dilLep1>=0){double mindR=1E100;  for(unsigned int l=0;l<selLeptons.size();l++){double dR = deltaR(leadingLep, selLeptons[l]); if(dR<mindR){mindR=dR; dilLep1=l;} }  }
                if(dilLep2>=0){double mindR=1E100;  for(unsigned int l=0;l<selLeptons.size();l++){double dR = deltaR(trailerLep, selLeptons[l]); if(dR<mindR){mindR=dR; dilLep2=l;} }  }

		//
		// HIGGS ANALYSIS
		//
		if(examineThisEvent) cout << "THE HIGGS ANALYSIS STARTS" << endl;

                float weightBeforeLepCorr = weight; //save the weight before the lepton corrections.
                                                    //necessary for the optimization.

		std::vector<TString> chTagsMain=chTags;


                //FIND the two highest pT leptons not coming from the Z, and with dR>0.1 from all other leptons in the event
                int higgsCandL1=-1, higgsCandL2=-1;
                for(int l=0   ;l<(int)selLeptons.size();l++){
                   if(l==dilLep1 || l==dilLep2)continue;
//                   printf("%+2i - pT=%+6.2f eta=%+6.2f phi=%+6.2f\n", selLeptons[l].id, selLeptons[l].pt(), selLeptons[l].eta(), selLeptons[l].phi());
                   if(deltaR(selLeptons[l],  selLeptons[dilLep1])<0.1)continue;
                   if(deltaR(selLeptons[l],  selLeptons[dilLep2])<0.1)continue;
                   if(higgsCandL1<0){higgsCandL1=l;continue;}
                   if(higgsCandL2<0 && deltaR(selLeptons[l],  selLeptons[higgsCandL1])>=0.1){higgsCandL2=l;break;}//ordered in pT, so all done
                }

              //Build the higgs candidate and determine the higgs category
 	      LorentzVector higgsCand(0,0,0,0);
              int HiggsShortId=-1;  int higgsCandId=0;
              string ChannelName = "none";   string signName = "";
              if(higgsCandL1>=0 && higgsCandL2>=0){
                   higgsCandId=selLeptons[higgsCandL1].id*selLeptons[higgsCandL2].id;
                   higgsCand = LorentzVector(selLeptons[higgsCandL1]+selLeptons[higgsCandL2]);
                   if(higgsCandId<0){signName="_OS";}else{signName="_SS";}
                   if(higgsCandId<0){HiggsShortId = 0;}else{HiggsShortId = 8;}
                   if(abs(selLeptons[dilLep1].id)==13){HiggsShortId += 0;}else{HiggsShortId += 4;}
                   switch(abs(higgsCandId)){
                      case 11*13:  ChannelName  = "elmu";  HiggsShortId+= 0; break;
                      case 11*15:  ChannelName  = "elha";  HiggsShortId+= 1; break;
                      case 13*15:  ChannelName  = "muha";  HiggsShortId+= 2; break;
                      case 15*15:  ChannelName  = "haha";  HiggsShortId+= 3; break;
                      default:     ChannelName  = "none";  HiggsShortId =-1; break;
                   }
              }               
   	      chTagsMain.push_back(chTagsMain[chTagsMain.size()-1] + signName + ChannelName); 

              //reweight the event to account for lept eff.
              if(isMC && higgsCandL1>=0 && abs(selLeptons[higgsCandL1].id)<15)weight *= lepEff.getLeptonEfficiency( selLeptons[higgsCandL1].pt(), selLeptons[higgsCandL1].eta(), abs(selLeptons[higgsCandL1].id), abs(selLeptons[higgsCandL1].id) ==11 ? "loose" : "loose" ).first;
              if(isMC && higgsCandL2>=0 && abs(selLeptons[higgsCandL2].id)<15)weight *= lepEff.getLeptonEfficiency( selLeptons[higgsCandL2].pt(), selLeptons[higgsCandL2].eta(), abs(selLeptons[higgsCandL2].id), abs(selLeptons[higgsCandL2].id) ==11 ? "loose" : "loose" ).first;

              //check if the pair pass Lepton Veto
              bool passLepVetoMain = true;
              for(int l=0;l<(int)selLeptons.size() && passLepVetoMain;l++){
                 if(l==dilLep1 || l==dilLep2 || l==higgsCandL1 || l==higgsCandL2) continue; //lepton already used in the dilepton pair or higgs candidate
                 if(abs(selLeptons[l].id)==15){
                    llvvTau    tau = selLeptons[l].tau;
                    if(tau.pt()<20) continue;
                    if(!tau.passId(llvvTAUID::againstElectronLoose) ||
                       !tau.passId(llvvTAUID::againstMuonLoose2)    ||
                       !tau.passId(llvvTAUID::byMediumCombinedIsolationDeltaBetaCorr3Hits)  ) continue;                    
                    passLepVetoMain = false; break;
                 }else{
                    passLepVetoMain = false; break;
                 }
              }

              //check if the pair pass b-jet veto
              bool passBJetVetoMain = true;
              for(int j1=0;j1<(int)selBJets.size();j1++){
                 if(dilLep1    !=-1 && deltaR(selBJets[j1]   , selLeptons[dilLep1 ])>0.4){passBJetVetoMain=false; break;}
                 if(dilLep2    !=-1 && deltaR(selBJets[j1]   , selLeptons[dilLep2 ])>0.4){passBJetVetoMain=false; break;}
                 if(higgsCandL1         !=-1 && deltaR(selBJets[j1]   , selLeptons[higgsCandL1      ])>0.4){passBJetVetoMain=false; break;}
                 if(higgsCandL2         !=-1 && deltaR(selBJets[j1]   , selLeptons[higgsCandL2      ])>0.4){passBJetVetoMain=false; break;}
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

              //check if the event pass loosest possible higgs selection
              bool passDPhiCut    =  (fabs(deltaPhi(zll.phi(), met.phi()))>1.5);
              bool passHiggsLoose = passHiggsCuts(selLeptons, rho, higgsCandId, higgsCandL1, higgsCandL2, 0.4, 0.4, llvvTAUID::decayModeFinding, 20); 
              bool passHiggsMain  = passHiggsLoose && passHiggsCuts(selLeptons, rho, higgsCandId, higgsCandL1, higgsCandL2, 0.2, 0.2, llvvTAUID::byMediumCombinedIsolationDeltaBetaCorr3Hits, 25);


		//SVFIT MASS
		LorentzVector higgsCand_SVFit = higgsCand;
		if(passZpt && passZeta && passDPhiCut && passHiggsLoose && passLepVetoMain && passBJetVetoMain){
//                  higgsCand_SVFit = getSVFit(met, selLeptons, higgsCandL1, higgsCandL2);  //compute svfit mass in a smart way
		}

                //build the higgs candH
		LorentzVector higgsCandH       = zll + higgsCand;
                LorentzVector higgsCandH_SVFit = zll + higgsCand_SVFit;
*/                


		//
		// NOW FOR THE CONTROL PLOTS
		//
/*
		mon.fillHisto("eventflow"      ,   chTagsMain,                 0, weight);
		if(selLeptons.size()>=2){
			//mon.fillHisto("nlep"           ,   chTags, selLeptons.size(), weight);
			mon.fillHisto("eventflow"      ,   chTagsMain,                 1, weight);
			if(passZmass){
				mon.fillHisto("eventflow"   ,   chTagsMain,                 2, weight);

				//pu control
				mon.fillHisto("nvtx"        ,   chTagsMain, nvtx,      weight);
				mon.fillHisto("nvtxraw"     ,   chTagsMain, nvtx,      weight/puWeight);
				mon.fillHisto("rho"         ,   chTagsMain, rho,       weight);
				mon.fillHisto("rho25"       ,   chTagsMain, rho25,     weight);

				//Z kinematics control
				mon.fillHisto("leadpt"      ,   chTagsMain, leadingLep.pt(), weight);      
				mon.fillHisto("trailerpt"   ,   chTagsMain, trailerLep.pt(), weight);      
				mon.fillHisto("leadeta"     ,   chTagsMain, leadingLep.eta(), weight);      
				mon.fillHisto("trailereta"  ,   chTagsMain, trailerLep.eta(), weight);      

				//analyze dilepton kinematics
				mon.fillHisto("zpt"         ,   chTagsMain, zll.pt(),      weight);      
				mon.fillHisto("zmass"       ,   chTagsMain, zll.mass(),    weight);  
				mon.fillHisto("zeta"        ,   chTagsMain, zll.eta(),     weight);
				mon.fillHisto("zy"          ,   chTagsMain, zll.Rapidity(),weight);

				if(passZpt && passZeta){
					mon.fillHisto("eventflow",   chTagsMain,                 3, weight);

					mon.fillHisto("ntaus"        ,  chTags, selTaus.size(), weight);
					mon.fillHisto("tauleadpt"    ,  chTagsMain, selTaus.size()>0?selTaus[0].pt():-1,  weight);
					mon.fillHisto("tauleadeta"   ,  chTagsMain, selTaus.size()>0?selTaus[0].eta():-10, weight);

					if(selLeptons.size()>=4){
						mon.fillHisto("eventflow",   chTagsMain,                 4, weight);

							if(passLepVetoMain){
								mon.fillHisto("eventflow",   chTagsMain,                 5, weight);

								if(passBJetVetoMain){
									mon.fillHisto("eventflow"	,   chTagsMain,                 6, weight);

                                                                        mon.fillHisto("dPhi_AZ"         , chTagsMain, fabs(deltaPhi(higgsCand.phi(), zll.phi())),    weight);
                                                                        mon.fillHisto("dPhi_AMet"       , chTagsMain, fabs(deltaPhi(higgsCand.phi(), met.phi())),    weight);
                                                                        mon.fillHisto("dPhi_ZMet"       , chTagsMain, fabs(deltaPhi(zll.phi(), met.phi())),    weight);
                                					mon.fillHisto("met"      	, chTagsMain, met.pt()         , weight);

				                                	if(passDPhiCut){
                          			       				mon.fillHisto("eventflow",   chTagsMain,                 7, weight);


			                        				if(passHiggsMain){
											mon.fillHisto("eventflow",   chTagsMain,                 8, weight);
											mon.fillHisto("eventflow"	,   chTagsMain,                10+HiggsShortId, weight);


											mon.fillHisto("Apt"       	, chTagsMain, higgsCand.pt(),    weight);
											mon.fillHisto("Amass"           , chTagsMain, higgsCand.mass(),  weight);
											mon.fillHisto("Amasssvfit"      , chTagsMain, higgsCand_SVFit.mass(),  weight);
											mon.fillHisto("Hmass"           , chTagsMain, higgsCandH.mass(),  weight);
											if(higgsCand.mass()>0 && higgsCand.mass()<=25)
												mon.fillHisto("Hmass_MA_0_25"           , chTagsMain, higgsCandH.mass(),  weight);
											if(higgsCand.mass()>25 && higgsCand.mass()<=50)
												mon.fillHisto("Hmass_MA_25_50"          , chTagsMain, higgsCandH.mass(),  weight);
											if(higgsCand.mass()>50 && higgsCand.mass()<=100)
												mon.fillHisto("Hmass_MA_50_100"         , chTagsMain, higgsCandH.mass(),  weight);
											if(higgsCand.mass()>100 && higgsCand.mass()<=200)
												mon.fillHisto("Hmass_MA_100_200"        , chTagsMain, higgsCandH.mass(),  weight);
											mon.fillHisto("Hpt"             , chTagsMain, higgsCandH.pt(),  weight);
											mon.fillHisto("Hmasssvfit"      , chTagsMain, higgsCandH_SVFit.mass(),  weight);
											if(higgsCand_SVFit.mass()>0 && higgsCand_SVFit.mass()<=25)
												mon.fillHisto("Hmasssvfit_MA_0_25"      , chTagsMain, higgsCandH_SVFit.mass(),  weight);
											if(higgsCand_SVFit.mass()>25 && higgsCand_SVFit.mass()<=50)
												mon.fillHisto("Hmasssvfit_MA_25_50"             , chTagsMain, higgsCandH_SVFit.mass(),  weight);
											if(higgsCand_SVFit.mass()>50 && higgsCand_SVFit.mass()<=100)
												mon.fillHisto("Hmasssvfit_MA_50_100"            , chTagsMain, higgsCandH_SVFit.mass(),  weight);
											if(higgsCand_SVFit.mass()>100 && higgsCand_SVFit.mass()<=200)
												mon.fillHisto("Hmasssvfit_MA_100_200"           , chTagsMain, higgsCandH_SVFit.mass(),  weight);
											mon.fillHisto("vismass2D"       ,  chTagsMain, higgsCand.mass(), higgsCandH.mass(), weight);
											mon.fillHisto("svfitmass2D"     ,  chTagsMain, higgsCand_SVFit.mass(), higgsCandH_SVFit.mass(), weight);

											mon.fillHisto("Anjets"    	, chTagsMain, NCleanedJetMain      , weight); 
											mon.fillHisto("Amet"      	, chTagsMain, met.pt()         , weight);
											fprintf(outTxtEvents, "%d %d %d\n",ev.eventAuxiliary().luminosityBlock(),ev.eventAuxiliary().run(),ev.eventAuxiliary().event());
								}//else{mon.fillHisto("failreason",chTags,11,weight);      } //BJETVETO 
							}//else{mon.fillHisto("failreason",chTags,10,weight);      } //LEPVETO
						}//else{mon.fillHisto("failreason",chTags,9,weight);      } //HIGGS
					}//else{mon.fillHisto("failreason",chTags,8,weight);      } //dPhiCut
                                   }//NLeptons+NTaus>=4
				}//else{mon.fillHisto("failreason",chTags,7,weight);      } //ZKin
			}//else{mon.fillHisto("failreason",chTags,6,weight);      } //ZMass
		}//else{mon.fillHisto("failreason",chTags,5,weight);      } //NLEP  
*/

/*
		//SYSTEMATIC STUDY on all events passing the basic preselection
                if(passZmass && passZpt && passZeta && selLeptons.size()>=4 && passLepVetoMain && passBJetVetoMain && passDPhiCut && passHiggsLoose){
			for(size_t ivar=0; ivar<nvarsToInclude; ivar++){
                        //independent on the index

                       std::vector<TString>& locTags =  chTagsMain;                           


                       int closestJetIndexL1=-1;  double pTL1=-1;  double dRL1=-1;
                       double dRminL1 = closestJet(selLeptons[higgsCandL1], selJets, closestJetIndexL1);
		       if(closestJetIndexL1>=0 && dRminL1<0.5){pTL1=selJets[closestJetIndexL1].pt(); dRL1=dRminL1;}else{pTL1=selLeptons[higgsCandL1].pt(); dRL1=-1;}
                       int closestJetIndexL2=-1;  double pTL2=-1;  double dRL2=-1;
                       double dRminL2 = closestJet(selLeptons[higgsCandL2], selJets, closestJetIndexL2);
		       if(closestJetIndexL2>=0 && dRminL2<0.5){pTL2=selJets[closestJetIndexL2].pt(); dRL2=dRminL2;}else{pTL2=selLeptons[higgsCandL2].pt(); dRL2=-1;}

                       treeEventId  = ev.eventAuxiliary().event();
                       treeLumiId   = ev.eventAuxiliary().luminosityBlock();
                       treeRunId    = ev.eventAuxiliary().run();
                       treeHiggsId  = higgsCandId;
                       treeVisMass  = higgsCandH.mass();
                       treeSVFMass  = higgsCandH_SVFit.mass();                                
                       treeLeg1Id   = selLeptons[higgsCandL1].id;
                       treeLeg1DR   = dRL1;
                       treeLeg1Pt   = pTL1;
                       treeLeg1Eta  = selLeptons[higgsCandL1].eta();
                       treeLeg1Iso  = abs(treeLeg1Id)==15?(float)(selLeptons[higgsCandL1].tau.idbits&0xFFFFFFFF):utils::cmssw::relIso(selLeptons[higgsCandL1].lep, rho);
                       treeLeg2Id   = selLeptons[higgsCandL2].id;
                       treeLeg2DR   = dRL2;
                       treeLeg2Pt   = pTL2;
                       treeLeg2Eta  = selLeptons[higgsCandL2].eta();
                       treeLeg2Iso  = abs(treeLeg2Id)==15?(float)(selLeptons[higgsCandL2].tau.idbits&0xFFFFFFFF):utils::cmssw::relIso(selLeptons[higgsCandL2].lep, rho);

		       for(unsigned int index=0; index<optim_Cuts_sumPt.size();index++){
		          float iweight = weightBeforeLepCorr; //nominal
			if(ivar==5)                        iweight *= TotalWeight_plus;        //pu up
			if(ivar==6)                        iweight *= TotalWeight_minus;       //pu down
			//                  if(ivar==7)                        iweight *= Q2Weight_plus;
			//                  if(ivar==8)                        iweight *= Q2Weight_down;
			//                  if(ivar==9)                        iweight *= PDFWeight_plus;
			//                  if(ivar==10)                       iweight *= PDFWeight_down;



                          bool passHiggs = passHiggsCuts(selLeptons, rho, higgsCandId, higgsCandL1, higgsCandL2, optim_Cuts_elIso[index], optim_Cuts_muIso[index], (float)optim_Cuts_taIso[index], optim_Cuts_sumPt[index]);
		          if(passHiggs){
                                treeCutIndex = index;
				mon.fillHisto(TString("svfit_shapes")+varNames[ivar],locTags,index,higgsCand_SVFit.mass(),iweight);
				mon.fillHisto(TString("Hsvfit_shapes")+varNames[ivar],locTags,index,higgsCandH_SVFit.mass(),iweight);
				mon.fillHisto(TString("FR_closestJetPt")+varNames[ivar],locTags,index,pTL1,iweight);
                                mon.fillHisto(TString("FR_closestJetPt")+varNames[ivar],locTags,index,pTL2,iweight);
                                if(ivar==0 && tree)tree->Fill();
 			 }		    
		       }//end of the loop on cutIndex
		}//end of the loop on the systematics
           }//end the IF condition
*/


	}//end of the event loop  
	printf("\n"); 

	//##############################################
	//########     SAVING HISTO TO FILE     ########
	//##############################################
	//save control plots to file
	outUrl += "/";
	outUrl += outFileUrl + ".root";
	printf("Results save in %s\n", outUrl.Data());

	//save all to the file
	TFile *ofile=TFile::Open(outUrl, "recreate");
	mon.Write();
        //if(tree){tree->SetDirectory(ofile); tree->Write();}
	ofile->Close();

	//save summary tuple
//	outUrl.ReplaceAll(".root","_summary.root");
//	ofile=TFile::Open(outUrl,"recreate");
//	tree->SetDirectory(ofile);
//	tree->Write();
//	ofile->Close();
	if(outTxtFile)fclose(outTxtFile);
//	fclose(outTxtEvents);

}  




