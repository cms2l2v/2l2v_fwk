#include <iostream>
#include <boost/shared_ptr.hpp>

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h" //for svfit

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


class SVFitBooster{
   public:
      std::map<uint32_t, double> precomputed;      
      SVFitBooster(){}
      ~SVFitBooster(){}            
   public:
      double getSVFit(llvvMet& met, llvvLeptonCollection& selLeptons, llvvTauCollection& selTaus, int higgsCandMu, int higgsCandEl, int higgsCandT1, int higgsCandT2){       
         //build key lookup index:
         if(higgsCandMu<0)higgsCandMu=255;           if(higgsCandEl<0)higgsCandEl=255;         if(higgsCandT1<0)higgsCandT1=255;         if(higgsCandT2<0)higgsCandT2=255;  //to make sure that Id=0 is OK (set to highest possible 8bit values)
         uint32_t key = ((higgsCandMu&0x0F)<<24) + ((higgsCandEl&0x0F)<<16) + ((higgsCandT1&0x0F)<<8) + ((higgsCandT2&0x0F)<<0);

         if(precomputed.find(key)==precomputed.end()){ //svfit mass not yet computed for this set of lepton, compute it...
            TMatrixD covMET(2, 2); // PFMET significance matrix
            covMET[0][0] = met.sigx2;
            covMET[0][1] = met.sigxy;
            covMET[1][0] = met.sigxy;
            covMET[1][1] = met.sigy2;

            std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
            if(higgsCandMu!=255)measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kLepDecay, svFitStandalone::LorentzVector(selLeptons[higgsCandMu].px(), selLeptons[higgsCandMu].py(), selLeptons[higgsCandMu].pz(), selLeptons[higgsCandMu].E()) ));
            if(higgsCandEl!=255)measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kLepDecay, svFitStandalone::LorentzVector(selLeptons[higgsCandEl].px(), selLeptons[higgsCandEl].py(), selLeptons[higgsCandEl].pz(), selLeptons[higgsCandEl].E()) ));
            if(higgsCandT1!=255)measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kHadDecay, svFitStandalone::LorentzVector(selTaus   [higgsCandT1].px(), selTaus   [higgsCandT1].py(), selTaus   [higgsCandT1].pz(), selTaus   [higgsCandT1].E()) ));
            if(higgsCandT2!=255)measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kHadDecay, svFitStandalone::LorentzVector(selTaus   [higgsCandT2].px(), selTaus   [higgsCandT2].py(), selTaus   [higgsCandT2].pz(), selTaus   [higgsCandT2].E()) ));
            SVfitStandaloneAlgorithm algo(measuredTauLeptons, svFitStandalone::Vector(met.px(), met.py(), 0) , covMET, 0);
            algo.addLogM(false);
//            algo.integrateMarkovChain();
            //algo.integrateVEGAS(); ////Use this instead for VEGAS integration
            algo.fit();
            if(algo.isValidSolution()){
               //double diTauMassErr = algo.massUncert();
               precomputed[key] = algo.getMass();
             }else{
               precomputed[key] = -1;
             }
         }
         return precomputed[key];
   }
};



       LorentzVector buildCandidates( llvvLeptonCollection& leptons, llvvTauCollection& taus,  
                                      llvvJetExtCollection& jets, llvvJetExtCollection& bjets,
                                      int& higgsCandId, int& higgsCandMu, int& higgsCandEl, int& higgsCandT1, int& higgsCandT2,
                                      bool& passHiggs, int& HiggsShortId, vector<TString>& chTags,
                                      bool& passLepVeto, bool& passBJetVeto, int& NCleanedJet,
                                      vector<float>& isoLep, vector<float>& sumPt, float& charge, 
                                      int dilLep1, int dilLep2, double rho, float weight, bool isMC, LeptonEfficiencySF lepEff){


	       LorentzVector higgsCand(0,0,0,0);
               HiggsShortId=0;
               higgsCandId=0;  higgsCandMu=-1; higgsCandEl=-1; higgsCandT1=-1; higgsCandT2=-1;
               llvvTAUID hpsIsoTau;

               passHiggs    = false;
               passBJetVeto = true;
               passLepVeto  = true;
               NCleanedJet  = 0;
               //**************************************
               // e - mu final state
               //**************************************

               for(int l1=0; l1<(int)leptons.size() && !higgsCandId; l1++){                       
                       float relIso1 = utils::cmssw::relIso(leptons[l1], rho);
                       if( relIso1 > isoLep.at(0) ) continue;
                       for(int l2=l1+1; l2<(int)leptons.size() && !higgsCandId; l2++){
                               float relIso2 = utils::cmssw::relIso(leptons[l2], rho);
                               cout << "lepton 1 id " << leptons[l1].id << " lepton 2 id: " << leptons[l2].id << endl;
                               cout << "iso 1 cut " << isoLep.at(0) << " iso lepton 1: " << relIso1 << endl;
                               cout << "iso 2 cut " << isoLep.at(0) << " iso lepton 2: " << relIso2 << endl;
                               if( relIso2 > isoLep.at(0) ) continue;
                               if(l1==dilLep1 || l1==dilLep2 || l2==dilLep1 || l2==dilLep2) continue; //lepton already used in the dilepton pair
			       if(leptons[l1].id*leptons[l2].id!=(charge*143)) continue;//Only consider opposite OR same sign objects, depending on your needs
                               if((leptons[l1].pt()+leptons[l2].pt()) < sumPt.at(0) ) continue;//Only consider pairs with LT>25 GeV

                               if(deltaR(leptons[l1], leptons[dilLep1])<0.1) continue;
                               if(deltaR(leptons[l1], leptons[dilLep2])<0.1) continue;
                               if(deltaR(leptons[l2], leptons[dilLep1])<0.1) continue;
                               if(deltaR(leptons[l2], leptons[dilLep2])<0.1) continue;
                               if(deltaR(leptons[l1], leptons[l2     ])<0.1) continue;
                               
                               int muId, elId;
                               if(abs(leptons[l1].id)==13){muId=l1; elId=l2;}else{muId=l2; elId=l1;}

                               higgsCandId=leptons[muId].id*leptons[elId].id;  higgsCandMu=muId; higgsCandEl=elId;
                               cout << "FOUND EM CANDIDATE of charge: " << higgsCandId << endl;
                               higgsCand = LorentzVector(leptons[muId]+leptons[elId]);
                               break;
                       }
               }//close e-mu selection
               for(int l1=0; l1<(int)leptons.size()&& !higgsCandId;l1++){
                       if(l1==dilLep1 || l1==dilLep2) continue; //lepton already used in the dilepton pair
                       for(int t1=0; t1<(int)taus   .size()&& !higgsCandId;t1++){
                               if(taus[t1].pt()<15 || fabs(taus[t1].eta())>2.3) continue;

                               if(deltaR(leptons[l1], leptons[dilLep1])<0.1) continue;
                               if(deltaR(leptons[l1], leptons[dilLep2])<0.1) continue;
                               if(deltaR(taus[t1]   , leptons[dilLep1])<0.1) continue;
                               if(deltaR(taus[t1]   , leptons[dilLep2])<0.1) continue;
                               if(deltaR(leptons[l1], taus[t1        ])<0.1) continue;
                               float relIso1 = utils::cmssw::relIso(leptons[l1], rho);
                               cout << "lepton 1 id " << leptons[l1].id << " tau id: " << taus[t1].id << endl;
                               cout << "iso ET cut " << isoLep.at(1) << " iso lepton : " << relIso1 << endl;
                               cout << "iso MT cut " << isoLep.at(2) << " iso lepton : " << relIso1 << endl;
                               cout << "tau iso loose? " << taus[t1].passId(llvvTAUID::byLooseCombinedIsolationDeltaBetaCorr3Hits) << endl;
                               cout << "tau against ele MVA5? " << taus[t1].passId(llvvTAUID::againstElectronTightMVA5) << endl;
                               cout << "tau against ele loose? " << taus[t1].passId(llvvTAUID::againstElectronLoose) << endl;
                               cout << "tau against muon tight2? " << taus[t1].passId(llvvTAUID::againstMuonTight2) << endl;

                               //e-tau
                               if(abs(leptons[l1].id)==11 &&
					       (leptons[l1].id*taus[t1].id!=(charge*165) || 
						relIso1>isoLep.at(1) || !taus[t1].passId(llvvTAUID::againstElectronTightMVA5) ||
						!taus[t1].passId(llvvTAUID::byLooseCombinedIsolationDeltaBetaCorr3Hits) ||
						(taus[t1].pt()+leptons[l1].pt()) < sumPt.at(1))) continue;
                               //mu-tau
                               if(abs(leptons[l1].id)!=11 &&
					       (leptons[l1].id*taus[t1].id!=(charge*195) ||
						relIso1>isoLep.at(2) || !taus[t1].passId(llvvTAUID::againstElectronLoose) ||
						!taus[t1].passId(llvvTAUID::againstMuonTight2) || !taus[t1].passId(llvvTAUID::byLooseCombinedIsolationDeltaBetaCorr3Hits) ||
						(taus[t1].pt()+leptons[l1].pt()) < sumPt.at(2))) continue;

                               higgsCandId=leptons[l1].id*taus[t1].id;  if(abs(leptons[l1].id)==11){higgsCandEl=l1;}else{higgsCandMu=l1;} higgsCandT1=t1;
                               cout << "FOUND LT CANDIDATE of charge " << higgsCandId << endl;
                               higgsCand = LorentzVector(leptons[l1]+taus[t1]);
                               break;
                       }
               }//close lep-tau selection
              for(int t1=0   ;dilLep2>=0 && t1<(int)taus   .size()&& !higgsCandId;t1++){
                       for(int t2=t1+1;dilLep2>=0 && t2<(int)taus   .size()&& !higgsCandId;t2++){

                               if(taus[t1].pt()<15 || fabs(taus[t1].eta())>2.3) continue;
                               if(taus[t2].pt()<15 || fabs(taus[t2].eta())>2.3) continue;
                               if(taus[t1].id*taus[t2].id!=(charge*225)) continue;
                               if((taus[t1].pt()+taus[t2].pt()) < sumPt.at(3) ) continue;

                               if(deltaR(taus[t1]   , leptons[dilLep1])<0.1) continue;
                               if(deltaR(taus[t1]   , leptons[dilLep2])<0.1) continue;
                               if(deltaR(taus[t2]   , leptons[dilLep1])<0.1) continue;
                               if(deltaR(taus[t2]   , leptons[dilLep2])<0.1) continue;
                               if(deltaR(taus[t1]   , taus[t2        ])<0.1) continue;
                               cout << "iso TT cut " << isoLep.at(3) << " iso tau 1: " << taus[t1].passId((float)isoLep.at(3)) << " iso tau 2: " << taus[t2].passId((float)isoLep.at(3)) << endl;
                               cout << "tau against ele loose? " << taus[t1].passId(llvvTAUID::againstElectronLoose) <<  " - " << taus[t2].passId(llvvTAUID::againstElectronLoose) << endl;

                               if(!taus[t1].passId(llvvTAUID::againstElectronLoose) || !taus[t1].passId((float)isoLep.at(3)) ) continue;
                               if(!taus[t2].passId(llvvTAUID::againstElectronLoose) || !taus[t2].passId((float)isoLep.at(3)) ) continue;
                               //if(!taus[t1].passId(llvvTAUID::againstElectronLoose) || !taus[t1].passId(llvvTAUID::byMediumCombinedIsolationDeltaBetaCorr3Hits) ) continue;
                               //if(!taus[t2].passId(llvvTAUID::againstElectronLoose) || !taus[t2].passId(llvvTAUID::byMediumCombinedIsolationDeltaBetaCorr3Hits) ) continue;

                               higgsCandId=taus[t1].id*taus[t2].id;  higgsCandT1=t1; higgsCandT2=t2;
                               cout << "FOUND TT CANDIDATE of charge " << higgsCandId << endl;
                               higgsCand = LorentzVector(taus[t1]+taus[t2]);
                               break;
                       }
               }//close tau-tau selection
	      if(higgsCandMu!=-1)weight *= isMC ? lepEff.getLeptonEfficiency( leptons[higgsCandMu].pt(), leptons[higgsCandMu].eta(), abs(leptons[higgsCandMu].id), abs(leptons[higgsCandMu].id) ==11 ? "loose" : "loose" ).first : 1.0;
	      if(higgsCandEl!=-1)weight *= isMC ? lepEff.getLeptonEfficiency( leptons[higgsCandEl].pt(), leptons[higgsCandEl].eta(), abs(leptons[higgsCandEl].id), abs(leptons[higgsCandEl].id) ==11 ? "loose" : "loose" ).first : 1.0;


                                higgsCandId = abs(higgsCandId);
				passHiggs = higgsCandId>100; 
 
				if( higgsCandId==143 ){ 
					chTags.push_back(chTags[chTags.size()-1] + string("_elmu")); 
                                	HiggsShortId=0+(abs(leptons[dilLep1].id)==13?0:4);
				}
				else if( higgsCandId==165 ){ 
					chTags.push_back(chTags[chTags.size()-1] + string("_elha")); 
					HiggsShortId=1+(abs(leptons[dilLep1].id)==13?0:4);
				}
				else if( higgsCandId==195 ){ 
					chTags.push_back(chTags[chTags.size()-1] + string("_muha")); 
					HiggsShortId=2+(abs(leptons[dilLep1].id)==13?0:4);
				}
				else if( higgsCandId==225 ){ 
					chTags.push_back(chTags[chTags.size()-1] + string("_haha")); 
					HiggsShortId=3+(abs(leptons[dilLep1].id)==13?0:4);
				}
				else if( higgsCandId== 15 ){ 
					chTags.push_back(chTags[chTags.size()-1] + string("_haCtrl"));
				}
				else      
					chTags.push_back(chTags[chTags.size()-1] + string("_none"));

                                //Lepton Veto
                                for(int l1=0;l1<(int)leptons.size();l1++){
                                        if(l1==dilLep1 || l1==dilLep2 || l1==higgsCandMu || l1==higgsCandEl) continue; //lepton already used in the dilepton pair or higgs candidate
                                        passLepVeto = false; break;
                                }
                                for(int t1=0;passLepVeto && t1<(int)taus   .size();t1++){
                                        if(t1==higgsCandT1 || t1==higgsCandT2) continue; //taus already used in the dilepton pair or higgs candidate
                                        if(taus[t1].pt()<20) continue;
                                        if(!taus[t1].passId(llvvTAUID::againstElectronLoose) || 
                                           !taus[t1].passId(llvvTAUID::againstMuonLoose2)    || 
                                           !taus[t1].passId(llvvTAUID::byMediumCombinedIsolationDeltaBetaCorr3Hits)  ) continue;
                                        passLepVeto = false; break;
                                }

                                //b-jet veto
                                for(int j1=0;j1<(int)bjets.size();j1++){
                                        if(dilLep1    !=-1 && deltaR(bjets[j1]   , leptons[dilLep1    ])>0.4){passBJetVeto=false; break;}
                                        if(dilLep2    !=-1 && deltaR(bjets[j1]   , leptons[dilLep2    ])>0.4){passBJetVeto=false; break;}
                                        if(higgsCandMu!=-1 && deltaR(bjets[j1]   , leptons[higgsCandMu])>0.4){passBJetVeto=false; break;}
                                        if(higgsCandEl!=-1 && deltaR(bjets[j1]   , leptons[higgsCandEl])>0.4){passBJetVeto=false; break;}
                                        if(higgsCandT1!=-1 && deltaR(bjets[j1]   , taus   [higgsCandT1])>0.4){passBJetVeto=false; break;}
                                        if(higgsCandT2!=-1 && deltaR(bjets[j1]   , taus   [higgsCandT2])>0.4){passBJetVeto=false; break;}
                                }

                                for(int j1=0;j1<(int)jets.size();j1++){
                                        if(dilLep1    !=-1 && deltaR(jets[j1]   , leptons[dilLep1    ])<0.4) continue;
                                        if(dilLep2    !=-1 && deltaR(jets[j1]   , leptons[dilLep2    ])<0.4) continue;
                                        if(higgsCandMu!=-1 && deltaR(jets[j1]   , leptons[higgsCandMu])<0.4) continue;
                                        if(higgsCandEl!=-1 && deltaR(jets[j1]   , leptons[higgsCandEl])<0.4) continue;
                                        if(higgsCandT1!=-1 && deltaR(jets[j1]   , taus   [higgsCandT1])<0.4) continue;
                                        if(higgsCandT2!=-1 && deltaR(jets[j1]   , taus   [higgsCandT2])<0.4) continue;
                                        NCleanedJet++;
                                }

				return higgsCand;

       }//close the function 



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
	//bool runLoosePhotonSelection(false); 
	//bool Cut_tautau_MVA_iso = true;
        ////***BOOLEAN variableS above not used for the moment,commented for sake of cleaning***////

	float minJetPtToApply(30);

	std::vector<int> jacknifeCfg=runProcess.getParameter<std::vector<int> >("jacknife");
	int jacknife(jacknifeCfg[0]), jacks(jacknifeCfg[1]);
	if(jacknife>0 && jacks>0) cout << "Jacknife will be applied to every " << jacknife << " out of " << jacks << " events" << endl;

	std::vector<std::string> urls=runProcess.getParameter<std::vector<std::string> >("input");
	TString url = TString(urls[0]);
	TString baseDir    = runProcess.getParameter<std::string>("dirName");
	TString outFileUrl(gSystem->BaseName(url));
	outFileUrl.ReplaceAll(".root","");
	if(mctruthmode!=0) { outFileUrl += "_filt"; outFileUrl += mctruthmode; }
	TString outdir=runProcess.getParameter<std::string>("outdir");
	TString outUrl( outdir );
	gSystem->Exec("mkdir -p " + outUrl);
	bool filterOnlyEE(false), filterOnlyMUMU(false);
	if(!isMC)
	{
		if(url.Contains("DoubleEle")) filterOnlyEE=true;
		if(url.Contains("DoubleMu"))  filterOnlyMUMU=true;
	}
	bool isSingleMuPD(!isMC && url.Contains("SingleMu"));  
	bool isV0JetsMC(isMC && (url.Contains("DYJetsToLL_50toInf") || url.Contains("WJets")));
	bool isSignal(isMC && (url.Contains("VBFNLO") || url.Contains("lljj")) );

	TString outTxtUrl= outUrl + "/" + outFileUrl + ".txt";
	FILE* outTxtFile = NULL;
	if(!isMC)outTxtFile = fopen(outTxtUrl.Data(), "w");
	printf("TextFile URL = %s\n",outTxtUrl.Data());


	//summary ntuple
	TString summaryTupleVarNames("ch:weight:nInitEvent:mjj:detajj:spt:setajj:dphijj:ystar:hardpt:fisher:llr:mva:ystar3:maxcjpt:ncjv:htcjv:ncjv15:htcjv15");
	TNtuple *summaryTuple = new TNtuple("ewkzp2j","ewkzp2j",summaryTupleVarNames);
	Float_t summaryTupleVars[summaryTupleVarNames.Tokenize(":")->GetEntriesFast()];
	summaryTuple->SetDirectory(0);

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

	TH1 *h=mon.addHistogram( new TH1F ("eventflow", ";;Events", 20,0,20) );
	h->GetXaxis()->SetBinLabel(1,"InitialEv");
	h->GetXaxis()->SetBinLabel(2,"Nlep#geq2");
	h->GetXaxis()->SetBinLabel(3,"Zmass");
	h->GetXaxis()->SetBinLabel(4,"Zkin");
	h->GetXaxis()->SetBinLabel(5,"Nlep+Ntau#geq4"); 
	h->GetXaxis()->SetBinLabel(6,"Higgs Cand");
	h->GetXaxis()->SetBinLabel(7,"Lep Veto");
	h->GetXaxis()->SetBinLabel(8,"Btag Veto");
	h->GetXaxis()->SetBinLabel(10,"mm_em");
	h->GetXaxis()->SetBinLabel(11,"mm_et");
	h->GetXaxis()->SetBinLabel(12,"mm_mt");
	h->GetXaxis()->SetBinLabel(13,"mm_tt");
	h->GetXaxis()->SetBinLabel(14,"ee_em");
	h->GetXaxis()->SetBinLabel(15,"ee_et");
	h->GetXaxis()->SetBinLabel(16,"ee_mt");
	h->GetXaxis()->SetBinLabel(17,"ee_tt");

	/*TH1 *h1=mon.addHistogram( new TH1F ("failreason", ";;Events", 20,0,20) );
	h1->GetXaxis()->SetBinLabel(1,"");
	h1->GetXaxis()->SetBinLabel(2,"");
	h1->GetXaxis()->SetBinLabel(3,"");
	h1->GetXaxis()->SetBinLabel(4,"");
	h1->GetXaxis()->SetBinLabel(6,"Nlep#geq2");
	h1->GetXaxis()->SetBinLabel(7,"Zmass");
	h1->GetXaxis()->SetBinLabel(8,"Zkin");
	h1->GetXaxis()->SetBinLabel(9,"Nlep+Ntau#geq4");
	h1->GetXaxis()->SetBinLabel(10,"HiggsCand"); 
	h1->GetXaxis()->SetBinLabel(11,"leptVeto");
	h1->GetXaxis()->SetBinLabel(12,"bjetVeto");*/

	//PLOTS ADDED FOR CHECK!
	/*mon.addHistogram(new TH1F ("drtaulep","#DeltaR(#tau,lep)",50,0,5));                         //filled before vetoing
	mon.addHistogram(new TH1F ("decaymode","decaymodefinding",2,0,2));
	mon.addHistogram(new TH1F ("antimuon2","antimuon2",2,0,2));
	mon.addHistogram(new TH1F ("chargeEM","id1*id2 of the tau decay products",452,-226,226));
	mon.addHistogram(new TH1F ("chargeLT","id1*id2 of the tau decay products",452,-226,226));
	mon.addHistogram(new TH1F ("chargeTT","id1*id2 of the tau decay products",452,-226,226));
	mon.addHistogram(new TH1F ("isomu","RelIso(#mu)",50,-0.5,1.5));
	mon.addHistogram(new TH1F("isoele","RelIso(ele)",50,-0.5,1.5));
        */

	mon.addHistogram( new TH1F("pthat",";#hat{p}_{T} [GeV];Events",50,0,1000) );
	mon.addHistogram( new TH1F("nup",";NUP;Events",10,0,10) );
	mon.addHistogram( new TH1F("nupfilt",";NUP;Events",10,0,10) );

	//pileup control
	mon.addHistogram( new TH1F( "nvtx",";Vertices;Events",50,-0.5,49.5) ); 
	mon.addHistogram( new TH1F( "nvtxraw",";Vertices;Events",50,-0.5,49.5) ); 
	mon.addHistogram( new TH1F( "rho",";#rho;Events",50,0,25) ); 
	mon.addHistogram( new TH1F( "rho25",";#rho(#eta<2.5);Events",50,0,25) ); 

	//lepton control
	mon.addHistogram( new TH1F( "nlep", ";nlep;Events", 10,0,10) );
	mon.addHistogram( new TH1F( "leadpt", ";p_{T}^{l};Events", 100,0,100) );
	mon.addHistogram( new TH1F( "leadeta", ";#eta^{l};Events", 50,-2.6,2.6) );
	mon.addHistogram( new TH1F( "trailerpt", ";p_{T}^{l};Events", 100,0,100) );
	mon.addHistogram( new TH1F( "trailereta", ";#eta^{l};Events", 50,-2.6,2.6) );
	mon.addHistogram( new TH1F( "leppt", ";p_{T}^{l};Events", 100,0,100) );

	//tau control
	mon.addHistogram( new TH1F( "ntaus",      ";ntaus;Events", 10,-0.5,9.5) );
	mon.addHistogram( new TH1F( "tauleadpt",  ";p_{T}^{#tau};Events", 100,0,100) );
	mon.addHistogram( new TH1F( "tauleadeta", ";#eta^{#tau};Events", 50,-2.6,2.6) );
	mon.addHistogram( new TH1F( "taupt",  ";p_{T}^{#tau};Events", 100,0,100) );
	mon.addHistogram( new TH1F( "taucharge",  ";p_{T}^{#tau};Events", 5,-2,2) );
	mon.addHistogram( new TH1F( "taudz",      ";dz^{#tau};Events", 50,0,10) );
	mon.addHistogram( new TH1F( "tauvz",      ";vz^{#tau};Events", 50,0,10) );

	//bjets control
	mon.addHistogram( new TH1F( "nbjets",      ";ntaus;Events", 6,0,6) );
	mon.addHistogram( new TH1F( "bjetpt",  ";p_{T}^{bjet};Events", 100,0,100) );
	mon.addHistogram( new TH1F( "bjetcsv", ";#eta^{#tau};Events", 50,0, 1) );

	//boson control
	mon.addHistogram( new TH1F( "qt",      ";p_{T}^{#gamma} [GeV];Events",500,0,1500));
	mon.addHistogram( new TH1F( "zpt",     ";p_{T}^{ll};Events", 100,0,100) );
	mon.addHistogram( new TH1F( "zptNM1",  ";p_{T}^{ll};Events", 100,0,100) );
	mon.addHistogram( new TH1F( "zeta",    ";#eta^{ll};Events", 50,-10,10) );
	mon.addHistogram( new TH1F( "zetaNM1", ";#eta^{ll};Events", 50,-10,10) );
	mon.addHistogram( new TH1F( "zy",      ";y^{ll};Events", 50,-6,6) );
	mon.addHistogram( new TH1F( "rapidity",";y^{ll};Events", 50,0,2) );
	mon.addHistogram( new TH1F( "zyNM1",   ";y^{ll};Events", 50,-6,6) );
	mon.addHistogram( new TH1F( "zmass",   ";M^{ll};Events", 60,60,120) );

	//higgs control
	mon.addHistogram( new TH1F( "Apt",      ";p_{T}^{A} [GeV];Events",25,0,100));
	mon.addHistogram( new TH1F( "Amass",    ";M^{A} [GeV];Events",25,0,300));
	mon.addHistogram( new TH1F( "Amasssvfit",    ";M^{A} [GeV];Events",25,0,300));
	mon.addHistogram( new TH1F( "Amet",    ";MET [GeV];Events",20,0,200));
	mon.addHistogram( new TH1F( "Anjets",   ";NJets;Events",10,-0.5,9.5));
	mon.addHistogram( new TH1F( "Hmass",    ";M^{H} [GeV];Events",50,0,600));
	mon.addHistogram( new TH1F( "Hpt",      ";p_{T}^{H} [GeV];Events",25,0,100));
	mon.addHistogram( new TH1F( "Hmasssvfit",    ";M^{H} [GeV];Events",50,0,600));
        double xbin[5]={0,25,50,100,200}; 
        (TH2F*)mon.addHistogram(new TProfile2D("vismass2D",      ";M_{A,VIS}; M_{H,VIS}", 4, xbin, 4, 0, 600));
        (TH2F*)mon.addHistogram(new TProfile2D("svfitmass2D",      ";M_{A,SVFit}; M_{H,SVFit}", 4, xbin, 4, 0, 600));
        

        std::vector<llvvTAUID> tauIDiso;
	tauIDiso.push_back(llvvTAUID::byLooseCombinedIsolationDeltaBetaCorr3Hits);
	tauIDiso.push_back(llvvTAUID::byMediumCombinedIsolationDeltaBetaCorr3Hits);
	tauIDiso.push_back(llvvTAUID::byLooseIsolationMVA3oldDMwLT);
	tauIDiso.push_back(llvvTAUID::byMediumIsolationMVA3oldDMwLT);

	std::vector<float>    optim_Cuts_sumptEM , optim_Cuts_sumptET , optim_Cuts_sumptMT , optim_Cuts_sumptTT ;
	std::vector<float>    optim_Cuts_lepIsoEM, optim_Cuts_lepIsoET, optim_Cuts_lepIsoMT, optim_Cuts_lepIsoTT;
	std::vector<float>    optim_Cuts_charge;

        for(float charge=-1;charge<=1;charge+=2){
	for(float lepIsoEM=0.30;lepIsoEM>=0.10;lepIsoEM-=0.1)
	{
		for(float lepIsoET=0.30;lepIsoET>=0.10;lepIsoET-=0.1)
		{
			for(float lepIsoMT=0.30;lepIsoMT>=0.10;lepIsoMT-=0.1)
			{
				for(float lepIsoTT=0;lepIsoTT<tauIDiso.size();lepIsoTT++)
				{
					for(float sumptEM=20;sumptEM<=100;sumptEM+=20)
					{
						for(float sumptET=sumptEM;sumptET<=100;sumptET+=20)
						{
							for(float sumptMT=sumptET;sumptMT<=100;sumptMT+=20)
							{
								for(float sumptTT=sumptMT+20;sumptTT<=100;sumptTT+=20)
								{
										optim_Cuts_charge.push_back(charge);									
										optim_Cuts_lepIsoEM.push_back(lepIsoEM);
										optim_Cuts_lepIsoET.push_back(lepIsoET);
										optim_Cuts_lepIsoMT.push_back(lepIsoMT);
										optim_Cuts_lepIsoTT.push_back(float(tauIDiso.at(lepIsoTT)));
										optim_Cuts_sumptEM.push_back(sumptEM);
										optim_Cuts_sumptET.push_back(sumptET);
										optim_Cuts_sumptMT.push_back(sumptMT);
										optim_Cuts_sumptTT.push_back(sumptTT);
								}
							}
						}
					}
				}
			}
		}
	}}

	TH2F* Hoptim_cuts  =(TH2F*)mon.addHistogram(new TProfile2D("optim_cut",      ";cut index;variable",       optim_Cuts_sumptEM.size(),0,optim_Cuts_sumptEM.size(), 10, 0, 10)) ;
	Hoptim_cuts->GetYaxis()->SetBinLabel(1, "SS/OS"); 
	Hoptim_cuts->GetYaxis()->SetBinLabel(2, "lepIso_EM<"); 
	Hoptim_cuts->GetYaxis()->SetBinLabel(3, "lepIso_ET<");
	Hoptim_cuts->GetYaxis()->SetBinLabel(4, "lepIso_MT<"); 
	Hoptim_cuts->GetYaxis()->SetBinLabel(5, "lepIso_TT<"); 
	Hoptim_cuts->GetYaxis()->SetBinLabel(6, "sumPt_EM>"); 
	Hoptim_cuts->GetYaxis()->SetBinLabel(7, "sumPt_ET>"); 
	Hoptim_cuts->GetYaxis()->SetBinLabel(8, "sumPt_MT>"); 
	Hoptim_cuts->GetYaxis()->SetBinLabel(9, "sumPt_TT>"); 

	for(unsigned int index=0;index<optim_Cuts_sumptEM.size();index++){
		Hoptim_cuts->Fill(index,0.0,optim_Cuts_charge[index]); 
		Hoptim_cuts->Fill(index,1.0,optim_Cuts_lepIsoEM[index]); 
		Hoptim_cuts->Fill(index,2.0,optim_Cuts_lepIsoET[index]); 
		Hoptim_cuts->Fill(index,3.0,optim_Cuts_lepIsoMT[index]); 
		Hoptim_cuts->Fill(index,4.0,(float)optim_Cuts_lepIsoTT[index]); 
		Hoptim_cuts->Fill(index,5.0,optim_Cuts_sumptEM[index]); 
		Hoptim_cuts->Fill(index,6.0,optim_Cuts_sumptET[index]); 
		Hoptim_cuts->Fill(index,7.0,optim_Cuts_sumptMT[index]); 
		Hoptim_cuts->Fill(index,8.0,optim_Cuts_sumptTT[index]); 
	}

	TH1F* Hoptim_systs     =  (TH1F*) mon.addHistogram( new TH1F ("optim_systs"    , ";syst;", nvarsToInclude,0,nvarsToInclude) ) ;
	for(size_t ivar=0; ivar<nvarsToInclude; ivar++)
	{
		Hoptim_systs->GetXaxis()->SetBinLabel(ivar+1, varNames[ivar]);
		mon.addHistogram( new TH2F (TString("svfit_shapes")+varNames[ivar],";cut index;|M_{A}|;Events",optim_Cuts_sumptEM.size(),0,optim_Cuts_sumptEM.size(),25,0,250) );
	}

	//##############################################
	//######## GET READY FOR THE EVENT LOOP ########
	//##############################################

	fwlite::ChainEvent ev(urls);
	const Int_t totalEntries= ev.size();

	//MC normalization (to 1/pb)
	float nInitEvent=1.0;
	if(isMC){
		nInitEvent = (float)utils::getMergeableCounterValue(urls, "startCounter");
	}
	double xsecWeight = xsec/nInitEvent;
	if(!isMC) xsecWeight=1.0;

	//jet energy scale and uncertainties 
	TString jecDir = runProcess.getParameter<std::string>("jecDir");
	gSystem->ExpandPathName(jecDir);
	FactorizedJetCorrector *jesCor        = utils::cmssw::getJetCorrector(jecDir,isMC);
	JetCorrectionUncertainty *totalJESUnc = new JetCorrectionUncertainty((jecDir+"/MC_Uncertainty_AK5PFchs.txt").Data());

	//muon energy scale and uncertainties
	MuScleFitCorrector *muCor=getMuonCorrector(jecDir,url);

	//pdf info
	//  PDFInfo *mPDFInfo=0;
	//  if(isMC)
	//    {
	//      TString pdfUrl(url);
	//      pdfUrl.ReplaceAll(".root","_pdf.root");
	//      pdfUrl.ReplaceAll("/MC","/pdf/MC");
	//      mPDFInfo=new PDFInfo(pdfUrl,"cteq66.LHgrid");
	//      for(int i=0; i<mPDFInfo->numberPDFs(); i++)
	//	{
	//	  TString var("_"); var+=i;
	//	  mon.addHistogram( new TH1F("vbfcandjetdeta"+var    , ";|#Delta #eta|;Jets",                        50,0,10) );
	//	  mon.addHistogram( new TH1F("vbfcandjet1eta"+var    , ";#eta;Jets",                                 50,0,5) );
	//	  mon.addHistogram( new TH1F("vbfcandjet1pt"+var     , ";p_{T} [GeV];Jets",                        50,0,500) );
	//	  mon.addHistogram( new TH1F("vbfcandjet2eta"+var    , ";#eta;Jets",                                 50,0,5) );
	//	  mon.addHistogram( new TH1F("vbfcandjet2pt"+var     , ";p_{T} [GeV];Jets",                        50,0,500) );
	//	}
	//    }

	//pileup weighting: based on vtx for now...
	edm::LumiReWeighting* LumiWeights = NULL;
	utils::cmssw::PuShifter_t PuShifters;
	double PUNorm[] = {1,1,1};
	if(isMC){
		std::vector<double> dataPileupDistributionDouble = runProcess.getParameter< std::vector<double> >("datapileup");
		std::vector<float> dataPileupDistribution; for(unsigned int i=0;i<dataPileupDistributionDouble.size();i++){dataPileupDistribution.push_back(dataPileupDistributionDouble[i]);}
		std::vector<float> mcPileupDistribution;
		utils::getMCPileupDistribution(ev,dataPileupDistribution.size(), mcPileupDistribution);
		while(mcPileupDistribution.size()<dataPileupDistribution.size())  mcPileupDistribution.push_back(0.0);
		while(mcPileupDistribution.size()>dataPileupDistribution.size())dataPileupDistribution.push_back(0.0);

		LumiWeights= new edm::LumiReWeighting(mcPileupDistribution,dataPileupDistribution);
		PuShifters=utils::cmssw::getPUshifters(dataPileupDistribution,0.05);
		utils::getPileupNormalization(mcPileupDistribution, PUNorm, LumiWeights, PuShifters);
	}

	gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE

	FILE* outTxtEvents = NULL;
	outTxtEvents = fopen(outTxtUrl.Data(), "w");
	bool examineThisEvent=false;

	//##############################################
	//########           EVENT LOOP         ########
	//##############################################
	//loop on all the events
	printf("Progressing Bar     :0%%       20%%       40%%       60%%       80%%       100%%\n");
	printf("Scanning the ntuple :");
	DuplicatesChecker duplicatesChecker;
	int nDuplicates(0);
	int step(totalEntries/50);
	for( int iev=0; iev<totalEntries; iev++ ){
		if(iev%step==0){printf(".");fflush(stdout);}
		if(!isMC && jacknife>0 && jacks>0 && iev%jacks==(uint)jacknife) continue;

		//##############################################   EVENT LOOP STARTS   ##############################################
		//load the event content from tree
		if(examineThisEvent) cout << "*** EVENT ***" << iev << endl;
		ev.to(iev);

		//int FaillingReason = 0;
		//    if(!isMC && duplicatesChecker.isDuplicate( ev.eventAuxiliary().run(),ev.eventAuxiliary().luminosityBlock(),ev.eventAuxiliary().event()) ) { nDuplicates++; continue; }
		//      edm::TriggerResultsByName tr = ev.triggerResultsByName("DataAna");      if(!tr.isValid()){printf("TR is invalid\n");continue;}
		//      bool passFilter=false;
		//      for(unsigned int i=0;i<tr.size()-1;i++){
		//         if(!tr.accept(i))continue;
		//         passFilter=true;
		//         printf("Path %3i %50s --> %1i\n",i, tr.triggerName(i).c_str(),tr.accept(i));
		//      }fflush(stdout);
		//      if(!passFilter){cout<<"rejected by the producer\n"; FaillingReason=1;}

		int nvtx = 0;
		fwlite::Handle< int > nvtxHandle;
		nvtxHandle.getByLabel(ev, "llvvObjectProducersUsed", "nvtx");
		if(nvtxHandle.isValid()){ nvtx = *nvtxHandle;}

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
		if(!tauCollHandle.isValid()){printf("llvvLeptonCollection Object NotFound\n");  continue;}
		llvvTauCollection taus = *tauCollHandle;

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

		//pileup weight
		float weight = 1.0;
		double TotalWeight_plus = 1.0;
		double TotalWeight_minus = 1.0;
		float puWeight(1.0);

                //**********************************************************************************************************************************
                //THIS MUST BE UNCOMMENTED FOR THE NORMAL ANALYSIS !!!
		if(isMC){
			puWeight          = LumiWeights->weight(genEv.ngenITpu) * PUNorm[0];
			weight            = xsecWeight*puWeight;
			TotalWeight_plus  = PuShifters[utils::cmssw::PUUP  ]->Eval(genEv.ngenITpu) * (PUNorm[2]/PUNorm[0]);
			TotalWeight_minus = PuShifters[utils::cmssw::PUDOWN]->Eval(genEv.ngenITpu) * (PUNorm[1]/PUNorm[0]);
		}
                //**********************************************************************************************************************************

		//
		// BELOW FOLLOWS THE ANALYSIS OF THE MAIN SELECTION WITH N-1 PLOTS
		//

		//
		// LEPTON ANALYSIS
		//
		llvvLeptonCollection selLeptons;
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
			selLeptons.push_back(leptons[ilep]);
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
			float relIso1 = utils::cmssw::relIso(selLeptons[l1], rho);
			if(examineThisEvent) cout << "ID/Pt/ISO lep1: " << selLeptons[l1].id <<"/"<<selLeptons[l1].pt()<<"/"<<relIso1<< endl;
			if( relIso1>0.30 ) continue;
			for(unsigned int l2=l1+1;l2<selLeptons.size();l2++){
				float relIso2 = utils::cmssw::relIso(selLeptons[l2], rho);
				if(examineThisEvent) cout << "ID/Pt/ISO lep2: " << selLeptons[l2].id <<"/"<<selLeptons[l2].pt()<<"/"<<relIso2<< endl;
				if(fabs(selLeptons[l1].id)!=fabs(selLeptons[l2].id)) continue; //only consider same flavor lepton pairs
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

                //**********************************************************************************************************************************
                //THIS MUST BE UNCOMMENTED FOR THE NORMAL ANALYSIS !!!
		//apply data/mc correction factors
		if(dilLep1>=0)weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[dilLep1].pt(), selLeptons[dilLep1].eta(), abs(selLeptons[dilLep1].id),  abs(selLeptons[dilLep1].id) ==11 ? "loose" : "loose" ).first : 1.0;
		if(dilLep2>=0)weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[dilLep2].pt(), selLeptons[dilLep2].eta(), abs(selLeptons[dilLep2].id),  abs(selLeptons[dilLep2].id) ==11 ? "loose" : "loose" ).first : 1.0;
                //**********************************************************************************************************************************

		std::vector<TString> chTags;
		bool isDileptonCandidate = false;
		chTags.push_back("all");
		if( abs(dilId)==121 && eeTrigger  ){ chTags.push_back("ee"); isDileptonCandidate=true; }
		if( abs(dilId)==169 && mumuTrigger){ chTags.push_back("mumu"); isDileptonCandidate=true; }
		if( !isDileptonCandidate           ) chTags.push_back("ct");
		
		bool passZpt = (zll.pt()>20);
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
				minDRlj = TMath::Min( minDRlj, deltaR(jets[ijet],selLeptons[ilep]) );
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
		llvvTauCollection selTaus;
		if(examineThisEvent) cout << "tau size " << taus.size() << endl;
		for(size_t itau=0; itau<taus.size(); itau++){
			llvvTau& tau = taus[itau];
			if(examineThisEvent) cout << "tau n: " << itau << " pt/eta" << tau.pt() << "/" << fabs(tau.eta()) << endl;
			//if(examineThisEvent) cout << "muonLoose/DM " << tau.passId(llvvTAUID::againstMuonLoose2)  << "/" << tau.passId(llvvTAUID::decayModeFindingNewDMs) << endl;
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

			selTaus.push_back(tau);         
		}
		if(examineThisEvent) cout << "Tau size is: " << selTaus.size() << endl;
		if(examineThisEvent) cout << "***TAU ANALYSIS IS FINISHED***" << endl;

		//
		// HIGGS ANALYSIS
		//
		if(examineThisEvent) cout << "THE HIGGS ANALYSIS STARTS" << endl;

                float weightBeforeLepCorr = weight; //save the weight before the lepton corrections.
                                                    //necessary for the optimization.




		LorentzVector higgsCand;
		LorentzVector higgsCandH;
                std::vector<float> sumPt;
                std::vector<float> isoLep;
                float charge = -1;
                sumPt.push_back(25); sumPt.push_back(30); sumPt.push_back(45); sumPt.push_back(70);
                isoLep.push_back(0.3); isoLep.push_back(0.3); isoLep.push_back(0.3); isoLep.push_back(llvvTAUID::byMediumCombinedIsolationDeltaBetaCorr3Hits);
                int higgsCandId, higgsCandMu, higgsCandEl, higgsCandT1, higgsCandT2, HiggsShortId;
                int NCleanedJet;
                bool passHiggs, passLepVeto, passBJetVeto;

                higgsCand = buildCandidates( selLeptons, selTaus,
                                             selJets, selBJets,      
                                             higgsCandId, higgsCandMu, higgsCandEl, higgsCandT1, higgsCandT2, 
                                             passHiggs, HiggsShortId, chTags,
                                             passLepVeto, passBJetVeto, NCleanedJet,
                                             isoLep, sumPt, charge, 
                                             dilLep1, dilLep2, rho, weight, isMC, lepEff);  
			
                higgsCandH = higgsCand + zll;


		//SVFIT MASS
                SVFitBooster svfitbooster; //initialize the svfit booster (need to be done once per event)
		double diTauMass = -1;
		if(passZpt && passZeta && passHiggs && passLepVeto && passBJetVeto){
                   diTauMass = svfitbooster.getSVFit(met, selLeptons, selTaus, higgsCandMu, higgsCandEl, higgsCandT1, higgsCandT2);  //compute svfit mass in a smart way
                  if(diTauMass<0)diTauMass = higgsCand.mass();
		}

		LorentzVector higgsCand_SVFit;
		LorentzVector higgsCandH_SVFit;

                higgsCand_SVFit = higgsCand;
                float SVFitEnergy = TMath::Sqrt(higgsCand.Px()*higgsCand.Px()+higgsCand.Py()*higgsCand.Py()+higgsCand.Pz()*higgsCand.Pz()+diTauMass*diTauMass);
                higgsCand_SVFit.SetE(SVFitEnergy); 
                higgsCandH_SVFit = higgsCand_SVFit + zll;

		//
		// NOW FOR THE CONTROL PLOTS
		//

		mon.fillHisto("eventflow"      ,   chTags,                 0, weight);
		if(selLeptons.size()>=2){
			//mon.fillHisto("nlep"           ,   chTags, selLeptons.size(), weight);
			mon.fillHisto("eventflow"      ,   chTags,                 1, weight);
			if(passZmass){
				mon.fillHisto("eventflow"   ,   chTags,                 2, weight);

				//pu control
				mon.fillHisto("nvtx"        ,   chTags, nvtx,      weight);
				mon.fillHisto("nvtxraw"     ,   chTags, nvtx,      weight/puWeight);
				mon.fillHisto("rho"         ,   chTags, rho,       weight);
				mon.fillHisto("rho25"       ,   chTags, rho25,     weight);

				//Z kinematics control
				mon.fillHisto("leadpt"      ,   chTags, leadingLep.pt(), weight);      
				mon.fillHisto("trailerpt"   ,   chTags, trailerLep.pt(), weight);      
				mon.fillHisto("leadeta"     ,   chTags, leadingLep.eta(), weight);      
				mon.fillHisto("trailereta"  ,   chTags, trailerLep.eta(), weight);      

				//analyze dilepton kinematics
				mon.fillHisto("zpt"         ,   chTags, zll.pt(),      weight);      
				mon.fillHisto("zmass"       ,   chTags, zll.mass(),    weight);  
				mon.fillHisto("zeta"        ,   chTags, zll.eta(),     weight);
				mon.fillHisto("zy"          ,   chTags, zll.Rapidity(),weight);

				if(passZpt && passZeta){
					mon.fillHisto("eventflow",   chTags,                 3, weight);

					//mon.fillHisto("ntaus"        ,  chTags, selTaus.size(), weight);
					mon.fillHisto("tauleadpt"    ,  chTags, selTaus.size()>0?selTaus[0].pt():-1,  weight);
					mon.fillHisto("tauleadeta"   ,  chTags, selTaus.size()>0?selTaus[0].eta():-10, weight);

					if(selLeptons.size()+selTaus.size()>=4){
						mon.fillHisto("eventflow",   chTags,                 4, weight);

						//SYSTEMATIC STUDY
						for(size_t ivar=0; ivar<nvarsToInclude; ivar++){
							float iweight = weightBeforeLepCorr; //nominal
							//if(ivar==5)                        iweight *= TotalWeight_plus;        //pu up
							//if(ivar==6)                        iweight *= TotalWeight_minus;       //pu down
							//                  if(ivar==7)                        iweight *= Q2Weight_plus;
							//                  if(ivar==8)                        iweight *= Q2Weight_down;
							//                  if(ivar==9)                        iweight *= PDFWeight_plus;
							//                  if(ivar==10)                       iweight *= PDFWeight_down;
							//re-assign the event category;
							for(unsigned int index=0; index<optim_Cuts_sumptEM.size();index++)
							{
								vector<float> sumPtCut;
								vector<float> isoLepCut;
								float charge_OS_SS;
								std::vector<TString> locTags = chTags;

								sumPtCut.push_back(optim_Cuts_sumptEM[index]);
								sumPtCut.push_back(optim_Cuts_sumptET[index]);
								sumPtCut.push_back(optim_Cuts_sumptMT[index]);
								sumPtCut.push_back(optim_Cuts_sumptTT[index]);
								isoLepCut.push_back(optim_Cuts_lepIsoEM[index]);
								isoLepCut.push_back(optim_Cuts_lepIsoET[index]);
								isoLepCut.push_back(optim_Cuts_lepIsoMT[index]);
								isoLepCut.push_back((float)optim_Cuts_lepIsoTT[index]);
								charge_OS_SS = optim_Cuts_charge[index];

								//build the optimal candidate
								LorentzVector higgsCandOpt = buildCandidates( selLeptons, selTaus,
										selJets, selBJets,
										higgsCandId, higgsCandMu, higgsCandEl, higgsCandT1, higgsCandT2, 
										passHiggs, HiggsShortId, locTags,
										passLepVeto, passBJetVeto, NCleanedJet,
										isoLepCut, sumPtCut, charge_OS_SS,
										dilLep1, dilLep2, rho, weight, isMC, lepEff);  
								//
								if(passHiggs && passLepVeto && passBJetVeto){
									SVFitBooster svfitboosterOpt; 
									double diTauMassOpt = -1;
									diTauMassOpt = svfitboosterOpt.getSVFit(met, selLeptons, selTaus, higgsCandMu, higgsCandEl, higgsCandT1, higgsCandT2);  
									if(diTauMassOpt<0) diTauMassOpt = higgsCandOpt.mass();
									mon.fillHisto(TString("svfit_shapes")+varNames[ivar],locTags,index,diTauMassOpt,iweight);}		    
							}
						}
						if(passHiggs){
							mon.fillHisto("eventflow",   chTags,                 5, weight);
							if(passLepVeto){
								mon.fillHisto("eventflow",   chTags,                 6, weight);
								if(passBJetVeto){
									mon.fillHisto("eventflow"	,   chTags,                 7, weight);
									mon.fillHisto("eventflow"	,   chTags,                 9+HiggsShortId, weight);
									mon.fillHisto("Apt"       	, chTags, higgsCand.pt(),    weight);
									mon.fillHisto("Amass"     	, chTags, higgsCand.mass(),  weight);
									mon.fillHisto("Amasssvfit"     	, chTags, higgsCand_SVFit.mass(),  weight);
									mon.fillHisto("Hmass"     	, chTags, higgsCandH.mass(),  weight);
									mon.fillHisto("Hpt"     	, chTags, higgsCandH.pt(),  weight);
									mon.fillHisto("Hmasssvfit"     	, chTags, higgsCandH_SVFit.mass(),  weight);
									mon.fillHisto("vismass2D"	,  chTags, higgsCand.mass(), higgsCandH.mass(), weight);
									mon.fillHisto("svfitmass2D"	,  chTags, higgsCand_SVFit.mass(), higgsCandH_SVFit.mass(), weight);
									mon.fillHisto("Anjets"    	, chTags, NCleanedJet      , weight); 
									mon.fillHisto("Amet"      	, chTags, met.pt()         , weight);
									fprintf(outTxtEvents, "%d %d %d\n",ev.eventAuxiliary().luminosityBlock(),ev.eventAuxiliary().run(),ev.eventAuxiliary().event());
								}//else{mon.fillHisto("failreason",chTags,11,weight);      } //BJETVETO 
							}//else{mon.fillHisto("failreason",chTags,10,weight);      } //LEPVETO
						}//else{mon.fillHisto("failreason",chTags,9,weight);      } //HIGGS
					}//else{mon.fillHisto("failreason",chTags,8,weight);      } //4Lep+Tau
				}//else{mon.fillHisto("failreason",chTags,7,weight);      } //ZKin
			}//else{mon.fillHisto("failreason",chTags,6,weight);      } //ZMass
		}//else{mon.fillHisto("failreason",chTags,5,weight);      } //NLEP  
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
	ofile->Close();

	//save summary tuple
	outUrl.ReplaceAll(".root","_summary.root");
	ofile=TFile::Open(outUrl,"recreate");
	summaryTuple->SetDirectory(ofile);
	summaryTuple->Write();
	ofile->Close();
	if(outTxtFile)fclose(outTxtFile);
	fclose(outTxtEvents);
}  




