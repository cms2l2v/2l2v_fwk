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
      std::map<uint32_t, LorentzVector> precomputed;      
      SVFitBooster(){}
      ~SVFitBooster(){}            
   public:
      LorentzVector getSVFit(llvvMet& met, llvvTauLeptonCollection& selLeptons, int higgsCandL1, int higgsCandL2){       
         if(higgsCandL1<0 || higgsCandL2<0) return LorentzVector(0,0,0,0);
         return LorentzVector(selLeptons[higgsCandL1]+selLeptons[higgsCandL2]); //DEBUG RETURN VIS MASS TO RUN FASTER WHILE DEBUGGING


         uint32_t key = ((higgsCandL1&0x0F)<<16) + (higgsCandL2&0x0F);

         if(precomputed.find(key)==precomputed.end()){ //svfit mass not yet computed for this set of lepton, compute it...
            TMatrixD covMET(2, 2); // PFMET significance matrix
            covMET[0][0] = met.sigx2;
            covMET[0][1] = met.sigxy;
            covMET[1][0] = met.sigxy;
            covMET[1][1] = met.sigy2;

            std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
            measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(abs(selLeptons[higgsCandL1].id)==15?svFitStandalone::kTauToHadDecay:abs(selLeptons[higgsCandL1].id)==11?svFitStandalone::kTauToElecDecay:svFitStandalone::kTauToMuDecay, svFitStandalone::LorentzVector(selLeptons[higgsCandL1].px(), selLeptons[higgsCandL1].py(), selLeptons[higgsCandL1].pz(), selLeptons[higgsCandL1].E()) ));
            measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(abs(selLeptons[higgsCandL2].id)==15?svFitStandalone::kTauToHadDecay:abs(selLeptons[higgsCandL1].id)==11?svFitStandalone::kTauToElecDecay:svFitStandalone::kTauToMuDecay, svFitStandalone::LorentzVector(selLeptons[higgsCandL2].px(), selLeptons[higgsCandL2].py(), selLeptons[higgsCandL2].pz(), selLeptons[higgsCandL2].E()) ));
            SVfitStandaloneAlgorithm algo(measuredTauLeptons, svFitStandalone::Vector(met.px(), met.py(), 0) , covMET, 0);
            algo.addLogM(false);
            algo.fit();
            if(algo.isValidSolution()){             precomputed[key] = algo.fittedDiTauSystem();
            }else{			            precomputed[key] = LorentzVector(0,0,0,0);		
            }
         }
         return precomputed[key];
   }
};


double closestJet(LorentzVectorF& obj, llvvJetExtCollection& selJets, int& closestJetIndex){
   double dRMin = 1E100;  closestJetIndex = -1;
   for(int j=0;j<(int)selJets.size();j++){
      double dR = deltaR(selJets[j], obj);
      if(dR<dRMin){dRMin=dR; closestJetIndex=j;}      
   }
   return dRMin;
}


       LorentzVector buildCandidates( llvvTauLeptonCollection& leptons,  
                                      llvvJetExtCollection& jets, llvvJetExtCollection& bjets,
                                      int& higgsCandId, int& L1, int& L2,
                                      bool& passHiggs, int& HiggsShortId, vector<TString>& chTags,
                                      bool& passLepVeto, bool& passBJetVeto, int& NCleanedJet,
                                      vector<float>& isoLep, vector<float>& sumPt, 
                                      int dilLep1, int dilLep2, double rho, float weight, bool isMC, LeptonEfficiencySF lepEff){


	       LorentzVector higgsCand(0,0,0,0);
               HiggsShortId=-1;
               higgsCandId=0;  

               passHiggs    = false;
               passBJetVeto = true;
               passLepVeto  = true;
               NCleanedJet  = 0;
               string ChannelName = "none";
               string signName = "";

                //FIND the two highest pT leptons not coming from the Z, and with dR>0.1 from all other leptons in the event
                L1=-1, L2=-1;
                for(int l=0   ;l<(int)leptons.size();l++){
                   if(l==dilLep1 || l==dilLep2)continue;
//                   printf("%+2i - pT=%+6.2f eta=%+6.2f phi=%+6.2f\n", leptons[l].id, leptons[l].pt(), leptons[l].eta(), leptons[l].phi());
                   if(deltaR(leptons[l],  leptons[dilLep1])<0.1)continue;
                   if(deltaR(leptons[l],  leptons[dilLep2])<0.1)continue;
                   if(L1<0){L1=l;continue;}
                   if(L2<0 && deltaR(leptons[l],  leptons[L1])>=0.1){L2=l;break;}//ordered in pT, so all done
                }

                if(L1>=0 && L2>=0){
                   higgsCandId=leptons[L1].id*leptons[L2].id;
                   higgsCand = LorentzVector(leptons[L1]+leptons[L2]);
                   if(higgsCandId<0){signName="_OS";}else{signName="_SS";}
                   if(higgsCandId<0){HiggsShortId = 0;}else{HiggsShortId = 8;}
                   if(abs(leptons[dilLep1].id)==13){HiggsShortId += 0;}else{HiggsShortId += 4;}
	
                   // e - mu final state	
                   if(abs(higgsCandId) == 11*13){
                      llvvLepton lep1 = leptons[L1].lep;
                      llvvLepton lep2 = leptons[L2].lep;

                      if((lep1.pt()+lep2.pt()) >= sumPt.at(0)
                      &&  utils::cmssw::relIso(lep1, rho) <= isoLep.at(abs(lep1.id)==11?0:1) 
                      &&  utils::cmssw::relIso(lep2, rho) <= isoLep.at(abs(lep2.id)==11?0:1)
                      ){
                         passHiggs    = true;  
                         ChannelName  = "elmu";
                         HiggsShortId+= 0;
                      }

		   // l - tau final state
                   }else if(abs(higgsCandId) == 11*15 || abs(higgsCandId) == 13*15){
                      llvvTau    tau = abs(leptons[L1].id)==15?leptons[L1].tau:leptons[L2].tau;
                      llvvLepton lep = abs(leptons[L1].id)==15?leptons[L2].lep:leptons[L1].lep;

                      if(tau.pt()>=15 && fabs(tau.eta())<=2.3){
                         float relIso1 = utils::cmssw::relIso(lep, rho);

                         //e-tau
                         if( abs(lep.id)==11 
                             && relIso1<=isoLep.at(abs(lep.id)==11?0:1) 
                             && tau.passId(llvvTAUID::againstElectronTightMVA5) 
                             && tau.passId((int)isoLep.at(2)) 
                             && (tau.pt()+leptons[L1].pt()) >= sumPt.at(0)){

                             passHiggs    = true;
                             ChannelName  = "elha";
                             HiggsShortId+= 1;
                          }

                         //mu-tau
			 if( abs(lep.id)==13 
                             && relIso1<=isoLep.at(abs(lep.id)==11?0:1)
                             && tau.passId(llvvTAUID::againstElectronLoose)
                             && tau.passId(llvvTAUID::againstMuonTight2)
                             && tau.passId((int)isoLep.at(2))
                             && (tau.pt()+lep.pt()) >= sumPt.at(0) ){

                             passHiggs    = true;
                             ChannelName  = "muha";
                             HiggsShortId+= 2;  
                         }    
                      }
                   // tau - tau final state
                   }else if(abs(higgsCandId) == 15*15){
                      llvvTau    tau1 = leptons[L1].tau;
                      llvvTau    tau2 = leptons[L2].tau;

                      if(tau1.pt()>=15 && fabs(tau1.eta())<=2.3
                      && tau2.pt()>=15 && fabs(tau2.eta())<=2.3
                      &&(tau1.pt()+tau2.pt()) >= sumPt.at(0)
                      && tau1.passId(llvvTAUID::againstElectronLoose) && tau1.passId((int)isoLep.at(2)) 
                      && tau2.passId(llvvTAUID::againstElectronLoose) && tau2.passId((int)isoLep.at(2))){

                      passHiggs    = true;  
                      ChannelName  = "haha";
                      HiggsShortId+= 3;
                   }
               }               

	      if(isMC && abs(leptons[L1].id)<15)weight *= lepEff.getLeptonEfficiency( leptons[L1].pt(), leptons[L1].eta(), abs(leptons[L1].id), abs(leptons[L1].id) ==11 ? "loose" : "loose" ).first;
              if(isMC && abs(leptons[L2].id)<15)weight *= lepEff.getLeptonEfficiency( leptons[L2].pt(), leptons[L2].eta(), abs(leptons[L2].id), abs(leptons[L2].id) ==11 ? "loose" : "loose" ).first;

   	      chTags.push_back(chTags[chTags.size()-1] + signName + ChannelName); 


              //Lepton Veto
              for(int l=0;l<(int)leptons.size() && passLepVeto;l++){
                 if(l==dilLep1 || l==dilLep2 || l==L1 || l==L2) continue; //lepton already used in the dilepton pair or higgs candidate
                 if(abs(leptons[l].id)==15){
                    llvvTau    tau = leptons[l].tau;
                    if(tau.pt()<20) continue;
                    if(!tau.passId(llvvTAUID::againstElectronLoose) ||
                       !tau.passId(llvvTAUID::againstMuonLoose2)    ||
                       !tau.passId(llvvTAUID::byMediumCombinedIsolationDeltaBetaCorr3Hits)  ) continue;                    
                    passLepVeto = false; break;
                 }else{
                    passLepVeto = false; break;
                 }
              }

              //b-jet veto
              for(int j1=0;j1<(int)bjets.size();j1++){
                 if(dilLep1    !=-1 && deltaR(bjets[j1]   , leptons[dilLep1 ])>0.4){passBJetVeto=false; break;}
                 if(dilLep2    !=-1 && deltaR(bjets[j1]   , leptons[dilLep2 ])>0.4){passBJetVeto=false; break;}
                 if(L1         !=-1 && deltaR(bjets[j1]   , leptons[L1      ])>0.4){passBJetVeto=false; break;}
                 if(L2         !=-1 && deltaR(bjets[j1]   , leptons[L2      ])>0.4){passBJetVeto=false; break;}
              }

              for(int j1=0;j1<(int)jets.size();j1++){
                 if(dilLep1    !=-1 && deltaR(jets[j1]   , leptons[dilLep1 ])<0.4) continue;
                 if(dilLep2    !=-1 && deltaR(jets[j1]   , leptons[dilLep2 ])<0.4) continue;
                 if(L1         !=-1 && deltaR(jets[j1]   , leptons[L1      ])<0.4) continue;
                 if(L2         !=-1 && deltaR(jets[j1]   , leptons[L2      ])<0.4) continue;
                 NCleanedJet++;
              }
           }else{
              passHiggs    = false;  
              ChannelName = "none";
              HiggsShortId = 0;
              chTags.push_back(chTags[chTags.size()-1] + "_" + ChannelName);
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

	TH1 *h=mon.addHistogram( new TH1F ("eventflow", ";;Events", 28,0,28) );
	h->GetXaxis()->SetBinLabel(1,"InitialEv");
	h->GetXaxis()->SetBinLabel(2,"Nlep#geq2");
	h->GetXaxis()->SetBinLabel(3,"Zmass");
	h->GetXaxis()->SetBinLabel(4,"Zkin");
	h->GetXaxis()->SetBinLabel(5,"Nlep+Ntau#geq4"); 
	h->GetXaxis()->SetBinLabel(6,"Higgs Cand");
	h->GetXaxis()->SetBinLabel(7,"Lep Veto");
	h->GetXaxis()->SetBinLabel(8,"Btag Veto");
	h->GetXaxis()->SetBinLabel(9,"");
	h->GetXaxis()->SetBinLabel(10,"OS em+mm");
	h->GetXaxis()->SetBinLabel(11,"OS et+mm");
	h->GetXaxis()->SetBinLabel(12,"OS mt+mm");
	h->GetXaxis()->SetBinLabel(13,"OS tt+mm");
	h->GetXaxis()->SetBinLabel(14,"OS em+ee");
	h->GetXaxis()->SetBinLabel(15,"OS et+ee");
	h->GetXaxis()->SetBinLabel(16,"OS mt+ee");
	h->GetXaxis()->SetBinLabel(17,"OS tt+ee");
        h->GetXaxis()->SetBinLabel(18,"SS em+mm");
        h->GetXaxis()->SetBinLabel(19,"SS et+mm");
        h->GetXaxis()->SetBinLabel(20,"SS mt+mm");
        h->GetXaxis()->SetBinLabel(21,"SS tt+mm");
        h->GetXaxis()->SetBinLabel(22,"SS em+ee");
        h->GetXaxis()->SetBinLabel(23,"SS et+ee");
        h->GetXaxis()->SetBinLabel(24,"SS mt+ee");
        h->GetXaxis()->SetBinLabel(25,"SS tt+ee");


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

	double xbin[5]={0,25,50,100,200};
	mon.addHistogram( new TH1F( "Apt",                   ";p_{T}^{A} [GeV];Events",25,0,100));
	mon.addHistogram( new TH1F( "Amass",                 ";M^{A} [GeV];Events",20,0,300));
	mon.addHistogram( new TH1F( "Amasssvfit",            ";M^{A} [GeV];Events",20,0,300));
	mon.addHistogram( new TH1F( "Amet",                  ";MET [GeV];Events",20,0,200));
	mon.addHistogram( new TH1F( "Anjets",                ";NJets;Events",10,-0.5,9.5));
	mon.addHistogram( new TH1F( "Hmass",                 ";M^{H} [GeV];Events",40,0,600));
	mon.addHistogram( new TH1F( "Hmass_MA_0_25",         ";M^{H} [GeV];Events",40,0,600));
	mon.addHistogram( new TH1F( "Hmass_MA_25_50",        ";M^{H} [GeV];Events",40,0,600));
	mon.addHistogram( new TH1F( "Hmass_MA_50_100",       ";M^{H} [GeV];Events",40,0,600));
	mon.addHistogram( new TH1F( "Hmass_MA_100_200",      ";M^{H} [GeV];Events",40,0,600));
	mon.addHistogram( new TH1F( "Hpt",                   ";p_{T}^{H} [GeV];Events",25,0,100));
	mon.addHistogram( new TH1F( "Hmasssvfit",            ";M^{H} [GeV];Events",40,0,600));
	mon.addHistogram( new TH1F( "Hmasssvfit_MA_0_25",    ";M^{H} [GeV];Events",40,0,600));
	mon.addHistogram( new TH1F( "Hmasssvfit_MA_25_50",   ";M^{H} [GeV];Events",40,0,600));
	mon.addHistogram( new TH1F( "Hmasssvfit_MA_50_100",  ";M^{H} [GeV];Events",40,0,600));
	mon.addHistogram( new TH1F( "Hmasssvfit_MA_100_200", ";M^{H} [GeV];Events",40,0,600));
	(TH2F*)mon.addHistogram(new TProfile2D("vismass2D",      ";M_{A,VIS}; M_{H,VIS}", 4, xbin, 4, 0, 600));
	(TH2F*)mon.addHistogram(new TProfile2D("svfitmass2D",      ";M_{A,SVFit}; M_{H,SVFit}", 4, xbin, 4, 0, 600));


        std::vector<llvvTAUID> tauIDiso;
        tauIDiso.push_back(llvvTAUID::decayModeFinding);
	tauIDiso.push_back(llvvTAUID::byLooseCombinedIsolationDeltaBetaCorr3Hits);
	tauIDiso.push_back(llvvTAUID::byMediumCombinedIsolationDeltaBetaCorr3Hits);
	tauIDiso.push_back(llvvTAUID::byLooseIsolationMVA3oldDMwLT);
	tauIDiso.push_back(llvvTAUID::byMediumIsolationMVA3oldDMwLT);

	std::vector<float>    optim_Cuts_sumPt;
	std::vector<float>    optim_Cuts_elIso, optim_Cuts_muIso, optim_Cuts_taIso;

 //DEBUG
	for(float elIso=0.40;elIso>=0.10;elIso-=0.1)
	{
		for(float muIso=0.40;muIso>=0.10;muIso-=0.1)
		{
				for(float taIso=0;taIso<tauIDiso.size();taIso++)
				{
					for(float sumPt=20;sumPt<=100;sumPt+=20)
					{
										optim_Cuts_elIso.push_back(elIso);
										optim_Cuts_muIso.push_back(muIso);
										optim_Cuts_taIso.push_back(float(tauIDiso.at(taIso)));
										optim_Cuts_sumPt.push_back(sumPt);
					}
				}
		}
	}
	TH2F* Hoptim_cuts  =(TH2F*)mon.addHistogram(new TProfile2D("optim_cut",      ";cut index;variable",       optim_Cuts_sumPt.size(),0,optim_Cuts_sumPt.size(), 4, 0, 4)) ;
	Hoptim_cuts->GetYaxis()->SetBinLabel(1, "eIso<"); 
	Hoptim_cuts->GetYaxis()->SetBinLabel(2, "muIso<");
	Hoptim_cuts->GetYaxis()->SetBinLabel(3, "tauIso<"); 
	Hoptim_cuts->GetYaxis()->SetBinLabel(4, "sumPt>"); 

	for(unsigned int index=0;index<optim_Cuts_sumPt.size();index++){
		Hoptim_cuts->Fill(index,0.0,optim_Cuts_elIso[index]); 
		Hoptim_cuts->Fill(index,1.0,optim_Cuts_muIso[index]); 
		Hoptim_cuts->Fill(index,2.0,(float)optim_Cuts_taIso[index]); 
		Hoptim_cuts->Fill(index,3.0,optim_Cuts_sumPt[index]); 
	}

	TH1F* Hoptim_systs     =  (TH1F*) mon.addHistogram( new TH1F ("optim_systs"    , ";syst;", nvarsToInclude,0,nvarsToInclude) ) ;
	for(size_t ivar=0; ivar<nvarsToInclude; ivar++)
	{
		Hoptim_systs->GetXaxis()->SetBinLabel(ivar+1, varNames[ivar]);
		mon.addHistogram( new TH2F (TString("svfit_shapes")+varNames[ivar],";cut index;|M_{A}|;Events",optim_Cuts_sumPt.size(),0,optim_Cuts_sumPt.size(),25,0,250) );
		mon.addHistogram( new TH2F (TString("FR_closestJetPt")+varNames[ivar],";cut index;jet p_{T} (GeV);Events",optim_Cuts_sumPt.size(),0,optim_Cuts_sumPt.size(),15,0,150) );
	}


        //create a tree and related variables to save higgs candidate info for each cutIndex values
        //this is done only on data 
        unsigned int treeEventId = 0;
        unsigned int treeLumiId  = 0;
        unsigned int treeRunId   = 0;
        unsigned int treeCutIndex= 0;        
        int          treeHiggsId = 0;
        float        treeVisMass = 0;        
        float        treeSVFMass = 0;
        int          treeLeg1Id  = 0;
        float        treeLeg1DR  = 0;
        float        treeLeg1Pt  = 0;
        float        treeLeg1Eta = 0;
        int          treeLeg2Id  = 0;
        float        treeLeg2DR  = 0;
        float        treeLeg2Pt  = 0;
        float        treeLeg2Eta = 0;

                
        TTree* tree = NULL;
        if(!isMC){
           tree = new TTree("CandTree","CandTree");
           tree->Branch("eventId", &treeEventId , string("eventId/i" ).c_str());
           tree->Branch("lumiId" , &treeLumiId  , string("lumiId/i"  ).c_str());
           tree->Branch("runId"  , &treeRunId   , string("runId/i"   ).c_str());
           tree->Branch("cutIndex",&treeCutIndex, string("cutIndex/i").c_str());
           tree->Branch("higgsId", &treeHiggsId , string("higgsId/I" ).c_str());
           tree->Branch("visMass", &treeVisMass , string("visMass/F" ).c_str());
           tree->Branch("svfMass", &treeSVFMass , string("svfMass/F" ).c_str());
           tree->Branch("leg1Id" , &treeLeg1Id  , string("leg1Id/I"  ).c_str());
           tree->Branch("leg1DR" , &treeLeg1DR  , string("leg1DR/F"  ).c_str());
           tree->Branch("leg1Pt" , &treeLeg1Pt  , string("leg1Pt/F"  ).c_str());
           tree->Branch("leg1Eta", &treeLeg1Eta , string("leg1Eta/F" ).c_str());
           tree->Branch("leg2Id" , &treeLeg2Id  , string("leg2Id/I"  ).c_str());
           tree->Branch("leg2DR" , &treeLeg2DR  , string("leg2DR/F"  ).c_str());
           tree->Branch("leg2Pt" , &treeLeg2Pt  , string("leg2Pt/F"  ).c_str());
           tree->Branch("leg2Eta", &treeLeg2Eta , string("leg2Eta/F" ).c_str());
           tree->SetDirectory(NULL);
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
//		cout << "*** EVENT ***" << iev << endl;
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
		if(!tauCollHandle.isValid()){printf("llvvTauCollection Object NotFound\n");  continue;}
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
		if(examineThisEvent) cout << "tau size " << taus.size() << endl;
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

		LorentzVector higgsCand;
		LorentzVector higgsCandH;
                std::vector<float> sumPt;
                std::vector<float> isoLep;
                sumPt.push_back(25); 
                isoLep.push_back(0.2); isoLep.push_back(0.2); isoLep.push_back(llvvTAUID::byMediumCombinedIsolationDeltaBetaCorr3Hits);
                int higgsCandId, higgsCandL1, higgsCandL2, HiggsShortId;
                int NCleanedJetMain;
                bool passHiggsMain, passLepVetoMain, passBJetVetoMain;
		std::vector<TString> chTagsMain=chTags;

                higgsCand = buildCandidates( selLeptons,
                                             selJets, selBJets,      
                                             higgsCandId, higgsCandL1, higgsCandL2, 
                                             passHiggsMain, HiggsShortId, chTagsMain,
                                             passLepVetoMain, passBJetVetoMain, NCleanedJetMain,
                                             isoLep, sumPt, 
                                             dilLep1, dilLep2, rho, weight, isMC, lepEff);  
			
                higgsCandH = higgsCand + zll;


		//SVFIT MASS
                SVFitBooster svfitbooster; //initialize the svfit booster (need to be done once per event)
		double diTauMass = -1;
                LorentzVector diTauSystem; 
		if(passZpt && passZeta && passHiggsMain && passLepVetoMain && passBJetVetoMain){
                  diTauSystem = svfitbooster.getSVFit(met, selLeptons, higgsCandL1, higgsCandL2);  //compute svfit mass in a smart way
                  if(diTauSystem.mass()<=0)diTauSystem = higgsCand;
                  diTauMass = diTauSystem.mass();
		}

		LorentzVector higgsCand_SVFit;
		LorentzVector higgsCandH_SVFit;

                higgsCand_SVFit = diTauSystem;
                higgsCandH_SVFit = higgsCand_SVFit + zll;
                

		//
		// NOW FOR THE CONTROL PLOTS
		//

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

//FIXME
					//mon.fillHisto("ntaus"        ,  chTags, selTaus.size(), weight);
//					mon.fillHisto("tauleadpt"    ,  chTagsMain, selTaus.size()>0?selTaus[0].pt():-1,  weight);
//					mon.fillHisto("tauleadeta"   ,  chTagsMain, selTaus.size()>0?selTaus[0].eta():-10, weight);

					if(selLeptons.size()>=4){
						mon.fillHisto("eventflow",   chTagsMain,                 4, weight);

						if(passHiggsMain){
//								for(uint st=0;st<chTagsMain.size();st++){
//								cout << "chTagsMain " << chTagsMain.at(st) << endl;
//                                                                }
							mon.fillHisto("eventflow",   chTagsMain,                 5, weight);
							if(passLepVetoMain){
								mon.fillHisto("eventflow",   chTagsMain,                 6, weight);
								if(passBJetVetoMain){
									mon.fillHisto("eventflow"	,   chTagsMain,                 7, weight);
									mon.fillHisto("eventflow"	,   chTagsMain,                 9+HiggsShortId, weight);
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

						//SYSTEMATIC STUDY on all events passing the basic preselection
						for(size_t ivar=0; ivar<nvarsToInclude; ivar++){
							//if(ivar==5)                        iweight *= TotalWeight_plus;        //pu up
							//if(ivar==6)                        iweight *= TotalWeight_minus;       //pu down
							//                  if(ivar==7)                        iweight *= Q2Weight_plus;
							//                  if(ivar==8)                        iweight *= Q2Weight_down;
							//                  if(ivar==9)                        iweight *= PDFWeight_plus;
							//                  if(ivar==10)                       iweight *= PDFWeight_down;
							//re-assign the event category;
							for(unsigned int index=0; index<optim_Cuts_sumPt.size();index++)
							{
  							        float iweight = weightBeforeLepCorr; //nominal

								vector<float> sumPtCut;
								vector<float> isoLepCut;
								float charge_OS_SS;

						                int NCleanedJet;
						                bool passHiggs, passLepVeto, passBJetVeto;
								std::vector<TString> locTags = chTags;
                                                                //	std::size_t found1 = (string)locTags[st].find("ee_");
								//	std::size_t found2 = (string)locTags[st].find("mumu_");
								//	if(found1!=std::string::npos)
								//		(string)locTags[st].erase((string)locTags[st].begin()+3, (string)locTags[st].end());
								//	if(found2!=std::string::npos)
								//		(string)locTags[st].erase((string)locTags[st].begin()+5, (string)locTags[st].end());

								sumPtCut.push_back(optim_Cuts_sumPt[index]);
								isoLepCut.push_back(optim_Cuts_elIso[index]);
								isoLepCut.push_back(optim_Cuts_muIso[index]);
								isoLepCut.push_back((float)optim_Cuts_taIso[index]);

								LorentzVector higgsCandOpt = buildCandidates( selLeptons,
										selJets, selBJets,
										higgsCandId, higgsCandL1, higgsCandL2, 
										passHiggs, HiggsShortId, locTags,
										passLepVeto, passBJetVeto, NCleanedJet,
										isoLepCut, sumPtCut,
										dilLep1, dilLep2, rho, iweight, isMC, lepEff);  
								

								if(passHiggs && passLepVeto && passBJetVeto){
									double diTauMassOpt = -1;
									diTauMassOpt = (svfitbooster.getSVFit(met, selLeptons, higgsCandL1, higgsCandL2)).mass();  
									if(diTauMassOpt<=0) diTauMassOpt = higgsCandOpt.mass();
									mon.fillHisto(TString("svfit_shapes")+varNames[ivar],locTags,index,diTauMassOpt,iweight);
          
                                                                        treeEventId  = ev.eventAuxiliary().event();
                                                                        treeLumiId   = ev.eventAuxiliary().luminosityBlock();
                                                                        treeRunId    = ev.eventAuxiliary().run();
                                                                        treeCutIndex = index;
                                                                        treeHiggsId  = higgsCandId;
                                                                        treeVisMass  = diTauMassOpt;
                                                                        treeSVFMass  = higgsCandOpt.mass();                                
                                                                        treeLeg1Id   = -1;
                                                                        treeLeg1DR   = -1;
                                                                        treeLeg1Pt   = -1;
                                                                        treeLeg1Eta  = -1;
                                                                        treeLeg2Id   = -1;
                                                                        treeLeg2DR   = -1;
                                                                        treeLeg2Pt   = -1;
                                                                        treeLeg2Eta  = -1;

                                                                        int closestJetIndex=-1;  double pT=-1;  double dR=-1;
									if(higgsCandL1!=-1){
                                                                           double dRmin = closestJet(selLeptons[higgsCandL1], selJets, closestJetIndex);
									   if(closestJetIndex>=0 && dRmin<0.5){pT=selJets[closestJetIndex].pt(); dR=dRmin;}else{pT=selLeptons[higgsCandL1].pt(); dR=-1;}
 									   mon.fillHisto(TString("FR_closestJetPt")+varNames[ivar],locTags,index,pT,iweight);
                                                                           if(treeLeg1Id==-1){treeLeg1Id=selLeptons[higgsCandL1].id; treeLeg1DR=dR; treeLeg1Pt=pT; treeLeg1Eta=selLeptons[higgsCandL1].eta();
                                                                           }else{             treeLeg2Id=selLeptons[higgsCandL1].id; treeLeg2DR=dR; treeLeg2Pt=pT; treeLeg2Eta=selLeptons[higgsCandL1].eta();}
									}
									if(higgsCandL2!=-1){
                                                                           double dRmin = closestJet(selLeptons[higgsCandL2], selJets, closestJetIndex);
                                                                           if(closestJetIndex>=0 && dRmin<0.5){pT=selJets[closestJetIndex].pt(); dR=dRmin;}else{pT=selLeptons[higgsCandL2].pt(); dR=-1;}
                                                                           mon.fillHisto(TString("FR_closestJetPt")+varNames[ivar],locTags,index,pT,iweight);
                                                                           if(treeLeg1Id==-1){treeLeg1Id=selLeptons[higgsCandL2].id; treeLeg1DR=dR; treeLeg1Pt=pT; treeLeg1Eta=selLeptons[higgsCandL2].eta();
                                                                           }else{             treeLeg2Id=selLeptons[higgsCandL2].id; treeLeg2DR=dR; treeLeg2Pt=pT; treeLeg2Eta=selLeptons[higgsCandL2].eta();}
									}

                                                                        if(ivar==0 && tree)tree->Fill();

								}		    
							}//end of the loop on cutIndex
						}//end of the loop on the systematics

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
        if(tree){tree->SetDirectory(ofile); tree->Write();}
	ofile->Close();
/*
	//save summary tuple
	outUrl.ReplaceAll(".root","_summary.root");
	ofile=TFile::Open(outUrl,"recreate");
	tree->SetDirectory(ofile);
	tree->Write();
	ofile->Close();
	if(outTxtFile)fclose(outTxtFile);
	fclose(outTxtEvents);
*/
}  




