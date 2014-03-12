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


#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
//#include "UserCode/llvv_fwk/interface/DataEventSummaryHandler.h"
#include "UserCode/llvv_fwk/interface/llvvObjects.h"
#include "UserCode/llvv_fwk/interface/TMVAUtils.h"
#include "UserCode/llvv_fwk/interface/LeptonEfficiencySF.h"
#include "UserCode/llvv_fwk/interface/PDFInfo.h"
#include "UserCode/llvv_fwk/interface/MuScleFitCorrector.h"
#include "UserCode/llvv_fwk/interface/GammaWeightsHandler.h"
#include "UserCode/llvv_fwk/interface/BtagUncertaintyComputer.h"
#include "UserCode/llvv_fwk/interface/TopPtWeighter.h"


#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
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

  bool debug = runProcess.getParameter<bool>("debug");
  bool isMC = runProcess.getParameter<bool>("isMC");
  double xsec = runProcess.getParameter<double>("xsec");
  int mctruthmode=runProcess.getParameter<int>("mctruthmode");
  bool runLoosePhotonSelection(false);

  bool Cut_tautau_MVA_iso = true;//currently hardcoded

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
  bool filterOnlyEE(false), filterOnlyEMU(false), filterOnlyMUMU(false), filterOnlySINGLEMU(false);
  if(!isMC)
    {
      if(url.Contains("DoubleEle")) filterOnlyEE=true;
      if(url.Contains("DoubleMu"))  filterOnlyMUMU=true;
      if(url.Contains("MuEG"))      filterOnlyEMU=true;
      if(url.Contains("SingleMu"))  filterOnlySINGLEMU=true;
    }
  bool isSingleMuPD(!isMC && url.Contains("SingleMu"));  
  bool isV0JetsMC(isMC && (url.Contains("DYJetsToLL_50toInf") || url.Contains("WJets")));
  bool isSignal(isMC && (url.Contains("VBFNLO") || url.Contains("lljj")) );

  TString outTxtUrl= outUrl + "/" + outFileUrl + ".txt";
  FILE* outTxtFile = NULL;
  if(!isMC)outTxtFile = fopen(outTxtUrl.Data(), "w");
  printf("TextFile URL = %s\n",outTxtUrl.Data());

  //residual corrections from Z+1 jet recoil
  std::vector<TGraph *> recoilResidualsGr;
  TFile *fIn=TFile::Open("/afs/cern.ch/user/p/psilva/work/CMSSW_5_3_3_patch2/src/CMGTools/HtoZZ2l2nu/data/recoilBalance.root");
  if(fIn){
    TGraph *gr=(TGraph *)fIn->Get("mumurecoilbalancevseta50toInfdydata2mc");
    if(gr) recoilResidualsGr.push_back( gr );
    gr=(TGraph *)fIn->Get("mumurecoilbalancevseta50toInfgdata2mc");
    if(gr) recoilResidualsGr.push_back( gr );
    fIn->Close();
    cout << "Read " << recoilResidualsGr.size() << " residual recoil systematics" << endl;
  }
 

  std::vector<string>  weightsFile = runProcess.getParameter<std::vector<string> >("weightsFile");
  
  TString wFile("");
  if(weightsFile.size()) wFile = TString(gSystem->ExpandPathName(weightsFile[0].c_str()));
  
  //b-tag efficiencies read b-tag efficiency map
  std::map<std::pair<TString,TString>, std::pair<TGraphErrors *,TGraphErrors *> > btagEffCorr;
  if(weightsFile.size() && isMC)
    {
      TString btagEffCorrUrl(wFile); btagEffCorrUrl += "/btagEff.root";
      gSystem->ExpandPathName(btagEffCorrUrl);
      TFile *btagF=TFile::Open(btagEffCorrUrl);
      if(btagF!=0 && !btagF->IsZombie())
	{
	  TList *dirs=btagF->GetListOfKeys();
	  for(int itagger=0; itagger<dirs->GetEntries(); itagger++)
	    {
	      TString iDir(dirs->At(itagger)->GetName());
	      btagEffCorr[ std::pair<TString,TString>(iDir,"b") ] 
		= std::pair<TGraphErrors *,TGraphErrors *>( (TGraphErrors *) btagF->Get(iDir+"/beff"),(TGraphErrors *) btagF->Get(iDir+"/sfb") );
	      btagEffCorr[ std::pair<TString,TString>(iDir,"c") ] 
		= std::pair<TGraphErrors *,TGraphErrors *>( (TGraphErrors *) btagF->Get(iDir+"/ceff"),(TGraphErrors *) btagF->Get(iDir+"/sfc") );
	      btagEffCorr[ std::pair<TString,TString>(iDir,"udsg") ] 
		= std::pair<TGraphErrors *,TGraphErrors *>( (TGraphErrors *) btagF->Get(iDir+"/udsgeff"),(TGraphErrors *) btagF->Get(iDir+"/sfudsg") );
	    }
	}
      cout << btagEffCorr.size() << " b-tag correction factors have been read" << endl;
    }
  
  std::vector<double> btagBins; btagBins.clear();
  btagBins.push_back(35  );
  btagBins.push_back(45  );
  btagBins.push_back(55  );
  btagBins.push_back(65  );
  btagBins.push_back(75  );
  btagBins.push_back(90  );
  btagBins.push_back(110 );
  btagBins.push_back(140 );
  btagBins.push_back(185 );
  btagBins.push_back(235 );
  btagBins.push_back(290 );
  btagBins.push_back(360 );
  btagBins.push_back(450 );
  btagBins.push_back(550 );
  btagBins.push_back(700 );


  

  //summary ntuple
  TString summaryTupleVarNames("ch:weight:nInitEvent:mjj:detajj:spt:setajj:dphijj:ystar:hardpt:fisher:llr:mva:ystar3:maxcjpt:ncjv:htcjv:ncjv15:htcjv15");
  TNtuple *summaryTuple = new TNtuple("ewkzp2j","ewkzp2j",summaryTupleVarNames);
  Float_t summaryTupleVars[summaryTupleVarNames.Tokenize(":")->GetEntriesFast()];
  summaryTuple->SetDirectory(0);

  //lepton efficiencies
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

  TH1 *h=mon.addHistogram( new TH1F ("eventflowsinglelepton", ";;Events", 5,0,5) );
  h->GetXaxis()->SetBinLabel(1,"1 lep, #geq 3 jets");
  h->GetXaxis()->SetBinLabel(2,"E_{T}^{miss} #geq 40 GeV");
  h->GetXaxis()->SetBinLabel(3,"#geq 1 btag");
  h->GetXaxis()->SetBinLabel(4,"1 #tau_{h} "); 
  h->GetXaxis()->SetBinLabel(5,"OS");

  h=mon.addHistogram( new TH1F ("eventflowdileptons", ";;Events", 6,0,6) );
  h->GetXaxis()->SetBinLabel(1,"2 leptons");
  h->GetXaxis()->SetBinLabel(2,"M>12 GeV #wedge |M-M_{Z}|>15 GeV");
  h->GetXaxis()->SetBinLabel(3,"#geq 2 jets");
  h->GetXaxis()->SetBinLabel(4,"op. sign"); 
  h->GetXaxis()->SetBinLabel(5,"OS");
  h->GetXaxis()->SetBinLabel(6,"#geq 2 b-tags");

  TString var("");
  h = mon.addHistogram( new TH1D("finalevtflowdileptons"+var,";Category;Events",6,0,6) ); 
  h->GetXaxis()->SetBinLabel(1,"=0 btags");
  h->GetXaxis()->SetBinLabel(2,"=1 btags");
  h->GetXaxis()->SetBinLabel(3,"=2 btags");
  h->GetXaxis()->SetBinLabel(4,"=3 btags");
  h->GetXaxis()->SetBinLabel(5,"=4 btags");
  h->GetXaxis()->SetBinLabel(6,"#geq5 btags");
  
  h = mon.addHistogram( new TH1D("finalevtflow2btagsdileptons"+var,";Category;Events",4,2,6) ); 
  h->GetXaxis()->SetBinLabel(1,"=2 btags");
  h->GetXaxis()->SetBinLabel(2,"=3 btags");
  h->GetXaxis()->SetBinLabel(3,"=4 btags");
  h->GetXaxis()->SetBinLabel(4,"#geq5 btags");
  
  //  h->GetXaxis()->SetBinLabel(1,"#geq 2 leptons");
  //  h->GetXaxis()->SetBinLabel(2,"|M-M_{Z}|<30");
  //  h->GetXaxis()->SetBinLabel(3,"Nlep+N#tau #geq 4");
  //  h->GetXaxis()->SetBinLabel(4,"HiggsCandidate"); 
  //  h->GetXaxis()->SetBinLabel(5,"lepton Veto");
  //  h->GetXaxis()->SetBinLabel(6,"Btag Veto");
  //  h->GetXaxis()->SetBinLabel(7," ");
  //  h->GetXaxis()->SetBinLabel(8,"mm_em");
  //  h->GetXaxis()->SetBinLabel(9,"mm_et");
  //  h->GetXaxis()->SetBinLabel(10,"mm_mt");
  //  h->GetXaxis()->SetBinLabel(11,"mm_tt");
  //  h->GetXaxis()->SetBinLabel(12,"ee_em");
  //  h->GetXaxis()->SetBinLabel(13,"ee_et");
  //  h->GetXaxis()->SetBinLabel(14,"ee_mt");
  //  h->GetXaxis()->SetBinLabel(15,"ee_tt");



  std::vector<TString> controlCats;
  controlCats.push_back("eq2leptons");
  controlCats.push_back("eq1jets");   
  controlCats.push_back("eq2jets");   
  controlCats.push_back(""); // FIXME: add DY reweighting
  controlCats.push_back("geq2btags");
  controlCats.push_back("lowmet");
  controlCats.push_back("eq1jetslowmet");
  controlCats.push_back("z");
  controlCats.push_back("zeq1jets");
  controlCats.push_back("zlowmet");
  controlCats.push_back("zeq1jetslowmet");
  
  for(size_t k=0; k<controlCats.size(); ++k){
    mon.addHistogram( new TH1F(controlCats[k]+"emva", "; e-id MVA; Electrons", 50, 0.95,1.0) );
    mon.addHistogram( new TH1F(controlCats[k]+"mll",";Dilepton invariant mass [GeV];Events",50,0,250) );
    mon.addHistogram( new TH1F(controlCats[k]+"ptll",";Dilepton transverse momentum [GeV];Events",50,0,250) );
    mon.addHistogram( new TH1F(controlCats[k]+"pte",";Electron transverse momentum [GeV];Events",50,0,500) );
    mon.addHistogram( new TH1F(controlCats[k]+"ptmu",";Muon transverse momentum [GeV];Events",50,0,500) );
    mon.addHistogram( new TH1F(controlCats[k]+"ptlep",";Lepton transverse momentum [GeV];Events",50,0,500) ); 
    mon.addHistogram( new TH1F(controlCats[k]+"sumpt",";Sum of lepton transverse momenta [GeV];Events",50,0,500) );
    mon.addHistogram( new TH1F(controlCats[k]+"ptmin",";Minimum lepton transverse momentum [GeV];Events",50,0,500) );
    mon.addHistogram( new TH1F(controlCats[k]+"dilarccosine",";#theta(l,l') [rad];Events",50,0,3.2) );
    mon.addHistogram( new TH1F(controlCats[k]+"mtsum",";M_{T}(l^{1},E_{T}^{miss})+M_{T}(l^{2},E_{T}^{miss}) [GeV];Events",100,0,1000) );
    mon.addHistogram( new TH1F(controlCats[k]+"met",";Missing transverse energy [GeV];Events",50,0,500) );
    mon.addHistogram( new TH1F(controlCats[k]+"metnotoppt",";Missing transverse energy [GeV];Events",50,0,500) );
    mon.addHistogram( new TH1F(controlCats[k]+"ht",";H_{T} [GeV];Events",50,0,1000) );
    mon.addHistogram( new TH1F(controlCats[k]+"htb",";H_{T} (bjets) [GeV];Events",50,0,1000) );
    mon.addHistogram( new TH1F(controlCats[k]+"htnol","; H_[T] (no leptons) [GeV];Events",50,0,1000) );
    mon.addHistogram( new TH1F(controlCats[k]+"htbnol","; H_[T] (bjets, no leptons) [GeV];Events",50,0,1000) );
  }
  
  mon.addHistogram( new TH1F("pthat",";#hat{p}_{T} [GeV];Events",50,0,1000) );
  mon.addHistogram( new TH1F("nup",";NUP;Events",10,0,10) );
  mon.addHistogram( new TH1F("nupfilt",";NUP;Events",10,0,10) );

  //pileup control
  mon.addHistogram( new TH1F( "nvtx",";Vertices;Events",50,0,50) ); 
  mon.addHistogram( new TH1F( "nvtxraw",";Vertices;Events",50,0,50) ); 
  mon.addHistogram( new TH1F( "rho",";#rho;Events",50,0,25) ); 
  mon.addHistogram( new TH1F( "rho25",";#rho(#eta<2.5);Events",50,0,25) ); 

  //lepton control
  mon.addHistogram( new TH1F( "leadpt", ";p_{T}^{l};Events", 50,0,500) );
  mon.addHistogram( new TH1F( "leadeta", ";#eta^{l};Events", 50,-2.6,2.6) );
  mon.addHistogram( new TH1F( "trailerpt", ";p_{T}^{l};Events", 50,0,500) );
  mon.addHistogram( new TH1F( "trailereta", ";#eta^{l};Events", 50,-2.6,2.6) );
  mon.addHistogram( new TH1F( "leppt", ";p_{T}^{l};Events", 50,0,500) );

  //jets control
  mon.addHistogram( new TH1F( "njets",      ";njets;Events", 6,0,6) );
  mon.addHistogram( new TH1F( "jetpt",  ";p_{T}^{jet};Events", 50,0,500) );
  mon.addHistogram( new TH1F( "jeteta", ";#eta^{jet};Events", 50,-2.6, 2.6) );

  //tau control
  mon.addHistogram( new TH1F( "ntaus",      ";ntaus;Events", 6,0,6) );
  mon.addHistogram( new TH1F( "tauleadpt",  ";p_{T}^{#tau};Events", 50,0,500) );
  mon.addHistogram( new TH1F( "tauleadeta", ";#eta^{#tau};Events", 50,-2.6,2.6) );
  mon.addHistogram( new TH1F( "taupt",  ";p_{T}^{#tau};Events", 50,0,500) );
  mon.addHistogram( new TH1F( "taueta",  ";#eta^{#tau};Events", 50,-2.6,2.6) );
  mon.addHistogram( new TH1F( "taucharge",  ";p_{T}^{#tau};Events", 5,-2,2) );
  mon.addHistogram( new TH1F( "taudz",      ";dz^{#tau};Events", 50,0,10) );
  mon.addHistogram( new TH1F( "tauvz",      ";vz^{#tau};Events", 50,0,10) );
  mon.addHistogram( new TH1F( "tauemfraction", ";emf^{#tau};Events", 50, 0., 5.) );
  mon.addHistogram( new TH1F( "taudizeta"    , ";dZ^{#tau};Events", 50, 0., 10.) );
  
  mon.addHistogram( new TH1F( "ntausos",      ";ntaus;Events", 6,0,6) );
  mon.addHistogram( new TH1F( "tauleadptos",  ";p_{T}^{#tau};Events", 50,0,500) );
  mon.addHistogram( new TH1F( "tauleadetaos", ";#eta^{#tau};Events", 50,-2.6,2.6) );
  mon.addHistogram( new TH1F( "tauptos",  ";p_{T}^{#tau};Events", 50,0,500) );
  mon.addHistogram( new TH1F( "tauetaos",  ";#eta^{#tau};Events", 50,-2.6,2.6) );
  mon.addHistogram( new TH1F( "tauchargeos",  ";p_{T}^{#tau};Events", 5,-2,2) );
  mon.addHistogram( new TH1F( "taudzos",      ";dz^{#tau};Events", 50,0,10) );
  mon.addHistogram( new TH1F( "tauvzos",      ";vz^{#tau};Events", 50,0,10) );
  mon.addHistogram( new TH1F( "tauemfractionos", ";emf^{#tau};Events", 50, 0., 5.) );
  mon.addHistogram( new TH1F( "taudizetaos"    , ";dZ^{#tau};Events", 50, 0., 10.) );
  

  //bjets control
  mon.addHistogram( new TH1F( "nbjets",      ";ntaus;Events", 6,0,6) );
  mon.addHistogram( new TH1F( "bjetpt",  ";p_{T}^{bjet};Events", 50,0,500) );
  mon.addHistogram( new TH1F( "bjetcsv", ";csv discriminator;Events", 50,0, 1) );

  //boson control
  mon.addHistogram( new TH1F( "qt",      ";p_{T}^{#gamma} [GeV];Events",500,0,1500));
  mon.addHistogram( new TH1F( "zpt",     ";p_{T}^{ll};Events", 50,0,500) );
  mon.addHistogram( new TH1F( "zptNM1",  ";p_{T}^{ll};Events", 50,0,500) );
  mon.addHistogram( new TH1F( "zeta",    ";#eta^{ll};Events", 50,-10,10) );
  mon.addHistogram( new TH1F( "zetaNM1", ";#eta^{ll};Events", 50,-10,10) );
  mon.addHistogram( new TH1F( "zy",      ";y^{ll};Events", 50,-6,6) );
  mon.addHistogram( new TH1F( "rapidity",";y^{ll};Events", 50,0,2) );
  mon.addHistogram( new TH1F( "zyNM1",   ";y^{ll};Events", 50,-6,6) );
  mon.addHistogram( new TH1F( "zmass",   ";M^{ll};Events", 100,40,250) );

  //higgs control
  mon.addHistogram( new TH1F( "higgspt",      ";p_{T}^{#higgs} [GeV];Events",500,0,1500));
  mon.addHistogram( new TH1F( "higgsmass",    ";M^{#higgs} [GeV];Events",100,0,500));
  mon.addHistogram( new TH1F( "higgsmasssvfit",    ";M^{#higgs} [GeV];Events",100,0,500));
  mon.addHistogram( new TH1F( "higgsmet",    ";MET [GeV];Events",100,0,500));
  mon.addHistogram( new TH1F( "higgsnjets",   ";NJets;Events",10,0,10));

  //statistical analysis
  std::vector<double> optim_Cuts_sumpt; 
  std::vector<double> optim_Cuts_lepIso; 
  for(double sumpt=10;sumpt<=100;sumpt+=5)
    {
      for(double lepIso=0.4;lepIso>=0.05;lepIso-=0.05)
	{
	  optim_Cuts_sumpt.push_back(sumpt);
	  optim_Cuts_lepIso.push_back(lepIso);
	} 
    }
  TH2F* Hoptim_cuts  =(TH2F*)mon.addHistogram(new TProfile2D("optim_cut",      ";cut index;variable",       optim_Cuts_sumpt.size(),0,optim_Cuts_sumpt.size(), 2, 0, 2)) ;
  Hoptim_cuts->GetYaxis()->SetBinLabel(1, "jpt1>");
  Hoptim_cuts->GetYaxis()->SetBinLabel(2, "jpt2>");
  for(unsigned int index=0;index<optim_Cuts_sumpt.size();index++){
    Hoptim_cuts->Fill(index,0.0,optim_Cuts_sumpt[index]); 
    Hoptim_cuts->Fill(index,1.0,optim_Cuts_lepIso[index]); 
  }

  TH1F* Hoptim_systs     =  (TH1F*) mon.addHistogram( new TH1F ("optim_systs"    , ";syst;", nvarsToInclude,0,nvarsToInclude) ) ;
  for(size_t ivar=0; ivar<nvarsToInclude; ivar++)
  {
    Hoptim_systs->GetXaxis()->SetBinLabel(ivar+1, varNames[ivar]);
    mon.addHistogram( new TH2F (TString("svfit_shapes")+varNames[ivar],";cut index;|#Delta #eta|;Events",optim_Cuts_sumpt.size(),0,optim_Cuts_sumpt.size(),25,0,250) );
  }


  //tau fakeRate
  mon.addHistogram( new TH2F( "taufakerate",    ";CutIndex;p_{T}^{#tau} [GeV]",optim_Cuts_sumpt.size(),0,optim_Cuts_sumpt.size(),20,0,100));

  
  //##############################################
  //######## GET READY FOR THE EVENT LOOP ########
  //##############################################
  
  //open the file and get events tree
  //  DataEventSummaryHandler evSummary;
  //  TFile *file = TFile::Open(url);
  //  printf("Looping on %s\n",url.Data());
  //  if(file==0) return -1;
  //  if(file->IsZombie()) return -1;
  //  if( !evSummary.attach( (TTree *) file->Get(baseDir+"/data") ) ) { file->Close();  return -1; }
  //  const Int_t totalEntries= evSummary.getEntries();
  
  fwlite::ChainEvent ev(urls);
  const Int_t totalEntries= ev.size();
 
  //MC normalization (to 1/pb)
  float nInitEvent(1.0);
  if(isMC){
     nInitEvent = (float)utils::getMergeableCounterValue(urls, "startCounter");
  }
  double xsecWeight(xsec/nInitEvent);
  if(!isMC) xsecWeight=1.0;

  bool isTTbarMC(isMC && (url.Contains("TTJets") || url.Contains("_TT_") ));
  // Top Pt weighter
  TopPtWeighter* topPtWgt=0;
  if(isTTbarMC){
    TString shapesDir("");
    if(weightsFile.size()) shapesDir=wFile;
    topPtWgt = new TopPtWeighter(outFileUrl, outUrl, shapesDir, ev );
  }

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
  edm::LumiReWeighting*          LumiWeights = NULL;
  edm::LumiReWeighting* singleLepLumiWeights = NULL;
  utils::cmssw::PuShifter_t PuShifters, singleLepPuShifters;
  double 
             PUNorm[] = {1,1,1},
    singleLepPUNorm[] = {1,1,1};

  if(isMC){
    std::vector<double> dataPileupDistributionDouble = runProcess.getParameter< std::vector<double> >("datapileup");
    std::vector<double> singleLepDataPileupDistributionDouble = runProcess.getParameter< std::vector<double> >("datapileupSingleLep"); 
    std::vector<float> dataPileupDistribution; for(unsigned int i=0;i<dataPileupDistributionDouble.size();i++){dataPileupDistribution.push_back(dataPileupDistributionDouble[i]);}
    std::vector<float> singleLepDataPileupDistribution; for(unsigned int i=0;i<singleLepDataPileupDistributionDouble.size();i++){singleLepDataPileupDistribution.push_back(singleLepDataPileupDistributionDouble[i]);}
    std::vector<float> mcPileupDistribution;
    std::vector<float> singleLepMcPileupDistribution;
    utils::getMCPileupDistribution(ev,dataPileupDistribution.size(), mcPileupDistribution);
    utils::getMCPileupDistribution(ev,singleLepDataPileupDistribution.size(), singleLepMcPileupDistribution);
    while(mcPileupDistribution.size()<dataPileupDistribution.size())  mcPileupDistribution.push_back(0.0);
    while(mcPileupDistribution.size()>dataPileupDistribution.size())dataPileupDistribution.push_back(0.0);
    while(singleLepMcPileupDistribution.size()<singleLepDataPileupDistribution.size())  singleLepMcPileupDistribution.push_back(0.0);
    while(singleLepMcPileupDistribution.size()>singleLepDataPileupDistribution.size())singleLepDataPileupDistribution.push_back(0.0);
    
    LumiWeights= new edm::LumiReWeighting(mcPileupDistribution,dataPileupDistribution);
    PuShifters=utils::cmssw::getPUshifters(dataPileupDistribution,0.05);
    utils::getPileupNormalization(mcPileupDistribution, PUNorm, LumiWeights, PuShifters);

    singleLepLumiWeights= new edm::LumiReWeighting(singleLepMcPileupDistribution,singleLepDataPileupDistribution);
    singleLepPuShifters=utils::cmssw::getPUshifters(singleLepDataPileupDistribution,0.05);
    utils::getPileupNormalization(singleLepMcPileupDistribution, singleLepPUNorm, singleLepLumiWeights, singleLepPuShifters);

  }
  
  gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE

  //##############################################
  //########           EVENT LOOP         ########
  //##############################################
  //loop on all the events
  printf("Progressing Bar     :0%%       20%%       40%%       60%%       80%%       100%%\n");
  printf("Scanning the ntuple :");
  DuplicatesChecker duplicatesChecker;
  int nDuplicates(0);
  int step(totalEntries/50); 
  for( int iev=0; iev<totalEntries; iev++ ) 
    {
      if(iev%step==0){printf(".");fflush(stdout);}
      if(!isMC && jacknife>0 && jacks>0 && iev%jacks==jacknife) continue;
      
      //##############################################   EVENT LOOP STARTS   ##############################################
      //load the event content from tree
      //      evSummary.getEntry(iev);
      ev.to(iev);
      //      DataEventSummary &ev = evSummary.getEvent();
      if(!isMC && duplicatesChecker.isDuplicate( ev.eventAuxiliary().run(),ev.eventAuxiliary().luminosityBlock(),ev.eventAuxiliary().event()) ) { nDuplicates++; continue; }

      //get the collection of generated Particles
      fwlite::Handle< llvvGenEvent > genEventHandle;
      genEventHandle.getByLabel(ev, "llvvObjectProducersUsed");
      if(!genEventHandle.isValid()){printf("llvvGenEvent Object NotFound\n");continue;}
      llvvGenEvent genEv = *genEventHandle;


      if(isV0JetsMC){
	mon.fillHisto("nup","",genEv.nup,1);
	if(genEv.nup>5) continue;
	mon.fillHisto("nupfilt","",genEv.nup,1);
      }
      
      fwlite::Handle< llvvGenParticleCollection > genPartCollHandle;
      genPartCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
      if(!genPartCollHandle.isValid()){printf("llvvGenParticleCollection Object NotFound\n");continue;}
      llvvGenParticleCollection gen = *genPartCollHandle;

      fwlite::Handle< llvvLeptonCollection > leptonCollHandle;
      leptonCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
      if(!leptonCollHandle.isValid()){printf("llvvLeptonCollection Object NotFound\n");continue;}
      llvvLeptonCollection leptons = *leptonCollHandle;

      fwlite::Handle< llvvElectronInfoCollection > electronInfoCollHandle;
      electronInfoCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
      if(!electronInfoCollHandle.isValid()){printf("llvvElectronInfoCollection Object NotFound\n");continue;}

      fwlite::Handle< llvvMuonInfoCollection > muonInfoCollHandle;
      muonInfoCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
      if(!muonInfoCollHandle.isValid()){printf("llvvMuonInfoCollection Object NotFound\n");continue;}

      fwlite::Handle< llvvTauCollection > tauCollHandle;
      tauCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
      if(!tauCollHandle.isValid()){printf("llvvLeptonCollection Object NotFound\n");continue;}
      llvvTauCollection taus = *tauCollHandle;

      fwlite::Handle< llvvJetCollection > jetCollHandle;
      jetCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
      if(!jetCollHandle.isValid()){printf("llvvJetCollection Object NotFound\n");continue;}
      llvvJetExtCollection jets;
      for(unsigned int i=0;i<jetCollHandle->size();i++){jets.push_back(llvvJetExt((*jetCollHandle)[i]));}

      fwlite::Handle< llvvMet > metHandle;
      metHandle.getByLabel(ev, "llvvObjectProducersUsed", "pfMETPFlow"); 
      //metHandle.getByLabel(ev, "llvvObjectProducersUsed", "pfType1CorrectedMet"); 
      if(!metHandle.isValid()){printf("llvvMet Object NotFound\n");continue;}
      llvvMet met = *metHandle;

      fwlite::Handle< std::vector<bool> > triggerBitsHandle;
      triggerBitsHandle.getByLabel(ev, "llvvObjectProducersUsed", "triggerBits");
      if(!triggerBitsHandle.isValid()){printf("triggerBits Object NotFound\n");continue;}
      std::vector<bool> triggerBits = *triggerBitsHandle;

      fwlite::Handle< std::vector<int> > triggerPrescalesHandle;
      triggerPrescalesHandle.getByLabel(ev, "llvvObjectProducersUsed", "triggerPrescales");
      if(!triggerPrescalesHandle.isValid()){printf("triggerPrescales Object NotFound\n");continue;}
      std::vector<int> triggerPrescales = *triggerPrescalesHandle;

      fwlite::Handle< double > rhoHandle;
      rhoHandle.getByLabel(ev, "kt6PFJets", "rho");
      if(!rhoHandle.isValid()){printf("rho Object NotFound\n");continue;}
      double rho(*rhoHandle);

      fwlite::Handle< double > rho25Handle;
      rho25Handle.getByLabel(ev, "kt6PFJetsCentral", "rho");
      if(!rho25Handle.isValid()){printf("rho25 Object NotFound\n");continue;}
      double rho25(*rho25Handle);

      fwlite::Handle< int > nvtxHandle;
      nvtxHandle.getByLabel(ev, "llvvObjectProducersUsed", "nvtx");
      if(!nvtxHandle.isValid()){printf("nvtx Object NotFound\n");continue;}
      int nvtx(*nvtxHandle);

      //require compatibilitity of the event with the PD
      bool eeTrigger          ( triggerBits[0]                   );
      bool emuTrigger         ( triggerBits[4] || triggerBits[5] );
      bool muTrigger          ( triggerBits[6]                   );
      //      bool eTrigger           ( triggerBits[12]                  ); // FIXME: must process singleElectron
      bool mumuTrigger        ( triggerBits[2] || triggerBits[3] );// || muTrigger;
      if(filterOnlyEE)       {                  emuTrigger=false; mumuTrigger=false; muTrigger=false; /*eTrigger=false;*/}
      if(filterOnlyEMU)      { eeTrigger=false;                   mumuTrigger=false; muTrigger=false; /*eTrigger=false;*/}
      if(filterOnlyMUMU)     { eeTrigger=false; emuTrigger=false;                    muTrigger=false; /*eTrigger=false;*/}
      if(filterOnlySINGLEMU) { eeTrigger=false; emuTrigger=false; mumuTrigger=false;                  /*eTrigger=false;*/}
      //      if(filterOnlySINGLEELE){ eeTrigger=false; emuTrigger=false; mumuTrigger=false; muTrigger=false;                }
      
      //      if(isSingleMuPD)   { eeTrigger=false; if( mumuTrigger || !muTrigger ) mumuTrigger= false;  }

      // PU weighting now follows channel selection, to take into account the different pileup distribution between singleMu and dileptons
      //pileup weight
      float weight(1.0);
      double TotalWeight_plus(1.0);
      double TotalWeight_minus(1.0);
      float puWeight(1.0);
      

      //
      //
      // BELOW FOLLOWS THE ANALYSIS OF THE MAIN SELECTION WITH N-1 PLOTS
      //
      //
      
      //
      // LEPTON ANALYSIS
      //
      llvvLeptonCollection selLeptons;
      llvvLeptonCollection selSingleLepLeptons;
      for(size_t ilep=0; ilep<leptons.size(); ilep++)
	{

	  bool 
	    passKin(true),         passId(true),         passIso(true),
	    passSingleLepKin(true),passSingleLepId(true),passSingleLepIso(true);

	  int lid(leptons[ilep].id);
	  //if(isSingleMuPD && abs(lid)==11) continue; // FIXME: raw veto
	  
	  //apply muon corrections
	  if(abs(lid)==13 && muCor){
            TLorentzVector p4(leptons[ilep].px(), leptons[ilep].py(), leptons[ilep].pz(), leptons[ilep].energy());
	    muCor->applyPtCorrection(p4 , lid<0 ? -1 :1 );
	    if(isMC) muCor->applyPtSmearing(p4, lid<0 ? -1 : 1, false);
            leptons[ilep].SetPxPyPzE(p4.Px(), p4.Py(), p4.Pz(), p4.Energy());
	  }
	  // no need for charge sign anymore
	  lid=abs(lid);
	  
	  TString lepStr( lid==13 ? "mu" : "e");
	  
	  //kinematics
	  //	  float leta = lid==11 ? leptons[ilep].getVal("sceta") : leptons[ilep].eta();
          float leta( lid==11 ? leptons[ilep].electronInfoRef->sceta : leptons[ilep].eta() );


	  // Dilepton kin
	  if(leptons[ilep].pt()< 20. )                   passKin=false;
	  if(abs(leta)> (lid==11 ? 2.5 : 2.4) )               passKin=false;
	  if(lid==11 && (abs(leta)>1.4442 && abs(leta)<1.5660))    passKin=false; // Crack veto

	  // SingleLepton kin
	  if(leptons[ilep].pt()< (lid==11 ? 35 :30 ) )  passSingleLepKin=false;
	  if(abs(leta)> (lid==11 ? 2.5 : 2.1) )              passSingleLepKin=false;
	  if(lid==11 && (abs(leta)>1.4442 && abs(leta)<1.5660))   passSingleLepKin=false; // Crack veto

	  //id
	  Int_t idbits(leptons[ilep].idbits);
	  if(lid==11){
	    if(leptons[ilep].electronInfoRef->isConv)              passId=false;
	    bool isLoose(((idbits >> 4) & 0x1));
	    if(!isLoose)                                           passId=false;

	    passSingleLepId=passId;
 	  }
	  else{
	    bool isLoose    = ((idbits >> 8) & 0x1);
	    bool isTight    (((idbits >> 10) & 0x1));
	    if(!isLoose)                                   passId=false;
	    if(!isTight)          	          passSingleLepId=false;
	  }
	  
	  //isolation
          float relIso( utils::cmssw::relIso(leptons[ilep], rho) );
          if( (lid==11 && relIso>0.15) || (lid==13 && relIso>0.20) ) passIso=false;
          if( (lid==11 && relIso>0.15)  || (lid==13 && relIso>0.12) ) passSingleLepIso=false;
	  
	  if(passId          && passIso          && passKin         ) selLeptons.push_back(leptons[ilep]);
	  if(passSingleLepId && passSingleLepIso && passSingleLepKin) selSingleLepLeptons.push_back(leptons[ilep]);
	}
      std::sort(selLeptons.begin()         , selLeptons.end()         , sort_llvvObjectByPt);
      std::sort(selSingleLepLeptons.begin(), selSingleLepLeptons.end(), sort_llvvObjectByPt);
      
      //at this point check if it's worth continuing
      //      if(selLeptons.size()<1) continue;


      //
      // SINGLE LEPTON ANALYSIS
      //
      
      // Leading lepton is our lepton of choice
      llvvLepton leadingSingleLep;
      int singleLeptonId(-1);
      if(selSingleLepLeptons.size()>0){
	leadingSingleLep=selSingleLepLeptons[0];
	singleLeptonId=leadingSingleLep.id;
      }
      
      // Veto additional leptons in the event
      size_t
	nVetoE(0),
	nVetoMu(0);
      
      if(selSingleLepLeptons.size()>0){
	for(size_t ilep=0; ilep<leptons.size(); ilep++)
	  {
	    if(leptons[ilep] == leadingSingleLep) continue; // Don't veto on the main lepton
	    int lid(abs(leptons[ilep].id));
	    //    	    if(isSingleMuPD && abs(lid)==11){ nVetoE++; continue;} // FIXME: temp Raw veto

	    bool passKin(true),passId(true),passIso(true);	    
	    // Muon correction already applied in the first loop. LoL.
	    
	    TString lepStr( lid==13 ? "mu" : "e");
	    
	    
	    //kinematics
	    //	  float leta = lid==11 ? leptons[ilep].getVal("sceta") : leptons[ilep].eta();
	    float leta( lid==11 ? leptons[ilep].electronInfoRef->sceta : leptons[ilep].eta() );
	    if(leptons[ilep].pt()< (lid==11 ? 20 : 10 ) )                   passKin=false; // Single lepton has a higher cut for muons
	    if(abs(leta)> 2.5 )                                             passKin=false; // Single lepton has 2.1 for eta. Double lepton 2.1 
	    if(lid==11 && (abs(leta)>1.4442 && abs(leta)<1.5660))           passKin=false; // Crack veto
	    
	    //id
	    Int_t idbits(leptons[ilep].idbits);
	    if(lid==11){
	      if(leptons[ilep].electronInfoRef->isConv)              passId=false;
	      bool isLoose ( ((idbits >> 4) & 0x1));
	      if(!isLoose)                                   passId=false;
	    }
	    else{
	      bool isLoose    ( ((idbits >> 8) & 0x1) ); // Veto is loose
	      //bool isTight    = ((idbits >> 10) & 0x1);
	      if(!isLoose)                                   passId=false;
	    }
	    
	    //isolation
	    float relIso( utils::cmssw::relIso(leptons[ilep], rho) );
	    if( (lid==11 && relIso>0.15) || (lid!=11 && relIso>0.2) ) passIso=false; // SingleLepton values
	    
	    if(!passId || !passIso || !passKin) continue;
	  
	    
	    if(lid==11) nVetoE++;
	    if(lid==13) nVetoMu++; 
	  }
      }
      
      //
      // DILEPTON ANALYSIS
      //

      int 
	dilLep1(-1), dilLep2(-1),
	dilId(-1);
      
      
      if(selLeptons.size()>=2){
	dilLep1=0;
	dilLep2=1;

	dilId = selLeptons[dilLep1].id * selLeptons[dilLep2].id;
	
      }

      bool isSameFlavour(abs(dilId)==121 || abs(dilId)==169);


      //check the channel
      //prepare the tag's vectors for histo filling
      std::vector<TString> chTags;
      chTags.push_back("ll"); // FIXME: legacy to remove.
      
      if( abs(dilId)==121 && eeTrigger  ) chTags.push_back("ee");
      else if( abs(singleLeptonId) == 13 && muTrigger /*&& selSingleLepLeptons.size()>0 */&& nVetoE==0 && nVetoMu==0 ) chTags.push_back("singlemu"); // selSingleLepLeptons.size() implicitly checked by abs(singleLeptonId)==13
      else if( abs(dilId)==143 && emuTrigger ) chTags.push_back("emu");
      else if( abs(dilId)==169 && mumuTrigger) chTags.push_back("mumu"); 
      
      
      if(chTags.size()<2) continue;

      if(debug && chTags.size() > 2 ) cout << "[DEBUG] ERROR. The same event has been chosen for inclusive final state and MORE THAN ONE exclusive final state" << endl;

      
      //apply data/mc correction factors
      
      // Dilepton full analysis
      // ----------------------
      if(chTags[1] == "ee" || chTags[1] == "emu" || chTags[1] == "mumu" /*|| chTags[0] == "ll" */ ){

	if(selLeptons.size()<2) continue; // 2 leptons
	
	if(isMC){
	  puWeight          = LumiWeights->weight(genEv.ngenITpu) * PUNorm[0];
	  weight            = xsecWeight*puWeight;
	  TotalWeight_plus  = PuShifters[utils::cmssw::PUUP  ]->Eval(genEv.ngenITpu) * (PUNorm[2]/PUNorm[0]);
	  TotalWeight_minus = PuShifters[utils::cmssw::PUDOWN]->Eval(genEv.ngenITpu) * (PUNorm[1]/PUNorm[0]);
	}
	
	weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[dilLep1].pt(), selLeptons[dilLep1].eta(), abs(selLeptons[dilLep1].id),  abs(selLeptons[dilLep1].id) ==11 ? "loose" : "loose" ).first : 1.0;
	weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[dilLep2].pt(), selLeptons[dilLep2].eta(), abs(selLeptons[dilLep2].id),  abs(selLeptons[dilLep2].id) ==11 ? "loose" : "loose" ).first : 1.0;
	
	//generator level
	LorentzVector genll;
	for(size_t ig=0; ig<gen.size(); ig++){
	  int pid(abs(gen[ig].id));
	  if(pid!=11 && pid!=13) continue;
	  genll += gen[ig];
	}
	
	//
	//JET/MET ANALYSIS
	//
	llvvJetExtCollection selJets, /*selJetsNoId,*/ selBJets;
	int njets(0), nbjets(0);
	
	for(size_t ijet=0; ijet<jets.size(); ijet++) 
	  {
	    //correct jet
	    float toRawSF(jets[ijet].torawsf);
	    LorentzVector rawJet(jets[ijet]*toRawSF);
	    jesCor->setJetEta(rawJet.eta());
	    jesCor->setJetPt(rawJet.pt());
	    jesCor->setJetA(jets[ijet].area);
	    jesCor->setRho(rho);
	    float newJECSF(jesCor->getCorrection());
	    jets[ijet].SetPxPyPzE(rawJet.px(),rawJet.py(),rawJet.pz(),rawJet.energy());
	    jets[ijet] *= newJECSF;
	    jets[ijet].torawsf = 1./newJECSF;
	    if(jets[ijet].pt()<15 || fabs(jets[ijet].eta())>4.7 ) continue;

	    //cross-clean with selected leptons and photons
	    double minDRlj(9999.),minDRlg(9999.);
	    for(size_t ilep=0; ilep<selLeptons.size(); ilep++)
	      minDRlj = TMath::Min( minDRlj, deltaR(jets[ijet],selLeptons[ilep]) );
	    if(minDRlj<0.4 || minDRlg<0.4) continue;
	    
	    
	    //jet id
	    // float pumva=jets[ijet].puMVA;
	    Int_t idbits(jets[ijet].idbits);
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
	  
	    double jetpt=jets[ijet].pt();
	    //selJetsNoId.push_back(jets[ijet]);
	    
	    if(!passPFloose || !passLooseSimplePuId || jets[ijet].pt()<minJetPtToApply || abs(jets[ijet].eta())>2.5) continue;
	    
	    //if(passPFloose && passLooseSimplePuId){
	    //  selJets.push_back(jets[ijet]);
	    //  if(jets[ijet].pt()>minJetPtToApply) njets++;
	    //}
	    selJets.push_back(jets[ijet]);
	    njets++;

	    
	    //bjets
	    mon.fillHisto("bjetpt"    ,  chTags, jets[ijet].pt(),  weight);
	    mon.fillHisto("bjetcsv"   ,  chTags, jets[ijet].origcsv,  weight);
	
	    bool hasCSVV1L(jets[ijet].csv >0.405);
	    bool hasBtagCorr(hasCSVV1L);
	    if(isMC){
	      // set a unique seed
	      double bseed_sin_phi = sin(jets[ijet].phi()*1000000);
	      double bseed = abs(static_cast<int>(bseed_sin_phi*100000));
	      
	      // get jet flavour
	      int bflavid= jets[ijet].genflav;
	      
	      //Initialize
	      BTagSFUtil btsfutil( bseed );
	      
	      TString flavKey("udsg");
	      if(abs(bflavid)==4) flavKey="c";
	      if(abs(bflavid)==5) flavKey="b";
	      std::pair<TString,TString> btagKey("csvL",flavKey);
	      if(btagEffCorr.find(btagKey)!=btagEffCorr.end())
		{
		  TGraphErrors* mceffGr=btagEffCorr[btagKey].first;
		  TGraphErrors* sfGr=btagEffCorr[btagKey].second;
		  if(mceffGr && sfGr){
		    float eff=mceffGr->Eval(jetpt);
		    float sf=sfGr->Eval(jetpt);	
		    if(var == "btagup" || var == "btagdown" || var == "unbtagup" || var == "unbtagdown"){
		      // Apply uncertainty shift
		      float sferr(0.), delta(1000000.);
		      for(size_t ind=0; ind<btagBins.size(); ++ind){
			float tempDelta( fabs(jetpt - btagBins[ind]) );
			if(tempDelta<delta){
			  delta=tempDelta;
			  if      (var == "btagup"   || var == "unbtagup")   sferr=fabs(sfGr->GetErrorYhigh(ind));  // Ensure positive number
			  else if (var == "btagdown" || var == "unbtagdown") sferr=0-fabs(sfGr->GetErrorYlow(ind)); // Ensure negative number
			}
		      }
		      
		      //		      for(int ind=0; ind<sfGr->GetN(); ++ind){
		      //			double xv(0.), yv(0.);
		      //			sfGr->GetPoint(ind,xv,yv);
		      //			float tempDelta(jetpt - xv );
		      //			if(tempDelta<delta){
		      //			  delta=tempDelta; 
		      //			  if      (var == "btagup"   || var == "unbtagup")   sferr=fabs(sfGr->GetErrorYhigh(ind));  // Ensure positive number
		      //			  else if (var == "btagdown" || var == "unbtagdown") sferr=0-fabs(sfGr->GetErrorYlow(ind)); // Ensure negative number
		      //			}
		      //		      }
		      sf+=sferr;
		    }
		    btsfutil.modifyBTagsWithSF(hasBtagCorr, sf, eff);	    
		  }
		}
	    }
	    //	    if(!hasCSVV1L) continue;
	    if(!hasBtagCorr) continue;
	    
	    selBJets.push_back(jets[ijet]);  
	    nbjets++;
	  }
	std::sort(selJets.begin(), selJets.end(), sort_llvvObjectByPt);
	std::sort(selBJets.begin(), selBJets.end(), sort_llvvObjectByPt);
	mon.fillHisto("nbjets"    ,  chTags, nbjets,  weight);
	
	
	// NOW FOR THE CONTROL PLOTS
	//
	//start analysis
	//      for(size_t ich=0; ich<chTags.size(); ich++)
	//	{
	//	  std::vector<TString> tags(1,chTags[ich]);
	
	//      bool otherTriggersVeto( mumuTrigger && eeTrigger ); 
	
	LorentzVector ll(selLeptons[dilLep1]+selLeptons[dilLep2]);
	             
	float mll(ll.mass());
	float thetall(utils::cmssw::getArcCos<llvvLepton>(selLeptons[dilLep1],selLeptons[dilLep2]));
	float mtsum(utils::cmssw::getMT<llvvLepton,llvvMet>(selLeptons[dilLep1],met)+utils::cmssw::getMT<llvvLepton,llvvMet>(selLeptons[dilLep2],met));
	bool isZcand( isSameFlavour && abs(mll-91.2)<15);
	bool passDileptonSelection(mll>12. && !isZcand);
	bool passJetSelection(selJets.size()>=2);
	bool passMetSelection(met.pt()>40); // !isSameFlavour || met.pt()>40); (charged Higgs case expects neutrinos in the ee,mumu final states too)
	bool passBtagSelection(selBJets.size()>=2);
	bool isOS(dilId<0 && dilId!=-1);
	if(debug && dilId==-1) cout << "[DEBUG] ERROR: in dilepton analysis after >=2 leptons cut there is an event with dilId not set (==-1)" << endl;

	// control distributions
	std::vector<TString> ctrlCategs;
	ctrlCategs.push_back("eq2leptons");
	if(        passDileptonSelection && selJets.size()==1                                                ) ctrlCategs.push_back("eq1jets");   
	if(        passDileptonSelection && passJetSelection                                                 ) ctrlCategs.push_back("eq2jets");   
	if(isOS && passDileptonSelection && passJetSelection  && passMetSelection                            ) ctrlCategs.push_back(""); // FIXME: add DY reweighting
	if(isOS && passDileptonSelection && passJetSelection  && passMetSelection  && passBtagSelection      ) ctrlCategs.push_back("geq2btags");
	if(isOS && passDileptonSelection && passJetSelection  && met.pt()<30       /*&& passBtagSelection*/  ) ctrlCategs.push_back("lowmet");
	if(isOS && passDileptonSelection && selJets.size()==1 && met.pt()<30       /*&& passBtagSelection*/  ) ctrlCategs.push_back("eq1jetslowmet");
	if(isOS && isZcand          && passJetSelection  && passMetSelection  /*&& passBtagSelection*/       ) ctrlCategs.push_back("z");
	if(isOS && isZcand          && selJets.size()==1 && passMetSelection  /*&& passBtagSelection*/       ) ctrlCategs.push_back("zeq1jets");
	if(isOS && isZcand          && passJetSelection  && met.pt()<30       /*&& passBtagSelection*/       ) ctrlCategs.push_back("zlowmet");
	if(isOS && isZcand          && selJets.size()==1 && met.pt()<30       /*&& passBtagSelection*/       ) ctrlCategs.push_back("zeq1jetslowmet");
	for(size_t icat=0; icat<ctrlCategs.size(); icat++)
	  {
	    double ptmin((selLeptons[dilLep1].pt()>selLeptons[dilLep2].pt()) ? selLeptons[dilLep2].pt() : selLeptons[dilLep1].pt() );
	    double sumpt(selLeptons[dilLep1].pt()+selLeptons[dilLep2].pt());
	    double 
	      ht    (selLeptons[dilLep1].pt()+selLeptons[dilLep2].pt()),
	      htb   (selLeptons[dilLep1].pt()+selLeptons[dilLep2].pt()),
	      htnol (0),
	      htbnol(0);
	    
	    std::vector<int> dileptons;
	    dileptons.push_back(dilLep1);
	    dileptons.push_back(dilLep2);
	     
	    for(size_t ilep=0; ilep<2; ++ilep)
	      {
		if(abs(selLeptons[dileptons[ilep]].id)==11){
		  //		  mon.fillHisto(ctrlCategs[icat]+"emva", chTags, selLeptons[dileptons[ilep]].electronInfoRef->mvatrigv0, weight); // "BadRefCore RefCore: Request to resolve a null or invalid reference to a product of type 'std::vector<llvvElectronInfo>' has been detected. Please modify the calling code to test validity before dereferencing."
		  mon.fillHisto(ctrlCategs[icat]+"pte" ,  chTags, selLeptons[dileptons[ilep]].pt(),        weight);
		  mon.fillHisto(ctrlCategs[icat]+"ptlep",chTags, selLeptons[dileptons[ilep]].pt(),        weight);
		}
		else if(abs(selLeptons[dileptons[ilep]].id)==13){
		  mon.fillHisto(ctrlCategs[icat]+"ptmu",  chTags, selLeptons[dileptons[ilep]].pt(),        weight);
		  mon.fillHisto(ctrlCategs[icat]+"ptlep", chTags, selLeptons[dileptons[ilep]].pt(),        weight);
		}
	      }
	    mon.fillHisto(ctrlCategs[icat]+"sumpt", chTags, sumpt, weight);
	    
	    for(size_t ijet=0; ijet<selJets.size(); ++ijet) // FIXME: am I sure that for HT I want to use only jets with pt>30, eta<2.5?
	      {
		ht+=selJets[ijet].pt();
		htnol+=selJets[ijet].pt();
		if(ijet<selBJets.size())
		  {
		    htb+=selBJets[ijet].pt();
		    htbnol+=selBJets[ijet].pt();
		  }
	      }
	    ht+=met.pt();
	    htb+=met.pt();
	    htnol+=met.pt();
	    htbnol+=met.pt();
	    
	    mon.fillHisto(ctrlCategs[icat]+"ptmin",        chTags, ptmin,           weight);
	    mon.fillHisto(ctrlCategs[icat]+"mll",          chTags, mll,             weight);
	    mon.fillHisto(ctrlCategs[icat]+"ptll",         chTags, ll.pt(),         weight);
	    mon.fillHisto(ctrlCategs[icat]+"mtsum",        chTags, mtsum,           weight);
	    mon.fillHisto(ctrlCategs[icat]+"dilarccosine", chTags, thetall,         weight);
	    mon.fillHisto(ctrlCategs[icat]+"met",          chTags, met.pt(),        weight);
	    //	    mon.fillHisto(ctrlCategs[icat]+"metnotoppt",   chTags, met.pt(),        weight/wgtTopPt); // FIXME: top Pt reweighing 
	    mon.fillHisto(ctrlCategs[icat]+"njets",        chTags, selJets.size(),  weight);
	    //      mon.fillHisto(ctrlCategs[icat]+"njetsnotoppt", chTags, selJets.size(),  weight/wgtTopPt); // FIXME: top Pt reweighing
	    mon.fillHisto(ctrlCategs[icat]+"nbjets",       chTags, selBJets.size(), weight);
	    mon.fillHisto(ctrlCategs[icat]+"ht",           chTags, ht,              weight);
	    mon.fillHisto(ctrlCategs[icat]+"htb",          chTags, htb,             weight);
	    mon.fillHisto(ctrlCategs[icat]+"htnol",        chTags, htnol,           weight);
	    mon.fillHisto(ctrlCategs[icat]+"htbnol",       chTags, htbnol,          weight);
	    
	    for(size_t ijet=0; ijet<selJets.size(); ijet++)
	      {
		if(selJets[ijet].pt()<30 || fabs(selJets[ijet].eta())>2.5) continue;
		TString label("jet"); label+=(ijet+1);
		//		if(ijet+1 >=3) continue; // unnecessary histos at the moment 
		int flavId(selJets[ijet].genflav);
		if(abs(flavId)==5 || abs(flavId)==4 ) flavId=abs(flavId)-1;
		else if(abs(flavId)>6)                flavId=1;
		else if(abs(flavId)==0)               flavId=0;
		else                                  flavId=2;
		mon.fillHisto(ctrlCategs[icat]+"pt"+label+"pt",        chTags, selJets[ijet].pt(), weight);
		mon.fillHisto(ctrlCategs[icat]+"pt"+label+"eta",       chTags, fabs(selJets[ijet].eta()), weight);
		mon.fillHisto(ctrlCategs[icat]+"pt"+label+"flav",      chTags, abs(flavId), weight);
		mon.fillHisto(ctrlCategs[icat]+"pt"+label+"nobsmearpt",chTags, abs(flavId)==5 ? selJets[ijet].pt() : selJets[ijet].jer, weight);
		mon.fillHisto(ctrlCategs[icat]+"pt"+label+"smearpt",   chTags,                                         selJets[ijet].jer, weight);
	      }
	    
	    for(size_t ijet=0; ijet<selBJets.size(); ijet++)
	      {
		TString label("jet"); label+=(ijet+1);
		//		if(ijet+1 >=3) continue; // unnecessary histos at the moment 
		//const data::PhysicsObject_t &genJet=selBJets[ijet].getObject("genJet");
		int flavId(selBJets[ijet].genflav);//genJet.info.find("id")->second;
		if(abs(flavId)==5 || abs(flavId)==4 ) flavId=abs(flavId)-1;
		else if(abs(flavId)>6)                flavId=1;
		else if(abs(flavId)==0)               flavId=0;
		else                                  flavId=2;
		mon.fillHisto(ctrlCategs[icat]+"btag"+label+"pt",        chTags, selBJets[ijet].pt(), weight,true);
		// mon.fillHisto(ctrlCategs[icat]+"btag"+label+"ptnotoppt", chTags, selBJets[ijet].pt(), weight/wgtTopPt,true);  // FIXME: top Pt reweighing
		mon.fillHisto(ctrlCategs[icat]+"btag"+label+"eta",       chTags, fabs(selBJets[ijet].eta()), weight);
		mon.fillHisto(ctrlCategs[icat]+"btag"+label+"flav",      chTags, abs(flavId), weight);
		mon.fillHisto(ctrlCategs[icat]+"btag"+label+"nobsmearpt",chTags, abs(flavId)==5 ? selBJets[ijet].pt() : selBJets[ijet].jer, weight);
		mon.fillHisto(ctrlCategs[icat]+"btag"+label+"smearpt",   chTags,                                       selBJets[ijet].jer, weight);
	      }
	  }
	
	mon.fillHisto("eventflowdileptons",chTags,0,weight);
	//pu control
	mon.fillHisto("nvtx"     ,   chTags, nvtx,      weight);
	mon.fillHisto("nvtxraw"  ,   chTags, nvtx,      weight/puWeight);
	
	//select the event
 	TString var(""); // FIXME: statistical analysis 
	if(passDileptonSelection){     mon.fillHisto("eventflowdileptons"+var, chTags, 1, weight);
	  if( passJetSelection){  mon.fillHisto("eventflowdileptons"+var, chTags, 2, weight);
	    // jets control
	    mon.fillHisto( "njets"  , chTags, selJets.size(), weight);
	    mon.fillHisto( "jetpt"  , chTags, selJets[0].pt(), weight);
	    mon.fillHisto( "jeteta" , chTags, selJets[0].eta(), weight);
	    if(passMetSelection){ mon.fillHisto("eventflowdileptons"+var, chTags, 3, weight);
	      if(isOS){	          mon.fillHisto("eventflowdileptons"+var, chTags, 4, weight);
		float nbtags(selBJets.size());
		if(nbtags>5)
		  mon.fillHisto("finaleventflowdileptons"+var, chTags, 5, weight);
		else
		  mon.fillHisto("finaleventflowdileptons"+var, chTags, nbtags, weight);
		
		if(nbtags>=2){ mon.fillHisto("eventflowdileptons"+var, chTags, 5, weight);
		  if(nbtags>5)
		    mon.fillHisto("finalevtflow2btagsdileptons"+var, chTags, 5, weight);	
		  else 
		    mon.fillHisto("finalevtflow2btagsdileptons"+var, chTags, nbtags, weight);
		}
	      }
	    }
	  }
	}

      }

      
      // Single lepton analysis
      // ----------------------
      if(chTags[1] == "singlemu"){

	if(selSingleLepLeptons.size()<1) continue;

	if(isMC){
	  puWeight          = singleLepLumiWeights->weight(genEv.ngenITpu) * singleLepPUNorm[0];
	  weight            = xsecWeight*puWeight;
	  TotalWeight_plus  = singleLepPuShifters[utils::cmssw::PUUP  ]->Eval(genEv.ngenITpu) * (singleLepPUNorm[2]/singleLepPUNorm[0]);
	  TotalWeight_minus = singleLepPuShifters[utils::cmssw::PUDOWN]->Eval(genEv.ngenITpu) * (singleLepPUNorm[1]/singleLepPUNorm[0]);
	}
	
	// 


	weight *= isMC ? lepEff.getLeptonEfficiency( selSingleLepLeptons[0].pt(), selSingleLepLeptons[0].eta(), abs(selSingleLepLeptons[0].id),  abs(selSingleLepLeptons[0].id) ==11 ? "loose" : "tight" ).first : 1.0;

	// Top Pt reweighting
	double tPt(0.), tbarPt(0.); 
	bool hasTop(false);
	int ngenLeptonsStatus3(0);
	float wgtTopPt(1.0), wgtTopPtUp(1.0), wgtTopPtDown(1.0);
	if(isMC)
	  {
	    for(size_t igen=0; igen<gen.size(); igen++){
	      if(gen[igen].status!=3) continue;
	      int absid=abs(gen[igen].id);
	      if(absid==6){
		hasTop=true;
		if(isTTbarMC){
		  if(gen[igen].id > 0) tPt=gen[igen].pt();
		  else              tbarPt=gen[igen].pt();
		}
	      }
	      if(absid!=11 && absid!=13 && absid!=15) continue;
	      ngenLeptonsStatus3++;
	    }
	    if(mctruthmode==1 && (ngenLeptonsStatus3!=2 || !hasTop)) continue;
	    if(mctruthmode==2 && (ngenLeptonsStatus3==2 || !hasTop)) continue;
	  }
	
	if(tPt>0 && tbarPt>0 && topPtWgt)
	  {
	    topPtWgt->computeWeight(tPt,tbarPt);
	    topPtWgt->getEventWeight(wgtTopPt, wgtTopPtUp, wgtTopPtDown);
	    wgtTopPtUp /= wgtTopPt;
	    wgtTopPtDown /= wgtTopPt;
	  }

	double muontrigeff(1.); 
	utils::cmssw::getSingleMuTrigEff(selSingleLepLeptons[0].pt(), abs(selSingleLepLeptons[0].eta()),muontrigeff);
	weight *= wgtTopPt*muontrigeff;
	// FIXME: add top pt weight syst


	//generator level
	LorentzVector genll;
	for(size_t ig=0; ig<gen.size(); ig++){
	  int pid(abs(gen[ig].id));
	  if(pid!=11 && pid!=13) continue;
	  genll += gen[ig];
	}
	
	//
	//JET/MET ANALYSIS
	//
	llvvJetExtCollection selJets, /*selJetsNoId,*/ selBJets;
	int njets(0), nbjets(0);
	int nTrueJets(0);
	int nTauJets(0);
	
	for(size_t ijet=0; ijet<jets.size(); ijet++) 
	  {
	    //correct jet
	    float toRawSF(jets[ijet].torawsf);
	    LorentzVector rawJet(jets[ijet]*toRawSF);
	    jesCor->setJetEta(rawJet.eta());
	    jesCor->setJetPt(rawJet.pt());
	    jesCor->setJetA(jets[ijet].area);
	    jesCor->setRho(rho);
	    float newJECSF(jesCor->getCorrection());
	    jets[ijet].SetPxPyPzE(rawJet.px(),rawJet.py(),rawJet.pz(),rawJet.energy());
	    jets[ijet] *= newJECSF;
	    jets[ijet].torawsf = 1./newJECSF;
	    if(jets[ijet].pt()<15 || fabs(jets[ijet].eta())>4.7 ) continue;
	    
	    
	    //cross-clean with selected leptons, photons and taus
	    double minDRlj(9999.),minDRlg(9999.), minDRtj(9999.);
	    for(size_t ilep=0; ilep<selSingleLepLeptons.size(); ilep++)
	      minDRlj = TMath::Min( minDRlj, deltaR(jets[ijet],selSingleLepLeptons[ilep]) );
	    for(size_t itau=0; itau<taus.size(); ++itau){
	      if( taus[itau].pt()<20.0 || fabs(taus[itau].eta()) > 2.3) continue;
	      minDRtj = TMath::Min( minDRtj, deltaR(jets[ijet],taus[itau] ) );
	    }
	    if(minDRlj<0.4 || minDRlg<0.4 || minDRtj<0.4) continue;
	    
	    
	    //jet id
	    // float pumva=jets[ijet].puMVA;
	    Int_t idbits(jets[ijet].idbits);
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
	    
	    //selJetsNoId.push_back(jets[ijet]);
	    double jetpt=jets[ijet].pt();
	    
	    if(jets[ijet].pt()>=30 && abs(jets[ijet].eta())<=2.5 && passPFloose && passLoosePuId)
	      nTrueJets++;
	    if(jets[ijet].pt()>=20 && abs(jets[ijet].eta())<=2.5 && passPFloose && passLoosePuId)
	      nTauJets++;
	    
	    if(!passPFloose || !passLoosePuId || jets[ijet].pt()<30 || abs(jets[ijet].eta())>2.5) continue;	    
	    //	    if(passPFloose && passLooseSimplePuId){
	    //	      selJets.push_back(jets[ijet]);
	    //	      if(jets[ijet].pt()>minJetPtToApply) njets++;
	    selJets.push_back(jets[ijet]);
	    njets++;

	    //bjets
	    mon.fillHisto("bjetpt"    ,  chTags, jets[ijet].pt(),  weight);
	    mon.fillHisto("bjetcsv"   ,  chTags, jets[ijet].origcsv,  weight);

	    bool hasCSVV1L(jets[ijet].csv >0.405);
	    bool hasBtagCorr(hasCSVV1L);
	    if(isMC){
	      // set a unique seed
	      double bseed_sin_phi = sin(jets[ijet].phi()*1000000);
	      double bseed = abs(static_cast<int>(bseed_sin_phi*100000));
	      
	      // get jet flavour
	      int bflavid= jets[ijet].genflav;
	      
	      //Initialize
	      BTagSFUtil btsfutil( bseed );
	      
	      TString flavKey("udsg");
	      if(abs(bflavid)==4) flavKey="c";
	      if(abs(bflavid)==5) flavKey="b";
	      std::pair<TString,TString> btagKey("csvL",flavKey);
	      if(btagEffCorr.find(btagKey)!=btagEffCorr.end())
		{
		  TGraphErrors* mceffGr=btagEffCorr[btagKey].first;
		  TGraphErrors* sfGr=btagEffCorr[btagKey].second;
		  if(mceffGr && sfGr){
		    float eff=mceffGr->Eval(jetpt);
		    float sf=sfGr->Eval(jetpt);	
		    if(var == "btagup" || var == "btagdown" || var == "unbtagup" || var == "unbtagdown"){
		      // Apply uncertainty shift
		      float sferr(0.), delta(1000000.);
		      for(size_t ind=0; ind<btagBins.size(); ++ind){
			float tempDelta( fabs(jetpt - btagBins[ind]) );
			if(tempDelta<delta){
			  delta=tempDelta;
			  if      (var == "btagup"   || var == "unbtagup")   sferr=fabs(sfGr->GetErrorYhigh(ind));  // Ensure positive number
			  else if (var == "btagdown" || var == "unbtagdown") sferr=0-fabs(sfGr->GetErrorYlow(ind)); // Ensure negative number
			}
		      }
		      
		      //		      for(int ind=0; ind<sfGr->GetN(); ++ind){
		      //			double xv(0.), yv(0.);
		      //			sfGr->GetPoint(ind,xv,yv);
		      //			float tempDelta(jetpt - xv );
		      //			if(tempDelta<delta){
		      //			  delta=tempDelta; 
		      //			  if      (var == "btagup"   || var == "unbtagup")   sferr=fabs(sfGr->GetErrorYhigh(ind));  // Ensure positive number
		      //			  else if (var == "btagdown" || var == "unbtagdown") sferr=0-fabs(sfGr->GetErrorYlow(ind)); // Ensure negative number
		      //			}
		      //		      }
		      sf+=sferr;
		    }
		    btsfutil.modifyBTagsWithSF(hasBtagCorr, sf, eff);	    
		  }
		}
	    }
	    //	    if(!hasCSVV1L) continue;
	    if(!hasBtagCorr) continue;

	      // if(jets[ijet].pt()>20 && fabs(jets[ijet].eta())<2.4 && jets[ijet].origcsv>0.898){
	    selBJets.push_back(jets[ijet]);  
	    nbjets++;
	  }
	std::sort(selJets.begin(), selJets.end(), sort_llvvObjectByPt);
	std::sort(selBJets.begin(), selBJets.end(), sort_llvvObjectByPt);
	mon.fillHisto("nbjets"    ,  chTags, nbjets,  weight);
	
	//
	// TAU ANALYSIS
	//
	llvvTauCollection selTaus;
	for(size_t itau=0; itau<taus.size(); itau++){
	  llvvTau& tau = taus[itau];
	  if(tau.pt()<20.0 || fabs(tau.eta()) > 2.3)continue; 
	  
	  
	  bool overalWithLepton(false);
	  for(int l1=0   ;l1<(int)selSingleLepLeptons.size();l1++){
	    if(deltaR(tau, selSingleLepLeptons[l1])<0.1){overalWithLepton=true; break;}
	  }
	  if(overalWithLepton)continue;
	  
	  //         printf("TauId: "); for(unsigned int i=0;i<64;i++){printf("%i ", (int) ((tau.idbits>>i)&1));}printf("\n");
	  
	  if(!tau.isPF) continue; // We want PF taus
	  if(abs(tau.dZ)>=0.5) continue;  
	  if( tau.emfraction >= 2. /*0.95*/ ) continue;
	  if(abs(tau.id/15.0) !=1) continue; // Non non-1 taus. Actually this should be always ok 
	  
	  
	  if(!tau.passId(llvvTAUID::againstElectronMediumMVA5))continue;
	  if(!tau.passId(llvvTAUID::againstMuonTight2))continue; 
	  if(!tau.passId(llvvTAUID::decayModeFinding))continue;
	  if(!tau.passId(llvvTAUID::byMediumCombinedIsolationDeltaBetaCorr3Hits))continue;
	  
	  selTaus.push_back(tau);         
	}
	
	
	// NOW FOR THE CONTROL PLOTS
	//
	//start analysis
	//      for(size_t ich=0; ich<chTags.size(); ich++)
	//	{
	//	  std::vector<TString> tags(1,chTags[ich]);
	
	bool otherTriggersVeto( eeTrigger || emuTrigger || mumuTrigger ); // It should be already excluded by the channel requirement  // No. emu and mumu trigger may fire but we select the event for mutau  
	bool additionalLeptonsVeto( nVetoE>0 || nVetoMu>0 );
	bool passMuonPlusJets( muTrigger && selSingleLepLeptons.size()==1 /*&& !otherTriggersVeto*/ && !additionalLeptonsVeto && nTrueJets>1 && nTauJets>2 );    
	bool passMet(met.pt()>40.);
	bool pass1bjet(nbjets>0);
	bool pass1tau(selTaus.size()==1);
	bool passOS(true);
	if(pass1tau)  passOS = ((selTaus[0].id)*(selSingleLepLeptons[0].id)<0);      
	
	if(passMuonPlusJets)                                               mon.fillHisto("eventflowsinglelepton",chTags,0,weight);
	if(passMuonPlusJets && passMet)                                    mon.fillHisto("eventflowsinglelepton",chTags,1,weight);
	if(passMuonPlusJets && passMet && pass1bjet)                       mon.fillHisto("eventflowsinglelepton",chTags,2,weight);
	if(passMuonPlusJets && passMet && pass1bjet && pass1tau)           mon.fillHisto("eventflowsinglelepton",chTags,3,weight);
	if(passMuonPlusJets && passMet && pass1bjet && pass1tau && passOS) mon.fillHisto("eventflowsinglelepton",chTags,4,weight);
	
	
	if(passMuonPlusJets){
	  //pu control
	  mon.fillHisto("nvtx"     ,   chTags, nvtx,      weight);
	  mon.fillHisto("nvtxraw"  ,   chTags, nvtx,      weight/puWeight);
	  
	  // jets control
	  mon.fillHisto( "njets"  , chTags, selJets.size(), weight);
	  mon.fillHisto( "jetpt"  , chTags, selJets[0].pt(), weight);
	  mon.fillHisto( "jeteta" , chTags, selJets[0].eta(), weight);
	  // FIXME: add exclusive plots
	}
	
	
	if(passMuonPlusJets && passMet && pass1bjet && pass1tau){
	  mon.fillHisto("ntaus_1tau"        ,  chTags, selTaus.size(), weight);
	  //mon.fillHisto("tauleadpt_1tau"    ,  chTags, selTaus.size()>0?selTaus[0].pt():-1,  weight);
	  //mon.fillHisto("tauleadeta_1tau"   ,  chTags, selTaus.size()>0?selTaus[0].eta():-10, weight);
	  mon.fillHisto("taupt"        , chTags, selTaus[0].pt(),  weight);
	  mon.fillHisto("taueta"       , chTags, selTaus[0].eta(),  weight);
	  mon.fillHisto("taucharge"    , chTags, selTaus[0].id/15.0, weight);
	  mon.fillHisto("taudz"        , chTags, selTaus[0].vz,  weight);
	  mon.fillHisto("tauvz"        , chTags, selTaus[0].z_expo,  weight);
	  mon.fillHisto("tauemfraction", chTags, selTaus[0].emfraction,  weight);
	  mon.fillHisto("taudizeta"    , chTags, selTaus[0].dZ,  weight);
	  
	  if(passOS){
	    mon.fillHisto("ntausos"        ,  chTags, selTaus.size(), weight);
	    //mon.fillHisto("tauleadptos"    ,  chTags, selTaus.size()>0?selTaus[0].pt():-1,  weight);
	    //mon.fillHisto("tauleadetaos"   ,  chTags, selTaus.size()>0?selTaus[0].eta():-10, weight);
	    mon.fillHisto("tauptos"        ,  chTags, selTaus[0].pt(),  weight);
	    mon.fillHisto("tauetaos"       , chTags, selTaus[0].eta(),  weight);
	    mon.fillHisto("tauchargeos"    ,  chTags, selTaus[0].id/15.0, weight);
	    mon.fillHisto("taudzos"        ,  chTags, selTaus[0].vz,  weight);
	    mon.fillHisto("tauvzos"        ,  chTags, selTaus[0].z_expo,  weight);
	    mon.fillHisto("tauemfractionos", chTags, selTaus[0].emfraction,  weight);
	    mon.fillHisto("taudizetaos"    , chTags, selTaus[0].dZ,  weight);
	    
	  }
	}
      }
      
      
      //       //
      //       // HIGGS ANALYSIS
      //       //
      // 
      //       LorentzVector higgsCand;
      //       int higgsCandId=0,  higgsCandMu=-1, higgsCandEl=-1, higgsCandT1=-1, higgsCandT2=-1;
      // 
      //       bool passNLep = (selLeptons.size()+selTaus.size())>=4;
      // 
      //       //Check if the event is compatible with a Mu-El candidate
      //       for(int l1=0   ;!higgsCandId && l1<(int)selLeptons.size();l1++){
      //       for(int l2=l1+1;!higgsCandId && l2<(int)selLeptons.size();l2++){
      //          if(l1==dilLep1 || l1==dilLep2 || l2==dilLep1 || l2==dilLep2)continue; //lepton already used in the dilepton pair
      //          if(selLeptons[l1].id*selLeptons[l2].id!=-143)continue;//Only consider opposite sign, opposite flavor pairs
      // 
      //          //check all deltaRs
      //          if(deltaR(selLeptons[l1], selLeptons[dilLep1])<0.1)continue;
      //          if(deltaR(selLeptons[l1], selLeptons[dilLep2])<0.1)continue;
      //          if(deltaR(selLeptons[l2], selLeptons[dilLep1])<0.1)continue;
      //          if(deltaR(selLeptons[l2], selLeptons[dilLep2])<0.1)continue;
      //          if(deltaR(selLeptons[l1], selLeptons[l2     ])<0.1)continue;
      // 
      // //         if(utils::cmssw::relIso(selLeptons[l1], rho)>0.25)continue;
      // //         if(utils::cmssw::relIso(selLeptons[l2], rho)>0.25)continue;
      // 
      // //         if((selLeptons[l1].pt() + selLeptons[l2].pt())<25)continue;
      //  
      //          int muId, elId;
      //          if(abs(selLeptons[l1].id)==13){muId=l1; elId=l2;}else{muId=l2; elId=l1;}
      // 
      //          higgsCand=selLeptons[muId]+selLeptons[elId]; higgsCandId=selLeptons[muId].id*selLeptons[elId].id;  higgsCandMu=muId; higgsCandEl=elId;
      //          break;//we found a candidate, stop the loop
      //       }}
      // 
      //       //Check if the event is compatible with a Lep-Tau candidate
      //       for(int l1=0;!higgsCandId && l1<(int)selLeptons.size();l1++){
      //       for(int t1=0;!higgsCandId && t1<(int)selTaus   .size();t1++){
      //          if(l1==dilLep1 || l1==dilLep2)continue; //lepton already used in the dilepton pair
      //          if(selLeptons[l1].id*selTaus[t1].id>=0)continue;//Only consider opposite sign pairs
      //          if(selTaus[t1].pt()<15 || fabs(selTaus[t1].eta())>2.3)continue;
      // 
      //          //check all deltaRs
      //          if(deltaR(selLeptons[l1], selLeptons[dilLep1])<0.1)continue;
      //          if(deltaR(selLeptons[l1], selLeptons[dilLep2])<0.1)continue;
      //          if(deltaR(selTaus[t1]   , selLeptons[dilLep1])<0.1)continue;
      //          if(deltaR(selTaus[t1]   , selLeptons[dilLep2])<0.1)continue;
      //          if(deltaR(selLeptons[l1], selTaus[t1        ])<0.1)continue;
      // 
      // //         if(utils::cmssw::relIso(selLeptons[l1], rho)>0.25)continue;
      // 
      //          if(abs(selLeptons[l1].id)==11 && (!selTaus[t1].passId(llvvTAUID::decayModeFinding) || !selTaus[t1].passId(llvvTAUID::againstElectronTightMVA3) || !selTaus[t1].passId(llvvTAUID::againstMuonLoose2) || !selTaus[t1].passId(llvvTAUID::byLooseCombinedIsolationDeltaBetaCorr3Hits) ))continue;
      //          if(abs(selLeptons[l1].id)!=11 && (!selTaus[t1].passId(llvvTAUID::decayModeFinding) || !selTaus[t1].passId(llvvTAUID::againstElectronLoose    ) || !selTaus[t1].passId(llvvTAUID::againstMuonTight2) || !selTaus[t1].passId(llvvTAUID::byLooseCombinedIsolationDeltaBetaCorr3Hits) ))continue;
      // 
      // //         if((selLeptons[l1].pt()+selTaus[t1].pt())<45 )continue;
      // 
      //          higgsCand=selLeptons[l1]+selTaus[t1]; higgsCandId=selLeptons[l1].id*selTaus[t1].id;  if(abs(selLeptons[l1].id)==11){higgsCandEl=l1;}else{higgsCandMu=l1;} higgsCandT1=t1;
      //          break;//we found a candidate, stop the loop         
      //       }}
      // 
      //       //Check if the event is compatible with a Tau-Tau candidate
      //       for(int t1=0   ;!higgsCandId && t1<(int)selTaus   .size();t1++){
      //       for(int t2=t1+1;!higgsCandId && t2<(int)selTaus   .size();t2++){
      // 
      //          //printf("Tau %6.2f %6.2f %+3i %i %i %i %i\n", selTaus[t1].pt(), fabs(selTaus[t1].eta()), selTaus[t1].id, (int)selTaus[t1].passId(llvvTAUID::decayModeFinding), (int)selTaus[t1].passId(llvvTAUID::againstElectronLoose), (int) selTaus[t1].passId(llvvTAUID::againstMuonLoose2), (int) selTaus[t1].passId(llvvTAUID::byLooseCombinedIsolationDeltaBetaCorr3Hits) );
      //          //printf("    %6.2f %6.2f %+3i %i %i %i %i\n", selTaus[t2].pt(), fabs(selTaus[t2].eta()), selTaus[t2].id, (int)selTaus[t2].passId(llvvTAUID::decayModeFinding), (int)selTaus[t2].passId(llvvTAUID::againstElectronLoose), (int) selTaus[t2].passId(llvvTAUID::againstMuonLoose2), (int) selTaus[t2].passId(llvvTAUID::byLooseCombinedIsolationDeltaBetaCorr3Hits) );
      // 
      //          if(selTaus[t1].pt()<15 || fabs(selTaus[t1].eta())>2.3)continue;
      //          if(selTaus[t2].pt()<15 || fabs(selTaus[t2].eta())>2.3)continue;
      // 
      //          //check all deltaRs
      //          if(deltaR(selTaus[t1]   , selLeptons[dilLep1])<0.1)continue;
      //          if(deltaR(selTaus[t1]   , selLeptons[dilLep2])<0.1)continue;
      //          if(deltaR(selTaus[t2]   , selLeptons[dilLep1])<0.1)continue;
      //          if(deltaR(selTaus[t2]   , selLeptons[dilLep2])<0.1)continue;
      //          if(deltaR(selTaus[t1]   , selTaus[t2        ])<0.1)continue;
      // 
      //          if(!selTaus[t1].passId(llvvTAUID::decayModeFinding) || !selTaus[t1].passId(llvvTAUID::againstElectronLoose) || !selTaus[t1].passId(llvvTAUID::againstMuonLoose2) || !selTaus[t1].passId(llvvTAUID::byLooseCombinedIsolationDeltaBetaCorr3Hits) )continue;
      //          if(!selTaus[t2].passId(llvvTAUID::decayModeFinding) || !selTaus[t2].passId(llvvTAUID::againstElectronLoose) || !selTaus[t2].passId(llvvTAUID::againstMuonLoose2) || !selTaus[t2].passId(llvvTAUID::byLooseCombinedIsolationDeltaBetaCorr3Hits) )continue;
      // 
      // //         if((selTaus[t1].pt()+selTaus[t2].pt())<75)continue;
      // 
      //          if(selTaus[t1].id*selTaus[t2].id<0){        
      //             higgsCand=selTaus[t1]+selTaus[t2]; higgsCandId=selTaus[t1].id*selTaus[t2].id;  higgsCandT1=t1; higgsCandT2=t2;
      //          }else{
      //             higgsCand=selTaus[t1]+selTaus[t2]; higgsCandId=selTaus[t1].id;  higgsCandT1=t1; higgsCandT2=t2;            
      //          }
      //          break;//we found a candidate, stop the loop
      //       }}
      // 
      // 
      //       //apply data/mc correction factors
      //       if(higgsCandMu!=-1)weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[higgsCandMu].pt(), selLeptons[higgsCandMu].eta(), abs(selLeptons[higgsCandMu].id),  abs(selLeptons[higgsCandMu].id) ==11 ? "loose" : "loose" ).first : 1.0;
      //       if(higgsCandEl!=-1)weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[higgsCandEl].pt(), selLeptons[higgsCandEl].eta(), abs(selLeptons[higgsCandEl].id),  abs(selLeptons[higgsCandEl].id) ==11 ? "loose" : "loose" ).first : 1.0;
      // 
      //       //check the channel
      //       //prepare the tag's vectors for histo filling
      //       bool passHiggs = abs(higgsCandId)>100;  int HiggsShortId = abs(selLeptons[dilLep1].id)==13?0:4;
      //            if( abs(higgsCandId)==143 ){ chTags.push_back(chTags[chTags.size()-1] + string("_elmu")); HiggsShortId+=0;}
      //       else if( abs(higgsCandId)==165 ){ chTags.push_back(chTags[chTags.size()-1] + string("_elha")); HiggsShortId+=1;}
      //       else if( abs(higgsCandId)==195 ){ chTags.push_back(chTags[chTags.size()-1] + string("_muha")); HiggsShortId+=2;}
      //       else if( abs(higgsCandId)==225 ){ chTags.push_back(chTags[chTags.size()-1] + string("_haha")); HiggsShortId+=3;}
      //       else if( abs(higgsCandId)== 15 ){ chTags.push_back(chTags[chTags.size()-1] + string("_haCtrl"));}
      //       else                              chTags.push_back(chTags[chTags.size()-1] + string("_none"));
      // 
      //       bool passLepVeto  = true;
      //       for(int l1=0;l1<(int)selLeptons.size();l1++){
      //          mon.fillHisto("leppt"    ,  chTags, selLeptons[l1].pt(),  weight);
      //          if(l1==dilLep1 || l1==dilLep2 || l1==higgsCandMu || l1==higgsCandEl)continue; //lepton already used in the dilepton pair or higgs candidate
      //          passLepVeto = false; break;
      //       }
      //       for(int t1=0;passLepVeto && t1<(int)selTaus   .size();t1++){         
      //          mon.fillHisto("taupt"    ,  chTags, selTaus[t1].pt(),  weight);
      //          mon.fillHisto("taucharge",  chTags, selTaus[t1].id/15.0, weight);
      //          mon.fillHisto("taudz"    ,  chTags, selTaus[t1].vz,  weight);
      //          mon.fillHisto("tauvz"    ,  chTags, selTaus[t1].z_expo,  weight);
      // 
      //          if(t1==higgsCandT1 || t1==higgsCandT2)continue; //lepton already used in the dilepton pair or higgs candidate
      //          if(selTaus[t1].pt()<20)continue;
      //          if(!selTaus[t1].passId(llvvTAUID::againstElectronLoose) || !selTaus[t1].passId(llvvTAUID::againstMuonLoose2) || !selTaus[t1].passId(llvvTAUID::byMediumCombinedIsolationDeltaBetaCorr3Hits)  )continue;
      //          passLepVeto = false; break;          
      //       } 
      // 
      //       bool passBJetVeto = true;;
      //       for(int j1=0;j1<(int)selBJets.size();j1++){
      //          if(                   deltaR(selBJets[j1]   , selLeptons[dilLep1    ])>0.4){passBJetVeto=false; break;}
      //          if(                   deltaR(selBJets[j1]   , selLeptons[dilLep2    ])>0.4){passBJetVeto=false; break;}
      //          if(higgsCandMu!=-1 && deltaR(selBJets[j1]   , selLeptons[higgsCandMu])>0.4){passBJetVeto=false; break;}
      //          if(higgsCandEl!=-1 && deltaR(selBJets[j1]   , selLeptons[higgsCandEl])>0.4){passBJetVeto=false; break;}
      //          if(higgsCandT1!=-1 && deltaR(selBJets[j1]   , selTaus   [higgsCandT1])>0.4){passBJetVeto=false; break;}
      //          if(higgsCandT2!=-1 && deltaR(selBJets[j1]   , selTaus   [higgsCandT2])>0.4){passBJetVeto=false; break;}
      //       }
      // 
      //       int NCleanedJet=0;
      //       for(int j1=0;j1<(int)selJets.size();j1++){
      //          if(                   deltaR(selJets[j1]   , selLeptons[dilLep1    ])<0.4)continue;
      //          if(                   deltaR(selJets[j1]   , selLeptons[dilLep2    ])<0.4)continue;
      //          if(higgsCandMu!=-1 && deltaR(selJets[j1]   , selLeptons[higgsCandMu])<0.4)continue;
      //          if(higgsCandEl!=-1 && deltaR(selJets[j1]   , selLeptons[higgsCandEl])<0.4)continue;
      //          if(higgsCandT1!=-1 && deltaR(selJets[j1]   , selTaus   [higgsCandT1])<0.4)continue;
      //          if(higgsCandT2!=-1 && deltaR(selJets[j1]   , selTaus   [higgsCandT2])<0.4)continue;
      //          NCleanedJet++;
      //       }
      // 
      // 
      // 
      //       //SVFIT MASS
      //       double diTauMass = -1;
      //       if(passHiggs && passLepVeto && passBJetVeto){
      //             //taken from https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorkingSummer2013
      //             TMatrixD covMET(2, 2); // PFMET significance matrix
      //             covMET[0][0] = met.sigx2;
      //             covMET[0][1] = met.sigxy;
      //             covMET[1][0] = met.sigxy;
      //             covMET[1][1] = met.sigy2;
      //             std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
      //             if(higgsCandMu!=-1)measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, NSVfitStandalone::LorentzVector(selLeptons[higgsCandMu].px(), selLeptons[higgsCandMu].py(), selLeptons[higgsCandMu].pz(), selLeptons[higgsCandMu].E()) ));
      //             if(higgsCandEl!=-1)measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, NSVfitStandalone::LorentzVector(selLeptons[higgsCandEl].px(), selLeptons[higgsCandEl].py(), selLeptons[higgsCandEl].pz(), selLeptons[higgsCandEl].E()) ));
      //             if(higgsCandT1!=-1)measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kHadDecay, NSVfitStandalone::LorentzVector(selTaus   [higgsCandT1].px(), selTaus   [higgsCandT1].py(), selTaus   [higgsCandT1].pz(), selTaus   [higgsCandT1].E()) ));
      //             if(higgsCandT2!=-1)measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kHadDecay, NSVfitStandalone::LorentzVector(selTaus   [higgsCandT2].px(), selTaus   [higgsCandT2].py(), selTaus   [higgsCandT2].pz(), selTaus   [higgsCandT2].E()) ));
      //             NSVfitStandaloneAlgorithm algo(measuredTauLeptons, NSVfitStandalone::Vector(met.px(), met.py(), 0) , covMET, 0);
      //             algo.addLogM(false);
      //             algo.integrateMarkovChain();
      //             //algo.integrateVEGAS(); ////Use this instead for VEGAS integration
      //             if(algo.isValidSolution()){
      //                diTauMass     = algo.getMass(); 
      //                double diTauMassErr = algo.massUncert();
      //             }
      //       }
      // 
      
      


      //      mon.fillHisto("eventflow",chTags,0,weight);
      //      if(passZmass)                                              mon.fillHisto("eventflow",chTags,1,weight);
      //      if(passZmass && passNLep)                                  mon.fillHisto("eventflow",chTags,2,weight);
      //      if(passZmass && passHiggs)                                 mon.fillHisto("eventflow",chTags,3,weight);
      //      if(passZmass && passHiggs && passLepVeto)                  mon.fillHisto("eventflow",chTags,4,weight);
      //      if(passZmass && passHiggs && passLepVeto && passBJetVeto)  mon.fillHisto("eventflow",chTags,5,weight);
      //      if(passZmass && passHiggs && passLepVeto && passBJetVeto)  mon.fillHisto("eventflow",chTags,7+HiggsShortId,weight);
      //      
      //      mon.fillHisto("zmass",    chTags, zll.mass(), weight);  
      //
      //        if(passZmass && passBJetVeto && passLepVeto && !passHiggs && abs(higgsCandId)== 15){
      //            mon.fillHisto("higgsmass"    , chTags, higgsCand.mass(),  weight);
      //           
      //            mon.fillHisto("taufakerate"     ,   chTags, 1,  selTaus   [higgsCandT1].jet.pt(),    weight);
      //            mon.fillHisto("taufakerate"     ,   chTags, 1,  selTaus   [higgsCandT2].jet.pt(),    weight);
      //
      //            for(size_t itau=0; itau<taus.size(); itau++){
      //               llvvTau& tau = taus[itau];
      //               if(tau.pt()<15.0 || fabs(tau.eta()) > 2.3)continue;
      //               mon.fillHisto("taufakerate"     ,   chTags, 0,  tau.jet.pt(),    weight);
      //            }
      //
      //        }
      //
      ////	  if(passZmass){
      //        if(passZmass && passBJetVeto && passLepVeto && passHiggs){
      //	
      //	    //pu control
      //            mon.fillHisto("nvtx"     ,   chTags, nvtx,      weight);
      //	    mon.fillHisto("nvtxraw"  ,   chTags, nvtx,      weight/puWeight);
      //	    mon.fillHisto("rho"      ,   chTags, rho,       weight);
      //            mon.fillHisto("rho25"    ,   chTags, rho25,     weight);
      //	
      //	    //Z kinematics control
      //	    mon.fillHisto("leadpt"      ,   chTags, leadingLep.pt(), weight);      
      //	    mon.fillHisto("trailerpt"   ,   chTags, trailerLep.pt(), weight);      
      //	    mon.fillHisto("leadeta"     ,   chTags, TMath::Max(fabs(leadingLep.eta()),fabs(trailerLep.eta())), weight);      
      //	    mon.fillHisto("trailereta"  ,   chTags, TMath::Min(fabs(leadingLep.eta()),fabs(trailerLep.eta())), weight);      
      //
      //	    mon.fillHisto("zpt"      , chTags, zll.pt(),  weight);      
      //	    mon.fillHisto("zeta"     , chTags, zll.eta(), weight);
      //	    mon.fillHisto("zy"       , chTags, zy,        weight);
      //	    
      ////	    if(passZpt && passZeta){
      //          if(passBJetVeto && passLepVeto){
      //	  
      //	      //analyze dilepton kinematics
      //	      mon.fillHisto("leadeta"   ,  chTags, leadingLep.eta(), weight);
      //	      mon.fillHisto("leadpt"    ,  chTags, leadingLep.pt(),  weight);
      //	      mon.fillHisto("trailereta",  chTags, trailerLep.eta(), weight);
      //	      mon.fillHisto("trailerpt",   chTags, trailerLep.pt(),  weight);
      // 
      //
      //              mon.fillHisto("ntaus"        ,  chTags, selTaus.size(), weight);
      //              mon.fillHisto("tauleadpt"    ,  chTags, selTaus.size()>0?selTaus[0].pt():-1,  weight);
      //              mon.fillHisto("tauleadeta"   ,  chTags, selTaus.size()>0?selTaus[0].eta():-10, weight);
      //
      //              if(passHiggs){
      //                 mon.fillHisto("higgspt"      , chTags, higgsCand.pt(),    weight);
      //                 mon.fillHisto("higgsmass"    , chTags, higgsCand.mass(),  weight);
      //                 mon.fillHisto("higgsmasssvfit", chTags, diTauMass,  weight);
      //                 mon.fillHisto("higgsnjets"   , chTags, NCleanedJet      , weight); 
      //                 mon.fillHisto("higgsmet"     , chTags, met.pt()         , weight);
      //              }
      //
      //
      //
      //
      
      //	      //STATISTICAL ANALYSIS
      //	      float Q2Weight_plus(1.0), Q2Weight_down(1.0);
      //	      float PDFWeight_plus(1.0), PDFWeight_down(1.0);
      //	      if(isSignal){
      ////		  if(mPDFInfo){
      ////		      std::vector<float> wgts=mPDFInfo->getWeights(iev);
      ////		      for(size_t ipw=0; ipw<wgts.size(); ipw++){
      ////			  PDFWeight_plus = TMath::Max(PDFWeight_plus,wgts[ipw]);
      ////			  PDFWeight_down = TMath::Min(PDFWeight_down,wgts[ipw]);
      ////			}
      ////		  }
      //	      }
      //
      //	      for(size_t ivar=0; ivar<nvarsToInclude; ivar++){
      //		float iweight = weight;                                                //nominal
      //		if(ivar==5)                        iweight *= TotalWeight_plus;        //pu up
      //		if(ivar==6)                        iweight *= TotalWeight_minus;       //pu down
      //		if(ivar==7)                        iweight *= Q2Weight_plus;
      //		if(ivar==8)                        iweight *= Q2Weight_down;
      //		if(ivar==9)                        iweight *= PDFWeight_plus;
      //		if(ivar==10)                       iweight *= PDFWeight_down;
      //
      ///*	    
      //		llvvJetExtCollection localSelJets;
      //		for(size_t ijet=0; ijet<jets.size(); ijet++){
      //	      
      //		  float rawpt=jets[ijet].pt();
      //		  float pt=rawpt;
      //		  if(ivar==1) pt=jets[ijet].jesup;
      //		  if(ivar==2) pt=jets[ijet].jesdown;
      //		  if(ivar==3) pt=jets[ijet].jerup;
      //		  if(ivar==4) pt=jets[ijet].jerdown;
      //		  if(pt<minJetPtToApply || fabs(jets[ijet].eta())>4.7) continue;
      //	      
      //		  Int_t idbits=jets[ijet].idbits;
      //		  bool passPFloose ( ((idbits>>0) & 0x1) );
      //		  //int puId((idbits>>3) & 0xf);
      //		  //bool passLoosePuId( ( puId >> 2) & 0x1);
      //		  int simplePuId( ( idbits >>7 ) & 0xf );
      //		  bool passLooseSimplePuId(  ( simplePuId >> 2) & 0x1);
      //		  if(!passPFloose || !passLooseSimplePuId) continue;
      //
      //		  llvvJetExt iSelJet(jets[ijet]);
      //		  iSelJet *= pt/rawpt;
      //		  localSelJets.push_back( iSelJet );
      //		}
      //		if(localSelJets.size()<2)  continue;
      //		std::sort(localSelJets.begin(), localSelJets.end(),  sort_llvvObjectByPt);
      //
      //		//recoil residuals uncertainty
      //		if( (ivar==11 || ivar==12) && recoilResidualsGr.size()==2)
      //		  {
      //		    for(size_t ijet=0; ijet<2; ijet++)
      //		      {
      //			float abseta=TMath::Abs(localSelJets[ijet].eta());
      //			float sfEnvelope=TMath::Max( TMath::Abs(1-recoilResidualsGr[0]->Eval(abseta)),
      //						     TMath::Abs(1-recoilResidualsGr[1]->Eval(abseta)) );
      //			localSelJets[ijet] *= (ivar==11 ? 1+sfEnvelope : 1-sfEnvelope);
      //		      }
      //		  }
      //*/
      //		//re-assign the event category;
      //		std::vector<TString> locTags = chTags;
      //
      //                double sumPt = 0;                  
      //                if(higgsCandMu!=-1) sumPt+= selLeptons[higgsCandMu].pt();
      //                if(higgsCandEl!=-1) sumPt+= selLeptons[higgsCandEl].pt();
      //                if(higgsCandT1!=-1) sumPt+= selTaus   [higgsCandT1].pt();
      //                if(higgsCandT2!=-1) sumPt+= selTaus   [higgsCandT2].pt();
      //
      //                double Lep1Iso=-1, Lep2Iso=-1;
      //                if(higgsCandMu!=-1)Lep1Iso = utils::cmssw::relIso(selLeptons[higgsCandMu], rho);
      //                if(higgsCandEl!=-1)Lep2Iso = utils::cmssw::relIso(selLeptons[higgsCandEl], rho);
      //
      //		for(unsigned int index=0; index<optim_Cuts_sumpt.size();index++)
      //		  {
      //                    bool passSumPt = sumPt>=optim_Cuts_sumpt[index];
      //                    bool passIso   = Lep1Iso<=optim_Cuts_lepIso[index] && Lep2Iso<=optim_Cuts_lepIso[index];
      //
      ////		    TString mjjCat("");
      ////		    std::vector<TString> localSelTags=getDijetCategories(mjj,detajj,locTags,mjjCat);
      //                    std::vector<TString> localSelTags=locTags;
      //                    if(passSumPt && passIso)mon.fillHisto(TString("svfit_shapes")+varNames[ivar],localSelTags,index,diTauMass,iweight);		    
      //		  }
      //
      //	      }
      //	    }//end passZpt && passZeta
      //
      //	  }//end passZmass
      //	}
      //
    }
  
  printf("\n"); 
  //  file->Close();
  
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
}  





