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

#include "FWCore/FWLite/interface/FWLiteEnabler.h"
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
#include "UserCode/llvv_fwk/interface/LeptonEfficiencySF.h"
#include "UserCode/llvv_fwk/interface/PDFInfo.h"
#include "UserCode/llvv_fwk/interface/rochcor2015.h"
#include "UserCode/llvv_fwk/interface/muresolution_run2.h"
#include "UserCode/llvv_fwk/interface/BTagCalibrationStandalone.h"
#include "UserCode/llvv_fwk/interface/BtagUncertaintyComputer.h"
#include "UserCode/llvv_fwk/interface/GammaWeightsHandler.h"

#include "UserCode/llvv_fwk/interface/PatUtils.h"
#include "UserCode/llvv_fwk/interface/TrigUtils.h"
#include "UserCode/llvv_fwk/interface/EwkCorrections.h"
#include "UserCode/llvv_fwk/interface/ZZatNNLO.h"

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

#include <time.h>

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
//  AutoLibraryLoader::enable();
  FWLiteEnabler::enable();

  
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

  bool filterOnlyEE(false), filterOnlyMUMU(false), filterOnlyEMU(false), filterOnlyPhoton(false), filterOnlyE(false), filterOnlyMU(false);
  if(!isMC){
      if(dtag.Contains("DoubleEle"))   filterOnlyEE=true;
      if(dtag.Contains("DoubleMu"))    filterOnlyMUMU=true;
      if(dtag.Contains("MuEG"))        filterOnlyEMU=true;
      if(dtag.Contains("SinglePhoton"))filterOnlyPhoton=true;     
      if(dtag.Contains("SingleMu"))    filterOnlyE=true;      
      if(dtag.Contains("SingleElectron"))filterOnlyMU=true;      
  }
  bool isV0JetsMC(false);//isMC && (dtag.Contains("DYJetsToLL_50toInf") || dtag.Contains("_WJets")));  #FIXME should be reactivated as soon as we have exclusive jet samples
  bool isWGmc(isMC && dtag.Contains("WG"));
  bool isZGmc(isMC && dtag.Contains("ZG"));
  bool isMC_GG  = isMC && ( string(dtag.Data()).find("GG" )  != string::npos);
  bool isMC_VBF = isMC && ( string(dtag.Data()).find("VBF")  != string::npos);
  bool isMC_125OnShell = isMC && (mctruthmode==521);
  if(isMC_125OnShell) mctruthmode=125;
  bool isMC_ZZ  = isMC && ( string(dtag.Data()).find("TeV_ZZ")  != string::npos);
  bool isMC_ZZ2l2nu  = isMC && ( string(dtag.Data()).find("TeV_ZZ2l2nu")  != string::npos);
  bool isMC_WZ  = isMC && ( string(dtag.Data()).find("TeV_WZ")  != string::npos);
  bool isMC_QCD = (isMC && dtag.Contains("QCD"));
  bool isMC_GJet = (isMC && dtag.Contains("GJet"));
  bool is2015data = (!isMC && dtag.Contains("2015")); 
  bool is2016data = (!isMC && dtag.Contains("2016")); 

  //tree info
  TString dirname = runProcess.getParameter<std::string>("dirName");

  //systematics
  bool runSystematics                        = runProcess.getParameter<bool>("runSystematics");
  std::vector<TString> varNames(1,"");

  std::vector<string> jetVarNames = {""};//, "_scale_jup","_scale_jdown", "_res_jup", "_res_jdown"};


  if(runSystematics){
     if(true){
        varNames.push_back("_scale_umetup"); varNames.push_back("_scale_umetdown");    //unclustered met
        varNames.push_back("_res_jup");      varNames.push_back("_res_jdown");    //jet energy resolution
        varNames.push_back("_scale_jup");    varNames.push_back("_scale_jdown");  //jet energy scale
        varNames.push_back("_scale_mup");    varNames.push_back("_scale_mdown");  //muon energy scale
        varNames.push_back("_scale_eup");    varNames.push_back("_scale_edown");  //electron energy scale
        varNames.push_back("_puup");         varNames.push_back("_pudown");      //pileup uncertainty 
        varNames.push_back("_eff_bup");      varNames.push_back("_eff_bdown");    //btag veto
        varNames.push_back("_lepveto");                                           //3rd lepton veto
        varNames.push_back("_th_factup");    varNames.push_back("_th_factdown"); //factorization and renormalization scales
        varNames.push_back("_th_pdf");                                           //pdf
        varNames.push_back("_th_alphas");                                         //alpha_s (QCD)
     }
     if(isMC_GG || isMC_VBF){
        varNames.push_back("_signal_lshapeup"); varNames.push_back("_signal_lshapedown"); //signal line shape (NNLO + interf)
        varNames.push_back("_signal_normup"); varNames.push_back("_signal_normdown"); //signal scale      (NNLO + interf)
     }
     if(isMC_ZZ){
        varNames.push_back("_th_ewkup"); varNames.push_back("_th_ewkdown"); //EWK+QCD corrections
     }
  }
  size_t nvarsToInclude=varNames.size();
  
  std::vector<std::string> allWeightsURL=runProcess.getParameter<std::vector<std::string> >("weightsFile");
  std::string weightsDir( allWeightsURL.size() ? allWeightsURL[0] : "");

  std::vector<std::string> gammaPtWeightsFiles =  runProcess.getParameter<std::vector<std::string> >("weightsFile");      
  GammaWeightsHandler* gammaWgtHandler = (gammaPtWeightsFiles.size()>0 && gammaPtWeightsFiles[0]!="") ? new GammaWeightsHandler(runProcess,"",true) : NULL;
  if(gammaWgtHandler)printf("gammaWgtHandler is activated\n");

  //HIGGS weights and uncertainties
  
  //narrow resonance    
  std::vector<std::pair<double, double> > NRparams;
  if(mctruthmode==125){
//    NRparams.push_back(std::make_pair<double,double>(5, -1));  //vary the width
//    NRparams.push_back(std::make_pair<double,double>(8, -1));
//    NRparams.push_back(std::make_pair<double,double>(10,-1));
//    NRparams.push_back(std::make_pair<double,double>(11,-1));
//    NRparams.push_back(std::make_pair<double,double>(12,-1));
//    NRparams.push_back(std::make_pair<double,double>(13,-1));
//    NRparams.push_back(std::make_pair<double,double>(14,-1));
//    NRparams.push_back(std::make_pair<double,double>(15,-1));
//    NRparams.push_back(std::make_pair<double,double>(16,-1));
//    NRparams.push_back(std::make_pair<double,double>(17,-1));
//    NRparams.push_back(std::make_pair<double,double>(18,-1));
//    NRparams.push_back(std::make_pair<double,double>(19,-1));
//    NRparams.push_back(std::make_pair<double,double>(20,-1));
//    NRparams.push_back(std::make_pair<double,double>(22,-1));
//    NRparams.push_back(std::make_pair<double,double>(25,-1));
//    NRparams.push_back(std::make_pair<double,double>(30,-1));
  }else if(suffix=="" && (isMC_GG || isMC_VBF)){ //consider the other points only when no suffix is being used    
      NRparams.push_back(std::make_pair<double,double>(1.0, 0.0) ); //cp, brnew
      NRparams.push_back(std::make_pair<double,double>(0.6, 0.0) ); //cp, brnew
      NRparams.push_back(std::make_pair<double,double>(0.3, 0.0) ); //cp, brnew
      NRparams.push_back(std::make_pair<double,double>(0.1, 0.0) ); //cp, brnew
  }
  if(NRparams.size()<=0)NRparams.push_back(std::make_pair<double,double>(-1.0, -1.0)); //no reweighting

  std::vector<TString>    NRsuffix; 
  std::vector<double[6] > lShapeWeights(NRparams.size());   //WEIGHT for LineShape (NNLO kFactors + Interf), split in shape and scale unc with the following format:  scaleNominal, shapeNominal, scaleDown, shapeDown, scaleUp, shapeUp
  for(unsigned int nri=0;nri<NRparams.size();nri++){
     char tmp[255];
     if(      NRparams[nri].first<0 && NRparams[nri].second<0){   sprintf(tmp,"%s", "");
     }else if(NRparams[nri].first>0 && NRparams[nri].second<0){   sprintf(tmp,"_width%6.2f"      , NRparams[nri].first);
     }else{                                                       sprintf(tmp,"_cp%3.2f_brn%3.2f",NRparams[nri].first, NRparams[nri].second); 
     } 
     NRsuffix.push_back(TString(tmp));
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
  

  //ELECTROWEAK CORRECTION WEIGHTS
  std::vector<std::vector<float>> ewkTable, ZZ_NNLOTable;
  if(isMC_ZZ2l2nu){
  	ewkTable = EwkCorrections::readFile_and_loadEwkTable(dtag);
		ZZ_NNLOTable = ZZatNNLO::readFile_and_loadTable(dtag);
	}

  //#######################################
  //####      LINE SHAPE WEIGHTS       ####
  //#######################################
  std::map<std::pair<double,double>, std::vector<std::pair<double, TGraph *> > > hLineShapeGrVec;  
  if(isMC_GG || isMC_VBF){
      TH1D* hGen=new TH1D("hGen", "hGen", 1000, 0, 4000);
      utils::getHiggsLineshapeFromMiniAOD(urls, hGen);
      printf("hGen integral = %f\n", hGen->Integral());

      TGraph* hLineShapeNominal= new TGraph(hGen);
      TFile* nrLineShapesFile=NULL;
      if(mctruthmode==125){
	TString nrLineShapesFileUrl(weightsDir+"/higgs_width_zz2l2nu.root");
	gSystem->ExpandPathName(nrLineShapesFileUrl);
	nrLineShapesFile=TFile::Open(nrLineShapesFileUrl);
      }else{
	TString nrLineShapesFileUrl(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/weights/NR_weightsRun2.root");
	gSystem->ExpandPathName(nrLineShapesFileUrl);
	nrLineShapesFile=TFile::Open(nrLineShapesFileUrl);
      }


      bool isGGZZContinuum = (dtag.Contains("ggZZ"));

      //loop over possible scenarios: SM or BSM
      for(size_t nri=0; nri<NRparams.size(); nri++){
         //recompute weights depending on the scenario (SM or BSM)
         TGraph* shapeWgtsGr      = new TGraph; shapeWgtsGr->SetName("shapeWgts_"+ NRsuffix[nri]);          float shapeNorm(0);     double signalNorm(1.0);
         TGraph* shapeWgts_upGr   = new TGraph; shapeWgts_upGr->SetName("shapeWgtsUp_"+ NRsuffix[nri]);     float shapeUpNorm(0);   double signalUpNorm(1.0);
         TGraph* shapeWgts_downGr = new TGraph; shapeWgts_downGr->SetName("shapeWgtsDown_"+ NRsuffix[nri]); float shapeDownNorm(0); double signalDownNorm(1.0);

         TGraph* nrWgtGr=NULL; TGraph* nrWgtUpGr=NULL; TGraph* nrWgtDownGr=NULL;
         if(isGGZZContinuum){
            nrWgtGr     = higgs::utils::weightGGZZContinuum(nrLineShapesFile,signalNorm    ,""  );
            nrWgtUpGr   = higgs::utils::weightGGZZContinuum(nrLineShapesFile,signalUpNorm  ,"Up");         
            nrWgtDownGr = higgs::utils::weightGGZZContinuum(nrLineShapesFile,signalDownNorm,"Dn");        
         }else if(mctruthmode!=125 && NRparams[nri].first>=0){
            nrWgtGr     = higgs::utils::weightNarrowResonnance(isMC_VBF,HiggsMass, NRparams[nri].first, NRparams[nri].second, nrLineShapesFile,signalNorm    ,""     );
            nrWgtUpGr   = higgs::utils::weightNarrowResonnance(isMC_VBF,HiggsMass, NRparams[nri].first, NRparams[nri].second, nrLineShapesFile,signalUpNorm  ,"_up"  );          
            nrWgtDownGr = higgs::utils::weightNarrowResonnance(isMC_VBF,HiggsMass, NRparams[nri].first, NRparams[nri].second, nrLineShapesFile,signalDownNorm,"_down"); 
         }
         if(!nrWgtUpGr){nrWgtUpGr=nrWgtGr;  signalUpNorm=signalNorm;}
         if(!nrWgtDownGr){nrWgtDownGr=nrWgtGr;  signalDownNorm=signalNorm;}

         double hySum=0;
         for(int ip=1; ip<=hGen->GetXaxis()->GetNbins(); ip++){
    	     Double_t hmass    = hGen->GetBinCenter(ip);
	     Double_t hy       = hGen->GetBinContent(ip);
             hySum            += hy;
		  
	     Double_t shapeWgt(1.0),shapeWgtUp(1.0),shapeWgtDown(1.0);
             if(mctruthmode==125){  //reweighting to SM higgs width various width 
                TString var("");
   	        if(dtag.Contains("ScaleUp"))   var="up";
	        if(dtag.Contains("ScaleDown")) var="down";
  	        shapeWgt       = higgs::utils::weightToH125Interference(hmass,NRparams[nri].first,nrLineShapesFile,var);
	        shapeWgtUp     = shapeWgt;
	        shapeWgtDown   = shapeWgt;
 	     }else if(NRparams[nri].first>=0 || isGGZZContinuum){ //reweighting to Narrow reasonnance or EWS model, or to ggZZ continuum
                shapeWgt       = nrWgtGr    ->Eval(hmass);
                shapeWgtUp     = nrWgtUpGr  ->Eval(hmass);
                shapeWgtDown   = nrWgtDownGr->Eval(hmass);                       
             }else{ //unknown case, do not reweight  (SM-Like)
     	        shapeWgt     = 1.0;
	        shapeWgtUp   = 1.0;
 	        shapeWgtDown = 1.0;
  	     }
	
             //if(ip==150)printf("weight for mZZ %f = %f x %f\n", hmass, shapeWgt, hy);
	     shapeWgtsGr     ->SetPoint(shapeWgtsGr     ->GetN(), hmass, shapeWgt);       shapeNorm     += shapeWgt*hy;
	     shapeWgts_upGr  ->SetPoint(shapeWgts_upGr  ->GetN(), hmass, shapeWgtUp);     shapeUpNorm   += shapeWgtUp*hy;
	     shapeWgts_downGr->SetPoint(shapeWgts_downGr->GetN(), hmass, shapeWgtDown);   shapeDownNorm += shapeWgtDown*hy;
         }
	      
         if(mctruthmode!=125){
	    if(hySum>0){
	       shapeNorm     /= hySum;
	       shapeUpNorm   /= hySum;
	       shapeDownNorm /= hySum;
     	    }

            //fix possible normalization issues
            printf("C'=%6.2f  BRNew=%6.2f shapeNorm = %f Up=%f Down=%f  signalNorm=%f Up=%f Down=%f\n", NRparams[nri].first, NRparams[nri].second, shapeNorm, shapeUpNorm, shapeDownNorm, signalNorm, signalUpNorm, signalDownNorm);
 	    for(Int_t ip=0; ip<shapeWgtsGr->GetN(); ip++){
  	         Double_t x,y;
	         shapeWgtsGr->GetPoint(ip,x,y);
	         shapeWgtsGr->SetPoint(ip,x,y/shapeNorm);   
 	         shapeWgts_upGr->GetPoint(ip,x,y);
	         shapeWgts_upGr->SetPoint(ip,x,y/shapeUpNorm);		    
  	         shapeWgts_downGr->GetPoint(ip,x,y);
	         shapeWgts_downGr->SetPoint(ip,x,y/shapeDownNorm);
	    }
         }

         //all done here...
	 std::vector<std::pair<double, TGraph*> > inrWgts = {std::make_pair(signalNorm/signalNorm, shapeWgtsGr), std::make_pair(signalUpNorm/signalNorm,shapeWgts_upGr), std::make_pair(signalDownNorm/signalNorm,shapeWgts_downGr)};
	 hLineShapeGrVec[ NRparams[nri] ] = inrWgts;
      }	  
      if(nrLineShapesFile){nrLineShapesFile->Close(); delete nrLineShapesFile;}
  }


  //##############################################
  //########    INITIATING HISTOGRAMS     ########
  //##############################################
  SmartSelectionMonitor mon;
  printf("Definition of plots");

  //generator level control : add an underflow entry to make sure the histo is kept
  ((TH1F*)mon.addHistogram( new TH1F( "higgsMass_raw",     ";Higgs Mass [GeV];Events", 500,0,1500) ))->Fill(-1.0,0.0001);
  for(unsigned int nri=0;nri<NRparams.size();nri++){ 
    ((TH1F*)mon.addHistogram( new TH1F( "higgsMass_shape"+NRsuffix[nri] , ";Higgs Mass;Events [GeV]", 500,0,1500) ))->Fill(-1.0,0.0001);
    ((TH1F*)mon.addHistogram( new TH1F( "higgsMass_shape&scale"+NRsuffix[nri] , ";Higgs Mass;Events [GeV]", 500,0,1500) ))->Fill(-1.0,0.0001);
  }

  mon.addHistogram( new TH1F( "wdecays",     ";W decay channel",5,0,5) );
  mon.addHistogram( new TH1F( "zdecays",     ";Z decay channel",6,0,6) );

  mon.addHistogram( new TH1F( "metFilter_eventflow",     ";metEventflow",20,0,20) );

  //event selection
  TH1F *h=(TH1F*) mon.addHistogram( new TH1F ("eventflow", ";;Events", 10,0,10) );
  h->GetXaxis()->SetBinLabel(1,"raw");
  h->GetXaxis()->SetBinLabel(2,"#geq 2 iso leptons");
  h->GetXaxis()->SetBinLabel(3,"|M-91|<15");
  h->GetXaxis()->SetBinLabel(4,"p_{T}>55");
  h->GetXaxis()->SetBinLabel(5,"3^{rd}-lepton veto");
  h->GetXaxis()->SetBinLabel(6,"b-veto"); 
  h->GetXaxis()->SetBinLabel(7,"#Delta #phi(jet,E_{T}^{miss})>0.5");
  h->GetXaxis()->SetBinLabel(8,"E_{T}^{miss}>80");
  h->GetXaxis()->SetBinLabel(9,"E_{T}^{miss}>125");


  h=(TH1F*) mon.addHistogram( new TH1F ("trigger", ";;Events", 10,0,10) );
  h->GetXaxis()->SetBinLabel(1,"#mu#mu");
  h->GetXaxis()->SetBinLabel(2,"#mu");
  h->GetXaxis()->SetBinLabel(3,"ee");
  h->GetXaxis()->SetBinLabel(4,"e");
  h->GetXaxis()->SetBinLabel(5,"e#mu");
  h->GetXaxis()->SetBinLabel(6,"#gamma"); 

  //pu control
  mon.addHistogram( new TH1F( "nvtx",";Vertices;Events",50,0,50) ); 
  mon.addHistogram( new TH1F( "nvtxraw",";Vertices;Events",50,0,50) ); 
  mon.addHistogram( new TH1F( "rho",";#rho;Events",50,0,25) ); 

  // photon control
  mon.addHistogram(new TH1F("npho",   ";Number of Photons;Events", 20, 0, 20) ); 
  mon.addHistogram(new TH1F("npho55", ";Number of Photons;Events", 20, 0, 20) ); 
  mon.addHistogram(new TH1F("npho100",";Number of Photons;Events", 20, 0, 20) ); 
  mon.addHistogram(new TH1F("photonpt", ";Photon pT [GeV];Events", 500, 0, 1000) ); 
  mon.addHistogram(new TH1F("phopt", ";Photon pT [GeV];Events", 500, 0, 1000) ); 
  mon.addHistogram(new TH1F("phoeta", ";Photon pseudo-rapidity;Events", 50, 0, 5) );
  mon.addHistogram(new TH1F("bosonnvtx", ";Photon #eta;Events", 50, 0, 50) );
  mon.addHistogram(new TH1F("bosoneta", ";Photon #eta;Events", 100, -5, 5) );
  mon.addHistogram(new TH1F("bosonphi", ";Photon #phi;Events", 80, -4, 4) );
  mon.addHistogram(new TH1F("bosonphiHG", ";Photon #phi;Events", 800, -4, 4) );
  mon.addHistogram(new TH1F("metphi", ";MET #phi;Events", 80, -4, 4) );
  mon.addHistogram(new TH1F("metphiUnCor", ";MET #phi;Events", 80, -4, 4) );
  mon.addHistogram(new TH1F("dphi_boson_met", ";#Delta #phi(#gamma,MET);Events", 40, 0, 4) );
  
  //lepton control
  mon.addHistogram( new TH1F( "nleptons",   ";Nleptons;Events",10,0,10) );
  mon.addHistogram( new TH1F( "leadpt",     ";Transverse momentum [GeV];Events", 50,0,500) );
  mon.addHistogram( new TH1F( "leadeta",    ";Pseudo-rapidity;Events", 50,0,2.6) );
  mon.addHistogram( new TH1F( "trailerpt",  ";Transverse momentum [GeV];Events", 50,0,500) );
  mon.addHistogram( new TH1F( "trailereta", ";Pseudo-rapidity;Events", 50,0,2.6) );
  mon.addHistogram( new TH1F( "zy",         ";Rapidity;Events", 50,0,3) );
  mon.addHistogram( new TH1F( "zmass",      ";Mass [GeV];Events / 2 GeV", 100,40,240) );
  mon.addHistogram( new TH1F( "zmass_btag50", ";Mass [GeV];Events / 2 GeV", 100,40,200) );
  mon.addHistogram( new TH1F( "zmass_bveto50",";Mass [GeV];Events / 2 GeV", 100,40,200) );
  mon.addHistogram( new TH1F( "zmass_btag80", ";Mass [GeV];Events / 2 GeV", 100,40,200) );
  mon.addHistogram( new TH1F( "zmass_bveto80",";Mass [GeV];Events / 2 GeV", 100,40,200) );
  mon.addHistogram( new TH1F( "zmass_btag125", ";Mass [GeV];Events / 2 GeV", 100,40,200) );
  mon.addHistogram( new TH1F( "zmass_bveto125",";Mass [GeV];Events / 2 GeV", 100,40,200) );
  mon.addHistogram( new TH1F( "zpt",        ";Transverse momentum [GeV];Events",100,0,1500));
  Double_t zptaxis[]= {0,15,30,45,60,75,90,105,120,135,150,165,180,195,210,225,240,255,270,285,300,315,330,345,360,375,390,405,435,465,495,525,555,585,615,675,735,795,855,975};
  Int_t nzptAxis=sizeof(zptaxis)/sizeof(Double_t);
  mon.addHistogram( new TH1F( "zpt_rebin",  ";Transverse momentum [GeV];Events / GeV", nzptAxis-1,zptaxis));
  mon.addHistogram( new TH1F( "zptMet125",        ";Transverse momentum [GeV];Events",100,0,1500));
  mon.addHistogram( new TH1F( "qmass",      ";Mass [GeV];Events / (1 GeV)",100,76,106));
  mon.addHistogram( new TH1F( "qt",         ";Transverse momentum [GeV];Events / GeV",1500,0,1500));
  mon.addHistogram( new TH1F( "qtraw",      ";Transverse momentum [GeV];Events / GeV",1500,0,1500));

  //extra leptons in the event
  mon.addHistogram( new TH1F( "nextraleptons", ";Extra leptons;Events",4,0,4) );
  mon.addHistogram( new TH1F( "thirdleptonpt", ";Transverse momentum;Events", 50,0,500) );
  mon.addHistogram( new TH1F( "thirdleptoneta", ";Pseudo-rapidity;Events", 50,0,2.6) );
  mon.addHistogram( new TH1F( "thirdleptonmt", ";Transverse mass(3^{rd} lepton,E_{T}^{miss}) [GeV];Events", 50,0,500) );


  mon.addHistogram( new TH1F("csv",      ";Combined Secondary Vertex;Jets",50,0.,1.) );
  mon.addHistogram( new TH1F("csvb",     ";Combined Secondary Vertex;Jets",50,0.,1.) );
  mon.addHistogram( new TH1F("csvc",     ";Combined Secondary Vertex;Jets",50,0.,1.) );
  mon.addHistogram( new TH1F("csvothers",";Combined Secondary Vertex;Jets",50,0.,1.) );
  mon.addHistogram( new TH1F("leadjetpt",    ";Transverse momentum [GeV];Events",50,0,1000) );
  mon.addHistogram( new TH1F("trailerjetpt", ";Transverse momentum [GeV];Events",50,0,1000) );
  mon.addHistogram( new TH1F("vbfjeteta",    ";Pseudo-rapidity;Events",25,0,5) );
  mon.addHistogram( new TH1F("fwdjeteta",    ";Pseudo-rapidity;Events",25,0,5) );
  mon.addHistogram( new TH1F("cenjeteta",       ";Pseudo-rapidity;Events",25,0,5) );
  mon.addHistogram( new TH1F("jetId",       ";Pseudo-rapidity;Events",25,0,5) );
  Double_t mjjaxis[32];
  mjjaxis[0]=0.01;
  for(size_t i=1; i<20; i++)  mjjaxis[i]   =50*i;        //0-1000
  for(size_t i=0; i<5; i++)   mjjaxis[20+i]=1000+100*i; //1000-1500
  for(size_t i=0; i<=5; i++)   mjjaxis[25+i]=1500+300*i; //1500-5000  
  mjjaxis[31]=5000;
  mon.addHistogram( new TH1F("vbfmjj"       , ";Dijet invariant mass [GeV];Events / GeV",31,mjjaxis) );
  mon.addHistogram( new TH1F("vbfdphijj"    , ";Azimuthal angle difference;Events",20,0,3.5) );
  mon.addHistogram( new TH1F("vbfdetajj"    , ";Pseudo-rapidity span;Events",20,0,10) );
  mon.addHistogram( new TH1F("vbfcjv"       , ";Central jet multiplicity;Events",5,0,5) );
  TH1F* hjetsfinal   = (TH1F*) mon.addHistogram( new TH1F("njetsfinal",   ";Jet multiplicity;Events",5,0,5) );
  TH1F* hjets   = (TH1F*) mon.addHistogram( new TH1F("njets",   ";Jet multiplicity;Events",5,0,5) );
  TH1F* hbtags  = (TH1F*) mon.addHistogram( new TH1F("nbtags",  ";b-tag multiplicity;Events",5,0,5) );
  for(int ibin=1; ibin<=hjets->GetXaxis()->GetNbins(); ibin++){
      TString label("");
      if(ibin==h->GetXaxis()->GetNbins()) label +="#geq";
      label += (ibin-1);
      hjets   ->GetXaxis()->SetBinLabel(ibin,label);
      hjetsfinal->GetXaxis()->SetBinLabel(ibin,label);      
      hbtags  ->GetXaxis()->SetBinLabel(ibin,label);
  } 

  mon.addHistogram( new TH1F( "mindphijmet",  ";min #Delta#phi(jet,E_{T}^{miss});Events",20,0,4) );
  mon.addHistogram( new TH1F( "mindphijmet25",  ";min #Delta#phi(jet,E_{T}^{miss});Events",20,0,4) );
  mon.addHistogram( new TH1F( "mindphijmet50",  ";min #Delta#phi(jet,E_{T}^{miss});Events",20,0,4) );
  mon.addHistogram( new TH1F( "mindphijmetNM1",  ";min #Delta#phi(jet,E_{T}^{miss});Events",20,0,4) );
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
  mon.addHistogram( new TH1F( "mtcheckpoint"  ,         ";Transverse mass [GeV];Events / GeV",160,150,1750) );
  mon.addHistogram( new TH1F( "metcheckpoint" ,         ";Missing transverse energy [GeV];Events / GeV",100,0,500) );

  mon.addHistogram( new TH1F( "met_Inbtag",          ";Missing transverse energy [GeV];Events / GeV",nmetAxis-1,metaxis) ); //50,0,1000) );
  mon.addHistogram( new TH1F( "met_Inbveto",          ";Missing transverse energy [GeV];Events / GeV",nmetAxis-1,metaxis) ); //50,0,1000) );
  mon.addHistogram( new TH1F( "met_Outbtag",          ";Missing transverse energy [GeV];Events / GeV",nmetAxis-1,metaxis) ); //50,0,1000) );
  mon.addHistogram( new TH1F( "met_Outbveto",          ";Missing transverse energy [GeV];Events / GeV",nmetAxis-1,metaxis) ); //50,0,1000) );

  for(size_t ivar=0; ivar<nvarsToInclude; ivar++){
     mon.addHistogram( new TH1F( (TString("metSyst")+varNames[ivar]),          ";Missing transverse energy [GeV];Events / GeV",nmetAxis-1,metaxis) ); //50,0,1000) );
     mon.addHistogram( new TH1F( (TString("mtSyst")+varNames[ivar]),         ";Transverse mass [GeV];Events / GeV",nmtAxis-1,mtaxis) );
  }


  mon.addHistogram( new TH1F( "mt_Inbtag50"  ,         ";Transverse mass [GeV];Events / GeV",nmtAxis-1,mtaxis) );
  mon.addHistogram( new TH1F( "mt_Inbveto50"  ,         ";Transverse mass [GeV];Events / GeV",nmtAxis-1,mtaxis) );
  mon.addHistogram( new TH1F( "mt_Inbtag80"  ,         ";Transverse mass [GeV];Events / GeV",nmtAxis-1,mtaxis) );
  mon.addHistogram( new TH1F( "mt_Inbveto80"  ,         ";Transverse mass [GeV];Events / GeV",nmtAxis-1,mtaxis) );
  mon.addHistogram( new TH1F( "mt_Inbtag125"  ,         ";Transverse mass [GeV];Events / GeV",nmtAxis-1,mtaxis) );
  mon.addHistogram( new TH1F( "mt_Inbveto125"  ,         ";Transverse mass [GeV];Events / GeV",nmtAxis-1,mtaxis) );
  mon.addHistogram( new TH1F( "mt_Outbtag50"  ,         ";Transverse mass [GeV];Events / GeV",nmtAxis-1,mtaxis) );
  mon.addHistogram( new TH1F( "mt_Outbveto50"  ,         ";Transverse mass [GeV];Events / GeV",nmtAxis-1,mtaxis) );
  mon.addHistogram( new TH1F( "mt_Outbtag80"  ,         ";Transverse mass [GeV];Events / GeV",nmtAxis-1,mtaxis) );
  mon.addHistogram( new TH1F( "mt_Outbveto80"  ,         ";Transverse mass [GeV];Events / GeV",nmtAxis-1,mtaxis) );
  mon.addHistogram( new TH1F( "mt_Outbtag125"  ,         ";Transverse mass [GeV];Events / GeV",nmtAxis-1,mtaxis) );
  mon.addHistogram( new TH1F( "mt_Outbveto125"  ,         ";Transverse mass [GeV];Events / GeV",nmtAxis-1,mtaxis) );


  mon.addHistogram( new TH1F( "mtfinal"  ,         ";Transverse mass [GeV];Events / GeV",nmtAxis-1,mtaxis) );
  mon.addHistogram( new TH1F( "metfinal",          ";Missing transverse energy [GeV];Events / GeV",nmetAxis-1,metaxis) ); //50,0,1000) );
  mon.addHistogram( new TH1F( "mindphijmetfinal",  ";min #Delta#phi(jet,E_{T}^{miss});Events",20,0,4) );
  mon.addHistogram( new TH1F("vbfmjjfinal"       , ";Dijet invariant mass [GeV];Events / GeV",31,mjjaxis) );
  mon.addHistogram( new TH1F("vbfdetajjfinal"    , ";Pseudo-rapidity span;Events",20,0,10) );
  
  mon.addHistogram( new TH1F( "mzz",   ";M_{ZZ} [GeV];Events / GeV", 150, 0, 1500) ); //The binning is the same than the one for the corrections.


  //
  // HISTOGRAMS FOR OPTIMIZATION and STATISTICAL ANALYSIS
  //
  //

  //NEED FOR ALL OPTIMIZATION
  TH1F* Hoptim_systs     =  (TH1F*) mon.addHistogram( new TH1F ("optim_systs"    , ";syst;;", nvarsToInclude,0,nvarsToInclude) ) ;
  for(size_t ivar=0; ivar<nvarsToInclude; ivar++){
      Hoptim_systs->GetXaxis()->SetBinLabel(ivar+1, varNames[ivar]);
  }

  std::vector<double> optim_Cuts1_met;
  for(double met=50;met<140;met+=5) {  optim_Cuts1_met    .push_back(met);  }
  TH2F* Hoptim_cuts  =(TH2F*)mon.addHistogram(new TProfile2D("optim_cut",      ";cut index;variable",       optim_Cuts1_met.size(),0,optim_Cuts1_met.size(), 1, 0, 1)) ;
  Hoptim_cuts->GetYaxis()->SetBinLabel(1, "met>");
  for(unsigned int index=0;index<optim_Cuts1_met.size();index++){ Hoptim_cuts    ->Fill(index, 0.0, optim_Cuts1_met[index]);  }
  for(size_t ivar=0; ivar<nvarsToInclude; ivar++){
      for(unsigned int nri=0;nri<NRparams.size();nri++){ 
	mon.addHistogram( new TH2F (TString("mt_shapes")+NRsuffix[nri]+varNames[ivar],";cut index;Transverse mass [GeV];Events",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(), 160,150,1750) );     
	mon.addHistogram( new TH2F (TString("met_shapes")+NRsuffix[nri]+varNames[ivar],";cut index;Missing transverse energy [GeV];Events",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),100 ,0,500) );     
	TH2F *h=(TH2F *) mon.addHistogram( new TH2F ("mt_shapes_NRBctrl"+NRsuffix[nri]+varNames[ivar],";cut index;Selection region;Events",optim_Cuts1_met.size(),0,optim_Cuts1_met.size(),6,0,6) );
	h->GetYaxis()->SetBinLabel(1,"M_{in}^{ll}/=0 b-tags");
	h->GetYaxis()->SetBinLabel(2,"M_{out}^{ll}/=0 b-tags");
	h->GetYaxis()->SetBinLabel(3,"M_{out+}^{ll}/=0 b-tags");
	h->GetYaxis()->SetBinLabel(4,"M_{in}^{ll}/#geq 1 b-tag");
	h->GetYaxis()->SetBinLabel(5,"M_{out}^{ll}/#geq 1 b-tag");
	h->GetYaxis()->SetBinLabel(6,"M_{out+}^{ll}/#geq 1 b-tag");
      }
    }

  std::vector<std::vector<double>> optim_Cuts_VBF; 
  for(double jet2Pt=20    ;jet2Pt<50;jet2Pt+=5) { 
     for(double jet1Pt=jet2Pt;jet1Pt<50;jet1Pt+=5) {
        for(double deta=2.0     ;deta<5.0;deta+=0.25) { 
           for(double mjj=100.0    ;mjj<1000;mjj+=50) {
              optim_Cuts_VBF.push_back( std::vector<double>{jet1Pt, jet2Pt, deta, mjj} );
  }}}}
   
  TH2F* Hoptim_cuts_VBF  =(TH2F*)mon.addHistogram(new TProfile2D("optim_cut_VBF",      ";cut index;variable",       optim_Cuts_VBF.size(),0,optim_Cuts_VBF.size(), 4, 0, 4)) ;
  Hoptim_cuts_VBF->GetYaxis()->SetBinLabel(1, "jet1 p_{T}>");
  Hoptim_cuts_VBF->GetYaxis()->SetBinLabel(2, "jet2 p_{T}>");
  Hoptim_cuts_VBF->GetYaxis()->SetBinLabel(3, "d#eta>");
  Hoptim_cuts_VBF->GetYaxis()->SetBinLabel(4, "M_{jj}>");
  for(unsigned int index=0;index<optim_Cuts_VBF.size();index++){ for(unsigned int cut=0;cut<4;cut++){ Hoptim_cuts_VBF    ->Fill(index, float(cut), optim_Cuts_VBF[index][cut]); }  }
  for(size_t ivar=0; ivar<nvarsToInclude; ivar++){
     if(ivar>0)continue;// do not fill for systematics in order to save time
     mon.addHistogram( new TH2F (TString("vbf_shapes")+varNames[ivar],";cut index;Transverse mass [GeV];Events",optim_Cuts_VBF.size(),0,optim_Cuts_VBF.size(), 1, 0, 2500) );     
  }

     
  //##############################################
  //######## GET READY FOR THE EVENT LOOP ########
  //##############################################

  //MET CORRection level
  pat::MET::METCorrectionLevel metcor = pat::MET::METCorrectionLevel::Type1XY;
  //jet energy scale and uncertainties 
  TString jecDir = runProcess.getParameter<std::string>("jecDir");
  gSystem->ExpandPathName(jecDir);
  FactorizedJetCorrector *jesCor = NULL;
  if(isMC || is2015data) jesCor        = utils::cmssw::getJetCorrector(jecDir,isMC);
  TString pf(isMC ? "MC" : "DATA");
  JetCorrectionUncertainty *totalJESUnc = NULL;
  if(isMC || is2015data) totalJESUnc = new JetCorrectionUncertainty((jecDir+"/"+pf+"_Uncertainty_AK4PFchs.txt").Data());
  //muon energy scale and uncertainties
  TString muscleDir = runProcess.getParameter<std::string>("muscleDir");
  gSystem->ExpandPathName(muscleDir);
  rochcor2015* muCor = NULL; //need to be updated for 2016
  if(isMC || is2015data) muCor = new rochcor2015();  //replace the MuScleFitCorrector we used at run1
  //photon and electron enerhy scale based on https://twiki.cern.ch/twiki/bin/viewauth/CMS/EGMSmearer    (adapted to the miniAOD/FWLite framework) 

	ElectronEnergyCalibratorRun2 ElectronEnCorrector;
	string EGammaEnergyCorrectionFile = "";
	PhotonEnergyCalibratorRun2 PhotonEnCorrector;
  	EpCombinationTool theEpCombinationTool;
	if(isMC || is2015data){
  	EGammaEnergyCorrectionFile = "EgammaAnalysis/ElectronTools/data/76X_16DecRereco_2015"; 
  	PhotonEnCorrector = PhotonEnergyCalibratorRun2(isMC, false, EGammaEnergyCorrectionFile);
  	PhotonEnCorrector.initPrivateRng(new TRandom(1234));
  	theEpCombinationTool.init((string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/weights/GBRForest_data_25ns.root").c_str(), "gedelectron_p4combination_25ns");  //got confirmation from Matteo Sani that this works for both data and MC 
  	ElectronEnCorrector = ElectronEnergyCalibratorRun2(theEpCombinationTool, isMC, false, EGammaEnergyCorrectionFile);
  	ElectronEnCorrector.initPrivateRng(new TRandom(1234));
	}
  //lepton efficiencies
  LeptonEfficiencySF lepEff;

  //b-tagging: beff and leff must be derived from the MC sample using the discriminator vs flavor
  //the scale factors are taken as average numbers from the pT dependent curves see:
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2012_Data_and_MC_EPS13_prescript
  BTagSFUtil btsfutil;
  float beff(0.68), sfb(0.99), sfbunc(0.015);
  float leff(0.13), sfl(1.05), sflunc(0.12);

  double btagLoose = 0.605; //https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation74X FIXME, I sent an email to Petra to know more (Hugo)
  // setup calibration readers
  BTagCalibration btagCalib("CSVv2", string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/weights/btagSF_CSVv2.csv");
  BTagCalibrationReader btagCal   (&btagCalib, BTagEntry::OP_LOOSE, "mujets", "central");  // calibration instance, operating point, measurement type, systematics type
  BTagCalibrationReader btagCalUp (&btagCalib, BTagEntry::OP_LOOSE, "mujets", "up"     );  // sys up
  BTagCalibrationReader btagCalDn (&btagCalib, BTagEntry::OP_LOOSE, "mujets", "down"   );  // sys down
  BTagCalibrationReader btagCalL  (&btagCalib, BTagEntry::OP_LOOSE, "comb", "central");  // calibration instance, operating point, measurement type, systematics type
  BTagCalibrationReader btagCalLUp(&btagCalib, BTagEntry::OP_LOOSE, "comb", "up"     );  // sys up
  BTagCalibrationReader btagCalLDn(&btagCalib, BTagEntry::OP_LOOSE, "comb", "down"   );  // sys down

  // from Btag SF and eff from https://indico.cern.ch/event/437675/#preview:1629681
  beff = 0.747; sfb = 0.899; //for Loose WP  //sfb is not actually used as it's taken from btagCal

  //pileup weighting
  edm::LumiReWeighting* LumiWeights = NULL;
  utils::cmssw::PuShifter_t PuShifters;
  double PUNorm[] = {1,1,1};

  //MC normalization (to 1/pb)
  double xsecWeight = 1.0;
  if(isMC){
          std::vector<double> dataPileupDistributionDouble = runProcess.getParameter< std::vector<double> >("datapileup");
          std::vector<float> dataPileupDistribution; for(unsigned int i=0;i<dataPileupDistributionDouble.size();i++){dataPileupDistribution.push_back(dataPileupDistributionDouble[i]);}
          std::vector<float> mcPileupDistribution;

          double totalNumEvent = utils::getMCPileupDistributionAndTotalEventFromMiniAOD(urls,dataPileupDistribution.size(), mcPileupDistribution);
          xsecWeight=xsec/totalNumEvent;

	  //utils::getMCPileupDistributionFromMiniAOD(urls,dataPileupDistribution.size(), mcPileupDistribution);
          while(mcPileupDistribution.size()<dataPileupDistribution.size())  mcPileupDistribution.push_back(0.0);
          while(mcPileupDistribution.size()>dataPileupDistribution.size())dataPileupDistribution.push_back(0.0);
          gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE
          LumiWeights = new edm::LumiReWeighting(mcPileupDistribution,dataPileupDistribution);
          PuShifters=utils::cmssw::getPUshifters(dataPileupDistribution,0.05);
          utils::getPileupNormalization(mcPileupDistribution, PUNorm, LumiWeights, PuShifters);
  }
 
  gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE


  higgs::utils::EventCategory eventCategoryInst(higgs::utils::EventCategory::EXCLUSIVE2JETSVBF); //jet(0,>=1)+vbf binning
  patUtils::MetFilter metFiler;
  if(!isMC){ 
		if(is2015data){
			metFiler.FillBadEvents(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/MetFilter/DoubleEG_RunD/DoubleEG_csc2015.txt");
			metFiler.FillBadEvents(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/MetFilter/DoubleEG_RunD/DoubleEG_ecalscn1043093.txt"); 
			metFiler.FillBadEvents(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/MetFilter/DoubleMuon_RunD/DoubleMuon_csc2015.txt");
			metFiler.FillBadEvents(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/MetFilter/DoubleMuon_RunD/DoubleMuon_ecalscn1043093.txt"); 
			metFiler.FillBadEvents(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/MetFilter/MuonEG_RunD/MuonEG_csc2015.txt");
			metFiler.FillBadEvents(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/MetFilter/MuonEG_RunD/MuonEG_ecalscn1043093.txt"); 
		}
		else if(is2016data){

		}
        //FIXME, we need to add here the single mu, single el, and gamma path
  }
  string debugText = "";

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
         double puWeightUp = 1.0;
         double puWeightDown = 1.0;
         float puWeight(1.0);

         //##############################################   EVENT LOOP STARTS   ##############################################
         //if(!isMC && duplicatesChecker.isDuplicate( ev.run, ev.lumi, ev.event) ) { nDuplicates++; continue; }

         //Skip bad lumi
         if(!isMC && !goodLumiFilter.isGoodLumi(ev.eventAuxiliary().run(),ev.eventAuxiliary().luminosityBlock()))continue;

         reco::GenParticleCollection gen;
         GenEventInfoProduct eventInfo;       
          if(isMC){
            fwlite::Handle< reco::GenParticleCollection > genHandle;
            genHandle.getByLabel(ev, "prunedGenParticles");
            if(genHandle.isValid()){ gen = *genHandle;}

            fwlite::Handle< GenEventInfoProduct > genEventInfoHandle;
            genEventInfoHandle.getByLabel(ev, "generator");        
            if(genEventInfoHandle.isValid()){ eventInfo = *genEventInfoHandle;}

            //WEIGHT for NLO negative interference
            //totalNumEvent+=eventInfo.weight();
            weight *= eventInfo.weight(); 
       

            //WEIGHT for Pileup
	    int ngenITpu = 0;
	    fwlite::Handle< std::vector<PileupSummaryInfo> > puInfoH;
            puInfoH.getByLabel(ev, "slimmedAddPileupInfo");
            for(std::vector<PileupSummaryInfo>::const_iterator it = puInfoH->begin(); it != puInfoH->end(); it++){
               if(it->getBunchCrossing()==0)      { ngenITpu += it->getTrueNumInteractions(); } //getPU_NumInteractions(); 
            }
            puWeight          = LumiWeights->weight(ngenITpu) * PUNorm[0];
            puWeightUp  = PuShifters[utils::cmssw::PUUP  ]->Eval(ngenITpu) * (PUNorm[2]/PUNorm[0]);
            puWeightDown = PuShifters[utils::cmssw::PUDOWN]->Eval(ngenITpu) * (PUNorm[1]/PUNorm[0]);
            weight *= puWeight;

            if(isMC && (mctruthmode==15 || mctruthmode==1113)){// && (string(dtag.Data()).find("Z#rightarrow")==0 || isMC_ZZ2l2nu))
                int prodId = 1;
                for( unsigned int k=0; k<gen.size(); ++k){	
                        if( gen[k].isHardProcess() && ( abs( gen[k].pdgId() ) == 11 || abs( gen[k].pdgId() ) == 13 || abs( gen[k].pdgId() )==15 ) ) prodId*=gen[k].pdgId(); 
                }
                if(mctruthmode==15   && abs(prodId)!=225)continue; //skip not tautau
                if(mctruthmode==1113 && abs(prodId)==225)continue; //skip tautau
            }

            if(isMC_VBF || isMC_GG || mctruthmode==125){
               LorentzVector higgs(0,0,0,0);
  	       for(unsigned int igen=0; igen<gen.size(); igen++){
                  if(!gen[igen].isHardProcess()) continue;
	          if(abs(gen[igen].pdgId())>=11 && abs(gen[igen].pdgId())<=16){ higgs += gen[igen].p4(); }
	       }
	       if(mctruthmode==125) {
	         if( isMC_125OnShell && higgs.mass()> 180) continue;
	         if(!isMC_125OnShell && higgs.mass()<=180) continue;
	       }
	  
               //Line shape weights 
               if(isMC_VBF || isMC_GG){    
                  mon.fillHisto("higgsMass_raw",    "all", higgs.mass(), weight);

                  //WEIGHT for LineShape (NNLO kFactors + Interf), split in shape and scale unc with the following format:  scaleNominal, shapeNominal, scaleDown, shapeDown, scaleUp, shapeUp
                  //compute weight correction for all NR shapes
                  for(unsigned int nri=0;nri<NRparams.size();nri++){ 
                     std::vector<std::pair<double, TGraph*> > shapeWgtGr = hLineShapeGrVec[NRparams[nri] ];
                     for(size_t iwgt=0; iwgt<shapeWgtGr.size(); iwgt++){ 
                        lShapeWeights[nri][iwgt*2+0]=shapeWgtGr[iwgt].first;
                        lShapeWeights[nri][iwgt*2+1]=shapeWgtGr[iwgt].second?shapeWgtGr[iwgt].second->Eval(higgs.mass()):1.0;
                     }
                     mon.fillHisto(TString("higgsMass_shape"      )+NRsuffix[nri], "all", higgs.mass(), weight*lShapeWeights[nri][1] );
                     mon.fillHisto(TString("higgsMass_shape&scale")+NRsuffix[nri], "all", higgs.mass(), weight*lShapeWeights[nri][1]*lShapeWeights[nri][0] );

                     //printf("NRI BW=%i --> %6.2e %6.2e %6.2e %6.2e %6.2e %6.2e\n", nri, lShapeWeights[nri][0], lShapeWeights[nri][1], lShapeWeights[nri][2], lShapeWeights[nri][3], lShapeWeights[nri][4], lShapeWeights[nri][5]);

                     //scale Up/Down by Nominal
                     if(lShapeWeights[nri][0]>0){lShapeWeights[nri][2]/=lShapeWeights[nri][0];  lShapeWeights[nri][4]/=lShapeWeights[nri][0];}
                     if(lShapeWeights[nri][1]>0){lShapeWeights[nri][3]/=lShapeWeights[nri][1];  lShapeWeights[nri][5]/=lShapeWeights[nri][1];}

                     //printf("NRI AW=%i --> %6.2e %6.2e %6.2e %6.2e %6.2e %6.2e\n", nri, lShapeWeights[nri][0], lShapeWeights[nri][1], lShapeWeights[nri][2], lShapeWeights[nri][3], lShapeWeights[nri][4], lShapeWeights[nri][5]);
                  }  
                  weight *= lShapeWeights[0][0]*lShapeWeights[0][1]; 
               }
            }
          }

          //apply trigger and require compatibilitiy of the event with the PD
          edm::TriggerResultsByName tr = ev.triggerResultsByName("HLT");
          if(!tr.isValid())return false;

          float triggerPrescale(1.0),triggerThreshold(0), triggerThresholdHigh(99999);
          char photonTriggerTreshName[255];
          bool mumuTrigger        = utils::passTriggerPatterns(tr, "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");                  
          bool muTrigger          = utils::passTriggerPatterns(tr, "HLT_IsoMu20_v*", "HLT_IsoTkMu20_v*", "HLT_IsoMu27_v*");                                               
          bool eeTrigger          = utils::passTriggerPatterns(tr, "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");       
          bool eTrigger           = utils::passTriggerPatterns(tr, "HLT_Ele23_WPLoose_Gsf_v*", "HLT_Ele22_eta2p1_WP75_Gsf_v*", "HLT_Ele22_eta2p1_WPLoose_Gsf_v*"); 
          bool emuTrigger         = utils::passTriggerPatterns(tr, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*", "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*");  
          bool photonTrigger      = patUtils::passPhotonTrigger(ev, triggerThreshold, triggerPrescale, triggerThresholdHigh);                                                                        
          bool passTrigger        = mumuTrigger||muTrigger||eeTrigger||eTrigger||emuTrigger||photonTrigger;

          if(  mumuTrigger)mon.fillHisto("trigger", "raw", 0 , weight);
          if(    muTrigger)mon.fillHisto("trigger", "raw", 1 , weight);
          if(    eeTrigger)mon.fillHisto("trigger", "raw", 2 , weight);
          if(     eTrigger)mon.fillHisto("trigger", "raw", 3 , weight);
          if(   emuTrigger)mon.fillHisto("trigger", "raw", 4 , weight);
          if(photonTrigger)mon.fillHisto("trigger", "raw", 5 , weight);

          if(!isMC && passTrigger){ //avoid double counting of events from different PD
             if(filterOnlyMUMU)     { passTrigger = mumuTrigger;}
             if(filterOnlyMU)       { passTrigger = muTrigger     && !mumuTrigger;}
             if(filterOnlyEE)       { passTrigger = eeTrigger     && !muTrigger  && !mumuTrigger;}
             if(filterOnlyE)        { passTrigger = eTrigger      && !eeTrigger  && !muTrigger && !mumuTrigger; }
             if(filterOnlyEMU)      { passTrigger = emuTrigger    && !eTrigger   && !eeTrigger && !muTrigger && !mumuTrigger; }
             if(filterOnlyPhoton)   { passTrigger = photonTrigger && !emuTrigger && !eTrigger  && !eeTrigger && !muTrigger && !mumuTrigger;}
          }

          if(passTrigger){
             if(  mumuTrigger)mon.fillHisto("trigger", "cleaned", 0 , weight);
             if(    muTrigger)mon.fillHisto("trigger", "cleaned", 1 , weight);
             if(    eeTrigger)mon.fillHisto("trigger", "cleaned", 2 , weight);
             if(     eTrigger)mon.fillHisto("trigger", "cleaned", 3 , weight);
             if(   emuTrigger)mon.fillHisto("trigger", "cleaned", 4 , weight);
             if(photonTrigger)mon.fillHisto("trigger", "cleaned", 5 , weight);
          }

          if(photonTrigger){sprintf(photonTriggerTreshName, "PhoTrg%i", int(triggerThreshold));}

          //ONLY RUN ON THE EVENTS THAT PASS OUR TRIGGERS
           if(!passTrigger && !photonTriggerStudy)continue;        


         //##############################################   EVENT PASSED THE TRIGGER   ######################################
          int metFilterValue = metFiler.passMetFilterInt( ev );
          mon.fillHisto("metFilter_eventflow", "", metFilterValue, weight);
          if( metFilterValue!=0 ) continue;	 //Note this must also be applied on MC
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
          pat::MET met = mets[0]; 

          pat::METCollection puppimets;
          fwlite::Handle< pat::METCollection > puppimetsHandle;
          puppimetsHandle.getByLabel(ev, "slimmedMETsPuppi");
          if(puppimetsHandle.isValid()){ puppimets = *puppimetsHandle;}
          LorentzVector puppimet = puppimets[0].p4(); 

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


         //MC crap for photon studies
         if(photonTrigger && (isWGmc || isZGmc)){
           int nge(0), ngm(0), ngt(0), ngj(0), ngnu(0);
           bool zFound(false), wFound(false);
           for(size_t ig=0; ig<gen.size(); ig++){
             if(gen[ig].status()!=3) continue;
             int id(abs(gen[ig].pdgId()));
             if(id==23) zFound=true;
             if(id==24) wFound=true;
             if(id==11) nge++;
             if(id==13) ngm++;
             if(id==15) ngt++;
             if(id==12 || id==14 || id==16) ngnu++;
             if((wFound || zFound) && id<6) ngj++;
           }
           if(zFound){
             int decBin=0;
             if(nge==2) decBin=1;
             if(ngm==2) decBin=2;
             if(ngt==2) decBin=3;
             if(ngj>=2) decBin=4;
             if(ngnu==2) decBin=5;
             mon.fillHisto("zdecays","",decBin,1);
           }
           if(wFound){
             int decBin=0;
             if(nge==1 && ngnu==1) decBin=1;
             if(ngm==1 && ngnu==1) decBin=2;
             if(ngt==1 && ngnu==1) decBin=3;
             if(ngj>=2) decBin=4;
             mon.fillHisto("wdecays","",decBin,1);
           }
         }

         //Resolve G+jet/QCD mixing (avoid double counting of photons)
         if (isMC_GJet || isMC_QCD ) {
           // iF GJet sample; accept only event with prompt photons                                                                 
           // if QCD sample; reject events with prompt photons in final state                                                       
             bool gPromptFound=false;
             for(size_t ig=0; ig<gen.size(); ig++){
               if((abs(gen[ig].pdgId())==22) && gen[ig].isPromptFinalState())  gPromptFound=true;               
             }
             if ( (isMC_GJet) && (!gPromptFound) ) continue; //reject event
             if ( (isMC_QCD) && gPromptFound ) continue; //reject event
         }

     		 //Electroweak corrections to ZZ and WZ(soon) simulations
     		 double ewkCorrectionsWeight = 1.;
     		 double ewkCorrections_error = 0.;
     		 if(isMC_ZZ2l2nu) ewkCorrectionsWeight = EwkCorrections::getEwkCorrections(dtag, gen, ewkTable, eventInfo, ewkCorrections_error);
     		 double ewkCorrections_up = (ewkCorrectionsWeight + ewkCorrections_error)/ewkCorrectionsWeight;
     	         double ewkCorrections_down = (ewkCorrectionsWeight - ewkCorrections_error)/ewkCorrectionsWeight;
     
       	 //final event weight
       	 weight *= ewkCorrectionsWeight;
    
    		 //NNLO corrections on ZZ2l2nu
    		 double ZZ_NNLOcorrectionsWeight =1.;
				 double mzz = - 404; // will be filled by getNNLOCorrections 
    		 if(isMC_ZZ2l2nu) ZZ_NNLOcorrectionsWeight = ZZatNNLO::getNNLOCorrections(dtag, gen, ZZ_NNLOTable, mzz);
				 
				 //final event weight
				 if(isMC_ZZ2l2nu) mon.fillHisto("mzz", "qqZZ_atNLO", mzz, weight);
				 weight *= ZZ_NNLOcorrectionsWeight;
				 if(isMC_ZZ2l2nu) mon.fillHisto("mzz", "qqZZ_atNNLO", mzz, weight);

         //
         //
         // BELOW FOLLOWS THE ANALYSIS OF THE MAIN SELECTION WITH N-1 PLOTS
         //
         //

         //
         // PHOTON ANALYSIS
         //
         pat::PhotonCollection selPhotons;	    
         int nPho55=0; int nPho100=0;
         for(size_t ipho=0; ipho<photons.size(); ipho++){
	    pat::Photon photon = photons[ipho]; 
 	    mon.fillHisto("phopt", "trg", photon.pt(), weight);
	    mon.fillHisto("phoeta", "trg", photon.eta(), weight);           
            //printf("photon pt=%6.2f eta=%+6.2f\n", photon.pt(), photon.eta());

//            if(photonTrigger && (photon.pt()<triggerThreshold || photon.pt()>triggerThresholdHigh))continue;
             if(photonTrigger && (photon.pt()<(triggerThreshold) || photon.pt()>(triggerThresholdHigh+10)))continue;           //add 5GeV to avoid trigger inefficiencies
            //printf("A\n");

            //calibrate photon energy
            //PhotonEnCorrector.calibrate(photon, ev.eventAuxiliary().run(), edm::StreamID::invalidStreamID()); 

            if(photon.pt()<55)continue;
            //printf("B\n");
            if(fabs(photon.superCluster()->eta())>1.4442 ) continue;
            //printf("C\n");

	    if(!patUtils::passId(photon, rho, patUtils::llvvPhotonId::Tight)) continue;
            //printf("D\n");

            selPhotons.push_back(photon);
            if(photon.pt()>55)nPho55++;
            if(photon.pt()>100)nPho100++;
         }           

         //
         // LEPTON ANALYSIS
         //
         
         //start by merging electrons and muons
         std::vector<patUtils::GenericLepton> leptons;
         for(size_t l=0;l<electrons.size();l++){leptons.push_back(patUtils::GenericLepton(electrons[l]));}      
         for(size_t l=0;l<muons    .size();l++){leptons.push_back(patUtils::GenericLepton(muons    [l]));}      
         std::sort(leptons.begin(),   leptons.end(), utils::sort_CandidatesByPt);

         std::vector<patUtils::GenericLepton> selLeptons, extraLeptons;
         LorentzVector muDiff(0,0,0,0);
         LorentzVector elDiff(0,0,0,0);
         for(size_t ilep=0; ilep<leptons.size(); ilep++){
             bool passKin(true),passId(true),passIso(true);
             bool passLooseLepton(true), passSoftMuon(true), passSoftElectron(true), passVetoElectron(true);
             int lid=leptons[ilep].pdgId();

             //no need for charge info any longer
             lid=abs(lid);
             TString lepStr( lid==13 ? "mu" : "e");



             //veto nearby photon (loose electrons are many times photons...)
             double minDRlg(9999.);
             for(size_t ipho=0; ipho<selPhotons.size(); ipho++){
               minDRlg=TMath::Min(minDRlg,deltaR(leptons[ilep].p4(),selPhotons[ipho].p4()));
             }
             //printf("lepton pt=%6.2f eta=%+6.2f %3i drtoPhoton%6.2f\n", leptons[ilep].pt(), leptons[ilep].eta(), lid, minDRlg);

             if(minDRlg<0.1) continue;

             //Cut based identification
             passId = lid==11?patUtils::passId(leptons[ilep].el, vtx[0], patUtils::llvvElecId::Tight) : patUtils::passId(leptons[ilep].mu, vtx[0], patUtils::llvvMuonId::Tight);
             passLooseLepton &= lid==11?patUtils::passId(leptons[ilep].el, vtx[0], patUtils::llvvElecId::Loose) : patUtils::passId(leptons[ilep].mu, vtx[0], patUtils::llvvMuonId::Loose);
             passSoftMuon &= lid==11? false : patUtils::passId(leptons[ilep].mu, vtx[0], patUtils::llvvMuonId::Soft);

             //isolation
             passIso = lid==11?patUtils::passIso(leptons[ilep].el,  patUtils::llvvElecIso::Tight) : patUtils::passIso(leptons[ilep].mu,  patUtils::llvvMuonIso::Tight);
             passLooseLepton &= lid==11?patUtils::passIso(leptons[ilep].el,  patUtils::llvvElecIso::Loose) : patUtils::passIso(leptons[ilep].mu,  patUtils::llvvMuonIso::Loose);

             //apply muon corrections
             if(abs(lid)==13 && passIso && passId){
                 passSoftMuon=false;
                 if(muCor){
                   float qter;
                   TLorentzVector p4(leptons[ilep].px(),leptons[ilep].py(),leptons[ilep].pz(),leptons[ilep].energy());
                   if(isMC){muCor->momcor_mc  (p4, lid<0 ? -1 :1, 0, qter);
                   }else if (is2015data){   muCor->momcor_data(p4, lid<0 ? -1 :1, 0, qter); 
                   }


                   if(ev.eventAuxiliary().event()==869607902 || ev.eventAuxiliary().event()==471854508){
                      printf("\nevent = %lli %i lepton pt, eta, phi %f %f %f --> %f %f %f\n", ev.eventAuxiliary().event(), int(ilep), leptons[ilep].pt(),leptons[ilep].eta(),leptons[ilep].phi(), p4.Pt(),p4.Eta(),p4.Phi() );
                   }

                   muDiff -= leptons[ilep].p4();
                   leptons[ilep].setP4(LorentzVector(p4.Px(),p4.Py(),p4.Pz(),p4.E() ) );
                   muDiff += leptons[ilep].p4();
                 }
               }

             //apply electron corrections             
             if(abs(lid)==11  && passIso && passId){
                elDiff -= leptons[ilep].p4();                   
                if (isMC || is2015data){
                	ElectronEnCorrector.calibrate(leptons[ilep].el, ev.eventAuxiliary().run(), edm::StreamID::invalidStreamID()); 
                	leptons[ilep] = patUtils::GenericLepton(leptons[ilep].el); //recreate the generic lepton to be sure that the p4 is ok
                }
                elDiff += leptons[ilep].p4();                 
             }

              //kinematics
             float leta = fabs(lid==11 ?  leptons[ilep].el.superCluster()->eta() : leptons[ilep].eta());
             if(leta> (lid==11 ? 2.5 : 2.4) )            passKin=false;
             if(lid==11 && (leta>1.4442 && leta<1.5660)) passKin=false;
             passLooseLepton &= passKin;
             passSoftMuon    &= passKin;
             if(lid==13){
               if(leptons[ilep].pt()<10) passLooseLepton=false;
               if(leptons[ilep].pt()<3)  passSoftMuon=false;
             }else if(lid==11){
               if(leptons[ilep].pt()<10) passLooseLepton=false;
             }
             if(leptons[ilep].pt()<25) passKin=false;
            


             if(passId && passIso && passKin)          selLeptons.push_back(leptons[ilep]); 
             else if(passLooseLepton || passSoftMuon)  extraLeptons.push_back(leptons[ilep]);

             if(ev.eventAuxiliary().event()==869607902 || ev.eventAuxiliary().event()==471854508){
                if(lid==13)printf("muon%i Id=%i Iso=%i LooseId=%i kin=%i\n", int(ilep), passId?1:0, passIso?1:0, passLooseLepton?1:0, passKin?1:0);
             }


           }
           std::sort(selLeptons.begin(),   selLeptons.end(), utils::sort_CandidatesByPt);
           std::sort(extraLeptons.begin(), extraLeptons.end(), utils::sort_CandidatesByPt);

           //update the met for lepton energy scales
           met.setP4(met.p4() - muDiff - elDiff); //note this also propagates to all MET uncertainties
           met.setUncShift(met.px() - muDiff.px()*0.01, met.py() - muDiff.py()*0.01, met.sumEt() - muDiff.pt()*0.01, pat::MET::METUncertainty::MuonEnUp);   //assume 1% uncertainty on muon rochester
           met.setUncShift(met.px() + muDiff.px()*0.01, met.py() + muDiff.py()*0.01, met.sumEt() + muDiff.pt()*0.01, pat::MET::METUncertainty::MuonEnDown); //assume 1% uncertainty on muon rochester
           met.setUncShift(met.px() - elDiff.px()*0.01, met.py() - elDiff.py()*0.01, met.sumEt() - elDiff.pt()*0.01, pat::MET::METUncertainty::ElectronEnUp);   //assume 1% uncertainty on electron scale correction
           met.setUncShift(met.px() + elDiff.px()*0.01, met.py() + elDiff.py()*0.01, met.sumEt() + elDiff.pt()*0.01, pat::MET::METUncertainty::ElectronEnDown); //assume 1% uncertainty on electron scale correction
         //
         //JET/MET ANALYSIS
         //
         //add scale/resolution uncertainties and propagate to the MET      
         if(isMC || is2015data) utils::cmssw::updateJEC(jets,jesCor,totalJESUnc,rho,vtx.size(),isMC); 

         //select the jets
         std::map<string, pat::JetCollection> selJetsVar;
         std::map<string, int   > njetsVar;
         std::map<string, int   > nbtagsVar;
         std::map<string, double> mindphijmetVar;
         for(unsigned int ivar=0;ivar<jetVarNames.size();ivar++){njetsVar[jetVarNames[ivar]] = 0;}  //initialize
         for(unsigned int ivar=0;ivar<jetVarNames.size();ivar++){mindphijmetVar[jetVarNames[ivar]] = 9999.0;}  //initialize
         nbtagsVar[""] = 0; nbtagsVar["_eff_bup"] = 0; nbtagsVar["_eff_bdown"] = 0;  //initialize

         for(size_t ijet=0; ijet<jets.size(); ijet++){
             pat::Jet jet = jets[ijet]; //copy the jet, such that we can update it

             if(jet.pt()<15 || fabs(jet.eta())>4.7 ) continue;

             //mc truth for this jet
             const reco::GenJet* genJet=jet.genJet();
             TString jetType( genJet && genJet->pt()>0 ? "truejetsid" : "pujetsid" );
             
             //cross-clean with selected leptons and photons
             double minDRlj(9999.); for(size_t ilep=0; ilep<selLeptons.size(); ilep++)  minDRlj = TMath::Min( minDRlj, deltaR(jet,selLeptons[ilep]) );
             double minDRlg(9999.); for(size_t ipho=0; ipho<selPhotons.size(); ipho++)  minDRlg = TMath::Min( minDRlg, deltaR(jet,selPhotons[ipho]) );
             if(minDRlj<0.4 || minDRlg<0.4) continue;
             
             //jet id
             bool passPFloose = patUtils::passPFJetID("Loose", jet);
             bool passLooseSimplePuId = true; //patUtils::passPUJetID(jet); //FIXME Broken in miniAOD V2 : waiting for JetMET fix. (Hugo)
             if(jet.pt()>30){
                 mon.fillHisto(jetType,"",fabs(jet.eta()),0);
                 if(passPFloose)                        mon.fillHisto("jetId", jetType,fabs(jet.eta()),1);
                 if(passLooseSimplePuId)                mon.fillHisto("jetId", jetType,fabs(jet.eta()),2);
                 if(passPFloose && passLooseSimplePuId) mon.fillHisto("jetId", jetType,fabs(jet.eta()),3);
             }
             if(!passPFloose || !passLooseSimplePuId) continue; 


            //check for btagging
            if(jet.pt()>30 && fabs(jet.eta())<2.5){
              bool hasCSVtag = (jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")>btagLoose);
              bool hasCSVtagUp = hasCSVtag;  
              bool hasCSVtagDown = hasCSVtag;
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
            }


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
         }
         //sort all jet collection by pT
         for(auto jetCollIt = selJetsVar.begin(); jetCollIt!=selJetsVar.end(); jetCollIt++){
            std::sort(jetCollIt->second.begin(), jetCollIt->second.end(), utils::sort_CandidatesByPt);
         }


         //save weight
         double initialWeight = weight;

         //compute scale uncertainty once and for all
         std::pair<double, double> scaleUncVar = patUtils::scaleVariation(ev);  //compute it only once          

         // LOOP ON SYSTEMATIC VARIATION FOR THE STATISTICAL ANALYSIS
         for(size_t ivar=0; ivar<nvarsToInclude; ivar++){
          if(!isMC && ivar>0 ) continue; //loop on variation only for MC samples

          //start from a nominal
          float weight = initialWeight;

           //Theoretical Uncertanties: PDF, Alpha and Scale
           if(varNames[ivar]=="_th_factup")     weight *= std::max(0.9, std::min(scaleUncVar.first , 1.1)); 
           if(varNames[ivar]=="_th_factdown")   weight *= std::max(0.9, std::min(scaleUncVar.second, 1.1));
           if(varNames[ivar]=="_th_alphas")     weight *= patUtils::alphaVariation(ev);
           if(varNames[ivar]=="_th_pdf")        weight *= patUtils::pdfVariation(ev);

           //EwkCorrections variation
           if ( varNames[ivar]=="_th_ewkup")    weight *= ewkCorrections_up;
           if ( varNames[ivar]=="_th_ewkdown")  weight *= ewkCorrections_down;

           //pileup variations
           if(varNames[ivar]=="_puup")          weight *= puWeightUp;
           if(varNames[ivar]=="_pudown")        weight *= puWeightDown;
          
           //recompute MET with variation
           LorentzVector imet = met.corP4(metcor);
           if(varNames[ivar]=="_scale_jup")      imet = met.shiftedP4(pat::MET::METUncertainty::JetEnUp           , metcor);
           if(varNames[ivar]=="_scale_jdown")    imet = met.shiftedP4(pat::MET::METUncertainty::JetEnDown         , metcor);
           if(varNames[ivar]=="_res_jup")        imet = met.shiftedP4(pat::MET::METUncertainty::JetResUp          , metcor);
           if(varNames[ivar]=="_res_jdown")      imet = met.shiftedP4(pat::MET::METUncertainty::JetResDown        , metcor);
           if(varNames[ivar]=="_scale_umetup")   imet = met.shiftedP4(pat::MET::METUncertainty::UnclusteredEnUp   , metcor);              
           if(varNames[ivar]=="_scale_umetdown") imet = met.shiftedP4(pat::MET::METUncertainty::UnclusteredEnDown , metcor);              
           if(varNames[ivar]=="_scale_mup")      imet = met.shiftedP4(pat::MET::METUncertainty::MuonEnUp          , metcor);
           if(varNames[ivar]=="_scale_mdown")    imet = met.shiftedP4(pat::MET::METUncertainty::MuonEnDown        , metcor);             
           if(varNames[ivar]=="_scale_eup")      imet = met.shiftedP4(pat::MET::METUncertainty::ElectronEnUp      , metcor); 
           if(varNames[ivar]=="_scale_edown")    imet = met.shiftedP4(pat::MET::METUncertainty::ElectronEnDown    , metcor);

           auto& selJets      = selJetsVar[""];        if(selJetsVar    .find(varNames[ivar].Data())!=selJetsVar    .end())selJets     = selJetsVar    [varNames[ivar].Data()];            
           auto& njets        = njetsVar [""];         if(njetsVar      .find(varNames[ivar].Data())!=njetsVar      .end())njets       = njetsVar      [varNames[ivar].Data()];
           auto& nbtags       = nbtagsVar[""];         if(nbtagsVar     .find(varNames[ivar].Data())!=nbtagsVar     .end())nbtags      = nbtagsVar     [varNames[ivar].Data()];
           auto& mindphijmet  = mindphijmetVar[""];    if(mindphijmetVar.find(varNames[ivar].Data())!=mindphijmetVar.end())mindphijmet = mindphijmetVar[varNames[ivar].Data()];

            //
            // ASSIGN CHANNEL
            //
            double weightBefLoop = weight;
            for(unsigned int L=0;L<3;L++){  //Loop to assign a Z-->ll channel to photons
               if(L>0 && !(photonTrigger && gammaWgtHandler) )continue; //run it only for photon reweighting
               weight = weightBefLoop;
               std::vector<TString> chTags;
               TString evCat;       
               int dilId(1);
               LorentzVector boson(0,0,0,0);
               if(selLeptons.size()==2  && !gammaWgtHandler){  //this is not run if photon reweighting is activated to avoid mixing           
                   for(size_t ilep=0; ilep<2; ilep++){
                       dilId *= selLeptons[ilep].pdgId();
                       int id(abs(selLeptons[ilep].pdgId()));
                       weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[ilep].pt(), selLeptons[ilep].eta(), id,  id ==11 ? "tight"    : "tight"    ).first : 1.0; //ID 
                       weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[ilep].pt(), selLeptons[ilep].eta(), id,  id ==11 ? "tightiso" : "tightiso" ).first : 1.0; //ISO w.r.t ID
                       boson += selLeptons[ilep].p4();
                   }        
                   //check the channel
                   if( abs(dilId)==121){  chTags.push_back("ee");   chTags.push_back("ll"); }
                   if( abs(dilId)==169){  chTags.push_back("mumu"); chTags.push_back("ll"); }
                   if( abs(dilId)==143){  chTags.push_back("emu");  }           
                   //if(isMC && abs(dilId)==169)weight *= lepEff.getTriggerEfficiencySF(selLeptons[0].pt(), selLeptons[0].eta(), selLeptons[1].pt(), selLeptons[1].eta(), dilId).first;  //commented for ee as inefficiencies should be covered by the singleMu/El triggers
                   evCat=eventCategoryInst.GetCategory(selJets,boson);            
               }else if(selPhotons.size()==1 && photonTrigger){
                   dilId=22;
                   if(L==0)                         {chTags.push_back("gamma");
                   }else if(L==1 && gammaWgtHandler){chTags.push_back("ee");   chTags.push_back("ll");
                   }else if(L==2 && gammaWgtHandler){chTags.push_back("mumu"); chTags.push_back("ll");
                   }else{ continue;
                   }
                   boson = selPhotons[0].p4();
                   evCat=eventCategoryInst.GetCategory(selJets,boson);            
                   if(L>0 && gammaWgtHandler)boson = gammaWgtHandler->getMassiveP4(boson, string(L==1?"ee":"mumu")+evCat);
                   std::vector<Float_t> photonVars;
                   photonVars.push_back(boson.pt());           
                   float photonWeightMain=1.0;
                   if(L>0 && gammaWgtHandler)photonWeightMain=gammaWgtHandler->getWeightFor(photonVars,string(L==1?"ee":"mumu")+evCat);
                   //if(L>0 && gammaWgtHandler)printf("Photon pT = %6.2f --> prescale=%6.2f weight=%6.2E forL=%i  cat=%s\n", boson.pt(), triggerPrescale, photonWeightMain, L, (string(L==1?"ee":"mumu")+evCat).Data());
                   weight *= triggerPrescale * photonWeightMain;
               }else{
                  continue;
               }

               std::vector<TString> tags(1,"all");
               for(size_t ich=0; ich<chTags.size(); ich++){
                 tags.push_back( chTags[ich] );
                 tags.push_back( chTags[ich]+evCat );
               }

               //////////////////////////
               //                      //
               //  BASELINE SELECTION  //
               //                      //
               //////////////////////////

               bool passMass(fabs(boson.mass()-91)<15);
               bool passQt(boson.pt()>55);
               bool passThirdLeptonVeto( selLeptons.size()==2 && extraLeptons.size()==0 );
               bool passBtags(nbtags==0); 
               bool passMinDphijmet( njets==0 || mindphijmet>0.5);

               if(dilId==22){
                   passMass=photonTrigger;
                   passThirdLeptonVeto=(selLeptons.size()==0 && extraLeptons.size()==0);
               }

              if(varNames[ivar]=="_lepveto" && !passThirdLeptonVeto){
                 int NExtraLep = std::max(0, int(selLeptons.size()) + int(extraLeptons.size()) - 2);
                 if(((rand()%1000)/1000.0) < pow(0.04, NExtraLep))passThirdLeptonVeto=true;  //4% Id uncertainty exponent Number of aditional leptons
              }

              double mt=0;
              if(passQt && passThirdLeptonVeto && passMinDphijmet && (boson.mass()>40 && boson.mass()<200)) mt=higgs::utils::transverseMass(boson,imet,true);             


               if(ivar==0){  //fill control plots only for the nominal systematic

                  mon.fillHisto("nleptons",tags,selLeptons.size(), weight);
                  mon.fillHisto("npho", tags, selPhotons.size(), weight);
                  mon.fillHisto("npho55", tags, nPho55, weight);
                  mon.fillHisto("npho100", tags, nPho100, weight);
                  if(photonTrigger && selPhotons.size()>0)mon.fillHisto("photonpt", tags[tags.size()-1]+photonTriggerTreshName,   selPhotons[0].pt(), weight);
                 

                  // Photon trigger efficiencies
                  // Must be run without the photonTrigger requirement on top of of the Analysis.
                  if (photonTriggerStudy && selPhotons.size() ){
                    TString tag="trigger";
                    pat::Photon iphoton = selPhotons[0];
               
                    mon.fillHisto("phopt", tag, iphoton.pt(),weight);
                    mon.fillHisto("phoeta", tag, iphoton.eta(), weight);
                    trigUtils::photonControlSample(ev, iphoton, mon, tag);
                    trigUtils::photonControlEff(ev, iphoton, mon, tag);

                    for(size_t itag=0; itag<tags.size(); itag++){		
                      //update the weight
                      TString icat=tags[itag];
                      mon.fillHisto("phopt", icat, iphoton.pt(),weight);
                      mon.fillHisto("phoeta", icat, iphoton.eta(),weight);
                      trigUtils::photonControlSample(ev, iphoton, mon, icat);
                      trigUtils::photonControlEff(ev, iphoton, mon, icat);
                    }
                  } // end Trigger efficiencies


    

                  mon.fillHisto("eventflow",  tags,0,weight);
                  mon.fillHisto("nvtxraw",  tags,vtx.size(),weight/puWeight);
                  mon.fillHisto("nvtx",  tags,vtx.size(),weight);
                  mon.fillHisto("rho",  tags,rho,weight);

                  if(chTags.size()==0) continue;
                  mon.fillHisto("eventflow",  tags,1,weight);
                  if(dilId!=22){
                    mon.fillHisto("leadpt",      tags,selLeptons[0].pt(),weight); 
                    mon.fillHisto("trailerpt",   tags,selLeptons[1].pt(),weight); 
                    mon.fillHisto("leadeta",     tags,fabs(selLeptons[0].eta()),weight); 
                    mon.fillHisto("trailereta",  tags,fabs(selLeptons[1].eta()),weight); 
                  }
           
                  mon.fillHisto("zmass", tags,boson.mass(),weight); 
                  mon.fillHisto("zy",    tags,fabs(boson.Rapidity()),weight); 

                  if(passMass){
                    mon.fillHisto("eventflow",tags, 2,weight);
                    mon.fillHisto("zpt",      tags, boson.pt(),weight);
                    mon.fillHisto("zpt_rebin",tags, boson.pt(),weight,true);
                    if(imet.pt()>125)mon.fillHisto("zptMet125",      tags, boson.pt(),weight);


                    //these two are used to reweight photon -> Z, the 3rd is a control
                    mon.fillHisto("qt",       tags, boson.pt(),weight,true); 
                    mon.fillHisto("qtraw",    tags, boson.pt(),weight/triggerPrescale,true); 

                    if(passQt){
                      mon.fillHisto("eventflow",tags,3,weight);
                      int nExtraLeptons((selLeptons.size()-2)+extraLeptons.size());
                      mon.fillHisto("nextraleptons",tags,nExtraLeptons,weight);
                      if(nExtraLeptons>0){
                        LorentzVector thirdLepton(selLeptons.size()>2 ?  selLeptons[1].p4() : extraLeptons[0].p4());
                        double dphi=fabs(deltaPhi(thirdLepton.phi(),imet.phi()));
                        double mt3rd=TMath::Sqrt(2*thirdLepton.pt()*imet.pt()*(1-TMath::Cos(dphi)));
                        mon.fillHisto("thirdleptonpt",tags,thirdLepton.pt(),weight);
                        mon.fillHisto("thirdleptoneta",tags,fabs(thirdLepton.eta()),weight);
                        mon.fillHisto("thirdleptonmt",tags,mt3rd,weight);
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
                        
                        if(passBtags){
                          mon.fillHisto("eventflow",tags,5,weight);

                          mon.fillHisto( "mindphijmet",tags,mindphijmet,weight);
                          if(imet.pt()>25)mon.fillHisto( "mindphijmet25",tags,mindphijmet,weight);
                          if(imet.pt()>50)mon.fillHisto( "mindphijmet50",tags,mindphijmet,weight);
                          if(imet.pt()>80)mon.fillHisto( "mindphijmetNM1",tags,mindphijmet,weight);
                          if(passMinDphijmet){
                            mon.fillHisto("eventflow",tags,6,weight);
                           
                            //this one is used to sample the boson mass: cuts may shape Z lineshape
                            mon.fillHisto("qmass",       tags, boson.mass(),weight); 
                            mon.fillHisto( "njets",tags,njets,weight);
       
                            double b_dphi=fabs(deltaPhi(boson.phi(),imet.phi()));
                            mon.fillHisto( "metphi",tags,imet.phi(),weight);
                            mon.fillHisto( "metphiUnCor",tags,met.corP4(pat::MET::METCorrectionLevel::Type1).phi(),weight);
                            mon.fillHisto( "bosonnvtx",tags,vtx.size(),weight);                                                               
                            mon.fillHisto( "bosoneta",tags,boson.eta(),weight);                                                                                
                            mon.fillHisto( "bosonphi",tags,boson.phi(),weight);                                                               
                            mon.fillHisto( "bosonphiHG",tags,boson.phi(),weight);
                            mon.fillHisto( "dphi_boson_met",tags,b_dphi,weight);
       
                            mon.fillHisto( "met",tags,imet.pt(),weight,true);
                            mon.fillHisto( "metpuppi",tags,puppimet.pt(),weight,true);
                            mon.fillHisto( "balance",tags,imet.pt()/boson.pt(),weight);

                            TVector2 met2(imet.px(),imet.py());
                            TVector2 boson2(boson.px(), boson.py());
                            double axialMet(boson2*met2); axialMet/=-boson.pt();
                            mon.fillHisto( "axialmet",tags,axialMet,weight);
       
                            mon.fillHisto( "mt",tags,mt,weight,true);

                            if(imet.pt()>optim_Cuts1_met[0]) {
                               mon.fillHisto( "mtcheckpoint",  tags, mt,       weight, true);
                               mon.fillHisto( "metcheckpoint", tags, imet.pt(), weight, true);
                            }

                            if(imet.pt()>80){
                              mon.fillHisto("eventflow",tags,7,weight);
                              mon.fillHisto( "mtNM1",tags,mt,weight,true);
                              mon.fillHisto( "balanceNM1",tags,imet.pt()/boson.pt(),weight);
                              mon.fillHisto( "axialmetNM1",tags,axialMet,weight);
                            }

                            if(imet.pt()>125){
                              mon.fillHisto("eventflow",tags,8,weight);

                              mon.fillHisto( "metfinal",tags,imet.pt(),weight,true);                      
                              mon.fillHisto( "mtfinal",tags,mt,weight,true);
                              mon.fillHisto( "mindphijmetfinal",tags,mindphijmet,weight);
                              mon.fillHisto( "njetsfinal",tags,njets,weight);
                              if(!isMC){
                                 char buffer[1024];
                                 sprintf(buffer, "\ncat=%s %9i:%6i:%9lli @ %50s\n",  tags[tags.size()-1].Data(), ev.eventAuxiliary().run(), ev.eventAuxiliary().luminosityBlock(), ev.eventAuxiliary().event(), urls[f].c_str() );  debugText+=buffer; 
                                 sprintf(buffer, " - nLep=%2i nSoftLept=%2i nPhotons=%2i  nJets=%2i\n", int(selLeptons.size()), int(extraLeptons.size()), int(selPhotons.size()), int(selJets.size())  ); debugText+=buffer;                               
                                 sprintf(buffer, " - MET=%8.2f mT=%8.2f nvtx=%3i\n", imet.pt(), mt, int(vtx.size()) ); debugText+=buffer;
                                 sprintf(buffer, " - MET type1XY=%8.2f type1=%8.2f uncorrected=%8.2f\n", met.corP4(pat::MET::METCorrectionLevel::Type1XY).pt(), met.corP4(pat::MET::METCorrectionLevel::Type1).pt(), met.corP4(pat::MET::METCorrectionLevel::Raw).pt() ); debugText+=buffer;
                                 sprintf(buffer, " - LeptonScale changes on MET mu=%8.2f  el=%8.2f\n", muDiff.pt(), elDiff.pt() ); debugText+=buffer;                        
                                 sprintf(buffer, " - Z pT=%6.2f eta=%+6.2f phi=%+6.2f\n", boson.pt(), boson.eta(), boson.phi() ); debugText+=buffer;
                                 if(selLeptons.size()>0)sprintf(buffer, " - lep0 Id=%+3i pT=%6.2f, eta=%+6.2f phi=%+6.2f\n", selLeptons[0].pdgId(), selLeptons[0].pt(), selLeptons[0].eta(), selLeptons[0].phi()  ); debugText+=buffer;
                                 if(selLeptons.size()>1)sprintf(buffer, " - lep1 Id=%+3i pT=%6.2f, eta=%+6.2f phi=%+6.2f\n", selLeptons[1].pdgId(), selLeptons[1].pt(), selLeptons[1].eta(), selLeptons[1].phi()  ); debugText+=buffer;                                
                               }
                            }

                            if(mt>500){
                              mon.fillHisto( "metNM1",tags,imet.pt(),weight,true);
                            }

                            //pre-VBF control
                            if(njets>=2){
                              LorentzVector dijet=selJets[0].p4()+selJets[1].p4();
                              float deta=fabs(selJets[0].eta()-selJets[1].eta());
                              float dphi=fabs(deltaPhi(selJets[0].phi(),selJets[1].phi()));
                              float pt1(selJets[0].pt()),pt2(selJets[1].pt());
                              mon.fillHisto( "leadjetpt",tags,pt1,weight);
                              mon.fillHisto( "trailerjetpt",tags,pt2,weight);
                              if(pt1>30 && pt2>30){
                                float eta1(selJets[0].eta()),eta2(selJets[1].eta());
                                float fwdEta( fabs(eta1)>fabs(eta2) ? eta1 : eta2);
                                float cenEta( fabs(eta1)>fabs(eta2) ? eta2 : eta1);
                                mon.fillHisto("vbfjeteta", tags,fabs(fwdEta),  weight);
                                mon.fillHisto("vbfjeteta", tags,fabs(cenEta),  weight);                          
                                mon.fillHisto("fwdjeteta",tags,fabs(fwdEta),  weight);
                                mon.fillHisto("cenjeteta",tags,fabs(cenEta),  weight);
                                mon.fillHisto("vbfdetajj",tags,deta,        weight);
                                if(deta>4.0){
                                  mon.fillHisto("vbfmjj",   tags,dijet.mass(),weight,true);
                                  if(dijet.mass()>500){
                                    mon.fillHisto("vbfdphijj",tags,dphi,        weight);
                                    int countJetVeto = 0;
                                    for(size_t ijet=2; ijet<selJets.size(); ijet++){
                                       if((selJets[ijet].eta()<selJets[0].eta() && selJets[ijet].eta()>selJets[1].eta()) ||
                                          (selJets[ijet].eta()>selJets[0].eta() && selJets[ijet].eta()<selJets[1].eta())){
                                          countJetVeto++;
                                       }
                                    }
                                    mon.fillHisto("vbfcjv",tags,countJetVeto,        weight);
                                  }
                                }

                                if(imet.pt()>125){

                                    mon.fillHisto("vbfdetajjfinal",tags,deta,        weight);                          
                                    if(deta>4.0){mon.fillHisto("vbfdphijjfinal",tags,dphi,        weight);}

                                    if(!isMC){
                                       char buffer[1024];                           
                                       sprintf(buffer, " - VBF mjj=%8.2f  dEta=%+6.2f dPhi=%+6.2f\n", dijet.mass(), deta, dphi   ); debugText+=buffer;
                                       sprintf(buffer, " - VBF jet1 pT=%6.2f eta=%+6.2f phi=%+6.2f\n", selJets[0].pt(), selJets[0].eta(), selJets[0].phi()   ); debugText+=buffer;                                
                                       sprintf(buffer, " - VBF jet2 pT=%6.2f eta=%+6.2f phi=%+6.2f\n", selJets[1].pt(), selJets[1].eta(), selJets[1].phi()   ); debugText+=buffer;                                                            
                                    }
                                }

                              }
                            }
                          }
                        }
                      }
                    }        
                  }

                  bool isZsideBand    ( (boson.mass()>40  && boson.mass()<70) || (boson.mass()>110 && boson.mass()<200) );              
                  if(passQt && passThirdLeptonVeto && passMinDphijmet && (boson.mass()>40 && boson.mass()<200)){
                     if(passBtags){                                      
                        if(imet.pt()>50 )mon.fillHisto("zmass_bveto50" , tags,boson.mass(),weight); 
                        if(imet.pt()>80 )mon.fillHisto("zmass_bveto80" , tags,boson.mass(),weight); 
                        if(imet.pt()>125)mon.fillHisto("zmass_bveto125", tags,boson.mass(),weight); 
                        if(passMass){
                            mon.fillHisto( "met_Inbveto",tags,imet.pt(),weight);
                           if(imet.pt()>50 )mon.fillHisto("mt_Inbveto50" , tags,mt,weight); 
                           if(imet.pt()>80 )mon.fillHisto("mt_Inbveto80" , tags,mt,weight); 
                           if(imet.pt()>125)mon.fillHisto("mt_Inbveto125", tags,mt,weight); 
                        }else if(isZsideBand){
                            mon.fillHisto( "met_Outbveto",tags,imet.pt(),weight);
                           if(imet.pt()>50 )mon.fillHisto("mt_Outbveto50" , tags,mt,weight); 
                           if(imet.pt()>80 )mon.fillHisto("mt_Outbveto80" , tags,mt,weight); 
                           if(imet.pt()>125)mon.fillHisto("mt_Outbveto125", tags,mt,weight); 
                        }
                     }else{
                        if(imet.pt()>50 )mon.fillHisto("zmass_btag50" , tags,boson.mass(),weight); 
                        if(imet.pt()>80 )mon.fillHisto("zmass_btag80" , tags,boson.mass(),weight); 
                        if(imet.pt()>125)mon.fillHisto("zmass_btag125", tags,boson.mass(),weight); 
                        if(passMass){
                            mon.fillHisto( "met_Inbtag",tags,imet.pt(),weight);
                           if(imet.pt()>50 )mon.fillHisto("mt_Inbtag50" , tags,mt,weight); 
                           if(imet.pt()>80 )mon.fillHisto("mt_Inbtag80" , tags,mt,weight); 
                           if(imet.pt()>125)mon.fillHisto("mt_Inbtag125", tags,mt,weight); 
                        }else if(isZsideBand){
                            mon.fillHisto( "met_Outbtag",tags,imet.pt(),weight);
                           if(imet.pt()>50 )mon.fillHisto("mt_Outbtag50" , tags,mt,weight); 
                           if(imet.pt()>80 )mon.fillHisto("mt_Outbtag80" , tags,mt,weight); 
                           if(imet.pt()>125)mon.fillHisto("mt_Outbtag125", tags,mt,weight); 
                        }
                     }
                  }

               }


               //if(ivar==0)printf("%9i:%9lli SYST:%30s  Met=%8.3f mT=%8.3f  Weight=%6.2E %i %i %i %i %i\n",  ev.eventAuxiliary().run(), ev.eventAuxiliary().event(), "NOSYST", imet.pt(), higgs::utils::transverseMass(boson,imet,true), weight, passBtags?1:0, passMass?1:0, passQt?1:0, passThirdLeptonVeto?1:0, passMinDphijmet?1:0 ); 

               if(passBtags && passMass && passQt && passThirdLeptonVeto && passMinDphijmet){

                  mon.fillHisto(TString("mtSyst")+varNames[ivar],tags, mt,weight);
                  mon.fillHisto(TString("metSyst")+varNames[ivar],tags, imet.pt(),weight);                    




                  //scan the MET cut and fill the shapes
                  for(unsigned int index=0;index<optim_Cuts1_met.size();index++){             
                     if(imet.pt()>optim_Cuts1_met[index]){
                       for(unsigned int nri=0;nri<NRparams.size();nri++){
                          //Higgs line shape
                          float shapeWeight = weight;   //used for shape dependent weights (avoid overwritting chWeights)
                          double weightToOtherNRI = ( (lShapeWeights[nri][0] * lShapeWeights[nri][1]) / (lShapeWeights[0][0] * lShapeWeights[0][1]) );  //remove weights form nri=0 as those are already in the nominal weight and apply the one for NRI!=0;
                          if(!std::isnan((double)weightToOtherNRI))shapeWeight *= weightToOtherNRI; 

                          if(varNames[ivar]=="_signal_normdown") shapeWeight*=lShapeWeights[nri][2];
                          if(varNames[ivar]=="_signal_lshapedown") shapeWeight*=lShapeWeights[nri][3];
                          if(varNames[ivar]=="_signal_normup"  ) shapeWeight*=lShapeWeights[nri][4];
                          if(varNames[ivar]=="_signal_lshapeup"  ) shapeWeight*=lShapeWeights[nri][5];

                          //if(nri==0 && index==0)printf("%9i:%9lli SYST:%30s  Met=%8.3f mT=%8.3f  Weight=%6.2E\n",  ev.eventAuxiliary().run(), ev.eventAuxiliary().event(), varNames[ivar].Data(), imet.pt(), mt, weight ); 

                          mon.fillHisto(TString("mt_shapes")+NRsuffix[nri]+varNames[ivar],tags,index, mt,shapeWeight);
                          mon.fillHisto(TString("met_shapes")+NRsuffix[nri]+varNames[ivar],tags,index, imet.pt(),shapeWeight);                    
                       }
                    }
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
  //
 
  //scale all events by 1/N to avoid the initial loop to stupidly count the events
  //mon.Scale(1.0/totalNumEvent);

  
  TString terminationCmd = "";
  //save control plots to file
  printf("Results save in local directory and moved to %s\n", outUrl.Data());
  
  //save all to the file
  terminationCmd += TString("mv out.root ") + outUrl + ";";
  TFile *ofile=TFile::Open("out.root", "recreate");
  mon.Write();
  ofile->Close();

  if(!isMC && debugText!=""){ 
     TString outTxtUrl= outUrl + ".txt";    
     terminationCmd += TString("mv out.txt ") + outTxtUrl + ";";
     FILE* outTxtFile = fopen("out.txt", "w");
     fprintf(outTxtFile, "%s", debugText.c_str());
     printf("TextFile URL = %s\n",outTxtUrl.Data());
     if(outTxtFile)fclose(outTxtFile);
  }

  //Now that everything is done, dump the list of lumiBlock that we processed in this job
  if(!isMC){
     terminationCmd += TString("mv out.json ") + ((outUrl.ReplaceAll(".root",""))+".json") + ";";
     goodLumiFilter.FindLumiInFiles(urls);
     goodLumiFilter.DumpToJson("out.json");
  }

  system(terminationCmd.Data());

}  
