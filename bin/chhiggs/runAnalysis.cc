#include <iostream>
#include <boost/shared_ptr.hpp>

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

//Load here all the dataformat that we will need
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
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

  TString suffix=runProcess.getParameter<std::string>("suffix");
  std::vector<std::string> urls=runProcess.getParameter<std::vector<std::string> >("input");
  TString baseDir    = runProcess.getParameter<std::string>("dirName");
  TString url = TString(argv[1]);
  TString outFileUrl(gSystem->BaseName(url));
  outFileUrl.ReplaceAll("_cfg.py","");
  if(mctruthmode!=0) { outFileUrl += "_filt"; outFileUrl += mctruthmode; }
  TString outdir=runProcess.getParameter<std::string>("outdir");
  TString outUrl( outdir );
  gSystem->Exec("mkdir -p " + outUrl);


  bool filterOnlyEE(false), filterOnlyMUMU(false), filterOnlyEMU(false);
  if(!isMC)
    {
      if(url.Contains("DoubleEle")) filterOnlyEE=true;
      if(url.Contains("DoubleMu"))  filterOnlyMUMU=true;
      if(url.Contains("MuEG"))      filterOnlyEMU=true;
    }
  bool isSingleMuPD(!isMC && url.Contains("SingleMu"));  
  bool isV0JetsMC(isMC && (url.Contains("DYJetsToLL_50toInf") || url.Contains("WJets")));
  bool isWGmc(isMC && url.Contains("WG"));
  bool isZGmc(isMC && url.Contains("ZG"));
  bool isMC_GG  = isMC && ( string(url.Data()).find("GG" )  != string::npos);
  bool isMC_VBF = isMC && ( string(url.Data()).find("VBF")  != string::npos);
  bool isMC_125OnShell = isMC && (mctruthmode==521);
  if(isMC_125OnShell) mctruthmode=125;
  bool isMC_ZZ  = isMC && ( string(url.Data()).find("TeV_ZZ")  != string::npos);
  bool isMC_WZ  = isMC && ( string(url.Data()).find("TeV_WZ")  != string::npos);

  TString outTxtUrl= outUrl + "/" + outFileUrl + ".txt";
  FILE* outTxtFile = NULL;
  if(!isMC)outTxtFile = fopen(outTxtUrl.Data(), "w");
  printf("TextFile URL = %s\n",outTxtUrl.Data());

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

  GammaWeightsHandler *gammaWgtHandler=0;
  if(mctruthmode==22 || mctruthmode==111) gammaWgtHandler=new GammaWeightsHandler(runProcess,"",true);

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

  //HIGGS weights and uncertainties
  
  //narrow resonance    
  double cprime = runProcess.getParameter<double>("cprime");
  double brnew  = runProcess.getParameter<double>("brnew");
  std::vector<std::pair<double, double> > NRparams;
  NRparams.push_back(std::make_pair<double,double>(1.0, 0));
//  NRparams.push_back(std::make_pair<double,double>(double(cprime),double(brnew)) );
//  if(mctruthmode==125){
//    NRparams.push_back(std::make_pair<double,double>(5, 0));
//    NRparams.push_back(std::make_pair<double,double>(8, 0));
//    NRparams.push_back(std::make_pair<double,double>(10,0));
//    NRparams.push_back(std::make_pair<double,double>(11,0));
//    NRparams.push_back(std::make_pair<double,double>(12,0));
//    NRparams.push_back(std::make_pair<double,double>(13,0));
//    NRparams.push_back(std::make_pair<double,double>(14,0));
//    NRparams.push_back(std::make_pair<double,double>(15,0));
//    NRparams.push_back(std::make_pair<double,double>(16,0));
//    NRparams.push_back(std::make_pair<double,double>(17,0));
//    NRparams.push_back(std::make_pair<double,double>(18,0));
//    NRparams.push_back(std::make_pair<double,double>(19,0));
//    NRparams.push_back(std::make_pair<double,double>(20,0));
//    NRparams.push_back(std::make_pair<double,double>(22,0));
//    NRparams.push_back(std::make_pair<double,double>(25,0));
//    NRparams.push_back(std::make_pair<double,double>(30,0));
//  }else if(suffix==""){ //consider the other points only when no suffix is being used    
//    for(double cp=0.1;cp<=1.0;cp+=0.1){
//       for(double brn=0.0; brn<=0.5;brn+=0.1){
//          NRparams.push_back(std::make_pair<double,double>((double)cp, (double)brn) );
//    }}
//  }


  std::vector<TGraph *> NRweightsGr;
  std::vector<double> NRweights(NRparams.size());
  std::vector<TString>NRsuffix; for(unsigned int nri=0;nri<NRparams.size();nri++){if(NRparams[nri].first<0 && NRparams[nri].second<0){NRsuffix.push_back(TString(""));}else{char tmp[255];sprintf(tmp,"_cp%3.2f_brn%3.2f",NRparams[nri].first, NRparams[nri].second); NRsuffix.push_back(TString(tmp));} }

  
  //STANDARD MODEL
  double HiggsMass=0; string VBFString = ""; string GGString("");
  TF1 *decayProbPdf=new TF1("relbw","(2*sqrt(2)*[0]*[1]*sqrt(pow([0],2)*(pow([0],2)+pow([1],2)))/(TMath::Pi()*sqrt(pow([0],2)+sqrt(pow([0],2)*(pow([0],2)+pow([1],2))))))/(pow(pow(x,2)-pow([0],2),2)+pow([0]*[1],2))",0,2000);
  if(isMC_GG){  
    size_t GGStringpos =  string(url.Data()).find("GG");
    string StringMass = string(url.Data()).substr(GGStringpos+5,4);  sscanf(StringMass.c_str(),"%lf",&HiggsMass);
    GGString = string(url.Data()).substr(GGStringpos);  
  }else if(isMC_VBF){
    size_t VBFStringpos =  string(url.Data()).find("VBF");
    string StringMass = string(url.Data()).substr(VBFStringpos+6,4);  sscanf(StringMass.c_str(),"%lf",&HiggsMass);
    VBFString = string(url.Data()).substr(VBFStringpos);
  }
  if(mctruthmode==125) HiggsMass=124;
  
  //#######################################
  //####      LINE SHAPE WEIGHTS       ####
  //#######################################
  bool useGenLineShape(true),useGenLineShapeForNR(true);
  TH1 *hGen=0;
  TGraph *hLineShapeNominal=0;
  std::map<std::pair<double,double>, std::vector<TGraph *> > hLineShapeGrVec;  
  TFile *nrLineShapesFile=0;
  if(isMC_GG || isMC_VBF)
    {
      if(mctruthmode==125){
	TString nrLineShapesFileUrl(weightsDir+"/higgs_width_zz2l2nu.root");
	gSystem->ExpandPathName(nrLineShapesFileUrl);
	nrLineShapesFile=TFile::Open(nrLineShapesFileUrl);
      }
      else if(useGenLineShapeForNR){
	TString nrLineShapesFileUrl(weightsDir+"/NR_weightsFromLoic.root");
	gSystem->ExpandPathName(nrLineShapesFileUrl);
	nrLineShapesFile=TFile::Open(nrLineShapesFileUrl);
      }

      TString lineShapeWeightsFileURL(weightsDir+"/");
      lineShapeWeightsFileURL += (isMC_VBF ? "VBFtoHtoZZLineShapes.root" : "GGtoHtoZZLineShapes.root");
      gSystem->ExpandPathName(lineShapeWeightsFileURL);
      TFile *fin=TFile::Open(lineShapeWeightsFileURL);     
      
      TString interferenceShapeWeightsFileUrl(weightsDir+"/");
      interferenceShapeWeightsFileUrl += (isMC_VBF ? "VBFtoHtoZZLineShapesInterference.root" : "GGtoHtoZZLineShapesInterference.root");
      TFile *fin_int=0;
      if(interferenceShapeWeightsFileUrl!="")// && (isMC_GG) 
	{
	  gSystem->ExpandPathName(interferenceShapeWeightsFileUrl);
	  fin_int=TFile::Open(interferenceShapeWeightsFileUrl);
	}
      if(fin) 
	{
	  if(!useGenLineShape) cout << "Line shape weights (and uncertainties) will be applied from " << fin->GetName() << endl;
	  if(fin_int)          cout << "Inteference terms (and uncertainties) will betaken from " << fin_int->GetName() << endl;
	  
	  char dirBuf[100];
	  sprintf(dirBuf,"H%d/",int(HiggsMass));
	  	  
	  hGen                   = (TH1 *) fin->Get(dirBuf+TString("gen")); hGen->SetDirectory(0); hGen->Scale(1./hGen->Integral()); 
	  if(!useGenLineShape)   hLineShapeNominal = new TGraph((TH1 *)fin->Get(dirBuf+TString("cps_shape")));	  
	  else                   hLineShapeNominal = new TGraph(hGen);


	  //CPS shape: if not required set all weights to 1.0
	  TGraph *cpsGr          = (TGraph *) fin->Get(dirBuf+TString("cps"));
	  if(useGenLineShape){
	    for(int ip=0; ip<cpsGr->GetN(); ip++){
	      Double_t x,y; 
	      cpsGr->GetPoint(ip,x,y);
	      cpsGr->SetPoint(ip,x,1.0);
	    }
	  }

	  //interference: if not found set all weights to 1.0, if only nominal then set weights for 100% variation 
	  TGraph *cpspintGr          = (TGraph *) (fin_int!=0? fin_int: fin)->Get(dirBuf+TString("nominal"));
	  if(cpspintGr==0) cpspintGr = (TGraph *) (fin_int!=0? fin_int: fin)->Get(dirBuf+TString("Ratio"));
	  TGraph *cpspint_upGr       = (TGraph *) (fin_int!=0? fin_int: fin)->Get(dirBuf+TString("up"));
	  TGraph *cpspint_downGr     = (TGraph *) (fin_int!=0? fin_int: fin)->Get(dirBuf+TString("down"));
	  if(cpspintGr==0)
	    {
	      cpspintGr = (TGraph *)cpsGr->Clone();
	      for(int ip=0; ip<cpspintGr->GetN(); ip++) { Double_t x,y; cpspintGr->GetPoint(ip,x,y); cpspintGr->SetPoint(ip,x,1); }
	      cpspint_upGr = (TGraph *) cpspintGr->Clone();
	      cpspint_downGr=(TGraph *) cpspintGr->Clone();
	    }
	  else if(cpspint_upGr==0 && cpspint_downGr==0)
	    {
	      cpspint_upGr=(TGraph *)cpspintGr->Clone();
	      cpspint_downGr=(TGraph *)cpspintGr->Clone();
	      for(int ip=0; ip<cpspintGr->GetN(); ip++) { 
		Double_t x,y; 
		cpspintGr->GetPoint(ip,x,y);
		cpspint_downGr->SetPoint(ip,x,1.0);
		float yDiff(fabs(1-y));
		float yMirror(1-2*yDiff);
		if(y>1)        yMirror=1+2*yDiff;
		if(yMirror<0)  yMirror=0;
		if(yMirror>10) yMirror=10;
		cpspint_upGr->SetPoint(ip,x,yMirror);
	      }
	    }
	  
	  //loop over possible scenarios: SM or BSM
	  for(size_t nri=0; nri<NRparams.size(); nri++)
	    {
	      //recompute weights depending on the scenario (SM or BSM)
	      TGraph *shapeWgtsGr      = new TGraph; shapeWgtsGr->SetName("shapeWgts_"+ NRsuffix[nri]);          float shapeNorm(0);
	      TGraph *shapeWgts_upGr   = new TGraph; shapeWgts_upGr->SetName("shapeWgtsUp_"+ NRsuffix[nri]);     float shapeUpNorm(0);
	      TGraph *shapeWgts_downGr = new TGraph; shapeWgts_downGr->SetName("shapeWgtsDown_"+ NRsuffix[nri]); float shapeDownNorm(0);

              TGraph* nrWgtGr    =NULL;
              TGraph* nrWgtUpGr  =NULL;
              TGraph* nrWgtDownGr=NULL;              

	      Float_t hySum(0);
	      for(int ip=1; ip<=hGen->GetXaxis()->GetNbins(); ip++)
		{
		  Double_t hmass    = hGen->GetBinCenter(ip);
		  Double_t hy       = hGen->GetBinContent(ip);
		  
		  Double_t shapeWgt(1.0),shapeWgtUp(1.0),shapeWgtDown(1.0);
		  if(NRparams[nri].first<0)
		    {
		      shapeWgt     = cpsGr->Eval(hmass) * cpspintGr->Eval(hmass);
		      shapeWgtUp   = cpsGr->Eval(hmass) * cpspint_upGr->Eval(hmass);
		      shapeWgtDown = cpsGr->Eval(hmass) * cpspint_downGr->Eval(hmass);
		    }
		  else if(mctruthmode==125){
		    TString var("");
		    if(url.Contains("ScaleUp"))   var="up";
		    if(url.Contains("ScaleDown")) var="down";
		    Double_t nrWgt = higgs::utils::weightToH125Interference(hmass,NRparams[nri].first,nrLineShapesFile,var);
		    shapeWgt       = cpsGr->Eval(hmass) * nrWgt;
		    shapeWgtUp     = shapeWgt;
		    shapeWgtDown   = shapeWgt;
		  }
		  else
		    {
                      if(!nrWgtGr)    nrWgtGr     = higgs::utils::weightNarrowResonnance(VBFString,HiggsMass, hmass, NRparams[nri].first, NRparams[nri].second, hLineShapeNominal,decayProbPdf,nrLineShapesFile);
                      if(!nrWgtUpGr)  nrWgtUpGr   = higgs::utils::weightNarrowResonnance(VBFString,HiggsMass, hmass, NRparams[nri].first, NRparams[nri].second, hLineShapeNominal,decayProbPdf,nrLineShapesFile,"_up");          
                      if(!nrWgtDownGr)nrWgtDownGr = higgs::utils::weightNarrowResonnance(VBFString,HiggsMass, hmass, NRparams[nri].first, NRparams[nri].second, hLineShapeNominal,decayProbPdf,nrLineShapesFile,"_down"); 

                    shapeWgt       = nrWgtGr    ->Eval(hmass);
                    shapeWgtUp     = nrWgtUpGr  ->Eval(hmass);
                    shapeWgtDown   = nrWgtDownGr->Eval(hmass);

//		      Double_t cpsWgt=cpsGr->Eval(hmass);
//		      shapeWgt       = cpsWgt * std::max(0.0, nrWgtGr    ->Eval(hmass));
//		      shapeWgtUp     = cpsWgt * std::max(0.0, nrWgtUpGr  ->Eval(hmass));
//		      shapeWgtDown   = cpsWgt * std::max(0.0, nrWgtDownGr->Eval(hmass));
		    }
		  
		  shapeWgtsGr->SetPoint(shapeWgtsGr->GetN(),           hmass, shapeWgt);       shapeNorm     += shapeWgt*hy;
		  shapeWgts_upGr->SetPoint(shapeWgts_upGr->GetN(),     hmass, shapeWgtUp);     shapeUpNorm   += shapeWgtUp*hy;
		  shapeWgts_downGr->SetPoint(shapeWgts_downGr->GetN(), hmass, shapeWgtDown);   shapeDownNorm += shapeWgtDown*hy;
		}
	      
	      if(mctruthmode!=125)
	      {
		cout << "C'=" << NRparams[nri].first << " " << hySum << " " << shapeNorm << endl;
		if(hySum>0){
		  shapeNorm     /= hySum;
		  shapeUpNorm   /= hySum;
		  shapeDownNorm /= hySum;
		}
		//fix possible normalization issues
		cout << "C'=" << NRparams[nri].first << " BRnew=" << NRparams[nri].second << " shape wgts will be re-scaled to preserve unitarity with: "
		     << " nominal=" << shapeNorm
		     << " up     =" << shapeUpNorm
		     << " down   =" << shapeDownNorm 
		     << endl;
		for(Int_t ip=0; ip<shapeWgtsGr->GetN(); ip++)
		  {
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
	      std::vector<TGraph *> inrWgts;
	      inrWgts.push_back( shapeWgtsGr      );
	      inrWgts.push_back( shapeWgts_upGr   );
	      inrWgts.push_back( shapeWgts_downGr );
	      hLineShapeGrVec[ NRparams[nri] ] = inrWgts;
	    }
	  
	  //close files
	  fin->Close();
	  delete fin;
	  if(fin_int){
	    fin_int->Close();
	    delete fin_int;
	  }
	}
    }


  //##############################################
  //########    INITIATING HISTOGRAMS     ########
  //##############################################
  SmartSelectionMonitor mon;

  //generator level control : add an underflow entry to make sure the histo is kept
  ((TH1F*)mon.addHistogram( new TH1F( "higgsMass_raw",     ";Higgs Mass [GeV];Events", 500,0,1500) ))->Fill(-1.0,0.0001);
  ((TH1F*)mon.addHistogram( new TH1F( "higgsMass_cpspint", ";Higgs Mass [GeV];Events", 500,0,1500) ))->Fill(-1.0,0.0001);
  for(unsigned int nri=0;nri<NRparams.size();nri++){ 
    ((TH1F*)mon.addHistogram( new TH1F( "higgsMass_4nr"+NRsuffix[nri] , ";Higgs Mass;Events [GeV]", 500,0,1500) ))->Fill(-1.0,0.0001);
  }

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

  //pu control
  mon.addHistogram( new TH1F( "nvtx",";Vertices;Events",50,0,50) ); 
  mon.addHistogram( new TH1F( "nvtxraw",";Vertices;Events",50,0,50) ); 
  mon.addHistogram( new TH1F( "rho",";#rho;Events",50,0,25) ); 

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

  //lepton control
  mon.addHistogram( new TH1F( "leadpt",     ";Transverse momentum [GeV];Events", 50,0,500) );
  mon.addHistogram( new TH1F( "leadeta",    ";Pseudo-rapidity;Events", 50,0,2.6) );
  mon.addHistogram( new TH1F( "trailerpt",  ";Transverse momentum [GeV];Events", 50,0,500) );
  mon.addHistogram( new TH1F( "trailereta", ";Pseudo-rapidity;Events", 50,0,2.6) );
  mon.addHistogram( new TH1F( "zy",         ";Rapidity;Events", 50,0,3) );
  mon.addHistogram( new TH1F( "zmass",      ";Mass [GeV];Events", 100,40,250) );
  mon.addHistogram( new TH1F( "zpt",        ";Transverse momentum [GeV];Events",100,0,1500));
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
  mon.addHistogram( new TH1F( "met",          ";Missing transverse energy [GeV];Events",nmetAxis-1,metaxis) ); //50,0,1000) );
  mon.addHistogram( new TH1F( "metNM1",        ";Missing transverse energy [GeV];Events",nmetAxis-1,metaxis) ); //50,0,1000) );
  Double_t mtaxis[]={100,120,140,160,180,200,220,240,260,280,300,325,350,375,400,450,500,600,700,800,900,1000,2000};
  Int_t nmtAxis=sizeof(mtaxis)/sizeof(Double_t);
  mon.addHistogram( new TH1F( "mt"  ,         ";Transverse mass;Events",nmtAxis-1,mtaxis) );
  mon.addHistogram( new TH1F( "mtNM1"  ,       ";Transverse mass;Events",nmtAxis-1,mtaxis) );
  mon.addHistogram( new TH1F( "mtresponse",   ";Transverse mass response;Events", 100,0,2) );
  mon.addHistogram( new TH1F( "mtcheckpoint"  ,         ";Transverse mass [GeV];Events",160,150,1750) );
  mon.addHistogram( new TH1F( "metcheckpoint" ,         ";Missing transverse energy [GeV];Events",100,0,500) );


  //
  // STATISTICAL ANALYSIS
  //
  std::vector<double> optim_Cuts1_met; 
  for(double met=50;met<140;met+=5) {  optim_Cuts1_met    .push_back(met);  }
  TProfile* Hoptim_cuts1_met     =  (TProfile*) mon.addHistogram( new TProfile ("optim_cut1_met"    , ";cut index;met"    ,optim_Cuts1_met.size(),0,optim_Cuts1_met.size()) ) ;
  mon.addHistogram( new TH1F ("metcount"    , ";E_{T}^{miss} cut [GeV];Total events"    ,optim_Cuts1_met.size(),0,optim_Cuts1_met.size()) );
  for(unsigned int index=0;index<optim_Cuts1_met.size();index++){ Hoptim_cuts1_met    ->Fill(index, optim_Cuts1_met[index]);  }
  TH1F* Hoptim_systs     =  (TH1F*) mon.addHistogram( new TH1F ("optim_systs"    , ";syst;", nvarsToInclude,0,nvarsToInclude) ) ;
  for(size_t ivar=0; ivar<nvarsToInclude; ivar++)
    {
      Hoptim_systs->GetXaxis()->SetBinLabel(ivar+1, varNames[ivar]);
      
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



     
  //##############################################
  //######## GET READY FOR THE EVENT LOOP ########
  //##############################################

  fwlite::ChainEvent ev(urls);
  const size_t totalEntries= ev.size();

  //MC normalization (to 1/pb)
  double xsecWeight = xsec/totalEntries;
  if(!isMC) xsecWeight=1.0;

  //jet energy scale and uncertainties 
  TString jecDir = runProcess.getParameter<std::string>("jecDir");
  gSystem->ExpandPathName(jecDir);
  FactorizedJetCorrector *jesCor        = utils::cmssw::getJetCorrector(jecDir,isMC);
  JetCorrectionUncertainty *totalJESUnc = new JetCorrectionUncertainty((jecDir+"/MC_Uncertainty_AK5PFchs.txt").Data());
  
  //muon energy scale and uncertainties
  MuScleFitCorrector *muCor=getMuonCorrector(jecDir,url);

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
          utils::getMCPileupDistributionFromMiniAOD(ev,dataPileupDistribution.size(), mcPileupDistribution);
          while(mcPileupDistribution.size()<dataPileupDistribution.size())  mcPileupDistribution.push_back(0.0);
          while(mcPileupDistribution.size()>dataPileupDistribution.size())dataPileupDistribution.push_back(0.0);
          gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE
          LumiWeights = new edm::LumiReWeighting(mcPileupDistribution,dataPileupDistribution);
          PuShifters=utils::cmssw::getPUshifters(dataPileupDistribution,0.05);
          utils::getPileupNormalization(mcPileupDistribution, PUNorm, LumiWeights, PuShifters);
  }
 
  gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE


  higgs::utils::EventCategory eventCategoryInst(higgs::utils::EventCategory::EXCLUSIVE2JETSVBF); //jet(0,>=1)+vbf binning


  //##############################################
  //########           EVENT LOOP         ########
  //##############################################
  //loop on all the events
  printf("Progressing Bar     :0%%       20%%       40%%       60%%       80%%       100%%\n");
  printf("Scanning the ntuple :");
  int treeStep(totalEntries/50);
  //DuplicatesChecker duplicatesChecker;
  //int nDuplicates(0);
  for( size_t iev=0; iev<totalEntries; iev++){
      if(iev%treeStep==0){printf(".");fflush(stdout);}

       //##############################################   EVENT LOOP STARTS   ##############################################
       ev.to(iev); //load the event content from the EDM file
       //if(!isMC && duplicatesChecker.isDuplicate( ev.run, ev.lumi, ev.event) ) { nDuplicates++; continue; }

       //apply trigger and require compatibilitiy of the event with the PD
       edm::TriggerResultsByName tr = ev.triggerResultsByName("HLT");
       if(!tr.isValid())return false;

      bool eeTrigger          = utils::passTriggerPatterns(tr, "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
      bool muTrigger          = utils::passTriggerPatterns(tr, "HLT_IsoMu24_eta2p1_v*");
      bool mumuTrigger        = utils::passTriggerPatterns(tr, "HLT_Mu17_Mu8_v*", "HLT_Mu17_TkMu8_v*"); 
      bool emuTrigger         = utils::passTriggerPatterns(tr, "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*", "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
      if(filterOnlyEE)   { mumuTrigger=false; emuTrigger=false;  }
      if(filterOnlyMUMU) { eeTrigger=false;   emuTrigger=false;  }
      if(isSingleMuPD)   { eeTrigger=false;   emuTrigger=false;  if( muTrigger && !mumuTrigger) mumuTrigger=true; else mumuTrigger=false; }
      if(filterOnlyEMU)  { eeTrigger=false;   mumuTrigger=false; }

      bool hasPhotonTrigger(false);
      float triggerPrescale(1.0),triggerThreshold(0);
      bool runPhotonSelection(mctruthmode==22 || mctruthmode==111);
      if(runPhotonSelection){
	  eeTrigger=false; mumuTrigger=false;

          std::string successfulPath="";
          if(     utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_v*")){ hasPhotonTrigger=true; triggerThreshold=92; }
          else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_v*")){ hasPhotonTrigger=true; triggerThreshold=77; }
          else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_v*")){ hasPhotonTrigger=true; triggerThreshold=50; }
          else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_v*")){ hasPhotonTrigger=true; triggerThreshold=36; }
          
          if(successfulPath!=""){ //get the prescale associated to it
             fwlite::Handle< pat::PackedTriggerPrescales > prescalesHandle;
             prescalesHandle.getByLabel(ev, "patTrigger");
             pat::PackedTriggerPrescales prescales = *prescalesHandle;
             const edm::TriggerResults& trResults =  prescales.triggerResults();
             prescales.setTriggerNames( ev.triggerNames(trResults) );
             triggerPrescale = prescales.getPrescaleForName(successfulPath);
          }
      }
      if(!(eeTrigger || muTrigger || mumuTrigger || emuTrigger || hasPhotonTrigger))continue;  //ONLY RUN ON THE EVENTS THAT PASS OUR TRIGGERS

       //##############################################   EVENT PASSED THE TRIGGER   #######################################


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

       pat::TauCollection taus;
       fwlite::Handle< pat::TauCollection > tausHandle;
       tausHandle.getByLabel(ev, "slimmedTaus");
       if(tausHandle.isValid()){ taus = *tausHandle;}

      if(isV0JetsMC){
         fwlite::Handle< LHEEventProduct > lheEPHandle;
         lheEPHandle.getByLabel(ev, "externalLHEProducer");

 	 mon.fillHisto("nup","",lheEPHandle->hepeup().NUP,1);
 	 if(lheEPHandle->hepeup().NUP>5) continue;
	 mon.fillHisto("nupfilt","",lheEPHandle->hepeup().NUP,1);
      }



      //MC crap for photon studies
      if(hasPhotonTrigger && (isWGmc || isZGmc)){
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


      //
      // DERIVE WEIGHTS TO APPLY TO SAMPLE
      //

       //pileup weight
       float weight = 1.0;
       double TotalWeight_plus = 1.0;
       double TotalWeight_minus = 1.0;
       float puWeight(1.0);

       if(isMC){          
          int ngenITpu = 0;

          fwlite::Handle< std::vector<PileupSummaryInfo> > puInfoH;
          puInfoH.getByLabel(ev, "addPileupInfo");
          for(std::vector<PileupSummaryInfo>::const_iterator it = puInfoH->begin(); it != puInfoH->end(); it++){
             if(it->getBunchCrossing()==0)      { ngenITpu += it->getPU_NumInteractions(); }
          }

          puWeight          = LumiWeights->weight(ngenITpu) * PUNorm[0];
          weight            = xsecWeight*puWeight;
          TotalWeight_plus  = PuShifters[utils::cmssw::PUUP  ]->Eval(ngenITpu) * (PUNorm[2]/PUNorm[0]);
          TotalWeight_minus = PuShifters[utils::cmssw::PUDOWN]->Eval(ngenITpu) * (PUNorm[1]/PUNorm[0]);
      }


      //Higgs specific weights
      float lShapeWeights[3]={1.0,1.0,1.0};
      for(unsigned int nri=0;nri<NRparams.size();nri++){NRweights[nri] = 1.0;}
      if(isMC){

	LorentzVector higgs(0,0,0,0);
	LorentzVector totLeptons(0,0,0,0);
	for(size_t igen=0; igen<gen.size(); igen++){
	  if(gen[igen].status()!=3) continue;
	  if(abs(gen[igen].pdgId())>=11 && abs(gen[igen].pdgId())<=16) totLeptons += gen[igen].p4();
	  if(gen[igen].pdgId()==25)                                      higgs=gen[igen].p4();
	}
	if(mctruthmode==125) {
	  higgs=totLeptons;
	  if(isMC_125OnShell && higgs.mass()>180) continue;
	  if(!isMC_125OnShell && higgs.mass()<=180) continue;
	}
	float shapeWeight(1.0);
        if((isMC_VBF || isMC_GG) && higgs.pt()>0){
	  {
	    //Line shape weights 
	    if(isMC_VBF || isMC_GG)
	      {
		std::vector<TGraph *> nominalShapeWgtGr=hLineShapeGrVec.begin()->second;
		for(size_t iwgt=0; iwgt<nominalShapeWgtGr.size(); iwgt++)
		  {
		    if(nominalShapeWgtGr[iwgt]==0) continue;
		    lShapeWeights[iwgt]=nominalShapeWgtGr[iwgt]->Eval(higgs.mass());
		  }
	      }
	    shapeWeight   = lShapeWeights[0];
	    
	    //control SM line shape
	    mon.fillHisto("higgsMass_raw",    "", higgs.mass(), puWeight);
	    mon.fillHisto("higgsMass_cpspint","", higgs.mass(), puWeight * shapeWeight);
	    
	    //compute weight correction for narrow resonnance
	    for(unsigned int nri=0;nri<NRparams.size();nri++){ 
	      if(NRparams[nri].first<0) continue;
	      std::vector<TGraph *> shapeWgtGr = hLineShapeGrVec[NRparams[nri] ];
	      NRweights[nri] = shapeWgtGr[0]->Eval(higgs.mass()); 
	      float iweight = puWeight * NRweights[nri];
	      mon.fillHisto(TString("higgsMass_4nr")+NRsuffix[nri], "", higgs.mass(), iweight );
	    }  
	  }
	}
  
	//final event weight
	weight = puWeight * shapeWeight;
      }

      //
      //
      // BELOW FOLLOWS THE ANALYSIS OF THE MAIN SELECTION WITH N-1 PLOTS
      //
      //

      //
      // photon selection
      //
      pat::PhotonCollection selPhotons;
      if(runPhotonSelection)
	{
	  //filter out number of prompt photons to avoid double counting
	 //  int ngenpho(0);
	  // 	  for(size_t igen=0; igen<gen.size(); igen++)
	  // 	    {
	  // 	      if(gen[igen].get("id")!=22 || gen[igen].get("status")!=1) continue;
	  // 	      float lxy=gen[igen].getVal("lxy");
	  // 	      if(lxy>0) continue;
	  // 	      ngenpho++;
	  //  }
	  //if(mctruthmode==111 && ngenpho>0) continue;
	  //if(mctruthmode==22 && ngenpho==0) continue;

	  //select the photons
	  for(size_t ipho=0; ipho<photons.size(); ipho++)
	    {
	      double pt=photons[ipho].pt();
	      double eta=photons[ipho].superCluster()->eta();

	      //if systematics are active loosen the selection to the medium working point
              //std::vector<std::pair<std::string, bool> > phid =  photons[ipho].photonIDs() ;   for(size_t PHIDI = 0 ; PHIDI<phid.size(); PHIDI++){ printf("%i   %s  %i\n", (int)PHIDI, phid[PHIDI].first.c_str(), phid[PHIDI].second?1:0);  }  //print available photon ID
              bool hasTightPhotonId = true;// photons[ipho].photonID("pidTight")>1.0;//FIXME not sure what the cut should be
	      bool passId = (photons[ipho].r9()>=0.9 && hasTightPhotonId); 

  	      //isolation 
	      //FIXME: Currently using straight the iso from PAT, not sure it's equivalent to what we were doing (just below)
	      bool passIso(true);    
	      float relIso = photons[ipho].particleIso() / photons[ipho].pt();  
              if(relIso<0.4)passIso = true;
//	      double gIso    = photons[ipho].getVal("gIso03");
//	      double gArea   = utils::cmssw::getEffectiveArea(22,eta,3,"gIso");	      
//	      double chIso   = photons[ipho].getVal("chIso03");
//	      double chArea  = utils::cmssw::getEffectiveArea(22,eta,3,"chIso");
//	      double nhIso   = photons[ipho].getVal("nhIso03");
//	      double nhArea  = utils::cmssw::getEffectiveArea(22,eta,3,"nhIso");
//	      passIso &= (TMath::Max(chIso-chArea*rho,0.0) < 0.7); 
//	      passIso &= (TMath::Max(nhIso-nhArea*rho,0.0) < 0.4+0.04*pt); 
//	      passIso &= (TMath::Max(gIso-gArea*rho,  0.0) < 0.5+0.005*pt); 
      
	      //select the photon
	      if(pt<triggerThreshold || fabs(eta)>1.4442 ) continue;
	      if(!passId) continue;
	      if(!passIso) continue; 
	      selPhotons.push_back(photons[ipho]);
	    }
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

	  //isolation
	  passIso = lid==11?patUtils::passIso(leptons[ilep].el,  patUtils::llvvElecIso::Tight) : patUtils::passIso(leptons[ilep].mu,  patUtils::llvvMuonIso::Tight);
	  
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
      int njets(0),nbtags(0),nbtagsJP(0);
      float mindphijmet(9999.);
      for(size_t ijet=0; ijet<jets.size(); ijet++){
	  if(jets[ijet].pt()<15 || fabs(jets[ijet].eta())>4.7 ) continue;

	  //mc truth for this jet
	  const reco::GenJet* genJet=jets[ijet].genJet();
	  TString jetType( genJet && genJet->pt()>0 ? "truejetsid" : "pujetsid" );
	  
	  //cross-clean with selected leptons and photons
	  float minDRlj(9999.),minDRlg(9999.);
          for(size_t ilep=0; ilep<selLeptons.size(); ilep++)
            minDRlj = TMath::Min( minDRlj, deltaR(jets[ijet],selLeptons[ilep]) );
	  for(size_t ipho=0; ipho<selPhotons.size(); ipho++)
	    minDRlg = TMath::Min( minDRlg, deltaR(jets[ijet],selPhotons[ipho]) );
	  if(minDRlj<0.4 || minDRlg<0.4) continue;
	  
	  //jet id
	  bool passPFloose = true; //FIXME --> Need to be updated according to te latest recipe;
	  float PUDiscriminant = jets[ijet].userFloat("pileupJetId:fullDiscriminant");
	  bool passLooseSimplePuId = true; //FIXME --> Need to be updated according to the latest recipe
	  if(jets[ijet].pt()>30){
	      mon.fillHisto(jetType,"",fabs(jets[ijet].eta()),0);
	      if(passPFloose)                        mon.fillHisto(jetType,"",fabs(jets[ijet].eta()),1);
	      if(passLooseSimplePuId)                mon.fillHisto(jetType,"",fabs(jets[ijet].eta()),2);
	      if(passPFloose && passLooseSimplePuId) mon.fillHisto(jetType,"",fabs(jets[ijet].eta()),3);
	  }
	  if(!passPFloose || !passLooseSimplePuId) continue;
	  selJets.push_back(jets[ijet]);
	  if(jets[ijet].pt()>30) {
	    njets++;
	    float dphijmet=fabs(deltaPhi(met.phi(), jets[ijet].phi()));
	    if(dphijmet<mindphijmet) mindphijmet=dphijmet;
	    if(fabs(jets[ijet].eta())<2.5){
	      bool hasCSVtag(jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags")>0.405);
	      //update according to the SF measured by BTV
	      if(isMC){
		  int flavId=jets[ijet].partonFlavour();
		  if(abs(flavId)==5)        btsfutil.modifyBTagsWithSF(hasCSVtag,sfb,beff);
		  else if(abs(flavId)==4)   btsfutil.modifyBTagsWithSF(hasCSVtag,sfb/5,beff);
		  else		            btsfutil.modifyBTagsWithSF(hasCSVtag,sfl,leff);
              }
	      nbtags   += hasCSVtag;
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
	ntaus++;
      }
      std::sort(selTaus.begin(), selTaus.end(), utils::sort_CandidatesByPt);
      

      //
      // ASSIGN CHANNEL
      //
      std::vector<TString> chTags;
      int dilId(1);
      LorentzVector boson(0,0,0,0);
      if(!runPhotonSelection && selLeptons.size()==2)
	{
 	  for(size_t ilep=0; ilep<2; ilep++)
	    {
	      dilId *= selLeptons[ilep].pdgId();
	      int id(abs(selLeptons[ilep].pdgId()));
	      weight *= isMC ? lepEff.getLeptonEfficiency( selLeptons[ilep].pt(), selLeptons[ilep].eta(), id,  id ==11 ? "loose" : "loose" ).first : 1.0;
	      boson += selLeptons[ilep].p4();
	    }
     
	  //check the channel
	  if( abs(dilId)==121 && eeTrigger)   chTags.push_back("ee");
	  if( abs(dilId)==169 && mumuTrigger) chTags.push_back("mumu"); 
	  if( abs(dilId)==143 && emuTrigger) chTags.push_back("emu"); 
	}
      else{
	if(hasPhotonTrigger && selPhotons.size()) {
	  dilId=22;
	  chTags.push_back("ee");
	  chTags.push_back("mumu");
	  boson = selPhotons[0].p4();
	  weight *= triggerPrescale;
	}
      }

      TString evCat=eventCategoryInst.GetCategory(selJets,boson);
      std::vector<TString> tags(1,"all");
      for(size_t ich=0; ich<chTags.size(); ich++){
	tags.push_back( chTags[ich] );
	tags.push_back( chTags[ich]+evCat );
      }

      mon.fillHisto("eventflow",  tags,0,weight);
      if(chTags.size()==0) continue;

      //
      // BASELINE SELECTION
      //
      bool passMass(fabs(boson.mass()-91)<15);
      bool passQt(boson.pt()>55);
      bool passThirdLeptonVeto( selLeptons.size()==2 && extraLeptons.size()==0 );
      bool passBtags(nbtags==0);
      bool passMinDphijmet( njets==0 || mindphijmet>0.5);
      if(runPhotonSelection)
	{
	  passMass=hasPhotonTrigger;
	  passThirdLeptonVeto=(selLeptons.size()==0 && extraLeptons.size()==0);
	}

      mon.fillHisto("eventflow",  tags,1,weight);
      mon.fillHisto("nvtxraw",  tags,vtx.size(),weight/puWeight);
      mon.fillHisto("nvtx",  tags,vtx.size(),weight);
      mon.fillHisto("rho",  tags,rho,weight);
      if(!runPhotonSelection){
	mon.fillHisto("leadpt",      tags,selLeptons[0].pt(),weight); 
	mon.fillHisto("trailerpt",   tags,selLeptons[1].pt(),weight); 
	mon.fillHisto("leadeta",     tags,fabs(selLeptons[0].eta()),weight); 
	mon.fillHisto("trailereta",  tags,fabs(selLeptons[1].eta()),weight); 

	mon.fillHisto("ntaus", tags, ntaus,weight);
	if(ntaus>0) mon.fillHisto("leadtaupt", tags, selTaus[0].pt(),weight);
      }
      mon.fillHisto("zmass", tags,boson.mass(),weight); 
      mon.fillHisto("zy",    tags,fabs(boson.Rapidity()),weight); 

      if(passMass){

	mon.fillHisto("eventflow",tags, 2,weight);
	mon.fillHisto("zpt",      tags, boson.pt(),weight);

	//these two are used to reweight photon -> Z, the 3rd is a control
	mon.fillHisto("qt",       tags, boson.pt(),weight,true); 
	mon.fillHisto("qtraw",    tags, boson.pt(),weight/triggerPrescale,true); 

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

	      float csv(selJets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"));
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

	      if(gammaWgtHandler!=0) {
		float lastPhotonWeight(1.0), totalPhotonWeight(0.0);
		LorentzVector lastMassiveBoson(boson);
		for(size_t itag=0; itag<tags.size(); itag++){
		  size_t idx(tags.size()-itag-1);
		  std::vector<Float_t> photonVars;
		  photonVars.push_back(boson.pt());
		  //photonVars.push_back(met.pt()/boson.pt());
		  float photonWeight=gammaWgtHandler->getWeightFor(photonVars,tags[idx]);
		  if(tags[idx]=="all")       { 
		    photonWeights[idx]=(totalPhotonWeight==0? 1.0:totalPhotonWeight); 
		  }
		  else if(photonWeight!=1.0) { 
		    photonWeights[idx]=photonWeight; 
		    massiveBoson[idx]=gammaWgtHandler->getMassiveP4(boson,tags[idx]);
		    totalPhotonWeight+=photonWeight; 
		    lastPhotonWeight=photonWeight; 
		    lastMassiveBoson=massiveBoson[idx];
		  }
		  else                       { 
		    photonWeights[idx]=lastPhotonWeight; 
		    massiveBoson[idx]=lastMassiveBoson;
		  }
		}
	      }


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
		  mon.fillHisto( "met",icat,met.pt(),iweight,true);
		  mon.fillHisto( "balance",icat,met.pt()/iboson.pt(),iweight);
		  TVector2 met2(met.px(),met.py());
		  TVector2 boson2(iboson.px(), iboson.py());
		  double axialMet(boson2*met2); axialMet/=-iboson.pt();
		  mon.fillHisto( "axialmet",icat,axialMet,iweight);
		  double mt=higgs::utils::transverseMass(iboson,met,true);
		  mon.fillHisto( "mt",icat,mt,iweight,true);
		  
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
	
	//Higgs line shape
	if(varNames[ivar].Contains("lshape"))
	  {
	    size_t idx( varNames[ivar].EndsWith("up") ? 1 : 2);
	    float shapeReWeight = lShapeWeights[idx];
	    if(lShapeWeights[0]==0) shapeReWeight=0;
	    else                    shapeReWeight /= lShapeWeights[0];
	    iweight                 *= shapeReWeight;
	  }	 

	//recompute MET/MT if JES/JER was varied
	LorentzVector    zvv = mets[0].p4();
//FIXME
//	if(varyJesUp)    zvv = mets[0].shiftedP4(pat::MET::METUncertainty::JetEnUp);
//	if(varyJesDown)  zvv = mets[0].shiftedP4(pat::MET::METUncertainty::JetEnDown);
//	if(varyJerUp)    zvv = mets[0].shiftedP4(pat::MET::METUncertainty::JetResUp);
//	if(varyJerDown)  zvv = mets[0].shiftedP4(pat::MET::METUncertainty::JetResDown);
//	if(varyUmetUp)   zvv = mets[0].shiftedP4(pat::MET::METUncertainty::UnclusteredEnUp);
//	if(varyUmetDown) zvv = mets[0].shiftedP4(pat::MET::METUncertainty::UnclusteredEnDown);
//	if(varyLesUp)    zvv = met[utils::cmssw::LESUP]; //FIXME  must vary all leptons separately: MuonEnUp/MuonEnDown/ElectronEnUp/ElectronEnDown/TauEnUp/TauEnDown
//	if(varyLesDown)  zvv = met[utils::cmssw::LESDOWN];

        pat::JetCollection tightVarJets;
	bool passLocalBveto(passBtags);
 	for(size_t ijet=0; ijet<jets.size(); ijet++){

	  float eta=jets[ijet].eta();
	  if( fabs(eta)>4.7 ) continue;
	  float pt=jets[ijet].pt();
          //FIXME
//	  if(varyJesUp)    pt=jets[ijet].getVal("jesup");
//	  if(varyJesDown)  pt=jets[ijet].getVal("jesdown");
//	  if(varyJerUp)    pt=jets[ijet].getVal("jerup");
//	  if(varyJerDown)  pt=jets[ijet].getVal("jerdown");
	  if(pt<30) continue;

	  //cross-clean with selected leptons and photons
	  double minDRlj(9999.),minDRlg(9999.);
          for(size_t ilep=0; ilep<selLeptons.size(); ilep++)
            minDRlj = TMath::Min( minDRlj, deltaR(jets[ijet].p4(),selLeptons[ilep].p4()) );
	  for(size_t ipho=0; ipho<selPhotons.size(); ipho++)
	    minDRlg = TMath::Min( minDRlg, deltaR(jets[ijet].p4(),selPhotons[ipho].p4()) );
	  if(minDRlj<0.4 || minDRlg<0.4) continue;
	  
	  //jet id
	  bool passPFloose = true;//FIXME
	  int simplePuId = true;//FIXME
	  bool passLooseSimplePuId = true;//FIXME
	  if(!passPFloose || !passLooseSimplePuId) continue;
	 
	  //jet is selected
	  tightVarJets.push_back(jets[ijet]);

	  //check b-tag
	  if(pt<30 || fabs(eta)>2.5) continue;
	  if(!isMC) continue;
	  if(!varyBtagUp && !varyBtagDown) continue;
	  int flavId=jets[ijet].partonFlavour();
	  bool hasCSVtag (jets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags")>0.405);
 	  if(varyBtagUp) {
	    if(abs(flavId)==5)        btsfutil.modifyBTagsWithSF(hasCSVtag,sfb+sfbunc,beff);
	    else if(abs(flavId)==4)   btsfutil.modifyBTagsWithSF(hasCSVtag,sfb/5+2*sfbunc,beff);
	    else		      btsfutil.modifyBTagsWithSF(hasCSVtag,sfl+sflunc,leff);
	  }
 	  else if(varyBtagDown) {
	    if(abs(flavId)==5)        btsfutil.modifyBTagsWithSF(hasCSVtag,sfb-sfbunc,beff);
	    else if(abs(flavId)==4)   btsfutil.modifyBTagsWithSF(hasCSVtag,sfb/5-2*sfbunc,beff);
	    else		      btsfutil.modifyBTagsWithSF(hasCSVtag,sfl-sflunc,leff);
 	  }
	  passLocalBveto |= hasCSVtag;
 	}
	
	bool isZsideBand    ( (boson.mass()>40  && boson.mass()<70) || (boson.mass()>110 && boson.mass()<200) );
	bool isZsideBandPlus( (boson.mass()>110 && boson.mass()<200) );
 	bool passPreselection                 (passMass && passQt && passThirdLeptonVeto && passMinDphijmet && passLocalBveto);
 	bool passPreselectionMbvetoMzmass     (            passQt && passThirdLeptonVeto && passMinDphijmet                  );          
	
 	//re-assign the event category to take migrations into account
 	TString evCat  = eventCategoryInst.GetCategory(tightVarJets,boson);
	for(size_t ich=0; ich<chTags.size(); ich++){
	  
	  TString tags_full=chTags[ich]+evCat;
	  float chWeight(iweight);

	  //update weight and mass for photons
	  LorentzVector iboson(boson);
	  if(gammaWgtHandler!=0)
	    {
	      std::vector<Float_t> photonVars;
	      photonVars.push_back(iboson.pt());
	      //photonVars.push_back(met[0].pt()/iboson.pt());
	      float photonWeight=gammaWgtHandler->getWeightFor(photonVars,tags_full);
	      chWeight *= photonWeight;
	      iboson   = gammaWgtHandler->getMassiveP4(iboson,tags_full);
	    }
	  
	  //updet the transverse mass
	  float mt =higgs::utils::transverseMass(iboson,zvv,true);

	  //scan the MET cut and fill the shapes
	  for(unsigned int index=0;index<optim_Cuts1_met.size();index++){             
	    
	    if(zvv.pt()>optim_Cuts1_met[index]){
	      for(unsigned int nri=0;nri<NRparams.size();nri++){
		
		float nrweight=chWeight*NRweights[nri];
		if(nri>0)
		  {
		    nrweight=chWeight*NRweights[nri];
		    if(lShapeWeights[0]==0) nrweight=0;
		    else                    nrweight/=lShapeWeights[0];
		  }
		
		if(passPreselection && ivar==0 && nri==0                                    )   mon.fillHisto("metcount", tags_full, index, nrweight);
		if(passPreselection                                                         )   mon.fillHisto(TString("mt_shapes")+NRsuffix[nri]+varNames[ivar],tags_full,index, mt,nrweight);
		if(passPreselection                                                         )   mon.fillHisto(TString("met_shapes")+NRsuffix[nri]+varNames[ivar],tags_full,index, zvv.pt(),nrweight);
		if(passPreselectionMbvetoMzmass && passMass          && passLocalBveto      )   mon.fillHisto("mt_shapes_NRBctrl"+NRsuffix[nri]+varNames[ivar],tags_full,index,0,nrweight);
		if(passPreselectionMbvetoMzmass && isZsideBand       && passLocalBveto      )   mon.fillHisto("mt_shapes_NRBctrl"+NRsuffix[nri]+varNames[ivar],tags_full,index,1,nrweight);
		if(passPreselectionMbvetoMzmass && isZsideBandPlus   && passLocalBveto      )   mon.fillHisto("mt_shapes_NRBctrl"+NRsuffix[nri]+varNames[ivar],tags_full,index,2,nrweight);
		if(passPreselectionMbvetoMzmass && passMass          && !passLocalBveto     )   mon.fillHisto("mt_shapes_NRBctrl"+NRsuffix[nri]+varNames[ivar],tags_full,index,3,nrweight);
		if(passPreselectionMbvetoMzmass && isZsideBand       && !passLocalBveto     )   mon.fillHisto("mt_shapes_NRBctrl"+NRsuffix[nri]+varNames[ivar],tags_full,index,4,nrweight);
		if(passPreselectionMbvetoMzmass && isZsideBandPlus   && !passLocalBveto     )   mon.fillHisto("mt_shapes_NRBctrl"+NRsuffix[nri]+varNames[ivar],tags_full,index,5,nrweight);
	      }
	    }
	  }
	}
      }
  }
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

  if(outTxtFile)fclose(outTxtFile);
}  







