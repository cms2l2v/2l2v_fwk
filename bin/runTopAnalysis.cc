#include <iostream>
#include <boost/shared_ptr.hpp>
#include <fstream>

#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/llvv_fwk/interface/DataEventSummaryHandler.h"
#include "UserCode/llvv_fwk/interface/LxyAnalysis.h"
#include "UserCode/llvv_fwk/interface/UEAnalysis.h"
#include "UserCode/llvv_fwk/interface/BTVAnalysis.h"
#include "UserCode/llvv_fwk/interface/RAnalysis.h"
#include "UserCode/llvv_fwk/interface/TopPtWeighter.h"
#include "UserCode/llvv_fwk/interface/LeptonEfficiencySF.h"
#include "UserCode/llvv_fwk/interface/MuScleFitCorrector.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "PhysicsTools/CondLiteIO/interface/RecordWriter.h"
#include "DataFormats/FWLite/interface/Record.h"
#include "DataFormats/FWLite/interface/EventSetup.h"
#include "DataFormats/FWLite/interface/ESHandle.h"
#include "CondFormats/PhysicsToolsObjects/interface/BinningPointByMap.h"
#include "RecoBTag/PerformanceDB/interface/BtagPerformance.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"

#include "TSystem.h"
#include "TFile.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TCanvas.h"
#include "TString.h"
#include "TDirectory.h"
#include "TEventList.h"
#include "TRandom.h"

#include <iostream>

using namespace std;

//
int main(int argc, char* argv[])
{
  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  //check arguments
  if ( argc < 2 ) {
    std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl;
    return 0;
  }
  
  //
  // configure
  //
  const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");
  std::vector<std::string> urls=runProcess.getParameter<std::vector<std::string> >("input");
  TString url = TString(urls[0]);
  TString baseDir     = runProcess.getParameter<std::string>("dirName");
  bool runSystematics = runProcess.getParameter<bool>("runSystematics");
  TString jecDir      = runProcess.getParameter<std::string>("jecDir");
  bool isMC           = runProcess.getParameter<bool>("isMC");
  int mcTruthMode     = runProcess.getParameter<int>("mctruthmode");
  double xsec         = runProcess.getParameter<double>("xsec");
  bool isV0JetsMC(isMC && (url.Contains("DYJetsToLL_50toInf") || url.Contains("WJets")));
  bool isTTbarMC(isMC && (url.Contains("TTJets") || url.Contains("_TT_")));
  TString out          = runProcess.getParameter<std::string>("outdir");
  bool saveSummaryTree = runProcess.getParameter<bool>("saveSummaryTree");
  std::vector<string>  weightsFile = runProcess.getParameter<std::vector<string> >("weightsFile");
   
  //jet energy scale uncertainties
  gSystem->ExpandPathName(jecDir);
  FactorizedJetCorrector *jesCor        = utils::cmssw::getJetCorrector(jecDir,isMC);
  JetCorrectionUncertainty *totalJESUnc = new JetCorrectionUncertainty((jecDir+"/MC_Uncertainty_AK5PFchs.txt").Data());

  //muon energy scale and uncertainties
  MuScleFitCorrector *muCor=getMuonCorrector(jecDir,url);

  //re-scale dy in the signal region
  std::map<TString,float> dySFmap;
  if(weightsFile.size() && url.Contains("DY") && isMC)
    {
      TString dyWgtsUrl(weightsFile[0].c_str()); dyWgtsUrl += "/top_dysf.root";
      gSystem->ExpandPathName(dyWgtsUrl);
      TFile *dyF=TFile::Open(dyWgtsUrl);
      if(dyF!=0 && !dyF->IsZombie())
	{
	  TH1* dysfH=(TH1 *)dyF->Get("dysf");
	  for(int ibin=1; ibin<=dysfH->GetXaxis()->GetNbins(); ibin++) { 
	    dySFmap[dysfH->GetXaxis()->GetBinLabel(ibin)]=dysfH->GetBinContent(ibin); 
	    cout << dysfH->GetXaxis()->GetBinLabel(ibin) << " " << dysfH->GetBinContent(ibin) << endl;
	  }
	  dyF->Close();
	}
    }

  //
  // check input file
  //
  TFile *inF = TFile::Open(url);
  if(inF==0) return -1;
  if(inF->IsZombie()) return -1;
  TString proctag=gSystem->BaseName(url);
  Ssiz_t pos=proctag.Index(".root");
  proctag.Remove(pos,proctag.Length());

  bool filterOnlyEE(false), filterOnlyEMU(false), filterOnlyMUMU(false);
  if(!isMC)
    {
      if(url.Contains("DoubleEle")) filterOnlyEE=true;
      if(url.Contains("DoubleMu"))  filterOnlyMUMU=true;
      if(url.Contains("MuEG"))      filterOnlyEMU=true;
    }
  
  //
  // pileup reweighter
  //
  std::vector<double> dataPileupDistributionDouble = runProcess.getParameter< std::vector<double> >("datapileup");
  std::vector<float> dataPileupDistribution; for(unsigned int i=0;i<dataPileupDistributionDouble.size();i++){dataPileupDistribution.push_back(dataPileupDistributionDouble[i]);}
  std::vector<float> mcPileupDistribution;
  if(isMC){
    TString puDist(baseDir+"/pileup");
    TH1F* histo = (TH1F *) inF->Get(puDist);
    if(!histo)std::cout<<"pileup histogram is null!!!\n";
    for(int i=1;i<=histo->GetNbinsX();i++){mcPileupDistribution.push_back(histo->GetBinContent(i));}
    delete histo;
  }
  while(mcPileupDistribution.size()<dataPileupDistribution.size())  mcPileupDistribution.push_back(0.0);
  while(mcPileupDistribution.size()>dataPileupDistribution.size())dataPileupDistribution.push_back(0.0);
  gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE
  edm::LumiReWeighting *LumiWeights= isMC ? new edm::LumiReWeighting(mcPileupDistribution,dataPileupDistribution): 0;
  utils::cmssw::PuShifter_t PuShifters;
  if(isMC) PuShifters=utils::cmssw::getPUshifters(dataPileupDistribution,0.05);
  
  //systematic variations for final selection
  std::vector<TString> systVars(1,"");
  if(isMC && runSystematics){
    systVars.push_back("jerdown");     systVars.push_back("jerup"); 
    systVars.push_back("jesdown");     systVars.push_back("jesup"); 
    systVars.push_back("umetdown");    systVars.push_back("umetup"); 
    systVars.push_back("pudown");      systVars.push_back("puup"); 
    if(isTTbarMC){
      systVars.push_back("topptdown"); systVars.push_back("topptup");
    }
  }

  //
  // control histograms
  //
  SmartSelectionMonitor controlHistos;
  TH1F* Hhepup        = (TH1F* )controlHistos.addHistogram(new TH1F ("heupnup"    , "hepupnup"    ,20,0,20) ) ;
  TH1F* Hcutflow      = (TH1F*) controlHistos.addHistogram(new TH1F ("cutflow"    , "cutflow"    ,5,0,5) ) ;
  TH1F* Hoptim_systs  = (TH1F*) controlHistos.addHistogram(new TH1F ("optim_systs"    , ";syst;", systVars.size(),0,systVars.size()) );
  

  //vertex multiplicity
  controlHistos.addHistogram( new TH1F ("nvertices", "; Vertex multiplicity; Events", 50, 0.,50.) );

  //event selection histogram
  TString labels[]={"2 leptons", "M>12 #wedge |M-M_{Z}|>15", "#geq 2 jets", "E_{T}^{miss}>40,0", "op. sign"};
  int nsteps=sizeof(labels)/sizeof(TString);
  TH1F *cutflowH = (TH1F *)controlHistos.addHistogram( new TH1F("evtflow",";Cutflow;Events",nsteps,0,nsteps) );
  for(int ibin=0; ibin<nsteps; ibin++) cutflowH->GetXaxis()->SetBinLabel(ibin+1,labels[ibin]);

  TString uelabels[]={"#geq 2 leptons", "dilepton", "#geq 2 jets", "#geq 2-btags"};
  int nuesteps=sizeof(uelabels)/sizeof(TString);
  TH1F *uecutflowH = (TH1F *)controlHistos.addHistogram( new TH1F("ueevtflow",";Cutflow;Events",nuesteps,0,nuesteps) );
  for(int ibin=0; ibin<nuesteps; ibin++) uecutflowH->GetXaxis()->SetBinLabel(ibin+1,uelabels[ibin]);
  
  for(size_t ivar=0;ivar<systVars.size(); ivar++) 
    {
      Hoptim_systs->GetXaxis()->SetBinLabel(ivar+1,systVars[ivar]);
      TH1D *finalCutflowH=new TH1D("finalevtflow"+systVars[ivar],";Category;Events",4,0,4);
      finalCutflowH->GetXaxis()->SetBinLabel(1,"=1 jets");
      finalCutflowH->GetXaxis()->SetBinLabel(2,"=2 jets");
      finalCutflowH->GetXaxis()->SetBinLabel(3,"=3 jets");
      finalCutflowH->GetXaxis()->SetBinLabel(4,"=4 jets");
      controlHistos.addHistogram( finalCutflowH );
    }


  TString ctrlCats[]={"","eq1jets","lowmet","eq1jetslowmet","osbtag","osbveto"};
  for(size_t k=0;k<sizeof(ctrlCats)/sizeof(TString); k++)
    {
      controlHistos.addHistogram( new TH1F("eeta", "; Pseudo-rapidity; Electrons", 50,0,2.5) );
      controlHistos.addHistogram( new TH1F("mueta", "; Pseudo-rapidity; Muons", 50,0,2.5) );

      controlHistos.addHistogram( new TH1F(ctrlCats[k]+"emva", "; e-id MVA; Electrons", 50, 0.95,1.0) );
      controlHistos.addHistogram( new TH1F(ctrlCats[k]+"mll",";Dilepton invariant mass [GeV];Events",50,0,250) );
      controlHistos.addHistogram( new TH1F(ctrlCats[k]+"ptll",";Dilepton transverse momentum [GeV];Events",50,0,250) );
      controlHistos.addHistogram( new TH1F(ctrlCats[k]+"met",";Missing transverse energy [GeV];Events",50,0,500) );
      controlHistos.addHistogram( new TH1F(ctrlCats[k]+"metnotoppt",";Missing transverse energy [GeV];Events",50,0,500) );
      TH1F *h=(TH1F *)controlHistos.addHistogram( new TH1F(ctrlCats[k]+"njets",";Jet multiplicity;Events",6,0,6) );
      TH1F *h2=(TH1F *)controlHistos.addHistogram( new TH1F(ctrlCats[k]+"njetsnotoppt",";Jet multiplicity;Events",6,0,6) );
      for(int ibin=1; ibin<=h->GetXaxis()->GetNbins(); ibin++)
	{
	  TString label( ibin==h->GetXaxis()->GetNbins() ? "#geq" : "=");
	  label += (ibin-1);
	  label += " jets";
	  h->GetXaxis()->SetBinLabel(ibin,label);
	  h2->GetXaxis()->SetBinLabel(ibin,label);
	  
	  if(ibin==1) continue;
	  label="jet"; label+=(ibin-1);
	  Float_t jetPtaxis[]={30,35,40,45,50,55,60,65,70,80,90,100,125,150,200,250,500};
	  const size_t nJetPtBins=sizeof(jetPtaxis)/sizeof(Float_t)-1;
	  controlHistos.addHistogram( new TH1F(ctrlCats[k]+"btag"+label+"pt",";Transverse momentum [GeV];Events",nJetPtBins,jetPtaxis) );
	  controlHistos.addHistogram( new TH1F(ctrlCats[k]+"btag"+label+"ptnotoppt",";Transverse momentum [GeV];Events",nJetPtBins,jetPtaxis) );
	  controlHistos.addHistogram( new TH1F(ctrlCats[k]+"btag"+label+"eta",";Pseudo-rapidity;Events",25,0,2.5) );
	  TH1F *flavH=(TH1F *)controlHistos.addHistogram( new TH1F(ctrlCats[k]+"btag"+label+"flav",";Flavor;Events",5,0,5) );
	  for(int ibin=1; ibin<=5; ibin++)
	    {
	      TString label("unmatched");
	      if(ibin==2) label="g";
	      if(ibin==3) label="uds";
	      if(ibin==4) label="c";
	      if(ibin==5) label="b";
	      flavH->GetXaxis()->SetBinLabel(ibin,label);
	    }
	}

      //for DY estimation
      for(size_t ivar=0;ivar<systVars.size(); ivar++) 
	{
	  controlHistos.addHistogram( new TH1F(ctrlCats[k]+"dilarccosine"+systVars[ivar],";#theta(l,l') [rad];Events",50,0,3.2) );
	  controlHistos.addHistogram( new TH1F(ctrlCats[k]+"mtsum"+systVars[ivar],";M_{T}(l^{1},E_{T}^{miss})+M_{T}(l^{2},E_{T}^{miss}) [GeV];Events",100,0,1000) );
	}
    }

  controlHistos.addHistogram( new TH1F("leadpt",";Transverse momentum [GeV];Events",    50,0,500) );
  controlHistos.addHistogram( new TH1F("trailerpt",";Transverse momentum [GeV];Events", 50,0,500) );
  controlHistos.addHistogram( new TH1F("dilpt",";|#Sigma #vec{p}_{T}| [GeV];Events",    50,0,500) );
  controlHistos.addHistogram( new TH1F("sumpt",";#Sigma p_{T} [GeV];Events",            50,0,500) );
  TH1F *hch=(TH1F *)controlHistos.addHistogram( new TH1F("dilcharge",";Charge;Events",2,0,2) );
  hch->GetXaxis()->SetBinLabel(1,"os");
  hch->GetXaxis()->SetBinLabel(2,"ss");
  
  //lepton efficiencies
  LeptonEfficiencySF lepEff;

  ///
  // process events file
  //
  DataEventSummaryHandler evSummary;
  if( !evSummary.attach( (TTree *) inF->Get(baseDir+"/data") ) )  { inF->Close();  return -1; }  
  const Int_t totalEntries=evSummary.getEntries();
  
  float cnorm=1.0;
  if(isMC){
    TH1F* cutflowH = (TH1F *) inF->Get(baseDir+"/cutflow");
    if(cutflowH) cnorm=cutflowH->GetBinContent(1);
  }
  Hcutflow->SetBinContent(1,cnorm);

  cout << "Processing: " << proctag << " @ " << url << endl
       << "Initial number of events: " << cnorm << endl
       << "Events in tree:           " << totalEntries << endl
       << " xSec x BR:               " << xsec << endl;

  //TOP analysis
  //UEAnalysis ueAn(controlHistos);
  //TNtuple *ueNtuple = ueAn.getSummaryTuple();
  //BTVAnalysis btvAn(controlHistos,runSystematics);
  //LxyAnalysis lxyAn(controlHistos,runSystematics);
  RAnalysis rAn(controlHistos,systVars);

  TopPtWeighter *topPtWgt=0;
  if(isTTbarMC ){
    TString shapesDir("");
    if(weightsFile.size()) shapesDir=weightsFile[0].c_str();
    topPtWgt = new TopPtWeighter( proctag, out, shapesDir, evSummary.getTree() );
  }

  //check if a summary should be saved
  Float_t evSummaryWeight(1.0);
  Float_t xsecWeight(isMC ? xsec/cnorm : 1.0);
  TFile *spyFile=0;
  TDirectory *spyDir=0;
  DataEventSummaryHandler *spyEvents=0;
  if(saveSummaryTree)
    {
      gSystem->Exec("mkdir -p " + out);
      gDirectory->SaveSelf();
      TString summaryName(out + "/" + proctag);
      if(mcTruthMode!=0) { summaryName += "_filt"; summaryName += mcTruthMode; } 
      summaryName += "_summary.root";
      gSystem->ExpandPathName(summaryName);
      cout << "Creating event summary file @ " << summaryName << endl;

      //open file
      spyEvents = new DataEventSummaryHandler;
      spyFile = TFile::Open(summaryName,"RECREATE");
      spyFile->rmdir(proctag);
      spyDir = spyFile->mkdir("dataAnalyzer");
      TTree *outT = evSummary.getTree()->CloneTree(0);
      outT->SetTitle("Event summary");
      outT->SetDirectory(spyDir);
      outT->SetAutoSave(1000000);
      outT->Branch("weight",&evSummaryWeight,"weight/F"); 
      spyEvents->init(outT,false);

      //add also other summary tuples
      //ueNtuple->SetDirectory(spyDir);
    }



  //
  // analyze (puf...)
  //
  DuplicatesChecker duplicatesChecker;
  int nDuplicates(0);
  int nChargeFlaws(0);
  for (int inum=0; inum < totalEntries; ++inum)
    {
      if(inum%500==0) { printf("\r [ %d/100 ]",int(100*float(inum)/float(totalEntries))); cout << flush; }
      evSummary.getEntry(inum);
      DataEventSummary &ev = evSummary.getEvent();
      if(!isMC && duplicatesChecker.isDuplicate( ev.run, ev.lumi, ev.event) ) { nDuplicates++; continue; }

      //pileup weight
      float weightNom(1.0),weightUp(1.0), weightDown(1.0);
      if(LumiWeights) {
	weightNom  = LumiWeights->weight(ev.ngenITpu);
	weightUp   = weightNom*PuShifters[utils::cmssw::PUUP]->Eval(ev.ngenITpu);
	weightDown = weightNom*PuShifters[utils::cmssw::PUDOWN]->Eval(ev.ngenITpu);
      }

      if(isV0JetsMC && ev.nup>5)                          continue;
      Hhepup->Fill(ev.nup,1);

      //MC truth
      data::PhysicsObjectCollection_t gen=evSummary.getPhysicsObject(DataEventSummaryHandler::GENPARTICLES);
      bool hasTop(false);
      int ngenLeptonsStatus3(0);
      float wgtTopPt(1.0), wgtTopPtUp(1.0), wgtTopPtDown(1.0);
      if(isMC)
	{
	  float pttop(0), ptantitop(0);
	  for(size_t igen=0; igen<gen.size(); igen++){
	    if(gen[igen].get("status")!=3) continue;
	    int absid=abs(gen[igen].get("id"));
	    if(absid==6) {
	      hasTop=true;
	      if(gen[igen].get("id")==6) pttop=gen[igen].pt();
	      else                       ptantitop=gen[igen].pt();
	    }
	    if(absid!=11 && absid!=13 && absid!=15) continue;
	    ngenLeptonsStatus3++;
	  }
	  if(mcTruthMode==1 && (ngenLeptonsStatus3!=2 || !hasTop)) continue;
	  if(mcTruthMode==2 && (ngenLeptonsStatus3==2 || !hasTop)) continue;
	  if(pttop>0 && ptantitop>0 && topPtWgt)
	    {
	      topPtWgt->computeWeight(pttop,ptantitop);
	      topPtWgt->getEventWeight(wgtTopPt, wgtTopPtUp, wgtTopPtDown );
	      wgtTopPtUp /= wgtTopPt;
	      wgtTopPtDown /= wgtTopPt;
	      //cout << wgtTopPt << " " << wgtTopPtUp << " " << wgtTopPtDown << endl;
	    }
	}


      Hcutflow->Fill(1,1);
      Hcutflow->Fill(2,weightNom);
      Hcutflow->Fill(3,weightUp);
      Hcutflow->Fill(4,weightDown);

      //trigger bits
      bool eeTrigger   = ev.t_bits[0];
      bool emuTrigger  = ev.t_bits[4] || ev.t_bits[5];
      bool mumuTrigger = ev.t_bits[2] || ev.t_bits[3];
      if(filterOnlyEE)   {                   emuTrigger=false;  mumuTrigger=false; }
      if(filterOnlyEMU)  { eeTrigger=false;                     mumuTrigger=false; }
      if(filterOnlyMUMU) { eeTrigger=false;  emuTrigger=false;                     }
 
      data::PhysicsObjectCollection_t leptons=evSummary.getPhysicsObject(DataEventSummaryHandler::LEPTONS);
      data::PhysicsObjectCollection_t selLeptons;
      for(size_t ilep=0; ilep<leptons.size(); ilep++)
	{
	  Int_t id=leptons[ilep].get("id");
	  bool passKin(true),passId(true),passIso(true);
	  if(abs(id)==11)
	    {
	      float sceta=leptons[ilep].getVal("sceta");
	      Float_t gIso    = leptons[ilep].getVal("gIso03");
	      Float_t chIso   = leptons[ilep].getVal("chIso03");
	      //Float_t puchIso = leptons[ilep].getVal("puchIso03");
	      Float_t nhIso   = leptons[ilep].getVal("nhIso03");
	      float relIso=(TMath::Max(nhIso+gIso-ev.rho*utils::cmssw::getEffectiveArea(11,sceta),Float_t(0.))+chIso)/leptons[ilep].pt();
	      if(leptons[ilep].pt()<20)                      passKin=false;
	      if(fabs(leptons[ilep].eta())>2.5)              passKin=false;
	      if(fabs(sceta)>1.4442 && fabs(sceta)<1.5660)   passKin=false;
	      if(leptons[ilep].getFlag("isconv"))            passId=false;
	      if(leptons[ilep].getVal("tk_d0")>0.4)          passId=false;
	      if(leptons[ilep].getVal("tk_lostInnerHits")>0) passId=false;
	      if(leptons[ilep].getVal("mvatrig")<0.5)        passId=false;
	      if(relIso>0.15)                                passIso=false;
	    }
	  else
	    {
	      if(muCor){
		TLorentzVector p4(leptons[ilep].px(),leptons[ilep].py(),leptons[ilep].pz(),leptons[ilep].energy());
		muCor->applyPtCorrection(p4 , id<0 ? -1 : 1 );
		if(isMC) muCor->applyPtSmearing(p4, id<0 ? -1 : 1, false);
		leptons[ilep].SetPxPyPzE(p4.Px(),p4.Py(),p4.Pz(),p4.E());
	      }
	      
	      Int_t idbits    = leptons[ilep].get("idbits");
	      bool isTight    = ((idbits >> 10) & 0x1);
	      Float_t gIso    = leptons[ilep].getVal("gIso04");
	      Float_t chIso   = leptons[ilep].getVal("chIso04");
	      Float_t puchIso = leptons[ilep].getVal("puchIso04");
	      Float_t nhIso   = leptons[ilep].getVal("nhIso04");
	      Float_t relIso=(TMath::Max(nhIso+gIso-0.5*puchIso,0.)+chIso)/leptons[ilep].pt();
	      if(leptons[ilep].pt()<20)                      passKin=false;
	      if(fabs(leptons[ilep].eta())>2.4)              passKin=false;
	      if(!isTight)                                   passId=false;
	      if(relIso>0.12)                                passIso=false;
	    }

	  if(!passKin || !passId || !passIso) continue;
	  selLeptons.push_back(leptons[ilep]);

	  const data::PhysicsObject_t &genLep=leptons[ilep].getObject("gen");
	  int genId=genLep.info.find("id")->second;
	  if(genId!=0 && id!=genId)  nChargeFlaws++;
	  
	}
      sort(selLeptons.begin(),selLeptons.end(),data::PhysicsObject_t::sortByPt);
     
      //select the leptons
      if(!eeTrigger && !emuTrigger && !mumuTrigger) continue;
      if(selLeptons.size()<2) continue;
      
      //determine the dilepton channel
      ev.cat=1;
      float llScaleFactor(1.0);
      for(size_t ilep=0; ilep<2; ilep++)
	{
	  ev.cat *= selLeptons[ilep].get("id");
	  int id(abs(selLeptons[ilep].get("id")));
	  llScaleFactor *= isMC ? lepEff.getLeptonEfficiency( selLeptons[ilep].pt(), selLeptons[ilep].eta(), id,  id ==11 ? "loose" : "tight" ).first : 1.0;
	}

      TString chName;
      bool isOS(ev.cat<0);
      bool isSameFlavor(false);
      if     (abs(ev.cat)==11*11 && eeTrigger)   { chName="ee";  isSameFlavor=true;  if(ngenLeptonsStatus3>=2) llScaleFactor*=0.972; }
      else if(abs(ev.cat)==11*13 && emuTrigger)  { chName="emu";                     if(ngenLeptonsStatus3>=2) llScaleFactor*=0.968; }
      else if(abs(ev.cat)==13*13 && mumuTrigger) { chName="mumu"; isSameFlavor=true; if(ngenLeptonsStatus3>=2) llScaleFactor*=0.955; }
      else                                       continue;
      std::vector<TString> ch(1,chName);
      if(isSameFlavor) ch.push_back("ll");
            
      //apply data/mc correction factors and update the event weight
      float weight = weightNom*llScaleFactor*wgtTopPt;

      //jet/met
      data::PhysicsObjectCollection_t recoMet=evSummary.getPhysicsObject(DataEventSummaryHandler::MET);     
      data::PhysicsObjectCollection_t jets=evSummary.getPhysicsObject(DataEventSummaryHandler::JETS);
      utils::cmssw::updateJEC(jets,jesCor,totalJESUnc,ev.rho,ev.nvtx,isMC);
      std::vector<LorentzVector> met=utils::cmssw::getMETvariations(recoMet[0],jets,selLeptons,isMC);

      //select the jets
      data::PhysicsObjectCollection_t looseJets,selJets;
      int nbtags(0);
      for(size_t ijet=0; ijet<jets.size(); ijet++)
	{
	  //cross-clean with selected leptons 
	  double minDRlj(9999.);
	  for(size_t ilep=0; ilep<selLeptons.size(); ilep++)
	    minDRlj = TMath::Min( minDRlj, deltaR(jets[ijet],selLeptons[ilep]) );
	  if(minDRlj<0.4) continue;
	  
	  //require to pass the loose id
	  Int_t idbits=jets[ijet].get("idbits");
	  bool passPFloose( ((idbits>>0) & 0x1));
	  if(!passPFloose) continue;

	  //top candidate jets
	  looseJets.push_back(jets[ijet]);
	  if(jets[ijet].pt()<30 || fabs(jets[ijet].eta())>2.5 ) continue;
	  selJets.push_back(jets[ijet]);
	  nbtags += (jets[ijet].getVal("csv")>0.405);
	}
      sort(looseJets.begin(),looseJets.end(),data::PhysicsObject_t::sortByCSV);
      sort(selJets.begin(),  selJets.end(),  data::PhysicsObject_t::sortByCSV);
      

      //select the event
      if(selLeptons.size()<2) continue;
      controlHistos.fillHisto("evtflow", ch, 0, weight);
      controlHistos.fillHisto("ueevtflow", ch, 0, weight);
      controlHistos.fillHisto("nvertices",  ch, ev.nvtx, weight);

      LorentzVector ll=selLeptons[0]+selLeptons[1];
      float mll=ll.mass();
      bool isZcand( isSameFlavor && fabs(mll-91)<15);
      bool passDilSelection(mll>12 && !isZcand);
      bool passJetSelection(selJets.size()>=2);
      bool passMetSelection( !isSameFlavor || met[0].pt()>40);


      //
      // NOMINAL SELECTION CONTROL
      //
      // define control category and define DY weight: to be only applied to events passing the MET cut (it's an efficiency correction)
      // for events with 0 jets with the DY sf from the 1 jet bin
      std::vector<TString> ctrlCategs;
      float dyWeight(1.0),ibtagdyWeight(1.0);
      if     (isOS && passDilSelection && passJetSelection  && passMetSelection)   { ctrlCategs.push_back("");                               if(dySFmap.find(chName)!=dySFmap.end()) dyWeight=dySFmap[chName]; }
      else if(isOS && passDilSelection                      && passMetSelection)   { if(selJets.size()==1) ctrlCategs.push_back("eq1jets");  if(dySFmap.find(chName+"eq1jets")!=dySFmap.end()) dyWeight=dySFmap[chName+"eq1jets"]; }
      else if(isOS && passDilSelection && passJetSelection  && !passMetSelection)  { ctrlCategs.push_back("lowmet"); }
      else if(isOS && passDilSelection && selJets.size()==1 && !passMetSelection)  { ctrlCategs.push_back("eq1jetslowmet"); }
      else if(isOS && passDilSelection && passJetSelection  && nbtags>=2)          { ctrlCategs.push_back("osbtag");                         if(dySFmap.find(chName+"osbtag")!=dySFmap.end())  ibtagdyWeight=dySFmap[chName+"osbtag"]; }
      else if(isOS && passDilSelection && passJetSelection  && nbtags==0)          { ctrlCategs.push_back("osbveto"); }

      //control distributions
      if(isOS && passDilSelection && passMetSelection) {
	controlHistos.fillHisto("njets",        ch, selJets.size(), weight*dyWeight);
	controlHistos.fillHisto("njetsnotoppt",  ch, selJets.size(), weight*dyWeight/wgtTopPt);
      }
      for(size_t icat=0; icat<ctrlCategs.size(); icat++)
	{
	  float iweight(weight);
	  if(passMetSelection) iweight *=dyWeight;
	  for(size_t ilep=0; ilep<2; ilep++)
	    {
	      if(abs(selLeptons[ilep].get("id"))!=11) {
		controlHistos.fillHisto(ctrlCategs[icat]+"mueta", ch, fabs(selLeptons[ilep].eta()), iweight);
	      }
	      else{
		controlHistos.fillHisto(ctrlCategs[icat]+"emva", ch, selLeptons[ilep].getVal("mvatrig"), iweight);
		controlHistos.fillHisto(ctrlCategs[icat]+"eeta", ch, fabs(selLeptons[ilep].eta()), iweight);
	      }
	    }

	  controlHistos.fillHisto(ctrlCategs[icat]+"mll",          ch, mll,            iweight);
	  controlHistos.fillHisto(ctrlCategs[icat]+"ptll",         ch, ll.pt(),        iweight);
	  controlHistos.fillHisto(ctrlCategs[icat]+"met",          ch, met[0].pt(),    iweight);
	  controlHistos.fillHisto(ctrlCategs[icat]+"metnotoppt",   ch, met[0].pt(),    iweight/wgtTopPt);
	  if(ctrlCategs[icat]!="") {
	    controlHistos.fillHisto(ctrlCategs[icat]+"njets",  ch, selJets.size(), iweight);
	    controlHistos.fillHisto(ctrlCategs[icat]+"njetsnotoppt",  ch, selJets.size(), iweight/wgtTopPt);
	  }

	  for(size_t ijet=0; ijet<selJets.size(); ijet++)
	    {
	      TString label("jet"); label+=(ijet+1);
	      const data::PhysicsObject_t &genJet=selJets[ijet].getObject("genJet");
	      int flavId=genJet.info.find("id")->second;
	      if(abs(flavId)==5 || abs(flavId)==4 ) flavId=abs(flavId)-1;
	      else if(abs(flavId)>6)                flavId=1;
	      else if(abs(flavId)==0)               flavId=0;
	      else                                  flavId=2;
	      controlHistos.fillHisto(ctrlCategs[icat]+"btag"+label+"pt",        ch, selJets[ijet].pt(), iweight, true);
	      controlHistos.fillHisto(ctrlCategs[icat]+"btag"+label+"ptnotoppt", ch, selJets[ijet].pt(), iweight/wgtTopPt, true);
	      controlHistos.fillHisto(ctrlCategs[icat]+"btag"+label+"eta",       ch, fabs(selJets[ijet].eta()), iweight);
	      controlHistos.fillHisto(ctrlCategs[icat]+"btag"+label+"flav",      ch, abs(flavId), iweight);
	    }
	}

      
      //if(passDilSelection &&                     passMetSelection && isOS) btvAn.analyze(selLeptons,looseJets,isMC,ev.nvtx,weightNom*llScaleFactor*dyWeight,weightUp*llScaleFactor*dyWeight,weightDown*llScaleFactor*dyWeight,hasTop);
      //if(passDilSelection && passJetSelection &&                     isOS) lxyAn.analyze(selLeptons,selJets,met[0],gen,weightNom*llScaleFactor*dyWeight);

      //select the event
      if(!passDilSelection) continue;
      controlHistos.fillHisto("evtflow", ch, 1, weight);
      if(isOS) controlHistos.fillHisto("ueevtflow", ch, 1, weight);
      
      if(passJetSelection) {

	controlHistos.fillHisto("evtflow", ch, 2, weight);

	//UE event analysis with PF candidates
	/*
	if(isOS){
	  controlHistos.fillHisto("ueevtflow", ch, 2, weight);
	  if(looseJets[0].pt()>30 && looseJets[1].pt()>30 && fabs(looseJets[0].eta())<2.5 && fabs(looseJets[1].eta())<2.5)
	    {
	      if(looseJets[0].getVal("csv")>0.405 && looseJets[1].getVal("csv")>0.405)
		{
		  float iweight( weight*ibtagdyWeight );
		  controlHistos.fillHisto("ueevtflow", ch, 3, iweight);
		  data::PhysicsObjectCollection_t pf = evSummary.getPhysicsObject(DataEventSummaryHandler::PFCANDIDATES);
		  ueAn.analyze(selLeptons,looseJets,met[0],pf,gen,ev.nvtx,iweight);
		  if(saveSummaryTree) ueAn.fillSummaryTuple(xsecWeight);
		}
	    }
	}
	*/
	
	//other analysis
	if(passMetSelection) {
	  
	  float iweight( weight*dyWeight );
	  controlHistos.fillHisto("evtflow", ch, 3, iweight);

	  float dilcharge( leptons[0].get("id")*leptons[1].get("id") ); 
	  controlHistos.fillHisto("dilcharge",     ch, (dilcharge<0 ? 0 : 1), iweight);

	  if(isOS) {
	    controlHistos.fillHisto("evtflow", ch, 4, iweight);

	    controlHistos.fillHisto("leadpt",    ch, selLeptons[0].pt(), iweight);
	    controlHistos.fillHisto("trailerpt", ch, selLeptons[1].pt(), iweight);
	    controlHistos.fillHisto("sumpt",     ch, selLeptons[0].pt()+selLeptons[1].pt(), iweight);
	    controlHistos.fillHisto("dilpt",     ch, ll.pt(), iweight);
    	    
	    //save selected event
	    if(spyEvents){
	      evSummaryWeight=xsecWeight*iweight;
	      spyEvents->getTree()->Fill();
	    }
	  }
	}
      }
    
      //
      // STATISTICAL ANALYSIS (with systs variations)
      //
      if(isOS)
	{
	  for(size_t ivar=0;ivar<systVars.size(); ivar++) 
	    {
	      TString var=systVars[ivar];

	      //re-select the jets
	      int nlocaljets(0),nlocalbtags(0);
	      data::PhysicsObjectCollection_t localSelJets;
	      for(size_t ijet=0; ijet<looseJets.size(); ijet++)
		{
		  float pt(looseJets[ijet].pt());
		  float eta(fabs(looseJets[ijet].eta()));
		  if(var=="jerup")     pt=looseJets[ijet].getVal("jerup");
		  if(var=="jerdown")   pt=looseJets[ijet].getVal("jerdown");
		  if(var=="jesup")     pt=looseJets[ijet].getVal("jesup");
		  if(var=="jesdown")   pt=looseJets[ijet].getVal("jesdown");
		  if(eta<2.5 && pt>30) {
		    nlocaljets++;
		    nlocalbtags += (looseJets[ijet].getVal("csv")>0.405);
		    localSelJets.push_back(looseJets[ijet]);
		  }
		}
	      if(nlocaljets<1) continue;
	      bool passLocalJetSelection(nlocaljets>1);
	
	      //re-select the MET
	      int metIdx(0);
	      if(var=="jerup")    metIdx=utils::cmssw::JERUP;   if(var=="jerdown")  metIdx=utils::cmssw::JERDOWN;
	      if(var=="jesup")    metIdx=utils::cmssw::JESUP;   if(var=="jesdown")  metIdx=utils::cmssw::JESDOWN;
	      if(var=="umetup")   metIdx=utils::cmssw::UMETUP;  if(var=="umetdown") metIdx=utils::cmssw::UMETDOWN;
	      LorentzVector iMet=met[metIdx];
	      bool passLocalMetSelection( !isSameFlavor || iMet.pt()>40 );

	      //event category and dy scale factor
	      TString localCtrlCateg(""), btagCtrlCateg("");
	      float idyWeight(1.0),ibtagdyWeight(1.0);
	      if     ( passLocalMetSelection &&  passLocalJetSelection)  { localCtrlCateg="";       if(dySFmap.find(chName)!=dySFmap.end()) idyWeight=dySFmap[chName];  }
	      else if(!passLocalMetSelection &&  passLocalJetSelection)  { localCtrlCateg="lowmet"; }
	      else if(!passLocalMetSelection && !passLocalJetSelection)  { localCtrlCateg="eq1jetslowmet"; }
	      else if( passLocalMetSelection && !passLocalJetSelection)  { localCtrlCateg="eq1jets"; if(dySFmap.find(chName+"eq1jets")!=dySFmap.end()) idyWeight=dySFmap[chName+"eq1jets"];    }
	      else if( passLocalJetSelection && nlocalbtags>=2)          { btagCtrlCateg="osbtag";   if(dySFmap.find(chName+"osbtag")!=dySFmap.end())  ibtagdyWeight=dySFmap[chName+"osbtag"]; }
	      else if( passLocalJetSelection && nlocalbtags==0)          { btagCtrlCateg="osbveto"; }

	      //re-assign event weight
	      float iweight(weightNom);
	      if(var=="puup")      iweight  = weightUp;
	      if(var=="pudown")    iweight  = weightDown;
	      if(var=="topptdown") iweight *= wgtTopPtUp;
	      if(var=="topptup")   iweight *= wgtTopPtDown;

	      iweight *= llScaleFactor;
	      if(passLocalMetSelection) iweight*=idyWeight;

	      float thetall=utils::cmssw::getArcCos<LorentzVector>(selLeptons[0],selLeptons[1]);
	      float mtsum=utils::cmssw::getMT<LorentzVector>(selLeptons[0],iMet)+utils::cmssw::getMT<LorentzVector>(selLeptons[1],iMet);
	      controlHistos.fillHisto(localCtrlCateg+"mtsum"+systVars[ivar],        ch, mtsum,          iweight);
	      controlHistos.fillHisto(localCtrlCateg+"dilarccosine"+systVars[ivar], ch, thetall,        iweight);
	      if(btagCtrlCateg!="")
		{
		  float ibtagweight(iweight*ibtagdyWeight/idyWeight);
		  controlHistos.fillHisto(btagCtrlCateg+"mtsum"+systVars[ivar],        ch, mtsum,          ibtagweight);
		  controlHistos.fillHisto(btagCtrlCateg+"dilarccosine"+systVars[ivar], ch, thetall,        ibtagweight);
		}

	      //final selection
	      if(!passLocalMetSelection) continue;
	      if(nlocaljets<1 || nlocaljets>4) continue;
	      controlHistos.fillHisto("finalevtflow"+systVars[ivar], ch, nlocaljets-1, iweight);

	      //R measurement
	      if(ivar==0) rAn.prepareAnalysis(selLeptons,localSelJets);
	      rAn.analyze(selLeptons,localSelJets,iweight,systVars[ivar],hasTop);

	    }
	}
    }
  if(nDuplicates) cout << "[Warning] found " << nDuplicates << " duplicate events in this ntuple" << endl;
  if(nChargeFlaws) cout << "[Warning] found " << nChargeFlaws << " charge flaws between reco and gen" << endl;

  //
  // close opened files
  // 
  inF->Close();
  if(spyFile){
    spyDir->cd(); 
    spyEvents->getTree()->Write();
    //ueNtuple->Write();
    spyFile->Close();
  }
    
  //
  // finally, save histos to local file
  //
  TString outUrl(out);
  gSystem->ExpandPathName(outUrl);
  gSystem->Exec("mkdir -p " + outUrl);
  outUrl += "/";
  outUrl += proctag;
  if(mcTruthMode!=0) { outUrl += "_filt"; outUrl += mcTruthMode; }
  outUrl += ".root";
  TFile *file=TFile::Open(outUrl, "recreate");
  controlHistos.Write();
  file->Close();
  
  //that's all folks!
}  
