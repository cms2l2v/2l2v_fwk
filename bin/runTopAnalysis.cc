#include <iostream>
#include <boost/shared_ptr.hpp>
#include <fstream>

#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/llvv_fwk/interface/DataEventSummaryHandler.h"
#include "UserCode/llvv_fwk/interface/LxyAnalysis.h"
#include "UserCode/llvv_fwk/interface/UEAnalysis.h"
#include "UserCode/llvv_fwk/interface/BTVAnalysis.h"
#include "UserCode/llvv_fwk/interface/LeptonEfficiencySF.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
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
  TString url         = runProcess.getParameter<std::string>("input");
  TString baseDir     = runProcess.getParameter<std::string>("dirName");
  bool runSystematics = runProcess.getParameter<bool>("runSystematics");
  TString jecDir      = runProcess.getParameter<std::string>("jecDir");
  bool isMC           = runProcess.getParameter<bool>("isMC");
  int mcTruthMode     = runProcess.getParameter<int>("mctruthmode");
  double xsec         = runProcess.getParameter<double>("xsec");
  bool isV0JetsMC(isMC && (url.Contains("DYJetsToLL_50toInf") || url.Contains("WJets")));
  bool hasTop          = (url.Contains("t#bar{t}") || url.Contains("TT") || url.Contains("SingleT"));
  TString out          = runProcess.getParameter<std::string>("outdir");
  bool saveSummaryTree = runProcess.getParameter<bool>("saveSummaryTree");
  std::vector<string>  weightsFile = runProcess.getParameter<std::vector<string> >("weightsFile");
   

  //jet energy scale uncertainties
  gSystem->ExpandPathName(jecDir);
  FactorizedJetCorrector *jesCor        = utils::cmssw::getJetCorrector(jecDir,isMC);
  JetCorrectionUncertainty *totalJESUnc = new JetCorrectionUncertainty((jecDir+"/MC_Uncertainty_AK5PFchs.txt").Data());

  //muon energy scale and uncertainties
  MuScleFitCorrector *muCor=utils::cmssw::getMuonCorrector(jecDir,url);

  //re-scale dy in the signal region
  std::map<TString,float> dySFmap;
  if(weightsFile.size() && url.Contains("DY") && isMC)
    {
      TFile *dyF=TFile::Open(weightsFile[0].c_str());
      TH1* dysfH=(TH1 *)dyF->Get("dysf");
      for(int ibin=1; ibin<=dysfH->GetXaxis()->GetNbins(); ibin++) { dySFmap[dysfH->GetXaxis()->GetBinLabel(ibin)]=dysfH->GetBinContent(ibin); 
	cout << dysfH->GetXaxis()->GetBinLabel(ibin) << " " << dysfH->GetBinContent(ibin) << endl;
      }
      dyF->Close();
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
  

  //
  // control histograms
  //
  SmartSelectionMonitor controlHistos;
  TH1F* Hhepup        = (TH1F* )controlHistos.addHistogram(new TH1F ("heupnup"    , "hepupnup"    ,20,0,20) ) ;
  TH1F* Hcutflow      = (TH1F*) controlHistos.addHistogram(new TH1F ("cutflow"    , "cutflow"    ,5,0,5) ) ;

  //vertex multiplicity
  controlHistos.addHistogram( new TH1F ("nvertices", "; Vertex multiplicity; Events", 50, 0.,50.) );

  //event selection histogram
  TString labels[]={"2 leptons", "M>12 #wedge |M-M_{Z}|>15", "#geq 2 jets", "E_{T}^{miss}>40,0", "op. sign"};
  int nsteps=sizeof(labels)/sizeof(TString);
  TH1F *cutflowH = (TH1F *)controlHistos.addHistogram( new TH1F("evtflow",";Cutflow;Events",nsteps,0,nsteps) );
  for(int ibin=0; ibin<nsteps; ibin++) cutflowH->GetXaxis()->SetBinLabel(ibin+1,labels[ibin]);

  TString ctrlCats[]={"","eq1jets","lowmet","eq1jetslowmet","zlowmet","zeq1jets","zeq1jetslowmet","z"};
  for(size_t k=0;k<sizeof(ctrlCats)/sizeof(TString); k++)
    {
      controlHistos.addHistogram( new TH1F(ctrlCats[k]+"emva", "; e-id MVA; Electrons", 50, 0.95,1.0) );
      controlHistos.addHistogram( new TH1F(ctrlCats[k]+"mll",";Dilepton invariant mass [GeV];Events",50,0,250) );
      controlHistos.addHistogram( new TH1F(ctrlCats[k]+"ptll",";Dilepton transverse momentum [GeV];Events",50,0,250) );
      controlHistos.addHistogram( new TH1F(ctrlCats[k]+"dilarccosine",";#theta(l,l') [rad];Events",50,0,3.2) );
      controlHistos.addHistogram( new TH1F(ctrlCats[k]+"mtsum",";M_{T}(l^{1},E_{T}^{miss})+M_{T}(l^{2},E_{T}^{miss}) [GeV];Events",100,0,1000) );
      controlHistos.addHistogram( new TH1F(ctrlCats[k]+"met",";Missing transverse energy [GeV];Events",100,0,1000) );
      TH1F *h=(TH1F *)controlHistos.addHistogram( new TH1F(ctrlCats[k]+"njets",";Jet multiplicity;Events",6,0,6) );
      for(int ibin=1; ibin<=h->GetXaxis()->GetNbins(); ibin++)
	{
	  TString label( ibin==h->GetXaxis()->GetNbins() ? "#geq" : "=");
	  label += (ibin-1);
	  label += " jets";
	  h->GetXaxis()->SetBinLabel(ibin,label);
	  
	  if(ibin==1) continue;
	  label="jet"; label+=(ibin-1);
	  controlHistos.addHistogram( new TH1F(ctrlCats[k]+"pt"+label+"pt",";Transverse momentum [GeV];Events",50,0,250) );
	  controlHistos.addHistogram( new TH1F(ctrlCats[k]+"pt"+label+"nobsmearpt",";Transverse momentum [GeV];Events",50,0,250) );
	  controlHistos.addHistogram( new TH1F(ctrlCats[k]+"pt"+label+"smearpt",";Transverse momentum [GeV];Events",50,0,250) );
	  controlHistos.addHistogram( new TH1F(ctrlCats[k]+"pt"+label+"eta",";Pseudo-rapidity;Events",50,0,2.5) );
	  TH1F *flav1H=(TH1F *)controlHistos.addHistogram( new TH1F(ctrlCats[k]+"pt"+label+"flav",";Flavor;Events",5,0,5) );
	  controlHistos.addHistogram( new TH1F(ctrlCats[k]+"btag"+label+"pt",";Transverse momentum [GeV];Events",50,0,250) );
	  controlHistos.addHistogram( new TH1F(ctrlCats[k]+"btag"+label+"eta",";Pseudo-rapidity;Events",50,0,2.5) );
	  TH1F *flav2H=(TH1F *)controlHistos.addHistogram( new TH1F(ctrlCats[k]+"btag"+label+"flav",";Flavor;Events",5,0,5) );
	  for(int ibin=1; ibin<=5; ibin++)
	    {
	      TString label("unmatched");
	      if(ibin==2) label="g";
	      if(ibin==3) label="uds";
	      if(ibin==4) label="c";
	      if(ibin==5) label="b";
	      flav1H->GetXaxis()->SetBinLabel(ibin,label);
	      flav2H->GetXaxis()->SetBinLabel(ibin,label);
	    }
	}
    }
  
  //lepton efficiencies
  LeptonEfficiencySF lepEff;

  //UEAnalysis ueAn(controlHistos);
  //BTVAnalysis btvAn(controlHistos,runSystematics);
  LxyAnalysis lxyAn(controlHistos,runSystematics);
  
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
      spyDir = spyFile->mkdir(proctag);
      TTree *outT = evSummary.getTree()->CloneTree(0);
      outT->SetTitle("Event summary");
      outT->SetDirectory(spyDir);
      outT->SetAutoSave(1000000);
      outT->Branch("weight",&evSummaryWeight,"weight/F"); 
      spyEvents->init(outT,false);
    }



  //
  // analyze (puf...)
  //
  DuplicatesChecker duplicatesChecker;
  int nDuplicates(0);
  for (int inum=0; inum < totalEntries; ++inum)
    {
      if(inum%500==0) { printf("\r [ %d/100 ]",int(100*float(inum)/float(totalEntries))); cout << flush; }
      evSummary.getEntry(inum);
      DataEventSummary &ev = evSummary.getEvent();
      if(!isMC && duplicatesChecker.isDuplicate( ev.run, ev.lumi, ev.event) ) { nDuplicates++; continue; }

      //pileup weight
      float weight(1.0),weightUp(1.0), weightDown(1.0);
      if(LumiWeights) {
	weight     = LumiWeights->weight(ev.ngenITpu);
	weightUp   = weight*PuShifters[utils::cmssw::PUUP]->Eval(ev.ngenITpu);
	weightDown = weight*PuShifters[utils::cmssw::PUDOWN]->Eval(ev.ngenITpu);
      }
      if(isV0JetsMC && ev.nup>5)                          continue;
      Hhepup->Fill(ev.nup,weight);

      Hcutflow->Fill(1,1);
      Hcutflow->Fill(2,weight);
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
	      //bool isTight    = ((idbits >> 10) & 0x1);
	      bool isLoose    = ((idbits >> 8) & 0x1);
	      Float_t gIso    = leptons[ilep].getVal("gIso04");
	      Float_t chIso   = leptons[ilep].getVal("chIso04");
	      Float_t puchIso = leptons[ilep].getVal("puchIso04");
	      Float_t nhIso   = leptons[ilep].getVal("nhIso04");
	      Float_t relIso=(TMath::Max(nhIso+gIso-0.5*puchIso,0.)+chIso)/leptons[ilep].pt();
	      if(leptons[ilep].pt()<20)                      passKin=false;
	      if(fabs(leptons[ilep].eta())>2.4)              passKin=false;
	      if(!isLoose)                                   passId=false;
	      if(relIso>0.20)                                passIso=false;
	    }

	  if(!passKin || !passId || !passIso) continue;
	  selLeptons.push_back(leptons[ilep]);
	}
      sort(selLeptons.begin(),selLeptons.end(),data::PhysicsObject_t::sortByPt);
     
      //select the leptons
      if(!eeTrigger && !emuTrigger && !mumuTrigger) continue;
      if(selLeptons.size()<2) continue;
      
      //apply data/mc correction factors
      int dilId(1);
      float llScaleFactor(1.0);
      for(size_t ilep=0; ilep<2; ilep++)
	{
	  dilId *= selLeptons[ilep].get("id");
	  int id(abs(selLeptons[ilep].get("id")));
	  llScaleFactor *= isMC ? lepEff.getLeptonEfficiency( selLeptons[ilep].pt(), selLeptons[ilep].eta(), id,  id ==11 ? "loose" : "tight" ).first : 1.0;
	}
      weight *= llScaleFactor;
      
      //set the channel
      TString chName;
      bool isOS(dilId<0);
      bool isSameFlavor(false);
      if     (abs(dilId)==11*11 && eeTrigger)   { chName="ee";  isSameFlavor=true; }
      else if(abs(dilId)==11*13 && emuTrigger)  { chName="emu"; }
      else if(abs(dilId)==13*13 && mumuTrigger) { chName="mumu"; isSameFlavor=true; }
      else                                       continue;
      std::vector<TString> ch(1,chName);
      if(isSameFlavor) ch.push_back("ll");
            
      //select the jets
      data::PhysicsObjectCollection_t jets=evSummary.getPhysicsObject(DataEventSummaryHandler::JETS);
      data::PhysicsObjectCollection_t looseJets,selJets;
      for(size_t ijet=0; ijet<jets.size(); ijet++)
	{
	  //correct jet
	  float toRawSF=jets[ijet].getVal("torawsf");
	  LorentzVector rawJet(jets[ijet]*toRawSF);
	  jesCor->setJetEta(rawJet.eta());
	  jesCor->setJetPt(rawJet.pt());
	  jesCor->setJetA(jets[ijet].getVal("area"));
	  jesCor->setRho(ev.rho);
	  jesCor->setNPV(ev.nvtx);
	  float newJECSF=jesCor->getCorrection();
	  jets[ijet].SetPxPyPzE(rawJet.px(),rawJet.py(),rawJet.pz(),rawJet.energy());
	  jets[ijet] *= newJECSF;
	  jets[ijet].setVal("torawsf",1./newJECSF);

	  //cross-clean with selected leptons 
	  double minDRlj(9999.);
	  for(size_t ilep=0; ilep<selLeptons.size(); ilep++)
	    minDRlj = TMath::Min( minDRlj, deltaR(jets[ijet],selLeptons[ilep]) );
	  if(minDRlj<0.4) continue;
	  
	  //require to pass the loose id
	  Int_t idbits=jets[ijet].get("idbits");
	  bool passPFloose( ((idbits>>0) & 0x1));
	  if(!passPFloose) continue;

	  //add scale/resolution uncertainties
	  const data::PhysicsObject_t &genJet=jets[ijet].getObject("genJet");
	  std::vector<float> smearPt=utils::cmssw::smearJER(jets[ijet].pt(),jets[ijet].eta(),genJet.pt());
	  jets[ijet].setVal("jer",     isMC ? smearPt[0] : jets[ijet].pt());
	  jets[ijet].setVal("jerup",   isMC ? smearPt[1] : jets[ijet].pt());
	  jets[ijet].setVal("jerdown", isMC ? smearPt[2] : jets[ijet].pt());
	  smearPt=utils::cmssw::smearJES(jets[ijet].pt(),jets[ijet].eta(), totalJESUnc);
	  jets[ijet].setVal("jesup",   isMC ? smearPt[0] : jets[ijet].pt());
	  jets[ijet].setVal("jesdown", isMC ? smearPt[1] : jets[ijet].pt());

	  //top candidate jets
	  looseJets.push_back(jets[ijet]);
	  if(jets[ijet].pt()<30 || fabs(jets[ijet].eta())>2.5 ) continue;
	  selJets.push_back(jets[ijet]);
	}
      sort(looseJets.begin(),looseJets.end(),data::PhysicsObject_t::sortByPt);
      sort(selJets.begin(),  selJets.end(),  data::PhysicsObject_t::sortByCSV);
      
      //the met
      data::PhysicsObjectCollection_t met=evSummary.getPhysicsObject(DataEventSummaryHandler::MET);
     
      //MC truth
      data::PhysicsObjectCollection_t gen=evSummary.getPhysicsObject(DataEventSummaryHandler::GENPARTICLES);
    
      //select the event
      if(selLeptons.size()<2) continue;
      controlHistos.fillHisto("evtflow", ch, 0, weight);
      controlHistos.fillHisto("nvertices",  ch, ev.nvtx, weight);

      LorentzVector ll=selLeptons[0]+selLeptons[1];
      float mll=ll.mass();
      float thetall=utils::cmssw::getArcCos<LorentzVector>(selLeptons[0],selLeptons[1]);
      float mtsum=utils::cmssw::getMT<LorentzVector>(selLeptons[0],met[0])+utils::cmssw::getMT<LorentzVector>(selLeptons[1],met[0]);
      bool isZcand( isSameFlavor && fabs(mll-91)<15);
      bool passDilSelection(mll>12 && !isZcand);
      bool passJetSelection(selJets.size()>=2);
      bool passMetSelection( !isSameFlavor || met[0].pt()>40);

      std::vector<TString> ctrlCategs;
      float dyWeight(1.0);
      if(isOS && passDilSelection && passJetSelection  && passMetSelection)   { ctrlCategs.push_back("");        if(dySFmap.find(chName)!=dySFmap.end()) dyWeight=dySFmap[chName]; }
      if(isOS && passDilSelection && selJets.size()==1 && passMetSelection)   { ctrlCategs.push_back("eq1jets"); if(dySFmap.find(chName+"eq1jets")!=dySFmap.end()) dyWeight=dySFmap[chName+"eq1jets"];}
      if(isOS && passDilSelection && passJetSelection  && !passMetSelection)  ctrlCategs.push_back("lowmet");
      if(isOS && passDilSelection && selJets.size()==1 && !passMetSelection)  ctrlCategs.push_back("eq1jetslowmet");
      if(isOS && isZcand          && passJetSelection  && passMetSelection)   ctrlCategs.push_back("z");
      if(isOS && isZcand          && selJets.size()==1 && passMetSelection)   ctrlCategs.push_back("zeq1jets");
      if(isOS && isZcand          && passJetSelection  && !passMetSelection)  ctrlCategs.push_back("zlowmet");
      if(isOS && isZcand          && selJets.size()==1 && !passMetSelection)  ctrlCategs.push_back("zeq1jetslowmet");

      //control distributions
      if(isOS && passDilSelection && passMetSelection) controlHistos.fillHisto("njets",        ch, selJets.size(), dyWeight*weight);
      for(size_t icat=0; icat<ctrlCategs.size(); icat++)
	{
	  float iweight(weight);
	  if(passMetSelection) iweight *= dyWeight;

	  for(size_t ilep=0; ilep<2; ilep++)
	    {
	      if(abs(selLeptons[ilep].get("id"))!=11) continue;
	      controlHistos.fillHisto(ctrlCategs[icat]+"emva", ch, selLeptons[ilep].getVal("mvatrig"), iweight);
	    }
	  controlHistos.fillHisto(ctrlCategs[icat]+"mll",          ch, mll,            iweight);
	  controlHistos.fillHisto(ctrlCategs[icat]+"ptll",         ch, ll.pt(),        iweight);
	  controlHistos.fillHisto(ctrlCategs[icat]+"mtsum",        ch, mtsum,          iweight);
	  controlHistos.fillHisto(ctrlCategs[icat]+"dilarccosine", ch, thetall,        iweight);
	  controlHistos.fillHisto(ctrlCategs[icat]+"met",          ch, met[0].pt(),    iweight);
	  if(ctrlCategs[icat]!="") controlHistos.fillHisto(ctrlCategs[icat]+"njets",        ch, selJets.size(), iweight);
	  
	  for(size_t ijet=0; ijet<looseJets.size(); ijet++)
	    {
	      if(looseJets[ijet].pt()<30 || fabs(looseJets[ijet].eta())>2.5) continue;
	      TString label("jet"); label+=(ijet+1);
	      const data::PhysicsObject_t &genJet=looseJets[ijet].getObject("genJet");
	      int flavId=genJet.info.find("id")->second;
	      if(abs(flavId)==5 || abs(flavId)==4 ) flavId=abs(flavId)-1;
	      else if(abs(flavId)>6)                flavId=1;
	      else if(abs(flavId)==0)               flavId=0;
	      else                                  flavId=2;
	      controlHistos.fillHisto(ctrlCategs[icat]+"pt"+label+"pt",        ch, looseJets[ijet].pt(), iweight);
	      controlHistos.fillHisto(ctrlCategs[icat]+"pt"+label+"eta",       ch, fabs(looseJets[ijet].eta()), iweight);
	      controlHistos.fillHisto(ctrlCategs[icat]+"pt"+label+"flav",      ch, abs(flavId), iweight);
	      controlHistos.fillHisto(ctrlCategs[icat]+"pt"+label+"nobsmearpt",ch, abs(flavId)==5 ? looseJets[ijet].pt() : looseJets[ijet].getVal("jer"), iweight);
	      controlHistos.fillHisto(ctrlCategs[icat]+"pt"+label+"smearpt",   ch,                                         looseJets[ijet].getVal("jer"), iweight);
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
	      controlHistos.fillHisto(ctrlCategs[icat]+"btag"+label+"pt",        ch, selJets[ijet].pt(), iweight);
	      controlHistos.fillHisto(ctrlCategs[icat]+"btag"+label+"eta",       ch, fabs(selJets[ijet].eta()), iweight);
	      controlHistos.fillHisto(ctrlCategs[icat]+"btag"+label+"flav",      ch, abs(flavId), iweight);
	      controlHistos.fillHisto(ctrlCategs[icat]+"btag"+label+"nobsmearpt",ch, abs(flavId)==5 ? selJets[ijet].pt() : selJets[ijet].getVal("jer"), iweight);
	      controlHistos.fillHisto(ctrlCategs[icat]+"btag"+label+"smearpt",   ch,                                       selJets[ijet].getVal("jer"), iweight);
	    }
	}

      
      //if(passDilSelection &&                     passMetSelection && isOS) btvAn.analyze(selLeptons,looseJets,isMC,ev.nvtx,weight*dyWeight,weightUp*dyWeight,weightDown*dyWeight,hasTop);
      if(passDilSelection && passJetSelection &&                     isOS) lxyAn.analyze(selLeptons,selJets,met[0],gen,weight);

      //select the event
      if(!passDilSelection) continue;
      controlHistos.fillHisto("evtflow", ch, 1, weight);
      
      if(!passJetSelection) continue;
      controlHistos.fillHisto("evtflow", ch, 2, weight);

      if(passMetSelection) {
	
	float iweight = weight*dyWeight;
	controlHistos.fillHisto("evtflow", ch, 3, iweight);

	if(isOS) {
	  controlHistos.fillHisto("evtflow", ch, 4, iweight);
	  if(spyEvents){
	    spyEvents->getEvent().cat=dilId;
	    evSummaryWeight=xsecWeight*iweight;
	    spyEvents->getTree()->Fill();
	  }
	}
      }
      
      //UE event analysis (no need to require MET, after 2-btags the events will be pure in ttbar)
      float nbtags(0);
      for(size_t ijet=0; ijet<selJets.size(); ijet++) nbtags += (selJets[ijet].getVal("supercsv")>0.531);
      if(!isOS || nbtags<2) continue;
      
      //PF candidates
      //data::PhysicsObjectCollection_t pf = evSummary.getPhysicsObject(DataEventSummaryHandler::PFCANDIDATES);
      //ueAn.analyze(selLeptons,selJets,met,pf,gen,weight);
    }
  if(nDuplicates) cout << "[Warning] found " << nDuplicates << " duplicate events in this ntuple" << endl;

  //
  // close opened files
  // 
  inF->Close();
  if(spyFile){
    spyDir->cd(); spyEvents->getTree()->Write();
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
