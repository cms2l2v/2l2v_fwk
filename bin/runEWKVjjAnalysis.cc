#include <iostream>
#include <boost/shared_ptr.hpp>

#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/llvv_fwk/interface/DataEventSummaryHandler.h"
#include "UserCode/llvv_fwk/interface/TMVAUtils.h"
#include "UserCode/llvv_fwk/interface/LeptonEfficiencySF.h"
#include "UserCode/llvv_fwk/interface/PDFInfo.h"
#include "UserCode/llvv_fwk/interface/MuScleFitCorrector.h"
#include "UserCode/llvv_fwk/interface/GammaWeightsHandler.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

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

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >::BetaVector BetaVector;

struct AnalysisBox_t{
  bool accept, acceptMass, acceptKin, acceptRap;
  LorentzVector vCand;
  data::PhysicsObject_t lead, trailer;
  TString cat;
};

//
AnalysisBox_t assignBox(data::PhysicsObjectCollection_t &leptons, data::PhysicsObjectCollection_t &photons, LorentzVector &met)
{
  std::vector<size_t> isoLeptons,nonIsoLeptons;
  for(size_t i=0; i<leptons.size(); i++){
    if(!leptons[i].getFlag("isIso")) nonIsoLeptons.push_back(i);
    else                             isoLeptons.push_back(i);
  }
  
  std::vector<size_t> isoPhotons,nonIsoPhotons;
  for(size_t i=0; i<photons.size(); i++){
    if(!photons[i].getFlag("isIso")) nonIsoPhotons.push_back(i);
    else                             isoPhotons.push_back(i);
  }

  //
  // ASSIGN THE BOX
  // 1. ISOLATED OBJECTS     : ll,    l,        photon
  // 2. NON-ISOLATED OBJECTS : ll,    l,        photon
  // THRESHOLDS                20/20, 30 or 26, 50
  //
  AnalysisBox_t box;
  if(isoLeptons.size()>=2){
    box.lead      = leptons[ isoLeptons[0] ] ;
    box.trailer   = leptons[ isoLeptons[1] ] ;
    if( abs(box.lead.get("id")*box.trailer.get("id"))==11*11 )      box.cat="ee";
    else if( abs(box.lead.get("id")*box.trailer.get("id"))==13*13 ) box.cat="mumu";
    else                                                    box.cat="";
    box.vCand      = box.lead+box.trailer;
    box.accept     = (box.lead.pt()>20 && box.trailer.pt()>20 && box.cat!="");
    box.acceptMass = (fabs(box.vCand.mass()-91)<15);
  }
  else if(isoLeptons.size()==1){
    box.lead      = leptons[ isoLeptons[0] ] ;
    box.trailer   = data::PhysicsObject_t(met.px(),met.py(),0,met.pt());
    box.vCand     = box.lead+box.trailer;
    box.accept    = true;
    if( abs(box.lead.get("id"))==11 )      { box.cat="e";  box.accept=(box.lead.pt()>30); }
    else if( abs(box.lead.get("id"))==13 ) { box.cat="mu"; box.accept=(box.lead.pt()>26); }
    box.accept    &= (box.trailer.pt()>30);
    double dphi=fabs(deltaPhi(box.lead.phi(),box.trailer.phi()));
    double mt=TMath::Sqrt(2*box.lead.pt()*box.trailer.pt()*(1-TMath::Cos(dphi)));
    box.acceptMass = (mt>30);
  }
  else if(isoPhotons.size()==1){
    box.lead    = photons[ isoPhotons[0] ]; 
    box.trailer = data::PhysicsObject_t(0,0,0,0);
    box.vCand   = photons[ isoPhotons[0] ]; 
    box.cat="g";
    box.accept     = true;
    box.acceptMass = true;
  }
  else if(nonIsoLeptons.size()>=2){
    box.lead      = leptons[ nonIsoLeptons[0] ] ;
    box.trailer   = leptons[ nonIsoLeptons[1] ] ;
    if( abs(box.lead.get("id")*box.trailer.get("id"))==11*11 )      box.cat="nonisoee";
    else if( abs(box.lead.get("id")*box.trailer.get("id"))==13*13 ) box.cat="nonisomumu";
    else                                                    box.cat="";
    box.vCand      = box.lead+box.trailer;
    box.accept     = (box.lead.pt()>20 && box.trailer.pt()>20 && box.cat!="");
    box.acceptMass = (fabs(box.vCand.mass()-91)<15);

  }
  else if(nonIsoLeptons.size()==1){
    box.lead      = leptons[ nonIsoLeptons[0] ] ;
    box.trailer   = data::PhysicsObject_t(met.px(),met.py(),0,met.pt());
    box.vCand     = box.lead+box.trailer;
    if( abs(box.lead.get("id"))==11 )      { box.cat="nonisoe";  box.accept=(box.lead.pt()>30); }
    else if( abs(box.lead.get("id"))==13 ) { box.cat="nonisomu"; box.accept=(box.lead.pt()>26); }
    box.accept    &= (box.trailer.pt()>30);
    double dphi=fabs(deltaPhi(box.lead.phi(),box.trailer.phi()));
    double mt=TMath::Sqrt(2*box.lead.pt()*box.trailer.pt()*(1-TMath::Cos(dphi)));
    box.acceptMass = (mt>30);
  }
  else if(nonIsoPhotons.size()==1){
    box.lead    = photons[ nonIsoPhotons[0] ]; 
    box.trailer = photons[ nonIsoPhotons[0] ]; 
    box.vCand   = data::PhysicsObject_t(0,0,0,0);
    box.cat="nonisog";
    box.accept=true;
    box.acceptMass=true;
  }

  //common requirements
  box.acceptKin = (box.vCand.pt()>50);
  box.acceptRap = (fabs(box.vCand.Rapidity())<1.4442);

  //fix me for W channels solve the pZ(neutrino) equation imposing mW and derive the full W kinematics

  return box;
}





//
bool passPhaseSpaceAcceptance(data::PhysicsObjectCollection_t &leptons, data::PhysicsObjectCollection_t &jets,bool isReco=false){
  if(leptons.size()<2 || jets.size()<2) return false;

  bool passLeptons(true);
  LorentzVector sumL(0,0,0,0);
  for(size_t i=0; i<2; i++){
    float pt=leptons[i].pt();
    float eta=leptons[i].eta();
    sumL+=leptons[i]; 
    if(isReco) passLeptons&=(fabs(eta)<2.5 && fabs(pt)>20); 
    else       passLeptons&=(fabs(eta)<2.5 && fabs(pt)>15); 
  }
  bool passDilepton( fabs(sumL.mass()-91)<30 );
  if(isReco) passDilepton=( fabs(sumL.mass()-91)<15);
  
  bool passJets(true);
  LorentzVector sumJ(0,0,0,0);
  for(size_t i=0; i<2; i++){
    float pt=jets[i].pt();
    float eta=jets[i].eta();
    sumJ+=jets[i]; 
    if(isReco) passJets&=(fabs(eta)<4.7 && fabs(pt)>30); 
    else       passJets&=(fabs(eta)<5.0 && fabs(pt)>25); 
  }
  bool passDijet( sumJ.mass()>120 );
  
  //build the final result
  //bool toReturn(passJets && passDijet);
  //bool toReturn(passJets && passDijet && passLeptons);
  bool toReturn(passJets && passDijet && passLeptons && passDilepton);
  return toReturn;
}


//
bool passGenAcceptance(data::PhysicsObjectCollection_t &gen,bool tauFilt=false){
  
  data::PhysicsObjectCollection_t jets,leptons;
  for(size_t i=0; i<gen.size(); i++){
    int status=gen[i].get("status");
    if(status!=3) continue;
    int pid=gen[i].get("id");
    if(fabs(pid)<6)                    jets.push_back(gen[i]);
    if(fabs(pid)==11 || fabs(pid)==13) leptons.push_back(gen[i]);
  }

  //select only ee or mumu events
  if(tauFilt) return (leptons.size()==2);

  //only last two quarks are outgoing
  std::reverse(jets.begin(),jets.end());
  return passPhaseSpaceAcceptance(leptons,jets);
}




//
float getAngle(LorentzVector &a, LorentzVector &b)
{
  TVector3 mom1(a.px(),a.py(),a.pz());
  TVector3 mom2(b.px(),b.py(),b.pz());
  double cosine = mom1.Dot(mom2)/(mom1.Mag()*mom2.Mag());
  return acos(cosine);
}

//
std::vector<TString> getDijetCategories(double mjj,double etajj,std::vector<TString> &curTags, TString &mjjCat)
{
  std::vector<TString> mjjCats;
  if(mjj<250)               mjjCats.push_back("mjjq016");
  if(mjj>=250 && mjj<350)   mjjCats.push_back("mjjq033");
  if(mjj>=350 && mjj<450)   mjjCats.push_back("mjjq049");
  if(mjj>=450 && mjj<550)   mjjCats.push_back("mjjq066");
  if(mjj>=550 && mjj<750)   mjjCats.push_back("mjjq083");
  if(mjj>=750 && mjj<1000)  mjjCats.push_back("mjjq092");
  if(mjj>=1000)             mjjCats.push_back("mjjq100");
  if(mjj>=1250)             mjjCats.push_back("highmjj");
  if(mjj>=750)              mjjCats.push_back("mjjgt092");  
  mjjCat=mjjCats[0];

  //include new tags
  std::vector<TString> selTags;
  for(size_t i=0; i<curTags.size(); i++)
    {
      TString itag=curTags[i];
      selTags.push_back(itag);
      for(size_t j=0; j<mjjCats.size(); j++)
	selTags.push_back(itag+mjjCats[j]);
    }

  //all done here
  return selTags;
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

  float minJetPtToApply(30);

  // ##############################
  // configure the process
  // ##############################
  const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");

  //input and output configuration
  std::vector<std::string> urls=runProcess.getParameter<std::vector<std::string> >("input");
  TString url = TString(urls[0]);
  TString baseDir    = runProcess.getParameter<std::string>("dirName");
  TString outFileUrl(gSystem->BaseName(url));
  outFileUrl.ReplaceAll(".root","");
  TString outdir=runProcess.getParameter<std::string>("outdir");
  TString outUrl( outdir );
  gSystem->Exec("mkdir -p " + outUrl);

  bool isMC = runProcess.getParameter<bool>("isMC");
  double xsec = runProcess.getParameter<double>("xsec");
  bool isPhotonPD(!isMC && url.Contains("Photon"));
  bool isDoubleElePD(!isMC && url.Contains("DoubleEle"));
  bool isDoubleMuPD(!isMC && url.Contains("DoubleMu"));
  bool isSingleElePD(!isMC && url.Contains("SingleEle"));  
  bool isSingleMuPD(!isMC && url.Contains("SingleMu"));  
  bool isV0JetsMC(isMC && (url.Contains("DYJetsToLL_50toInf") || url.Contains("WJets")));
  //bool isSignal(isMC && (url.Contains("VBFNLO") || url.Contains("lljj") || url.Contains("ajj") || url.Contains("lvjj") ) );

  //jet energy scale and uncertainties 
  TString jecDir = runProcess.getParameter<std::string>("jecDir");
  gSystem->ExpandPathName(jecDir);
  FactorizedJetCorrector *jesCor        = utils::cmssw::getJetCorrector(jecDir,isMC);
  JetCorrectionUncertainty *totalJESUnc = new JetCorrectionUncertainty((jecDir+"/MC_Uncertainty_AK5PFchs.txt").Data());

  //lepton efficiencies and corrections
  LeptonEfficiencySF lepEff;
  MuScleFitCorrector *muCor=getMuonCorrector(jecDir,url);

  //photon efficiencies


  //##############################################
  //########    INITIATING HISTOGRAMS     ########
  //##############################################
  SmartSelectionMonitor mon;
  
  TH1F* Hcutflow  = (TH1F*) mon.addHistogram(  new TH1F ("cutflow"    , "cutflow"    ,6,0,6) ) ;

  mon.addHistogram( new TH1F("nup",";NUP;Events",10,0,10) );
  mon.addHistogram( new TH1F("nupfilt",";NUP;Events",10,0,10) );
 
  TH1 *h=mon.addHistogram( new TH1F ("eventflow", ";;Events", 5,0,5) );
  h->GetXaxis()->SetBinLabel(1,"Box assign.");
  h->GetXaxis()->SetBinLabel(2,"M_{V}");
  h->GetXaxis()->SetBinLabel(3,"p_{T}^{V}>50");
  h->GetXaxis()->SetBinLabel(4,"#eta^{V}<1.4442");
  h->GetXaxis()->SetBinLabel(5,"#geq 2 jets"); 

  //pileup control
  mon.addHistogram( new TH1F( "nvtx",";Vertices;Events",50,0,50) ); 
  mon.addHistogram( new TH1F( "nvtxraw",";Vertices;Events",50,0,50) ); 

  //lepton control
  mon.addHistogram( new TH1F( "leadpt",    ";Transverse momentum [GeV];Events", 50,0,500) );
  mon.addHistogram( new TH1F( "trailerpt", ";Transverse momentum [GeV];Events", 50,0,500) );

  //boson control
  mon.addHistogram( new TH1F( "qt",      ";Transverse momentum [GeV];Events",500,0,1500));
  mon.addHistogram( new TH1F( "qeta",    ";Pseudo-rapidity;Events", 50,0,10) );
  mon.addHistogram( new TH1F( "qy",      ";Rapidity;Events", 50,0,6) );
  mon.addHistogram( new TH1F( "qmass",   ";Mass or transverse mass [GeV];Events", 100,0,500) );


  //balance histograms
  for(size_t ireg=0; ireg<3; ireg++){
    TString regStr("lt15collinear");
    if(ireg==1) regStr="gt50back2back";
    if(ireg==2) regStr="lt50back2back";
    mon.addHistogram( new TH1F("recoilbalance"+regStr, "; p_{T}(jet)/p_{T}; Jets", 100,0,5) );
    mon.addHistogram( new TH2F("recoilbalancevseta"+regStr, "; #eta(jet); <p_{T}(jet)/p_{T}>", 50,0,5,100,0,5) );
    TH2 *idH=(TH2 *)mon.addHistogram( new TH2F("recoilbalanceid"+regStr, "; Pseudo-rapidity; ID", 50,0,5,4,0,4) );
    idH->GetYaxis()->SetBinLabel(1,"no id");  
    idH->GetYaxis()->SetBinLabel(2,"PF");
    idH->GetYaxis()->SetBinLabel(3,"PU");
    idH->GetYaxis()->SetBinLabel(4,"PF+PU");
    if(ireg==0)
      {
	mon.addHistogram( (TH2 *)idH->Clone("truejetsid") );
	mon.addHistogram( (TH2 *)idH->Clone("pujetsid") );
      }
  }
  
  //jet control
  mon.addHistogram( new TH1F("jetpt"       , ";p_{T} [GeV];Events",50,0,400) );
  mon.addHistogram( new TH1F("jeteta"       , ";|#eta|;Events",25,0,5) );
  h=mon.addHistogram( new TH1F ("njets", ";Jet multiplicity;Events", 5,0,5) );
  h->GetXaxis()->SetBinLabel(1,"=0 jets");
  h->GetXaxis()->SetBinLabel(2,"=1 jets");
  h->GetXaxis()->SetBinLabel(3,"=2 jets");
  h->GetXaxis()->SetBinLabel(4,"=3 jets");
  h->GetXaxis()->SetBinLabel(5,"#geq 4 jets"); 

  //vbf control
  mon.addHistogram( new TH1F("vbfcandjet1pt"     , ";Leading jet p_{T} [GeV];Events",     50,0,1000) );
  mon.addHistogram( new TH1F("vbfcandjet2pt"     , ";Trailer jet p_{T} [GeV];Events",     50,0,500) );
  mon.addHistogram( new TH1F("vbfcandjet1eta"    , ";Forward jet #eta;Events",                                 25,0,5) );
  mon.addHistogram( new TH1F("vbfcandjet2eta"    , ";Central jet #eta;Events",                                 25,0,5) );
  mon.addHistogram( new TH1F("vbfcandjetdeta"    , ";Dijet pseudo-rapidity distance (#Delta#eta);Events",   50,0,10) );  
  Double_t mjjaxis[32];
  mjjaxis[0]=0.01;
  for(size_t i=1; i<20; i++)  mjjaxis[i]   =50*i;        //0-1000
  for(size_t i=0; i<5; i++)   mjjaxis[20+i]=1000+100*i; //1000-1500
  for(size_t i=0; i<=5; i++)   mjjaxis[25+i]=1500+300*i; //1500-5000  
  mjjaxis[31]=5000;
  mon.addHistogram( new TH1F("vbfmjj"            , ";Dijet invariant mass [GeV];Events",               31,mjjaxis) );
  mon.addHistogram( new TH1F("vbfdphijj"         , ";Dijet azimuthal opening (#Delta#phi_{jj});Events", 20,0,3.5) );
  mon.addHistogram( new TH1F("vbfystar"          , ";y^{*}_{ll}=#eta_{ll}-(#eta_{j1}+#eta_{j2})/2;Events",       10,0,5) );
  mon.addHistogram( new TH1F("vbfpt"             , ";Dijet p_{T} [GeV];Dijets",       50,0,500) );
  mon.addHistogram( new TH1F("vbfspt"            , ";Relative balance (#Delta^{rel}_{p_{T}});Events",                   50,0,1) );
  mon.addHistogram( new TH1F("vbfhardpt"         , ";Hard process p_{T} [GeV];Events",                   25,0,250) );

  
  //##############################################
  //######## GET READY FOR THE EVENT LOOP ########
  //##############################################

  //open the file and get events tree
  DataEventSummaryHandler evSummary;
  TFile *file = TFile::Open(url);
  printf("Looping on %s\n",url.Data());
  if(file==0) return -1;
  if(file->IsZombie()) return -1;
  if( !evSummary.attach( (TTree *) file->Get(baseDir+"/data") , false) ) { file->Close();  return -1; }
  const Int_t totalEntries= evSummary.getEntries();
 
  //MC normalization (to 1/pb)
  float cnorm=1.0;
  if(isMC){
    TH1F* cutflowH = (TH1F *) file->Get(baseDir+"/cutflow");
    if(cutflowH) cnorm=cutflowH->GetBinContent(1);
    printf("cnorm = %f\n",cnorm);
  }
  Hcutflow->SetBinContent(1,cnorm);
  //double xsecWeight = xsec;
  //if(!isMC) xsecWeight=1.0;

  //pileup weighting
  std::vector<double> dataPileupDistributionDouble = runProcess.getParameter< std::vector<double> >("datapileup");
  std::vector<float> dataPileupDistribution; for(unsigned int i=0;i<dataPileupDistributionDouble.size();i++){dataPileupDistribution.push_back(dataPileupDistributionDouble[i]);}
  std::vector<float> mcPileupDistribution;
  if(isMC){
    TString puDist(baseDir+"/pileup");
    TH1F* histo = (TH1F *) file->Get(puDist);
    if(!histo) std::cout<<"pileup histogram is null!!!\n";
    for(int i=1;i<=histo->GetNbinsX();i++){mcPileupDistribution.push_back(histo->GetBinContent(i));}
    delete histo;
  }
  while(mcPileupDistribution.size()<dataPileupDistribution.size())  mcPileupDistribution.push_back(0.0);
  while(mcPileupDistribution.size()>dataPileupDistribution.size())dataPileupDistribution.push_back(0.0);

  gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE
  edm::LumiReWeighting *LumiWeights= isMC ? new edm::LumiReWeighting(mcPileupDistribution,dataPileupDistribution): 0;
  utils::cmssw::PuShifter_t PuShifters;
  if(isMC) { PuShifters=utils::cmssw::getPUshifters(dataPileupDistribution,0.05); }

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
  
      //load the event content from tree
      evSummary.getEntry(iev);
      DataEventSummary &ev = evSummary.getEvent();
      if(!isMC && duplicatesChecker.isDuplicate( ev.run, ev.lumi, ev.event) ) { nDuplicates++; continue; }

      //gen level filter for V+njets
      if(isV0JetsMC){
	mon.fillHisto("nup","",ev.nup,1);
	if(ev.nup>5) continue;
	mon.fillHisto("nupfilt","",ev.nup,1);
      }

      //pileup weight
      float weight = 1.0;
      double TotalWeight_plus = 1.0;
      double TotalWeight_minus = 1.0;
      float puWeight(1.0);
      if(isMC){
        puWeight          = LumiWeights->weight(ev.ngenITpu);
	weight            = puWeight;
        TotalWeight_plus  = PuShifters[utils::cmssw::PUUP]->Eval(ev.ngenITpu);
        TotalWeight_minus = PuShifters[utils::cmssw::PUDOWN]->Eval(ev.ngenITpu);
      }
      Hcutflow->Fill(1,1);
      Hcutflow->Fill(2,weight);
      Hcutflow->Fill(3,weight*TotalWeight_minus);
      Hcutflow->Fill(4,weight*TotalWeight_plus);
      Hcutflow->Fill(5,weight);

      //check trigger bits and prescale
      float triggerPrescale(1.0),triggerThreshold(0);
      bool eeTrigger          = ev.t_bits[0];
      bool mumuTrigger        = (ev.t_bits[2] || ev.t_bits[3]);
      bool muTrigger          = ev.t_bits[6];
      bool eTrigger           = ev.t_bits[13];
      bool gTrigger(false);
      for(size_t itrig=10; itrig>=7; itrig--)
	{
	  if(!ev.t_bits[itrig]) continue;
	  gTrigger=true;
	  triggerPrescale=ev.t_prescale[itrig];
	  if(itrig==10) triggerThreshold=90;
	  if(itrig==9)  triggerThreshold=75;
	  if(itrig==8)  triggerThreshold=50;
	  if(itrig==7)  triggerThreshold=36;
	  break;
	}
      if(!isMC)
	{
	  if(!isDoubleElePD && !isPhotonPD) gTrigger=false;
	  if(!isDoubleElePD)                eeTrigger=false;
	  if(!isDoubleMuPD)                 mumuTrigger=false;
	  if(!isSingleMuPD)                 muTrigger=false;
	  if(!isSingleElePD)                eTrigger=false;
	}
   
      //
      // photon selection
      //
      data::PhysicsObjectCollection_t photons=evSummary.getPhysicsObject(DataEventSummaryHandler::PHOTONS);
      data::PhysicsObjectCollection_t selPhotons;
      for(size_t ipho=0; ipho<photons.size(); ipho++)
	{
	  double pt=photons[ipho].pt();
	  double eta=photons[ipho].getVal("sceta");
	  
	  //if systematics are active loosen the selection to the medium working point
	  Int_t idbits( photons[ipho].get("id") );
	  bool hasTightPhotonId( (idbits >> 2 ) & 0x1 );
	  double gIso    = photons[ipho].getVal("gIso03");
	  double gArea   = utils::cmssw::getEffectiveArea(22,eta,3,"gIso");	      
	  double chIso   = photons[ipho].getVal("chIso03");
	  double chArea  = utils::cmssw::getEffectiveArea(22,eta,3,"chIso");
	  double nhIso   = photons[ipho].getVal("nhIso03");
	  double nhArea  = utils::cmssw::getEffectiveArea(22,eta,3,"nhIso");
	  
	  //select the photon
	  if(pt<triggerThreshold || fabs(eta)>1.4442 ) continue;
	  bool passId(true);
	  if( photons[ipho].getVal("r9")<0.9 ) passId=false;
	  if(!hasTightPhotonId) passId=false;
	  if(!passId) continue;
	  bool passIso(true);
	  passIso &= (TMath::Max(chIso-chArea*ev.rho,0.0) < 0.7); 
	  passIso &= (TMath::Max(nhIso-nhArea*ev.rho,0.0) < 0.4+0.04*pt); 
	  passIso &= (TMath::Max(gIso-gArea*ev.rho,  0.0) < 0.5+0.005*pt); 
	  photons[ipho].setFlag("isIso",passIso);
	  selPhotons.push_back(photons[ipho]);
	}
    

      //
      // LEPTON ANALYSIS
      //
      data::PhysicsObjectCollection_t leptons=evSummary.getPhysicsObject(DataEventSummaryHandler::LEPTONS);
      data::PhysicsObjectCollection_t selLeptons;
      for(size_t ilep=0; ilep<leptons.size(); ilep++)
	{
	  bool passKin(true),passId(true),passIso(true);
	  int lid=leptons[ilep].get("id");

	  //apply muon corrections
	  if(lid==13 && muCor){
	    TLorentzVector p4(leptons[ilep].px(),leptons[ilep].py(),leptons[ilep].pz(),leptons[ilep].energy());
	    muCor->applyPtCorrection(p4 , lid<0 ? -1 :1 );
	    if(isMC) muCor->applyPtSmearing(p4, lid<0 ? -1 : 1, false);
	    leptons[ilep].SetPxPyPzE(p4.Px(),p4.Py(),p4.Pz(),p4.E());
	  }

	  //no need for charge info any longer
	  lid=abs(lid);
	  TString lepStr( lid==13 ? "mu" : "e");

	  //veto nearby photon (loose electrons are many times photons...)
	  double minDRlg(9999.);
	  for(size_t ipho=0; ipho<selPhotons.size(); ipho++)
	    minDRlg=TMath::Min(minDRlg,deltaR(leptons[ilep],selPhotons[ipho]));
	  if(minDRlg<0.1) continue;
	  
	  //kinematics
	  float leta = lid==11 ? leptons[ilep].getVal("sceta") : leptons[ilep].eta();
	  if(leptons[ilep].pt()<15)                   passKin=false;
	  if(leta> (lid==11 ? 2.5 : 2.4) )            passKin=false;
	  if(lid==11 && (leta>1.4442 && leta<1.5660)) passKin=false;

	  //id
	  Int_t idbits = leptons[ilep].get("idbits");
	  if(lid==11){
	    if(leptons[ilep].getFlag("isconv"))              passId=false;
	    bool isLoose = ((idbits >> 4) & 0x1);
	    if(!isLoose)                                     passId=false;
 	  }
	  else{
	    bool isLoose    = ((idbits >> 8) & 0x1);
	    if(!isLoose)                                     passId=false;
	  }

	  //isolation
	  Float_t gIso    = leptons[ilep].getVal(lid==11 ? "gIso03"    : "gIso04");
	  Float_t chIso   = leptons[ilep].getVal(lid==11 ? "chIso03"   : "chIso04");
	  Float_t puchIso = leptons[ilep].getVal(lid==11 ? "puchIso03" : "puchIso04");  
	  Float_t nhIso   = leptons[ilep].getVal(lid==11 ? "nhIso03"   : "nhIso04");
	  float relIso= lid==11 ?
	    (TMath::Max(nhIso+gIso-ev.rho*utils::cmssw::getEffectiveArea(11,leptons[ilep].getVal("sceta")),Float_t(0.))+chIso)/leptons[ilep].pt() :
	    (TMath::Max(nhIso+gIso-0.5*puchIso,0.)+chIso)/leptons[ilep].pt()
	    ;
	  if(lid==11){
	    if(relIso>0.15)                                passIso=false;
	  }
	  else{
	    if(relIso>0.20)                                passIso=false;
	  }
	  leptons[ilep].setFlag("isIso",passIso);

	  if(!passId || !passKin) continue;
	  selLeptons.push_back(leptons[ilep]);
	}
      std::sort(selLeptons.begin(), selLeptons.end(), data::PhysicsObject_t::sortByPt);

      //get JET/MET objects
      data::PhysicsObjectCollection_t jets=evSummary.getPhysicsObject(DataEventSummaryHandler::JETS);
      utils::cmssw::updateJEC(jets,jesCor,totalJESUnc,ev.rho,ev.nvtx,isMC);
      data::PhysicsObjectCollection_t recoMet=evSummary.getPhysicsObject(DataEventSummaryHandler::MET);
      std::vector<LorentzVector> met=utils::cmssw::getMETvariations(recoMet[0],jets,selLeptons,isMC);

      
      //get the analysis box and check trigger compatibility
      AnalysisBox_t box=assignBox(selLeptons,selPhotons,met[0]);
      if( ( box.cat=="ee"   || box.cat=="nonisoee")   && !eeTrigger )   box.accept=false;
      if( ( box.cat=="e"    || box.cat=="nonisoe")    && !eTrigger )    box.accept=false;
      if( ( box.cat=="mumu" || box.cat=="nonisomumu") && !mumuTrigger ) box.accept=false;
      if( ( box.cat=="mu"   || box.cat=="nonisomu")   && !muTrigger )   box.accept=false;
      if( ( box.cat=="g"    || box.cat=="nonisog")    && !gTrigger )    box.accept=false;
      if(!box.accept) continue;


      //apply data/mc correction factors
      bool boxHasTrailer(box.cat.EndsWith("ee") || box.cat.EndsWith("mumu"));
      std::vector<TString> chTags(1,box.cat);
      if(!box.cat.EndsWith("g"))
	{
	  int leadId(abs(box.lead.get("id")));
	  weight *= isMC ? lepEff.getLeptonEfficiency( box.lead.pt(), box.lead.eta(), leadId,  leadId ==11 ? "loose" : "loose" ).first : 1.0;	  
	  if(boxHasTrailer){
	    int trailerId(abs(box.trailer.get("id")));
	    weight *= isMC ? lepEff.getLeptonEfficiency( box.trailer.pt(), box.trailer.eta(), trailerId,  trailerId ==11 ? "loose" : "loose" ).first : 1.0;	  
	  }
	}
      else{
	weight *= triggerPrescale;
      }

      
      //select the jets for the box
      int njets(0),njetsNoId(0); 
      data::PhysicsObjectCollection_t selJets, selJetsNoId;
      for(size_t ijet=0; ijet<jets.size(); ijet++) 
	{
	  double jeta=fabs(jets[ijet].eta());
	  if(jets[ijet].pt()<15 || jeta>4.7 ) continue;

	  //mc truth for this jet
	  const data::PhysicsObject_t &genJet=jets[ijet].getObject("genJet");
	  TString jetType( genJet.pt()>0 ? "truejetsid" : "pujetsid" );
	  
	  //cross-clean with objects in the box
	  double drLeadJet( deltaR(jets[ijet], box.lead) );
	  double drTrailerJet( boxHasTrailer ? deltaR(jets[ijet], box.trailer) : 9999. );
	  if( drLeadJet<0.4 || drTrailerJet<0.4) continue;
	  
	  //jet id
	  Int_t idbits=jets[ijet].get("idbits");
	  bool passPFloose( ((idbits>>0) & 0x1));
	  int simplePuId( ( idbits >>7 ) & 0xf );
	  bool passLooseSimplePuId(  ( simplePuId >> 2) & 0x1);
	  if(jets[ijet].pt()>30)
	    {
	      mon.fillHisto(jetType,"",jeta,0);
	      if(passPFloose)                        mon.fillHisto(jetType,"",jeta,1);
	      if(passLooseSimplePuId)                mon.fillHisto(jetType,"",jeta,2);
	      if(passPFloose && passLooseSimplePuId) mon.fillHisto(jetType,"",jeta,3);
	    }
	  
	  selJetsNoId.push_back(jets[ijet]);
	  if(jets[ijet].pt()>minJetPtToApply) njetsNoId++;
	  if(passPFloose && passLooseSimplePuId){
	    selJets.push_back(jets[ijet]);
	    if(jets[ijet].pt()>minJetPtToApply) njets++;
	  }
	}
      std::sort(selJets.begin(), selJets.end(), data::PhysicsObject_t::sortByPt);
      
      //analyze dijets
      float maxPt(0), minPt(0), maxEta(0), minEta(0), maxAbsEta(0), minAbsEta(0);
      float dphijj(0), detajj(0), mjj(0), ptjj(0), spt(0);
      float hardpt(0), ystar(0);
      if(njets>=2)
	{
	  LorentzVector jet1=selJets[0];
	  LorentzVector jet2=selJets[1];
	  maxPt=max(jet1.pt(),jet2.pt());
	  minPt=min(jet1.pt(),jet2.pt());
	  maxAbsEta=max(fabs(jet1.eta()),fabs(jet2.eta()));
	  minAbsEta=min(fabs(jet1.eta()),fabs(jet2.eta()));
	  maxEta=max(jet1.eta(),jet2.eta());
	  minEta=min(jet1.eta(),jet2.eta());
	  detajj=fabs(maxEta-minEta);
	  dphijj=deltaPhi(jet1.phi(),jet2.phi());
	    
	  LorentzVector vbfSyst=jet1+jet2;
	  mjj=vbfSyst.mass();
	  ptjj=vbfSyst.pt();
	  spt=ptjj/(jet1.pt()+jet2.pt());

	  LorentzVector hardSyst=vbfSyst+box.vCand;
	  hardpt=hardSyst.pt();
	  ystar=box.vCand.Rapidity()-0.5*(jet1.Rapidity()+jet2.Rapidity());
	}

            
      //
      // NOW FOR THE CONTROL PLOTS
      //
      //start analysis
      for(size_t ich=0; ich<chTags.size(); ich++)
	{
	  std::vector<TString> tags(1,chTags[ich]);
	  
	  mon.fillHisto("eventflow",tags,0,weight);
	  if(box.acceptMass){
	    mon.fillHisto("eventflow",tags,1,weight);
	    if(box.acceptKin){
	      mon.fillHisto("eventflow",tags,2,weight);
	      if(box.acceptRap){
		mon.fillHisto("eventflow",tags,3,weight);
		if(njets>1){
		  mon.fillHisto("eventflow",tags,4,weight);
		}
	      }
	    }
	  }
	
	  mon.fillHisto("leadpt",    tags, box.lead.pt(), weight);  
	  mon.fillHisto("trailerpt", tags, box.trailer.pt(), weight);  
	  
	  mon.fillHisto("qmass",    tags, box.vCand.mass(), weight);  
	  if(box.acceptMass){
	
	    mon.fillHisto("nvtx"     ,   tags, ev.nvtx,      weight);
	    mon.fillHisto("nvtxraw"  ,   tags, ev.nvtx,      weight/puWeight);
    
	    //balance control
	    if(njetsNoId==1)
	      {
		//set as pu if no matched gen jet
		//bool isPUjet( selJetsNoId[0].getObject("genJet").pt()==0 ); 
		
		//kinematics
		float balance = selJetsNoId[0].pt()/box.vCand.pt();
		float dphi    = fabs( deltaPhi(selJetsNoId[0].phi(),box.vCand.phi()) );
		TString regStr("");
		if(dphi<1   && box.vCand.pt()<15) regStr="lt15collinear";
		if(dphi>2.7 && box.vCand.pt()>50) regStr="gt50back2back";
		if(dphi>2.7 && box.vCand.pt()<=50) regStr="lt50back2back";
		if(regStr!="")
		  {
		    //ids
		    Int_t idbits=selJetsNoId[0].get("idbits");
		    bool passPFloose ( ((idbits>>0) & 0x1) );
		    int simplePuId( ( idbits >>7 ) & 0xf );
		    bool passLooseSimplePuId(  ( simplePuId >> 2) & 0x1);

		    mon.fillHisto("recoilbalance"+regStr,      tags,balance, weight);
		    mon.fillHisto("recoilbalancevseta"+regStr, tags,fabs(selJetsNoId[0].eta()), balance, weight);
		    mon.fillHisto("recoilbalanceid"+regStr,    tags,fabs(selJetsNoId[0].eta()),0, weight);
		    if(passPFloose)                         mon.fillHisto("recoilbalanceid"+regStr,tags,fabs(selJetsNoId[0].eta()),1, weight);
		    if(passLooseSimplePuId)                 mon.fillHisto("recoilbalanceid"+regStr,tags,fabs(selJetsNoId[0].eta()),2, weight);
		    if(passPFloose && passLooseSimplePuId)  mon.fillHisto("recoilbalanceid"+regStr,tags,fabs(selJetsNoId[0].eta()),3, weight);
		  }
	      }
	    //end balance control

	    if(box.acceptKin){
	      mon.fillHisto("qeta"    , tags, fabs(box.vCand.eta()),      weight);
	      mon.fillHisto("qy"      , tags, fabs(box.vCand.Rapidity()), weight);
	    }

	    if(box.acceptKin && box.acceptRap){

	      mon.fillHisto("njets",tags, njets, weight);
	      for(size_t ijet=0; ijet<selJets.size(); ijet++)
		{
		  float pt=selJets[ijet].pt();
		  if(pt<minJetPtToApply) continue;
		  mon.fillHisto("jetpt",  tags, pt, weight);
		  mon.fillHisto("jeteta", tags, fabs(selJets[ijet].eta()), weight);
		}
	      
	      //signal region
	      if(njets>=2)
		{
		  TString mjjCat("");
		  std::vector<TString> selTags;
		  selTags = getDijetCategories(mjj,detajj,tags,mjjCat);

		  mon.fillHisto("qt"          ,       selTags, box.vCand.pt(),      weight);
		  mon.fillHisto("vbfcandjetpt",       selTags, maxPt,weight);
		  mon.fillHisto("vbfcandjetpt",       selTags, minPt,weight);
		  mon.fillHisto("vbfcandjet1pt",      selTags, maxPt,weight);
		  mon.fillHisto("vbfcandjet2pt",      selTags, minPt,weight);
		  mon.fillHisto("vbfcandjet1eta",     selTags, maxAbsEta, weight);
		  mon.fillHisto("vbfcandjet2eta",     selTags, minAbsEta, weight);
		  mon.fillHisto("vbfcandjeteta",      selTags, maxAbsEta, weight);
		  mon.fillHisto("vbfcandjeteta",      selTags, minAbsEta, weight);
		  mon.fillHisto("vbfcandjetdeta",     selTags, detajj,weight);
		  mon.fillHisto("vbfmjj",             selTags, mjj,weight,true);
		  mon.fillHisto("vbfpt",              selTags, ptjj,weight);
		  mon.fillHisto("vbfspt",             selTags, spt,weight);
		  mon.fillHisto("vbfdphijj",          selTags, fabs(dphijj),weight);
		  mon.fillHisto("vbfhardpt",          selTags, hardpt,weight);
		  mon.fillHisto("vbfystar",           selTags, fabs(ystar),weight);
		}
	    }
	  }
	}
    }	
  printf("\n"); 
  file->Close();
  
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
}  





