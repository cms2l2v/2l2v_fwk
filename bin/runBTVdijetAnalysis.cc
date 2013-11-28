#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/llvv_fwk/interface/DataEventSummaryHandler.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h" 
#include "TGraph2D.h"

#include <fstream>
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
  const edm::ParameterSet &runProcess       = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");
  std::vector<std::string> urls=runProcess.getParameter<std::vector<std::string> >("input");
  TString url = TString(urls[0]);
  TString baseDir     = runProcess.getParameter<std::string>("dirName");
  TString jecDir      = runProcess.getParameter<std::string>("jecDir");
  bool isMC           = runProcess.getParameter<bool>("isMC");
  double xsec         = runProcess.getParameter<double>("xsec");
  TString out          = runProcess.getParameter<std::string>("outdir");
  std::vector<string>  weightsFile = runProcess.getParameter<std::vector<string> >("weightsFile");

 //jet energy scale uncertainties
  gSystem->ExpandPathName(jecDir);
  FactorizedJetCorrector *jesCor        = utils::cmssw::getJetCorrector(jecDir,isMC);
  JetCorrectionUncertainty *totalJESUnc = new JetCorrectionUncertainty((jecDir+"/MC_Uncertainty_AK5PFchs.txt").Data());
  
  //
  // check input file
  //
  TFile *inF = TFile::Open(url);
  if(inF==0) return -1;
  if(inF->IsZombie()) return -1;
  TString proctag=gSystem->BaseName(url);
  Ssiz_t pos=proctag.Index(".root");
  proctag.Remove(pos,proctag.Length());

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

  //x_b weights
  std::vector<TString> systVars(1,"");
  std::map<TString,TGraph *> xbWeights;
  if(isMC && weightsFile.size()){
    TString xbFileUrl(weightsFile[0]+"/xb_weights.root");
    gSystem->ExpandPathName(xbFileUrl);
    cout << "Reading x_b weights from " << xbFileUrl << endl;
    TFile *xbFile=TFile::Open(xbFileUrl);
    TH1 *nominal  = (TH1 *) xbFile->Get("pythia_Z2_nominal");
    TH1 *corcella = (TH1 *) xbFile->Get("pythia_Z2_corcella");
    TH1 *p11      = (TH1 *)  xbFile->Get("pythia_Z2_p11");
    p11->Divide(nominal);        xbWeights["p11"]=new TGraph(p11);            systVars.push_back("p11");
    corcella->Divide(nominal);   xbWeights["corcella"]=new TGraph(corcella);  systVars.push_back("corcella");
    xbFile->Close();
  }

  //book histograms
  SmartSelectionMonitor controlHistos;
  TH1 *Hcutflow=(TH1 *)controlHistos.addHistogram(  new TH1F ("cutflow"    , "cutflow"    ,6,0,6) ) ;
  for(int ibin=1; ibin<=6; ibin++) Hcutflow->SetBinContent(ibin,1);

  controlHistos.addHistogram( new TH1F("pthat", ";#hat{p}_{T} [GeV]; Events",100,0,1500) );
  controlHistos.addHistogram( new TH1F ("nvertices", "; Vertex multiplicity; Events", 50, 0.,50.) );

  //jet control
  TString jetFlavors[]={"","b","udsg","c","unmatched"};
  size_t nJetFlavors=sizeof(jetFlavors)/sizeof(TString);
  for(size_t ijet=1; ijet<=2; ijet++)
    {
      TString jetctr(""); jetctr += ijet;
      controlHistos.addHistogram( new TH1F ("jet"+jetctr, "; Jet #"+jetctr+" transverse momentum [GeV]; Jets", 100, 0.,1000.) );
      controlHistos.addHistogram( new TH1F ("jet"+jetctr+"eta", "; Jet #"+jetctr+" pseudo-rapidity; Jets", 30, 0.,3.) );
      TH1 *hflav=controlHistos.addHistogram( new TH1F ("jet"+jetctr+"flav", "; Jet #"+jetctr+" flavour; Jets", nJetFlavors, 0.,nJetFlavors) );
      for(int xbin=1; xbin<=hflav->GetXaxis()->GetNbins(); xbin++) hflav->GetXaxis()->SetBinLabel(xbin,jetFlavors[xbin-1]);
    }
  
  //tag and probe analysis for Sec Vtx
  const Double_t ptBins[]={30,35,40,45,50,55,60,65,70,80,90,100,120,140,160,180,200,250,350,400,500,750,1000};
  Int_t nPtbins=sizeof(ptBins)/sizeof(Double_t)-1;
  TString svxAlgos[]={"svx","ivf"};
  for(size_t ialgo=0; ialgo<2; ialgo++)
    {
      for(size_t iflav=0; iflav<nJetFlavors; iflav++)
	{
	  size_t nVars(1);
	  if(jetFlavors[iflav]=="b") nVars=systVars.size();
	  for(size_t ivar=0; ivar<nVars; ivar++)
	    {
	      TString prefix(systVars[ivar]+jetFlavors[iflav]);
	      controlHistos.addHistogram( new TH2F (prefix+"recoil"+svxAlgos[ialgo]+"mass",       "; SecVtx Mass [GeV]; Jets",                            50, 0.,10.,  nPtbins,ptBins) );
	      controlHistos.addHistogram( new TH2F (prefix+"recoil"+svxAlgos[ialgo]+"lxy",        "; SecVtx L_{xy} [cm]; Jets",                           100, 0.,10., nPtbins,ptBins) );	 
	      controlHistos.addHistogram( new TH2F (prefix+"recoil"+svxAlgos[ialgo]+"lxysig",     "; #sigma/L_{xy} ; Jets",                               100, 0.,5.,  nPtbins,ptBins) );	 
	      controlHistos.addHistogram( new TH2F (prefix+"recoil"+svxAlgos[ialgo]+"dr",         "; #Delta R(jet,SecVtx L_{xy}); Jets",       50, 0.,1.0,  nPtbins,ptBins) );	 
	      controlHistos.addHistogram( new TH2F (prefix+"recoil"+svxAlgos[ialgo]+"ptfrac",     "; p_{T}(SecVtx L_{xy})/p_{T}(jet); Jets",   50, 0.,2.0,  nPtbins,ptBins) );	 
	      controlHistos.addHistogram( new TH2F (prefix+"recoil"+svxAlgos[ialgo]+"neutemfrac", "; Neutral EM fraction; Jets",               50, 0.,1.0,  nPtbins,ptBins) );	 
	      controlHistos.addHistogram( new TH2F (prefix+"recoil"+svxAlgos[ialgo]+"chfrac",     "; Charged fraction; Jets",                  50, 0.,1.0,  nPtbins,ptBins) );	 
	      controlHistos.addHistogram( new TH2F (prefix+"recoil"+svxAlgos[ialgo]+"neuthadfrac","; Neutral Had fraction; Jets",              50, 0.,1.0,  nPtbins,ptBins) );	 
	      controlHistos.addHistogram( new TH2F (prefix+"recoil"+svxAlgos[ialgo]+"mufrac",     "; Muon fraction; Jets",                     50, 0.,1.0,  nPtbins,ptBins) );	 
	      controlHistos.addHistogram( new TH2F (prefix+"recoil"+svxAlgos[ialgo]+"ntracks", "; SecVtx track multiplicity; Jets",                    10, 0.,10,   nPtbins,ptBins) );	 
	      controlHistos.addHistogram( new TH2F (prefix+"recoil"+svxAlgos[ialgo]+"nmuons",  "; Soft muon multiplicity; Jets",                       3, 0.,3,     nPtbins,ptBins) );	 
	      controlHistos.addHistogram( new TH2F (prefix+"recoil"+svxAlgos[ialgo]+"dimuon",  "; Soft di-muon mass [GeV]; Jets",                      100, 0.,10,  nPtbins,ptBins) );	 
	      controlHistos.addHistogram( new TH2F (prefix+"recoil"+svxAlgos[ialgo]+"dimuondr",      "; #Delta R(jet,#mu#mu); Jets",                   50, 0.,1.0,  nPtbins,ptBins) );	 
	      controlHistos.addHistogram( new TH2F (prefix+"recoil"+svxAlgos[ialgo]+"dimuonptfrac",  "; p_{T}(#mu#mu)/p_{T}(jet); Jets",               50, 0.,2.0,  nPtbins,ptBins) );	 
	    }
	}
    }

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
      float weightNom(1.0);
      if(LumiWeights) 	weightNom  = LumiWeights->weight(ev.ngenITpu);

      bool hasTrigger(false);
      float triggerPrescale(1.0);
      // int triggerIdx(-1);
      for(int itrig=0; itrig<=4; itrig++)
	{
	  if(!ev.t_bits[itrig]) continue;
	  hasTrigger=true;
	  //triggerPrescale=ev.t_prescale[itrig];
	  //  triggerIdx=itrig;
	  break;
	}
      if(isMC) hasTrigger=true; 
      if(!hasTrigger) continue;
      
      //weight for the event
      float weight(weightNom*triggerPrescale);

      //require 2 jets in the event
      data::PhysicsObjectCollection_t jets=evSummary.getPhysicsObject(DataEventSummaryHandler::JETS);
      utils::cmssw::updateJEC(jets,jesCor,totalJESUnc,ev.rho,ev.nvtx,isMC);
      data::PhysicsObjectCollection_t selJets;
      for(size_t ijet=0; ijet<jets.size(); ijet++)
	{
	  //require to pass the loose id
	  Int_t idbits=jets[ijet].get("idbits");
	  bool passPFloose( ((idbits>>0) & 0x1));
	  if(!passPFloose) continue;
	  
	  //top candidate jets
	  if(jets[ijet].pt()<30 || fabs(jets[ijet].eta())>2.5 ) continue;
	  selJets.push_back(jets[ijet]);
	}
      sort(selJets.begin(),  selJets.end(),  data::PhysicsObject_t::sortByPt);
     if(selJets.size()!=2) continue;
      
     //require at least one soft muon
     data::PhysicsObjectCollection_t pf = evSummary.getPhysicsObject(DataEventSummaryHandler::PFCANDIDATES);
     std::vector<int> nTriggerSoftMuons(selJets.size(),0), nSoftMuons(selJets.size(),0);
     std::vector<LorentzVector> muonsP4(selJets.size(),LorentzVector(0,0,0,0));
     for(size_t ijet=0; ijet<selJets.size(); ijet++){
       size_t pfstart=selJets[ijet].get("pfstart");
       size_t pfend=selJets[ijet].get("pfend");
       if(pf.size()<pfstart || pf.size()<pfend-1) continue;
       for(size_t ipfn=pfstart; ipfn<=pfend; ipfn++)
	 {
	   int id=pf[ipfn].get("id");
	   if(abs(id)!=13) continue;
	   if(pf[ipfn].pt()<1) continue;
	   nSoftMuons[ijet]++;
	   muonsP4[ijet] += pf[ipfn];
	   if(pf[ipfn].pt()<5) continue;
	   nTriggerSoftMuons[ijet]++;
	 }
     }
     if(nTriggerSoftMuons[0]+nTriggerSoftMuons[1]==0) continue;
     
     //azimuthal angle
     float dphijj=deltaPhi(selJets[0].phi(),selJets[1].phi());
     if(fabs(dphijj)<2.7) continue;
     
     //balancing variable
     float balance=selJets[1].pt()/selJets[0].pt();
     if(balance<0.9) continue;

     //now decide which one is the tag and which one is the probe (give preference to higher pT jet for tag) 
     int tagJetIdx(0), probeJetIdx(1);
     if( nTriggerSoftMuons[0]==0 ) { tagJetIdx=1; probeJetIdx=0; }

     //tag categories
     bool tagHasCSVL(selJets[tagJetIdx].getVal("csv")>0.405), tagHasCSVM(selJets[tagJetIdx].getVal("csv")>0.783), tagHasCSVT(selJets[tagJetIdx].getVal("csv")>0.920);
     std::vector<TString> catsToFill;
     catsToFill.push_back("");
     if(tagHasCSVL) catsToFill.push_back("Ltag");
     if(tagHasCSVM) catsToFill.push_back("Mtag");
     if(tagHasCSVT) catsToFill.push_back("Ttag");

     //gen control
     controlHistos.fillHisto("pthat",catsToFill,ev.pthat,weight);
     
     //vertices control
     controlHistos.fillHisto("nvertices",catsToFill,ev.nvtx,weight);

     //jet control
     std::vector<int> jetFlav(selJets.size(),0);
     for(size_t ijet=0; ijet<2; ijet++)
       {
	 int idx(tagJetIdx);
	 if(ijet==1) idx=probeJetIdx;
	 
	 //kinematics
	 float pt( selJets[ idx ].pt() );
	 float eta( fabs( selJets[ idx ].eta() ) );
	  
	 //mc truth, if available 
	 int iflav(0);
	 if(isMC)
	   {
	     const data::PhysicsObject_t &genJet=selJets[idx].getObject("genJet");
	     if(genJet.pt()>0) {
	       jetFlav[ijet]=genJet.info.find("id")->second;
	       if(abs(jetFlav[ijet])==5)      iflav=1;
	       else if(abs(jetFlav[ijet])==4) iflav=3;
	       else                            iflav=2;
	     }
	     else                              { iflav=4; jetFlav[ijet]=9999; }
	   }
 
	 //fill the histograms
	 TString jetctr(""); jetctr+=(ijet+1);
	 controlHistos.fillHisto("jet"+jetctr,        catsToFill, pt,        weight);
	 controlHistos.fillHisto("jet"+jetctr+"eta",  catsToFill, fabs(eta), weight);
	 controlHistos.fillHisto("jet"+jetctr+"flav", catsToFill, iflav,     weight);
       }
           
     //secondary vertex characteristics evaluated for the probe
     float recoilPtNorm(TMath::Min(selJets[probeJetIdx].pt(),ptBins[nPtbins-1]));
     int probeFlav( jetFlav[probeJetIdx] );
     float probeXb(-1);
     TString probeFlavStr("udsg");
     //for b quarks match the B-hadron
     if(abs(probeFlav)==5)         {
       probeFlavStr="b";
       const data::PhysicsObject_t &genParton=selJets[probeJetIdx].getObject("gen");
       int genId=genParton.info.find("id")->second;
       if(abs(genId)==5 && genParton.pt()>0){
	 data::PhysicsObjectCollection_t mctruth=evSummary.getPhysicsObject(DataEventSummaryHandler::GENPARTICLES);
	 for(size_t imc=0; imc<mctruth.size(); imc++)
	   {
	     int id=mctruth[imc].get("id");
	     if(abs(id)<500) continue;
	     if(deltaR(mctruth[imc],genParton)>0.5) continue;
	     probeXb=mctruth[imc].pt()/genParton.pt();
	     break;
	   }
       }
     }
     else if(abs(probeFlav)==4)    probeFlavStr="c";
     else if(abs(probeFlav)==9999) probeFlavStr="unmatched";
     float probeEta = selJets[probeJetIdx].eta();
     if(fabs(probeEta)>1.1) continue;     
     int nmuons=nSoftMuons[probeJetIdx];
     LorentzVector mumu=muonsP4[probeJetIdx];
     for(size_t ialgo=0; ialgo<2; ialgo++)
       {
	 const data::PhysicsObject_t &svx = selJets[probeJetIdx].getObject(svxAlgos[ialgo]);
	 float lxy=svx.vals.find("lxy")->second;
	 if(lxy<=0) continue;
	 float lxyErr=svx.vals.find("lxyErr")->second;
	 int ntrk=svx.info.find("ntrk")->second;

	 controlHistos.fillHisto("recoil"+svxAlgos[ialgo]+"mass",    catsToFill, svx.mass(),                      recoilPtNorm, weight);
	 controlHistos.fillHisto("recoil"+svxAlgos[ialgo]+"lxy",     catsToFill, lxy,                             recoilPtNorm, weight);
	 controlHistos.fillHisto("recoil"+svxAlgos[ialgo]+"lxysig",  catsToFill, lxyErr/lxy,                      recoilPtNorm, weight);
	 controlHistos.fillHisto("recoil"+svxAlgos[ialgo]+"dr",      catsToFill, deltaR(selJets[probeJetIdx],svx),   recoilPtNorm, weight);
	 controlHistos.fillHisto("recoil"+svxAlgos[ialgo]+"ptfrac",  catsToFill, svx.pt()/selJets[probeJetIdx].pt(), recoilPtNorm, weight);
	 controlHistos.fillHisto("recoil"+svxAlgos[ialgo]+"neutemfrac",  catsToFill, selJets[probeJetIdx].getVal("neutEmFrac"),   recoilPtNorm, weight);
	 controlHistos.fillHisto("recoil"+svxAlgos[ialgo]+"chfrac",  catsToFill, selJets[probeJetIdx].getVal("chHadFrac"),        recoilPtNorm, weight);
	 controlHistos.fillHisto("recoil"+svxAlgos[ialgo]+"neuthadfrac",  catsToFill, selJets[probeJetIdx].getVal("neutHadFrac"), recoilPtNorm, weight);
	 controlHistos.fillHisto("recoil"+svxAlgos[ialgo]+"mufrac",  catsToFill, selJets[probeJetIdx].getVal("muFrac"),           recoilPtNorm, weight);
	 controlHistos.fillHisto("recoil"+svxAlgos[ialgo]+"ntracks", catsToFill, ntrk,                            recoilPtNorm, weight);
	 controlHistos.fillHisto("recoil"+svxAlgos[ialgo]+"nmuons",  catsToFill, nmuons,                          recoilPtNorm, weight);
	 if(nmuons==2)
	   {
	     controlHistos.fillHisto("recoil"+svxAlgos[ialgo]+"dimuon",       catsToFill, mumu.mass(),                        recoilPtNorm, weight);
	     controlHistos.fillHisto("recoil"+svxAlgos[ialgo]+"dimuonptfrac", catsToFill, mumu.pt()/selJets[probeJetIdx].pt(),   recoilPtNorm, weight);
	     controlHistos.fillHisto("recoil"+svxAlgos[ialgo]+"dimuondr",     catsToFill, deltaR(mumu,selJets[probeJetIdx]),     recoilPtNorm, weight);
	   }

	 if(!isMC) continue;
	 size_t nVars(1);
	 if(probeFlavStr=="b") nVars=systVars.size();
	 for(size_t ivar=0; ivar<nVars; ivar++)
	   {
	     TString prefix(systVars[ivar]+probeFlavStr);
	     float iweight(weight);
	     if(xbWeights.find( systVars[ivar] ) != xbWeights.end() && probeXb>0 ) iweight *= xbWeights[systVars[ivar]]->Eval(probeXb);

	     controlHistos.fillHisto(prefix+"recoil"+svxAlgos[ialgo]+"mass",   catsToFill, svx.mass(),                         recoilPtNorm,iweight);
	     controlHistos.fillHisto(prefix+"recoil"+svxAlgos[ialgo]+"lxy",    catsToFill, lxy,                                recoilPtNorm,iweight);
	     controlHistos.fillHisto(prefix+"recoil"+svxAlgos[ialgo]+"lxysig",  catsToFill, lxyErr/lxy,                      recoilPtNorm, weight);
	     controlHistos.fillHisto(prefix+"recoil"+svxAlgos[ialgo]+"dr",     catsToFill, deltaR(selJets[probeJetIdx],svx),   recoilPtNorm,iweight);
	     controlHistos.fillHisto(prefix+"recoil"+svxAlgos[ialgo]+"ptfrac", catsToFill, svx.pt()/selJets[probeJetIdx].pt(), recoilPtNorm,iweight);
	     controlHistos.fillHisto(prefix+"recoil"+svxAlgos[ialgo]+"ntracks", catsToFill, ntrk,                              recoilPtNorm, iweight);
	     controlHistos.fillHisto(prefix+"recoil"+svxAlgos[ialgo]+"nmuons",  catsToFill, nmuons,                            recoilPtNorm, iweight);
	     controlHistos.fillHisto(prefix+"recoil"+svxAlgos[ialgo]+"neutemfrac",  catsToFill, selJets[probeJetIdx].getVal("neutEmFrac"),   recoilPtNorm, iweight);
	     controlHistos.fillHisto(prefix+"recoil"+svxAlgos[ialgo]+"chfrac",  catsToFill, selJets[probeJetIdx].getVal("chHadFrac"),        recoilPtNorm, iweight);
	     controlHistos.fillHisto(prefix+"recoil"+svxAlgos[ialgo]+"neuthadfrac",  catsToFill, selJets[probeJetIdx].getVal("neutHadFrac"), recoilPtNorm, iweight);
	     controlHistos.fillHisto(prefix+"recoil"+svxAlgos[ialgo]+"mufrac",  catsToFill, selJets[probeJetIdx].getVal("muFrac"),           recoilPtNorm, iweight);
	     if(nmuons==2)
	       {
		 controlHistos.fillHisto(prefix+"recoil"+svxAlgos[ialgo]+"dimuon",       catsToFill, mumu.mass(),                        recoilPtNorm, iweight);
		 controlHistos.fillHisto(prefix+"recoil"+svxAlgos[ialgo]+"dimuonptfrac", catsToFill, mumu.pt()/selJets[probeJetIdx].pt(),   recoilPtNorm, iweight);
		 controlHistos.fillHisto(prefix+"recoil"+svxAlgos[ialgo]+"dimuondr",     catsToFill, deltaR(mumu,selJets[probeJetIdx]),     recoilPtNorm, iweight);
	       }
	   }
       }
    }
  inF->Close();
  
  //
  // save histos to local file
  //
  TString outUrl(out);
  gSystem->ExpandPathName(outUrl);
  gSystem->Exec("mkdir -p " + outUrl);
  outUrl += "/";
  outUrl += proctag + ".root";
  TFile *file=TFile::Open(outUrl, "recreate");
  controlHistos.Write();
  file->Close();
 }
