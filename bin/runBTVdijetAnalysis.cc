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


  //book histograms
  SmartSelectionMonitor controlHistos;
  TH1 *Hcutflow=(TH1 *)controlHistos.addHistogram(  new TH1F ("cutflow"    , "cutflow"    ,6,0,6) ) ;
  for(int ibin=1; ibin<=6; ibin++) Hcutflow->SetBinContent(ibin,1);

  controlHistos.addHistogram( new TH1F ("triggeridx", "; Trigger index; Events", 5, 0.,5.) );
  controlHistos.addHistogram( new TH1F("pthat", ";#hat{p}_{T} [GeV]; Events",100,0,1500) );
  controlHistos.addHistogram( new TH1F ("nvertices", "; Vertex multiplicity; Events", 50, 0.,50.) );
  controlHistos.addHistogram( new TH1F ("nvertices_unwgt", "; Vertex multiplicity; Events", 50, 0.,50.) );
 
  TH1 * h = controlHistos.addHistogram( new TH1F ("njets", "; Jet multiplicity; Events", 6, 0.,6.) );
  h->GetXaxis()->SetBinLabel(1,"=0 jets");
  h->GetXaxis()->SetBinLabel(2,"=1 jets");
  h->GetXaxis()->SetBinLabel(3,"=2 jets");
  h->GetXaxis()->SetBinLabel(4,"=3 jets");
  h->GetXaxis()->SetBinLabel(5,"=4 jets");
  h->GetXaxis()->SetBinLabel(6,"#geq 5 jets");

  //jet control
  controlHistos.addHistogram( new TH1F ("rawjetpt", ";Transverse momentum [GeV]; Jets", 100, 0.,1000.) );
  controlHistos.addHistogram( new TH1F ("jetpt", ";Transverse momentum [GeV]; Jets",    100, 0.,1000.) );
  controlHistos.addHistogram( new TH1F ("jeteta", ";Pseudo-rapidity; Events", 30, 0.,3.) );
  for(size_t ijet=1; ijet<=2; ijet++)
    {
      TString jetctr(""); jetctr += ijet;
      controlHistos.addHistogram( new TH1F ("jet"+jetctr, "; Jet #"+jetctr+" transverse momentum [GeV]; Events", 100, 0.,1000.) );
      controlHistos.addHistogram( new TH1F ("jet"+jetctr+"eta", "; Jet #"+jetctr+" pseudo-rapidity; Events", 30, 0.,3.) );
    }
  
  //tag and probe analysis for Sec Vtx
  TString jetFlavors[]={"","b","udsg","c","unmatched"};
  size_t nJetFlavors=sizeof(jetFlavors)/sizeof(TString);
  const Double_t ptBins[]={30,35,40,45,50,55,60,65,70,80,90,100,120,140,160,180,200,250,350,400,500,750,1000};
  Int_t nPtbins=sizeof(ptBins)/sizeof(Double_t)-1;
  TString svxAlgos[]={"svx","ivf"};
  for(size_t ialgo=0; ialgo<2; ialgo++)
    {
      for(size_t iflav=0; iflav<nJetFlavors; iflav++)
	{
	  for(size_t itag=0; itag<4; itag++)
	    {
	      TString pf("");
	      if(itag==1) pf="Ltag";
	      if(itag==2) pf="Mtag";
	      if(itag==3) pf="Ttag";	  
	      if(ialgo==0 && iflav==0) {
		controlHistos.addHistogram( new TH1F ("deltapt"+pf, ";Relative balance; Events", 100, 0.,1.) );
		controlHistos.addHistogram( new TH1F("balance"+pf,";p_{T}^{1}/p_{T}^{2};Events",100,0.,1.) );
		controlHistos.addHistogram( new TH1F("dphijj"+pf,";#Delta#phi(j^{(1)},j^{(2)});Events",100,-3.2,3.2) );
	      }
	      controlHistos.addHistogram( new TH2F (jetFlavors[iflav]+"recoil"+svxAlgos[ialgo]+"mass"+pf,   "; SecVtx Mass [GeV]; Jets / (0.2 GeV)",                50, 0.,10., nPtbins,ptBins) );
	      controlHistos.addHistogram( new TH2F (jetFlavors[iflav]+"recoil"+svxAlgos[ialgo]+"lxy"+pf,    "; SecVtx L_{xy} [cm]; Jets / (0.1 cm)",                100, 0.,10., nPtbins,ptBins) );	 
	      controlHistos.addHistogram( new TH2F (jetFlavors[iflav]+"recoil"+svxAlgos[ialgo]+"dr"+pf,     "; #Delta R(jet,SecVtx L_{xy}); Jets / (0.1 cm)",       50, 0.,1.0, nPtbins,ptBins) );	 
	      controlHistos.addHistogram( new TH2F (jetFlavors[iflav]+"recoil"+svxAlgos[ialgo]+"ptfrac"+pf, "; p_{T}(SecVtx L_{xy})/p_{T}(jet); Jets / (0.1 cm)",   50, 0.,2.0, nPtbins,ptBins) );	 
	      if(iflav==0){
		controlHistos.addHistogram(  new TProfile("recoil"+svxAlgos[ialgo]+"lxyvsphi"+pf,"#phi [rad]",100,0,3.2) );
		controlHistos.addHistogram(  new TProfile("recoil"+svxAlgos[ialgo]+"lxyvseta"+pf,"#eta",100,0,2.5) );
		controlHistos.addHistogram(  new TProfile("recoil"+svxAlgos[ialgo]+"lxyvspt"+pf,"p_{T} [GeV]",nPtbins,ptBins) );
		controlHistos.addHistogram(  new TProfile("recoil"+svxAlgos[ialgo]+"lxyvsrun"+pf, "Run number",19000,190000,209000) ) ;
		controlHistos.addHistogram(  new TProfile("recoil"+svxAlgos[ialgo]+"lxyvsnvtx"+pf, "Vertices",15,0,30) );
		
		
		for(size_t ireg=0; ireg<5; ireg++)
		  {
		    TString reg("inc");
		    if(ireg==1) reg="0to0p9";
		    if(ireg==2) reg="0p9to1p1";
		    if(ireg==3) reg="1p1to1p5";
		    if(ireg==4) reg="1p5to2p5";
		    controlHistos.addHistogram( new TH1F (reg+"recoil"+svxAlgos[ialgo]+"lxy"+pf, ";SecVtx L_{xy} [cm]; Jets", 100, 0.,10.) );
		  }
	      }
	    }
	}
    }

  controlHistos.addHistogram(  new TH1F ("leadmm"   , ";Di-muon invariant mass [GeV]; Events;"    ,100,0,10) ) ; 
  controlHistos.addHistogram(  new TH1F ("trailermm", ";Di-muon invariant mass [GeV]; Events;"    ,100,0,10) ) ; 
  
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
      int triggerIdx(-1);
      for(int itrig=0; itrig<=4; itrig++)
	{
	  if(!ev.t_bits[itrig]) continue;
	  hasTrigger=true;
	  triggerPrescale=ev.t_prescale[itrig];
	  triggerIdx=itrig;
	  break;
	}
      if(!hasTrigger) continue;

      //weight for the event
      float weight(weightNom*triggerPrescale);

      //select the jets
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


      std::vector<TString> catsToFill;
      catsToFill.push_back("all");

      controlHistos.fillHisto("triggeridx",catsToFill,triggerIdx,weight);
      controlHistos.fillHisto("njets",catsToFill,selJets.size(),weight);

      //require 2 jets in the event
      if(selJets.size()!=2) continue;

      controlHistos.fillHisto("pthat",catsToFill,ev.pthat,weight);
      controlHistos.fillHisto("nvertices_unwgt",catsToFill,ev.nvtx,weight/weightNom);
      controlHistos.fillHisto("nvertices",catsToFill,ev.nvtx,weight);



      //dijet analysis: check for different b-tagging working points applied on the leading pT jet
      int nSoftMu(false);
      size_t pfstart=jets[0].get("pfstart");
      size_t pfend=jets[0].get("pfend");
      LorentzVector leadMM(0,0,0,0);
      data::PhysicsObjectCollection_t pf = evSummary.getPhysicsObject(DataEventSummaryHandler::PFCANDIDATES);
      for(size_t ipfn=pfstart; ipfn<=pfend; ipfn++)
	{
	  int id=pf[ipfn].get("id");
	  if(abs(id)!=13) continue;
	  nSoftMu++;
	  leadMM += pf[ipfn];
	}
      bool tagHasCSVL(selJets[0].getVal("csv")>0.405), tagHasCSVM(selJets[0].getVal("csv")>0.783), tagHasCSVT(selJets[0].getVal("csv")>0.920);
      std::vector<TString> tagBtags;
      tagBtags.push_back("");
      if(tagHasCSVL) tagBtags.push_back("Ltag");
      if(tagHasCSVM) tagBtags.push_back("Mtag");
      if(tagHasCSVT) tagBtags.push_back("Ttag");
      if(nSoftMu==0) continue;

      //azimuthal angle
      float dphijj=deltaPhi(selJets[0].phi(),selJets[1].phi());
      for(size_t ipf=0; ipf<tagBtags.size(); ipf++) controlHistos.fillHisto("dphijj"+tagBtags[ipf],catsToFill,dphijj,weight);
      if(fabs(dphijj)<2.7) continue;
      
      //balancing variables
      float balance=selJets[1].pt()/selJets[0].pt();
      LorentzVector dijet=selJets[0]+selJets[1];
      double deltaPt=dijet.pt()/(selJets[0].pt()+selJets[1].pt());
      for(size_t ipf=0; ipf<tagBtags.size(); ipf++) {
	controlHistos.fillHisto("balance"+tagBtags[ipf],catsToFill,balance,weight);
	controlHistos.fillHisto("deltapt"+tagBtags[ipf],catsToFill,deltaPt,weight);
      }
      if(balance<0.9) continue;

      //selected jet kinematics
      for(size_t ijet=0; ijet<selJets.size(); ijet++)
	{
	  //kinematics
	  float pt=selJets[ijet].pt();
	  float eta=fabs(selJets[ijet].eta());
	  
	  TString jetctr(""); jetctr+=(ijet+1);
	  controlHistos.fillHisto("jet"+jetctr,catsToFill,pt,weight);
	  controlHistos.fillHisto("jet"+jetctr+"eta",catsToFill,fabs(eta),weight);
	  controlHistos.fillHisto("rawjetpt", catsToFill,pt,weight/triggerPrescale);
	  controlHistos.fillHisto("jetpt", catsToFill,pt,weight);
	  controlHistos.fillHisto("jeteta",catsToFill,fabs(eta),weight);
	}

      //probe characteristics
      float recoilPtNorm(TMath::Min(selJets[1].pt(),ptBins[nPtbins-1]));
      int flavId(9999);
      if(isMC){
	const data::PhysicsObject_t &genJet=selJets[1].getObject("genJet");
	if(genJet.pt()>0) flavId=genJet.info.find("id")->second;
      }
      TString jetFlav("unmatched");
      if(fabs(flavId)==5)      jetFlav="b";
      else if(fabs(flavId)==4) jetFlav="c";
      else if(fabs(flavId)<4)  jetFlav="udsg";
      
      //secondary vertex characteristics evaluated for the probe
      for(size_t ialgo=0; ialgo<2; ialgo++)
	{
	  const data::PhysicsObject_t &svx=selJets[1].getObject(svxAlgos[ialgo]);
	  float lxy=svx.vals.find("lxy")->second;
	  if(lxy<=0) continue;
	  
	  for(size_t ipf=0; ipf<tagBtags.size(); ipf++)
	    {
	      TString pf=tagBtags[ipf];

	      controlHistos.fillProfile("recoil"+svxAlgos[ialgo]+"lxyvseta"+pf,   catsToFill, fabs(selJets[1].eta()), lxy, weight);
	      if(fabs(selJets[1].eta())<1.1)
		{
		  controlHistos.fillProfile("recoil"+svxAlgos[ialgo]+"lxyvsphi"+pf,   catsToFill, fabs(selJets[1].phi()), lxy, weight);
		  controlHistos.fillProfile("recoil"+svxAlgos[ialgo]+"lxyvspt"+pf,    catsToFill, fabs(selJets[1].pt()),  lxy, weight);
		  controlHistos.fillProfile("recoil"+svxAlgos[ialgo]+"lxyvsrun"+pf,   catsToFill, ev.run,              lxy, weight);
		  controlHistos.fillProfile("recoil"+svxAlgos[ialgo]+"lxyvsnvtx"+pf,  catsToFill, ev.nvtx,             lxy, weight);

		  controlHistos.fillHisto("recoil"+svxAlgos[ialgo]+"mass"+pf,   catsToFill, svx.mass(),                 recoilPtNorm,weight);
		  controlHistos.fillHisto("recoil"+svxAlgos[ialgo]+"lxy"+pf,    catsToFill, lxy,                    recoilPtNorm,weight);
		  controlHistos.fillHisto("recoil"+svxAlgos[ialgo]+"dr"+pf,     catsToFill, deltaR(selJets[1],svx),   recoilPtNorm,weight);
		  controlHistos.fillHisto("recoil"+svxAlgos[ialgo]+"ptfrac"+pf, catsToFill, svx.pt()/selJets[1].pt(), recoilPtNorm,weight);
		  if(isMC){
		    controlHistos.fillHisto(jetFlav+"recoil"+svxAlgos[ialgo]+"mass"+pf,   catsToFill, svx.mass(),               recoilPtNorm,weight);
		    controlHistos.fillHisto(jetFlav+"recoil"+svxAlgos[ialgo]+"lxy"+pf,    catsToFill, lxy,                      recoilPtNorm,weight);
		    controlHistos.fillHisto(jetFlav+"recoil"+svxAlgos[ialgo]+"dr"+pf,     catsToFill, deltaR(selJets[1],svx),   recoilPtNorm,weight);
		    controlHistos.fillHisto(jetFlav+"recoil"+svxAlgos[ialgo]+"ptfrac"+pf, catsToFill, svx.pt()/selJets[1].pt(), recoilPtNorm,weight);
		  }
		}
	      
	      TString reg("1p5to2p5");
	      if(fabs(selJets[1].eta())<0.9) reg="0to0p9";
	      else if(fabs(selJets[1].eta())<1.1) reg="0p9to1p1";
	      else if(fabs(selJets[1].eta())<1.5) reg="1p1to1p5";
	      controlHistos.fillHisto(reg+"recoil"+svxAlgos[ialgo]+"lxy"+pf, catsToFill, lxy, weight);
	      controlHistos.fillHisto("increcoil"+svxAlgos[ialgo]+"lxy"+pf, catsToFill, lxy, weight);
	    }
	}
      
      //J/Psi
      int ntrailerSoftMu(0);
      pfstart=jets[1].get("pfstart");
      pfend=jets[1].get("pfend");
      LorentzVector trailerMM(0,0,0,0);
      for(size_t ipfn=pfstart; ipfn<=pfend; ipfn++)
	{
	  int id=pf[ipfn].get("id");
	  if(abs(id)!=13) continue;
	  ntrailerSoftMu++;
	  trailerMM += pf[ipfn];
	}

      if(nSoftMu==2) controlHistos.fillHisto("leadmm",catsToFill,leadMM.mass(),weight);
      if(ntrailerSoftMu==2) controlHistos.fillHisto("trailermm",catsToFill,trailerMM.mass(),weight);

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
