#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
//#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"

#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

#include "EgammaAnalysis/ElectronTools/interface/EGammaCutBasedEleId.h"
#include "EGamma/EGammaAnalysisTools/interface/PFIsolationEstimator.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "CMGTools/External/interface/PileupJetIdAlgo.h"


//Tau stuff... somecleaning needed
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/BaseTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauTagInfo.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "RecoTauTag/TauTagTools/interface/PFTauElementsOperators.h"
#include "RecoTauTag/RecoTau/interface/TauDiscriminationProducerBase.h"

#include "DataFormats/PatCandidates/interface/Tau.h"



#include "TH1D.h"

#include "UserCode/llvv_fwk/interface/llvvObjects.h"
#include "UserCode/llvv_fwk/interface/MacroUtils.h"

using namespace std;
using namespace edm;
using namespace reco;

//
class llvvObjectProducers : public edm::EDFilter 
{

	public:
		llvvObjectProducers(const edm::ParameterSet &iConfig);
		virtual bool filter(edm::Event& iEvent, const edm::EventSetup &iSetup) ;
                virtual bool beginRun(edm::Run & iRun, edm::EventSetup const & iSetup);

	private:
		//selection configuration
		edm::ParameterSet analysisCfg_;

		//pf isolation for e/g objects
		PFIsolationEstimator eIsolator03_, eIsolator04_, gIsolator03_, gIsolator04_;

		//pileup jet id
		PileupJetIdAlgo *puJetIdAlgo_,*cutBasedPuJetIdAlgo_;


               //tool to trace prescale changes
               HLTConfigProvider hltConfig_;

                //pfParticleFlow
                int keepPfCandidates;

                //  keep all GEN
                bool keepFullGenInfo_;

                double   tauPtCut;
                uint64_t tauIdMask;
                double   jetPtCut;
                double   pfCandPtCut;
                double   pfCandDzMax;
};



using namespace std;


//
llvvObjectProducers::llvvObjectProducers(const edm::ParameterSet &iConfig)
{
	//configure selection
	analysisCfg_ = iConfig; 

	//init e/g PF isolation tools
	eIsolator03_.initializeElectronIsolation(true); eIsolator03_.setConeSize(0.3);
	eIsolator04_.initializeElectronIsolation(true); eIsolator04_.setConeSize(0.4);
	gIsolator03_.initializePhotonIsolation(true);   gIsolator03_.setConeSize(0.3);
	gIsolator04_.initializePhotonIsolation(true);   gIsolator04_.setConeSize(0.4);

	//pileup jet id
	puJetIdAlgo_ = new PileupJetIdAlgo( analysisCfg_.getParameter< std::vector<edm::ParameterSet> >( "pujetidAlgo")[0] );
	cutBasedPuJetIdAlgo_ = new PileupJetIdAlgo( analysisCfg_.getParameter< std::vector<edm::ParameterSet> >( "pujetidAlgo")[1] );

        //pfParticleFlow
        keepPfCandidates =  analysisCfg_.getParameter<int>("keepPfCandidates");

         //  keep all GEN
        keepFullGenInfo_ = analysisCfg_.getParameter<bool>("keepFullGenInfo");

        produces<llvvGenEvent>();
        produces<llvvGenParticleCollection>();
        produces<llvvLeptonCollection>();
        produces<llvvMuonInfoCollection>();
        produces<llvvElectronInfoCollection>();
        produces<llvvTauCollection>();
        produces<llvvPhotonCollection>();
        produces<llvvJetCollection>();
        produces<llvvPFParticleCollection>();
        std::vector<edm::InputTag> metSources=analysisCfg_.getParameter<std::vector<edm::InputTag> >("metSource");
        for(unsigned int i=0;i<metSources.size();i++){
           produces<llvvMet>(metSources[i].label());
        }
        produces<std::vector<std::string>, edm::InRun>("triggerPaths");
        produces<std::vector<bool> >("triggerBits");
        produces<std::vector<int> >("triggerPrescales");
        produces< int >("nvtx");

        tauPtCut    = 15;
        tauIdMask   =  ((uint64_t) 1 << llvvTAUID::decayModeFinding) + ((uint64_t) 1 << llvvTAUID::againstMuonLoose) + ((uint64_t) 1 << llvvTAUID::againstElectronLoose);

        jetPtCut    = 10;
        pfCandPtCut = 0.5;
        pfCandDzMax = 0.3;


}

//Save trigger names once per run (instead of once per event)
bool llvvObjectProducers::beginRun(edm::Run & iRun, edm::EventSetup const & iSetup){
  std::vector<std::string> triggerPaths = analysisCfg_.getParameter<std::vector<std::string> >("triggerPaths");
  std::auto_ptr<std::vector<std::string> > triggerPathNameOut (new std::vector<std::string>());
  std::vector<std::string>& triggerPathName = *triggerPathNameOut;
  for(unsigned int i=0;i<triggerPaths.size();i++){triggerPathName.push_back(triggerPaths[i]);}  
  iRun.put(triggerPathNameOut,"triggerPaths");

  bool changed(true);
  hltConfig_.init(iRun, iSetup,"HLT", changed);

  return true;
}

//
bool llvvObjectProducers::filter(edm::Event& iEvent, const edm::EventSetup &iSetup) 
{
  //////////////////////////////////   //////////////////////////////////   //////////////////////////////////   //////////////////////////////////
  //Prepare objects to be saved in the EDM file
  std::auto_ptr<llvvGenEvent> genEvOut(new llvvGenEvent());
  llvvGenEvent& genEv = *genEvOut;

  std::auto_ptr<llvvGenParticleCollection> genPartCollOut(new llvvGenParticleCollection());
  llvvGenParticleCollection& genPartColl = *genPartCollOut;

  std::auto_ptr<llvvLeptonCollection> lepCollOut(new llvvLeptonCollection());
  llvvLeptonCollection& lepColl = *lepCollOut;

  std::auto_ptr<llvvMuonInfoCollection> muInfoCollOut(new llvvMuonInfoCollection());
  llvvMuonInfoCollection& muInfoColl = *muInfoCollOut;

  std::auto_ptr<llvvElectronInfoCollection> elInfoCollOut(new llvvElectronInfoCollection());
  llvvElectronInfoCollection& elInfoColl = *elInfoCollOut;

  std::auto_ptr<llvvTauCollection> tauCollOut(new llvvTauCollection());
  llvvTauCollection& tauColl = *tauCollOut;

  std::auto_ptr<llvvPhotonCollection> phoCollOut(new llvvPhotonCollection());
  llvvPhotonCollection& phoColl = *phoCollOut;

  std::auto_ptr<llvvJetCollection> jetCollOut(new llvvJetCollection());
  llvvJetCollection& jetColl = *jetCollOut;

  std::auto_ptr<llvvPFParticleCollection> pfCollOut(new llvvPFParticleCollection());
  llvvPFParticleCollection& pfColl = *pfCollOut;

  std::auto_ptr<std::vector<bool> > triggerBitsOut(new std::vector<bool>());
  std::vector<bool>& triggerBits = *triggerBitsOut;

  std::auto_ptr<std::vector<int> > triggerPrescalesOut(new std::vector<int>());
  std::vector<int>& triggerPrescales = *triggerPrescalesOut;

  std::auto_ptr<int > nvtxOut(new int());
  int& nvtx = *nvtxOut;

  //////////////////////////////////   //////////////////////////////////   //////////////////////////////////   //////////////////////////////////

  //
  // Global objects
  //
  ESHandle<MagneticField> B;
  iSetup.get<IdealMagneticFieldRecord > ().get(B);
  const MagneticField* magField = B.product();

  ESHandle<GlobalTrackingGeometry> geomHandle;
  iSetup.get<GlobalTrackingGeometryRecord > ().get(geomHandle);


  //
  // Gen Info
  //
  bool isData=iEvent.isRealData();
  int nHardProcGenLeptons(0),nHardProcGenBosons(0);
  if(!isData){

    //pileup
    edm::Handle<std::vector<PileupSummaryInfo> > puInfoH;
    iEvent.getByType(puInfoH);
    genEv.ngenITpu    = 0;
    genEv.ngenOOTpu   = 0;
    genEv.ngenOOTpum1 = 0;
    genEv.ngenTruepu  = 0;
    if(puInfoH.isValid()){
        for(std::vector<PileupSummaryInfo>::const_iterator it = puInfoH->begin(); it != puInfoH->end(); it++){
            if(it->getBunchCrossing()==0)      { genEv.ngenITpu += it->getPU_NumInteractions();   genEv.ngenTruepu  = it->getTrueNumInteractions(); }
            else if(it->getBunchCrossing()<0)  { genEv.ngenOOTpum1 += it->getPU_NumInteractions(); }
            else                               { genEv.ngenOOTpu   += it->getPU_NumInteractions(); }
          }
      }

    //pdf info
    edm::Handle<GenEventInfoProduct> genEventInfoProd;
    iEvent.getByType( genEventInfoProd );
    if(genEventInfoProd.isValid())
      {
        genEv.genWeight = genEventInfoProd->weight();
        genEv.qscale = genEventInfoProd->qScale();
        if(genEventInfoProd->pdf()){
            genEv.qscale = genEventInfoProd->pdf()->scalePDF;
            genEv.x1  = genEventInfoProd->pdf()->x.first;
            genEv.x2  = genEventInfoProd->pdf()->x.second;
            genEv.id1 = genEventInfoProd->pdf()->id.first;
            genEv.id2 = genEventInfoProd->pdf()->id.second;
          }
        if(genEventInfoProd->binningValues().size()>0) genEv.pthat = genEventInfoProd->binningValues()[0];
      }

    //matrix element info
    Handle<LHEEventProduct> lheH;
    iEvent.getByType(lheH);
    if(lheH.isValid()) genEv.nup=lheH->hepeup().NUP;


     Handle<View<Candidate> > genParticlesH;
     iEvent.getByLabel(analysisCfg_.getParameter<edm::InputTag>("genSource"), genParticlesH);

     //analyze first hard process
     bool isSherpa(false);
     for(size_t i = 0; i < genParticlesH->size(); ++ i){
        const reco::GenParticle & p = dynamic_cast<const reco::GenParticle &>( (*genParticlesH)[i] );
        if(abs(p.pdgId())==2212 && p.status()==4) isSherpa=true; //this is only used by sherpa
        bool isHardProc(p.status()==3);
        bool isStableOfInterest(keepFullGenInfo_ && p.status()==1 && ((abs(p.pdgId())==22 && p.pt()>20) || (p.charge()!=0 && p.pt()>0.5 && fabs(p.eta())<3.0 )) );
        if(!isHardProc && !isStableOfInterest) continue;

        //check if lepton is ok to accept (for unfolding purposes)
        bool passLepAccCut( p.charge()!=0 &&  p.pt()>0.5 && fabs(p.eta())<3.0 );
      bool leptonIsOkToAccept(isHardProc && abs(p.pdgId())>=11 && abs(p.pdgId())<=14);
      if( leptonIsOkToAccept )
        {
          if( !keepFullGenInfo_ ) leptonIsOkToAccept=true;
          else
            {
           leptonIsOkToAccept &= passLepAccCut;
      //madgraph and pythia-based
           if(  p.numberOfMothers()  == 1 &&  p.mother()->pdgId() != p.pdgId() && p.mother()->pdgId() != 23 && abs (p.mother()->pdgId() ) != 24 ) leptonIsOkToAccept=false;
           
           //sherpa-like
           if(  p.numberOfMothers()  == 2 && (abs(p.mother(0)->pdgId()) != abs(p.pdgId()) ||  abs(p.mother(1)->pdgId()) != abs(p.pdgId())) )  leptonIsOkToAccept=false;
            }
        }
      nHardProcGenLeptons += leptonIsOkToAccept;
      nHardProcGenBosons  += (isHardProc && (abs(p.pdgId())==24 || abs(p.pdgId())==23));
              
        llvvGenParticle part;
        part.SetPxPyPzE(p.px(), p.py(), p.pz(), p.energy());
        part.id     =p.pdgId();
        part.status =p.status();
        part.lxy    =0;
           
        //check if photon is prompt or radiated from quark/line
        if(fabs(p.pdgId())==22){
           bool isPrompt(false);
           for(size_t b = 0; b < p.numberOfMothers(); ++ b){
              const reco::Candidate *p_m = p.mother(b);
              if(abs(p_m->pdgId())>25) continue;
              isPrompt=true;
           }
          if(!isPrompt) part.lxy=99999.;
        } 

     genPartColl.push_back(part);
     }
       
       // FSR photons (if full gen info is set to true)
       if(keepFullGenInfo_)
         {
           int NGenPart = genPartColl.size() ;
           for(int j = 0; j < NGenPart; j++ )
             {
               if ( fabs(genPartColl[j].status) != 1 && fabs(genPartColl[j].status) != 3 ) continue; 
               if(fabs(genPartColl[j].id) != 11 && fabs(genPartColl[j].id) != 13 ) continue;
               for(size_t i = 0; i < genParticlesH->size(); ++ i)
                 {
                   const reco::GenParticle & p = dynamic_cast<const reco::GenParticle &>( (*genParticlesH)[i] );
                   if (!(abs(p.pdgId()) == 22 && p.pt() <= 20 &&  p.pt() > 1e-6 ) ) continue ;
                   if( deltaR( genPartColl[j].eta(), genPartColl[j].phi(), p.eta(), p.phi()) > 0.15) continue;
           
                    llvvGenParticle part;
                    part.SetPxPyPzE(p.px(), p.py(), p.pz(), p.energy());
                    part.id     =p.pdgId();
                    part.status =p.status();
                    part.lxy    =0;
                    genPartColl.push_back(part);
                 }
             }
        }
       
       //heavy flavors
       for (size_t i=0; i<genParticlesH->size(); i++)
         {
           const reco::GenParticle & p = dynamic_cast<const reco::GenParticle &>( (*genParticlesH)[i] );
           if(p.status()!=3 || abs(p.pdgId())!=5) continue;

           const reco::Candidate *fs = utils::cmssw::getGeneratorFinalStateFor(&p,isSherpa);
           if(fs->numberOfDaughters()==0) continue;
           fs = utils::cmssw::getGeneratorFinalStateFor( fs->daughter(0), isSherpa );
           for(size_t j=0; j<fs->numberOfDaughters(); j++)
             {
               const reco::Candidate *d=fs->daughter(j);
               if(d==0) continue;
               int absid=abs(d->pdgId());
               if(!utils::cmssw::isBhadron(absid)) continue;
               
               //find first stable particle to trace decay length
               const reco::Candidate *gd=d;
               while(1)
                 {
                   if(gd->status()==1) break;
                   if(gd==0 || gd->numberOfDaughters()==0) break;
                   const reco::Candidate *newGd=gd; 
                   for(size_t k=0; k<gd->numberOfDaughters(); k++)
                     {
                       if(gd->daughter(k)->pdgId()==22) continue;
                       if(gd->daughter(k)==gd) continue;
                       newGd=gd->daughter(k);
                     }
                   if(gd==newGd) break;
                   gd=newGd;
                 }

               float vxMother(p.vx()),     vyMother(p.vy());
               float vxDaughter(gd->vx()), vyDaughter(gd->vy());

               //save B-hadron information
               llvvGenParticle part;
               part.SetPxPyPzE(d->px(), d->py(), d->pz(), d->energy());
               part.id     =d->pdgId();
               part.status =d->status();
               part.lxy    =sqrt(pow(vxMother-vxDaughter,2)+pow(vyMother-vyDaughter,2));
               genPartColl.push_back(part);
             }
         }
   }












  //////////////////////////////////   //////////////////////////////////   //////////////////////////////////   //////////////////////////////////
  // process triggers
        edm::InputTag trigSource              = analysisCfg_.getParameter<edm::InputTag>("triggerSource");
        std::vector<std::string> triggerPaths = analysisCfg_.getParameter<std::vector<std::string> >("triggerPaths");
        std::vector<int> triggerCats          = analysisCfg_.getParameter<std::vector<int> >("triggerCats");

        //ev.tn = triggerPaths.size();  
        edm::Handle<edm::TriggerResults> triggerBitsH;
        iEvent.getByLabel( trigSource, triggerBitsH);
        const edm::TriggerNames &triggerNames = iEvent.triggerNames( *triggerBitsH );

        for(size_t it=0; it<triggerPaths.size(); it++){
           triggerBits.push_back(false);
           triggerPrescales.push_back(0);
        }
        for (size_t itrig = 0; itrig != triggerBitsH->size(); ++itrig)
        {
                if( !triggerBitsH->wasrun(itrig) || triggerBitsH->error(itrig) || !triggerBitsH->accept(itrig) )continue;
                std::string trigName = triggerNames.triggerName(itrig);
                for(size_t it=0; it<triggerPaths.size(); it++)
                {
                        if(trigName.find(triggerPaths[it]) == std::string::npos) continue;
                        triggerBits[it] = true;
                        triggerPrescales[it] = hltConfig_.prescaleValue(iEvent, iSetup, trigName);
                        break;
                }
        }

  //////////////////////////////////   //////////////////////////////////   //////////////////////////////////   //////////////////////////////////
  //process vertex and beam spot
	edm::Handle<reco::BeamSpot> beamSpotH;
	iEvent.getByLabel( analysisCfg_.getParameter<edm::InputTag>("beamSpotSource"), beamSpotH);
	edm::Handle<reco::VertexCollection> vtxH;
	iEvent.getByLabel( analysisCfg_.getParameter<edm::InputTag>("vtxSource"), vtxH);
	if(vtxH->size()==0) return false;
        nvtx = vtxH->size();
	reco::VertexRef primVtx(vtxH,0);

	edm::Handle< double > rho, rho25;
	iEvent.getByLabel( analysisCfg_.getParameter<edm::InputTag>("rhoSource"),rho);
	iEvent.getByLabel( analysisCfg_.getParameter<edm::InputTag>("rho25Source"), rho25);


  //////////////////////////////////   //////////////////////////////////   //////////////////////////////////   //////////////////////////////////
	//
	// charged leptons
	// https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId
	// https://twiki.cern.ch/twiki/bin/view/CMS/EgammaCutBasedIdentification
	// https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPFBasedIsolation
	//
        int nElecs=0, nMuons=0, nPhotons=0;
  //////////////////////////////////   //////////////////////////////////   //////////////////////////////////   //////////////////////////////////
  //Muons
	edm::Handle<View<Candidate> > muH, eH;
	iEvent.getByLabel( analysisCfg_.getParameter<edm::InputTag>("muonSource"),     muH);
	iEvent.getByLabel( analysisCfg_.getParameter<edm::InputTag>("electronSource"), eH);

	edm::Handle<reco::ConversionCollection> convH;
	iEvent.getByLabel( analysisCfg_.getParameter<edm::InputTag>("conversionSource"), convH);
	EcalClusterLazyTools egLazyTool( iEvent, iSetup, analysisCfg_.getParameter<edm::InputTag>("ebrechitsSource"), analysisCfg_.getParameter<edm::InputTag>("eerechitsSource") );

	//particle flow candidates
	edm::Handle<reco::PFCandidateCollection> pfH;
	iEvent.getByLabel( analysisCfg_.getParameter<edm::InputTag>("pfSource"),pfH);

	for(size_t imu=0; imu< muH->size(); ++imu)
	{
		reco::CandidatePtr muonPtr    = muH->ptrAt(imu);
		const pat::Muon *muon         = dynamic_cast<const pat::Muon *>( muonPtr.get() );
		const reco::Candidate *genLep = muon->genLepton();

		//apply a loose pre-selection to our muon candidates 
		bool isPF( muon->isPFMuon() );
		bool isGlobal( muon->isGlobalMuon() );
		bool isTracker( muon->isTrackerMuon() );
		bool isLoose( isPF && (isGlobal || isTracker) );
		if(!isLoose) continue;
		if(muon->pt()<3 || fabs(muon->eta())>2.5) continue;

                llvvLepton lep;
                llvvMuonInfo muInfo;

		//store information
                lep.SetPxPyPzE(muon->px(), muon->py(), muon->pz(), muon->energy()); 
		lep.id                         = -13*muon->charge();
//		lep.pid                        = ev.mn;
		lep.isPF                       = isPF;                
		lep.genid                      = genLep ? genLep->pdgId() :0;
                if(genLep)lep.gen.SetPxPyPzE(genLep->px(), genLep->py(), genLep->pz(), genLep->energy());
                if(isGlobal){lep.trk.SetPxPyPzE(muon->globalTrack()->px(), muon->globalTrack()->py(), muon->globalTrack()->pz(), muon->globalTrack()->p());
                }else{       lep.trk.SetPxPyPzE(muon->innerTrack()->px(), muon->innerTrack()->py(), muon->innerTrack()->pz(), muon->innerTrack()->p()); }
//		lep.trkpt                      = isGlobal ? muon->globalTrack()->pt()                                    : muon->innerTrack()->pt();
//		lep.trketa                     = isGlobal ? muon->globalTrack()->eta()                                   : muon->innerTrack()->eta();
//		lep.trkphi                     = isGlobal ? muon->globalTrack()->phi()                                   : muon->innerTrack()->phi();
		lep.trkchi2                    = isGlobal ? muon->globalTrack()->normalizedChi2()                        : muon->innerTrack()->normalizedChi2();	  
		lep.trkValidPixelHits          = isGlobal ? muon->globalTrack()->hitPattern().numberOfValidPixelHits()   : muon->innerTrack()->hitPattern().numberOfValidPixelHits();
		lep.trkValidTrackerHits        = isGlobal ? muon->globalTrack()->hitPattern().numberOfValidTrackerHits() : muon->innerTrack()->hitPattern().numberOfValidTrackerHits();
		lep.trkLostInnerHits           = muon->innerTrack()->trackerExpectedHitsInner().numberOfLostHits();
		lep.trkPtErr                   = fabs(muon->innerTrack()->ptError()/muon->innerTrack()->pt());
		lep.d0                         = fabs(muon->muonBestTrack()->dxy(primVtx->position()));
		lep.dZ                         = fabs(muon->muonBestTrack()->dz(primVtx->position()));
		std::pair<bool,Measurement1D> ip3dRes   = utils::cmssw::getImpactParameter<reco::TrackRef>(muon->innerTrack(), primVtx, iSetup, true);
		lep.ip3d                       = ip3dRes.second.value();
		lep.ip3dsig                    = ip3dRes.second.significance();
		lep.ecalIso03                  = muon->isolationR03().emEt;
		lep.hcalIso03                  = muon->isolationR03().hadEt;
		lep.trkIso03                   = muon->isolationR03().sumPt;
		lep.ecalIso04                  = muon->isolationR05().emEt; //oh well....
		lep.hcalIso04                  = muon->isolationR05().hadEt;
		lep.trkIso04                   = muon->isolationR05().sumPt;
		lep.gIso03                     = muon->pfIsolationR03().sumPhotonEt;
		lep.chIso03                    = muon->pfIsolationR03().sumChargedHadronPt;
		lep.puchIso03                  = muon->pfIsolationR03().sumPUPt;
		lep.nhIso03                    = muon->pfIsolationR03().sumNeutralHadronEt;
		lep.gIso04                     = muon->pfIsolationR04().sumPhotonEt;
		lep.chIso04                    = muon->pfIsolationR04().sumChargedHadronPt;
		lep.puchIso04                  = muon->pfIsolationR04().sumPUPt;
		lep.nhIso04                    = muon->pfIsolationR04().sumNeutralHadronEt;
                //muon specific info
		muInfo.nMatches                   = muon->numberOfMatches();
		muInfo.nMatchedStations           = muon->numberOfMatchedStations();
		muInfo.validMuonHits              = isGlobal ? muon->globalTrack()->hitPattern().numberOfValidMuonHits() : 0.;
		muInfo.innerTrackChi2             = isTracker ? muon->innerTrack()->normalizedChi2() : 0.;
		muInfo.trkLayersWithMeasurement   = muon->track()->hitPattern().trackerLayersWithMeasurement();
		muInfo.pixelLayersWithMeasurement = isTracker ? muon->innerTrack()->hitPattern().pixelLayersWithMeasurement() : 0.;

		bool isTight( isPF                                 && isGlobal              
				&& fabs(lep.d0)<0.2 	 && fabs(lep.dZ)<0.5         && lep.trkValidPixelHits>0         && lep.trkchi2<10. 
				&& muInfo.validMuonHits>0.     && muInfo.nMatchedStations>1   && muInfo.trkLayersWithMeasurement>5 );
		bool isSoft(isTracker && muon->muonID("TMOneStationTight") 
				&& fabs(lep.d0)<3.  && fabs(lep.dZ)<30.
				&& muInfo.trkLayersWithMeasurement>5 && muInfo.pixelLayersWithMeasurement>1  && muInfo.innerTrackChi2 < 1.8 );
		bool isHighNew = true;//LoicQ Commented out muon::isHighPtMuon(dynamic_cast<const reco::Muon &>(*muon), dynamic_cast<const reco::Vertex &> (*primVtx)) ;

		//save id summary
		lep.idbits                     = 
			( (int(muon->muonID("GlobalMuonPromptTight")) & 0x1)   << 0)
			| ( (int(muon->muonID("TMLastStationLoose")) & 0x1)    << 1)
			| ( (int(muon->muonID("TMLastStationTight")) & 0x1)    << 2)
			| ( (int(muon->muonID("TMLastStationAngTight")) & 0x1) << 3)
			| ( (int(muon->muonID("TMOneStationTight")) & 0x1)     << 4)
			| ( isTracker                                          << 5)
			| ( isGlobal                                           << 6)
			| ( isPF                                               << 7)
			| ( isLoose                                            << 8)
			| ( isSoft                                             << 9)
			| ( isTight                                            << 10)
			| ( isHighNew                                          << 11);

		//add trigger match
		int TrigSum(0);
		for(size_t it=0; it<triggerPaths.size(); it++)
		{
			string tempTrigName = triggerPaths[it] + "*";
			if ( muon->triggerObjectMatchesByPath(tempTrigName).size() > 0 ) TrigSum |= (1<<it); 
		}
		lep.Tbits = TrigSum ;

               muInfoColl.push_back(muInfo);
               lepColl.push_back(lep);
               if(muon->pt()>18) nMuons++;
	}

  //////////////////////////////////   //////////////////////////////////   //////////////////////////////////   //////////////////////////////////
  //Electrons

	for(size_t iele=0; iele< eH->size(); ++iele)
	{
		reco::CandidatePtr elePtr       = eH->ptrAt(iele);
		const pat::Electron *ele        = dynamic_cast<const pat::Electron *>( elePtr.get() );
		const reco::Candidate *genLep   = ele->genLepton();
		const reco::GsfElectron *gsfEle = dynamic_cast<const reco::GsfElectron *>(ele);

		//pre-selection
		if(ele->gsfTrack().isNull() || ele->superCluster().isNull() || gsfEle==0) continue;
		if(ele->pt()<10 || !(ele->isEB() || ele->isEE()) )                        continue;
		bool overlapFound(false);
		for(unsigned int ilep=0; ilep<lepColl.size(); ilep++)
		{
			if( deltaR( lepColl[ilep].eta(), lepColl[ilep].phi(), ele->eta(), ele->phi()) > 0.1) continue;
			overlapFound=true;
			break;
		}
		if(overlapFound) continue;


                llvvLepton lep;
                llvvElectronInfo elInfo;

		//store information
                lep.SetPxPyPzE(ele->px(), ele->py(), ele->pz(),  ele->energy());
		lep.id                         = -11*ele->charge();
//		lep.pid                        = ev.egn;
		lep.isPF                       = ele->isPF(); 
		lep.genid                      = genLep ? genLep->pdgId() :0;
                if(genLep)lep.gen.SetPxPyPzE(genLep->px(), genLep->py(), genLep->pz(), genLep->energy());
//		lep.genpx                      = genLep ? genLep->px() : 0;
//		lep.genpy                      = genLep ? genLep->py() : 0;
//		lep.genpz                      = genLep ? genLep->pz() : 0;
//		lep.genen                      = genLep ? genLep->energy() : 0;
                lep.trk.SetPxPyPzE(ele->gsfTrack()->px(), ele->gsfTrack()->py(), ele->gsfTrack()->pz(), ele->gsfTrack()->p());
//		lep.trkpt                      = ele->gsfTrack()->pt();
//		lep.trketa                     = ele->gsfTrack()->eta();
//		lep.trkphi                     = ele->gsfTrack()->phi();
		lep.trkchi2                    = ele->gsfTrack()->normalizedChi2();
		lep.trkValidPixelHits          = ele->gsfTrack()->hitPattern().numberOfValidPixelHits();
		lep.trkValidTrackerHits        = ele->gsfTrack()->hitPattern().numberOfValidTrackerHits();
		lep.trkLostInnerHits           = ele->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits();
		lep.trkPtErr                   = fabs(ele->gsfTrack()->ptError()/ele->gsfTrack()->pt());
		lep.d0                         = fabs(ele->gsfTrack()->dxy(primVtx->position()));
		lep.dZ                         = fabs(ele->gsfTrack()->dz(primVtx->position()));
		std::pair<bool,Measurement1D> ip3dRes = utils::cmssw::getImpactParameter<reco::GsfTrackRef>(ele->gsfTrack(), primVtx, iSetup, true);
		lep.ip3d                       = ip3dRes.second.value();
		lep.ip3dsig                    = ip3dRes.second.significance();
		lep.ecalIso03                  = ele->dr03EcalRecHitSumEt();
		lep.hcalIso03                  = ele->dr03HcalTowerSumEt();
		lep.trkIso03                   = ele->dr03TkSumPt();
		lep.ecalIso04                  = ele->dr04EcalRecHitSumEt();
		lep.hcalIso04                  = ele->dr04HcalTowerSumEt();
		lep.trkIso04                   = ele->dr04TkSumPt();
		eIsolator03_.fGetIsolation(gsfEle, &(*pfH), primVtx, vtxH);
		lep.gIso03                     = eIsolator03_.getIsolationPhoton();
		lep.chIso03                    = eIsolator03_.getIsolationCharged();
		lep.puchIso03                  = 0;
		lep.nhIso03                    = eIsolator03_.getIsolationNeutral();
		eIsolator04_.fGetIsolation(gsfEle, &(*pfH), primVtx, vtxH);
		lep.gIso04                     = eIsolator04_.getIsolationPhoton();
		lep.chIso04                    = eIsolator04_.getIsolationCharged();
		lep.puchIso04                  = 0;
		lep.nhIso04                    = eIsolator04_.getIsolationNeutral();


		elInfo.isConv                   = ConversionTools::hasMatchedConversion(*gsfEle,convH,beamSpotH->position());
		elInfo.eopin                    = ele->eSuperClusterOverP(); 
		elInfo.eopout                   = ele->eEleClusterOverPout();
		elInfo.sce                      = ele->superCluster()->energy();
		elInfo.sceta                    = ele->superCluster()->eta();
		elInfo.scphi                    = ele->superCluster()->phi();
		elInfo.fbrem                    = ele->fbrem();
		elInfo.sihih                    = ele->sigmaIetaIeta();
		vector<float> cov                       = egLazyTool.localCovariances(*ele->superCluster()->seed()); 
		elInfo.sipip                    = sqrt(cov[2]); 
		elInfo.sihip                    = cov[1];
		elInfo.r9                       = egLazyTool.e3x3(*ele->superCluster()->seed())/ele->superCluster()->rawEnergy();
		elInfo.mvatrigv0                = ele->electronID("mvaTrigV0");
		elInfo.mvanontrigv0             = ele->electronID("mvaNonTrigV0");
		elInfo.hoe                      = ele->hadronicOverEm();
		elInfo.dphiin                   = ele->deltaPhiSuperClusterTrackAtVtx();
		elInfo.detain                   = ele->deltaEtaSuperClusterTrackAtVtx();
		elInfo.ooemoop                  = (1.0/ele->ecalEnergy() - ele->eSuperClusterOverP()/ele->ecalEnergy());

		//save id summary
		bool hasVetoId(false);
		lep.idbits = 
			(ele->ecalDrivenSeed()                                                                                  << 0)
			| ( ele->trackerDrivenSeed()                                                                            << 1)
			| ( EgammaCutBasedEleId::PassEoverPCuts(elInfo.sceta,elInfo.eopin,elInfo.fbrem) << 2);

		for(size_t iid=0; iid<4; iid++)
		{
			int id(EgammaCutBasedEleId::VETO);
			if(iid==1) id=EgammaCutBasedEleId::LOOSE;
			if(iid==2) id=EgammaCutBasedEleId::MEDIUM;
			if(iid==3) id=EgammaCutBasedEleId::TIGHT;
			bool hasId=EgammaCutBasedEleId::PassWP( EgammaCutBasedEleId::WorkingPoint(id), ele->isEB(), ele->pt(), ele->eta(), 
					elInfo.detain, elInfo.dphiin, elInfo.sihih, elInfo.hoe, elInfo.ooemoop, 
					lep.d0, lep.dZ, 
					0., 0., 0., elInfo.isConv, lep.trkLostInnerHits, *rho);
			if(iid==0) hasVetoId=hasId;
			lep.idbits |=  (hasId << (3+iid));
		}
		for(size_t iid=0; iid<2; iid++)
		{
			int id(EgammaCutBasedEleId::TRIGGERTIGHT);
			if(iid==1) id=EgammaCutBasedEleId::TRIGGERWP70;
			bool hasId = EgammaCutBasedEleId::PassTriggerCuts(EgammaCutBasedEleId::TriggerWorkingPoint(id),
					ele->isEB(), ele->pt(),
					elInfo.detain, elInfo.dphiin, elInfo.sihih, elInfo.hoe,
					lep.trkIso03, lep.ecalIso03, lep.hcalIso03);
			lep.idbits |= (hasId << (7+iid));
		}

		// add heep selector , they have quite a few 
		bool boolHeep( ele -> userInt("HEEPId") < 1 ) ;
		lep.idbits |= ( boolHeep <<  9 ) ;

		//add trigger match
		int TrigSum(0);
		for(size_t it=0; it<triggerPaths.size(); it++)
		{
			string tempTrigName = triggerPaths[it] + "*";
			if ( ele->triggerObjectMatchesByPath(tempTrigName).size() > 0 ) TrigSum |= (1<<it); 
		}
		lep.Tbits = TrigSum ;
		//require a very loose baseline id
		if(!hasVetoId) continue;

               elInfoColl.push_back(elInfo);
               lepColl.push_back(lep);
               if(ele->pt()>18) nElecs++;
	}


  //////////////////////////////////   //////////////////////////////////   //////////////////////////////////   //////////////////////////////////
  //Taus
        if(analysisCfg_.getParameter<edm::InputTag>("tauSource").label()!=""){
		Handle<pat::TauCollection> tauH;
		iEvent.getByLabel(analysisCfg_.getParameter<edm::InputTag>("tauSource"), tauH);
		for(size_t itau=0; itau<tauH->size(); itau++)
		{
		   const pat::Tau* tau = &((*tauH)[itau]);
		   const reco::Candidate *genLep   = tau->genLepton();

		   llvvTau tauInfo;

		   tauInfo.SetPxPyPzE(tau->px(), tau->py(), tau->pz(),  tau->energy());
		   tauInfo.id                         = -15*tau->charge();
		   tauInfo.isPF                       = tau->isPFTau();
		   tauInfo.genid                      = genLep ? genLep->pdgId() :0;
		   if(genLep)tauInfo.gen.SetPxPyPzE(genLep->px(), genLep->py(), genLep->pz(), genLep->energy());
		   if(tau->leadPFChargedHadrCand().isNonnull() &&  tau->leadPFChargedHadrCand()->trackRef().isNonnull()){
		   tauInfo.trkchi2                    = tau->leadPFChargedHadrCand()->trackRef()->normalizedChi2();
		   tauInfo.trkValidPixelHits          = tau->leadPFChargedHadrCand()->trackRef()->hitPattern().numberOfValidPixelHits();
		   tauInfo.trkValidTrackerHits        = tau->leadPFChargedHadrCand()->trackRef()->hitPattern().numberOfValidTrackerHits();
		   tauInfo.trkLostInnerHits           = tau->leadPFChargedHadrCand()->trackRef()->trackerExpectedHitsInner().numberOfLostHits();
		   tauInfo.trkPtErr                   = fabs(tau->leadPFChargedHadrCand()->trackRef()->ptError()/tau->leadPFChargedHadrCand()->trackRef()->pt());
		   tauInfo.d0                         = fabs(tau->leadPFChargedHadrCand()->trackRef()->dxy(primVtx->position()));
		   tauInfo.dZ                         = fabs(tau->leadPFChargedHadrCand()->trackRef()->dz(primVtx->position()));
		   std::pair<bool,Measurement1D> ip3dRes = utils::cmssw::getImpactParameter<reco::TrackRef>(tau->leadPFChargedHadrCand()->trackRef(), primVtx, iSetup, true);
		   tauInfo.ip3d                       = ip3dRes.second.value();
		   tauInfo.ip3dsig                    = ip3dRes.second.significance();
		   }
		   tauInfo.vz      = tau->vz();
		   tauInfo.z_expo = 0;
		   if (tau->leadPFChargedHadrCand().isNonnull() && tau->leadPFChargedHadrCand()->trackRef().isNonnull()) {
		     reco::TransientTrack track(tau->leadPFChargedHadrCand()->trackRef(), magField, geomHandle);
		     TransverseImpactPointExtrapolator extrapolator(magField);
		     TrajectoryStateOnSurface closestOnTransversePlaneState = extrapolator.extrapolate(track.impactPointState(), GlobalPoint(beamSpotH->position().x(), beamSpotH->position().y(), 0.0));
		     tauInfo.z_expo = closestOnTransversePlaneState.globalPosition().z();
		   }

		   tauInfo.jet.SetPxPyPzE(tau->pfJetRef().get()->px(), tau->pfJetRef().get()->py(), tau->pfJetRef().get()->pz(), tau->pfJetRef().get()->energy() );
		   tauInfo.numChargedParticlesSigCone   = tau->signalPFChargedHadrCands().size();
		   tauInfo.numNeutralHadronsSigCone     = tau->signalPFNeutrHadrCands().size();
		   tauInfo.numPhotonsSigCone            = tau->signalPFGammaCands().size();
		   tauInfo.numParticlesSigCone          = tau->signalPFCands().size();
		   tauInfo.numPiZeroSigCone             = tau->signalPiZeroCandidates().size();

		   tauInfo.numChargedParticlesIsoCone   = tau->isolationPFChargedHadrCands().size();
		   tauInfo.numNeutralHadronsIsoCone     = tau->isolationPFNeutrHadrCands().size();
		   tauInfo.numPhotonsIsoCone            = tau->isolationPFGammaCands().size();
		   tauInfo.numParticlesIsoCone          = tau->isolationPFCands().size();
		   tauInfo.ptSumChargedParticlesIsoCone = tau->isolationPFChargedHadrCandsPtSum();
		   tauInfo.ptSumPhotonsIsoCone          = tau->isolationPFGammaCandsEtSum();
	 
		   tauInfo.mva_e_pi                     = (tau->leadPFChargedHadrCand().isNonnull() ? tau->leadPFChargedHadrCand()->mva_e_pi() : 0.);
		   tauInfo.mva_pi_mu                    = (tau->leadPFChargedHadrCand().isNonnull() ? tau->leadPFChargedHadrCand()->mva_pi_mu() : 0.);
		   tauInfo.mva_e_mu                     = (tau->leadPFChargedHadrCand().isNonnull() ? tau->leadPFChargedHadrCand()->mva_e_mu() : 0.);
		   tauInfo.hcalEnergy                   = (tau->leadPFChargedHadrCand().isNonnull() ? tau->leadPFChargedHadrCand()->hcalEnergy() : 0.);
		   tauInfo.ecalEnergy                   = (tau->leadPFChargedHadrCand().isNonnull() ? tau->leadPFChargedHadrCand()->ecalEnergy() : 0.);
		   tauInfo.emfraction                   = tau->emFraction();

                   tauInfo.idbits  = 0;
		   tauInfo.idbits += tau->tauID("decayModeFinding"                           )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::decayModeFinding;
		   tauInfo.idbits += tau->tauID("byVLooseCombinedIsolationDeltaBetaCorr"     )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::byVLooseCombinedIsolationDeltaBetaCorr;
		   tauInfo.idbits += tau->tauID("byLooseCombinedIsolationDeltaBetaCorr"      )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::byLooseCombinedIsolationDeltaBetaCorr;
		   tauInfo.idbits += tau->tauID("byMediumCombinedIsolationDeltaBetaCorr"     )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::byMediumCombinedIsolationDeltaBetaCorr;
		   tauInfo.idbits += tau->tauID("byTightCombinedIsolationDeltaBetaCorr"      )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::byTightCombinedIsolationDeltaBetaCorr;
		   tauInfo.idbits += tau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits" )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::byLooseCombinedIsolationDeltaBetaCorr3Hits;
		   tauInfo.idbits += tau->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits")<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::byMediumCombinedIsolationDeltaBetaCorr3Hits;
		   tauInfo.idbits += tau->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits" )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::byTightCombinedIsolationDeltaBetaCorr3Hits;
		   tauInfo.idbits += tau->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"   )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::byCombinedIsolationDeltaBetaCorrRaw3Hits;
		   tauInfo.idbits += tau->tauID("againstElectronLoose"                       )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::againstElectronLoose;                    
		   tauInfo.idbits += tau->tauID("againstElectronMedium"                      )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::againstElectronMedium;
		   tauInfo.idbits += tau->tauID("againstElectronTight"                       )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::againstElectronTight;
		   tauInfo.idbits += tau->tauID("againstElectronMVA3category"                )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::againstElectronMVA3category;
                   tauInfo.idbits += tau->tauID("againstElectronMVA3raw"                     )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::againstElectronMVA3raw;
		   tauInfo.idbits += tau->tauID("againstElectronLooseMVA3"                   )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::againstElectronLooseMVA3;
		   tauInfo.idbits += tau->tauID("againstElectronMediumMVA3"                  )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::againstElectronMediumMVA3;
		   tauInfo.idbits += tau->tauID("againstElectronTightMVA3"                   )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::againstElectronTightMVA3;
		   tauInfo.idbits += tau->tauID("againstElectronVTightMVA3"                  )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::againstElectronVTightMVA3;
		   tauInfo.idbits += tau->tauID("againstMuonLoose"                           )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::againstMuonLoose;
		   tauInfo.idbits += tau->tauID("againstMuonMedium"                          )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::againstMuonMedium;
		   tauInfo.idbits += tau->tauID("againstMuonTight"                           )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::againstMuonTight;
		   tauInfo.idbits += tau->tauID("againstMuonLoose2"                          )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::againstMuonLoose2;
		   tauInfo.idbits += tau->tauID("againstMuonMedium2"                         )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::againstMuonMedium2;
		   tauInfo.idbits += tau->tauID("againstMuonTight2"                          )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::againstMuonTight2;
		   tauInfo.idbits += tau->tauID("againstMuonLoose3"                          )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::againstMuonLoose3;
		   tauInfo.idbits += tau->tauID("againstMuonTight3"                          )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::againstMuonTight3;
		   tauInfo.idbits += tau->tauID("byIsolationMVAraw"                          )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::byIsolationMVAraw;
		   tauInfo.idbits += tau->tauID("byLooseIsolationMVA"                        )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::byLooseIsolationMVA;
		   tauInfo.idbits += tau->tauID("byMediumIsolationMVA"                       )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::byMediumIsolationMVA;
		   tauInfo.idbits += tau->tauID("byTightIsolationMVA"                        )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::byTightIsolationMVA;
		   tauInfo.idbits += tau->tauID("byIsolationMVA2raw"                         )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::byIsolationMVA2raw;
		   tauInfo.idbits += tau->tauID("byLooseIsolationMVA2"                       )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::byLooseIsolationMVA2;
		   tauInfo.idbits += tau->tauID("byMediumIsolationMVA2"                      )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::byMediumIsolationMVA2;
		   tauInfo.idbits += tau->tauID("byTightIsolationMVA2"                       )<=0.5 ? 0 : (uint64_t) 1 << llvvTAUID::byTightIsolationMVA2;


                   //save charged hadron information
       	           for(unsigned int iCharged=0; iCharged < tau->signalPFChargedHadrCands().size() && iCharged<3; iCharged++){
		      const reco::PFCandidateRef& cand = tau->signalPFChargedHadrCands().at(iCharged);
	   	      math::XYZTLorentzVector candP4 = cand->p4();
                      tauInfo.tracks.push_back(LorentzVectorF(candP4.px(), candP4.py(), candP4.pz(), candP4.energy()));
	           }

                   //save neutral hadron information
 	           for(unsigned int iPi0=0; iPi0 < tau->signalPiZeroCandidates().size() && iPi0<2; iPi0++){
		      const reco::RecoTauPiZero& cand = tau->signalPiZeroCandidates().at(iPi0);
		      math::XYZTLorentzVector candP4 = cand.p4();
                      tauInfo.pi0s.push_back(LorentzVectorF(candP4.px(), candP4.py(), candP4.pz(), candP4.energy()));
 	           }

                   if(tauInfo.pt()>tauPtCut && (tauInfo.idbits&tauIdMask)>0 ){
                      tauColl.push_back(tauInfo);
                   }
		}
        }


  //////////////////////////////////   //////////////////////////////////   //////////////////////////////////   //////////////////////////////////
  //Photons
	//
	// photon selection
	// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonID2012
	//
	edm::Handle<edm::View<reco::Candidate> > photonH;
	iEvent.getByLabel( analysisCfg_.getParameter<edm::InputTag>("photonSource"), photonH );
	edm::Handle<reco::GsfElectronCollection> gsfEleH;
	iEvent.getByLabel( "gsfElectrons", gsfEleH );
	for(size_t ipho=0; ipho<photonH->size(); ipho++)
	{
		const reco::Photon *pho = dynamic_cast<const reco::Photon *>( photonH->ptrAt(ipho).get() );
		if(pho==0) continue;
		if(pho->pt()<20 || !(pho->isEB() || pho->isEE())) continue;
		bool matchesElectron(ConversionTools::hasMatchedPromptElectron(pho->superCluster(), gsfEleH, convH, beamSpotH->position()));
		bool matchesMuon(false);
                for(unsigned int ilep=0; ilep<lepColl.size(); ilep++)
                {
			if(fabs(lepColl[ilep].id)==11) continue;
                        if( deltaR( lepColl[ilep].eta(), lepColl[ilep].phi(), pho->eta(), pho->phi()) > 0.1) continue;
			matchesMuon=true;
			break;
		}
		if(matchesElectron || matchesMuon) continue;


                llvvPhoton gamma;

                //store information
                gamma.SetPxPyPzE(pho->px(), pho->py(), pho->pz(),  pho->energy());
//		gamma.pid                        = ev.egn;
		gamma.ecalIso03                  = pho->ecalRecHitSumEtConeDR03();
		gamma.hcalIso03                  = pho->hcalTowerSumEtConeDR03();
		gamma.trkIso03                   = pho->trkSumPtHollowConeDR03();
		gamma.ecalIso04                  = pho->ecalRecHitSumEtConeDR04();
		gamma.hcalIso04                  = pho->hcalTowerSumEtConeDR04();
		gamma.trkIso04                   = pho->trkSumPtHollowConeDR04();
		gIsolator03_.fGetIsolation(pho, &(*pfH), primVtx, vtxH);
		gamma.gIso03                     = gIsolator03_.getIsolationPhoton();
		gamma.chIso03                    = gIsolator03_.getIsolationCharged();
		gamma.puchIso03                  = 0;
		gamma.nhIso03                    = gIsolator03_.getIsolationNeutral();
		gIsolator04_.fGetIsolation(pho, &(*pfH), primVtx, vtxH);
		gamma.gIso04                     = gIsolator04_.getIsolationPhoton();
		gamma.chIso04                    = gIsolator04_.getIsolationCharged();
		gamma.puchIso04                  = 0;
		gamma.nhIso04                    = gIsolator04_.getIsolationNeutral();
		gamma.isConv                   = !( ConversionTools::matchedConversion(*(pho->superCluster()),convH,beamSpotH->position()).isNull() );
		gamma.sce                      = pho->superCluster()->energy();
		gamma.sceta                    = pho->superCluster()->eta();
		gamma.scphi                    = pho->superCluster()->phi();
		gamma.sihih                    = pho->sigmaIetaIeta();
		vector<float> cov                       = egLazyTool.localCovariances(*pho->superCluster()->seed()); 
		gamma.sipip                    = sqrt(cov[2]); 
		gamma.sihip                    = cov[1];
		gamma.r9                       = egLazyTool.e3x3(*pho->superCluster()->seed())/pho->superCluster()->rawEnergy();
		gamma.hoe                      = pho->hadTowOverEm();

		bool isLoose ( gamma.sihih<0.012 && gamma.hoe<0.05 );
		bool isMedium( gamma.sihih<0.011 && gamma.hoe<0.05 );
		bool isTight ( gamma.sihih<0.011 && gamma.hoe<0.05 );
		if(pho->isEE())
		{
			isLoose  = (gamma.sihih<0.034 && gamma.hoe<0.05 );
			isMedium = (gamma.sihih<0.033 && gamma.hoe<0.05 );
			isTight  = (gamma.sihih<0.031 && gamma.hoe<0.05 );
		}

		gamma.idbits    = (isLoose << 0) | (isMedium << 1 ) | (isTight << 2);
		if(isLoose && pho->isEB() && gamma.r9>0.9 && pho->pt()>20) nPhotons++;

                phoColl.push_back(gamma);
	}


  //////////////////////////////////   //////////////////////////////////   //////////////////////////////////   //////////////////////////////////
  //Filter out uninteresting events

  //now check if at least one trigger condition is fullfilled
  bool filterOut(true);
  for(unsigned int itrig=0; itrig<triggerCats.size(); itrig++){
      if(!triggerBits[itrig]) continue;
      int cat=triggerCats[itrig];
      if     (cat==11   && nElecs==0)                continue;
      else if(cat==13   && nMuons==0)                continue;
      else if(cat==22   && nPhotons==0)              continue;
      else if(cat==1111 && nElecs<2)                 continue;
      else if(cat==1113 && (nMuons==0 || nElecs==0)) continue;
      else if(cat==1313 && nMuons<2)                 continue;
      filterOut=false;
      break;
    }
  if(!isData && nHardProcGenLeptons>0 && nHardProcGenBosons>0 && keepFullGenInfo_) filterOut=false;
  if(filterOut) return false;

  //////////////////////////////////   //////////////////////////////////   //////////////////////////////////   //////////////////////////////////
  //Jets
	//
	// jets
	// https://twiki.cern.ch/twiki/bin/view/CMS/JetID
	// https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetID
	// https://twiki.cern.ch/twiki/bin/view/CMS/GluonTag
	// 
        if(analysisCfg_.getParameter<edm::InputTag>("jetSource").label()!=""){

		Handle<pat::JetCollection> jetH;
		iEvent.getByLabel( analysisCfg_.getParameter<edm::InputTag>("jetSource"), jetH);
		edm::Handle<edm::ValueMap<float> >  qgTaggerH;
		iEvent.getByLabel("QGTagger","qgLikelihood", qgTaggerH);
		edm::Handle<reco::JetTagCollection> tchpTags,   jpTags,    ssvheTags,    ivfTags,    origcsvTags,    csvTags,    jpcsvTags,    slcsvTags, supercsvTags;
		iEvent.getByLabel("trackCountingHighPurBJetTags",                  tchpTags);          
		iEvent.getByLabel("jetProbabilityBJetTags",                        jpTags);            
		iEvent.getByLabel("simpleSecondaryVertexHighEffBJetTags",          ssvheTags);  
		iEvent.getByLabel("simpleInclusiveSecondaryVertexHighEffBJetTags", ivfTags); 
		iEvent.getByLabel("combinedSecondaryVertexBJetTags",               origcsvTags);       
		iEvent.getByLabel("combinedSecondaryVertexRetrainedBJetTags",      csvTags);           
		iEvent.getByLabel("combinedCSVJPBJetTags",                         jpcsvTags);         
		iEvent.getByLabel("combinedCSVSLBJetTags",                         slcsvTags);         
		iEvent.getByLabel("combinedCSVJPSLBJetTags",                 supercsvTags);      
		edm::Handle<std::vector<reco::SecondaryVertexTagInfo> > svTagInfo, ivfTagInfo;
		iEvent.getByLabel("secondaryVertexTagInfos",                        svTagInfo);
		iEvent.getByLabel("inclusiveSecondaryVertexFinderTagInfosFiltered", ivfTagInfo);

		for(unsigned int ijet=0; ijet<jetH->size(); ++ijet)
		{
			edm::RefToBase<pat::Jet> jetRef(edm::Ref<pat::JetCollection>(jetH,ijet));
			const pat::Jet *patjet              = &((*jetH)[ijet]);
			const reco::Candidate *genParton = patjet->genParton();
			const reco::GenJet *genJet       = patjet->genJet();

			//pre-selection (note: raw jet energy must be used otherwise you'll have large inefficiencies for |eta|>3!!!!)
			float rawJetEn( patjet->correctedJet("Uncorrected").energy() );
			float nhf( (patjet->neutralHadronEnergy() + patjet->HFHadronEnergy())/rawJetEn );
			float nef( patjet->neutralEmEnergy()/rawJetEn );
			float cef( patjet->chargedEmEnergy()/rawJetEn );
			float chf( patjet->chargedHadronEnergy()/rawJetEn );
			float nch    = patjet->chargedMultiplicity();
			float nconst = patjet->numberOfDaughters();
			bool passLooseId(nhf<0.99  && nef<0.99 && nconst>1);
			bool passMediumId(nhf<0.95 && nef<0.95 && nconst>1);
			bool passTightId(nhf<0.90  && nef<0.90 && nconst>1);
			if(fabs(patjet->eta())<2.4) {
				passLooseId  &= (chf>0 && nch>0 && cef<0.99);
				passMediumId &= (chf>0 && nch>0 && cef<0.99);
				passTightId  &= (chf>0 && nch>0 && cef<0.99);
			}
			if(patjet->pt()<10 || fabs(patjet->eta())>4.7 /*|| !passLooseId*/) continue;

			//save information
			llvvJet jet;
			jet.SetPxPyPzE(patjet->px(), patjet->py(), patjet->pz(), patjet->energy());
			jet.torawsf     = patjet->correctedJet("Uncorrected").pt()/patjet->pt();
			jet.genflav     = patjet->partonFlavour();
			jet.genid       = genParton ? genParton->pdgId() : 0;
			if(genParton)jet.gen.SetPxPyPzE(genParton->px(), genParton->py(), genParton->pz(), genParton->energy());
	//		jet.genpx       = genParton ? genParton->px()    : 0;
	//		jet.genpy       = genParton ? genParton->py()    : 0;
	//		jet.genpz       = genParton ? genParton->pz()    : 0;
	//		jet.genen       = genParton ? genParton->energy(): 0;
			if(genJet)jet.genj.SetPxPyPzE(genJet->px(), genJet->py(), genJet->pz(), genJet->energy());
	//		jet.genjpx      = genJet    ? genJet->px()       : 0;
	//		jet.genjpy      = genJet    ? genJet->py()       : 0;
	//		jet.genjpz      = genJet    ? genJet->pz()       : 0;
	//		jet.genjen      = genJet    ? genJet->energy()   : 0;
			jet.neutHadFrac = patjet->neutralHadronEnergyFraction();
			jet.neutEmFrac  = patjet->neutralEmEnergyFraction();
			jet.chHadFrac   = patjet->chargedHadronEnergyFraction();
			jet.muFrac      = patjet->muonEnergyFraction();
			jet.area        = patjet->jetArea();

			jet.tchp        = (*tchpTags)[ijet].second;
			jet.jp          = (*jpTags)[ijet].second;
			jet.ssvhe       = (*ssvheTags)[ijet].second;
			jet.ivf         = (*ivfTags)[ijet].second;
			jet.origcsv     = (*origcsvTags)[ijet].second;
			jet.csv         = (*csvTags)[ijet].second;
			jet.jpcsv       = (*jpcsvTags)[ijet].second;
			jet.slcsv       = (*slcsvTags)[ijet].second;
			jet.supercsv    = (*supercsvTags)[ijet].second;

			//secondary vertex from associated tracks
			if(svTagInfo.isValid() && svTagInfo->size()>ijet);
			{
				const reco::SecondaryVertexTagInfo &sv= (*svTagInfo)[ijet];
				int nsvtx=sv.nVertices();
				jet.svxNtrk=0;
				jet.svxLxy=0;
				jet.svxLxyErr=0;
				jet.svxM=0;
				jet.svxPx=0;
				jet.svxPy=0;
				jet.svxPz=0;
				if(nsvtx)
				{  
					for (reco::Vertex::trackRef_iterator titt = sv.secondaryVertex(0).tracks_begin(); titt != sv.secondaryVertex(0).tracks_end(); titt++) jet.svxNtrk++;
					jet.svxLxy    = sv.flightDistance(0).value();
					jet.svxLxyErr = sv.flightDistance(0).error();
					jet.svxM      = sv.secondaryVertex(0).p4().mass();
					jet.svxPx     = sv.secondaryVertex(0).p4().px();
					jet.svxPy     = sv.secondaryVertex(0).p4().py();
					jet.svxPz     = sv.secondaryVertex(0).p4().pz();
				}
			}

			//secondary vertex from inclusive tracks
			if(ivfTagInfo.isValid() && ivfTagInfo->size()>ijet)
			{
				const reco::SecondaryVertexTagInfo &sv= (*ivfTagInfo)[ijet];
				int nsvtx=sv.nVertices();
				jet.ivfNtrk=0;
				jet.ivfLxy=0;
				jet.ivfLxyErr=0;
				jet.ivfM=0;
				jet.ivfPx=0;
				jet.ivfPy=0;
				jet.ivfPz=0;
				if(nsvtx)
				{
					for (reco::Vertex::trackRef_iterator titt = sv.secondaryVertex(0).tracks_begin(); titt != sv.secondaryVertex(0).tracks_end(); titt++) jet.ivfNtrk++;
					jet.ivfLxy    = sv.flightDistance(0).value();
					jet.ivfLxyErr = sv.flightDistance(0).error();
					jet.ivfM      = sv.secondaryVertex(0).p4().mass();
					jet.ivfPx     = sv.secondaryVertex(0).p4().px();
					jet.ivfPy     = sv.secondaryVertex(0).p4().py();
					jet.ivfPz     = sv.secondaryVertex(0).p4().pz();
				}
			}

			PileupJetIdentifier cutBasedPuIdentifier = cutBasedPuJetIdAlgo_->computeIdVariables(dynamic_cast<const reco::Jet*>(patjet), 0, primVtx.get(), *vtxH.product(), true);
			PileupJetIdentifier puIdentifier         = puJetIdAlgo_->computeIdVariables(dynamic_cast<const reco::Jet*>(patjet), 0, primVtx.get(), *vtxH.product(), true);
			jet.beta        = puIdentifier.beta();
			jet.betaStar    = puIdentifier.betaStar();
			jet.dRMean      = puIdentifier.dRMean();
			jet.dR2Mean     = puIdentifier.dR2Mean();
			jet.ptRMS       = puIdentifier.ptRMS();
			jet.ptD         = puIdentifier.ptD();
			jet.etaW        = puIdentifier.etaW();
			jet.phiW        = puIdentifier.phiW();
			jet.puMVA       = puIdentifier.mva();
			jet.qgMVA       = qgTaggerH.isValid() ? (*qgTaggerH)[jetRef] : 0;

			//save pf constituents (only for jets with pT>20 in the central region)
			jet.pfstart=-1;
			jet.pfend=-1;
			if(keepPfCandidates>0 && patjet->pt()>20 && fabs(patjet->eta())<3)
			{
				const std::vector<reco::PFCandidatePtr> pfConst = patjet->getPFConstituents();
				jet.pfstart=pfColl.size();
				for(unsigned int ipf=0; ipf<pfConst.size(); ipf++)
				{
					llvvPFParticle pfPart;
					pfPart.SetPxPyPzE(pfConst[ipf]->px(), pfConst[ipf]->py(), pfConst[ipf]->pz(), pfConst[ipf]->energy());
					pfPart.id     = pfConst[ipf]->pdgId();
					pfPart.charge = pfConst[ipf]->charge();
					pfColl.push_back(pfPart);
				}
				jet.pfend=pfColl.size()-1;
			}

			//a summary of the id bits
			jet.idbits =
				(passLooseId << 0)
				| (passMediumId << 1)
				| (passTightId << 2)
				| ( ( uint(puIdentifier.idFlag()) & 0xf ) << 3 )
				| ( ( uint(cutBasedPuIdentifier.idFlag()) & 0xf ) << 7 );

                        if(jet.pt()>jetPtCut){
   		          jetColl.push_back(jet);
                        }
		}
        }

  //////////////////////////////////   //////////////////////////////////   //////////////////////////////////   //////////////////////////////////
  //PF Candidates 
        //
        // charged PF candidates which haven't been clustered
        //
        if(keepPfCandidates>1){
           for(size_t ipf=0; ipf<pfH->size(); ipf++)
           {
                   const reco::PFCandidate &cand=(*pfH)[ipf];

                   //require charged and with track associated
                   if(cand.charge()==0) continue;
                   reco::TrackBaseRef trackBaseRef( cand.trackRef() );
                   if(trackBaseRef.isNull()) continue;

                   //minimum pT of 300 MeV
                   if(cand.pt()<0.3) continue;

                   //check for overlaps
                   bool matches(false);
                   for(unsigned int jpf=0; jpf<pfColl.size(); jpf++)
                   {
                           if( deltaR( pfColl[jpf].eta(), pfColl[jpf].phi(), cand.eta(), cand.phi()) > 0.1) continue;
                           matches=true;
                           break;
                   }
                   if(matches) continue;

                   //require it to be associated to the primary vertex
                   int bestVtx(-1);
                   float bestDz(9999.);
                   if(trackBaseRef.isAvailable())
                   {
                           for(size_t jVtx=0; jVtx<vtxH->size(); jVtx++)
                           {
                                   const reco::VertexRef vtxref(vtxH,jVtx);
                                   float vtxDz( fabs( trackBaseRef->dz( vtxref->position()) ) );
                                   if(vtxDz > bestDz) continue;
                                   bestDz=vtxDz;
                                   bestVtx=jVtx;
                           }
                   }
                   if(bestVtx!=0) continue;

                   llvvPFParticle pfPart;
                   pfPart.SetPxPyPzE(cand.px(), cand.py(), cand.pz(), cand.energy());
                   pfPart.id     = cand.pdgId();
                   pfPart.charge = cand.charge();
                   if(pfPart.pt()>pfCandPtCut && bestDz<pfCandDzMax)pfColl.push_back(pfPart);
           }
       }


  //////////////////////////////////   //////////////////////////////////   //////////////////////////////////   //////////////////////////////////
  //MET
        //
        // MET
        //

     std::vector<edm::InputTag> metSources=analysisCfg_.getParameter<std::vector<edm::InputTag> >("metSource");
     for(unsigned int imet=0;imet<metSources.size();imet++){
        Handle<View<reco::PFMET> > metH;
        iEvent.getByLabel(metSources[imet], metH);
        if(!metH.isValid())continue;

        //Output object in EDM format
        std::auto_ptr<llvvMet> metOut(new llvvMet());
        llvvMet& met = *metOut;

        //////////////////////////////////

          met.SetPxPyPzE(metH->ptrAt(0)->px(), metH->ptrAt(0)->py(), metH->ptrAt(0)->pz(), metH->ptrAt(0)->energy());
          met.sigx2 = metH->ptrAt(0)->getSignificanceMatrix()(0,0);
          met.sigxy = metH->ptrAt(0)->getSignificanceMatrix()(0,1);
          met.sigy2 = metH->ptrAt(0)->getSignificanceMatrix()(1,1);
          met.sig   = metH->ptrAt(0)->significance();

         iEvent.put(metOut, metSources[imet].label()); //save the object to the event here, to keep it in the loop
       }


  //////////////////////////////////   //////////////////////////////////   //////////////////////////////////   //////////////////////////////////
  //Save the object to the EDM EVENT


       // adding llvvMuonInfo to the EVENT and add Ref to MuonInfo Object to the Lepton Object
       edm::OrphanHandle<llvvMuonInfoCollection> muInfoCollHandle= iEvent.put(muInfoCollOut);
       int NMuons=0;
       for(int i=0;i<(int)lepColl.size();i++) {
          if(abs(lepColl[i].id)!=13)continue;
          lepColl[i].muonInfoRef = llvvMuonInfoRef(muInfoCollHandle,NMuons);
          NMuons++;
       }

       // adding llvvElectronInfo to the EVENT and add Ref to ElectronInfo Object to the Lepton Object
       edm::OrphanHandle<llvvElectronInfoCollection> elInfoCollHandle= iEvent.put(elInfoCollOut);
       int NElectrons=0;
       for(int i=0;i<(int)lepColl.size();i++) {
          if(abs(lepColl[i].id)!=11)continue;
          lepColl[i].electronInfoRef = llvvElectronInfoRef(elInfoCollHandle,NElectrons);
          NElectrons++;
       }

       //adding the lepton collection (including all references) to the EVENT
       iEvent.put(genEvOut);
       iEvent.put(lepCollOut);         
       iEvent.put(tauCollOut);
       iEvent.put(phoCollOut);
       iEvent.put(jetCollOut);
       iEvent.put(pfCollOut);
       iEvent.put(genPartCollOut);
       iEvent.put(triggerBitsOut, "triggerBits");
       iEvent.put(triggerPrescalesOut, "triggerPrescales");
       iEvent.put(nvtxOut, "nvtx");
       return true;
}




DEFINE_FWK_MODULE(llvvObjectProducers);

