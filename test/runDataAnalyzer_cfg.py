import FWCore.ParameterSet.Config as cms

process = cms.Process("DataAna")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = gtag

process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True),
                                        SkipEvent = cms.untracked.vstring('ProductNotFound')
                                        ) 

#the source and output
try:
    print 'Processing %d inputs'%len(inputList)
except:
    if isMC : inputList = cms.untracked.vstring('/store/relval/CMSSW_5_3_6-START53_V14/RelValTTbar/GEN-SIM-RECO/v2/00000/16D5D599-F129-E211-AB60-00261894390B.root')
    else    : inputList = cms.untracked.vstring('/store/data//Run2012A/DoubleMu/AOD//22Jan2013-v1/20000/F4C34C30-B581-E211-8269-003048FFD7A2.root') 
process.source = cms.Source("PoolSource",
                            fileNames = inputList
                            )

try:
    print 'Will store all PF candidates?' + str(storeAllPF)
except:
    storeAllPF=True
    print 'Will store all PF candidates?' + str(storeAllPF)
    
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

try:
    print 'EDM output set to %s'%outFile
except:
    outFile='Events.root'
process.out = cms.OutputModule("PoolOutputModule",
                               outputCommands = cms.untracked.vstring('keep *'),
                               fileName = cms.untracked.string(outFile)
                               )

if(isMC) : process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

try:
    print 'TFileService output set to %s'%(tfsOutputFile)
except:
    tfsOutputFile="DataAnalysis.root"
process.TFileService = cms.Service("TFileService", fileName = cms.string(tfsOutputFile))



##-------------------- Import the JEC services -----------------------
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

##-------------------- Import the Jet RECO modules -----------------------
process.load('RecoJets.Configuration.RecoPFJets_cff')

process.kt6PFJets.doRhoFastjet = True
process.ak5PFJets.doAreaFastjet = True

#apply a good vertex selector and filter out scraping
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter( "PrimaryVertexObjectFilter",
                                                   filterParams = pvSelector.clone( minNdof = cms.double(4.0),
                                                                                    maxZ = cms.double(24.0),
                                                                                    maxd0 = cms.double(2.0)
                                                                                    ),
                                                   src=cms.InputTag('offlinePrimaryVertices')
                                                   )
process.goodVertexFilter = cms.EDFilter("GoodVertexFilter",
                                        vertexCollection = cms.InputTag('goodOfflinePrimaryVertices'),
                                        minimumNDOF = cms.uint32(4),
                                        maxAbsZ = cms.double(24),
                                        maxd0 = cms.double(2)
                                        )

process.noscraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
                                  debugOn = cms.untracked.bool(False),
                                  numtrack = cms.untracked.uint32(10),
                                  thresh = cms.untracked.double(0.25)
                               )

# optional MET filters
# cf.https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFilters
process.load('RecoMET.METFilters.metFilters_cff')
process.hcalLaserEventFilter.taggingMode   = cms.bool(True)
process.EcalDeadCellTriggerPrimitiveFilter.taggingMode=cms.bool(True)
process.eeBadScFilter.taggingMode           = cms.bool(True)
process.ecalLaserCorrFilter.taggingMode     = cms.bool(True)
process.trackingFailureFilter.VertexSource  = cms.InputTag('goodOfflinePrimaryVertices')
process.trackingFailureFilter.taggingMode   = cms.bool(True)
process.manystripclus53X.taggedMode         = cms.untracked.bool(True)
process.manystripclus53X.forcedValue        = cms.untracked.bool(False)
process.toomanystripclus53X.taggedMode      = cms.untracked.bool(True)
process.toomanystripclus53X.forcedValue     = cms.untracked.bool(False)
process.logErrorTooManyClusters.taggedMode  = cms.untracked.bool(True)
process.logErrorTooManyClusters.forcedValue = cms.untracked.bool(False)  

process.metFilteringTaggers = cms.Sequence(process.HBHENoiseFilter*
                                           process.hcalLaserEventFilter *
                                           process.EcalDeadCellTriggerPrimitiveFilter *
                                           process.eeBadScFilter *
                                           process.ecalLaserCorrFilter *
                                           process.trackingFailureFilter *
                                           process.trkPOGFilters)

#PF2PAT
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.tools.pfTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *

postfix = "PFlow"
jetAlgo="AK5"
jecLevels=['L1FastJet', 'L2Relative', 'L3Absolute']
if(not isMC): jecLevels.append('L2L3Residual')

usePF2PAT(process,
          runPF2PAT=True,
          jetAlgo=jetAlgo,
          runOnMC=isMC,
          postfix=postfix,
          jetCorrections=('AK5PFchs', jecLevels),
          pvCollection=cms.InputTag('goodOfflinePrimaryVertices'),
          typeIMetCorrections=False)

#setup trigger matching
from UserCode.llvv_fwk.triggerMatching_cfg import *
addTriggerMatchingTo(process)

#custom electrons
useGsfElectrons(process,postfix=postfix,dR="03")
process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')
process.eidMVASequence = cms.Sequence(  process.mvaTrigV0 + process.mvaNonTrigV0 )
process.patElectronsPFlow.electronIDSources.mvaTrigV0    = cms.InputTag("mvaTrigV0")
process.patElectronsPFlow.electronIDSources.mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")
from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import *
process.selectedPatElectronsPFlowHeep = cms.EDProducer("HEEPAttStatusToPAT",
                                                       eleLabel = cms.InputTag("selectedPatElectronsWithTrigger"),
                                                       barrelCuts = cms.PSet(heepBarrelCuts),
                                                       endcapCuts = cms.PSet(heepEndcapCuts),
                                                       applyRhoCorrToEleIsol = cms.bool(True),
                                                       eleIsolEffectiveAreas = cms.PSet (heepEffectiveAreas),
                                                       eleRhoCorrLabel = cms.InputTag("kt6PFJets:rho"),
                                                       verticesLabel = cms.InputTag("goodOfflinePrimaryVertices"),
                                                       )

#custom muons
process.patMuonsPFlow.pfMuonSource = cms.InputTag("pfSelectedMuonsPFlow")
process.muonMatchPFlow.src = cms.InputTag("pfSelectedMuonsPFlow")

#custom jets for CHS
process.pfPileUpPFlow.checkClosestZVertex = cms.bool(False)
process.pfPileUpIsoPFlow.checkClosestZVertex = cms.bool(False)
getattr(process,"pfNoPileUp"+postfix).enable = True
getattr(process,"pfNoMuon"+postfix).enable = False     # to use muon-cleaned electron collection set to True (check iso)
getattr(process,"pfNoElectron"+postfix).enable = False # to use electron-cleaned tau collection set to True (check iso)
getattr(process,"pfNoTau"+postfix).enable = False      # to use tau-cleaned jet collection set to True (check what is a tau)
getattr(process,"pfNoJet"+postfix).enable = True       # this i guess it's for photons...      

#add q/g discriminator
process.load('QuarkGluonTagger.EightTeV.QGTagger_RecoJets_cff')
process.QGTagger.srcJets    = cms.InputTag("selectedPatJets"+postfix)
process.QGTagger.isPatJet  = cms.untracked.bool(True)
process.QGTagger.useCHS    = cms.untracked.bool(True) 
process.QGTagger.srcRho    = cms.InputTag('kt6PFJets','rho')
process.QGTagger.srcRhoIso = cms.InputTag('kt6PFJetsCentral','rho')
process.qgSequence=cms.Sequence(process.goodOfflinePrimaryVerticesQG+process.QGTagger)

#compute rho from central pf candidates only
from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
process.kt6PFJetsCentral = kt4PFJets.clone( rParam = cms.double(0.6),
                                            doAreaFastjet = cms.bool(True),
                                            doRhoFastjet = cms.bool(True),
                                            Rho_EtaMax = cms.double(2.5),
                                            Ghost_EtaMax = cms.double(2.5) )

from UserCode.llvv_fwk.btvDefaultSequence_cff import *
btvDefaultSequence(process,isMC,"selectedPatJets"+postfix,"goodOfflinePrimaryVertices")

# cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMetAnalysis
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
process.load("JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi")
process.pfType1CorrectedMet.applyType0Corrections = cms.bool(False)
process.pfType1CorrectedMet.srcType1Corrections = cms.VInputTag( cms.InputTag('pfMETcorrType0'),
                                                                 cms.InputTag('pfJetMETcorr', 'type1')
                                                                 )


#the analyzer
from UserCode.llvv_fwk.dataAnalyzer_cfi import *
try:
    if runDijetsAnalysis :
        process.dataAnalyzer = dijetAnalyzer.clone()
        print 'Running dijet version of the dataAnalyzer'
except:
    process.dataAnalyzer = dataAnalyzer.clone()
    print 'Running standard dataAnalyzer'
    if isTauEmbed:
        process.dataAnalyzer.cfg.triggerCats[2]=1113
        process.dataAnalyzer.cfg.triggerCats[3]=1113
        print 'Tweaking for tau embedded samples'
        
process.dataAnalyzer.cfg.storeAllPF=storeAllPF

try:
    if doUnfold:
        process.dataAnalyzer.cfg.keepFullGenInfo = cms.bool(True)
        print 'Warning running unfolding: ntuples will save lots of information on gen level and won`t use trigger'
except:
    print 'Basic generator level information will be stored (only status=3 + photons status 1)'
                
#counters for specific filters
process.startCounter = cms.EDProducer("EventCountProducer")
process.scrapCounter = process.startCounter.clone()
process.vtxCounter   = process.startCounter.clone()
process.metCounter   = process.startCounter.clone() 
process.p = cms.Path( process.startCounter
                      *process.noscraping
                      *process.scrapCounter
                      *process.goodOfflinePrimaryVertices
                      *process.goodVertexFilter
                      *process.vtxCounter
                      *process.metFilteringTaggers
		      *process.metCounter
                      *process.eidMVASequence
                      *getattr(process,"patPF2PATSequence"+postfix)
                      *process.btvSequence
                      *process.kt6PFJetsCentral
                      *process.qgSequence
                      *process.type0PFMEtCorrection*process.producePFMETCorrections
                      *process.selectedPatElectronsWithTrigger*process.selectedPatElectronsPFlowHeep
                      *process.selectedPatMuonsTriggerMatch 
                      *process.dataAnalyzer
                      )

	




