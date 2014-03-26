import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.PileupJetID_cfi  import pileupJetIdProducerChs
#pileupJetIdProducerChs.algos[0].tmvaWeights=cms.string("RecoJets/JetProducers/data/TMVAClassificationCategory_JetID_53X_chs_Dec2012.weights.xml")  

llvvGenParticleProducer = cms.EDFilter( "llvvGenParticleProducer",
   genSource       = cms.InputTag("genParticles"),
)


llvvObjectProducers = cms.EDFilter( "llvvObjectProducers",
                     genSource        = cms.InputTag("genParticles"),
                     keepFullGenInfo  = cms.bool(False),
                     vtxSource        = cms.InputTag("goodOfflinePrimaryVertices"),
                     beamSpotSource   = cms.InputTag("offlineBeamSpot"),
                     pfSource         = cms.InputTag("particleFlow"),
                     tauSource        = cms.InputTag("selectedPatTausPFlow"),
                     boostedTauSource = cms.InputTag("patTausBoost"),
                     muonSource       = cms.InputTag("selectedPatMuonsTriggerMatch"),
                     electronSource   = cms.InputTag("selectedPatElectronsWithTrigger"),
                     photonSource     = cms.InputTag("photons"),
                     conversionSource = cms.InputTag("allConversions"),
                     ebrechitsSource  = cms.InputTag("reducedEcalRecHitsEB"),
                     eerechitsSource  = cms.InputTag("reducedEcalRecHitsEE"),
                     rhoSource        = cms.InputTag("kt6PFJets:rho"),
                     rho25Source      = cms.InputTag("kt6PFJetsCentral:rho"),
                     jetSource        = cms.InputTag("selectedPatJetsPFlow"),
                     pujetidAlgo      = pileupJetIdProducerChs.algos,
                     keepPfCandidates = cms.int32(0), #0 PFCandidates are not saved, #1 save PF candidates in Jets, #2 Save all with pT>0.3
                     metSource        = cms.VInputTag("pfMETPFlow","pfMet","pfType1CorrectedMet","pfType1p2CorrectedMet", "pfMEtMVA"),
                     triggerSource    = cms.InputTag("TriggerResults::HLT"),
                     triggerPaths     = cms.vstring(
                                             'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v',
                                             'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v',
                                             'HLT_Mu17_Mu8_v',
                                             'HLT_Mu17_TkMu8_v',
                                             'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v',
                                             'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v',
                                             'HLT_IsoMu24_eta2p1_v',
                                             'HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_v',
                                             'HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_v',
                                             'HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_v',
                                             'HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_v',
                                             'HLT_Photon250_NoHE_v1_v',
                                             'HLT_Photon300_NoHE_v1_v',
                                             'HLT_Ele27_WP80_v',
                                             'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v',
                                             'HLT_IsoMu24_v'
                                             ),
                     triggerCats      = cms.vint32(
                                            1111,
                                            1111,
                                            1313,
                                            1313,
                                            1113,
                                            1113,
                                            13,
                                            22,
                                            22,
                                            22,
                                            22,
                                            22,
                                            22,
                                            11,
                                            11,
                                            13
                                            ),

)

## configure specifically for a dijet analysis
dijetObjectProducers = llvvObjectProducers.clone()
dijetObjectProducers.triggerPaths=cms.vstring("BTagMu_DiJet20","BTagMu_DiJet40","BTagMu_DiJet70","BTagMu_DiJet110","BTagMu_Jet300")
dijetObjectProducers.triggerCats=cms.vint32(1,1,1,1,1)
