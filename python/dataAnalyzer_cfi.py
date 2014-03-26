import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.PileupJetID_cfi  import pileupJetIdProducerChs
#pileupJetIdProducerChs.algos[0].tmvaWeights=cms.string("RecoJets/JetProducers/data/TMVAClassificationCategory_JetID_53X_chs_Dec2012.weights.xml")



dataAnalyzer = cms.EDAnalyzer( "DataAnalyzer",
                               cfg=cms.PSet( metFilters=cms.vstring(
                                                                    'HBHENoiseFilter',
                                                                    'hcalLaserEventFilter',
                                                                    'EcalDeadCellTriggerPrimitiveFilter',
                                                                    'eeBadScFilter',
                                                                    'ecalLaserCorrFilter',
                                                                    'trackingFailureFilter'),
                                             triggerSource = cms.InputTag("TriggerResults::HLT"),
                                             triggerPaths = cms.vstring('HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v',
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
								                                                        'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v'
                                                                        ),
                                             triggerCats  = cms.vint32(1111,
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
								                                                       11
                                                                       ),
                                             genSource           = cms.InputTag("genParticles"),
                                             keepFullGenInfo     = cms.bool(False),
                                             vtxSource           = cms.InputTag("goodOfflinePrimaryVertices"),
                                             beamSpotSource      = cms.InputTag("offlineBeamSpot"),
                                             pfSource            = cms.InputTag("particleFlow"),
                                             storeAllPF          = cms.bool(True),
                                             skipMCTrigSelection = cms.bool(True),
                                             muonSource          = cms.InputTag("selectedPatMuonsTriggerMatch"),
                                             electronSource      = cms.InputTag("selectedPatElectronsWithTrigger"),
                                             # electronSource      = cms.InputTag("selectedPatElectronsPFlowHeep"),
                                             photonSource        = cms.InputTag("photons"),
                                             conversionSource    = cms.InputTag("allConversions"),
                                             ebrechitsSource     = cms.InputTag("reducedEcalRecHitsEB"),
                                             eerechitsSource     = cms.InputTag("reducedEcalRecHitsEE"),
                                             jetSource           = cms.InputTag("selectedPatJetsPFlow"),
                                             pujetidAlgo         = pileupJetIdProducerChs.algos,
                                             metSource           = cms.VInputTag("pfMETPFlow","pfMet","pfType1CorrectedMet","pfType1p2CorrectedMet"),
                                             rhoSource           = cms.InputTag("kt6PFJets:rho"),
                                             rho25Source         = cms.InputTag("kt6PFJetsCentral:rho")
                                             )
                               )


#######################################
# 7 TeV trigger pathcs and categories #
#######################################
triggerPaths7TeV = cms.vstring('HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v',                                       #ee: 0-1
                               'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v',
                               'HLT_DoubleMu7_v',                                                                            #mumu: 2-5
                               'HLT_Mu13_Mu8_v',
                               'HLT_Mu17_Mu8_v',
                               'HLT_Mu17_TkMu8_v',
                               'HLT_Mu17_Ele8_CaloIdL_v',                                                                    #emu:6-9
                               'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v',
                               'HLT_Mu8_Ele17_CaloIdL_v',
                               'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v',
                               'HLT_IsoMu17_v',                                                                              #mu:10-12
                               'HLT_IsoMu24_v',
                               'HLT_IsoMu30_eta2p1_v',
                               'HLT_Photon30_CaloIdVL_IsoL_v',                                                               #gamma:13-20
                               'HLT_Photon50_CaloIdVL_IsoL_v',
                               'HLT_Photon75_CaloIdVL_IsoL_v',
                               'HLT_Photon90_CaloIdVL_IsoL_v',
                               'HLT_Photon125_v',
                               'HLT_Photon125_NoSpikeFilter_v',
                               'HLT_Photon135_v',
                               'HLT_Photon200_NoHE_v'
                               ),
triggerCats7TeV  = cms.vint32(1111,
                              1111,
                              1313,
                              1313,
                              1313,
                              1313,
                              1113,
                              1113,
                              1113,
                              1113,
                              13,
                              13,
                              13,
                              22,
                              22,
                              22,
                              22,
                              22,
                              22,
                              22,
                              22
                              )


## configure specifically for a dijet analysis
dijetAnalyzer = dataAnalyzer.clone()
dijetAnalyzer.cfg.triggerPaths=cms.vstring("BTagMu_DiJet20","BTagMu_DiJet40","BTagMu_DiJet70","BTagMu_DiJet110","BTagMu_Jet300")
dijetAnalyzer.cfg.triggerCats=cms.vint32(1,1,1,1,1)



