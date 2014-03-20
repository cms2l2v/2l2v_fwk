import FWCore.ParameterSet.Config as cms

from SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi import *
from SimTracker.TrackAssociation.TrackAssociatorByHits_cfi import *
from RecoBTag.ImpactParameter.impactParameter_cfi import *

def btvDefaultSequence(process, isMC=True, jetCollection="selectedPatJetsPFlow",vtxCollection="goodOfflinePrimaryVertices") :

    process.load("PhysicsTools.JetMCAlgos.CaloJetsMCFlavour_cfi")
    process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
    process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
    process.load("SimTracker.TrackHistory.TrackHistory_cff")
    process.load("SimTracker.TrackHistory.TrackClassifier_cff")
    process.load("RecoBTau.JetTagComputer.jetTagRecord_cfi")
    
    ############################################################################
    # for Impact Parameter based taggers
    ############################################################################
    process.load("RecoBTag.ImpactParameter.negativeOnlyJetBProbabilityComputer_cfi")
    process.load("RecoBTag.ImpactParameter.negativeOnlyJetProbabilityComputer_cfi")
    process.load("RecoBTag.ImpactParameter.positiveOnlyJetProbabilityComputer_cfi")
    process.load("RecoBTag.ImpactParameter.positiveOnlyJetBProbabilityComputer_cfi")
    process.load("RecoBTag.ImpactParameter.negativeTrackCounting3D2ndComputer_cfi")
    process.load("RecoBTag.ImpactParameter.negativeTrackCounting3D3rdComputer_cfi")
    process.load("RecoBTag.Configuration.RecoBTag_cff")
    process.load("RecoJets.JetAssociationProducers.ak5JTA_cff")
    process.load("RecoBTag.ImpactParameter.negativeOnlyJetBProbabilityJetTags_cfi")
    process.load("RecoBTag.ImpactParameter.negativeOnlyJetProbabilityJetTags_cfi")
    process.load("RecoBTag.ImpactParameter.positiveOnlyJetProbabilityJetTags_cfi")
    process.load("RecoBTag.ImpactParameter.positiveOnlyJetBProbabilityJetTags_cfi")
    process.load("RecoBTag.ImpactParameter.negativeTrackCountingHighPur_cfi")
    process.load("RecoBTag.ImpactParameter.negativeTrackCountingHighEffJetTags_cfi")
    process.load("RecoBTag.ImpactParameter.jetProbabilityBJetTags_cfi")
    process.load("RecoBTag.ImpactParameter.jetBProbabilityBJetTags_cfi")

    ############################################################################
    # for Secondary Vertex taggers
    ############################################################################
    process.load("RecoBTag.SecondaryVertex.secondaryVertexTagInfos_cfi")
    process.load("RecoBTag.SecondaryVertex.secondaryVertexNegativeTagInfos_cfi")
    process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexHighEffBJetTags_cfi")
    process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexHighPurBJetTags_cfi")
    process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexNegativeHighEffBJetTags_cfi")
    process.load("RecoBTag.SecondaryVertex.simpleSecondaryVertexNegativeHighPurBJetTags_cfi")
    process.load("RecoBTag.SecondaryVertex.combinedSecondaryVertexNegativeBJetTags_cfi")
    process.load("RecoBTag.SecondaryVertex.combinedSecondaryVertexNegativeES_cfi")
    process.load("RecoBTag.SecondaryVertex.combinedSecondaryVertexPositiveBJetTags_cfi")
    process.load("RecoBTag.SecondaryVertex.combinedSecondaryVertexPositiveES_cfi")

    ############################################################################
    # for Inclusive Secondary Vertexing
    ############################################################################
    process.load("RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff")
    process.load("RecoBTag.SecondaryVertex.bToCharmDecayVertexMerger_cfi")
    process.load("RecoBTag.SecondaryVertex.inclusiveSecondaryVertexFinderTagInfos_cfi")
    process.load("RecoBTag.SecondaryVertex.simpleInclusiveSecondaryVertexBJetTags_cfi")

    process.combinedInclusiveSecondaryVertexPositiveBJetTags = process.combinedInclusiveSecondaryVertexBJetTags.clone()
    process.combinedInclusiveSecondaryVertexPositiveBJetTags.jetTagComputer = cms.string('combinedSecondaryVertexPositive')

    ############################################################################
    # For the Retrained CSV
    ############################################################################
    process.combinedSecondaryVertexRetrained = process.combinedSecondaryVertex.clone(
      calibrationRecords = cms.vstring(
      'CombinedSVRetrainRecoVertex',
      'CombinedSVRetrainPseudoVertex',
      'CombinedSVRetrainNoVertex'
      )
      )
    process.combinedSecondaryVertexRetrainedBJetTags = process.combinedSecondaryVertexBJetTags.clone(
      jetTagComputer = cms.string('combinedSecondaryVertexRetrained')
      )

    ##############################################
    # Get calibrations for the CSV taggers
    ##############################################
    process.load("CondCore.DBCommon.CondDBSetup_cfi")
    process.BTauMVAJetTagComputerRecord = cms.ESSource("PoolDBESSource",
                                                       process.CondDBSetup,
                                                       timetype = cms.string('runnumber'),
                                                       toGet = cms.VPSet(cms.PSet(
        record = cms.string('BTauGenericMVAJetTagComputerRcd'),
        tag = cms.string('MVAComputerContainer_Retrained53X_JetTags_v2')
        )),
                                                       connect = cms.string('frontier://FrontierProd/CMS_COND_PAT_000'),
                                                       BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService')
                                                       )
    
    process.es_prefer_BTauMVAJetTagComputerRecord = cms.ESPrefer("PoolDBESSource","BTauMVAJetTagComputerRecord") 

    ### JP calibration for MC only 
#    if(isMC) :
#        process.GlobalTag.toGet = cms.VPSet(
#            cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
#                     tag = cms.string("TrackProbabilityCalibration_2D_MC53X_v2"),
#                     connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
#            cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
#                     tag = cms.string("TrackProbabilityCalibration_3D_MC53X_v2"),
#                     connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
#            )
#    else:
#        process.GlobalTag.toGet = cms.VPSet(
#            cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
#                     tag = cms.string("TrackProbabilityCalibration_2D_Data53X_v2"),
#                     connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
#            cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
#                     tag = cms.string("TrackProbabilityCalibration_3D_Data53X_v2"),
#                     connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
#            )
        
    #user specific configurations
    process.AK5byRef.jets                       = jetCollection
    process.ak5JetTracksAssociatorAtVertex.jets = jetCollection
    process.softMuonTagInfos.jets               = jetCollection 
    process.softPFMuonsTagInfos.primaryVertex = cms.InputTag(vtxCollection)
    process.softPFMuonsTagInfos.jets            = jetCollection 
    process.softPFElectronsTagInfos.primaryVertex = cms.InputTag(vtxCollection)
    process.softPFElectronsTagInfos.jets        = jetCollection

    


    #the sequence
    process.mainBtvSequence=cms.Sequence(process.inclusiveVertexing*process.inclusiveMergedVerticesFiltered*process.bToCharmDecayVertexMerged
                                         *process.ak5JetTracksAssociatorAtVertex
                                         *process.btagging
                                         *process.negativeTrackCountingHighEffJetTags*process.negativeTrackCountingHighPur
                                         *process.secondaryVertexTagInfos*process.simpleSecondaryVertexHighEffBJetTags*process.simpleSecondaryVertexHighPurBJetTags
                                         *process.secondaryVertexNegativeTagInfos*process.simpleSecondaryVertexNegativeHighEffBJetTags*process.simpleSecondaryVertexNegativeHighPurBJetTags
                                         *process.combinedSecondaryVertexRetrainedBJetTags
                                         *process.inclusiveSecondaryVertexFinderFilteredTagInfos
                                         *process.inclusiveSecondaryVertexFinderTagInfos
                                         *process.simpleInclusiveSecondaryVertexHighEffBJetTags
                                         *process.simpleInclusiveSecondaryVertexHighPurBJetTags
                                         #*process.doubleSecondaryVertexHighEffBJetTags
                                         *process.combinedInclusiveSecondaryVertexBJetTags)
                                         
    if(isMC) : process.btvSequence=cms.Sequence(process.myPartons*process.AK5Flavour*process.mainBtvSequence)
    else     : process.btvSequence=cms.Sequence(process.mainBtvSequence)
