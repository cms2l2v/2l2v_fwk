import FWCore.ParameterSet.Config as cms

process = cms.Process("DataAna")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "START53_V23::All"

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
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('/store/relval/CMSSW_5_3_6-START53_V14/RelValTTbar/GEN-SIM-RECO/v2/00000/16D5D599-F129-E211-AB60-00261894390B.root')
                            )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.out = cms.OutputModule("PoolOutputModule",
                               outputCommands = cms.untracked.vstring('keep *'),
                               fileName = cms.untracked.string("Events.root")
                               )


process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.TFileService = cms.Service("TFileService", fileName = cms.string("DataAnalysis.root"))

process.accAnalyzer = cms.EDAnalyzer("GeneratorLevelAcceptanceAnalyzer")

process.p = cms.Path( process.accAnalyzer )

