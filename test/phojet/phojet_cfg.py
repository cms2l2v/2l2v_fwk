import FWCore.ParameterSet.Config as cms

process = cms.Process("PhoJet")

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('file:input.root'))
#import PSet

#print PSet.process.source.fileNames

runProcess = cms.PSet(
    #input = cms.untracked.vstring("file:input.root"),
    #input = PSet.process.source.fileNames, 
    input = process.source.fileNames, 
    outdir = cms.string("results"),
    debug = cms.bool(True),
    isMC = cms.bool(True),
    xsec = cms.double(9999.99),
    mctruthmode = cms.int32(22)
)
