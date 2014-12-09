import FWCore.ParameterSet.Config as cms

runProcess = cms.PSet(
    input = cms.vstring("file:input.root"),
    outdir = cms.string("results"),
    debug = cms.bool(True),
    isMC = cms.bool(True),
    xsec = cms.double(9999.99),
    mctruthmode = cms.int32(22)
)
