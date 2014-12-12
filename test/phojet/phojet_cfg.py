import FWCore.ParameterSet.Config as cms

process = cms.Process("PhoJet")

inputFiles = cms.untracked.vstring("file:input.root")
try:
    import PSet
    print "loading from PSet..."
    inputFiles =  PSet.process.source.fileNames
except:
    print "not able to import" 
    pass

print inputFiles

runProcess = cms.PSet(
    #input = cms.untracked.vstring("file:input.root"),
    input = inputFiles, 
    #outdir = cms.string("results"),
    output = cms.string("output.root"),
    debug = cms.bool(True),
    isMC = cms.bool(True),
    xsec = cms.double(9999.99),
    mctruthmode = cms.int32(22)
)
