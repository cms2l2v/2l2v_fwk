import FWCore.ParameterSet.Config as cms

ewkzp2j = cms.PSet(  weightsDir = cms.string("${CMSSW_BASE}/src/UserCode/llvv_fwk/data/weights"),
                     methodList = cms.vstring('Fisher','LikelihoodD','BDTD', 'FisherCat'),
                     varsList   = cms.vstring("mjj","detajj","spt")
                     )
