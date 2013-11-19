import FWCore.ParameterSet.Config as cms

ewkzp2jBase = cms.PSet( weightsDir = cms.string("/afs/cern.ch/user/p/psilva/work/ewkzp2j_5311/ll/base_weights"),
                        methodList = cms.vstring('Fisher','BDTD','MLP'),
                        varsList   = cms.vstring("mjj","detajj","spt")
                        )

ewkzp2jFull = cms.PSet( weightsDir = cms.string("/afs/cern.ch/user/p/psilva/work/ewkzp2j_5311/ll/full_weights"),
                        methodList = cms.vstring('Fisher','BDTD', 'MLP'),
                        varsList   = cms.vstring("mjj","detajj","setajj","eta1","eta2","pt1","pt2","spt","qg1","qg2")
                        )

ewkzp2jFullNoQG = cms.PSet( weightsDir = cms.string("/afs/cern.ch/user/p/psilva/work/ewkzp2j_5311/ll/weights"),
                           methodList = cms.vstring('Fisher','BDTD', 'MLP'),
                           varsList   = cms.vstring("mjj","detajj","setajj","eta1","eta2","pt1","pt2","spt")                           
                           )
