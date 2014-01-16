import FWCore.ParameterSet.Config as cms
import os,sys

is7TeV=True
isMC=True
isTauEmbed=False
storeAllPF=False
gtag="START53_V23::All"

inputList=cms.untracked.vstring('/store/mc/Summer11dr53X/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV3-v1/10000/00A3A090-76EE-E211-9D7D-002618943935.root')
tfsOutputFile='Events.root'
outFile='edm_Events.root'

execfile( os.path.expandvars('${CMSSW_BASE}/src/UserCode/llvv_fwk/test/runDataAnalyzer_cfg.py'))

