import os,sys

isMC=True
isTauEmbed=False
storeAllPF=False
gtag="START53_V23::All"

from UserCode.llvv_fwk.storeTools_cff import configureSourceFromCommandLine
outFile, inputListArray = configureSourceFromCommandLine()
inputList=cms.untracked.vstring(inputListArray)
inputList=cms.untracked.vstring('/store/caf/user/tjkim/mc/Summer12_DR53X/TTJets_FullLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v2/00000/000C5D15-AB1A-E211-8BDE-00215E22053A.root',
                                '/store/caf/user/tjkim/mc/Summer12_DR53X/TTJets_FullLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v2/00000/0010005A-421A-E211-9E3C-E41F13181DA4.root')
tfsOutputFile=outFile
outFile=os.path.dirname(outFile)+'/edm_'+os.path.basename(outFile)

execfile( os.path.expandvars('${CMSSW_BASE}/src/UserCode/llvv_fwk/test/runDataAnalyzer_cfg.py'))

