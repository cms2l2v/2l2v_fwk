import os,sys

isMC=True
isTauEmbed=False
storeAllPF=False
gtag="START53_V23::All"

from UserCode.llvv_fwk.storeTools_cff import configureSourceFromCommandLine
outFile, inputList = configureSourceFromCommandLine()
tfsOutputFile=outFile
outFile=os.path.dirname(outFile)+'/edm_'+os.path.basename(outFile)

execfile( os.path.expandvars('${CMSSW_BASE}/src/UserCode/llvv_fwk/test/runDataAnalyzer_cfg.py'))

