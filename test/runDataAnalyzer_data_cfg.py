import os,sys

isMC=False
isTauEmbed=False
gtag="FT_53_V21_AN4::All"

from UserCode.2l2v_fwk.storeTools_cff import configureSourceFromCommandLine
outFile, inputList = configureSourceFromCommandLine()
tfsOutputFile=outFile
outFile=os.path.dirname(outFile)+'/edm_'+os.path.basename(outFile)

execfile( os.path.expandvars('${CMSSW_BASE}/src/UserCode/2l2v_fwk/test/runDataAnalyzer_cfg.py'))

