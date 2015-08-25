#rm $CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v/results/*.root
#in the command bellow, you can replace '8nh' by 'crab' to run using crab instead of LSF/Condor
runAnalysisOverSamples.py -e runPhoJetAnalysis -j $CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v/phojet/samples.json -o $CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v/phojet/results -c $CMSSW_BASE/src/UserCode/llvv_fwk/test/runAnalysis_cfg.py.templ -p "@useMVA=False @trig=False @saveSummaryTree=False @runSystematics=False @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0" -s 8nh --report True
