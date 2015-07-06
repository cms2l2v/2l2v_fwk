#rm $CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v/results/*.root
#in the command bellow, you can replace '8nh' by 'crab' to run using crab instead of LSF/Condor

#used to take additional arguments from the command line
arguments=''; for var in "$@"; do arguments=$arguments" "$var; done
if [ $# -ge 1 ]; then echo "Additional arguments will be considered: "$arguments ;fi 

runAnalysisOverSamples.py -e runHZZ2l2vAnalysis -j $CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v/samples.json -o $CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v/results -d  /store/group/phys_higgs/cmshzz2l2v/2014_03_20/ -c $CMSSW_BASE/src/UserCode/llvv_fwk/test/runAnalysis_cfg.py.templ -p "@useMVA=True @saveSummaryTree=True @runSystematics=True @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0" -s 8nh --report True $arguments
