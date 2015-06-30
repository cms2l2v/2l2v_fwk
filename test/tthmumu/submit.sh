#!/bin/bash

# test JSONFILE=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_test.json
JSONFILE=$CMSSW_BASE/src/UserCode/llvv_fwk/data/tthmumu/phys14_samples.json
#JSONFILE=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_ttbar.json
# test OUTDIR=$CMSSW_BASE/src/UserCode/llvv_fwk/test/chhiggs/results_ttbar3
OUTDIR=$CMSSW_BASE/src/UserCode/llvv_fwk/test/tthmumu/results

mkdir -p ${OUTDIR}

runLocalAnalysisOverSamples.py -e runTtHMuMuAnalysis -j ${JSONFILE} -o ${OUTDIR} -d  /dummy/ -c $CMSSW_BASE/src/UserCode/llvv_fwk/test/runAnalysis_cfg.py.templ -p "@useMVA=False @saveSummaryTree=False @runSystematics=False @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0" -s 8nh

#runLocalAnalysisOverSamples.py -e runChHiggsAnalysis -j $CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_samples.json -o $CMSSW_BASE/src/UserCode/llvv_fwk/test/chhiggs/results_nonexclusive -d  /dummy/ -c $CMSSW_BASE/src/UserCode/llvv_fwk/test/runAnalysis_cfg.py.templ -p "@useMVA=False @saveSummaryTree=False @runSystematics=False @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0" -s 8nh

#runLocalAnalysisOverSamples.py -e runChHiggsAnalysis -j $CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_ttbar.json -o $CMSSW_BASE/src/UserCode/llvv_fwk/test/chhiggs/results -d  /store/group/phys_higgs/cmshzz2l2v/2014_03_20/ -c $CMSSW_BASE/src/UserCode/llvv_fwk/test/runAnalysis_cfg.py.templ -p "@useMVA=False @saveSummaryTree=False @runSystematics=False @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0" -s 8nh
