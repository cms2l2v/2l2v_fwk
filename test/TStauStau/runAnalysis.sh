#!/bin/bash

# Inclusive
runLocalAnalysisOverSamples.py -e runTStauStauAnalysisFWLite -j tstaustau_samples_full.json -d /store/group/phys_higgs/cmshzz2l2v/2014_04_20/ -o Processed_Inclusive -c runAnalysis_cfg.py.templ -p "@useMVA=True @saveSummaryTree=False @runSystematics=False @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0 @exclusiveRun=False" -s 8nh
# Exclusive
#runLocalAnalysisOverSamples.py -e runTStauStauAnalysisFWLite -j tstaustau_samples_full.json -d /store/group/phys_higgs/cmshzz2l2v/2014_04_20/ -o Processed_Exclusive -c runAnalysis_cfg.py.templ -p "@useMVA=True @saveSummaryTree=False @runSystematics=False @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0 @exclusiveRun=True" -s 8nh
