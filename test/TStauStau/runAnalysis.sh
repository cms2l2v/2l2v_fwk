#!/bin/bash

## New stuff, Single lepton datasets aren't finished
runLocalAnalysisOverSamples.py -e runTStauStauAnalysisFWLite -j tstaustau_samples_full.json -d /store/group/phys_higgs/cmshzz2l2v/2014_04_20/ -o Processed -c runAnalysis_cfg.py.templ -p "@useMVA=True @saveSummaryTree=False @runSystematics=False @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0" -s 8nh
