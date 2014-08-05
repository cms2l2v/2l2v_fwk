#!/bin/bash

# Inclusive
#runLocalAnalysisOverSamples.py -e runTStauStauAnalysisFWLite -j boosted_quick_fix.json -d /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/mergedTuples/ -o /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/QuickFix/ -c runAnalysis_cfg.py.templ -p "@useMVA=True @saveSummaryTree=True @runSystematics=False @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0 @exclusiveRun=False" -s 8nh
# Exclusive
runLocalAnalysisOverSamples.py -e runTStauStauAnalysisFWLite -j boosted_quick_fix.json -d /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/mergedTuples/ -o /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/QuickFix/ -c runAnalysis_cfg.py.templ -p "@useMVA=True @saveSummaryTree=True @runSystematics=False @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0 @exclusiveRun=True" -s 8nh
