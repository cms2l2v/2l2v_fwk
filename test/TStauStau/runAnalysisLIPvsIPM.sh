#!/bin/bash

# Exclusive without quick fix
## LIP
runLocalAnalysisOverSamples.py -e runTStauStauAnalysisFWLite -j tstaustau_samples_full.json -d /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/mergedTuples/ -o /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/LIP_selection/ -c runAnalysis_cfg.py.templ -p "@useMVA=True @saveSummaryTree=True @runSystematics=False @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0 @exclusiveRun=True @selection=LIP @doQuickFix=False" -s 8nh
## IPM
runLocalAnalysisOverSamples.py -e runTStauStauAnalysisFWLite -j tstaustau_samples_full.json -d /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/mergedTuples/ -o /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/IPM_selection/ -c runAnalysis_cfg.py.templ -p "@useMVA=True @saveSummaryTree=True @runSystematics=False @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0 @exclusiveRun=True @selection=IPM @doQuickFix=False" -s 8nh

# Exclusive with quick fix
## LIP
runLocalAnalysisOverSamples.py -e runTStauStauAnalysisFWLite -j tstaustau_samples_full.json -d /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/mergedTuples/ -o /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/LIP_selection_QF/ -c runAnalysis_cfg.py.templ -p "@useMVA=True @saveSummaryTree=True @runSystematics=False @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0 @exclusiveRun=True @selection=LIP @doQuickFix=True" -s 8nh
## IPM
runLocalAnalysisOverSamples.py -e runTStauStauAnalysisFWLite -j tstaustau_samples_full.json -d /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/mergedTuples/ -o /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/IPM_selection_QF/ -c runAnalysis_cfg.py.templ -p "@useMVA=True @saveSummaryTree=True @runSystematics=False @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0 @exclusiveRun=True @selection=IPM @doQuickFix=True" -s 8nh
