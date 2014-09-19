#!/bin/bash

# Exclusive without quick fix
## LIP with 2012D HLT and lower lepton pT thresholds
runLocalAnalysisOverSamples.py -e runTStauStauAnalysisFWLite -j tstaustau_samples_full.json -d /lustre/ncg.ingrid.pt/cmst3/store/user/cbeiraod/14_08_06_2l2nu_EDMtuples_merged/ -o /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/LIP_2012D_selection/ -c runAnalysis_cfg.py.templ -p "@useMVA=True @saveSummaryTree=True @runSystematics=False @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0 @exclusiveRun=True @selection=LIP @doLIPTauID=True @doLIPBVeto=True @periodHLT=2012D @doCombined=False" -s 8nh
