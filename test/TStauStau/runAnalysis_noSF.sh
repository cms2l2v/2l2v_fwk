#!/bin/bash

runLocalAnalysisOverSamples.py -e runTStauStauAnalysisFWLite -j tstaustau_samples_full.json -d  /lustre/ncg.ingrid.pt/cmst3/store/user/cbeiraod/14_08_06_2l2nu_EDMtuples_merged/ -o /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/Results_noSF/ -c runAnalysis_cfg.py.templ -p "@saveSummaryTree=True @runSystematics=False @exclusiveRun=True @applyScaleFactors=False @doSVfit=False" -s 8nh
