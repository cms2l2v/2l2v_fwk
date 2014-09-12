#!/bin/bash

# Exclusive without quick fix
## IPM with LIP Tau ID
runLocalAnalysisOverSamples.py -e runTStauStauAnalysisFWLite -j tstaustau_samples_full.json -d /lustre/ncg.ingrid.pt/cmst3/store/user/cbeiraod/14_08_06_2l2nu_EDMtuples_merged/ -o /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/IPM_selection_LIPTauID/ -c runAnalysis_cfg.py.templ -p "@useMVA=True @saveSummaryTree=True @runSystematics=False @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0 @exclusiveRun=True @selection=IPM @doLIPTauID=True @doLIPBVeto=False @periodHLT=2012B" -s 8nh
## IPM with LIP B-Veto
runLocalAnalysisOverSamples.py -e runTStauStauAnalysisFWLite -j tstaustau_samples_full.json -d /lustre/ncg.ingrid.pt/cmst3/store/user/cbeiraod/14_08_06_2l2nu_EDMtuples_merged/ -o /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/IPM_selection_LIPBVeto/ -c runAnalysis_cfg.py.templ -p "@useMVA=True @saveSummaryTree=True @runSystematics=False @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0 @exclusiveRun=True @selection=IPM @doLIPTauID=False @doLIPBVeto=True @periodHLT=2012B" -s 8nh
