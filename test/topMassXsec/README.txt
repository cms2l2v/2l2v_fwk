###
Charged Higgs Analysis
###

Mailto: pietro.vischia at gmail.com

Analyze
-------
cd $CMSSW_BASE/src/UserCode/llvv_fwk/


### Run analysis
runLocalAnalysisOverSamples.py -e runTopMassAnalysis -j $CMSSW_BASE/src/UserCode/llvv_fwk/data/zhtautau_samples.json -o /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_fwlite/ -d  /store/group/phys_higgs/cmshzz2l2v/2013_08_30/ -c $CMSSW_BASE/src/UserCode/llvv_fwk/test/runAnalysis_cfg.py.templ -p "@useMVA=True @saveSummaryTree=True @runSystematics=True @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0" -s 8nh


Plots & tables
--------------


### FWLite test
runPlotterFWLite --iEcm 8 --iLumi 19577 --inDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_fwlite/ --outDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_fwlite/plots/ --outFile /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_fwlite/plotter.root  --json $CMSSW_BASE/src/UserCode/llvv_fwk/data/zhtautau_samples.json


Tables only
-----------

Datacards
---------
