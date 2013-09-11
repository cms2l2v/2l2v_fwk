runLocalAnalysisOverSamples.py -e runZHTauTauAnalysisFWLite -j $CMSSW_BASE/src/UserCode/llvv_fwk/data/zhtautau_samples.json -o $CMSSW_BASE/src/UserCode/llvv_fwk/test/zhtautau/results -d  /storage/data/cms/users/quertenmont/Higgs/Ntuples/13_08_30/ -c $CMSSW_BASE/src/UserCode/llvv_fwk/test/runAnalysis_cfg.py.templ -p "@useMVA=True @saveSummaryTree=True @runSystematics=True @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0" -s 8nh
#runPlotterFWLite --iEcm 8 --iLumi 19577 --inDir $CMSSW_BASE/src/UserCode/llvv_fwk/test/zhtautau/results/ --outDir $CMSSW_BASE/src/UserCode/llvv_fwk/test/zhtautau/plots/ --outFile $CMSSW_BASE/src/UserCode/llvv_fwk/test/zhtautau/plotter.root  --json $CMSSW_BASE/src/UserCode/llvv_fwk/data/zhtautau_samples.json


