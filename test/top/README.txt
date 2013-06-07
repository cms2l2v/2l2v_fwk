###
TOP
###
runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/top_samples.json -d /store/cmst3/user/psilva/539_ntuples -o ~/work/top_539/ -c test/runAnalysis_cfg.py.templ -p "@runSystematics=False @saveSummaryTree=True" -s 8nh
runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/top_syst_samples.json -d /store/cmst3/user/psilva/539_ntuples -o ~/work/top_539/ -c test/runAnalysis_cfg.py.templ -p "@runSystematics=False @saveSummaryTree=True" -s 8nh
runPlotter --iLumi 19683 --inDir ~/work/top_539/ --json data/top_samples.json --outFile ~/work/top_539/plotter_raw.root --noPlot
runPlotter --iLumi 19683 --inDir ~/work/top_539/ --json data/top_syst_samples.json --outFile ~/work/top_539/plotter_syst.root --noPlot

runPlotter --iLumi 19683 --inDir ~/work/top_539/ --outDir ~/work/top_539/plots --json data/top_samples.json --outFile ~/work/top_539/plotter.root --noLog  --showUnc  --plotExt .pdf

