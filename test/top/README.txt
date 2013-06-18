###
TOP
###

### run base
runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/top_samples.json -d /store/cmst3/user/psilva/Summer13_ntuples -o ~/work/top_539/ -c test/runAnalysis_cfg.py.templ -p "@runSystematics=False @saveSummaryTree=True" -s 8nh
runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/top_syst_samples.json -d /store/cmst3/user/psilva/Summer13_ntuples -o ~/work/top_539/ -c test/runAnalysis_cfg.py.templ -p "@runSystematics=False @saveSummaryTree=True" -s 8nh

### fit dy
runPlotter --iLumi 19683 --inDir ~/work/top_539/ --json data/top_samples.json --outFile ~/work/top_539/plotter_dy_raw.root --noPlot --only mtsum --only dilarc
runPlotter --iLumi 19683 --inDir ~/work/top_539/ --json data/top_syst_samples.json --outFile ~/work/top_539/plotter_syst_fordyfit.root --noPlot --only mtsum --only dilarc
fitDYforTop --in ~/work/top_539/plotter_dy_raw.root  --ttrep ~/work/top_539/plotter_syst_fordyfit.root --smooth --out dyFit
mv top_dysf.root data/weights/

### re-run re-scaling DY
runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/top_samples.json -d /store/cmst3/user/psilva/Summer13_ntuples -o ~/work/top_539/ -c test/runAnalysis_cfg.py.templ -p "@runSystematics=False @saveSummaryTree=True @weightsFile='data/weights/top_dysf.root'" -t DY -s 8nh

runPlotter --iLumi 19683 --inDir ~/work/top_539/ --json data/top_samples.json --outFile ~/work/top_539/plotter.root --noPlot
runPlotter --iLumi 19683 --inDir ~/work/top_539/ --json data/top_syst_samples.json --outFile ~/work/top_539/plotter_syst.root --noPlot

