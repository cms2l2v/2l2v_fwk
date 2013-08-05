###
TOP
###

### run base
runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/top_samples.json      -d /store/cmst3/user/psilva/5311_ntuples -o ~/work/top_539/raw/  -c test/runAnalysis_cfg.py.templ -p "@runSystematics=True @saveSummaryTree=False" -s 8nh
runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/top_syst_samples.json -d /store/cmst3/user/psilva/5311_ntuples -o ~/work/top_539/syst/ -c test/runAnalysis_cfg.py.templ -p "@runSystematics=False @saveSummaryTree=True"  -s 8nh

### fit dy and re-run with new scale
runPlotter --iLumi 19683 --inDir ~/work/top_539/raw/  --json data/top_samples.json      --outFile ~/work/top_539/plotter_dy_raw.root        --noPlot --only mtsum --only dilarc
runPlotter --iLumi 19683 --inDir ~/work/top_539/syst/ --json data/top_syst_samples.json --outFile ~/work/top_539/plotter_syst_fordyfit.root --noPlot --only mtsum --only dilarc
fitDYforTop --in ~/work/top_539/plotter_dy_raw.root  --ttrep ~/work/top_539/plotter_syst_fordyfit.root --smooth --out dyFit --syst
mv top_dysf.root data/weights/
runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/top_samples.json      -d /store/cmst3/user/psilva/Summer13_ntuples -o ~/work/top_539/raw/  -c test/runAnalysis_cfg.py.templ -p "@runSystematics=True @saveSummaryTree=False @weightsFile='data/weights/top_dysf.root'" -s 8nh -t DY


### fit ttbar signal strength and re-run 
runPlotter --iLumi 19683 --inDir ~/work/top_539/raw/  --json data/top_samples.json      --outFile ~/work/top_539/plotter_forxsec.root      --noPlot --only finalevt
runPlotter --iLumi 19683 --inDir ~/work/top_539/syst/ --json data/top_syst_samples.json --outFile ~/work/top_539/plotter_syst_forxsec.root --noPlot --only finalevt

fitTTbarCrossSection --in ~/work/top_539/plotter_forxsec.root --json data/top_samples.json --syst ~/work/top_539/plotter_syst_forxsec.root --bins 2,3,4 --out xsec      > xsec_result.txt
fitTTbarCrossSection --in ~/work/top_539/plotter_forxsec.root --json data/top_samples.json --syst ~/work/top_539/plotter_syst_forxsec.root --bins 1     --out xsec/1jet > xsec1j_result.txt
fitTTbarCrossSection --in ~/work/top_539/plotter_forxsec.root --json data/top_samples.json --syst ~/work/top_539/plotter_syst_forxsec.root --bins 2     --out xsec/2jet > xsec2j_result.txt
fitTTbarCrossSection --in ~/work/top_539/plotter_forxsec.root --json data/top_samples.json --syst ~/work/top_539/plotter_syst_forxsec.root --bins 3     --out xsec/3jet > xsec3j_result.txt
fitTTbarCrossSection --in ~/work/top_539/plotter_forxsec.root --json data/top_samples.json --syst ~/work/top_539/plotter_syst_forxsec.root --bins 4     --out xsec/4jet > xsec4j_result.txt



runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/top_samples.json      -d /store/cmst3/user/psilva/Summer13_ntuples  -o ~/work/top_539/final -c test/runAnalysis_cfg.py.templ -p "@runSystematics=False @saveSummaryTree=True @weightsFile='data/weights/top_dysf.root'" -s 8nh

runPlotter --iLumi 19683 --inDir ~/work/top_539/final/ --json data/top_samples.json --outFile ~/work/top_539/plotter.root --noLog --showUnc
runPlotter --iLumi 19683 --inDir ~/work/top_539/syst/ --json data/top_syst_samples.json --outFile ~/work/top_539/plotter_syst.root --noPlot

