###
TOP
###

### run base
runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/top_samples.json      -d /store/cmst3/user/psilva/5311_ntuples -o ~/work/top_5311/nom/  -c test/runAnalysis_cfg.py.templ -p "@runSystematics=False @saveSummaryTree=False @weightsFile='data/weights/'" -s 8nh
runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/top_syst_samples.json -d /store/cmst3/user/psilva/5311_ntuples -o ~/work/top_5311/syst/ -c test/runAnalysis_cfg.py.templ -p "@runSystematics=False @saveSummaryTree=True @weightsFile='data/weights/'"  -s 8nh

### fit dy and re-run with new scale
runPlotter --iLumi 19701 --inDir ~/work/top_5311/nom/  --json data/top_samples.json      --outFile ~/work/top_5311/plotter_dy_nom.root        --noPlot --only mtsum --only dilarc
runPlotter --iLumi 19701 --inDir ~/work/top_5311/syst/ --json data/top_syst_samples.json --outFile ~/work/top_5311/plotter_syst.root          --noPlot
fitDYforTop --in ~/work/top_5311/plotter_dy_nom.root  --ttrep ~/work/top_5311/plotter_syst.root --smooth --out dyFit 
mv top_dysf.root data/weights/

### re-run final selection for cross section measurement
runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/top_samples.json      -d /store/cmst3/user/psilva/5311_ntuples -o ~/work/top_5311/nom/  -c test/runAnalysis_cfg.py.templ -p "@runSystematics=True @saveSummaryTree=False @weightsFile='data/weights/'" -s 8nh
runPlotter --iLumi 19701 --inDir ~/work/top_5311/nom/  --json data/top_samples.json      --outFile ~/work/top_5311/plotter_nom.root      --noLog --plotExt .pdf --showUnc
fitDYforTop --in ~/work/top_5311/plotter_nom.root  --ttrep ~/work/top_5311/plotter_syst.root --smooth --out dyFit --syst

fitTTbarCrossSection --in ~/work/top_5311/plotter_nom.root --json data/top_samples.json --syst ~/work/top_5311/plotter_syst.root --bins 1,2,3,4 --out xsec/1inc > xsec/xsec1jinc_result.txt
fitTTbarCrossSection --in ~/work/top_5311/plotter_nom.root --json data/top_samples.json --syst ~/work/top_5311/plotter_syst.root --bins 2,3,4   --out xsec/2inc > xsec/xsec2jinc_result.txt
fitTTbarCrossSection --in ~/work/top_5311/plotter_nom.root --json data/top_samples.json --syst ~/work/top_5311/plotter_syst.root --bins 1       --out xsec/1jet > xsec/xsec1jexc_result.txt
fitTTbarCrossSection --in ~/work/top_5311/plotter_nom.root --json data/top_samples.json --syst ~/work/top_5311/plotter_syst.root --bins 2       --out xsec/2jet > xsec/xsec2jexc_result.txt
fitTTbarCrossSection --in ~/work/top_5311/plotter_nom.root --json data/top_samples.json --syst ~/work/top_5311/plotter_syst.root --bins 3       --out xsec/3jet > xsec/xsec3jexc_result.txt
fitTTbarCrossSection --in ~/work/top_5311/plotter_nom.root --json data/top_samples.json --syst ~/work/top_5311/plotter_syst.root --bins 4       --out xsec/4jet > xsec/xsec4jexc_result.txt

### lepton-jet invariant mass analysis
fitMljSpectrum --in ~/work/top_5311/plotter_nom.root --syst ~/work/top_5311/plotter_syst.root --json data/top_samples.json --systJson data/top_syst_samples.json

runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/top_samples.json      -d /store/cmst3/user/psilva/Summer13_ntuples  -o ~/work/top_5311/final -c test/runAnalysis_cfg.py.templ -p "@runSystematics=False @saveSummaryTree=True @weightsFile='data/weights/top_dysf.root'" -s 8nh

runPlotter --iLumi 19683 --inDir ~/work/top_5311/final/ --json data/top_samples.json --outFile ~/work/top_5311/plotter.root --noLog --showUnc
runPlotter --iLumi 19683 --inDir ~/work/top_5311/syst/ --json data/top_syst_samples.json --outFile ~/work/top_5311/plotter_syst.root --noPlot


### UE studies

runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/top_samples.json      -d /store/cmst3/user/psilva/5311_ntuples -o ~/work/top_5311/data/ -c test/runAnalysis_cfg.py.templ -p "@runSystematics=False @saveSummaryTree=True @weightsFile='data/weights/top_dysf.root'" -s 8nh
runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/top_syst_samples.json -d /store/cmst3/user/psilva/5311_ntuples -o ~/work/top_5311/syst/ -c test/runAnalysis_cfg.py.templ -p "@runSystematics=False @saveSummaryTree=True"  -s 8nh

runPlotter --iLumi 19701 --inDir ~/work/top_5311/syst/ --json data/top_syst_samples.json --outFile ~/work/top_5311/plotter_syst_ue.root --noPlot
runPlotter --iLumi 19701 --inDir ~/work/top_5311/data/ --json data/top_samples.json --outFile ~/work/top_5311/plotter_ue.root --noPlot