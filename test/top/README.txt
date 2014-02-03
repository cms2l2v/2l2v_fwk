###
TOP
###

### l+jets
runLocalAnalysisOverSamples.py -e runLjetsAnalysis -j data/top_lj_samples.json  -d /store/cmst3/user/psilva/5311_ntuples -o ~/work/top_5311/lj/  -c test/runAnalysis_cfg.py.templ -p "@runSystematics=False @saveSummaryTree=False @weightsFile='data/weights/'" -s 8nh
runPlotter --iLumi 19701 --inDir ~/work/top_5311/lj/  --json data/top_lj_samples.json      --outFile ~/work/top_5311/plotter_lj.root    --noLogy --showUnc

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
runPlotter --iLumi 19701 --inDir ~/work/top_5311/nom/  --json data/top_samples.json   --outFile ~/work/top_5311/plotter_nom.root      --noLog --plotExt .pdf --showUnc
runPlotter --iLumi 19701 --inDir ~/work/top_5311/nom/  --json data/dataperiod_comp.json     --noRoot --noLog --plotExt .pdf --showUnc --outDir stability
fitDYforTop --in ~/work/top_5311/plotter_nom.root  --ttrep ~/work/top_5311/plotter_syst.root --smooth --out dyFit --syst

fitTTbarCrossSection --in ~/work/top_5311/plotter_nom.root --json data/top_samples.json --syst ~/work/top_5311/plotter_syst.root --bins 1,2,3,4 --out xsec/1inc > xsec/xsec1jinc_result.txt
fitTTbarCrossSection --in ~/work/top_5311/plotter_nom.root --json data/top_samples.json --syst ~/work/top_5311/plotter_syst.root --bins 2,3,4   --out xsec/2inc > xsec/xsec2jinc_result.txt
fitTTbarCrossSection --in ~/work/top_5311/plotter_nom.root --json data/top_samples.json --syst ~/work/top_5311/plotter_syst.root --bins 1       --out xsec/1jet > xsec/xsec1jexc_result.txt
fitTTbarCrossSection --in ~/work/top_5311/plotter_nom.root --json data/top_samples.json --syst ~/work/top_5311/plotter_syst.root --bins 2       --out xsec/2jet > xsec/xsec2jexc_result.txt
fitTTbarCrossSection --in ~/work/top_5311/plotter_nom.root --json data/top_samples.json --syst ~/work/top_5311/plotter_syst.root --bins 3       --out xsec/3jet > xsec/xsec3jexc_result.txt
fitTTbarCrossSection --in ~/work/top_5311/plotter_nom.root --json data/top_samples.json --syst ~/work/top_5311/plotter_syst.root --bins 4       --out xsec/4jet > xsec/xsec4jexc_result.txt

## TOP-12-035

### lepton-jet invariant mass analysis
fitMljSpectrum --in ~/work/top_5311/plotter_nom.root --syst ~/work/top_5311/plotter_syst.root --json data/top_samples.json --systJson data/top_syst_samples.json

### heavy flavor fits
for r in 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99 1.0; do
#for r in -1.0; do
for i in `seq 1 200`; do
    submit2batch.sh -q8nh -JHFC${i} runHFCClosureTests.sh ${r} ${i} ~/work/top_5311/hfc_closure_${r}/ 
done
done

fitHeavyFlavorContent --par test/top/hfcParams_2012_data_cfg.json --btag test/top/csvL_2012_data_cfg.json --in ~/work/top_5311/plotter_nom.root --fit 0
fitHeavyFlavorContent --par test/top/hfcParams_2012_data_cfg.json --btag test/top/csvL_2012_data_cfg.json --in ~/work/top_5311/plotter_nom.root --fit 3
fitHeavyFlavorContent --par test/top/hfcParams_2012_data_cfg.json --btag test/top/csvL_2012_data_cfg.json --in ~/work/top_5311/plotter_nom.root --fit 4

## BTV-13-001
for i in `seq 0 4`; do rundFtM --flav ~/www/BTV-13-001/dFtM/csv${i}/flavbreakup.json --btag ~/www/BTV-13-001/dFtM/csv${i}/csv${i}_eff.json --in ~/www/BTV-13-001/dFtM/csv${i}/csv${i}_btags.root --fit 0; done
a=(`ls ~/www/BTV-13-001/dFtM/csv*/*workspace.root`)
allWS=""
for i in ${a[@]}; do allWS="${allWS},${i}"; done 
rundFtM  --ws ${allWS}

## TOP-13-007

runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/top_samples.json      -d /store/cmst3/user/psilva/5311_ntuples -o ~/work/top_5311/data/ -c test/runAnalysis_cfg.py.templ -p "@runSystematics=False @saveSummaryTree=True @weightsFile='data/weights/top_dysf.root'" -s 8nh
runLocalAnalysisOverSamples.py -e runTopAnalysis -j data/top_syst_samples.json -d /store/cmst3/user/psilva/5311_ntuples -o ~/work/top_5311/syst/ -c test/runAnalysis_cfg.py.templ -p "@runSystematics=False @saveSummaryTree=True"  -s 8nh

runPlotter --iLumi 19701 --inDir ~/work/top_5311/syst/ --json data/top_syst_samples.json --outFile ~/work/top_5311/plotter_syst_ue.root --noPlot
runPlotter --iLumi 19701 --inDir ~/work/top_5311/data/ --json data/top_samples.json --outFile ~/work/top_5311/plotter_ue.root --noPlot

## TOP-12-030

runLocalAnalysisOverSamples.py -e runBTVdijetAnalysis -j test/top/qcd_samples.json -d /store/cmst3/user/psilva/5311_qcd_ntuples -o ~/work/top_5311/qcd/ -c test/runAnalysis_cfg.py.templ -p "@runSystematics=False @saveSummaryTree=False @weightsFile='data/weights/'"  -s 8nh
runPlotter --iLumi 19407 --inDir ~/work/top_5311/qcd/ --json test/top/qcd_samples.json --outFile ~/work/top_5311/plotter_qcd.root 

