###
PDF variations
###
runLocalAnalysisOverSamples.py -e computePDFvariations -j data/top_samples.json -d /store/cmst3/user/psilva/539_ntuples -o ~/work/top_539/ -c test/runAnalysis_cfg.py.templ -p "@runSystematics=False @saveSummaryTree=True" -t MC -s 8nh

####
EWK Z+2jets
####
runLocalAnalysisOverSamples.py -e runVBFZAnalysis -j data/vbfz_samples.json -d /store/cmst3/user/psilva/539_ntuples -o ~/work/ewkzp2j_539 -c test/runAnalysis_cfg.py.templ -p "@runSystematics=False @useMVA=True" -s 8nh
runPlotter --iLumi 19683 --inDir ~/work/ewkzp2j_539/ --json data/vbfz_samples.json --outFile ~/work/ewkzp2j_539/plotter.root 

runLocalAnalysisOverSamples.py -e runVBFZAnalysis -j data/vbfz_photon_samples.json -d /store/cmst3/user/psilva/539_ntuples -o ~/work/ewkzp2j_539 -c test/runAnalysis_cfg.py.templ -p "@runSystematics=False @useMVA=True" -s 8nh
runPlotter --iLumi 19683 --inDir ~/work/ewkzp2j_539/ --json data/vbfz_photon_samples.json --outFile ~/work/ewkzp2j_539/plotter_g_raw.root 





