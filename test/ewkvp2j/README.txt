####
CLOSURE TESTS
####
sh test/ewkvp2j/steerClosureTest.sh 0
sh test/ewkvp2j/steerClosureTest.sh 1
sh test/ewkvp2j/steerClosureTest.sh 2
sh test/ewkvp2j/steerClosureTest.sh 3
sh test/ewkvp2j/steerClosureTest.sh 4

####
EWK Z+2jets ANALYSIS
####
runLocalAnalysisOverSamples.py -e runVBFZAnalysis -j data/vbfz_samples.json -d /store/cmst3/user/psilva/539_ntuples -o ~/work/ewkzp2j_539/ll/ -c test/runAnalysis_cfg.py.templ -p "@runSystematics=False @useMVA=True" -s 8nh
runLocalAnalysisOverSamples.py -e runVBFZAnalysis -j data/vbfz_syst_samples.json -d /store/cmst3/user/psilva/539_ntuples -o ~/work/ewkzp2j_539/ll/ -c test/runAnalysis_cfg.py.templ -p "@runSystematics=False @useMVA=True" -s 8nh -t jj
runPlotter --iLumi 19736 --inDir ~/work/ewkzp2j_539/ll/ --json data/vbfz_samples.json --outFile ~/work/ewkzp2j_539/plotter.root 
runPlotter --iLumi 19736 --inDir ~/work/ewkzp2j_539/ll/ --json data/vbfz_syst_samples.json --outFile ~/work/ewkzp2j_539/plotter_syst.root 

runLocalAnalysisOverSamples.py -e runVBFZAnalysis -j data/vbfz_photon_samples.json -d /store/cmst3/user/psilva/539_ntuples -o ~/work/ewkzp2j_539/g/data/ -c test/runAnalysis_cfg.py.templ -p "@runSystematics=False @useMVA=True" -s 8nh
runPlotter --iLumi 19736 --inDir ~/work/ewkzp2j_539/g/data/ --json data/vbfz_photon_samples.json --outFile ~/work/ewkzp2j_539/plotter_g_raw.root 


####
BOXED ANALYSIS
####
runLocalAnalysisOverSamples.py -e  runEWKVjjAnalysis -j data/ewkvjj_samples.json -d /store/cmst3/user/psilva/5311_ntuples -o ~/work/ewkzp2j_5311/ -c test/runAnalysis_cfg.py.templ -p "@runSystematics=False @useMVA=False" -s 8nh
runPlotter --iLumi 19736 --inDir ~/work/ewkzp2j_5311/ --json data/ewkvjj_samples.json --outFile ~/work/ewkzp2j_5311/plotter.root

