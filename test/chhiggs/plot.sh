mkdir -p ~/www/13TeV_topsel_nopu
mkdir -p ~/www/13TeV_topsel_nopu_noqcd
cp ~/www/HIG-13-026/index.php ~/www/13TeV_topsel_nopu
cp ~/www/HIG-13-026/index.php ~/www/13TeV_topsel_nopu_noqcd

runPlotter --iEcm 13 --iLumi 1000 --inDir $CMSSW_BASE/src/UserCode/llvv_fwk/test/chhiggs/results/ --outDir ~/www/13TeV_topsel_nopu --outFile ~/www/13TeV_topsel_nopu/plotter.root  --json $CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_samples.json --no2D --plotExt .png --plotExt .pdf

runPlotter --iEcm 13 --iLumi 1000 --inDir $CMSSW_BASE/src/UserCode/llvv_fwk/test/chhiggs/results/ --outDir ~/www/13TeV_topsel_nopu_noqcd --outFile ~/www/13TeV_topsel_nopu_noqcd/plotter.root  --json $CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_plot.json --no2D --plotExt .png --plotExt .pdf
