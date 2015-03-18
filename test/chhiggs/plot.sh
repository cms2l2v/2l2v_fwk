#!/bin/bash

#DIR=~/www/13TeV_topsel_fix/
DIR=~/www/13TeV_topsel_fix_nonexclusive2/
mkdir -p ${DIR}
rm -rf ${DIR}*
cp ~/www/HIG-13-026/index.php ${DIR}

#runPlotterFWLite --iEcm 13 --iLumi 19700 --inDir $CMSSW_BASE/src/UserCode/llvv_fwk/test/chhiggs/results/ --outDir ${DIR} --outFile ${DIR}plotter.root  --json $CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_plot.json --no2D --noPowers --plotExt .png --plotExt .pdf --only eventflow
#runPlotterFWLite --iEcm 13 --iLumi 19700 --inDir $CMSSW_BASE/src/UserCode/llvv_fwk/test/chhiggs/results/ --outDir ${DIR} --outFile ${DIR}plotter.root  --json $CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_plot.json --no2D --noPowers --plotExt .png --plotExt .pdf

runFixedPlotter --iEcm 13 --iLumi 19700 --inDir $CMSSW_BASE/src/UserCode/llvv_fwk/test/chhiggs/results_ttbar/ --outDir ${DIR} --outFile ${DIR}plotter.root  --json $CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_ttbar.json --cutflow all_initNorm --no2D --noPowers --plotExt .png  --only eventflow 


# Lessen the burden on the web browser

#for CHAN in emu ee mumu singlemu singlee
#do
#    mkdir ${DIR}temp${CHAN}/
#    mv ${DIR}${CHAN}* ${DIR}temp${CHAN}/
#    mv ${DIR}temp${CHAN}/ ${DIR}${CHAN}/
#    cp ~/www/HIG-13-026/index.php ${DIR}${CHAN}/
#done