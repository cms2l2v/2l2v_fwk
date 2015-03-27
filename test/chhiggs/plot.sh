#!/bin/bash

#DIR=~/www/13TeV_topsel_fix/
DIR=~/www/13TeV_topsel_fix_nonexclusive/
mkdir -p ${DIR}
rm -rf ${DIR}*
cp ~/www/HIG-13-026/index.php ${DIR}

#JSONFILE=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_test.json
JSONFILEDILEPTON=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_plot_dileptons.json
JSONFILELEPTAU=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_plot_leptontau.json
INDIR=$CMSSW_BASE/src/UserCode/llvv_fwk/test/chhiggs/results_ttbar2/


#runPlotterFWLite --iEcm 13 --iLumi 19700 --inDir $CMSSW_BASE/src/UserCode/llvv_fwk/test/chhiggs/results/ --outDir ${DIR} --outFile ${DIR}plotter.root  --json $CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_plot.json --no2D --noPowers --plotExt .png --plotExt .pdf --only eventflow
#runPlotterFWLite --iEcm 13 --iLumi 19700 --inDir $CMSSW_BASE/src/UserCode/llvv_fwk/test/chhiggs/results/ --outDir ${DIR} --outFile ${DIR}plotter.root  --json $CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_plot.json --no2D --noPowers --plotExt .png --plotExt .pdf

# Dilepton
runFixedPlotter --iEcm 13 --iLumi 5822 --inDir ${INDIR} --outDir ${DIR} --outFile ${DIR}plotter.root  --json ${JSONFILEDILEPTON} --cutflow all_initNorm --forceMerge --no2D --noPowers --plotExt .png  --only ee_eventflow --only emu_eventflow --only mumu_eventflow

# Leptontau
runFixedPlotter --iEcm 13 --iLumi 5822 --inDir ${INDIR} --outDir ${DIR} --outFile ${DIR}plotter.root  --json ${JSONFILELEPTAU} --cutflow all_initNorm --forceMerge --no2D --noPowers --plotExt .png  --only singlee_eventflow --only singlemu_eventflow


# Lessen the burden on the web browser

#for CHAN in emu ee mumu singlemu singlee
#do
#    mkdir ${DIR}temp${CHAN}/
#    mv ${DIR}${CHAN}* ${DIR}temp${CHAN}/
#    mv ${DIR}temp${CHAN}/ ${DIR}${CHAN}/
#    cp ~/www/HIG-13-026/index.php ${DIR}${CHAN}/
#done