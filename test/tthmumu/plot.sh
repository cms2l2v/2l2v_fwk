#!/bin/bash

DIR=~/www/13TeV_ttHMuMu/
mkdir -p ${DIR}
rm -rf ${DIR}*
cp ~/www/HIG-13-026/index.php ${DIR}

#JSONFILEDILEPTON=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_plot_dileptons.json
#JSONFILELEPTAU=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_plot_leptontau.json

# test LUMI=5822
# test JSONFILEDILEPTON=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_test.json
# test JSONFILELEPTAU=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_test.json
# test INDIR=$CMSSW_BASE/src/UserCode/llvv_fwk/test/chhiggs/results_ttbar3/
# test ONLYDILEPTON="--only ee_eventflow --only emu_eventflow --only mumu_eventflow"
# test ONLYLEPTAU="--only singlee_eventflow --only singlemu_eventflow"


LUMI=1000
JSONFILE=$CMSSW_BASE/src/UserCode/llvv_fwk/data/tthmumu/phys14_samples.json
INDIR=$CMSSW_BASE/src/UserCode/llvv_fwk/test/tthmumu/results/
ONLY=""

runFixedPlotter --iEcm 13 --iLumi ${LUMI} --inDir ${INDIR} --outDir ${DIR} --outFile ${DIR}plotter.root  --json ${JSONFILE} --cutflow all_initNorm --forceMerge --no2D --noPowers --plotExt .pdf --generatePseudoData ${ONLY}
runFixedPlotter --iEcm 13 --iLumi ${LUMI} --inDir ${INDIR} --outDir ${DIR} --outFile ${DIR}plotter.root  --json ${JSONFILE} --cutflow all_initNorm --forceMerge --no2D --noPowers --plotExt .png --generatePseudoData ${ONLY}

# Lessen the burden on the web browser

##for CHAN in emu ee mumu singlemu singlee
##do
##    mkdir ${DIR}temp${CHAN}/
##    mv ${DIR}${CHAN}* ${DIR}temp${CHAN}/
##    mv ${DIR}temp${CHAN}/ ${DIR}${CHAN}/
##    cp ~/www/HIG-13-026/index.php ${DIR}${CHAN}/
##done