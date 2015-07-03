#!/bin/bash

#DIR=~/www/13TeV_topsel_fix/
#DIR=~/www/13TeV_topsel_fix_ala8/
DIR=~/www/13TeV_topsel_fix_ala9/
DIR=~/www/13TeV_topsel_fix_talk/
DIR=~/www/13TeV_topsel_fix_posttalk/
#DIR=~/www/13TeV_forinclusivecheck/
mkdir -p ${DIR}
rm -rf ${DIR}*
cp ~/www/HIG-13-026/index.php ${DIR}


# test LUMI=5822
# test JSONFILEDILEPTON=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_test.json
# test JSONFILELEPTAU=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_test.json
# test INDIR=$CMSSW_BASE/src/UserCode/llvv_fwk/test/chhiggs/results_ttbar3/
# test ONLYDILEPTON="--only ee_eventflow --only emu_eventflow --only mumu_eventflow"
# test ONLYLEPTAU="--only singlee_eventflow --only singlemu_eventflow"


LUMI=1000
JSONFILEDILEPTON=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_plot_dileptons.json
JSONFILELEPTAU=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_plot_leptontau.json
#JSONFILEDILEPTON=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_samples.json
#JSONFILELEPTAU=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_samples.json
#INDIR=$CMSSW_BASE/src/UserCode/llvv_fwk/test/chhiggs/results_ttbar4/
INDIR=$CMSSW_BASE/src/UserCode/llvv_fwk/test/chhiggs/results_ttbar/
#ONLYDILEPTON=" --onlyStartWith emu --onlyStartWith ee --onlyStartWith mumu "
#ONLYLEPTAU=" --onlyStartWith singlemu --onlyStartWith singlee "
ONLYDILEPTON=" --onlyStartWith emu_eventflow --onlyStartWith ee_eventflow --onlyStartWith mumu_eventflow --onlyStartWith emu_step1leadpt --onlyStartWith emu_step1leadeta  --onlyStartWith emu_step3leadjetpt  --onlyStartWith emu_step3leadjeteta  --onlyStartWith emu_step3met  --onlyStartWith emu_step3nbjets --onlyStartWith ee_step1leadpt --onlyStartWith ee_step1leadeta  --onlyStartWith ee_step3leadjetpt  --onlyStartWith ee_step3leadjeteta  --onlyStartWith ee_step3met  --onlyStartWith ee_step3nbjets --onlyStartWith mumu_step1leadpt --onlyStartWith mumu_step1leadeta  --onlyStartWith mumu_step3leadjetpt  --onlyStartWith mumu_step3leadjeteta  --onlyStartWith mumu_step3met  --onlyStartWith mumu_step3nbjets"

ONLYLEPTAU=" --onlyStartWith singlemu_step6met --onlyStartWith singlee_step6met  --onlyStartWith singlemu_step6tauleadpt --onlyStartWith singlee_step6tauleadpt --onlyStartWith singlemu_final --onlyStartWith singlee_final --onlyStartWith singlee_eventflowslep --onlyStartWith singlemu_eventflowslep  --onlyStartWith singlemu_step1leadpt  --onlyStartWith singlemu_step1leadeta --onlyStartWith singlemu_step2leadjetpt --onlyStartWith singlemu_step2leadjeteta --onlyStartWith singlemu_step3met --onlyStartWith singlemu_step3nbjets --onlyStartWith singlee_step1leadpt  --onlyStartWith singlee_step1leadeta --onlyStartWith singlee_step2leadjetpt --onlyStartWith singlee_step2leadjeteta --onlyStartWith singlee_step3met --onlyStartWith singlee_step3nbjets"

#ONLYDILEPTON=" --onlyContain emu_eventflow --onlyContain ee_eventflow --onlyContain mumu_eventflow "
#ONLYLEPTAU=" --onlyContain singlemu_eventflowslep --onlyContain singlee_eventflowslep "

#runPlotterFWLite --iEcm 13 --iLumi 19700 --inDir $CMSSW_BASE/src/UserCode/llvv_fwk/test/chhiggs/results/ --outDir ${DIR} --outFile ${DIR}plotter.root  --json $CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_plot.json --no2D --noPowers --plotExt .png --plotExt .pdf --only eventflow
#runPlotterFWLite --iEcm 13 --iLumi 19700 --inDir $CMSSW_BASE/src/UserCode/llvv_fwk/test/chhiggs/results/ --outDir ${DIR} --outFile ${DIR}plotter.root  --json $CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_plot.json --no2D --noPowers --plotExt .png --plotExt .pdf

# Dilepton
runFixedPlotter --iEcm 13 --iLumi ${LUMI} --inDir ${INDIR} --outDir ${DIR} --outFile ${DIR}plotter_dilepton.root  --json ${JSONFILEDILEPTON} --cutflow all_initNorm --forceMerge --no2D --noPowers --plotExt .pdf --generatePseudoData --onlyStartWith optim_systs ${ONLYDILEPTON}
runFixedPlotter --iEcm 13 --iLumi ${LUMI} --inDir ${INDIR} --outDir ${DIR} --outFile ${DIR}plotter_dilepton.root  --json ${JSONFILEDILEPTON} --cutflow all_initNorm --useMerged  --no2D --noPowers --plotExt .png --generatePseudoData --onlyStartWith optim_systs ${ONLYDILEPTON}
runFixedPlotter --iEcm 13 --iLumi ${LUMI} --inDir ${INDIR} --outDir ${DIR} --outFile ${DIR}plotter_dilepton.root  --json ${JSONFILEDILEPTON} --cutflow all_initNorm --useMerged  --no2D --noPowers --plotExt .C   --generatePseudoData --onlyStartWith optim_systs ${ONLYDILEPTON}

# Leptontau
runFixedPlotter --iEcm 13 --iLumi ${LUMI} --inDir ${INDIR} --outDir ${DIR} --outFile ${DIR}plotter_ltau.root  --json ${JSONFILELEPTAU} --cutflow all_initNorm --forceMerge --no2D --noPowers --plotExt .pdf --generatePseudoData --onlyStartWith optim_systs ${ONLYLEPTAU}
runFixedPlotter --iEcm 13 --iLumi ${LUMI} --inDir ${INDIR} --outDir ${DIR} --outFile ${DIR}plotter_ltau.root  --json ${JSONFILELEPTAU} --cutflow all_initNorm --useMerged  --no2D --noPowers --plotExt .png --generatePseudoData --onlyStartWith optim_systs ${ONLYLEPTAU}
runFixedPlotter --iEcm 13 --iLumi ${LUMI} --inDir ${INDIR} --outDir ${DIR} --outFile ${DIR}plotter_ltau.root  --json ${JSONFILELEPTAU} --cutflow all_initNorm --useMerged  --no2D --noPowers --plotExt .C   --generatePseudoData --onlyStartWith optim_systs ${ONLYLEPTAU}

# Lessen the burden on the web browser

#for CHAN in emu ee mumu singlemu singlee
#do
#    mkdir ${DIR}temp${CHAN}/
#    mv ${DIR}${CHAN}* ${DIR}temp${CHAN}/
#    mv ${DIR}temp${CHAN}/ ${DIR}${CHAN}/
#    cp ~/www/HIG-13-026/index.php ${DIR}${CHAN}/
#done