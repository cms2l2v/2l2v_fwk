#!/bin/bash

#DIR=~/www/13TeV_topsel_fix_ala9/
#DIR=~/www/13TeV_topsel_fix_talk/
BASEDIR=~/www/13TeV_xsec_plots/
mkdir -p ${BASEDIR}
rm -rf ${BASEDIR}*
cp ~/www/HIG-13-026/index.php ${BASEDIR}

JSONFILEDILEPTON=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/test/phys14_plot_dileptons.json
JSONFILELEPTAU=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/test/phys14_plot_leptontau.json

JSONFILEDILEPTON=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_plot_dileptons.json
JSONFILELEPTAU=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_plot_leptontau.json

INDIR=$CMSSW_BASE/src/UserCode/llvv_fwk/test/chhiggs/results_ttbar/
ONLYDILEPTON=" --onlyStartWith emu_eventflow --onlyStartWith ee_eventflow --onlyStartWith mumu_eventflow --onlyStartWith emu_step1leadpt --onlyStartWith emu_step6leadpt  --onlyStartWith emu_step6met --onlyStartWith emu_step1leadeta  --onlyStartWith emu_step3leadjetpt  --onlyStartWith emu_step3leadjeteta  --onlyStartWith emu_step3met  --onlyStartWith emu_step3nbjets --onlyStartWith ee_step1leadpt --onlyStartWith ee_step1leadeta  --onlyStartWith ee_step3leadjetpt  --onlyStartWith emu_step6leadpt  --onlyStartWith emu_step6met --onlyStartWith ee_step3leadjeteta  --onlyStartWith ee_step3met  --onlyStartWith ee_step3nbjets --onlyStartWith mumu_step1leadpt  --onlyStartWith emu_step6leadpt  --onlyStartWith emu_step6met --onlyStartWith mumu_step1leadeta  --onlyStartWith mumu_step3leadjetpt  --onlyStartWith mumu_step3leadjeteta  --onlyStartWith mumu_step3met  --onlyStartWith mumu_step3nbjets "

ONLYLEPTAU=" --onlyStartWith singlemu_step6met --onlyStartWith singlee_step6met  --onlyStartWith singlemu_step6tauleadpt --onlyStartWith singlee_step6tauleadpt --onlyStartWith singlemu_final --onlyStartWith singlee_final --onlyStartWith singlee_eventflowslep --onlyStartWith singlemu_eventflowslep  --onlyStartWith singlemu_step1leadpt  --onlyStartWith singlemu_step1leadeta --onlyStartWith singlemu_step2leadjetpt --onlyStartWith singlemu_step2leadjeteta --onlyStartWith singlemu_step3met --onlyStartWith singlemu_step3nbjets --onlyStartWith singlee_step1leadpt  --onlyStartWith singlee_step1leadeta --onlyStartWith singlee_step2leadjetpt --onlyStartWith singlee_step2leadjeteta --onlyStartWith singlee_step3met --onlyStartWith singlee_step3nbjets"

PSEUDODATA=" --generatePseudoData "
PSEUDODATA=" "

#for LUMI in 500 1000 5000 10000
for LUMI in 5.7
do
    DIR="${BASEDIR}${LUMI}/"
    mkdir -p ${DIR}
    cp ~/www/HIG-13-026/index.php ${DIR}

    # Dilepton
    runFixedPlotter --iEcm 13 --iLumi ${LUMI} --inDir ${INDIR} --outDir ${DIR} --outFile ${DIR}plotter_dilepton.root  --json ${JSONFILEDILEPTON} --cutflow all_initNorm --no2D --noPowers --plotExt .pdf ${PSEUDODATA} --onlyStartWith optim_systs ${ONLYDILEPTON}
    runFixedPlotter --iEcm 13 --iLumi ${LUMI} --inDir ${INDIR} --outDir ${DIR} --outFile ${DIR}plotter_dilepton.root  --json ${JSONFILEDILEPTON} --cutflow all_initNorm --no2D --noPowers --plotExt .png ${PSEUDODATA} --onlyStartWith optim_systs ${ONLYDILEPTON}
    runFixedPlotter --iEcm 13 --iLumi ${LUMI} --inDir ${INDIR} --outDir ${DIR} --outFile ${DIR}plotter_dilepton.root  --json ${JSONFILEDILEPTON} --cutflow all_initNorm --no2D --noPowers --plotExt .C   ${PSEUDODATA} --onlyStartWith optim_systs ${ONLYDILEPTON}
    
    # Leptontau
    runFixedPlotter --iEcm 13 --iLumi ${LUMI} --inDir ${INDIR} --outDir ${DIR} --outFile ${DIR}plotter_ltau.root  --json ${JSONFILELEPTAU} --cutflow all_initNorm --no2D --noPowers --plotExt .pdf ${PSEUDODATA} --onlyStartWith optim_systs ${ONLYLEPTAU}
    runFixedPlotter --iEcm 13 --iLumi ${LUMI} --inDir ${INDIR} --outDir ${DIR} --outFile ${DIR}plotter_ltau.root  --json ${JSONFILELEPTAU} --cutflow all_initNorm --no2D --noPowers --plotExt .png ${PSEUDODATA} --onlyStartWith optim_systs ${ONLYLEPTAU}
    runFixedPlotter --iEcm 13 --iLumi ${LUMI} --inDir ${INDIR} --outDir ${DIR} --outFile ${DIR}plotter_ltau.root  --json ${JSONFILELEPTAU} --cutflow all_initNorm --no2D --noPowers --plotExt .C   ${PSEUDODATA} --onlyStartWith optim_systs ${ONLYLEPTAU}
    
done

exit 0
# Lessen the burden on the web browser

#for CHAN in emu ee mumu singlemu singlee
#do
#    mkdir ${DIR}temp${CHAN}/
#    mv ${DIR}${CHAN}* ${DIR}temp${CHAN}/
#    mv ${DIR}temp${CHAN}/ ${DIR}${CHAN}/
#    cp ~/www/HIG-13-026/index.php ${DIR}${CHAN}/
#done