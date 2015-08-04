#!/bin/bash

#JSONFILE=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/xsec_samples.json
JSONFILE=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/ttbar_samples.json

#OUTDIR=$CMSSW_BASE/src/UserCode/llvv_fwk/test/chhiggs/results_ttbar/
OUTDIR=$CMSSW_BASE/src/UserCode/llvv_fwk/test/chhiggs/results_ttbar_2/

QUEUE="8nh"
#QUEUE="crab"

mkdir -p ${OUTDIR}

if [ "${1}" = "submit" ]; then 

    if [ "${2}" = "all" ]; then 
        JSONFILE=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/xsec_samples.json
    elif [ "${2}" = "data" ]; then 
        JSONFILE=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/data_samples.json
    elif [ "${2}" = "mc" ]; then 
        JSONFILE=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/mc_samples.json
    elif [ "${2}" = "ttbar" ]; then 
        JSONFILE=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/ttbar_samples.json
    fi
    
    runAnalysisOverSamples.py -e runChHiggsAnalysis -j ${JSONFILE} -o ${OUTDIR} -d  /dummy/ -c $CMSSW_BASE/src/UserCode/llvv_fwk/test/runAnalysis_cfg.py.templ -p "@useMVA=False @saveSummaryTree=False @runSystematics=False @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0" -s ${QUEUE}

elif [ "${1}" = "lumi" ]; then 
    rm myjson.json
    cat test/chhiggs/results_ttbar/*SingleMuon*.json > myjson.json
    echo "myjson.json has been recreated."
    echo "Now please edit myjson.json to remove the additional \"}{\" characters"
    echo "After that, please run:"
    echo "python scripts/lumi/lcr2/lcr2.py -i myjson.json"

elif [ "${1}" = "plot" ]; then 

    BASEDIR=~/www/13TeV_xsec_plots/
    mkdir -p ${BASEDIR}
    cp ~/www/HIG-13-026/index.php ${BASEDIR}

    JSONFILEDILEPTON=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/test/phys14_plot_dileptons.json
    JSONFILELEPTAU=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/test/phys14_plot_leptontau.json
    
    JSONFILEDILEPTON=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_plot_dileptons.json
    JSONFILELEPTAU=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_plot_leptontau.json

    #JSONFILEDILEPTON=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_plot_dileptons_inclusivetop.json
    #JSONFILELEPTAU=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_plot_leptontau_inclusivetop.json
    
    INDIR=${OUTDIR}
    ONLYDILEPTON=" --onlyStartWith emu_eventflow --onlyStartWith ee_eventflow --onlyStartWith mumu_eventflow --onlyStartWith emu_step1leadpt --onlyStartWith emu_step6leadpt  --onlyStartWith emu_step6met --onlyStartWith emu_step1leadeta  --onlyStartWith emu_step3leadjetpt  --onlyStartWith emu_step3leadjeteta  --onlyStartWith emu_step3met  --onlyStartWith emu_step3nbjets --onlyStartWith ee_step1leadpt --onlyStartWith ee_step1leadeta  --onlyStartWith ee_step3leadjetpt  --onlyStartWith emu_step6leadpt  --onlyStartWith emu_step6met --onlyStartWith ee_step3leadjeteta  --onlyStartWith ee_step3met  --onlyStartWith ee_step3nbjets --onlyStartWith mumu_step1leadpt  --onlyStartWith emu_step6leadpt  --onlyStartWith emu_step6met --onlyStartWith mumu_step1leadeta  --onlyStartWith mumu_step3leadjetpt  --onlyStartWith mumu_step3leadjeteta  --onlyStartWith mumu_step3met  --onlyStartWith mumu_step3nbjets "
    
#    ONLYLEPTAU=" --onlyStartWith singlemu_step6met --onlyStartWith singlee_step6met  --onlyStartWith singlemu_step6tauleadpt  --onlyStartWith singlemu_step6tauleadeta --onlyStartWith singlee_step6tauleadpt --onlyStartWith singlee_step6tauleadeta --onlyStartWith singlemu_final --onlyStartWith singlee_final --onlyStartWith singlee_eventflowslep --onlyStartWith singlemu_eventflowslep  --onlyStartWith singlemu_step1leadpt  --onlyStartWith singlemu_step1leadeta --onlyStartWith singlemu_step2leadjetpt --onlyStartWith singlemu_step2leadjeteta --onlyStartWith singlemu_step3met --onlyStartWith singlemu_step3nbjets --onlyStartWith singlee_step1leadpt  --onlyStartWith singlee_step1leadeta --onlyStartWith singlee_step2leadjetpt --onlyStartWith singlee_step2leadjeteta --onlyStartWith singlee_step3met --onlyStartWith singlee_step3nbjets"
    ONLYLEPTAU=" --onlyStartWith singlee_eventflowslep --onlyStartWith singlemu_eventflowslep --onlyStartWith singlemu_step1nvtx --onlyStartWith singlemu_step1nvtxraw --onlyStartWith singlee_step1nvtx --onlyStartWith singlee_step1nvtxraw "
    PSEUDODATA=" --generatePseudoData "
    PSEUDODATA=" "
    
#for LUMI in 500 1000 5000 10000
#   29.539 16.670 
    for LUMI in 40.240
    do
        DIR="${BASEDIR}${LUMI}_2/"
        mkdir -p ${DIR}
        cp ~/www/HIG-13-026/index.php ${DIR}

        # Dilepton
        #runFixedPlotter --iEcm 13 --iLumi ${LUMI} --inDir ${INDIR} --outDir ${DIR} --outFile ${DIR}plotter_dilepton.root  --json ${JSONFILEDILEPTON} --cutflow all_initNorm --no2D --noPowers --plotExt .pdf --plotExt .png --plotExt .C ${PSEUDODATA} --onlyStartWith optim_systs ${ONLYDILEPTON}
        
        # Leptontau
        runFixedPlotter --iEcm 13 --iLumi ${LUMI} --inDir ${INDIR} --outDir ${DIR} --outFile ${DIR}plotter_ltau.root  --json ${JSONFILELEPTAU} --cutflow all_initNorm --no2D --noPowers --plotExt .pdf --plotExt .png --plotExt .C ${PSEUDODATA} --onlyStartWith optim_systs ${ONLYLEPTAU}
        
    done
   
    # Lessen the burden on the web browser
    
    #for CHAN in emu ee mumu singlemu singlee
    #do
    #    mkdir ${DIR}temp${CHAN}/
    #    mv ${DIR}${CHAN}* ${DIR}temp${CHAN}/
    #    mv ${DIR}temp${CHAN}/ ${DIR}${CHAN}/
    #    cp ~/www/HIG-13-026/index.php ${DIR}${CHAN}/
    #done

fi

exit 0