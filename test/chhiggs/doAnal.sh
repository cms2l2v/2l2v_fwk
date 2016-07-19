#!/bin/bash

#JSONFILE=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/xsec_samples.json
JSONFILE=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/ttbar_samples.json

OUTDIR=$CMSSW_BASE/src/UserCode/llvv_fwk/test/chhiggs/results_ttbar/

QUEUE="1nh"
#QUEUE="criminal"
#QUEUE="8nm"
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
    elif [ "${2}" = "chhiggs" ]; then
        JSONFILE=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/chhiggs_samples.json
    elif [ "${2}" = "ntuplizer" ]; then
        JSONFILE=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/xsec_samples.json
    fi

    RESUBMIT=""
    if [ "${3}" = "resubmit" ]; then
        RESUBMIT=" -F "
    fi
    
    runAnalysisOverSamples.py -e runChHiggsAnalysis -j ${JSONFILE} -o ${OUTDIR} -d  /dummy/ -c $CMSSW_BASE/src/UserCode/llvv_fwk/test/runAnalysis_cfg.py.templ -p "@useMVA=False @saveSummaryTree=False @runSystematics=False @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0" ${RESUBMIT} -f 8 --NFile 2 -s ${QUEUE}
    
    
elif [ "${1}" = "lumi" ]; then 
    rm myjson.json
    # Normtag from: /afs/cern.ch/user/l/lumipro/public/normtag_file/moriond16_normtag.json
    NORMTAG="/afs/cern.ch/user/l/lumipro/public/normtag_file/moriond16_normtag.json"
    #cat ${OUTDIR}/*SingleMuon*.json > myjson.json
    cat ${OUTDIR}/Data13TeV_SingleMuon2015D_*.json > myjson.json
    STARTINGJSON="data/chhiggs/Cert_246908-258750_13TeV_PromptReco_Collisions15_25ns_JSON.txt"
    #STARTINGJSON="data/json/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt"
    STARTINGJSON="data/json/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt"
    sed -i -e "s#}{#, #g" myjson.json; 
    sed -i -e "s#, ,#, #g" myjson.json;
    echo "myjson.json has been recreated and the additional \"}{\" have been fixed."
    echo "Now running brilcalc according to the luminosity group recommendation:"
    echo "brilcalc lumi --normtag /afs/cern.ch/user/c/cmsbril/public/normtag_json/OfflineNormtagV1.json -i myjson.json"
    export PATH=$HOME/.local/bin:/opt/brilconda/bin:$PATH    
    brilcalc lumi --normtag ${NORMTAG} -i myjson.json
    echo "To be compared with the output of the full json:"
    echo "brilcalc lumi --normtag ${STARTINGJSON}"
    #brilcalc lumi -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-255031_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt -n 0.962
    brilcalc lumi --normtag ${NORMTAG} -i ${STARTINGJSON}
elif [ "${1}" = "pileup" ]; then
    echo "Computing pileup:"
    echo "pileupCalc.py -i data/json/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 69000 --maxPileupBin 50 --numPileupBins 50 MyDataPileupHistogram.root"
    pileupCalc.py -i data/json/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 69000 --maxPileupBin 50 --numPileupBins 50 MyDataPileupHistogram.root

elif [ "${1}" = "plot" ]; then 

    BASEDIR=~/www/13TeV_xsec_plots/
    mkdir -p ${BASEDIR}
    cp ~/www/HIG-13-026/index.php ${BASEDIR}

    JSONFILEDILEPTON=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/plot_xsec_dileptons.json
    JSONFILELEPTAU=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/plot_xsec_leptontau.json

    INDIR=${OUTDIR}
    ONLYDILEPTON=" --onlyStartWith emu_eventflow --onlyStartWith ee_eventflow --onlyStartWith mumu_eventflow --onlyStartWith emu_alteventflow --onlyStartWith ee_alteventflow --onlyStartWith mumu_alteventflow --onlyStartWith emu_step1leadpt --onlyStartWith emu_step6leadpt  --onlyStartWith emu_step6met --onlyStartWith emu_step1leadeta  --onlyStartWith emu_step3leadjetpt  --onlyStartWith emu_step3leadjeteta  --onlyStartWith emu_step3met  --onlyStartWith emu_step3nbjets --onlyStartWith ee_step1leadpt --onlyStartWith ee_step1leadeta  --onlyStartWith ee_step3leadjetpt  --onlyStartWith emu_step6leadpt  --onlyStartWith emu_step6met --onlyStartWith ee_step3leadjeteta  --onlyStartWith ee_step3met  --onlyStartWith ee_step3nbjets --onlyStartWith mumu_step1leadpt  --onlyStartWith emu_step6leadpt  --onlyStartWith emu_step6met --onlyStartWith mumu_step1leadeta  --onlyStartWith mumu_step3leadjetpt  --onlyStartWith mumu_step3leadjeteta  --onlyStartWith mumu_step3met  --onlyStartWith mumu_step3nbjets "
    
    ONLYLEPTAU=" --onlyStartWith singlemu_alteventflowslep --onlyStartWith singlee_alteventflowslep --onlyStartWith singlemu_altstep5met --onlyStartWith singlemu_altstep6met --onlyStartWith singlemu_altstep7met --onlyStartWith singlemu_altstep5nbtags --onlyStartWith singlemu_altstep6nbtags --onlyStartWith singlemu_altstep7nbtags --onlyStartWith singlemu_step3nbtags --onlyStartWith singlee_step3nbtags --onlyStartWith singlemu_step6nbtags --onlyStartWith singlee_step6nbtags --onlyStartWith singlemu_step3njets --onlyStartWith singlee_step3njets  --onlyStartWith singlemu_step6met --onlyStartWith singlee_step6met  --onlyStartWith singlemu_step6tauleadpt  --onlyStartWith singlemu_step6tauleadeta --onlyStartWith singlee_step6tauleadpt --onlyStartWith singlee_step6tauleadeta --onlyStartWith singlemu_final --onlyStartWith singlee_final --onlyStartWith singlee_eventflowslep --onlyStartWith singlemu_eventflowslep  --onlyStartWith singlemu_step1leadpt  --onlyStartWith singlemu_step1leadeta --onlyStartWith singlemu_step2leadjetpt --onlyStartWith singlemu_step2leadjeteta --onlyStartWith singlemu_step3met --onlyStartWith singlemu_step3nbjets --onlyStartWith singlee_step1leadpt  --onlyStartWith singlee_step1leadeta --onlyStartWith singlee_step2leadjetpt --onlyStartWith singlee_step2leadjeteta --onlyStartWith singlee_step3met --onlyStartWith singlee_step3nbjets --onlyStartWith singlemu_step1nvtx --onlyStartWith singlemu_step1nvtxraw --onlyStartWith singlee_step1nvtx --onlyStartWith singlee_step1nvtxraw  --onlyStartWith singlemu_step2leadpt  --onlyStartWith singlemu_step2leadeta   --onlyStartWith singlemu_step3leadpt  --onlyStartWith singlemu_step3leadeta  --onlyStartWith singlemu_step4leadpt  --onlyStartWith singlemu_step4leadeta  --onlyStartWith singlemu_step5leadpt  --onlyStartWith singlemu_step5leadeta  --onlyStartWith singlemu_step6leadpt  --onlyStartWith singlemu_step6leadeta  --onlyStartWith singlee_step2leadpt  --onlyStartWith singlee_step2leadeta   --onlyStartWith singlee_step3leadpt  --onlyStartWith singlee_step3leadeta  --onlyStartWith singlee_step4leadpt  --onlyStartWith singlee_step4leadeta  --onlyStartWith singlee_step5leadpt  --onlyStartWith singlee_step5leadeta  --onlyStartWith singlee_step6leadpt  --onlyStartWith singlee_step6leadeta --onlyStartWith singlemu_step4met --onlyStartWith singlemu_step5met --onlyStartWith singlee_step4met --onlyStartWith singlee_step5met "

    ONLYSYSLIST=" --onlyStartWith optim_systs "
    if [ "${2}" = "all" ]; then 
        ONLYDILEPTON=" --onlyStartWith ee --onlyStartWith mumu --onlyStartWith emu "
        ONLYLEPTAU=" --onlyStartWith singlemu --onlyStartWith singlee "
        JSONFILEDILEPTON=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/plot_chhiggs_dileptons.json
        JSONFILELEPTAU=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/plot_chhiggs_leptontau.json

    elif [ "${2}" = "evtflows" ]; then
        ONLYDILEPTON=" --onlyStartWith ee_xseceventflowdilep --onlyStartWith ee_chhiggseventflowdilep --onlyStartWith emu_xseceventflowdilep --onlyStartWith emu_chhiggseventflowdilep --onlyStartWith mumu_xseceventflowdilep --onlyStartWith mumu_chhiggseventflowdilep  "
        ONLYLEPTAU=" --onlyStartWith singlee_xseceventflowslep --onlyStartWidh singlee_chhiggseventflowslep --onlyStartWith singlemu_xseceventflowslep --onlyStartWith singlemu_chhiggseventflowslep "
    elif [ "${2}" = "chhiggs" ]; then
        ONLYDILEPTON=" --onlyStartWith ee_xseceventflowdilep --onlyStartWith ee_chhiggseventflowdilep --onlyStartWith emu_xseceventflowdilep --onlyStartWith emu_chhiggseventflowdilep --onlyStartWith mumu_xseceventflowdilep --onlyStartWith mumu_chhiggseventflowdilep  "
        ONLYLEPTAU=" --onlyStartWith singlee_xseceventflowslep --onlyStartWidh singlee_chhiggseventflowslep --onlyStartWith singlemu_xseceventflowslep --onlyStartWith singlemu_chhiggseventflowslep "
        
        JSONFILEDILEPTON=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/chhiggs_samples.json
        JSONFILELEPTAU=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/chhiggs_samples.json

        JSONFILEDILEPTON=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/plot_chhiggs_dileptons.json
        JSONFILELEPTAU=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/plot_chhiggs_leptontau.json

    fi
    

#    ONLYLEPTAU=" --onlyStartWith singlee_eventflowslep --onlyStartWith singlemu_eventflowslep"
    PSEUDODATA=" --generatePseudoData "
    PSEUDODATA=" "
    
#for LUMI in 500 1000 5000 10000
    for LUMI in 1280
    do
        DIR="${BASEDIR}${LUMI}/"
        #if [ "${2}" = "chhiggs" ]; then
        #    DIR="${DIR}/chhiggs/"
        #    LUMI=45000
        #fi
        #LUMI=${3}
        mkdir -p ${DIR}
        cp ~/www/HIG-13-026/index.php ${DIR}

        echo ${LUMI}
	echo ${INDIR}
	echo ${DIR}
	echo ${JSONFILEDILEPTON}
	echo ${JSONFILELEPTAU}
	echo ${PSEUDODATA}
	echo ${ONLYSYSLIST}
	echo ${ONLYDILEPTON}
	echo ${ONLYLEPTAU}
        # Dilepton
        runFixedPlotter --iEcm 13 --iLumi ${LUMI} --inDir ${INDIR} --outDir ${DIR} --outFile ${DIR}plotter_dilepton.root  --json ${JSONFILEDILEPTON} --cutflow all_initNorm --no2D --noPowers --plotExt .pdf --plotExt .png --plotExt .C ${PSEUDODATA} --onlyStartWith optim_systs ${ONLYDILEPTON}
        
        # Leptontau
        runFixedPlotter --iEcm 13 --iLumi ${LUMI} --inDir ${INDIR} --outDir ${DIR} --outFile ${DIR}plotter_ltau.root  --json ${JSONFILELEPTAU} --cutflow all_initNorm --no2D --noPowers --plotExt .pdf --plotExt .png --plotExt .C ${PSEUDODATA} ${ONLYSYSLIST} ${ONLYLEPTAU}
        
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