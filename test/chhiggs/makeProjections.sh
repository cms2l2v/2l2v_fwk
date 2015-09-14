#!/bin/bash


for LUMI in 500 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6000 6500 7000 7500 8000 8500 9000 9500 10000 10500 11000 11500 12000 12500 13000 13500 14000 14500 15000 15500 16000 16500 17000 17500 18000 18500 19000 20000 20500 21000 21500 22000 22500 23000 23500 24000 24500 25000 25500 26000 26500 27000 27500 28000 28500 29000 30000
do
    echo "LUMINOSITY ${LUMI}" >> makeEfficiencies.C
    sh test/chhiggs/doAnal.sh plot evtflows ${LUMI}
    sh test/chhiggs/doAnal.sh plot chhiggs ${LUMI}
    for FINALSTATE in emu ee mumu singlemu singlee
    do
        SUFFIX=""
        TOGREP="dileptons"
        if [ "${FINALSTATE}" = "singlemu" ] || [ "${FINALSTATE}" = "singlee" ]; then
            SUFFIX="slep"
            TOGREP="ltau"
        fi
        rm ${FINALSTATE}_eventflow${SUFFIX}.tex
        echo "${FINALSTATE}" >> makeEfficiencies.C
        wget http://vischia.web.cern.ch/vischia/13TeV_xsec_plots/16.1/chhiggs/${FINALSTATE}_eventflow${SUFFIX}.tex 
        cat ${FINALSTATE}_eventflow${SUFFIX}.tex | grep M-150 >> makeEfficiencies.C
        cat ${FINALSTATE}_eventflow${SUFFIX}.tex | grep M-500 >> makeEfficiencies.C
        rm ${FINALSTATE}_eventflow${SUFFIX}.tex
        wget http://vischia.web.cern.ch/vischia/13TeV_xsec_plots/16.1/${FINALSTATE}_eventflow${SUFFIX}.tex 
        cat ${FINALSTATE}_eventflow${SUFFIX}.tex | grep ${TOGREP} >> makeEfficiencies.C
        cat ${FINALSTATE}_eventflow${SUFFIX}.tex | grep Total >> makeEfficiencies.C
        rm ${FINALSTATE}_eventflow${SUFFIX}.tex
    done
done


exit 0