#!/bin/bash

OUTRATIO=~/www/13TeV_xsec/ratio/
DEMU="${OUTRATIO}eventflow/DataCard_emu.dat"
DMUTAU="${OUTRATIO}eventflowslep/DataCard_singlemu.dat"
DRATIO="${OUTRATIO}DataCard_ratio.dat"
WRATIO="${OUTRATIO}DataCard_ratio.root"
RATIOAREA="${OUTRATIO}fit"

mkdir -p ${RATIOAREA}

combineCards.py Name1=${DEMU} Name2=${DMUTAU} > ${DRATIO}

text2workspace.py -P UserCode.llvv_fwk.XsecRatioModel:xsecRatioModel ${DRATIO} -o ${WRATIO}

# runPLRanalysis --in ${WRATIO} --lumi 1000 --ecm 13 --cl 0.68 
# mv PLR* ${OUTRATIO}
# mv textable.tex ${OUTRATIO}


exit 0
