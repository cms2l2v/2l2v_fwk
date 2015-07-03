#!/bin/bash

OUTRATIO=~/www/13TeV_xsec/ratio/
DEMU="${OUTRATIO}eventflow/DataCard_emu.dat"
DMUTAU="${OUTRATIO}eventflowslep/DataCard_singlemu.dat"
DRATIO="${OUTRATIO}DataCard_ratio.dat"
WRATIO="${OUTRATIO}DataCard_ratio.root"
RATIOAREA="${OUTRATIO}fit"

mkdir -p ${RATIOAREA}

echo "COMBINING CARDS"
combineCards.py Name1=${DEMU} Name2=${DMUTAU} > ${DRATIO}

echo "BUILDING WORKSPACE"
text2workspace.py -v 3 -P UserCode.llvv_fwk.XsecRatioModel:xsecRatioModel ${DRATIO} -o ${WRATIO}

dot -Tpng -o ${WRATIO}_directedgraph.png ${WRATIO}.dot
fdp -Tpng -o ${WRATIO}_springmodel.png ${WRATIO}.dot

convert ${WRATIO}_directedgraph.png ${WRATIO}_directedgraph.pdf
convert ${WRATIO}_springmodel.png   ${WRATIO}_springmodel.pdf
# runPLRanalysis --in ${WRATIO} --lumi 1000 --ecm 13 --cl 0.68 
# mv PLR* ${OUTRATIO}
# mv textable.tex ${OUTRATIO}


exit 0
