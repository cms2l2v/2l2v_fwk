#!/bin/bash

#500 1000 5000 10000
for LUMI in 16.1
do
    OUTRATIO=~/www/13TeV_xsec/${LUMI}/ratio/
    DEMU="${OUTRATIO}eventflow/DataCard_emu.dat"
    DMUTAU="${OUTRATIO}eventflowslep/DataCard_singlemu.dat"
    DRATIO="${OUTRATIO}DataCard_ratio.dat"
    WRATIO="${OUTRATIO}DataCard_ratio.root"
    RATIOAREA="${OUTRATIO}fit"
    
    mkdir -p ${RATIOAREA}
    cp ~/www/HIG-13-026/index.php ${RATIOAREA}

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
done

exit 0
