#!/bin/bash

RTOGEN=${1}
STEP=${2}
OUTPUT=${3}
INPUT=/afs/cern.ch/user/p/psilva/work/top_5311/plotter_nom.root
SYSTINPUT=/afs/cern.ch/user/p/psilva/work/top_5311/plotter_syst.root
MYCMSSW=/afs/cern.ch/user/p/psilva/work/CMSSW_5_3_11
#PARFILE=${MYCMSSW}/src/UserCode/llvv_fwk/test/top/hfcParams_2012_mc_cfg.json
PARFILE=${MYCMSSW}/src/UserCode/llvv_fwk/test/top/hfcParams_2012_pythiamc_cfg.json
#BTAGFILE=${MYCMSSW}/src/UserCode/llvv_fwk/test/top/csvL_2012_mc_cfg.json
BTAGFILE=${MYCMSSW}/src/UserCode/llvv_fwk/test/top/csvL_2012_pythiamc_cfg.json

JOB=${RTOGEN}${STEP}

mkdir -p ${OUTPUT}
cd ${MYCMSSW}/src
export SCRAM_ARCH=slc5_amd64_gcc462
eval `scram r -sh`
mkdir /tmp/hfc_${JOB}
cd /tmp/hfc_${JOB}
fitHeavyFlavorContent --par ${PARFILE}  --btag ${BTAGFILE} --in ${INPUT} --fit 0 --npe 5 --seed ${STEP} --syst ${SYSTINPUT} --rtogen ${RTOGEN} --fast
mv HFCClosureTest.root ${OUTPUT}/HFCClosureTest_${JOB}.root
mv *.pdf ${OUTPUT}/
cd -
rm -rf /tmp/hfc_${JOB}
