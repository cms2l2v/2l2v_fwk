#!/bin/bash

#export SCRAM_ARCH=slc6_amd64_gcc472
#export BUILD_ARCH=slc5_amd64_gcc462
#export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
#cd /exper-sw/cmst3/cmssw/users/cbeiraod/SLC6/CMSSW_5_3_15/src/
#eval `scramv1 runtime -sh`

#cd /exper-sw/cmst3/cmssw/users/cbeiraod/SLC6/CMSSW_5_3_15/src/UserCode/llvv_fwk/test/TStauStau/CutOptimization/
#runCutOptimizer --json cutOptim.json --outDir ./Results/ --plotExt .png --plotExt .root

qsub runRound0.sh
qsub runRound1.sh
qsub runRound2.sh
qsub runRound3.sh
qsub runRound4.sh
qsub runRound5.sh
qsub runRound6.sh
qsub runRound7.sh
qsub runRound8.sh
qsub runRound9.sh
qsub runRound10.sh
qsub runRound11.sh
qsub runRound12.sh
qsub runRound13.sh
qsub runRound14.sh
qsub runRound15.sh
qsub runRound16.sh
qsub runRound17.sh
qsub runRound18.sh
qsub runRound19.sh
qsub runRound20.sh
qsub runRound21.sh
qsub runRound22.sh
qsub runRound23.sh
qsub runRound24.sh
qsub runRound25.sh
qsub runRound26.sh
qsub runRound27.sh
qsub runRound28.sh
qsub runRound29.sh
qsub runRound30.sh
qsub runRound31.sh
qsub runRound32.sh
qsub runRound33.sh
qsub runRound34.sh
qsub runRound35.sh
qsub runRound36.sh
qsub runRound37.sh
qsub runRound38.sh
qsub runRound39.sh
qsub runRound40.sh
qsub runRound41.sh
