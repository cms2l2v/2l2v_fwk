#! /bin/sh

pwd
export SCRAM_ARCH=slc6_amd64_gcc472
export BUILD_ARCH=slc5_amd64_gcc462
export VO_CMS_SW_DIR=/nfs/soft/cms
cd /afs/cern.ch/work/v/vischia/private/code/CMSSW_5_3_15/src/UserCode/llvv_fwk
eval `scramv1 runtime -sh`
cd test/chhiggs/
sh tmvaClassifier.sh
