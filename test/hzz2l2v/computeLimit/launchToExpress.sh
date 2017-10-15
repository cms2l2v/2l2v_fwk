#Should be launched with qsub -j oe -q express -F "0800 100.0" launchToExpress.sh
source $VO_CMS_SW_DIR/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc491
export BUILD_ARCH=slc6_amd64_gcc491
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
export XRD_NETWORKSTACK=IPv4
cd /storage_mnt/storage/user/delannoy/HZZ2l2v/2_October_LimitsImprovement/CMSSW_7_4_7/src/UserCode/llvv_fwk/test/hzz2l2v/computeLimit
eval `scramv1 runtime -sh`
sh produceLimitsAndPlotsFromDatacards.sh $1 $2
