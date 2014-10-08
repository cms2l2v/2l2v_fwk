#!/bin/bash

cd /exper-sw/cmst3/cmssw/users/cbeiraod/Stau_Production/CMSSW_5_3_15/src
#cmsenv
eval `scramv1 runtime -sh`
cd -

addSVfit --inFile {inFile} --outFile {outFile}
