#!/bin/bash

#INPUTDIR="~/work/ewkzp2j_5311/ll/"
INPUTDIR="/afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_5315_mva/"

echo "Running minimal training"
#root -b -q "${CMSSW_BASE}/src/UserCode/llvv_fwk/test/chhiggs/tmvaClassifier.C+(\"Fisher,BDTD,MLP\",\"${INPUTDIR}\",true,false)"
root -b -q "${CMSSW_BASE}/src/UserCode/llvv_fwk/test/chhiggs/tmvaClassifier.C+(\"BDTD\",\"${INPUTDIR}\",true,false)"
mv TMVA.root TMVA_minimal.root

# Does not matter
#echo "Running full training but no q/g discriminator"
#root -b -q "${CMSSW_BASE}/src/UserCode/llvv_fwk/test/chhiggs/tmvaClassifier.C+(\"Fisher,BDTD,MLP\",\"${INPUTDIR}\",false,false)"
#mv TMVA.root TMVA_noqg_full.root

echo "Running full training"
#root -b -q "${CMSSW_BASE}/src/UserCode/llvv_fwk/test/chhiggs/tmvaClassifier.C+(\"Fisher,BDTD,MLP\",\"${INPUTDIR}\",false,true)"
root -b -q "${CMSSW_BASE}/src/UserCode/llvv_fwk/test/chhiggs/tmvaClassifier.C+(\"BDTD\",\"${INPUTDIR}\",false,true)"
mv TMVA.root TMVA_full.root

