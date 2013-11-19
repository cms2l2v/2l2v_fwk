#!/bin/bash

INPUTDIR="~/work/ewkzp2j_5311/ll/"

echo "Running minimal training"
root -b -q "${CMSSW_BASE}/src/UserCode/llvv_fwk/test/ewkvp2j/tmvaClassifier.C+(\"Fisher,BDTD,MLP\",\"${INPUTDIR}\",true,false)"
mv TMVA.root TMVA_minimal.root

echo "Running full training but no q/g discriminator"
root -b -q "${CMSSW_BASE}/src/UserCode/llvv_fwk/test/ewkvp2j/tmvaClassifier.C+(\"Fisher,BDTD,MLP\",\"${INPUTDIR}\",false,false)"
mv TMVA.root TMVA_noqg_full.root

echo "Running full training"
root -b -q "${CMSSW_BASE}/src/UserCode/llvv_fwk/test/ewkvp2j/tmvaClassifier.C+(\"Fisher,BDTD,MLP\",\"${INPUTDIR}\",false,true)"
mv TMVA.root TMVA_full.root

