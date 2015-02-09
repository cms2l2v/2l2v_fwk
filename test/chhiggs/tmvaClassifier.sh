#!/bin/bash

#INPUTDIR="~/work/ewkzp2j_5311/ll/"
INPUTDIR="/afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_5315_mva/"

for method in Fisher MLP BDTD
do
    echo "Running minimal training"
    root -b -q "${CMSSW_BASE}/src/UserCode/llvv_fwk/test/chhiggs/tmvaClassifier.C+(\"${method}\",\"${INPUTDIR}\",true,false)"
    mv TMVA.root TMVA_${method}_minimal.root

    # Does not matter
    #echo "Running full training but no q/g discriminator"
    #root -b -q "${CMSSW_BASE}/src/UserCode/llvv_fwk/test/chhiggs/tmvaClassifier.C+(\"Fisher,BDTD,MLP\",\"${INPUTDIR}\",false,false)"
    #mv TMVA.root TMVA_noqg_full.root
    
    echo "Running full training"
    root -b -q "${CMSSW_BASE}/src/UserCode/llvv_fwk/test/chhiggs/tmvaClassifier.C+(\"${method}\",\"${INPUTDIR}\",false,true)"
    mv TMVA.root TMVA_${method}_full.root
done

exit 0
