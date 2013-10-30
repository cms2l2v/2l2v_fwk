#!/bin/bash

INPUTDIR="~/work/ewkzp2j_5311/ll/"

echo "Running minimal training"
root -l -b -q tmvaClassifier.C++'("Fisher,BDTD,MLP","${INPUTDIR}",false)'
mv TMVA.root TMVA_minimal.root

echo "Running full training"
root -l -b -q tmvaClassifier.C++'("Fisher,BDTD,MLP","${INPUTDIR}",true)'
mv TMVA.root TMVA_full.root

