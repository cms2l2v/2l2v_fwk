#!/bin/bash

#rm ~/www/TStauStau/MCvsData/2012{B,C,D,BCD}/{exclusive,inclusive}{/,/DY/}*.{png,pdf,root,C,tex}
OUTDIR=/lustre/ncg.ingrid.pt/cmslocal/cbeiraod/Plots_noSF
INDIR=/lustre/ncg.ingrid.pt/cmslocal/cbeiraod/Results_noSF

## 2012ABCD
runPlotterFWLite --noPowers --iEcm 8 --iLumi 19672 --inDir $INDIR/ --outDir $OUTDIR/ --outFile $OUTDIR/plotter.root --json tstaustau_samples_full.json --plotExt .png --plotExt .root

## 2012A
runPlotterFWLite --noPowers --iEcm 8 --iLumi 876.225 --inDir $INDIR/ --outDir "$OUTDIR""_2012A/" --outFile "$OUTDIR""_2012A/plotter.root" --json tstaustau_samples_2012A.json --plotExt .png --plotExt .root

## 2012B
runPlotterFWLite --noPowers --iEcm 8 --iLumi 4412 --inDir $INDIR/ --outDir "$OUTDIR""_2012B/" --outFile "$OUTDIR""_2012B/plotter.root" --json tstaustau_samples_2012B.json --plotExt .png --plotExt .root

## 2012C
runPlotterFWLite --noPowers --iEcm 8 --iLumi 7055 --inDir $INDIR/ --outDir "$OUTDIR""_2012C/" --outFile "$OUTDIR""_2012C/plotter.root" --json tstaustau_samples_2012C.json --plotExt .png --plotExt .root

## 2012D
runPlotterFWLite --noPowers --iEcm 8 --iLumi 7329 --inDir $INDIR/ --outDir "$OUTDIR""_2012D/" --outFile "$OUTDIR""_2012D/plotter.root" --json tstaustau_samples_2012D.json --plotExt .png --plotExt .root
