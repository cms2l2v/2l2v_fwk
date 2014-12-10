#!/bin/bash

#rm ~/www/TStauStau/MCvsData/2012{B,C,D,BCD}/{exclusive,inclusive}{/,/DY/}*.{png,pdf,root,C,tex}
OUTDIR=./OUT
INDIR=/lustre/ncg.ingrid.pt/cmslocal/cbeiraod/Results
PLOTDIR=./Plots

if [[ -d $PLOTDIR/ ]]; then
  rm -Rf $PLOTDIR/
fi
mkdir $PLOTDIR/

runPlotterFWLite --noPowers --iEcm 8 --iLumi 19672 --inDir $INDIR/ --outDir $OUTDIR/ --outFile $OUTDIR/plotter.root --json tstaustau_separate_wjets.json --plotExt .png --plotExt .root
cp $OUTDIR/Q* $PLOTDIR/
rm -Rf $OUTDIR


makeDistributions --json tstaustau_samples_full.json --inDir $INDIR/ --outDir $OUTDIR/ --variables variables.json --signalSelection "(stauMass-neutralinoMass)==100"
cp $OUTDIR/cosPhi.png                      $PLOTDIR/DeltaM100_cosPhi.png
cp $OUTDIR/Q80.png                         $PLOTDIR/DeltaM100_Q80.png
cp $OUTDIR/Q80VsCosPhi.png                 $PLOTDIR/DeltaM100_Q80PlusCosPhi.png
cp $OUTDIR/TStauStau_120_20_cosPhi_Q80.png $PLOTDIR/DeltaM100_cosPhiVsQ80.png
rm -Rf $OUTDIR


makeDistributions --json tstaustau_samples_full.json --inDir $INDIR/ --outDir $OUTDIR/ --variables variables.json --signalSelection "(stauMass-neutralinoMass)==20"
cp $OUTDIR/cosPhi.png                      $PLOTDIR/DeltaM20_cosPhi.png
cp $OUTDIR/Q80.png                         $PLOTDIR/DeltaM20_Q80.png
cp $OUTDIR/Q80VsCosPhi.png                 $PLOTDIR/DeltaM20_Q80PlusCosPhi.png
cp $OUTDIR/TStauStau_120_20_cosPhi_Q80.png $PLOTDIR/DeltaM20_cosPhiVsQ80.png
rm -Rf $OUTDIR


makeDistributions --json tstaustau_samples_full.json --inDir $INDIR/ --outDir $OUTDIR/ --variables variables.json --signalSelection "(stauMass-neutralinoMass)==50"
cp $OUTDIR/cosPhi.png                      $PLOTDIR/DeltaM50_cosPhi.png
cp $OUTDIR/Q80.png                         $PLOTDIR/DeltaM50_Q80.png
cp $OUTDIR/Q80VsCosPhi.png                 $PLOTDIR/DeltaM50_Q80PlusCosPhi.png
cp $OUTDIR/TStauStau_120_20_cosPhi_Q80.png $PLOTDIR/DeltaM50_cosPhiVsQ80.png
rm -Rf $OUTDIR


makeDistributions --json tstaustau_samples_full.json --inDir $INDIR/ --outDir $OUTDIR/ --variables variables.json --signalSelection "(stauMass-neutralinoMass)==200"
cp $OUTDIR/cosPhi.png                      $PLOTDIR/DeltaM200_cosPhi.png
cp $OUTDIR/Q80.png                         $PLOTDIR/DeltaM200_Q80.png
cp $OUTDIR/Q80VsCosPhi.png                 $PLOTDIR/DeltaM200_Q80PlusCosPhi.png
cp $OUTDIR/TStauStau_120_20_cosPhi_Q80.png $PLOTDIR/DeltaM200_cosPhiVsQ80.png
rm -Rf $OUTDIR


makeDistributions --json tstaustau_samples_full.json --inDir $INDIR/ --outDir $OUTDIR/ --variables variables.json --signalSelection "(stauMass==50) &&(neutralinoMass==20)"
cp $OUTDIR/cosPhi.png                      $PLOTDIR/50_20_cosPhi.png
cp $OUTDIR/Q80.png                         $PLOTDIR/50_20_Q80.png
cp $OUTDIR/Q80VsCosPhi.png                 $PLOTDIR/50_20_Q80PlusCosPhi.png
cp $OUTDIR/TStauStau_120_20_cosPhi_Q80.png $PLOTDIR/50_20_cosPhiVsQ80.png
rm -Rf $OUTDIR


makeDistributions --json tstaustau_samples_full.json --inDir $INDIR/ --outDir $OUTDIR/ --variables variables.json --signalSelection "(stauMass==120) &&(neutralinoMass==20)"
cp $OUTDIR/cosPhi.png                      $PLOTDIR/120_20_cosPhi.png
cp $OUTDIR/Q80.png                         $PLOTDIR/120_20_Q80.png
cp $OUTDIR/Q80VsCosPhi.png                 $PLOTDIR/120_20_Q80PlusCosPhi.png
cp $OUTDIR/TStauStau_120_20_cosPhi_Q80.png $PLOTDIR/120_20_cosPhiVsQ80.png
rm -Rf $OUTDIR


makeDistributions --json tstaustau_samples_full.json --inDir $INDIR/ --outDir $OUTDIR/ --variables variables.json --signalSelection "(stauMass==120) &&(neutralinoMass==90)"
cp $OUTDIR/cosPhi.png                      $PLOTDIR/120_90_cosPhi.png
cp $OUTDIR/Q80.png                         $PLOTDIR/120_90_Q80.png
cp $OUTDIR/Q80VsCosPhi.png                 $PLOTDIR/120_90_Q80PlusCosPhi.png
cp $OUTDIR/TStauStau_120_20_cosPhi_Q80.png $PLOTDIR/120_90_cosPhiVsQ80.png
rm -Rf $OUTDIR


makeDistributions --json tstaustau_samples_full.json --inDir $INDIR/ --outDir $OUTDIR/ --variables variables.json --signalSelection "(stauMass==240) &&(neutralinoMass==20)"
cp $OUTDIR/cosPhi.png                      $PLOTDIR/240_20_cosPhi.png
cp $OUTDIR/Q80.png                         $PLOTDIR/240_20_Q80.png
cp $OUTDIR/Q80VsCosPhi.png                 $PLOTDIR/240_20_Q80PlusCosPhi.png
cp $OUTDIR/TStauStau_120_20_cosPhi_Q80.png $PLOTDIR/240_20_cosPhiVsQ80.png
rm -Rf $OUTDIR


makeDistributions --json tstaustau_samples_full.json --inDir $INDIR/ --outDir $OUTDIR/ --variables variables.json --signalSelection "(stauMass==240) &&(neutralinoMass==90)"
cp $OUTDIR/cosPhi.png                      $PLOTDIR/240_90_cosPhi.png
cp $OUTDIR/Q80.png                         $PLOTDIR/240_90_Q80.png
cp $OUTDIR/Q80VsCosPhi.png                 $PLOTDIR/240_90_Q80PlusCosPhi.png
cp $OUTDIR/TStauStau_120_20_cosPhi_Q80.png $PLOTDIR/240_90_cosPhiVsQ80.png
rm -Rf $OUTDIR


makeDistributions --json tstaustau_samples_full.json --inDir $INDIR/ --outDir $OUTDIR/ --variables variables.json --signalSelection "(stauMass==240) &&(neutralinoMass==210)"
cp $OUTDIR/cosPhi.png                      $PLOTDIR/240_210_cosPhi.png
cp $OUTDIR/Q80.png                         $PLOTDIR/240_210_Q80.png
cp $OUTDIR/Q80VsCosPhi.png                 $PLOTDIR/240_210_Q80PlusCosPhi.png
cp $OUTDIR/TStauStau_120_20_cosPhi_Q80.png $PLOTDIR/240_210_cosPhiVsQ80.png
rm -Rf $OUTDIR
