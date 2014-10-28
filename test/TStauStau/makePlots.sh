#!/bin/bash

#rm ~/www/TStauStau/MCvsData/2012{B,C,D,BCD}/{exclusive,inclusive}{/,/DY/}*.{png,pdf,root,C,tex}

## 2012ABCD
runPlotterFWLite --noPowers --iEcm 8 --iLumi 19672 --inDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/Results/ --outDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/Plots/ --outFile /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/Plots/plotter.root --json tstaustau_samples_full.json --plotExt .png --plotExt .root

## 2012A
runPlotterFWLite --noPowers --iEcm 8 --iLumi 876.225 --inDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/Results/ --outDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/Plots_2012A/ --outFile /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/Plots_2012A/plotter.root --json tstaustau_samples_2012A.json --plotExt .png --plotExt .root

## 2012B
runPlotterFWLite --noPowers --iEcm 8 --iLumi 4412 --inDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/Results/ --outDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/Plots_2012B/ --outFile /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/Plots_2012B/plotter.root --json tstaustau_samples_2012B.json --plotExt .png --plotExt .root

## 2012C
runPlotterFWLite --noPowers --iEcm 8 --iLumi 7055 --inDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/Results/ --outDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/Plots_2012C/ --outFile /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/Plots_2012C/plotter.root --json tstaustau_samples_2012C.json --plotExt .png --plotExt .root

## 2012D
runPlotterFWLite --noPowers --iEcm 8 --iLumi 7329 --inDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/Results/ --outDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/Plots_2012D/ --outFile /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/Plots_2012D/plotter.root --json tstaustau_samples_2012D.json --plotExt .png --plotExt .root
