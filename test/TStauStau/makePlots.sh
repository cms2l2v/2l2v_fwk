#!/bin/bash

#rm ~/www/TStauStau/MCvsData/2012{B,C,D,BCD}/{exclusive,inclusive}{/,/DY/}*.{png,pdf,root,C,tex}

# Only Run2012B
## LIP
runPlotterFWLite --noPowers --iEcm 8 --iLumi 4429 --inDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/LIP_selection/ --outDir ./LIP_Selection_2012B/ --outFile plotter_2012B_Exclusive_LIP.root --json tstaustau_samples_exclusive_2012B.json --plotExt .png --plotExt .pdf --plotExt .C --plotExt .root
## IPM
runPlotterFWLite --noPowers --iEcm 8 --iLumi 4429 --inDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/IPM_selection/ --outDir ./IPM_Selection_2012B/ --outFile plotter_2012B_Exclusive_IPM.root --json tstaustau_samples_exclusive_2012B.json --plotExt .png --plotExt .pdf --plotExt .C --plotExt .root

# Only Run2012C
## LIP
runPlotterFWLite --noPowers --iEcm 8 --iLumi 7147 --inDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/LIP_selection/ --outDir ./LIP_Selection_2012C/ --outFile plotter_2012C_Exclusive_LIP.root --json tstaustau_samples_exclusive_2012C.json --plotExt .png --plotExt .pdf --plotExt .C --plotExt .root
## IPM
runPlotterFWLite --noPowers --iEcm 8 --iLumi 7147 --inDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/IPM_selection/ --outDir ./IPM_Selection_2012C/ --outFile plotter_2012C_Exclusive_IPM.root --json tstaustau_samples_exclusive_2012C.json --plotExt .png --plotExt .pdf --plotExt .C --plotExt .root

# Only Run2012BC
## LIP
runPlotterFWLite --noPowers --iEcm 8 --iLumi 11576 --inDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/LIP_selection/ --outDir ./LIP_Selection_2012BC/ --outFile plotter_2012BC_Exclusive_LIP.root --json tstaustau_samples_exclusive_2012BC.json --plotExt .png --plotExt .pdf --plotExt .C --plotExt .root
## IPM
runPlotterFWLite --noPowers --iEcm 8 --iLumi 11576 --inDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/IPM_selection/ --outDir ./IPM_Selection_2012BC/ --outFile plotter_2012BC_Exclusive_IPM.root --json tstaustau_samples_exclusive_2012BC.json --plotExt .png --plotExt .pdf --plotExt .C --plotExt .root

# Only Run2012D
## LIP
runPlotterFWLite --noPowers --iEcm 8 --iLumi 7317 --inDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/LIP_selection/ --outDir ./LIP_Selection_2012D/ --outFile plotter_2012D_Exclusive_LIP.root --json tstaustau_samples_exclusive_2012D.json --plotExt .png --plotExt .pdf --plotExt .C --plotExt .root
## IPM
runPlotterFWLite --noPowers --iEcm 8 --iLumi 7317 --inDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/IPM_selection/ --outDir ./IPM_Selection_2012D/ --outFile plotter_2012D_Exclusive_IPM.root --json tstaustau_samples_exclusive_2012D.json --plotExt .png --plotExt .pdf --plotExt .C --plotExt .root

# Runs 2012{B,C,D}
## LIP
runPlotterFWLite --noPowers --iEcm 8 --iLumi 18893 --inDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/LIP_selection/ --outDir ./LIP_Selection_2012BCD/ --outFile plotter_2012BCD_Exclusive_LIP.root --json tstaustau_samples_exclusive_2012BCD.json --plotExt .png --plotExt .pdf --plotExt .C --plotExt .root
## IPM
runPlotterFWLite --noPowers --iEcm 8 --iLumi 18893 --inDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/IPM_selection/ --outDir ./IPM_Selection_2012BCD/ --outFile plotter_2012BCD_Exclusive_IPM.root --json tstaustau_samples_exclusive_2012BCD.json --plotExt .png --plotExt .pdf --plotExt .C --plotExt .root
