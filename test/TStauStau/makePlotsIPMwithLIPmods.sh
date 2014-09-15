#!/bin/bash

# Only Run2012B
## IPM with LIP Tau ID
runPlotterFWLite --noPowers --iEcm 8 --iLumi 4429 --inDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/IPM_selection_LIPTauID/ --outDir ./IPM_Selection_LIPTauID_2012B/ --outFile plotter_2012B_Exclusive_IPM_LIPTauID.root --json tstaustau_samples_exclusive_2012B.json --plotExt .png --plotExt .pdf --plotExt .C --plotExt .root
## IPM with LIP b-veto
runPlotterFWLite --noPowers --iEcm 8 --iLumi 4429 --inDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/IPM_selection_LIPBVeto/ --outDir ./IPM_Selection_LIPBVeto_2012B/ --outFile plotter_2012B_Exclusive_IPM_LIPBVeto.root --json tstaustau_samples_exclusive_2012B.json --plotExt .png --plotExt .pdf --plotExt .C --plotExt .root

# Only Run2012C
## IPM with LIP Tau ID
runPlotterFWLite --noPowers --iEcm 8 --iLumi 7147 --inDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/IPM_selection_LIPTauID/ --outDir ./IPM_Selection_LIPTauID_2012C/ --outFile plotter_2012C_Exclusive_IPM_LIPTauID.root --json tstaustau_samples_exclusive_2012C.json --plotExt .png --plotExt .pdf --plotExt .C --plotExt .root
## IPM with LIP b-veto
runPlotterFWLite --noPowers --iEcm 8 --iLumi 7147 --inDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/IPM_selection_LIPBVeto/ --outDir ./IPM_Selection_LIPBVeto_2012C/ --outFile plotter_2012C_Exclusive_IPM_LIPBVeto.root --json tstaustau_samples_exclusive_2012C.json --plotExt .png --plotExt .pdf --plotExt .C --plotExt .root

# Only Run2012BC
## IPM with LIP Tau ID
runPlotterFWLite --noPowers --iEcm 8 --iLumi 11576 --inDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/IPM_selection_LIPTauID/ --outDir ./IPM_Selection_LIPTauID_2012BC/ --outFile plotter_2012BC_Exclusive_IPM_LIPTauID.root --json tstaustau_samples_exclusive_2012BC.json --plotExt .png --plotExt .pdf --plotExt .C --plotExt .root
## IPM with LIP b-veto
runPlotterFWLite --noPowers --iEcm 8 --iLumi 11576 --inDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/IPM_selection_LIPBVeto/ --outDir ./IPM_Selection_LIPBVeto_2012BC/ --outFile plotter_2012BC_Exclusive_IPM_LIPBVeto.root --json tstaustau_samples_exclusive_2012BC.json --plotExt .png --plotExt .pdf --plotExt .C --plotExt .root

# Only Run2012D
## IPM with LIP Tau ID
runPlotterFWLite --noPowers --iEcm 8 --iLumi 7317 --inDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/IPM_selection_LIPTauID/ --outDir ./IPM_Selection_LIPTauID_2012D/ --outFile plotter_2012D_Exclusive_IPM_LIPTauID.root --json tstaustau_samples_exclusive_2012D.json --plotExt .png --plotExt .pdf --plotExt .C --plotExt .root
## IPM with LIP b-veto
runPlotterFWLite --noPowers --iEcm 8 --iLumi 7317 --inDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/IPM_selection_LIPBVeto/ --outDir ./IPM_Selection_LIPBVeto_2012D/ --outFile plotter_2012D_Exclusive_IPM_LIPBVeto.root --json tstaustau_samples_exclusive_2012D.json --plotExt .png --plotExt .pdf --plotExt .C --plotExt .root

# Only Run2012BCD
## IPM with LIP Tau ID
runPlotterFWLite --noPowers --iEcm 8 --iLumi 18893 --inDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/IPM_selection_LIPTauID/ --outDir ./IPM_Selection_LIPTauID_2012BCD/ --outFile plotter_2012BCD_Exclusive_IPM_LIPTauID.root --json tstaustau_samples_exclusive_2012BCD.json --plotExt .png --plotExt .pdf --plotExt .C --plotExt .root
## IPM with LIP b-veto
runPlotterFWLite --noPowers --iEcm 8 --iLumi 18893 --inDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/IPM_selection_LIPBVeto/ --outDir ./IPM_Selection_LIPBVeto_2012BCD/ --outFile plotter_2012BCD_Exclusive_IPM_LIPBVeto.root --json tstaustau_samples_exclusive_2012BCD.json --plotExt .png --plotExt .pdf --plotExt .C --plotExt .root
