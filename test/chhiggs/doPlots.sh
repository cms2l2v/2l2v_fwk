#!/bin/bash


if [ "${1}" = "past" ]; then
    runPlotter --iLumi 19702 --inDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_b/ --outDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_d/plots --json data/plot-ch-higgs_samples.json --outFile /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_d/plotter-with-ch-higgs_new.root --showUnc --plotExt .pdf --noPowers --onlyStartWith emu
    runPlotter --iLumi 19702 --inDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_b/ --outDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_d/plots --json data/plot-ch-higgs_samples.json --outFile /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_d/plotter-with-ch-higgs_new.root --showUnc --plotExt .png --noPowers --onlyStartWith emu
	#
    #
    #Tables only
    #-----------
    #runPlotter --iLumi 19702 --inDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311/ --outDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311/plots --json data/ch-higgs_samples.json --outFile /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311/plotter-all-ch-higgs.root --noPlot --noPowers 
    #runPlotter --iLumi 19702 --inDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311/ --outDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311/plots --json data/top_samples.json --outFile /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311/plotter-sm-only.root --noPlot --noPowers 
    #
    #
    #


elif [ "${1}" = "current" ]; then
# Fixed run 

    if [ "${2}" = "plots" ]; then
	# Plots
	runPlotter --iLumi 19702 --inDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_e/ --outDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_e/plots --json data/plot-ch-higgs_samples.json --outFile /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_e/plotter-forPlotting.root --showUnc --plotExt .pdf --noPowers --onlyStartWith emu
	runPlotter --iLumi 19702 --inDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_e/ --outDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_e/plots --json data/plot-ch-higgs_samples.json --outFile /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_e/plotter-forPlotting.root --showUnc --plotExt .png --noPowers --onlyStartWith emu
    elif [ "${2}" = "tables" ]; then
	# Tables
	runPlotter --iLumi 19702 --inDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_e/ --outDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_e/plots --json data/chhiggs/all-samples.json --outFile /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_e/plotter-forTables.root --showUnc --noPlots --noPowers --onlyStartWith emu
    elif [ "${2}" = "datacards" ]; then
	# Datacards
    runPlotter --iLumi 19702 --inDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_e/ --outDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_e/plots --json data/chhiggs/all-samples_higgs1pb.json --outFile /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_e/plotter-forDatacards.root --showUnc --noPlots --noPowers --onlyStartWith emu

    fi
fi
