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
    #    BASEDIR=/afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_e/
    BASEDIR=/afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_g/
    if [ "${2}" = "anal_sus" ]; then
	runLocalAnalysisOverSamples.py -e runChHiggsAnalysis -j data/chhiggs/ch-higgs_samples.json -d /afs/cern.ch/work/v/vischia/private/store/5311_ntuples/ -o ${BASEDIR} -c test/runAnalysis_cfg.py.templ -p "@runSystematics=True @saveSummaryTree=False" -s 8nh
    elif [ "${2}" = "anal_sm" ]; then
       runLocalAnalysisOverSamples.py -e runChHiggsAnalysis -j data/top_samples.json      -d /store/cmst3/user/psilva/5311_ntuples/             -o ${BASEDIR}     -c test/runAnalysis_cfg.py.templ -p "@runSystematics=True @saveSummaryTree=False" -s 8nh

    elif [ "${2}" = "plots" ]; then
	# Plots
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_evtflow_pdf.root             --showUnc --plotExt .pdf --noPowers --onlyStartWith emu_evtflow                                 &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_geq2btagsmet_pdf.root        --showUnc --plotExt .pdf --noPowers --onlyStartWith emu_geq2btagsmet		     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_geq2btagsnbjets_pdf.root     --showUnc --plotExt .pdf --noPowers --onlyStartWith emu_geq2btagsnbjets		     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_geq2btagsptlep_pdf.root      --showUnc --plotExt .pdf --noPowers --onlyStartWith emu_geq2btagsptlep		     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_geq2btagssumpt_pdf.root      --showUnc --plotExt .pdf --noPowers --onlyStartWith emu_geq2btagssumpt		     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_met_pdf.root                 --showUnc --plotExt .pdf --noPowers --onlyStartWith emu_met					     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_mll_pdf.root                 --showUnc --plotExt .pdf --noPowers --onlyStartWith emu_mll					     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_mtsum_pdf.root               --showUnc --plotExt .pdf --noPowers --onlyStartWith emu_mtsum				     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_nbjets_pdf.root              --showUnc --plotExt .pdf --noPowers --onlyStartWith emu_nbjets				     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_njets_pdf.root               --showUnc --plotExt .pdf --noPowers --onlyStartWith emu_njets				     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_nvertices_pdf.root           --showUnc --plotExt .pdf --noPowers --onlyStartWith emu_nvertices			     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_nverticesUnweighted_pdf.root --showUnc --plotExt .pdf --noPowers --onlyStartWith emu_nverticesUnweighted	     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_pte_pdf.root                 --showUnc --plotExt .pdf --noPowers --onlyStartWith emu_pte					     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_ptjet1eta_pdf.root           --showUnc --plotExt .pdf --noPowers --onlyStartWith emu_ptjet1eta			     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_ptjet1pt_pdf.root            --showUnc --plotExt .pdf --noPowers --onlyStartWith emu_ptjet1pt				     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_ptjet2eta_pdf.root           --showUnc --plotExt .pdf --noPowers --onlyStartWith emu_ptjet2eta			     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_ptjet2pt_pdf.root            --showUnc --plotExt .pdf --noPowers --onlyStartWith emu_ptjet2pt				     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_ptmin_pdf.root               --showUnc --plotExt .pdf --noPowers --onlyStartWith emu_ptmin				     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_ptmu_pdf.root                --showUnc --plotExt .pdf --noPowers --onlyStartWith emu_ptmu					     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_sumpt_pdf.root               --showUnc --plotExt .pdf --noPowers --onlyStartWith emu_sumpt                                     &

	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_evtflow_png.root             --showUnc --plotExt .png --noPowers --onlyStartWith emu_evtflow                                 &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_geq2btagsmet_png.root        --showUnc --plotExt .png --noPowers --onlyStartWith emu_geq2btagsmet		     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_geq2btagsnbjets_png.root     --showUnc --plotExt .png --noPowers --onlyStartWith emu_geq2btagsnbjets		     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_geq2btagsptlep_png.root      --showUnc --plotExt .png --noPowers --onlyStartWith emu_geq2btagsptlep		     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_geq2btagssumpt_png.root      --showUnc --plotExt .png --noPowers --onlyStartWith emu_geq2btagssumpt		     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_met_png.root                 --showUnc --plotExt .png --noPowers --onlyStartWith emu_met					     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_mll_png.root                 --showUnc --plotExt .png --noPowers --onlyStartWith emu_mll					     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_mtsum_png.root               --showUnc --plotExt .png --noPowers --onlyStartWith emu_mtsum				     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_nbjets_png.root              --showUnc --plotExt .png --noPowers --onlyStartWith emu_nbjets				     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_njets_png.root               --showUnc --plotExt .png --noPowers --onlyStartWith emu_njets				     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_nvertices_png.root           --showUnc --plotExt .png --noPowers --onlyStartWith emu_nvertices			     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_nverticesUnweighted_png.root --showUnc --plotExt .png --noPowers --onlyStartWith emu_nverticesUnweighted	     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_pte_png.root                 --showUnc --plotExt .png --noPowers --onlyStartWith emu_pte					     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_ptjet1eta_png.root           --showUnc --plotExt .png --noPowers --onlyStartWith emu_ptjet1eta			     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_ptjet1pt_png.root            --showUnc --plotExt .png --noPowers --onlyStartWith emu_ptjet1pt				     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_ptjet2eta_png.root           --showUnc --plotExt .png --noPowers --onlyStartWith emu_ptjet2eta			     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_ptjet2pt_png.root            --showUnc --plotExt .png --noPowers --onlyStartWith emu_ptjet2pt				     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_ptmin_png.root               --showUnc --plotExt .png --noPowers --onlyStartWith emu_ptmin				     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_ptmu_png.root                --showUnc --plotExt .png --noPowers --onlyStartWith emu_ptmu					     &
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/plot-ch-higgs_samples.json --outFile ${BASEDIR}plotter-forPlotting_sumpt_png.root               --showUnc --plotExt .png --noPowers --onlyStartWith emu_sumpt                                     &


#	mkdir ${BASEDIR}plots/tabs
#	mv ${BASEDIR}plotse* ${BASEDIR}plots/tabs/
    elif [ "${2}" = "tables" ]; then
	# Tables
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/chhiggs/all-samples.json --outFile ${BASEDIR}plotter-forTables.root --showUnc --noPlots --noPowers --onlyStartWith emu_evtflow
    elif [ "${2}" = "datacards" ]; then
	# Datacards
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}datacards/ --json data/chhiggs/all-samples.json --outFile ${BASEDIR}plotter-forSystTableInPAS.root --showUnc --noPlots --noPowers --onlyStartWith emu_evtflow
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}datacards/ --json data/chhiggs/all-samples.json --outFile ${BASEDIR}plotter-forSystTable.root --showUnc --noPlots --noPowers --onlyStartWith emu_finalevtflow2btags
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/chhiggs/all-samples_higgs1pb.json --outFile ${BASEDIR}plotter-all-samplesForDatacards_finalevtflow_norm.root --noPlot --noPowers  --onlyStartWith emu_finalevtflow2btags
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}datacards/ --json data/chhiggs/all-samples.json --outFile ${BASEDIR}plotter-forSystTableInPAS_optim.root --showUnc --noPlots --noPowers --onlyStartWith all_optim
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}datacards/ --json data/chhiggs/all-samples.json --outFile ${BASEDIR}plotter-forSystTable_optim.root --showUnc --noPlots --noPowers --onlyStartWith all_optim
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/chhiggs/all-samples_higgs1pb.json --outFile ${BASEDIR}plotter-all-samplesForDatacards_finalevtflow_norm_optim.root --noPlot --noPowers  --onlyStartWith all_optim


hadd -f ${BASEDIR}plotter-forSystTableInPAS_def.root			        ${BASEDIR}plotter-forSystTableInPAS.root			 	${BASEDIR}plotter-forSystTableInPAS_optim.root
hadd -f ${BASEDIR}plotter-forSystTable_def.root                              ${BASEDIR}plotter-forSystTable.root                             	${BASEDIR}plotter-forSystTable_optim.root                             
hadd -f ${BASEDIR}plotter-all-samplesForDatacards_finalevtflow_norm_def.root ${BASEDIR}plotter-all-samplesForDatacards_finalevtflow_norm.root	${BASEDIR}plotter-all-samplesForDatacards_finalevtflow_norm_optim.root

	sh test/chhiggs/doAllDatacardsByDecay.sh
    fi
fi
