#!/bin/bash

#
# full analysis
#

step=$1
optim_step=$2
outdir="/afs/cern.ch/user/p/psilva/work/ewkzp2j_5311"
indir="/store/cmst3/user/psilva/5311_ntuples"
cfg="$CMSSW_BASE/src/UserCode/llvv_fwk/test/runAnalysis_cfg.py.templ"


#prepare output directories
mkdir -p ${outdir}/ll/
mkdir -p ${outdir}/g/data/raw_loose
mkdir -p ${outdir}/g/data/raw_tight
mkdir -p ${outdir}/g/data/qt_loose
mkdir -p ${outdir}/g/data/qt_tight
mkdir -p ${outdir}/g/data/qt_pure

if [ "$step" == "1" ]; then
    echo "Submitting first pass"
    runLocalAnalysisOverSamples.py -e runVBFZAnalysis -j data/vbfz_samples.json        -d ${indir} -o ${outdir}/ll/              -c ${cfg} -p "@runSystematics=False @useMVA=True"  -s 1nd
    runLocalAnalysisOverSamples.py -e runVBFZAnalysis -j data/vbfz_photon_samples.json -d ${indir} -o ${outdir}/g/data/raw_tight -c ${cfg} -p "@runSystematics=False @useMVA=True" -s 1nd
    runLocalAnalysisOverSamples.py -e runVBFZAnalysis -j data/vbfz_photon_samples.json -d ${indir} -o ${outdir}/g/data/raw_loose -c ${cfg} -p "@runSystematics=True @useMVA=True"  -s 1nd
fi

if [ "$step" == "2" ]; then
    echo "Computing weights"
    runPlotter --iLumi 19736 --inDir ${outdir}/ll/               --json data/vbfz_samples.json        --outFile ${outdir}/plotter.root             --forceMerged --noPlot
    runPlotter --iLumi 19736 --inDir ${outdir}/g/data/raw_tight/ --json data/vbfz_photon_samples.json --outFile ${outdir}/plotter_g_tight_raw.root --forceMerged --noPlot
    runPlotter --iLumi 19736 --inDir ${outdir}/g/data/raw_loose/ --json data/vbfz_photon_samples.json --outFile ${outdir}/plotter_g_loose_raw.root --forceMerged --noPlot

    root -b -q "${CMSSW_BASE}/src/UserCode/llvv_fwk/test/ewkvp2j/FitQtSpectrum.C+(\"${outdir}/plotter.root\",\"${outdir}/plotter_g_loose_raw.root\",ALL)" 
    mv gammawgts.root ${outdir}/loose_gamma_weights.root

    root -b -q "${CMSSW_BASE}/src/UserCode/llvv_fwk/test/ewkvp2j/FitQtSpectrum.C+(\"${outdir}/plotter.root\",\"${outdir}/plotter_g_tight_raw.root\",ALL)" 
    mv gammawgts.root ${outdir}/tight_gamma_weights.root

fi

if [ "$step" == "3" ]; then
    echo "Submitting second pass, weighting in q_T"
    runLocalAnalysisOverSamples.py -e runVBFZAnalysis -j data/vbfz_photon_samples.json -d ${indir} -o ${outdir}/g/data/qt_tight -c ${cfg} -p "@runSystematics=False @useMVA=True @weightsFile='${outdir}/tight_gamma_weights.root'" -s 8nh
    runLocalAnalysisOverSamples.py -e runVBFZAnalysis -j data/vbfz_photon_samples.json -d ${indir} -o ${outdir}/g/data/qt_loose -c ${cfg} -p "@runSystematics=True @useMVA=True @weightsFile='${outdir}/loose_gamma_weights.root'"  -s 8nh
fi

if [ "$step" == "4" ]; then
    echo "Running final plotter"
    runPlotter --iLumi 19736 --inDir ${outdir}/g/data/qt_tight/ --json data/vbfz_photon_samples.json --outFile ${outdir}/plotter_g_tight_qt.root --forceMerged --noPlot
    runPlotter --iLumi 19736 --inDir ${outdir}/g/data/qt_loose/ --json data/vbfz_photon_samples.json --outFile ${outdir}/plotter_g_loose_qt.root --forceMerged --noPlot
fi

if [ "$step" == "5" ]; then
    echo "Running optimization"
    python test/ewkvp2j/optimize_VBFZ.py -p ${optim_step} -b True -s MLP_shapes -o MLP_analysis
fi

if [ "$step" == "6" ]; then
    echo "Running final plotter"
fi
