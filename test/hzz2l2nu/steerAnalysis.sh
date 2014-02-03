#!/bin/bash

#
# full analysis
#

step=$1
optim_step=$2
outdir="/afs/cern.ch/user/p/psilva/work/hzz_5311/"
indir="/store/cmst3/user/psilva/5311_ntuples"
cfg="$CMSSW_BASE/src/UserCode/llvv_fwk/test/runAnalysis_cfg.py.templ"

#prepare output directories
mkdir -p ${outdir}/ll/plots
mkdir -p ${outdir}/g/raw/plots
mkdir -p ${outdir}/g/qt/plots

lumi=19736

if [ "$step" == "1" ]; then
    echo "Submitting first pass (warning no systs yet)"
#    runLocalAnalysisOverSamples.py -e runHZZ2l2nuAnalysis -j data/htozz_samples.json        -d ${indir} -o ${outdir}/ll    -c ${cfg} -p "@runSystematics=True @useMVA=False @weightsFile='${CMSSW_BASE}/src/UserCode/llvv_fwk/data/weights/'" -s 8nh
#    runLocalAnalysisOverSamples.py -e runHZZ2l2nuAnalysis -j data/bulkg_samples.json        -d ${indir} -o ${outdir}/ll    -c ${cfg} -p "@runSystematics=True @useMVA=False @weightsFile='${CMSSW_BASE}/src/UserCode/llvv_fwk/data/weights/'" -s 8nh -t Bulk
#    runLocalAnalysisOverSamples.py -e runHZZ2l2nuAnalysis -j data/htozz_photon_samples.json -d ${indir} -o ${outdir}/g/raw -c ${cfg} -p "@runSystematics=False @useMVA=True"                                                                   -s 1nd
    runLocalAnalysisOverSamples.py -e runHZZ2l2nuAnalysis -j data/htozz_int_samples.json    -d ${indir} -o ${outdir}/ll    -c ${cfg} -p "@runSystematics=True @useMVA=False @weightsFile='${CMSSW_BASE}/src/UserCode/llvv_fwk/data/weights/'" -s 8nh
fi

if [ "$step" == "2" ]; then
    echo "Computing weights"
    runPlotter --iLumi ${lumi} --inDir ${outdir}/ll/    --json data/htozz_samples.json        --outFile ${outdir}/plotter.root       --forceMerged --outDir ${outdir}/ll/plots
#    runPlotter --iLumi ${lumi} --inDir ${outdir}/ll/    --json data/bulkg_samples.json        --outFile ${outdir}/plotter_bulkg.root --useMerged   --outDir ${outdir}/ll/plots
#    runPlotter --iLumi ${lumi} --inDir ${outdir}/g/raw/ --json data/htozz_photon_samples.json --outFile ${outdir}/plotter_g_raw.root --forceMerged --outDir ${outdir}/g/raw/plots
    #root -b -q "${CMSSW_BASE}/src/UserCode/llvv_fwk/test/ewkvp2j/FitQtSpectrum.C+(\"${outdir}/plotter.root\",\"${outdir}/plotter_g_raw.root\",ALL,HZZ)" 
    #mv gammawgts.root ${outdir}/gamma_weights.root
#    runPlotter --iLumi ${lumi} --inDir ${outdir}/ll/ --json data/htozz_int_samples.json --outFile ${outdir}/plotter_int.root --outDir ${outdir}/ll/plots --noPlot
fi

if [ "$step" == "3" ]; then
    echo "Submitting second pass, weighting in q_T"
    runLocalAnalysisOverSamples.py -e runHZZ2l2nuAnalysis -j data/htozz_photon_samples.json -d ${indir} -o ${outdir}/g/qt -c ${cfg} -p "@runSystematics=False @useMVA=True @weightsFile='${outdir}/gamma_weights.root'" -s 8nh
fi

if [ "$step" == "4" ]; then
    echo "Running final plotter"
    runPlotter --iLumi ${lumi} --inDir ${outdir}/g/qt/ --json data/htozz_photon_samples.json --outFile ${outdir}/plotter_g_qt.root --forceMerged --outDir ${outdir}/g/qt/plots
fi

if [ "$step" == "5" ]; then
    echo "Running DY substitution"
    root -b -q "${CMSSW_BASE}/src/UserCode/llvv_fwk/test/ewkvp2j/runDYprediction.C+(\"${outdir}/plotter.root\",\"${outdir}/plotter_g_qt.root\",\"\",\"\",SEARCH)"
fi

if [ "$step" == "6" ]; then
    echo "FIXME"
    echo "Running optimization"
    #python test/ewkvp2j/optimize_VBFZ.py -p ${optim_step} -b True -s Fisher_shapes -o Fisher_analysis  
fi



