#!/bin/bash

#
# closure test
#

step=$1
outdir="/afs/cern.ch/user/p/psilva/work/hzz_5311/"
indir="/store/cmst3/user/psilva/5311_ntuples"
cfg="$CMSSW_BASE/src/UserCode/llvv_fwk/test/runAnalysis_cfg.py.templ"


#prepare output directories
mkdir -p ${outdir}/dy/
mkdir -p ${outdir}/g/raw
mkdir -p ${outdir}/g/qt_pure

if [ "$step" == "1" ]; then
    echo "Submitting first pass"
    runLocalAnalysisOverSamples.py -e runHZZ2l2nuAnalysis -j data/vbfz_syst_samples.json   -d ${indir} -o ${outdir}/dy          -c ${cfg} -p "@runSystematics=False @weightsFile='${CMSSW_BASE}/src/UserCode/llvv_fwk/data/weights/'" -s 8nh -t DY;
    runLocalAnalysisOverSamples.py -e runHZZ2l2nuAnalysis -j data/htozz_photon_samples.json -d ${indir} -o ${outdir}/g/raw -c ${cfg} -p "@runSystematics=False @useMVA=True" -t MC -s 8nh
fi

if [ "$step" == "2" ]; then
    echo "Computing weights"
    runPlotter --json data/vbfz_syst_samples.json   --inDir ${outdir}/dy/     --iLumi 19800 --iEcm 8 --outFile ${outdir}/plotter_dy_closure.root             --noPlot --forceMerged;
    runPlotter --json data/htozz_photon_samples.json --inDir ${outdir}/g/raw/ --iLumi 19800 --iEcm 8 --outFile ${outdir}/plotter_dy_closure_g_raw.root --noPlot --forceMerged;
    root -b -q "${CMSSW_BASE}/src/UserCode/llvv_fwk/test/ewkvp2j/FitQtSpectrum.C+(\"${outdir}/plotter_dy_closure.root\",\"${outdir}/plotter_dy_closure_g_raw.root\",PUREG,HZZ)";
    mv gammawgts.root ${outdir}/pure_gamma_mcweights.root;
fi

if [ "$step" == "3" ]; then
    echo "Submitting second pass, weighting in q_T"
    runLocalAnalysisOverSamples.py -e runHZZ2l2nuAnalysis -j data/htozz_photon_samples.json -d ${indir} -o ${outdir}/g/qt_pure  -c ${cfg} -p "@runSystematics=False @useMVA=True @weightsFile='${outdir}/pure_gamma_mcweights.root'"                       -s 8nh -t MC;
fi

if [ "$step" == "4" ]; then
    echo "Running final plotter"
    runPlotter --json data/htozz_photon_samples.json --inDir ${outdir}/g/qt_pure/  --iLumi 19800 --iEcm 8 --outFile ${outdir}/plotter_dy_closure_g_qt_pure.root  --noPlot --forceMerged;

fi

if [ "$step" == "5" ]; then
    echo "Producing the closure tests report"
fi
