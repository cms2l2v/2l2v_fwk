#!/bin/bash

#
# closure test
#

step=$1
outdir="/afs/cern.ch/user/p/psilva/work/ewkzp2j_5311"
indir="/store/cmst3/user/psilva/5311_ntuples"
cfg="$CMSSW_BASE/src/UserCode/llvv_fwk/test/runAnalysis_cfg.py.templ"


#prepare output directories
mkdir -p ${outdir}/dy/
mkdir -p ${outdir}/g/raw_loose
mkdir -p ${outdir}/g/raw_tight
mkdir -p ${outdir}/g/qt_loose
mkdir -p ${outdir}/g/qt_tight
mkdir -p ${outdir}/g/qt_pure

if [ "$step" == "0" ]; then
    echo "Submitting PDF variations"
    runLocalAnalysisOverSamples.py -e computePDFvariations -j data/vbfz_samples.json -o ${outdir}/dy -d ${indir} -c ${cfg} -s 1nd;
    runLocalAnalysisOverSamples.py -e computePDFvariations -j data/vbfz_syst_samples.json -o ${outdir}/dy -d ${indir} -c ${cfg} -s 1nd;
    runLocalAnalysisOverSamples.py -e computePDFvariations -j data/vbfz_photon_samples.json -d ${indir} -o ${outdir}/g/raw_loose  -c ${cfg} -s 1nd;
fi

if [ "$step" == "1" ]; then
    echo "Submitting first pass"
    runLocalAnalysisOverSamples.py -e runVBFZAnalysis -j data/vbfz_syst_samples.json   -d ${indir} -o ${outdir}/dy          -c ${cfg} -p "@runSystematics=False @useMVA=True" -s 8nh -t DY;
    runLocalAnalysisOverSamples.py -e runVBFZAnalysis -j data/vbfz_photon_samples.json -d ${indir} -o ${outdir}/g/raw_tight -c ${cfg} -p "@useMVA=True"                       -s 8nh -t MC;
    runLocalAnalysisOverSamples.py -e runVBFZAnalysis -j data/vbfz_photon_samples.json -d ${indir} -o ${outdir}/g/raw_loose -c ${cfg} -p "@runSystematics=True @useMVA=True"  -s 8nh -t MC;
fi

if [ "$step" == "2" ]; then
    echo "Computing weights"
    runPlotter --json data/vbfz_syst_samples.json   --inDir ${outdir}/dy/          --iLumi 19800 --iEcm 8 --outFile ${outdir}/plotter_dy_closure.root             --noPlot --forceMerged;
    runPlotter --json data/vbfz_photon_samples.json --inDir ${outdir}/g/raw_tight/ --iLumi 19800 --iEcm 8 --outFile ${outdir}/plotter_dy_closure_g_raw_tight.root --noPlot --forceMerged;
    runPlotter --json data/vbfz_photon_samples.json --inDir ${outdir}/g/raw_loose/ --iLumi 19800 --iEcm 8 --outFile ${outdir}/plotter_dy_closure_g_raw_loose.root --noPlot --forceMerged;

    root -b -q "${CMSSW_BASE}/src/UserCode/llvv_fwk/test/ewkvp2j/FitQtSpectrum.C+(\"${outdir}/plotter_dy_closure.root\",\"${outdir}/plotter_dy_closure_g_raw_loose.root\",ALL)";
    mv gammawgts.root ${outdir}/loose_gamma_mcweights.root;
   
    root -b -q "${CMSSW_BASE}/src/UserCode/llvv_fwk/test/ewkvp2j/FitQtSpectrum.C+(\"${outdir}/plotter_dy_closure.root\",\"${outdir}/plotter_dy_closure_g_raw_loose.root\",PUREG)";
    mv gammawgts.root ${outdir}/pure_gamma_mcweights.root;
   
    root -b -q "${CMSSW_BASE}/src/UserCode/llvv_fwk/test/ewkvp2j/FitQtSpectrum.C+(\"${outdir}/plotter_dy_closure.root\",\"${outdir}/plotter_dy_closure_g_raw_tight.root\",ALL)";
    mv gammawgts.root ${outdir}/tight_gamma_mcweights.root;
fi

if [ "$step" == "3" ]; then
    echo "Submitting second pass, weighting in q_T"
    runLocalAnalysisOverSamples.py -e runVBFZAnalysis -j data/vbfz_photon_samples.json -d ${indir} -o ${outdir}/g/qt_tight -c ${cfg} -p "@useMVA=True @weightsFile='${outdir}/tight_gamma_mcweights.root'"                      -s 8nh -t MC;
    runLocalAnalysisOverSamples.py -e runVBFZAnalysis -j data/vbfz_photon_samples.json -d ${indir} -o ${outdir}/g/qt_loose -c ${cfg} -p "@runSystematics=True @useMVA=True @weightsFile='${outdir}/loose_gamma_mcweights.root'" -s 8nh -t MC;
    runLocalAnalysisOverSamples.py -e runVBFZAnalysis -j data/vbfz_photon_samples.json -d ${indir} -o ${outdir}/g/qt_pure  -c ${cfg} -p "@useMVA=True @weightsFile='${outdir}/pure_gamma_mcweights.root'"                       -s 8nh -t MC;
fi

if [ "$step" == "4" ]; then
    echo "Running final plotter"
    runPlotter --json data/vbfz_photon_samples.json --inDir ${outdir}/g/qt_tight/ --iLumi 19800 --iEcm 8 --outFile ${outdir}/plotter_dy_closure_g_qt_tight.root --noPlot --forceMerged;
    runPlotter --json data/vbfz_photon_samples.json --inDir ${outdir}/g/qt_loose/ --iLumi 19800 --iEcm 8 --outFile ${outdir}/plotter_dy_closure_g_qt_loose.root --noPlot --forceMerged;
    runPlotter --json data/vbfz_photon_samples.json --inDir ${outdir}/g/qt_pure/  --iLumi 19800 --iEcm 8 --outFile ${outdir}/plotter_dy_closure_g_qt_pure.root  --noPlot --forceMerged;
fi

if [ "$step" == "5" ]; then
    echo "Producing the closure tests report"
fi
