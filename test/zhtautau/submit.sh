#!/usr/bin/env bash

#--------------------------------------------------
# Global Code 
#--------------------------------------------------

if [[ $# -eq 0 ]]; then 
    printf "NAME\n\tsubmit.sh - Main driver to submit jobs\n"
    printf "\nSYNOPSIS\n"
    printf "\n\t%-5s\n" "./submit.sh [OPTION]" 
    printf "\nOPTIONS\n" 
    printf "\n\t%-5s  %-40s\n"  "0"  "completely clean up the directory" 
    printf "\n\t%-5s  %-40s\n"  "1"  "run 'runZHTauTauAnalysis' on samples.json" 
    printf "\n\t%-5s  %-40s\n"  "1.1"  "run 'runZHTauTauAnalysis' on photon_samples.json" 
    printf "\n\t%-5s  %-40s\n"  "2"  "compute integrated luminosity from processed samples" 
    printf "\n\t%-5s  %-40s\n"  "3"  "make plots and combine root files" 
    printf "\n\t%-5s  %-40s\n"  "3.1"  "make plots for photon_smaples" 
fi

step=$1   #variable that store the analysis step to run

#Additional arguments to take into account
arguments=''; for var in "${@:2}"; do arguments=$arguments" "$var; done
if [[ $# -ge 4 ]]; then echo "Additional arguments will be considered: "$arguments ;fi 

#--------------------------------------------------
# Global Variables
#--------------------------------------------------
SUFFIX=_2016_03_21
#SUFFIX=$(date +"_%Y_%m_%d") 
MAINDIR=$CMSSW_BASE/src/UserCode/llvv_fwk/test/zhtautau
JSON=$MAINDIR/samples.json
RESULTSDIR=$MAINDIR/results$SUFFIX
PLOTSDIR=$MAINDIR/plots${SUFFIX}
PLOTTER=$MAINDIR/plotter${SUFFIX}

#printf "Result dir is set as: \n\t%s\n" "$RESULTSDIR"


case $step in
    
    0)  #analysis cleanup
	echo "ALL DATA WILL BE LOST! [N/y]?"
	read answer
	if [[ $answer == "y" ]];
	then
	    echo "CLEANING UP..."
	    rm -rdf $RESULTSDIR $PLOTSDIR LSFJOB_* core.* *.sh.e* *.sh.o*
	fi
	;;

    1)  #submit jobs for 2l2v analysis
	echo "JOB SUBMISSION"
	echo "Input: " $JSON
	echo "Output: " $RESULTSDIR

	queue='8nh'
	#IF CRAB3 is provided in argument, use crab submissiong instead of condor/lsf
	if [[ $arguments == *"crab3"* ]]; then queue='crab3' ;fi 
	runAnalysisOverSamples.py -e runZHTauTauAnalysis -j $JSON -o $RESULTSDIR  -c $MAINDIR/../runAnalysis_cfg.py.templ -p "@useMVA=True @saveSummaryTree=True @runSystematics=True @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0" -s $queue --report True $arguments
	;;

    1.1)  #submit jobs for 2l2v photon jet analysis
	echo "JOB SUBMISSION for Photon + Jet analysis"
	queue='8nh'
	JSON=$MAINDIR/photon_samples.json
	echo "Input: " $JSON
	echo "Output: " $RESULTSDIR
	if [[ $arguments == *"crab3"* ]]; then queue='crab3' ;fi 
	runAnalysisOverSamples.py -e runZHTauTauAnalysis -j $JSON -o $RESULTSDIR -c $MAINDIR/../runAnalysis_cfg.py.templ -p "@useMVA=True @saveSummaryTree=True @runSystematics=True @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0" -s $queue --report True $arguments
	;;

    2)  #extract integrated luminosity of the processed lumi blocks
	echo "MISSING LUMI WILL APPEAR AS DIFFERENCE LUMI ONLY IN in.json"
	mergeJSON.py --output=$RESULTSDIR/json_all.json        $RESULTSDIR/Data*.json
	mergeJSON.py --output=$RESULTSDIR/json_doubleMu.json   $RESULTSDIR/Data*_DoubleMu*.json
	mergeJSON.py --output=$RESULTSDIR/json_doubleEl.json   $RESULTSDIR/Data*_DoubleElectron*.json
	mergeJSON.py --output=$RESULTSDIR/json_muEG.json   $RESULTSDIR/Data*_MuEG*.json
	mergeJSON.py --output=$RESULTSDIR/json_in.json  Cert_*.txt
	echo "MISSING LUMI BLOCKS IN DOUBLE MU DATASET"
	compareJSON.py --diff $RESULTSDIR/json_in.json $RESULTSDIR/json_doubleMu.json 
	echo "MISSING LUMI BLOCKS IN DOUBLE ELECTRON DATASET"
	compareJSON.py --diff $RESULTSDIR/json_in.json $RESULTSDIR/json_doubleEl.json 
	echo "MISSING LUMI BLOCKS IN MUON EGAMMA DATASET"
	compareJSON.py --diff $RESULTSDIR/json_in.json $RESULTSDIR/json_muEG.json 
	
	echo "COMPUTE INTEGRATED LUMINOSITY"
	export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.0.3/bin:$PATH
	pip install --install-option="--prefix=$HOME/.local" brilws &> /dev/null #will be installed only the first time
	brilcalc lumi --normtag ~lumipro/public/normtag_file/OfflineNormtagV2.json -i $RESULTSDIR/json_all.json -u /pb -o $RESULTSDIR/LUMI.txt 
	tail -n 3 $RESULTSDIR/LUMI.txt  
	;;

    3)  # make plots and combined root files
        if [ -f $RESULTSDIR/LUMI.txt ]; then
           INTLUMI=`tail -n 1 $RESULTSDIR/LUMI.txt | cut -d ',' -f 6`
        else
           INTLUMI=2268.759 #correspond to the value from DoubleMu OR DoubleEl OR MuEG without jobs failling and golden JSON
           echo "WARNING: $RESULTSDIR/LUMI.txt file is missing so use fixed integrated luminosity value, this might be different than the dataset you ran on"
        fi
	echo "MAKE PLOTS AND SUMMARY ROOT FILE, BASED ON AN INTEGRATED LUMINOSITY OF $INTLUMI"
	runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outDir $PLOTSDIR/ --outFile $PLOTTER.root  --json $JSON --no2D $arguments
	ln -s -f $PLOTTER.root $MAINDIR/plotter.root
	;;

    4)  # make plots and combined root files
        extractFRWeights --inFile $PLOTTER.root --outFile FR_Weights.root --outDir $PLOTSDIR/FakeRateWeights

esac

