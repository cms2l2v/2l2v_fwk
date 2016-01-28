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
    printf "\n\t%-5s  %-40s\n"  "1"  "run 'runHZZ2l2vAnalysis' on samples.json" 
    printf "\n\t%-5s  %-40s\n"  "1.1"  "run 'runHZZ2l2vAnalysis' on photon_samples.json" 
    printf "\n\t%-5s  %-40s\n"  "1.2"  "run 'runHZZ2l2vAnalysis' on photon_samples.json with photon re-weighting"
    printf "\n\t%-5s  %-40s\n"  "2"  "compute integrated luminosity from processed samples" 
    printf "\n\t%-5s  %-40s\n"  "2.1"  "compute integrated luminosity from processed photon samples" 
    printf "\n\t%-5s  %-40s\n"  "3"  "make plots and combine root files" 
    printf "\n\t%-5s  %-40s\n"  "3.1"  "make plots for photon_samples" 
    printf "\n\t%-5s  %-40s\n"  "3.2"  "make plots for photon_samples with photon re-weighting" 
fi

step=$1   #variable that store the analysis step to run

#Additional arguments to take into account
arguments=''; for var in "${@:2}"; do arguments=$arguments" "$var; done
if [[ $# -ge 4 ]]; then echo "Additional arguments will be considered: "$arguments ;fi 

#--------------------------------------------------
# Global Variables
#--------------------------------------------------
SUFFIX=_2016_01_27
#SUFFIX=_debug
#SUFFIX=$(date +"_%Y_%m_%d") 
MAINDIR=$CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v
JSON=$MAINDIR/samples.json
RESULTSDIR=$MAINDIR/results$SUFFIX
PLOTSDIR=$MAINDIR/plots${SUFFIX}
PLOTTER=$MAINDIR/plotter${SUFFIX}

echo "Input: " $JSON
echo "Output: " $RESULTSDIR

queue='8nh'
#IF CRAB3 is provided in argument, use crab submissiong instead of condor/lsf
if [[ $arguments == *"crab3"* ]]; then queue='crab3' ;fi 


################################################# STEPS between 0 and 1
if [[ $step == 0 ]]; then   
        #analysis cleanup
	echo "ALL DATA WILL BE LOST! [N/y]?"
	read answer
	if [[ $answer == "y" ]];
	then
	    echo "CLEANING UP..."
	    rm -rdf $RESULTSDIR $PLOTSDIR LSFJOB_* core.* *.sh.e* *.sh.o*
	fi
fi #end of step0

###  ############################################## STEPS between 1 and 2
if [[ $step > 0.999 &&  $step < 2 ]]; then
   if [[ $step == 1 || $step == 1.0 ]]; then        #submit jobs for 2l2v analysis
	echo "JOB SUBMISSION"
	runAnalysisOverSamples.py -e runHZZ2l2vAnalysis -j $JSON -o $RESULTSDIR  -c $MAINDIR/../runAnalysis_cfg.py.templ -p "@useMVA=True @saveSummaryTree=True @runSystematics=True @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0" -s $queue --report True --key 2l2v_mcbased $arguments
   fi

   if [[ $step == 1 || $step == 1.1 ]]; then        #submit jobs for 2l2v photon jet analysis
	echo "JOB SUBMISSION for Photon + Jet analysis"
	runAnalysisOverSamples.py -e runHZZ2l2vAnalysis -j $JSON -o $RESULTSDIR -c $MAINDIR/../runAnalysis_cfg.py.templ -p "@useMVA=True @saveSummaryTree=True @runSystematics=True @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0" -s $queue --report True --key 2l2v_photoncontrol $arguments 
   fi

   if [[ $step == 1 || $step == 1.2 ]]; then        #submit jobs for 2l2v analysis with photon re-weighting
	echo "JOB SUBMISSION for Photon + Jet analysis with photon re-weighting"                                                                                   
        runAnalysisOverSamples.py -e runHZZ2l2vAnalysis -j $JSON -o $RESULTSDIR -c $MAINDIR/../runAnalysis_cfg.py.templ -p "@useMVA=True @saveSummaryTree=True @weightsFile=$PWD/photonWeights_RunD.root @runSystematics=True @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0" -s $queue --report True --key 2l2v_photonsOnly $arguments 
   fi                              
fi

###  ############################################## STEPS between 2 and 3
if [[ $step > 1.999 && $step < 3 ]]; then
   if [[ $step == 2 ]]; then    #extract integrated luminosity of the processed lumi blocks
	echo "MISSING LUMI WILL APPEAR AS DIFFERENCE LUMI ONLY IN in.json"
	mergeJSON.py --output=$RESULTSDIR/json_all.json        $RESULTSDIR/Data*.json
	mergeJSON.py --output=$RESULTSDIR/json_doubleMu.json   $RESULTSDIR/Data*_DoubleMu*.json
	mergeJSON.py --output=$RESULTSDIR/json_doubleEl.json   $RESULTSDIR/Data*_DoubleElectron*.json
	mergeJSON.py --output=$RESULTSDIR/json_muEG.json       $RESULTSDIR/Data*_MuEG*.json
        mergeJSON.py --output=$RESULTSDIR/json_gamma.json      $RESULTSDIR/Data13TeV_SinglePhoton*.json
	mergeJSON.py --output=$RESULTSDIR/json_in.json  Cert_*.txt
	echo "MISSING LUMI BLOCKS IN DOUBLE MU DATASET"
	compareJSON.py --diff $RESULTSDIR/json_in.json $RESULTSDIR/json_doubleMu.json 
	echo "MISSING LUMI BLOCKS IN DOUBLE ELECTRON DATASET"
	compareJSON.py --diff $RESULTSDIR/json_in.json $RESULTSDIR/json_doubleEl.json 
	echo "MISSING LUMI BLOCKS IN MUON EGAMMA DATASET"
	compareJSON.py --diff $RESULTSDIR/json_in.json $RESULTSDIR/json_muEG.json 
	echo "MISSING LUMI BLOCKS IN Single Photon DATASET"
	compareJSON.py --diff $RESULTSDIR/json_in.json $RESULTSDIR/json_gamma.json 

	echo "COMPUTE INTEGRATED LUMINOSITY"
	export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.0.3/bin:$PATH
	pip install --upgrade --install-option="--prefix=$HOME/.local" brilws &> /dev/null #will be installed only the first time
	brilcalc lumi --normtag ~lumipro/public/normtag_file/OfflineNormtagV2.json -i $RESULTSDIR/json_all.json -u /pb -o $RESULTSDIR/LUMI.txt 
	tail -n 3 $RESULTSDIR/LUMI.txt  
     fi
  fi     

###  ############################################## STEPS between 3 and 4
if [[ $step > 2.999 && $step < 4 ]]; then
    if [ -f $RESULTSDIR/LUMI.txt ]; then
      INTLUMI=`tail -n 1 $RESULTSDIR/LUMI.txt | cut -d ',' -f 6`
    else
      INTLUMI=2215.182 #correspond to the value from DoubleMu OR DoubleEl OR MuEG without jobs failling and golden JSON
      echo "WARNING: $RESULTSDIR/LUMI.txt file is missing so use fixed integrated luminosity value, this might be different than the dataset you ran on"
    fi

     if [[ $step == 3 || $step == 3.0 ]]; then  # make plots and combined root files
	echo "MAKE SUMMARY ROOT FILE, BASED ON AN INTEGRATED LUMINOSITY OF $INTLUMI"
#        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outFile ${PLOTTER}.root  --json $JSON --noPlot --fileOption RECREATE --key 2l2v_mcbased                     $arguments        
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outFile ${PLOTTER}.root  --json $JSON --noPlot --fileOption UPDATE   --key 2l2v_mcbased --doInterpollation  $arguments        
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outFile ${PLOTTER}.root  --json $JSON --noPlot --fileOption UPDATE   --key 2l2v_photoncontrol               $arguments 
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outFile ${PLOTTER}.root  --json $JSON --noPlot --fileOption UPDATE   --key 2l2v_photonsOnly                 $arguments
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outFile ${PLOTTER}.root  --json $JSON --noPlot --fileOption UPDATE   --key 2l2v_datadriven                  $arguments 
	ln -s -f ${PLOTTER}.root plotter.root 
    fi        
   

    if [[ $step == 3 || $step == 3.1 ]]; then  # make plots and combine root files for photon + jet study    
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outDir $PLOTSDIR/mcbased/       --outFile ${PLOTTER}.root  --json $JSON --no2D --key 2l2v_mcbased --fileOption READ $arguments
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outDir $PLOTSDIR/mcbased_blind/ --outFile ${PLOTTER}.root  --json $JSON --no2D --key 2l2v_mcbased --fileOption READ --blind 100 --only "(all|ll|mumu|ee|emu)(|eq0jets|geq1jets|vbf)_(met|metpuppi)" $arguments
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outDir $PLOTSDIR/mcbased_blind/ --outFile ${PLOTTER}.root  --json $JSON --no2D --key 2l2v_mcbased --fileOption READ --blind 325 --only "(all|ll|mumu|ee|emu)(|eq0jets|geq1jets|vbf)_mt" $arguments
      
    fi

    if [[ $step == 3 || $step == 3.15 ]]; then  # make plots and combine root files for photon + jet study    
	echo "MAKE PLOTS AND SUMMARY ROOT FILE for Photon sample"
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outDir $PLOTSDIR/photons/ --outFile ${PLOTTER}.root  --json $JSON --no2D --key 2l2v_photoncontrol --fileOption READ --only "(all|ll|mumu|ee|emu|gamma)(_qt|_qmass|_met|_mt)"  $arguments 
    fi

    if [[ $step == 3 || $step == 3.2 ]]; then # make plots and combine root files for photon + jet study    
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outDir $PLOTSDIR/datadriven/       --outFile ${PLOTTER}.root  --json $JSON --no2D --key 2l2v_datadriven --fileOption READ $arguments
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outDir $PLOTSDIR/datadriven_blind/ --outFile ${PLOTTER}.root  --json $JSON --no2D --key 2l2v_datadriven --fileOption READ --blind 100 --only "(all|ll|mumu|ee|emu)(|eq0jets|geq1jets|vbf)_(met|metpuppi)" $arguments
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outDir $PLOTSDIR/datadriven_blind/ --outFile ${PLOTTER}.root  --json $JSON --no2D --key 2l2v_datadriven --fileOption READ --blind 325 --only "(all|ll|mumu|ee|emu)(|eq0jets|geq1jets|vbf)_mt" $arguments
    fi
fi

