if [[ $# -eq 0 ]]; then 
   echo "no argument provided, rerun the command with the step you wish to run"
   echo "possible options are:"
   echo "   " $0 "0; # completely clean up the directory (all data will be lost)"
   echo "   " $0 "1; # run 'runHZZ2l2vAnalysis' on all samples defined in samples.json'"
   echo "   " $0 "2; # compute the integrated luminosity from the processed data samples"
   echo "   " $0 "3; # make plots and combined root files"
   exit 1
fi
step=$1   #variable that store the analysis step to run

#Additional arguments to take into account
arguments=''; for var in "${@:2}"; do arguments=$arguments" "$var; done
if [[ $# -ge 2 ]]; then echo "Additional arguments will be considered: "$arguments ;fi 


#GLOBAL VARIABLES
SUFFIX=_2015_09_17
MAINDIR=$CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v
JSON=$MAINDIR/samples.json
RESULTSDIR=$MAINDIR/results$SUFFIX
PLOTSDIR=$MAINDIR/plots$SUFFIX
PLOTTER=$MAINDIR/plotter$SUFFIX
echo $RESULTSDIR

#GLOBAL CODE

#analysis cleanup
if [[ $step -eq 0 ]]; then
   echo "CLEANING UP"
   rm -rdf $RESULTSDIR $PLOTSDIR LSFJOB_* core.* *.sh.e* *.sh.o*
fi   

#submit jobs for 2l2v analysis
if [[ $step -eq 1 ]]; then
   echo "JOB SUBMISSION"
   queue='8nh'
   if [[ $arguments == *"crab3"* ]]; then queue='crab3' ;fi #IF CRAB3 is provided in argument, use crab submissiong instead of condor/lsf
   runAnalysisOverSamples.py -e runHZZ2l2vAnalysis -j $JSON -o $RESULTSDIR  -c $MAINDIR/../runAnalysis_cfg.py.templ -p "@useMVA=True @saveSummaryTree=True @runSystematics=True @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0" -s $queue --report True $arguments
fi

#extract integrated luminosity of the processed lumi blocks
if [[ $step -eq 2 ]]; then
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
   brilcalc lumi -i $RESULTSDIR/json_all.json -n 0.962 -u /pb -o $RESULTSDIR/LUMI.txt
   tail -n 3 $RESULTSDIR/LUMI.txt  
fi

if [[ $step -eq 3 ]]; then
   INTLUMI=`tail -n 1 $RESULTSDIR/LUMI.txt | cut -d ',' -f 6`
   echo "MAKE PLOTS AND SUMMARY ROOT FILE, BASED ON AN INTEGRATED LUMINOSITY OF $INTLUMI"
   runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outDir $PLOTSDIR/ --outFile $PLOTTER.root  --json $JSON --no2D $arguments
   ln -s $PLOTTER.root $MAINDIR/plotter.root
fi
