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


#analysis cleanup
if [[ $step -eq 0 ]]; then
   echo "CLEANING UP"
   rm -rdf results plots LSFJOB_* core.* *.sh.e* *.sh.o*
fi   

#submit jobs for 2l2v analysis
if [[ $step -eq 1 ]]; then
   echo "JOB SUBMISSION"
   queue='8nh'
   if [[ $arguments == *"crab3"* ]]; then queue='crab3' ;fi #IF CRAB3 is provided in argument, use crab submissiong instead of condor/lsf
   runAnalysisOverSamples.py -e runHZZ2l2vAnalysis -j $CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v/samples.json -o $CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v/results -d  /store/group/phys_higgs/cmshzz2l2v/2014_03_20/ -c $CMSSW_BASE/src/UserCode/llvv_fwk/test/runAnalysis_cfg.py.templ -p "@useMVA=True @saveSummaryTree=True @runSystematics=True @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0" -s $queue --report True $arguments
fi

#extract integrated luminosity of the processed lumi blocks
if [[ $step -eq 2 ]]; then
   echo "MISSING LUMI WILL APPEAR AS DIFFERENCE LUMI ONLY IN in.json"
   mergeJSON.py --output=processedData.json   results/Data*.json
   mergeJSON.py --output=in.json  Cert_*.txt
   compareJSON.py --diff in.json processedData.json 

   echo "COMPUTE INTEGRATED LUMINOSITY"
   export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.0.3/bin:$PATH
   pip install --install-option="--prefix=$HOME/.local" brilws &> /dev/null #will be installed only the first time
   brilcalc lumi -i processedData.json -n 0.962 -u /pb -o results/LUMI.txt
   tail -n 3 results/LUMI.txt  
fi

if [[ $step -eq 3 ]]; then
   INTLUMI=`tail -n 1 results/LUMI.txt | cut -d ',' -f 6`
   echo "MAKE PLOTS AND SUMMARY ROOT FILE, BASED ON AN INTEGRATED LUMINOSITY OF $INTLUMI"
   runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $PWD/results/ --outDir $PWD/plots/ --outFile $PWD/plotter.root  --json $PWD/samples.json --no2D $arguments
fi
