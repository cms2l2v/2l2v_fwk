###
Installation
###

export SCRAM_ARCH=slc5_amd64_gcc462
scramv1 project CMSSW CMSSW_5_3_9
cd CMSSW_5_3_9/src/
curl -O https://raw.github.com/pfs/usercode/master/TAGS.txt | sh

###
create ntuple
###
Crab files are found in test/grid
Create locally at CERN:
runOverSamples.py -j data/vbfz_samples.json -R "tmp>5000" -n 1 -d aoddir -p "-castor=/afs/cern.ch/user/p/psilva/work/ntuples -cfg=$CMSSW_BASE/src/UserCode/EWKV/test/runDataAnalyzer_mc_cfg.py" -t lljj -s 8nh

###
analysis specific
###
if you need set your stuff under test/analysis_subdir
in general to run an executable or script you can do:
runLocalAnalysisOverSamples.py -e my_exe -j data/my_samples.json -d my_input_dir -o my_output_dir -c test/runAnalysis_cfg.py.templ -p "@par1=val1 @par2=val2" -s queue

