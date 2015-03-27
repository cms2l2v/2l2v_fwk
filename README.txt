## Installation

# Phys14
cmsrel CMSSW_7_2_2 #released used for PHYS14 miniAOD samples  #OUTDATED
cd CMSSW_7_2_2/src/
cmsenv
git clone git@github.com:veelken/SVfit_standalone TauAnalysis/SVfitStandalone
git clone git@github.com:quertenmont/2l2v_fwk.git UserCode/llvv_fwk
cd UserCode/llvv_fwk
git checkout remotes/origin/72Xfwk #checkout the branch we are interested in
git checkout -b 72Xfwk_modified #copy the branch to a new one to host future modifications (ease pullrequest and code merging)
cd ../..


export SCRAM_ARCH=slc5_amd64_gcc462
# For CVS access
# export CVSROOT :ext:YOURUSERNAME@lxplus.cern.ch:/afs/cern.ch/user/c/cvscmssw/public/CMSSW
scramv1 project CMSSW CMSSW_5_3_15
cd CMSSW_5_3_15/src/
cmsenv

wget -q -O - --no-check-certificate https://raw.github.com/cms2l2v/2l2v_fwk/master/TAGS.txt | sh

## Ntuples

### Creating

Crab files are found in test/grid. To create locally at CERN:

runOverSamples.py -j data/vbfz_samples.json -R "tmp>5000" -n 1 -d aoddir -p "-castor=/afs/cern.ch/user/p/psilva/work/ntuples -cfg=$CMSSW_BASE/src/UserCode/llvv_fwk/test/runDataAnalyzer_mc_cfg.py" -t lljj -s 8nh

### Analysing

5311: /store/cmst3/user/psilva/5311_ntuples


## Analysis specific

If you need set your stuff under test/analysis_subdir

In general to run an executable or script you can do:

runLocalAnalysisOverSamples.py -e my_exe -j data/my_samples.json -d my_input_dir -o my_output_dir -c test/runAnalysis_cfg.py.templ -p "@par1=val1 @par2=val2" -s queue

