# Installation for 76X (2015)
```bash 
export SCRAM_ARCH=slc6_amd64_gcc493
scramv1 project CMSSW CMSSW_7_6_3_patch2
cd CMSSW_7_6_3_patch2/src/
cmsenv
wget -O - --no-check-certificate https://raw.githubusercontent.com/cms2l2v/2l2v_fwk/master/TAGS.txt | sh
```

# Installation for 80X (2016)
```bash 
export SCRAM_ARCH=slc6_amd64_gcc530
cmsrel CMSSW_8_0_14
cd CMSSW_8_0_14/src/
cmsenv

# The following fail for the moment:
#git cms-merge-topic -u matteosan1:smearer_76X

# This fails also: (see https://hypernews.cern.ch/HyperNews/CMS/get/sw-develtools/2414.html )
#git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
#cd HiggsAnalysis/CombinedLimit
#git checkout remotes/origin/74x-root6 
#cd ../..

git clone -b svFit_2015Apr03 https://github.com/veelken/SVfit_standalone.git TauAnalysis/SVfitStandalone

git clone https://github.com/cms2l2v/2l2v_fwk.git UserCode/llvv_fwk
cd UserCode/llvv_fwk
git checkout -b modified #copy the branch to a new one to host future modifications (ease pull request and code merging)
cd ../..

#since HiggsAnalysis is broken in 80X, copy one of the missing file to allow compilation and change some paths
wget https://raw.githubusercontent.com/cms-analysis/HiggsAnalysis-CombinedLimit/74x-root6/src/th1fmorph.cc -P UserCode/llvv_fwk/src/
wget https://raw.githubusercontent.com/cms-analysis/HiggsAnalysis-CombinedLimit/74x-root6/interface/th1fmorph.h -P UserCode/llvv_fwk/interface/
find UserCode/llvv_fwk/ -type f -name '*.cc' -exec sed -i -e 's/HiggsAnalysis\/CombinedLimit\/interface\/th1fmorph.h/UserCode\/llvv_fwk\/interface\/th1fmorph.h/g' {} \;

scramv1 b -j 16
```

# An important note about PR in 80X (2016)
Please before doing your PR, be sure to:
 - test your changes
 - reverse the changes due to th1fmorph thanks to the following command (and be sure to not add those new files in your git add):
```bash 
#TO DO before a PR
cd $CMSSW_BASE/src
find UserCode/llvv_fwk/ -type f -name '*.cc' -exec sed -i -e 's/UserCode\/llvv_fwk\/interface\/th1fmorph.h/HiggsAnalysis\/CombinedLimit\/interface\/th1fmorph.h/g' {} \;
```
 - then you can make your PR
 - if you want to continue developing, don't forget to reverse this once again:
```bash 
#TO DO after a PR to be able to continue to work
cd $CMSSW_BASE/src
find UserCode/llvv_fwk/ -type f -name '*.cc' -exec sed -i -e 's/HiggsAnalysis\/CombinedLimit\/interface\/th1fmorph.h/UserCode\/llvv_fwk\/interface\/th1fmorph.h/g' {} \;
```


# For developers

We have decided to use pull-request mode for the master development.

- Fork the code with your personal github ID. See [details](https://help.github.com/articles/fork-a-repo/)
- Make a clean git clone in the UserCode directory
```
cd $CMSSW_BASE/src/UserCode 
git clone git@github.com:yourgithubid/2l2v_fwk.git llvv_fwk
cd llvv_fwk
git remote add upstream git@github.com:cms2l2v/2l2v_fwk.git
git remote update
git merge upstream/master
```
- Make your own change and commit
```
git commit -a -m "Added feature A, B, C"
git push
```
- Make a pull request agains the cms2l2v. See [details](https://help.github.com/articles/using-pull-requests/)


# For creating private samples
Just use official recipe and produce MINIAOD samples in 76X


# Analysing
See example files in UserCode/llvv_fwk/test/hzz2l2nu/ Executing the
file 'submit.sh' will run the analysis on all samples defined in
'samples.json'.

Please see the full usage of 'submit.sh' with

> ./submit.sh

# Analysis specific
If you need set your stuff under test/analysis_subdir

In general to run an executable or script you can do:
```
runLocalAnalysisOverSamples.py -e my_exe -j data/my_samples.json -d my_input_dir -o my_output_dir -c test/runAnalysis_cfg.py.templ -p "@par1=val1 @par2=val2" -s queue
```

 The easiest is however to create a submit.sh script in the directory
 of your analysis.  The executable can also be added in the
 UserCode/llvv_fwk/test/hzz2l2nu/bin directory (also add it to the
 Buildfile there)

