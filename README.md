# Installation
```bash 
export SCRAM_ARCH=slc6_amd64_gcc493
scramv1 project CMSSW CMSSW_7_6_3_patch2
cd CMSSW_7_6_3_patch2/src/
cmsenv
wget -O - --no-check-certificate https://raw.githubusercontent.com/cms2l2v/2l2v_fwk/master/TAGS.txt | sh
```

# For developers

We have decided to use pull-request mode for the master development.

- Fork the code with your personal github ID. See [details](https://help.github.com/articles/fork-a-repo/)
- Make a clean git clone in the UserCode directory
```
cd path/to/CMSSW_7_6_3_patch2/src/UserCode 
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

