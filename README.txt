### Installation
export SCRAM_ARCH=slc6_amd64_gcc491
scramv1 project CMSSW CMSSW_7_4_2
cd CMSSW_7_4_2/src/
cmsenv

wget -O - --no-check-certificate https://raw.githubusercontent.com/quertenmont/2l2v_fwk/74Xfwk/TAGS.txt | sh

### For creating private samples
# just use official recipe and produce MINIAOD samples in 74X

### Analysing
# see example files in UserCode/llvv_fwk/test/hzz2l2nu/
# Executing the file 'submit.sh' will run the analysis on all samples defined
# in 'samples.json'.  Output files can be recombined in one unique summary
# root file using 'plot.sh'.  The same script will also produce all
# meaningfull plots for the histogram saved by your analysis


### Analysis specific
# If you need set your stuff under test/analysis_subdir

# In general to run an executable or script you can do:
# runLocalAnalysisOverSamples.py -e my_exe -j data/my_samples.json -d my_input_dir -o my_output_dir -c test/runAnalysis_cfg.py.templ -p "@par1=val1 @par2=val2" -s queue
# The easiest is however to create a submit.sh script in the directory of your
# analysis
# Nex executable can also be added in the UserCode/llvv_fwk/test/hzz2l2nu/bin
# directory (also add it to the Buildfile there)
