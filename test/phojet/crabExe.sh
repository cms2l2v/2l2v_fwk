#!/bin/sh

#  Submit PhoJetAnalysis Jobs
#  Author: Xin Shi <Xin.Shi@cern.ch> 
#  2014.12.12


# Check the current dir and list the files
pwd
ls -lR

# Copy files to lib and bin dir
cp ./lib/$SCRAM_ARCH/* $CMSSW_BASE/lib/$SCRAM_ARCH
cp runPhoJetAnalysis $CMSSW_BASE/bin/$SCRAM_ARCH

# Create the JobRepport
cmsRun -j FrameworkJobReport.xml phojet_cfg.py

# Run the analysis
runPhoJetAnalysis phojet_cfg.py

