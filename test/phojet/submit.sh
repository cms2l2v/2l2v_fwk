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

# # Setup grid proxy # no need after crab 3.3.12 
# cp x509_proxy $CMSSW_BASE/
# export X509_USER_PROXY=$CMSSW_BASE/x509_proxy

# Now run the analysis
runPhoJetAnalysis phojet_cfg.py


