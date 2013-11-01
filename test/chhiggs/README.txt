###
Charged Higgs Analysis
###

Mailto: pietro.vischia at gmail.com

Analyze
-------
cd $CMSSW_BASE/src/UserCode/llvv_fwk/

# Process ntuples
sh test/chhiggs/doAnal.sh current anal_sus
sh test/chhiggs/doAnal.sh current anal_sm
# Produce plots, tables and datacards
sh test/chhiggs/doAnal.sh current plots
sh test/chhiggs/doAnal.sh current tables
sh test/chhiggs/doAnal.sh current datacards
# Copy plots for putting in AN
sh test/chhiggs/doAnal.sh current put



### FWLite test
runLocalAnalysisOverSamples.py -e runChHiggsAnalysisFWLite -j $CMSSW_BASE/src/UserCode/llvv_fwk/data/zhtautau_samples.json -o /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_fwlite/ -d  /store/group/phys_higgs/cmshzz2l2v/2013_08_30/ -c $CMSSW_BASE/src/UserCode/llvv_fwk/test/runAnalysis_cfg.py.templ -p "@useMVA=True @saveSummaryTree=True @runSystematics=True @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0" -s 8nh


Plots & tables
--------------


### FWLite test
runPlotterFWLite --iEcm 8 --iLumi 19577 --inDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_fwlite/ --outDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_fwlite/plots/ --outFile /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_fwlite/plotter.root  --json $CMSSW_BASE/src/UserCode/llvv_fwk/data/zhtautau_samples.json


Tables only
-----------

Datacards
---------

Options
--out       --> output director
--in        --> input file from plotter
--syst      --> input file with syst shapes
--json      --> json file with the sample descriptor
--histo     --> name of histogram to be used
--noPowers --> Do not use powers of 10 for numbers in tables
--bins      --> list of bins to be used (they must be comma separated without space)
--ch        --> list of channels to be used (they must be comma separated without space)



###
CHANGELOG
###
2013-09-26: Datacards for likelihood fit
2013-09-17: Extended mass range (350, 400, 500, 600, 700) 
2013-08-31: Top pt reweighting and so on
2013-07-22: finished adding samples
2013-07-18: fixed nasty bug
2013-06-21: implemented leff and pu systematics
2013-06-11: added datacards production code. 
	    started implementing systs
	    switched to loose lepton ID 
	    switched to >=2 btags (signal has large multiplicity) 
2013-06-07: decoupled runChHiggsAnalysis from runTopAnalysis, improved cutflow




### 
TODO list
###
- Reimplement LandSShapesProducer (github.com/vischia/TopTaus) w/ different outputs (add a switch for the input interface)
- DY reweighting
- add smart datacards building for different signal masses
- finish adding taus to ntuplizer (singleLepton and singleJet triggers, trigger matching)
- import tauDilepton analysis from https://github.com/vischia/TopTaus
