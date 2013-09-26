###
Charged Higgs Analysis
###

Mailto: pietro.vischia at gmail.com

Analyze
-------

runLocalAnalysisOverSamples.py -e runChHiggsAnalysis -j data/ch-higgs_samples.json -d /afs/cern.ch/work/v/vischia/private/store/539_ntuples/ -o /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_b/ -c test/runAnalysis_cfg.py.templ -p "@runSystematics=True @saveSummaryTree=False" -s 8nh

runLocalAnalysisOverSamples.py -e runChHiggsAnalysis -j data/ch-higgs_samples.json -d /afs/cern.ch/work/v/vischia/private/store/5311_ntuples/ -o /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_b/ -c test/runAnalysis_cfg.py.templ -p "@runSystematics=True @saveSummaryTree=False" -s 8nh
runLocalAnalysisOverSamples.py -e runChHiggsAnalysis -j data/top_samples.json      -d /store/cmst3/user/psilva/5311_ntuples/             -o /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_b/     -c test/runAnalysis_cfg.py.templ -p "@runSystematics=True @saveSummaryTree=False" -s 8nh


Plots & tables
--------------

runPlotter --iLumi 19702 --inDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311/ --outDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311/plots --json data/plot-ch-higgs_samples.json --outFile /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311/plotter-with-ch-higgs.root --showUnc --plotExt .png
runPlotter --iLumi 19702 --inDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311/ --outDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311/plots --json data/top_samples.json --outFile /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311/plotter-smonly.root --showUnc --plotExt .png  

--onlyStartWith emu_evtflow

Tables only
-----------
runPlotter --iLumi 19702 --inDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311/ --outDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311/plots --json data/ch-higgs_samples.json --outFile /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311/plotter-all-ch-higgs.root --noPlot --noPowers 
runPlotter --iLumi 19702 --inDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311/ --outDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311/plots --json data/top_samples.json --outFile /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311/plotter-sm-only.root --noPlot --noPowers 



Datacards
---------
prepareChHiggsDatacards  --in /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311/plotter-with-ch-higgs.root --out /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311/datacards/ --json data/plot-ch-higgs_samples.json --noPowers --histo evtflow --bin 5  


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
