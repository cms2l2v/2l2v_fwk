###
Charged Higgs Analysis
###



Analyze
-------
NOTE: run with @runSystematics=False, 'cause the systs are not fully implemented yet, so you would produce N times the same base plots :)


runLocalAnalysisOverSamples.py -e runChHiggsAnalysis -j data/ch-higgs_samples.json -d /afs/cern.ch/work/v/vischia/private/store/539_ntuples/ -o /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs/ -c test/runAnalysis_cfg.py.templ -p "@runSystematics=True @saveSummaryTree=False" -s 8nh
runLocalAnalysisOverSamples.py -e runChHiggsAnalysis -j data/top_samples.json      -d /store/cmst3/user/psilva/Summer13_ntuples/             -o /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs/ -c test/runAnalysis_cfg.py.templ -p "@runSystematics=True @saveSummaryTree=False" -s 8nh


Plots & tables
--------------
runPlotter --iLumi 19683 --inDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs/ --outDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs/plots --json data/plot-ch-higgs_samples.json --outFile /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs/plotter-with-ch-higgs.root --showUnc --plotExt .pdf 
runPlotter --iLumi 19683 --inDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs/ --outDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs/plots --json data/plot-ch-higgs_samples.json --outFile /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs/plotter-with-ch-higgs.root --showUnc --plotExt .png  

--onlyStartWith emu_evtflow


Tables only
-----------
runPlotter --iLumi 19683 --inDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs/ --outDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs/plots --json data/ch-higgs_samples.json --outFile /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs/plotter-all-ch-higgs.root --noPlot --noPowers 
runPlotter --iLumi 19683 --inDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs/ --outDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs/plots --json data/top_samples.json --outFile /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs/plotter-sm-only.root --noPlot --noPowers 



Datacards
---------
prepareChHiggsDatacards  --in /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs/plotter-with-ch-higgs.root --out /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs/datacards/ --json data/plot-ch-higgs_samples.json --noPowers --histo evtflow --bin 5  

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
