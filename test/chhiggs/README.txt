###
Charged Higgs Analysis
###



Analyze
-------
runLocalAnalysisOverSamples.py -e runChHiggsAnalysis -j data/top_samples.json -d /store/cmst3/user/psilva/539_ntuples -o /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/top_539/ -c test/runAnalysis_cfg.py.templ -p "@runSystematics=True @saveSummaryTree=False" -s 8nh
runLocalAnalysisOverSamples.py -e runChHiggsAnalysis -j data/ch-higgs_samples.json -d /afs/cern.ch/work/v/vischia/private/store/539_ntuples -o /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/top_539/ -c test/runAnalysis_cfg.py.templ -p "@runSystematics=True @saveSummaryTree=False" -s 8nh
runLocalAnalysisOverSamples.py -e runChHiggsAnalysis -j data/top_syst_samples.json -d /store/cmst3/user/psilva/539_ntuples -o /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/top_539/ -c test/runAnalysis_cfg.py.templ -p "@runSystematics=True @saveSummaryTree=False"  -s 8nh

Plots & tables
--------------
runPlotter --iLumi 19683 --inDir  /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/top_539/ --outDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/top_539/plots --json data/plot-ch-higgs_samples.json --outFile /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/top_539/plotter-with-ch-higgs.root  --showUnc   --plotExt .pdf 
runPlotter --iLumi 19683 --inDir  /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/top_539/ --outDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/top_539/plots --json data/plot-ch-higgs_samples.json --outFile /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/top_539/plotter-with-ch-higgs.root  --showUnc   --plotExt .png 

Datacards
---------
prepareChHiggsDatacards  --in /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/top_539/plotter-with-ch-higgs.root --out /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/top_539/datacards --syst /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/top_539/plotter-with-ch-higgs.root --json data/plot-ch-higgs_samples.json --noPowers --histo evtflow --bin 5  
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
2013-06-07: decoupled runChHiggsAnalysis from runTopAnalysis, improved cutflow




### 
TODO list
###
- create datacards builder
- finish adding taus to ntuplizer (singleLepton and singleJet triggers, trigger matching)
- import tauDilepton analysis from https://github.com/vischia/TopTaus
