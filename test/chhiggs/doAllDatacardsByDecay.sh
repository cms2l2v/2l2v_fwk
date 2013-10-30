#!/bin/bash


# Do input 
# runPlotter --iLumi 19702 --inDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_b/ --outDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_b/plots --json data/allChsamples.json --outFile /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_b/plotter-all-samplesForDatacards_evtflow.root --noPlot --noPowers  --onlyStartWith emu_evtflow &

# runPlotter --iLumi 19702 --inDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_b/ --outDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_b/plots --json data/allChsamples.json --outFile /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_b/plotter-all-samplesForDatacards_finalevtflow.root --noPlot --noPowers  --onlyStartWith emu_finalevtflow &


### runPlotter --iLumi 19702 --inDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_b/ --outDir /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_d/plots --json data/chhiggs/all-samples_higgs1pb.json --outFile /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_d/plotter-all-samplesForDatacards_finalevtflow_norm.root --noPlot --noPowers  --onlyStartWith emu_finalevtflow 


mv datacardsByDecayNew datacardsByDecayNew_bak


mkdir -p datacardsByDecayNew/180
mkdir -p datacardsByDecayNew/200 
mkdir -p datacardsByDecayNew/220 
mkdir -p datacardsByDecayNew/250 
mkdir -p datacardsByDecayNew/300
mkdir -p datacardsByDecayNew/350
mkdir -p datacardsByDecayNew/400 
mkdir -p datacardsByDecayNew/500 
mkdir -p datacardsByDecayNew/600 
mkdir -p datacardsByDecayNew/700

mkdir -p datacardsByDecaySyst/180
mkdir -p datacardsByDecaySyst/200 
mkdir -p datacardsByDecaySyst/220 
mkdir -p datacardsByDecaySyst/250 
mkdir -p datacardsByDecaySyst/300
mkdir -p datacardsByDecaySyst/350
mkdir -p datacardsByDecaySyst/400 
mkdir -p datacardsByDecaySyst/500 
mkdir -p datacardsByDecaySyst/600 
mkdir -p datacardsByDecaySyst/700

mkdir -p datacardsByDecaySystPAS/180
mkdir -p datacardsByDecaySystPAS/200 
mkdir -p datacardsByDecaySystPAS/220 
mkdir -p datacardsByDecaySystPAS/250 
mkdir -p datacardsByDecaySystPAS/300
mkdir -p datacardsByDecaySystPAS/350
mkdir -p datacardsByDecaySystPAS/400 
mkdir -p datacardsByDecaySystPAS/500 
mkdir -p datacardsByDecaySystPAS/600 
mkdir -p datacardsByDecaySystPAS/700

# latest is in chhiggs directory, plotter-all-samplesForDatacards.root 
# 
for i in 180 200 220 250 300 350 500 600 700
	do
#  prepareChHiggsDatacards --in /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_b/plotter-all-samplesForDatacards_finalevtflow.root --out datacards/${i}/ --json tempjson/${i}samples.json --noPowers --histo finalevtflow1 --bin 1 --ch emu &
#  prepareChHiggsDatacards --in /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_b/plotter-all-samplesForDatacards_finalevtflow.root --out datacards/${i}/ --json tempjson/${i}samples.json --noPowers --histo finalevtflow2 --bin 1 --ch emu &
#  prepareChHiggsDatacards --in /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_b/plotter-all-samplesForDatacards_finalevtflow.root --out datacards/${i}/ --json tempjson/${i}samples.json --noPowers --histo finalevtflow3 --bin 1 --ch emu &
#  prepareChHiggsDatacards --in /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_b/plotter-all-samplesForDatacards_finalevtflow.root --out datacards/${i}/ --json tempjson/${i}samples.json --noPowers --histo finalevtflow4 --bin 1 --ch emu &
  

  prepareChHiggsDatacards --in /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_e/plotter-all-samplesForDatacards_finalevtflow_norm_def.root --out datacardsByDecayNew/${i}/ --suffix tb --json /afs/cern.ch/work/v/vischia/private/results/HIG-13-026/tempjsonByFinalState/${i}_tb.json --noPowers --histo finalevtflow2btags --bin 1 --ch emu & 
  prepareChHiggsDatacards --in /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_e/plotter-all-samplesForDatacards_finalevtflow_norm_def.root --out datacardsByDecayNew/${i}/ --suffix taunu --json /afs/cern.ch/work/v/vischia/private/results/HIG-13-026/tempjsonByFinalState/${i}_taunu.json --noPowers --histo finalevtflow2btags --bin 1 --ch emu & 

  prepareChHiggsDatacards --in /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_e/plotter-forSystTable_def.root --out datacardsByDecaySyst/${i}/ --suffix tb --json /afs/cern.ch/work/v/vischia/private/results/HIG-13-026/tempjsonByFinalState/${i}_tb.json --noPowers --histo finalevtflow2btags --bin 1 --ch emu & 
  prepareChHiggsDatacards --in /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_e/plotter-forSystTable_def.root --out datacardsByDecaySyst/${i}/ --suffix taunu --json /afs/cern.ch/work/v/vischia/private/results/HIG-13-026/tempjsonByFinalState/${i}_taunu.json --noPowers --histo finalevtflow2btags --bin 1 --ch emu & 

  prepareChHiggsDatacards --in /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_e/plotter-forSystTableInPAS_def.root --out datacardsByDecaySystPAS/${i}/ --suffix tb --json /afs/cern.ch/work/v/vischia/private/results/HIG-13-026/tempjsonByFinalState/${i}_tb.json --noPowers --histo evtflow --bin 1 --ch emu & 
  prepareChHiggsDatacards --in /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_e/plotter-forSystTableInPAS_def.root --out datacardsByDecaySystPAS/${i}/ --suffix taunu --json /afs/cern.ch/work/v/vischia/private/results/HIG-13-026/tempjsonByFinalState/${i}_taunu.json --noPowers --histo evtflow --bin 1 --ch emu & 


   #  prepareChHiggsDatacards --in /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_c/plotter-datacards.root --out datacardsByDecayNew/${i}/ --suffix taunu --json tempjsonByFinalState/${i}_taunu.json --noPowers --histo finalevtflow2btags --bin 1 --ch emu &


#  prepareChHiggsDatacards --in /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_b/plotter-all-samplesForDatacards_finalevtflow_norm.root --out datacardsByDecay/${i}/ --suffix tb --json tempjsonByFinalState/${i}_tb.json --noPowers --histo finalevtflow2 --bin 1 --ch emu & 
#  prepareChHiggsDatacards --in /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_b/plotter-all-samplesForDatacards_finalevtflow_norm.root --out datacardsByDecay/${i}/ --suffix tb --json tempjsonByFinalState/${i}_tb.json --noPowers --histo finalevtflow3 --bin 1 --ch emu & 
#  prepareChHiggsDatacards --in /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_b/plotter-all-samplesForDatacards_finalevtflow_norm.root --out datacardsByDecay/${i}/ --suffix tb --json tempjsonByFinalState/${i}_tb.json --noPowers --histo finalevtflow4 --bin 1 --ch emu & 
#  prepareChHiggsDatacards --in /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_b/plotter-all-samplesForDatacards_finalevtflow_norm.root --out datacardsByDecay/${i}/ --suffix tb --json tempjsonByFinalState/${i}_tb.json --noPowers --histo finalevtflow5 --bin 1 --ch emu & 
#
#  prepareChHiggsDatacards --in /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_b/plotter-all-samplesForDatacards_finalevtflow_norm.root --out datacardsByDecay/${i}/ --suffix taunu --json tempjsonByFinalState/${i}_taunu.json --noPowers --histo finalevtflow2 --bin 1 --ch emu &
#  prepareChHiggsDatacards --in /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_b/plotter-all-samplesForDatacards_finalevtflow_norm.root --out datacardsByDecay/${i}/ --suffix taunu --json tempjsonByFinalState/${i}_taunu.json --noPowers --histo finalevtflow3 --bin 1 --ch emu &
#  prepareChHiggsDatacards --in /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_b/plotter-all-samplesForDatacards_finalevtflow_norm.root --out datacardsByDecay/${i}/ --suffix taunu --json tempjsonByFinalState/${i}_taunu.json --noPowers --histo finalevtflow4 --bin 1 --ch emu &
#  prepareChHiggsDatacards --in /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_b/plotter-all-samplesForDatacards_finalevtflow_norm.root --out datacardsByDecay/${i}/ --suffix taunu --json tempjsonByFinalState/${i}_taunu.json --noPowers --histo finalevtflow5 --bin 1 --ch emu &
#
done
exit




