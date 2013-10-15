#!/bin/bash

#mkdir -p datacardsByDecayNew/180
#mkdir -p datacardsByDecayNew/200 
#mkdir -p datacardsByDecayNew/220 
#mkdir -p datacardsByDecayNew/250 
#mkdir -p datacardsByDecayNew/300

mkdir -p datacardsByDecayNew/350
mkdir -p datacardsByDecayNew/400 
mkdir -p datacardsByDecayNew/500 
mkdir -p datacardsByDecayNew/600 
mkdir -p datacardsByDecayNew/700

# latest is in chhiggs directory, plotter-all-samplesForDatacards.root 
# 180 200 220 250 300
for i in 350 400 500 600 700
	do
#  prepareChHiggsDatacards --in /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_b/plotter-all-samplesForDatacards_finalevtflow.root --out datacards/${i}/ --json tempjson/${i}samples.json --noPowers --histo finalevtflow1 --bin 1 --ch emu &
#  prepareChHiggsDatacards --in /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_b/plotter-all-samplesForDatacards_finalevtflow.root --out datacards/${i}/ --json tempjson/${i}samples.json --noPowers --histo finalevtflow2 --bin 1 --ch emu &
#  prepareChHiggsDatacards --in /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_b/plotter-all-samplesForDatacards_finalevtflow.root --out datacards/${i}/ --json tempjson/${i}samples.json --noPowers --histo finalevtflow3 --bin 1 --ch emu &
#  prepareChHiggsDatacards --in /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_b/plotter-all-samplesForDatacards_finalevtflow.root --out datacards/${i}/ --json tempjson/${i}samples.json --noPowers --histo finalevtflow4 --bin 1 --ch emu &
  

  prepareChHiggsDatacards --in /afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_c/plotter-datacards.root --out datacardsByDecayNew/${i}/ --suffix tb --json tempjsonByFinalState/${i}_tb.json --noPowers --histo finalevtflow2btags --bin 1 --ch emu & 

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




