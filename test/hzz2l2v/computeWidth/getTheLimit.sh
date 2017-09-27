#!/bin/sh

computeLimit --m 125 --in /eos/cms/store/user/hbrun/analysis/plotters/plotter_width_21sept.root --syst --index 17,17,17 --bins eq0jets,geq1jets,vbf --json $CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v/samples_widthLimits.json --key 2l2v_datadriven_SM --BackExtrapol --statBinByBin 0.00001 --dropBckgBelow 0.00001 --subNRB --histo mt_shapes --histoVBF mt_shapes --systpostfix _13TeV --shape --skipQQH --signalSufix "_cp1.00_brn0.00" --blindWithSignal

sh combineCards.sh

text2workspace.py card_combined.dat -o workspace.root -P UserCode.llvv_fwk.HiggsWidth2017:higgswidth --PO verbose

combine -n Exp -M MultiDimFit workspace.root -t -1 --expectSignal=1 --algo=grid --points 300 --setPhysicsModelParameterRanges r=0.0,15
python plotMuScan.py
