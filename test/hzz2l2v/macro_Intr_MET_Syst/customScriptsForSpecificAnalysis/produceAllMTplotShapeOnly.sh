mkdir -p AllMTPlots
#ALL
cp SystBySyst/ALL/InstrMET_systematics.root $CMSSW_BASE/src/UserCode/llvv_fwk/data/InstrMET_systematics/.
mkdir -p AllMTPlots/ALL
cd AllMTPlots/ALL
computeLimit --m 800 --in $CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v/plotter_2017_03_21_forLimits.root --syst --index  17,17,17 --bins eq0jets,geq1jets,vbf --json $CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v/samples_full2016_GGH.json --key 2l2v_datadriven_SM    --BackExtrapol --statBinByBin 0.00001 --dropBckgBelow 0.00001  --subNRB --blind  --histo mt_shapes --histoVBF mt_shapes  --systpostfix _13TeV --shape --skipQQH  --signalSufix "_cp100.00_brn0.00"
cd ../..
#WGamma
cp SystBySyst/WGAMMA_ONLY/InstrMET_systematics.root $CMSSW_BASE/src/UserCode/llvv_fwk/data/InstrMET_systematics/.
mkdir -p AllMTPlots/WGAMMA_ONLY
cd AllMTPlots/WGAMMA_ONLY
computeLimit --m 800 --in $CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v/plotter_2017_03_21_forLimits.root --syst --index  17,17,17 --bins eq0jets,geq1jets,vbf --json $CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v/samples_full2016_GGH.json --key 2l2v_datadriven_SM    --BackExtrapol --statBinByBin 0.00001 --dropBckgBelow 0.00001  --subNRB --blind  --histo mt_shapes --histoVBF mt_shapes  --systpostfix _13TeV --shape --skipQQH  --signalSufix "_cp100.00_brn0.00"
cd ../..
#ZGamma
cp SystBySyst/ZGAMMA_ONLY/InstrMET_systematics.root $CMSSW_BASE/src/UserCode/llvv_fwk/data/InstrMET_systematics/.
mkdir -p AllMTPlots/ZGAMMA_ONLY
cd AllMTPlots/ZGAMMA_ONLY
computeLimit --m 800 --in $CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v/plotter_2017_03_21_forLimits.root --syst --index  17,17,17 --bins eq0jets,geq1jets,vbf --json $CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v/samples_full2016_GGH.json --key 2l2v_datadriven_SM    --BackExtrapol --statBinByBin 0.00001 --dropBckgBelow 0.00001  --subNRB --blind  --histo mt_shapes --histoVBF mt_shapes  --systpostfix _13TeV --shape --skipQQH  --signalSufix "_cp100.00_brn0.00"
cd ../..
#Wjets
cp SystBySyst/WJETS_ONLY/InstrMET_systematics.root $CMSSW_BASE/src/UserCode/llvv_fwk/data/InstrMET_systematics/.
mkdir -p AllMTPlots/WJETS_ONLY
cd AllMTPlots/WJETS_ONLY
computeLimit --m 800 --in $CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v/plotter_2017_03_21_forLimits.root --syst --index  17,17,17 --bins eq0jets,geq1jets,vbf --json $CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v/samples_full2016_GGH.json --key 2l2v_datadriven_SM    --BackExtrapol --statBinByBin 0.00001 --dropBckgBelow 0.00001  --subNRB --blind  --histo mt_shapes --histoVBF mt_shapes  --systpostfix _13TeV --shape --skipQQH  --signalSufix "_cp100.00_brn0.00"
cd ../..
#Method
cp SystBySyst/METHOD_ONLY/InstrMET_systematics.root $CMSSW_BASE/src/UserCode/llvv_fwk/data/InstrMET_systematics/.
mkdir -p AllMTPlots/METHOD_ONLY
cd AllMTPlots/METHOD_ONLY
computeLimit --m 800 --in $CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v/plotter_2017_03_21_forLimits.root --syst --index  17,17,17 --bins eq0jets,geq1jets,vbf --json $CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v/samples_full2016_GGH.json --key 2l2v_datadriven_SM    --BackExtrapol --statBinByBin 0.00001 --dropBckgBelow 0.00001  --subNRB --blind  --histo mt_shapes --histoVBF mt_shapes  --systpostfix _13TeV --shape --skipQQH  --signalSufix "_cp100.00_brn0.00"
cd ../..
#Gamma stats
cp SystBySyst/GAMMASTATS_ONLY/InstrMET_systematics.root $CMSSW_BASE/src/UserCode/llvv_fwk/data/InstrMET_systematics/.
mkdir -p AllMTPlots/GAMMASTATS_ONLY
cd AllMTPlots/GAMMASTATS_ONLY
computeLimit --m 800 --in $CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v/plotter_2017_03_21_forLimits.root --syst --index  17,17,17 --bins eq0jets,geq1jets,vbf --json $CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v/samples_full2016_GGH.json --key 2l2v_datadriven_SM    --BackExtrapol --statBinByBin 0.00001 --dropBckgBelow 0.00001  --subNRB --blind  --histo mt_shapes --histoVBF mt_shapes  --systpostfix _13TeV --shape --skipQQH  --signalSufix "_cp100.00_brn0.00"
cd ../..


#Copy important files
cp AllMTPlots/ALL/Yields.tex AllMTPlots/Yields_ALL.tex
cp AllMTPlots/WGAMMA_ONLY/Yields.tex AllMTPlots/Yields_WGAMMA_ONLY.tex
cp AllMTPlots/ZGAMMA_ONLY/Yields.tex AllMTPlots/Yields_ZGAMMA_ONLY.tex
cp AllMTPlots/WJETS_ONLY/Yields.tex AllMTPlots/Yields_WJETS_ONLY.tex
cp AllMTPlots/METHOD_ONLY/Yields.tex AllMTPlots/Yields_METHOD_ONLY.tex
cp AllMTPlots/GAMMASTATS_ONLY/Yields.tex AllMTPlots/Yields_GAMMASTATS_ONLY.tex

cp AllMTPlots/ALL/plot_Shape.png AllMTPlots/plot_Shape_ALL.png
cp AllMTPlots/WGAMMA_ONLY/plot_Shape.png AllMTPlots/plot_Shape_WGAMMA_ONLY.png
cp AllMTPlots/ZGAMMA_ONLY/plot_Shape.png AllMTPlots/plot_Shape_ZGAMMA_ONLY.png
cp AllMTPlots/WJETS_ONLY/plot_Shape.png AllMTPlots/plot_Shape_WJETS_ONLY.png
cp AllMTPlots/METHOD_ONLY/plot_Shape.png AllMTPlots/plot_Shape_METHOD_ONLY.png
cp AllMTPlots/GAMMASTATS_ONLY/plot_Shape.png AllMTPlots/plot_Shape_GAMMASTATS_ONLY.png

