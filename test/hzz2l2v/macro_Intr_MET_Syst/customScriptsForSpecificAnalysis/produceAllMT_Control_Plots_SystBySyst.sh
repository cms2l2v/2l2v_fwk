#ALL
cp SystBySyst/ALL/InstrMET_systematics.root $CMSSW_BASE/src/UserCode/llvv_fwk/data/InstrMET_systematics/.
cd ..
sed 's/$PLOTSDIR\/datadriven\//$PLOTSDIR\/datadriven\/ALL\//' submit.sh | sed 's/$PLOTSDIR\/datadriven_blind\//$PLOTSDIR\/datadriven_blind\/ALL\//' > tmp_submit.sh
sh tmp_submit.sh 3.222
cd -
#WGamma
cp SystBySyst/WGAMMA_ONLY/InstrMET_systematics.root $CMSSW_BASE/src/UserCode/llvv_fwk/data/InstrMET_systematics/.
cd ..
sed 's/$PLOTSDIR\/datadriven\//$PLOTSDIR\/datadriven\/WGAMMA_ONLY\//' submit.sh | sed 's/$PLOTSDIR\/datadriven_blind\//$PLOTSDIR\/datadriven_blind\/WGAMMA_ONLY\//' > tmp_submit.sh
sh tmp_submit.sh 3.222
cd -
#ZGamma
cp SystBySyst/ZGAMMA_ONLY/InstrMET_systematics.root $CMSSW_BASE/src/UserCode/llvv_fwk/data/InstrMET_systematics/.
cd ..
sed 's/$PLOTSDIR\/datadriven\//$PLOTSDIR\/datadriven\/ZGAMMA_ONLY\//' submit.sh | sed 's/$PLOTSDIR\/datadriven_blind\//$PLOTSDIR\/datadriven_blind\/ZGAMMA_ONLY\//' > tmp_submit.sh
sh tmp_submit.sh 3.222
cd -
#Wjets
cp SystBySyst/WJETS_ONLY/InstrMET_systematics.root $CMSSW_BASE/src/UserCode/llvv_fwk/data/InstrMET_systematics/.
cd ..
sed 's/$PLOTSDIR\/datadriven\//$PLOTSDIR\/datadriven\/WJETS_ONLY\//' submit.sh | sed 's/$PLOTSDIR\/datadriven_blind\//$PLOTSDIR\/datadriven_blind\/WJETS_ONLY\//' > tmp_submit.sh
sh tmp_submit.sh 3.222
cd -
#Method
cp SystBySyst/METHOD_ONLY/InstrMET_systematics.root $CMSSW_BASE/src/UserCode/llvv_fwk/data/InstrMET_systematics/.
cd ..
sed 's/$PLOTSDIR\/datadriven\//$PLOTSDIR\/datadriven\/METHOD_ONLY\//' submit.sh | sed 's/$PLOTSDIR\/datadriven_blind\//$PLOTSDIR\/datadriven_blind\/METHOD_ONLY\//' > tmp_submit.sh
sh tmp_submit.sh 3.222
cd -
#Gamma stats
cp SystBySyst/GAMMASTATS_ONLY/InstrMET_systematics.root $CMSSW_BASE/src/UserCode/llvv_fwk/data/InstrMET_systematics/.
cd ..
sed 's/$PLOTSDIR\/datadriven\//$PLOTSDIR\/datadriven\/GAMMASTATS_ONLY\//' submit.sh | sed 's/$PLOTSDIR\/datadriven_blind\//$PLOTSDIR\/datadriven_blind\/GAMMASTATS_ONLY\//' > tmp_submit.sh
sh tmp_submit.sh 3.222
cd -


#Remove tmp files
cd ..
rm -f tmp_submit.sh
cd -

