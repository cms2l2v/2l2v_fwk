cp MakeSyst_custom_forMTandMET.C test.C 
rm -f InstrMET_systematics.root
#all except gamma stats and genuine stats
sed 's/ALL true/ALL false/' test.C | sed 's/GAMMASTATS_ONLY true/GAMMASTATS_ONLY false/' | sed 's/GENUINESTATS_ONLY true/GENUINESTATS_ONLY false/' | sed 's/DO_CORRECTIONS true/DO_CORRECTIONS false/' | sed 's/void MakeSyst_custom_forMTandMET()/void tmpMakeSyst()/' > tmpMakeSyst.C
root -l -q -b tmpMakeSyst.C
mv InstrMET_systematics.root $CMSSW_BASE/src/UserCode/llvv_fwk/data/InstrMET_systematics/InstrMET_systematics_ALL_EXCEPT_GAMMASTATS.root 
rm -f tmpMakeSyst.C
#all except gamma stats and ZnunuG
#sed 's/ALL true/ALL false/' test.C | sed 's/GAMMASTATS_ONLY true/GAMMASTATS_ONLY false/' | sed 's/ZGAMMA_ONLY true/ZGAMMA_ONLY false/' | sed 's/void MakeSyst_custom_forMTandMET()/void tmpMakeSyst()/' > tmpMakeSyst.C
#root -l -q -b tmpMakeSyst.C
#mv InstrMET_systematics.root $CMSSW_BASE/src/UserCode/llvv_fwk/data/InstrMET_systematics/InstrMET_systematics_ALL_EXCEPT_GAMMASTATS.root 
#rm -f tmpMakeSyst.C
#Gamma stats only
sed 's/ALL true/ALL false/' test.C | sed 's/WGAMMA_ONLY true/WGAMMA_ONLY false/' | sed 's/ZGAMMA_ONLY true/ZGAMMA_ONLY false/' | sed 's/WJETS_ONLY true/WJETS_ONLY false/' | sed 's/METHOD_ONLY true/METHOD_ONLY false/' | sed 's/DO_CORRECTIONS true/DO_CORRECTIONS false/' | sed 's/void MakeSyst_custom_forMTandMET()/void tmpMakeSyst()/' > tmpMakeSyst.C
root -l -q -b tmpMakeSyst.C
mv InstrMET_systematics.root $CMSSW_BASE/src/UserCode/llvv_fwk/data/InstrMET_systematics/InstrMET_systematics_GAMMASTATS.root
rm -f tmpMakeSyst.C
rm -f test.C
