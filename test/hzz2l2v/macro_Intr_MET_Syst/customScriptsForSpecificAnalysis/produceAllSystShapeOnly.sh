#sed 's/ALL true/ALL false/' test.C | sed 's/WGAMMA_ONLY true/WGAMMA_ONLY false/' | sed 's/ZGAMMA_ONLY true/ZGAMMA_ONLY false/' | sed 's/WJETS_ONLY true/WJETS_ONLY false/' | sed 's/METHOD_ONLY true/METHOD_ONLY false/' | sed 's/GAMMASTATS_ONLY true/GAMMASTATS_ONLY false/' | sed 's/void MakeSyst_custom_forMTandMET()/void tmpMakeSyst()/' > tmpMakeSyst.C
mkdir -p SystBySyst
rm -f InstrMET_systematics.root
#all
sed 's/WGAMMA_ONLY true/WGAMMA_ONLY false/' test.C | sed 's/ZGAMMA_ONLY true/ZGAMMA_ONLY false/' | sed 's/WJETS_ONLY true/WJETS_ONLY false/' | sed 's/METHOD_ONLY true/METHOD_ONLY false/' | sed 's/GAMMASTATS_ONLY true/GAMMASTATS_ONLY false/' | sed 's/void MakeSyst_custom_forMTandMET()/void tmpMakeSyst()/' > tmpMakeSyst.C
root -l -q -b tmpMakeSyst.C
mkdir -p SystBySyst/ALL
mv InstrMET_systematics.root SystBySyst/ALL/.
rm -f tmpMakeSyst.C
#Wgamma only
sed 's/ALL true/ALL false/' test.C | sed 's/ZGAMMA_ONLY true/ZGAMMA_ONLY false/' | sed 's/WJETS_ONLY true/WJETS_ONLY false/' | sed 's/METHOD_ONLY true/METHOD_ONLY false/' | sed 's/GAMMASTATS_ONLY true/GAMMASTATS_ONLY false/' | sed 's/void MakeSyst_custom_forMTandMET()/void tmpMakeSyst()/' > tmpMakeSyst.C
root -l -q -b tmpMakeSyst.C
mkdir -p SystBySyst/WGAMMA_ONLY
mv InstrMET_systematics.root SystBySyst/WGAMMA_ONLY/.
rm -f tmpMakeSyst.C
#Zgamma only
sed 's/ALL true/ALL false/' test.C | sed 's/WGAMMA_ONLY true/WGAMMA_ONLY false/' | sed 's/WJETS_ONLY true/WJETS_ONLY false/' | sed 's/METHOD_ONLY true/METHOD_ONLY false/' | sed 's/GAMMASTATS_ONLY true/GAMMASTATS_ONLY false/' | sed 's/void MakeSyst_custom_forMTandMET()/void tmpMakeSyst()/' > tmpMakeSyst.C
root -l -q -b tmpMakeSyst.C
mkdir -p SystBySyst/ZGAMMA_ONLY
mv InstrMET_systematics.root SystBySyst/ZGAMMA_ONLY/.
rm -f tmpMakeSyst.C
#Wjets only
sed 's/ALL true/ALL false/' test.C | sed 's/WGAMMA_ONLY true/WGAMMA_ONLY false/' | sed 's/ZGAMMA_ONLY true/ZGAMMA_ONLY false/' | sed 's/METHOD_ONLY true/METHOD_ONLY false/' | sed 's/GAMMASTATS_ONLY true/GAMMASTATS_ONLY false/' | sed 's/void MakeSyst_custom_forMTandMET()/void tmpMakeSyst()/' > tmpMakeSyst.C
root -l -q -b tmpMakeSyst.C
mkdir -p SystBySyst/WJETS_ONLY
mv InstrMET_systematics.root SystBySyst/WJETS_ONLY/.
rm -f tmpMakeSyst.C
#Method only
sed 's/ALL true/ALL false/' test.C | sed 's/WGAMMA_ONLY true/WGAMMA_ONLY false/' | sed 's/ZGAMMA_ONLY true/ZGAMMA_ONLY false/' | sed 's/WJETS_ONLY true/WJETS_ONLY false/' | sed 's/GAMMASTATS_ONLY true/GAMMASTATS_ONLY false/' | sed 's/void MakeSyst_custom_forMTandMET()/void tmpMakeSyst()/' > tmpMakeSyst.C
root -l -q -b tmpMakeSyst.C
mkdir -p SystBySyst/METHOD_ONLY
mv InstrMET_systematics.root SystBySyst/METHOD_ONLY/.
rm -f tmpMakeSyst.C
#Gamma stats only
sed 's/ALL true/ALL false/' test.C | sed 's/WGAMMA_ONLY true/WGAMMA_ONLY false/' | sed 's/ZGAMMA_ONLY true/ZGAMMA_ONLY false/' | sed 's/WJETS_ONLY true/WJETS_ONLY false/' | sed 's/METHOD_ONLY true/METHOD_ONLY false/' | sed 's/void MakeSyst_custom_forMTandMET()/void tmpMakeSyst()/' > tmpMakeSyst.C
root -l -q -b tmpMakeSyst.C
mkdir -p SystBySyst/GAMMASTATS_ONLY
mv InstrMET_systematics.root SystBySyst/GAMMASTATS_ONLY/.
rm -f tmpMakeSyst.C

