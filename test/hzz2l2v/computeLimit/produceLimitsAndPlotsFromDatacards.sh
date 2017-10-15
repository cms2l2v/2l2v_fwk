#Before launching this script you need to have launch optimize_WideWidth.py -p 4 to have produced all datacards.
#Then you can modify those cards by hand.... or change the computeLimit code would be even more straitghforward.
#NB: no need to change the combined one, this will be done by this script

if [[ $# -ne 2 ]]; then
	echo "Wrong number of parameters! This command asks for mass (0100, 0800, 1500...) and width (5.0, 10.0 or 100.0)"
  exit 0;
fi


massWithZero=$1
width=$2
mass=$((10#$massWithZero)) #small trick to remove leading zeros if mass is given in the format 0800


echo "Producing mT plots and limits for mass of $mass and width of $width GeV"
echo "Warning: width should be:  5.0, 10.0 or 100.0"

cd /storage_mnt/storage/user/delannoy/HZZ2l2v/2_October_LimitsImprovement/CMSSW_7_4_7/src/UserCode/llvv_fwk/test/hzz2l2v/computeLimit/cards_SB13TeV_SM_cp${width}0_brn0.00/$massWithZero;
sh combineCards.sh;
text2workspace.py card_combined.dat -o workspace.root -P UserCode.llvv_fwk.HeavyScalarMod:heavyscalarmod --PO verbose --PO 'is2l2nu' --PO m="${mass}" --PO w="${width}"
combine -M ProfileLikelihood --signif --pvalue -m $mass  workspace.root --setPhysicsModelParameters fvbf=0 --freezeNuisances fvbf > COMB.log;
mv higgsCombineTest.ProfileLikelihood.mH${mass}.root higgsCombineTest.ProfileLikelihood.mH${mass}_ggH.root;
combine -M ProfileLikelihood --signif --pvalue -m $mass  workspace.root --setPhysicsModelParameters fvbf=1 --freezeNuisances fvbf > COMB_VBF.log;
mv higgsCombineTest.ProfileLikelihood.mH${mass}.root higgsCombineTest.ProfileLikelihood.mH${mass}_qqH.root;
combine -M ProfileLikelihood --signif --pvalue -m $mass  workspace.root --setPhysicsModelParameters fvbf=1.0 > COMB_float.log;
mv higgsCombineTest.ProfileLikelihood.mH${mass}.root higgsCombineTest.ProfileLikelihood.mH${mass}_ppH.root;
combine -M MaxLikelihoodFit --saveNormalizations --saveShapes --saveWithUncertainties -m $mass workspace.root  --setPhysicsModelParameters fvbf=0.0 --freezeNuisances fvbf
mv mlfit.root mlfit_ggH.root;
python /storage_mnt/storage/user/delannoy/HZZ2l2v/2_October_LimitsImprovement/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py  mlfit_ggH.root -g Nuisance_CrossCheck_ggH.root
sh /storage_mnt/storage/user/delannoy/HZZ2l2v/2_October_LimitsImprovement/CMSSW_7_4_7/src/UserCode/llvv_fwk/test/hzz2l2v/computeLimit/produceThePerCategoryPlot.sh $mass ggH
combine -M MaxLikelihoodFit --saveNormalizations --saveShapes --saveWithUncertainties -m $mass workspace.root  --setPhysicsModelParameters fvbf=1.0 --freezeNuisances fvbf
mv mlfit.root mlfit_qqH.root;
python /storage_mnt/storage/user/delannoy/HZZ2l2v/2_October_LimitsImprovement/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py  mlfit_qqH.root -g Nuisance_CrossCheck_qqH.root
sh /storage_mnt/storage/user/delannoy/HZZ2l2v/2_October_LimitsImprovement/CMSSW_7_4_7/src/UserCode/llvv_fwk/test/hzz2l2v/computeLimit/produceThePerCategoryPlot.sh $mass qqH
combine -M MaxLikelihoodFit --saveNormalizations --saveShapes --saveWithUncertainties -m $mass workspace.root  --setPhysicsModelParameters fvbf=1.0
mv mlfit.root mlfit_ppH.root;
python /storage_mnt/storage/user/delannoy/HZZ2l2v/2_October_LimitsImprovement/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py  mlfit_ppH.root -g Nuisance_CrossCheck_ppH.root
sh /storage_mnt/storage/user/delannoy/HZZ2l2v/2_October_LimitsImprovement/CMSSW_7_4_7/src/UserCode/llvv_fwk/test/hzz2l2v/computeLimit/produceThePerCategoryPlot.sh $mass ppH
combine -M Asymptotic --picky -m $mass workspace.root --setPhysicsModelParameters fvbf=0 --freezeNuisances fvbf  >>  COMB.log;
mv higgsCombineTest.Asymptotic.mH${mass}.root higgsCombineTest.Asymptotic.mH${mass}_ggH.root;
combine -M Asymptotic --picky -m $mass workspace.root --setPhysicsModelParameters fvbf=1 --freezeNuisances fvbf  >>  COMB_VBF.log;
mv higgsCombineTest.Asymptotic.mH${mass}.root higgsCombineTest.Asymptotic.mH${mass}_qqH.root;
combine -M Asymptotic --picky -m $mass workspace.root --setPhysicsModelParameters fvbf=1  >>  COMB_float.log;
mv higgsCombineTest.Asymptotic.mH${mass}.root higgsCombineTest.Asymptotic.mH${mass}_ppH.root;
combine -M Asymptotic --picky -m $mass workspace.root --setPhysicsModelParameters fvbf=0 --freezeNuisances fvbf --run blind >>  COMB_blind.log;
mv higgsCombineTest.Asymptotic.mH${mass}.root higgsCombineTest.Asymptotic.mH${mass}_ggH_blinded.root;
combine -M Asymptotic --picky -m $mass workspace.root --setPhysicsModelParameters fvbf=1 --freezeNuisances fvbf --run blind  >>  COMB_VBF_blind.log;
mv higgsCombineTest.Asymptotic.mH${mass}.root higgsCombineTest.Asymptotic.mH${mass}_qqH_blinded.root;
combine -M Asymptotic --picky -m $mass workspace.root --setPhysicsModelParameters fvbf=1 --run blind  >>  COMB_float_blind.log;
mv higgsCombineTest.Asymptotic.mH${mass}.root higgsCombineTest.Asymptotic.mH${mass}_ppH_blinded.root;
cd -;
echo "Done!"
exit 1;
