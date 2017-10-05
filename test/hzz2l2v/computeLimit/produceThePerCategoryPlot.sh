#!/usr/bin/env bash

#help
if [[ "$#" == "0" ]]; then
  echo "This command produces the final MT plots per category with uncertainties and for all cases (datacard_input, prefit, postfitB and postfitSandB).";
    echo "Usage: 'produceThePerCategoryPlot.sh [mass] [ggH|qqH|floating]'";
    echo "#Mass  = the mass";
    echo "#Mode  = the production mode (ggH only, qqH only, or a floating ratio between the two)           ";
    echo "";
    echo "Example: 'sh produceThePerCategoryPlot.sh 0800 ggH'";
    exit 0;
fi

#check CMSSW
if [ "$CMSSW_BASE" == "" ]; then
    echo "Please set up a cmsenv"
    exit 0
fi

#store arguments
mass_old=$1  #the mass, has to be given with 4 numbers (0100, 0200, ... 1000, 1500...)
mode=$2  #the production mode (ggH only, qqH only or a floating ratio between the two)

mass=$((10#$mass_old)) #small trick to remove leading zeros if mass is given in the format 0800


echo "Producing MT plots for a mass of ${mass} for the ${mode} production mode"
mkdir -p mt_plots

#The code will produce the 6 MT plots for 4 cases:
##1. Normalisation and Error directly from inputs to datacards
sed s/TEMPLATE_MASS/${mass}/g ${CMSSW_BASE}/src/UserCode/llvv_fwk/test/hzz2l2v/computeLimit/drawThePerCategoryPlot_TEMPLATE.C | \
sed s/TEMPLATE_MODE/${mode}/g | \
sed s/TEMPLATE_USE_MLFIT/false/g | \
sed s/TEMPLATE_MACRO_NAME/drawThePerCategoryPlot_M${mass}_${mode}_inputsDatacard/g \
> drawThePerCategoryPlot_M${mass}_${mode}_inputsDatacard.C

root -l -q -b drawThePerCategoryPlot_M${mass}_${mode}_inputsDatacard.C 
mkdir -p mt_plots/M${mass}_${mode}_inputsDatacard
mv mtplots_* mt_plots/M${mass}_${mode}_inputsDatacard/.

##2. Normalisation and Error directly from prefit
sed s/TEMPLATE_MASS/${mass}/g ${CMSSW_BASE}/src/UserCode/llvv_fwk/test/hzz2l2v/computeLimit/drawThePerCategoryPlot_TEMPLATE.C | \
sed s/TEMPLATE_MODE/${mode}/g | \
sed s/TEMPLATE_USE_MLFIT/true/g | \
sed s/TEMPLATE_ERROR/shapes_prefit/g | \
sed s/TEMPLATE_MACRO_NAME/drawThePerCategoryPlot_M${mass}_${mode}_prefit/g \
> drawThePerCategoryPlot_M${mass}_${mode}_prefit.C

root -l -q -b drawThePerCategoryPlot_M${mass}_${mode}_prefit.C 
mkdir -p mt_plots/M${mass}_${mode}_prefit
mv mtplots_* mt_plots/M${mass}_${mode}_prefit/.


##3. Normalisation and Error directly from postfit with background only
sed s/TEMPLATE_MASS/${mass}/g ${CMSSW_BASE}/src/UserCode/llvv_fwk/test/hzz2l2v/computeLimit/drawThePerCategoryPlot_TEMPLATE.C | \
sed s/TEMPLATE_MODE/${mode}/g | \
sed s/TEMPLATE_USE_MLFIT/true/g | \
sed s/TEMPLATE_ERROR/shapes_fit_b/g | \
sed s/TEMPLATE_MACRO_NAME/drawThePerCategoryPlot_M${mass}_${mode}_postfitB/g \
> drawThePerCategoryPlot_M${mass}_${mode}_postfitB.C

root -l -q -b drawThePerCategoryPlot_M${mass}_${mode}_postfitB.C 
mkdir -p mt_plots/M${mass}_${mode}_postfitB
mv mtplots_* mt_plots/M${mass}_${mode}_postfitB/.


##4. Normalisation and Error directly from postfit with backgrond+signal
sed s/TEMPLATE_MASS/${mass}/g ${CMSSW_BASE}/src/UserCode/llvv_fwk/test/hzz2l2v/computeLimit/drawThePerCategoryPlot_TEMPLATE.C | \
sed s/TEMPLATE_MODE/${mode}/g | \
sed s/TEMPLATE_USE_MLFIT/true/g | \
sed s/TEMPLATE_ERROR/shapes_fit_s/g | \
sed s/TEMPLATE_MACRO_NAME/drawThePerCategoryPlot_M${mass}_${mode}_postfitSandB/g \
> drawThePerCategoryPlot_M${mass}_${mode}_postfitSandB.C

root -l -q -b drawThePerCategoryPlot_M${mass}_${mode}_postfitSandB.C 
mkdir -p mt_plots/M${mass}_${mode}_postfitSandB
mv mtplots_* mt_plots/M${mass}_${mode}_postfitSandB/.

#Cleaning
rm -f drawThePerCategoryPlot_M${mass}_${mode}_inputsDatacard.C
rm -f drawThePerCategoryPlot_M${mass}_${mode}_prefit.C
rm -f drawThePerCategoryPlot_M${mass}_${mode}_postfitB.C
rm -f drawThePerCategoryPlot_M${mass}_${mode}_postfitSandB.C
