#!/usr/bin/env bash

#--------------------------------------------------
# Global Code 
#--------------------------------------------------

if [[ $# -eq 0 ]]; then 
    printf "NAME\n\tsubmit.sh - Main driver to submit jobs\n"
    printf "\nSYNOPSIS\n"
    printf "\n\t%-5s\n" "./submit.sh [OPTION]" 
    printf "\nOPTIONS\n" 
## Run 2l2vAnalysis over samples
    printf "\n\t%-5s  %-40s\n"  "0"  "completely clean up the directory" 
    printf "\n\t%-5s  %-40s\n"  "1"  "run 'runHZZ2l2vAnalysis' on all samples" 
    printf "\n\t%-5s  %-40s\n"  "1.0"  "make plots for mcbased analysis"
    printf "\n\t%-5s  %-40s\n"  "1.1"  "run 'runHZZ2l2vAnalysis' on photon samples" 
    printf "\n\t%-5s  %-40s\n"  "1.2"  "run 'runHZZ2l2vAnalysis' on photon samples with photon re-weighting"

## Run photonZ closure tests:
    printf "\n\t%-5s  %-40s\n"  "1.3"  "run 'runHZZ2l2vAnalysis' on photon and DY samples in MC only"
    printf "\n\t%-5s  %-40s\n"  "1.4"  "run 'runHZZ2l2vAnalysis' on photon samples with photon re-weighting with MC weights"  

## Merge Results
    printf "\n\t%-5s  %-40s\n"  "2"  "compute integrated luminosity from processed samples" 
    printf "\n\t%-5s  %-40s\n"  "2.1"  "compute integrated luminosity from processed photon samples" 
    printf "\n\t%-5s  %-40s\n"  "3.0"  "make plots and combine root files" 

## Extract photon weights
    printf "\n\t%-5s  %-40s\n"  "3.01" "extract photon weights using bins in data or MC" 

## Make plots in mcbased(_blind), photons, datadriven(_blind) cases
    printf "\n\t%-5s  %-40s\n"  "3.1"  "make plots for mcbased analysis"  
    printf "\n\t%-5s  %-40s\n"  "3.15"  "make plots for photon_samples" 
    printf "\n\t%-5s  %-40s\n"  "3.16"  "make plots for photonZ analysis MC closure tests"
    printf "\n\t%-5s  %-40s\n"  "3.17"  "make plots for instr. MET and genuine Met"
    printf "\n\t%-5s  %-40s\n"  "3.2"  "make plots with data-driven MET"
fi

step=$1   #variable that store the analysis step to run

#Additional arguments to take into account
arguments=''; for var in "${@:2}"; do arguments=$arguments" "$var; done
if [[ $# -ge 4 ]]; then echo "Additional arguments will be considered: "$arguments ;fi 

#--------------------------------------------------
# Global Variables
#--------------------------------------------------

SUFFIX=_2017_01_25

#SUFFIX=$(date +"_%Y_%m_%d") 
MAINDIR=$CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v
JSON=$MAINDIR/samples_full2016_GGH_VBF.json
GOLDENJSON=$CMSSW_BASE/src/UserCode/llvv_fwk/data/json/
RESULTSDIR=$MAINDIR/results$SUFFIX
PLOTSDIR=$MAINDIR/plots${SUFFIX}
PLOTTER=$MAINDIR/plotter${SUFFIX}

echo "Input: " $JSON
echo "Output: " $RESULTSDIR

queue='8nh'

#IF CRAB3 is provided in argument, use crab submissiong instead of condor/lsf
if [[ $arguments == *"crab3"* ]]; then queue='crab3' ;fi 
 
################################################# STEPS between 0 and 1
if [[ $step == 0 ]]; then   
        #analysis cleanup
	echo "ALL DATA WILL BE LOST! [N/y]?"
	read answer
	if [[ $answer == "y" ]];
	then
	    echo "CLEANING UP..."
	    rm -rdf $RESULTSDIR $PLOTSDIR LSFJOB_* core.* *.sh.e* *.sh.o*
	fi
fi #end of step0

###  ############################################## STEPS between 1 and 2
if [[ $step > 0.999 &&  $step < 2 ]]; then
   if [[ $step == 1 || $step == 1.0 ]]; then        #submit jobs for 2l2v analysis
	echo "JOB SUBMISSION"
	runAnalysisOverSamples.py -e runHZZ2l2vAnalysis -j $JSON -o $RESULTSDIR  -c $MAINDIR/../runAnalysis_cfg.py.templ -p "@data_pileup=datapileup_latest @jacknife=0 @saveSummaryTree=True @runSystematics=False @useMVA=True @jacks=0" -s $queue --report True --key 2l2v_mcbased $arguments
   fi

   if [[ $step == 1 || $step == 1.1 ]]; then        #submit jobs for 2l2v photon jet analysis
	echo "JOB SUBMISSION for Photon + Jet analysis"
	runAnalysisOverSamples.py -e runHZZ2l2vAnalysis -j $JSON -o $RESULTSDIR -c $MAINDIR/../runAnalysis_cfg.py.templ -p "@data_pileup=datapileup_latest @jacknife=0 @saveSummaryTree=True @runSystematics=False @useMVA=True @jacks=0"  -s $queue --report True --key 2l2v_photoncontrol --skipkey 2l2v_mcbased $arguments 
   fi

   if [[ $step == 1 || $step == 1.2 ]]; then        #submit jobs for 2l2v analysis with photon re-weighting
	echo "JOB SUBMISSION for Photon + Jet analysis with photon re-weighting"                                                                        
        runAnalysisOverSamples.py -e runHZZ2l2vAnalysis -j $JSON -o $RESULTSDIR -c $MAINDIR/../runAnalysis_cfg.py.templ -p "@data_pileup=datapileup_latest @useMVA=True @saveSummaryTree=True @weightsFile=$PWD/photonWeights_run2016.root @puWeightsFile=$PWD/get_photon_pileupwgt.root  @runSystematics=False @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0" -s $queue --report True --key 2l2v_photonsOnly $arguments 
   fi    
                     
   if [[ $step == 1.3 ]]; then        #submit jobs for 2l2v analysis with photon re-weighting in MC 
       echo "JOB SUBMISSION for Photon + Jet analysis with photon re-weighting in MC : step1"  
       runAnalysisOverSamples.py -e runHZZ2l2vAnalysis -j $JSON -o $RESULTSDIR -c $MAINDIR/../runAnalysis_cfg.py.templ -p "@data_pileup=datapileup_latest @useMVA=True @saveSummaryTree=True @runSystematics=False @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0" -s $queue --report True --key 2l2v_photonZ $arguments
   fi

   if [[ $step == 1.4 ]]; then        #submit jobs for 2l2v analysis with photon re-weighting in MC                                              
        echo "JOB SUBMISSION for Photon + Jet analysis with photon re-weighting in MC : step2"                                          
        runAnalysisOverSamples.py -e runHZZ2l2vAnalysis -j $JSON -o $RESULTSDIR -c $MAINDIR/../runAnalysis_cfg.py.templ -p "@data_pileup=datapileup_latest @useMVA=True  @saveSummaryTree=True @weightsFile=$PWD/photonWeights_run2016MC.root @runSystematics=False @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0" -s $queue --report True --key 2l2v_mcphotonsOnly $arguments                                                                                                             
   fi    

   if [[ $HOSTNAME =~ "iihe" ]]; then yes | big-submission $RESULTSDIR/FARM/inputs/big.cmd; fi

 
fi

###  ############################################## STEPS between 2 and 3
if [[ $step > 1.999 && $step < 3 ]]; then
   if [[ $step == 2 ]]; then    #extract integrated luminosity of the processed lumi blocks
	echo "MISSING LUMI WILL APPEAR AS DIFFERENCE LUMI ONLY IN in.json"
	mergeJSON.py --output=$RESULTSDIR/json_all.json        $RESULTSDIR/Data*.json
	mergeJSON.py --output=$RESULTSDIR/json_doubleMu.json   $RESULTSDIR/Data*_DoubleMu*.json
	mergeJSON.py --output=$RESULTSDIR/json_doubleEl.json   $RESULTSDIR/Data*_DoubleElectron*.json
	mergeJSON.py --output=$RESULTSDIR/json_muEG.json       $RESULTSDIR/Data*_MuEG*.json
        mergeJSON.py --output=$RESULTSDIR/json_gamma.json      $RESULTSDIR/Data*_SinglePhoton*.json
	mergeJSON.py --output=$RESULTSDIR/json_in.json  $GOLDENJSON/Cert_*.txt
	echo "MISSING LUMI BLOCKS IN DOUBLE MU DATASET"
	compareJSON.py --diff $RESULTSDIR/json_in.json $RESULTSDIR/json_doubleMu.json 
	echo "MISSING LUMI BLOCKS IN DOUBLE ELECTRON DATASET"
	compareJSON.py --diff $RESULTSDIR/json_in.json $RESULTSDIR/json_doubleEl.json 
	echo "MISSING LUMI BLOCKS IN MUON EGAMMA DATASET"
	compareJSON.py --diff $RESULTSDIR/json_in.json $RESULTSDIR/json_muEG.json 
	echo "MISSING LUMI BLOCKS IN Single Photon DATASET"
	compareJSON.py --diff $RESULTSDIR/json_in.json $RESULTSDIR/json_gamma.json 

	echo "COMPUTE INTEGRATED LUMINOSITY"
	export LD_LIBRARY_PATH=/afs/cern.ch/cms/lumi/brilconda-1.1.7/root/lib
        export PYTHONPATH=/afs/cern.ch/cms/lumi/brilconda-1.1.7/root/lib 
        export PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH 
        export ROOTSYS=/afs/cern.ch/cms/lumi/brilconda-1.1.7/root 
        export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.1.7/bin:$PATH 
#	export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.0.3/bin:$PATH

	pip uninstall brilws -y 
	pip install --upgrade --install-option="--prefix=$HOME/.local" brilws &> /dev/null #will be installed only the first time

	if [[ $JSON =~ "2016" ]]; then
	    brilcalc lumi -b "STABLE BEAMS" --normtag /afs/cern.ch/user/l/lumipro/public/Normtags/normtag_DATACERT.json -i $RESULTSDIR/json_all.json -u /pb -o $RESULTSDIR/LUMI.txt 
	else
	    brilcalc lumi --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/moriond16_normtag.json -i $RESULTSDIR/json_all.json -u /pb -o $RESULTSDIR/LUMI.txt 
	fi
	tail -n 3 $RESULTSDIR/LUMI.txt  
     fi
  fi     

###  ############################################## STEPS between 3 and 4
if [[ $step > 2.999 && $step < 4 ]]; then
    if [ -f $RESULTSDIR/LUMI.txt ]; then
      INTLUMI=`tail -n 3 $RESULTSDIR/LUMI.txt | cut -d ',' -f 6`
    else
	if [[ $JSON =~ "full2016" ]]; then  
	    INTLUMI=35866.932
            echo "Please run step==2 above to calculate int. luminosity for 2016 data!" 
        else
	    echo "Please run step==2 above to calculate int. luminosity!"
#            INTLUMI=2268.759 #correspond to the value from DoubleMu OR DoubleEl OR MuEG without jobs failling and golden JSON 
        fi                                                                                                                   
	echo "WARNING: $RESULTSDIR/LUMI.txt file is missing so use fixed integrated luminosity value, this might be different than the dataset you ran on"
    fi
    
    if [[ $step == 3 || $step == 3.0 ]]; then  # make plots and combined root files
	echo "MAKE SUMMARY ROOT FILE, BASED ON AN INTEGRATED LUMINOSITY OF $INTLUMI"
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outFile ${PLOTTER}.root  --json $JSON --noPlot --fileOption RECREATE --key 2l2v_mcbased $arguments        
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outFile ${PLOTTER}.root  --json $JSON --noPlot --fileOption UPDATE   --key 2l2v_mcbased --doInterpollation  $arguments       
        cp  ${PLOTTER}.root  ${PLOTTER}_MCOnly.root

        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outFile ${PLOTTER}.root  --json $JSON --noPlot --fileOption UPDATE   --key 2l2v_photoncontrol               $arguments 
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outFile ${PLOTTER}.root  --json $JSON --noPlot --fileOption UPDATE   --key 2l2v_photonsOnly                 $arguments
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outFile ${PLOTTER}.root  --json $JSON --noPlot --fileOption UPDATE   --key genuineMet                       $arguments 
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outFile ${PLOTTER}.root  --json $JSON --noPlot --fileOption UPDATE   --key 2l2v_datadriven                  $arguments 
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outFile ${PLOTTER}.root  --json $JSON --noPlot --fileOption UPDATE   --key 2l2v_datadrivenplot              $arguments 
	ln -s -f ${PLOTTER}.root plotter.root 
    fi        
   
    if [[ $step == 3 || $step == 3.01 ]]; then  # make plots and combine root files for photon + jet study    
#        extractPhotonWeights --inFile plotter.root --outFile photonWeights_RunDNew.root --outDir $PLOTSDIR/photonWeights --fitf true
#	extractPhotonWeights_UsingBins --inFile ${PLOTTER}.root --outFile photonWeights_run2016MC.root --outDir $PLOTSDIR/photonMCWeights --mode MC --purity noQCD
	extractPhotonWeights_UsingBins --inFile ${PLOTTER}.root --outFile photonWeights_run2016.root --outDir $PLOTSDIR/photonWeights --mode DATA 
    fi


    if [[ $step == 3 || $step == 3.1 ]]; then  # make plots and combine root files for mcbased study    
	runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outFile ${PLOTTER}.root  --json $JSON --noPlot --fileOption RECREATE --key 2l2v_mcbased $arguments 
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outDir $PLOTSDIR/mcbased/ --outFile ${PLOTTER}.root  --json $JSON --plotExt .png --plotExt .pdf  --key 2l2v_mcbased --fileOption READ $arguments
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outDir $PLOTSDIR/mcbased_blind/ --outFile ${PLOTTER}.root  --json $JSON --no2D --plotExt .png --plotExt .pdf  --key 2l2v_mcbased --fileOption READ --blind 125 --only "(all|ll|mumu|ee|emu)(|eq0jets|geq1jets|vbf)_(met|metpuppi)" $arguments
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outDir $PLOTSDIR/mcbased_blind/ --outFile ${PLOTTER}.root  --json $JSON --no2D --plotExt .png --plotExt .pdf  --key 2l2v_mcbased --fileOption READ --blind 400 --only "(all|ll|mumu|ee|emu)(|eq0jets|geq1jets|vbf)_mt" $arguments      
    fi


    if [[ $step == 3 || $step == 3.11 ]]; then  # make plots and combine root files for mcbased study    
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outDir $PLOTSDIR/mcbased/ --outFile ${PLOTTER}.root  --json samples.json --no2D --plotExt .png --plotExt .pdf  --key 2l2v_mcbased --fileOption READ --showUnc 0.0 --only "emu(|eq0jets|geq1jets|vbf)_mt_(In|Out)b.*" --rebin 3  --noLog --removeRatioPlot  $arguments 
    fi


    if [[ $step == 3 || $step == 3.15 ]]; then  # make plots and combine root files for photon + jet study    
	echo "MAKE PLOTS AND SUMMARY ROOT FILE for Photon sample"
	runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outFile ${PLOTTER}.root  --json $JSON --noPlot --fileOption UPDATE   --key 2l2v_photoncontrol $arguments
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outDir $PLOTSDIR/photons/ --outFile ${PLOTTER}.root  --json $JSON --no2D --plotExt .png --plotExt .pdf  --key 2l2v_photoncontrol --fileOption READ --only "(all|ll|mumu|ee|emu|gamma)(_qt|_qmass|_met|_mt|npho.*)"  $arguments 
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outDir $PLOTSDIR/photons/ --outFile ${PLOTTER}.root  --json $JSON --no2D --plotExt .png --plotExt .pdf  --key 2l2v_photoncontrol --fileOption READ --only "gamma.*"  $arguments 
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outDir $PLOTSDIR/photons/ --outFile ${PLOTTER}.root  --json $JSON --no2D --plotExt .png --plotExt .pdf  --key 2l2v_photoncontrol --fileOption READ --only "trg.*"  $arguments 
    fi

    if [[ $step == 3 || $step == 3.16 ]]; then # make photonZ MC closure tests 
	# make sure you have at least DY MC and gamma+jets_reweighted in plotter.root
	echo "MAKE PLOTS for photonZ analysis closure test"
	runPhotonZClosure --inFile ${PLOTTER}.root --outDir $PLOTSDIR/photonZclosure --mode MC --channel ll --purity noQCD
    fi

    if [[ $step == 3 || $step == 3.17 ]]; then  # make plots and combine root files for photon + jet study    
	echo "MAKE PLOTS AND SUMMARY ROOT FILE for GenuineMet in photon sample"
	runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outFile ${PLOTTER}.root  --json $JSON --noPlot --fileOption UPDATE   --key 2l2v_photonsOnly $arguments
	runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outFile ${PLOTTER}.root  --json $JSON --noPlot --fileOption UPDATE   --key genuineMet $arguments
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outDir $PLOTSDIR/photonsGM/ --outFile ${PLOTTER}.root  --json $JSON --no2D --plotExt .png --plotExt .pdf  --key genuineMet --fileOption READ --only "(ll|gamma)(|eq0jets|geq1jets|vbf)(_.*)" $arguments 
    fi

    if [[ $step == 3 || $step == 3.2 ]]; then # make plots and combine root files for photon + jet study   
	runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outFile ${PLOTTER}.root  --json $JSON --noPlot --fileOption UPDATE   --key 2l2v_datadriven $arguments
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outDir $PLOTSDIR/datadriven/       --outFile ${PLOTTER}.root  --json $JSON --no2D --plotExt .png --plotExt .pdf  --key 2l2v_datadriven --fileOption READ $arguments
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outDir $PLOTSDIR/datadriven_blind/ --outFile ${PLOTTER}.root  --json $JSON --no2D --plotExt .png --plotExt .pdf  --key 2l2v_datadriven --fileOption READ --blind 125 --only "(all|ll|mumu|ee|emu)(|eq0jets|geq1jets|vbf)_(met|metpuppi)" $arguments
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outDir $PLOTSDIR/datadriven_blind/ --outFile ${PLOTTER}.root  --json $JSON --no2D --plotExt .png --plotExt .pdf  --key 2l2v_datadriven --fileOption READ --blind 400 --only "(all|ll|mumu|ee|emu)(|eq0jets|geq1jets|vbf)_mt" $arguments
    fi

    if [[ $step == 3 || $step == 3.21 ]]; then  # make plots and combine root files for photon + jet study    
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outDir $PLOTSDIR/datadriven/       --outFile ${PLOTTER}.root  --json $JSON --no2D --plotExt .png --plotExt .pdf  --key 2l2v_datadrivenplot --fileOption READ --showUnc 0.06 --rebin -1 --only "(all|ll|mumu|ee|emu)(|eq0jets|geq1jets|vbf)_(mtSyst|metSyst)" $arguments 
    fi

    if [[ $step == 3.30 ]]; then  # make plots and combined root files
        echo "MAKE SUMMARY ROOT FILE, BASED ON AN INTEGRATED LUMINOSITY OF $INTLUMI FOR LIMIT STUDIES"
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outFile ${PLOTTER}_forLimits.root  --json $JSON --noPlot --fileOption UPDATE --key 2l2v_mcbased --only "(all_optim_systs|all_optim_cut|(mumu|ee|emu)(|eq0jets|geq1jets|vbf)_(mt)(_shapes|_shapes_NRBctrl)(|_cp100.00_brn0.00|_cp10.00_brn0.00|_cp5.00_brn0.00)(|_scale_umetup|_scale_umetdown|_res_jup|_res_jdown|_scale_jup|_scale_jdown|_scale_mup|_scale_mdown|_stat_eup|_stat_edown|_stat_eup|_stat_edown|_sys_eup|_sys_edown|_GS_eup|_GS_edown|_resRho_eup|_resRho_edown|_puup|_pudown|_eff_bup|_eff_bdown|_lepveto|_th_pdf|_th_alphas|_th_ewkup|_th_ewkdown))" $arguments

        cp  ${PLOTTER}_forLimits.root  ${PLOTTER}_forLimits_MCOnly.root

    fi

    if [[ $step == 3.31 ]]; then
        echo ""

        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outFile ${PLOTTER}_forLimits.root  --json $JSON --noPlot --fileOption UPDATE   --key 2l2v_photoncontrol --only "(all_optim_systs|all_optim_cut|(mumu|ee|emu)(|eq0jets|geq1jets|vbf)_(mt)(_shapes|_shapes_NRBctrl)(|_cp100.00_brn0.00|_cp10.00_brn0.00|_cp5.00_brn0.00)(|_scale_umetup|_scale_umetdown|_res_jup|_res_jdown|_scale_jup|_scale_jdown|_scale_mup|_scale_mdown|_stat_eup|_stat_edown|_stat_eup|_stat_edown|_sys_eup|_sys_edown|_GS_eup|_GS_edown|_resRho_eup|_resRho_edown|_puup|_pudown|_eff_bup|_eff_bdown|_lepveto|_th_pdf|_th_alphas|_th_ewkup|_th_ewkdown))" $arguments
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outFile ${PLOTTER}_forLimits.root  --json $JSON --noPlot --fileOption UPDATE   --key 2l2v_photonsOnly  --only "(all_optim_systs|all_optim_cut|(mumu|ee|emu)(|eq0jets|geq1jets|vbf)_(mt)(_shapes|_shapes_NRBctrl)(|_cp100.00_brn0.00|_cp10.00_brn0.00|_cp5.00_brn0.00)(|_scale_umetup|_scale_umetdown|_res_jup|_res_jdown|_scale_jup|_scale_jdown|_scale_mup|_scale_mdown|_stat_eup|_stat_edown|_stat_eup|_stat_edown|_sys_eup|_sys_edown|_GS_eup|_GS_edown|_resRho_eup|_resRho_edown|_puup|_pudown|_eff_bup|_eff_bdown|_lepveto|_th_pdf|_th_alphas|_th_ewkup|_th_ewkdown))" $arguments
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outFile ${PLOTTER}_forLimits.root  --json $JSON --noPlot --fileOption UPDATE   --key genuineMet   --only "(all_optim_systs|all_optim_cut|(mumu|ee|emu)(|eq0jets|geq1jets|vbf)_(mt)(_shapes|_shapes_NRBctrl)(|_cp100.00_brn0.00|_cp10.00_brn0.00|_cp5.00_brn0.00)(|_scale_umetup|_scale_umetdown|_res_jup|_res_jdown|_scale_jup|_scale_jdown|_scale_mup|_scale_mdown|_stat_eup|_stat_edown|_stat_eup|_stat_edown|_sys_eup|_sys_edown|_GS_eup|_GS_edown|_resRho_eup|_resRho_edown|_puup|_pudown|_eff_bup|_eff_bdown|_lepveto|_th_pdf|_th_alphas|_th_ewkup|_th_ewkdown))" $arguments

        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outFile ${PLOTTER}_forLimits.root  --json $JSON --noPlot --fileOption UPDATE   --key 2l2v_datadriven  --only "(all_optim_systs|all_optim_cut|(mumu|ee|emu)(|eq0jets|geq1jets|vbf)_(mt)(_shapes|_shapes_NRBctrl)(|_cp100.00_brn0.00|_cp10.00_brn0.00|_cp5.00_brn0.00)(|_scale_umetup|_scale_umetdown|_res_jup|_res_jdown|_scale_jup|_scale_jdown|_scale_mup|_scale_mdown|_stat_eup|_stat_edown|_stat_eup|_stat_edown|_sys_eup|_sys_edown|_GS_eup|_GS_edown|_resRho_eup|_resRho_edown|_puup|_pudown|_eff_bup|_eff_bdown|_lepveto|_th_pdf|_th_alphas|_th_ewkup|_th_ewkdown))" $arguments
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outFile ${PLOTTER}_forLimits.root  --json $JSON --noPlot --fileOption UPDATE   --key 2l2v_datadrivenplot  --only "(all_optim_systs|all_optim_cut|(mumu|ee|emu)(|eq0jets|geq1jets|vbf)_(mt)(_shapes|_shapes_NRBctrl)(|_cp100.00_brn0.00|_cp10.00_brn0.00|_cp5.00_brn0.00)(|_scale_umetup|_scale_umetdown|_res_jup|_res_jdown|_scale_jup|_scale_jdown|_scale_mup|_scale_mdown|_stat_eup|_stat_edown|_stat_eup|_stat_edown|_sys_eup|_sys_edown|_GS_eup|_GS_edown|_resRho_eup|_resRho_edown|_puup|_pudown|_eff_bup|_eff_bdown|_lepveto|_th_pdf|_th_alphas|_th_ewkup|_th_ewkdown))" $arguments
    fi

    if [[ $step == 3.32 ]]; then  # make plots and combined root files
        echo "MAKE SUMMARY ROOT FILE, BASED ON AN INTEGRATED LUMINOSITY OF $INTLUMI TO COMPUTE MELA XSEC"
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outFile ${PLOTTER}_forMELAXSec.root  --json $JSON --noPlot --fileOption UPDATE --key 2l2v_mcbased --only "(all_higgsMass_shape)(|_cp100.00_brn0.00|_cp10.00_brn0.00|_cp5.00_brn0.00)" $arguments

    fi


fi

