#!/bin/bash

if [ "${1}" = "fwlite" ]; then
    BASEDIR=/afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_fwlite/
    if [ "${2}" = "anal_sus" ]; then
	echo "Not ready yet" 
    elif [ "${2}" = "anal_sm" ]; then
	runLocalAnalysisOverSamples.py -e runChHiggsAnalysisFWLite -j $CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/all-samples_fwlite.json -o ${BASEDIR} -d  /store/group/phys_higgs/cmshzz2l2v/2013_08_30/ -c $CMSSW_BASE/src/UserCode/llvv_fwk/test/runAnalysis_cfg.py.templ -p "@useMVA=True @saveSummaryTree=True @runSystematics=True @automaticSwitch=False @is2011=False @jacknife=0 @jacks=0 @weightsFile='${CMSSW_BASE}/src/UserCode/llvv_fwk/data/weights/'" -s 8nh
    elif [ "${2}" = "plots" ]; then
	runPlotterFWLite --iEcm 8 --iLumi 19782 --inDir ${BASEDIR} --outDir ${BASEDIR}plots/ --outFile ${BASEDIR}plotter.root  --json $CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/plot-samples_fwlite.json --noPowers --showUnc
   elif [ "${2}" = "display" ]; then	
	cp ${BASEDIR}/plots/ee_*       ~/www/newAnal/ee/
	cp ${BASEDIR}/plots/emu_*      ~/www/newAnal/emu/
	cp ${BASEDIR}/plots/mumu_*     ~/www/newAnal/mumu/
	cp ${BASEDIR}/plots/singlemu_* ~/www/newAnal/singlemu/
   elif [ "${2}" = "cleanDisplay" ]; then	
	rm ~/www/newAnal/ee/*png
	rm ~/www/newAnal/emu/*png
	rm ~/www/newAnal/mumu/*png
	rm ~/www/newAnal/singlemu/*png

	rm ~/www/newAnal/ee/*tex
	rm ~/www/newAnal/emu/*tex
	rm ~/www/newAnal/mumu/*tex
	rm ~/www/newAnal/singlemu/*tex
   elif [ "${2}" = "cleanArea" ]; then	
	rm ${BASEDIR}/plots/*
	rm -r ${BASEDIR}/FARM
	rm -r ${BASEDIR}/*py
	rm -r ${BASEDIR}/*txt
	rm -r ${BASEDIR}/*root

    fi
elif [ "${1}" = "current" ]; then
# Fixed run 
    #    BASEDIR=/afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_e/
#    BASEDIR=/afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_dio/
    BASEDIR=/afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311_5315/
# run in 5315
    if [ "${2}" = "anal_sus" ]; then
	runLocalAnalysisOverSamples.py -e runChHiggsAnalysis -j data/chhiggs/ch-higgs_samples.json -d /afs/cern.ch/work/v/vischia/private/store/5311_ntuples/ -o ${BASEDIR} -c test/runAnalysis_cfg.py.templ -p "@runSystematics=True @saveSummaryTree=False @weightsFile='${CMSSW_BASE}/src/UserCode/llvv_fwk/data/weights/'" -s 8nh
    elif [ "${2}" = "anal_sm" ]; then
       runLocalAnalysisOverSamples.py -e runChHiggsAnalysis -j data/top_samples_pre.json -d /store/cmst3/user/psilva/5311_ntuples/             -o ${BASEDIR}     -c test/runAnalysis_cfg.py.templ -p "@runSystematics=True @saveSummaryTree=False @weightsFile='${CMSSW_BASE}/src/UserCode/llvv_fwk/data/weights/'" -s 8nh

    elif [ "${2}" = "plots" ]; then
	# Plots
	JSONFILE=data/chhiggs/plot-ch-higgs_tanb30_samples.json
	#JSONFILE=data/chhiggs/plot-ch-higgs_samples.json
	for plotList in evtflow geq2btagsmet geq2btagsnbjets geq2btagsptlep geq2btagssumpt met mll mtsum nbjets njets nvertices nverticesUnweighted pte ptjet1eta ptjet1pt ptjet2eta ptjet2pt ptmin ptmu sumpt
	  do
	  for chanList in emu ee mumu
	    do
	    for formatList in pdf png C
	      do
	      runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json ${JSONFILE} --outFile ${BASEDIR}plotter-forPlotting_${chanList}_${plotList}_${formatList}.root             --showUnc --plotExt .${formatList} --noPowers --onlyStartWith ${chanList}_${plotList}                                 &
	    done
	  done
	done

    elif [ "${2}" = "tables" ]; then
	# Tables
	runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/chhiggs/all-samples.json --outFile ${BASEDIR}plotter-forTables.root --showUnc --noPlots --noPowers --onlyStartWith emu_evtflow
	mkdir -p ${BASEDIR}tables
	mv ${BASEDIR}plotsemu* ${BASEDIR}tables/
	mv ${BASEDIR}plotsee* ${BASEDIR}tables/
	mv ${BASEDIR}plotsmumu* ${BASEDIR}tables/

    elif [ "${2}" = "datacards" ]; then
	
	for chanList in emu ee mumu
	  do
	  runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}datacards/ --json data/chhiggs/all-samples.json     --outFile ${BASEDIR}plotter_${chanList}_forSystTableInPAS.root --showUnc --noPlots --noPowers --onlyStartWith ${chanList}_evtflow
	  runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}datacards/ --json data/chhiggs/all-samples.json     --outFile ${BASEDIR}plotter_${chanList}_forSystTable.root      --showUnc --noPlots --noPowers --onlyStartWith ${chanList}_finalevtflow2btags
	  runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/chhiggs/all-samples_higgs1pb.json --outFile ${BASEDIR}plotter_${chanList}_all-samplesForDatacards_finalevtflow_norm.root --noPlot --noPowers  --onlyStartWith ${chanList}_finalevtflow2btags
	  runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}datacards/ --json data/chhiggs/all-samples.json     --outFile ${BASEDIR}plotter_${chanList}_forSystTableInPAS_optim.root --showUnc --noPlots --noPowers --onlyStartWith all_optim
	  runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}datacards/ --json data/chhiggs/all-samples.json     --outFile ${BASEDIR}plotter_${chanList}_forSystTable_optim.root --showUnc --noPlots --noPowers --onlyStartWith all_optim
	  runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}plots --json data/chhiggs/all-samples_higgs1pb.json --outFile ${BASEDIR}plotter_${chanList}_all-samplesForDatacards_finalevtflow_norm_optim.root --noPlot --noPowers  --onlyStartWith all_optim
	  
	  # Merge syst components list 
	  hadd -f ${BASEDIR}plotter_${chanList}_forSystTableInPAS_def.root		             ${BASEDIR}plotter_${chanList}_forSystTableInPAS.root			 	${BASEDIR}plotter_${chanList}_forSystTableInPAS_optim.root
	  hadd -f ${BASEDIR}plotter_${chanList}_forSystTable_def.root                              ${BASEDIR}plotter_${chanList}_forSystTable.root                             	${BASEDIR}plotter_${chanList}_forSystTable_optim.root                             
	  hadd -f ${BASEDIR}plotter_${chanList}_all-samplesForDatacards_finalevtflow_norm_def.root ${BASEDIR}plotter_${chanList}_all-samplesForDatacards_finalevtflow_norm.root	${BASEDIR}plotter_${chanList}_all-samplesForDatacards_finalevtflow_norm_optim.root
	  
	done
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

	for i in 180 200 220 250 300 350 500 600 700
	  do
	  for chanList in emu ee mumu  
	    do
	    prepareChHiggsDatacards --in ${BASEDIR}plotter_${chanList}_all-samplesForDatacards_finalevtflow_norm_def.root --out datacardsByDecayNew/${i}/ --suffix tb --json /afs/cern.ch/work/v/vischia/private/results/HIG-13-026/tempjsonByFinalState/${i}_tb.json --noPowers --histo finalevtflow2btags --bin 1 --ch ${chanList} & 
	    prepareChHiggsDatacards --in ${BASEDIR}plotter_${chanList}_all-samplesForDatacards_finalevtflow_norm_def.root --out datacardsByDecayNew/${i}/ --suffix taunu --json /afs/cern.ch/work/v/vischia/private/results/HIG-13-026/tempjsonByFinalState/${i}_taunu.json --noPowers --histo finalevtflow2btags --bin 1 --ch ${chanList} & 
	    
	    prepareChHiggsDatacards --in ${BASEDIR}plotter_${chanList}_forSystTable_def.root --out datacardsByDecaySyst/${i}/ --suffix tb --json /afs/cern.ch/work/v/vischia/private/results/HIG-13-026/tempjsonByFinalState/${i}_tb.json --noPowers --histo finalevtflow2btags --bin 1 --ch ${chanList} & 
	    prepareChHiggsDatacards --in ${BASEDIR}plotter_${chanList}_forSystTable_def.root --out datacardsByDecaySyst/${i}/ --suffix taunu --json /afs/cern.ch/work/v/vischia/private/results/HIG-13-026/tempjsonByFinalState/${i}_taunu.json --noPowers --histo finalevtflow2btags --bin 1 --ch ${chanList} & 
	    
	    prepareChHiggsDatacards --in ${BASEDIR}plotter_${chanList}_forSystTableInPAS_def.root --out datacardsByDecaySystPAS/${i}/ --suffix tb --json /afs/cern.ch/work/v/vischia/private/results/HIG-13-026/tempjsonByFinalState/${i}_tb.json --noPowers --histo evtflow --bin 1 --ch ${chanList} & 
	    prepareChHiggsDatacards --in ${BASEDIR}plotter_${chanList}_forSystTableInPAS_def.root --out datacardsByDecaySystPAS/${i}/ --suffix taunu --json /afs/cern.ch/work/v/vischia/private/results/HIG-13-026/tempjsonByFinalState/${i}_taunu.json --noPowers --histo evtflow --bin 1 --ch ${chanList} & 
	  done
	done
	
    elif [ "${2}" = "mhmax" ]; then
	
	mkdir -p datacardsByDecayMhmax/180
	mkdir -p datacardsByDecayMhmax/200 
	mkdir -p datacardsByDecayMhmax/220 
	mkdir -p datacardsByDecayMhmax/250 
	mkdir -p datacardsByDecayMhmax/300

	mkdir -p datacardsByDecayScan/180
	mkdir -p datacardsByDecayScan/200 
	mkdir -p datacardsByDecayScan/220 
	mkdir -p datacardsByDecayScan/250 
	mkdir -p datacardsByDecayScan/300
	
	#JSONDIR=/afs/cern.ch/work/v/vischia/private/results/HIG-13-026/tempjsonByFinalState/
	JSONDIR=/afs/cern.ch/work/v/vischia/private/results/HIG-13-026/tempjsonByFinalState_5315/
	
	for i in 180 200 220 250 300
	  do
	   for chanList in emu ee mumu
	     do
	     prepareChHiggsDatacards --in ${BASEDIR}plotter_${chanList}_forSystTable_def.root --out datacardsByDecayMhmax/${i}/ --suffix mhmax --json ${JSONDIR}${i}_tb.json --noPowers --histo finalevtflow2btags --bin 1 --ch ${chanList} & 
	     prepareChHiggsDatacards --in ${BASEDIR}plotter_${chanList}_all-samplesForDatacards_finalevtflow_norm_def.root --out datacardsByDecayScan/${i}/ --suffix scan --json ${JSONDIR}${i}_tb.json --noPowers --histo finalevtflow2btags --bin 1 --ch ${chanList} & 
	   done
	done
	
    elif [ "${2}" = "mhmodp" ]; then
	mkdir -p ${BASEDIR}outputs/
	mkdir -p ${BASEDIR}outputs/tables_tanb5/
	mkdir -p ${BASEDIR}outputs/tables_tanb30/
	
	for chanList in emu ee mumu
	  do
	  runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}outputs/ --json data/chhiggs/signal_mhmodp_tanb5_fy_x10.json --outFile ${BASEDIR}outputs/plotter_${chanList}_all-samples_mhmodp_finalevtflow_tanb5.root --noPlot --noPowers  --onlyStartWith ${chanList}_evtflow
	  mv ${BASEDIR}outputs/*tex ${BASEDIR}outputs/tables_tanb5/ 
	  runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}outputs/ --json data/chhiggs/signal_mhmodp_tanb30_fy_x10.json --outFile ${BASEDIR}outputs/plotter_${chanList}_all-samples_mhmodp_finalevtflow_tanb30.root --noPlot --noPowers  --onlyStartWith ${chanList}_evtflow
	  mv ${BASEDIR}outputs/*tex ${BASEDIR}outputs/tables_tanb30/ 
	  runPlotter --iLumi 19702 --inDir ${BASEDIR} --outDir ${BASEDIR}outputs/ --json data/chhiggs/all-samples_mhmodp.json --outFile ${BASEDIR}outputs/plotter_${chanList}_all-samples_mhmodp_finalevtflow.root --noPlot --noPowers  --onlyStartWith ${chanList}_finalevtflow2btags
	## # 
	## # mkdir -p ${BASEDIR}outputs/datacardsMhmodp_tanb5/180
	## # mkdir -p ${BASEDIR}outputs/datacardsMhmodp_tanb5/200 
	## # mkdir -p ${BASEDIR}outputs/datacardsMhmodp_tanb5/220 
	## # mkdir -p ${BASEDIR}outputs/datacardsMhmodp_tanb5/250 
	## # mkdir -p ${BASEDIR}outputs/datacardsMhmodp_tanb5/300
	## # 
	  mkdir -p ${BASEDIR}outputs/datacardsMhmodp_tanb30/180
	  mkdir -p ${BASEDIR}outputs/datacardsMhmodp_tanb30/200 
	  mkdir -p ${BASEDIR}outputs/datacardsMhmodp_tanb30/220 
	  mkdir -p ${BASEDIR}outputs/datacardsMhmodp_tanb30/250 
	  mkdir -p ${BASEDIR}outputs/datacardsMhmodp_tanb30/300
	## # 
	  for i in 180 200 220 250 300
	    do
	## #   JSONDIR=/afs/cern.ch/work/v/vischia/private/results/HIG-13-026/tempjsonByFinalState_5315_tanb5/
	## #   prepareChHiggsDatacards --in ${BASEDIR}outputs/plotter_{chanList}_all-samples_mhmodp_finalevtflow.root --out ${BASEDIR}outputs/datacardsMhmodp_tanb5/${i}/ --suffix mhmodp_tanb5 --json ${JSONDIR}${i}_tb.json --noPowers --histo finalevtflow2btags --bin 1 --ch ${chanList} & 
	## # 
	    JSONDIR=/afs/cern.ch/work/v/vischia/private/results/HIG-13-026/tempjsonByFinalState_5315_tanb30/
	    prepareChHiggsDatacards --in ${BASEDIR}outputs/plotter_${chanList}_all-samples_mhmodp_finalevtflow.root --out ${BASEDIR}outputs/datacardsMhmodp_tanb30/${i}/ --suffix mhmodp_tanb30 --json ${JSONDIR}${i}_tb.json --noPowers --histo finalevtflow2btags --bin 1 --ch ${chanList} & 
	  done
	done
	
    elif [ "${2}" = "put" ]; then
	outputdir=tempDirForNotePlots/
	
	mkdir -p ${outputdir}
	
	for plotList in evtflow met mll mtsum nbjets njets nvertices nverticesUnweighted pte ptjet1eta ptjet1pt ptjet2eta ptjet2pt ptmin ptmu sumpt geq2btagsmet geq2btagsnbjets geq2btagsptlep geq2btagssumpt 
	  do
	  for chanList in emu ee mumu
	    do
	    for formatList in pdf png C
	    if [ "${plotList}" != "evtflow" ]; then
		if [ "${formatList}" == "C" ]; then
		    continue
		fi
	    fi
	    done
	  done
	done
	cp ${BASEDIR}plots/${chanList}_${plotList}.${formatList} ${outputdir}
	
	scp -r ${outputdir} lnlip02.lip.pt:~/
	
	rm -rf ${outputdir} 
	
    fi

fi

exit 0