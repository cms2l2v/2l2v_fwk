#!/bin/bash

# Fit the ttbar cross section using shape
BASEDIR=~/www/13TeV_xsec_plots/
JSONFILEDILEPTON=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_xsec_dileptons.json
JSONFILELEPTAU=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_xsec_leptontau.json

# histo choices for dilepton: finalnbjets, finalmt, eventflow(bin6)
# histo choices for leptau: finaltaur, finaltaupolarization, finaldphilepmet, finaldphitaumet, finaldphileptau, finaltaupt, eventflowslep(bin6)
# Remember to blank the BINS* variables when not using eventflow histograms ;)

###500 1000 5000 10000
for LUMI in 16.1
do
    DIR=${BASEDIR}${LUMI}/
    for HISTODILEPTON in  eventflow
    do
        OUTDILEPTON=~/www/13TeV_xsec/${LUMI}/dileptons/${HISTODILEPTON}/
        mkdir -p ${OUTDILEPTON}
        cp ~/www/HIG-13-026/index.php ${OUTDILEPTON}
        BINSDILEPTON=""
        BINSLEPTAU=""
        if [ "${HISTODILEPTON}" == "eventflow" ]
        then
            BINSDILEPTON=" --bins 6 "
            BINSLEPTAU=" --bins 6 "
        fi
        # Dileptons
        fitTTbarCrossSection --out ${OUTDILEPTON} --in ${DIR}plotter_dilepton.root --json ${JSONFILEDILEPTON} --histo ${HISTODILEPTON} ${BINSDILEPTON} --ch emu,ee,mumu
        cp ~/www/HIG-13-026/index.php ${OUTDILEPTON}
        
    done
    
    #for HISTOLEPTAU in finaltaur finaltaupolarization finaldphilepmet finaldphitaumet finaldphileptau finaltaupt eventflowslep finalnbjets
    for HISTOLEPTAU in finaltaupolarization finaldphilepmet eventflowslep finalnbjets
    do
        OUTLEPTAU=~/www/13TeV_xsec/${LUMI}/leptau/${HISTOLEPTAU}/
        mkdir -p ${OUTLEPTAU}
        cp ~/www/HIG-13-026/index.php ${OUTDILEPTAU}
        BINSDILEPTON=""
        BINSLEPTAU=""
        if [ "${HISTOLEPTAU}" == "eventflowslep" ]
        then
            BINSDILEPTON=" --bins 6 "
            BINSLEPTAU=" --bins 6 "
        fi
    # Leptontau
        fitTTbarCrossSection --out ${OUTLEPTAU}   --in ${DIR}plotter_ltau.root     --json ${JSONFILELEPTAU}   --histo ${HISTOLEPTAU} ${BINSLEPTAU} --ch singlemu,singlee
        cp ~/www/HIG-13-026/index.php ${OUTLEPTAU}
    done

    
    for HISTOFORRATIO in eventflowslep eventflow
    do
        OUTRATIO=~/www/13TeV_xsec/${LUMI}/ratio/${HISTOFORRATIO}/
        mkdir -p ${OUTRATIO}
        cp ~/www/HIG-13-026/index.php ${OUTRATIO}
        BINS=""
        INFILE=""
        JSON=""
        CH=""
        if [ "${HISTOFORRATIO}" == "eventflowslep" ]
        then
            BINS=" --bins 6 "
            INFILE="plotter_ltau.root"
            JSON=${JSONFILELEPTAU}
            CH=" --ch singlemu,singlee"
        fi
        if [ "${HISTOFORRATIO}" == "eventflow" ]
        then
            BINS=" --bins 6 "
            INFILE="plotter_dilepton.root"
            JSON=${JSONFILEDILEPTON}
            CH=" --ch emu,ee,mumu"
        fi
        fitTTbarCrossSection --out ${OUTRATIO} --in ${DIR}${INFILE} --json ${JSON} --histo ${HISTOFORRATIO} ${BINS} ${CH}
    done
done