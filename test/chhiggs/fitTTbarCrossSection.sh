#!/bin/bash

# Fit the ttbar cross section using shape
DIR=~/www/13TeV_topsel_fix_ala9/
JSONFILEDILEPTON=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_xsec_dileptons.json
JSONFILELEPTAU=$CMSSW_BASE/src/UserCode/llvv_fwk/data/chhiggs/phys14_xsec_leptontau.json

# histo choices for dilepton: finalnbjets, finalmt, eventflow(bin6)
# histo choices for leptau: finaltaur, finaltaupolarization, finaldphilepmet, finaldphitaumet, finaldphileptau, finaltaupt, eventflowslep(bin6)
# Remember to blank the BINS* variables when not using eventflow histograms ;)
HISTODILEPTON=finalmt
HISTOLEPTAU=finaltaupolarization
#HISTODILEPTON=eventflow
#HISTOLEPTAU=eventflowslep
BINSDILEPTON=""
BINSLEPTAU=""
#BINSDILEPTON=" --bins 6 "
#BINSLEPTAU=" --bins 6 "

OUTDILEPTON=~/www/13TeV_xsec/dileptons/
OUTLEPTAU=~/www/13TeV_xsec/leptau/



# Dileptons
fitTTbarCrossSection --out ${OUTDILEPTON} --in ${DIR}plotter_dilepton.root --json ${JSONFILEDILEPTON} --histo ${HISTODILEPTON} ${BINSDILEPTON} --ch emu,ee,mumu

# Leptontau
fitTTbarCrossSection --out ${OUTLEPTAU}   --in ${DIR}plotter_ltau.root     --json ${JSONFILELEPTAU}   --histo ${HISTOLEPTAU} ${BINSLEPTAU} --ch singlemu,singlee

cp ~/www/HIG-13-026/index.php ${OUTDILEPTON}
cp ~/www/HIG-13-026/index.php ${OUTLEPTAU}