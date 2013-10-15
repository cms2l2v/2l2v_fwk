#!/bin/bash

inputdir=/afs/cern.ch/work/v/vischia/private/code/tau_dilepton/chhiggs_5311/plots/
outputdir=tempDirForNotePlots/

mkdir -p ${outputdir}


cp ${inputdir}emu_evtflow.pdf                  ${outputdir}
cp ${inputdir}emu_met.pdf		       ${outputdir}
cp ${inputdir}emu_mll.pdf		       ${outputdir}
cp ${inputdir}emu_mtsum.pdf		       ${outputdir}
cp ${inputdir}emu_nbjets.pdf		       ${outputdir}
cp ${inputdir}emu_njets.pdf		       ${outputdir}
cp ${inputdir}emu_nvertices.pdf		       ${outputdir}
cp ${inputdir}emu_nverticesUnweighted.pdf      ${outputdir}
cp ${inputdir}emu_pte.pdf		       ${outputdir}
cp ${inputdir}emu_ptjet1eta.pdf		       ${outputdir}
cp ${inputdir}emu_ptjet1pt.pdf		       ${outputdir}
cp ${inputdir}emu_ptjet2eta.pdf		       ${outputdir}
cp ${inputdir}emu_ptjet2pt.pdf		       ${outputdir}
cp ${inputdir}emu_ptmin.pdf		       ${outputdir}
cp ${inputdir}emu_ptmu.pdf		       ${outputdir}
cp ${inputdir}emu_sumpt.pdf		       ${outputdir}
cp ${inputdir}emu_evtflow.png		       ${outputdir}
cp ${inputdir}emu_met.png		       ${outputdir}
cp ${inputdir}emu_mll.png		       ${outputdir}
cp ${inputdir}emu_mtsum.png		       ${outputdir}
cp ${inputdir}emu_nbjets.png		       ${outputdir}
cp ${inputdir}emu_njets.png		       ${outputdir}
cp ${inputdir}emu_nvertices.png		       ${outputdir}
cp ${inputdir}emu_nverticesUnweighted.png      ${outputdir}
cp ${inputdir}emu_pte.png		       ${outputdir}
cp ${inputdir}emu_ptjet1eta.png		       ${outputdir}
cp ${inputdir}emu_ptjet1pt.png		       ${outputdir}
cp ${inputdir}emu_ptjet2eta.png		       ${outputdir}
cp ${inputdir}emu_ptjet2pt.png		       ${outputdir}
cp ${inputdir}emu_ptmin.png		       ${outputdir}
cp ${inputdir}emu_ptmu.png		       ${outputdir}
cp ${inputdir}emu_sumpt.png                    ${outputdir}       


scp -r ${outputdir} lnlip02.lip.pt:~/

rm -rf ${outputdir} 
