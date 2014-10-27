#!/bin/bash

# Only Run2012D
## LIP
runPlotterFWLite --noPowers --iEcm 8 --iLumi 7329 --inDir /lustre/ncg.ingrid.pt/cmslocal/cbeiraod/LIP_2012D_selection/ --outDir ./LIP_2012DHLT_Selection/ --outFile plotter_Exclusive_LIP_2012DHLT.root --json tstaustau_samples_exclusive_2012D.json --plotExt .png --plotExt .pdf --plotExt .C --plotExt .root
