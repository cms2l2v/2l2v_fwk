#!/bin/bash

# Only Run2012B
# Inclusive
runPlotterFWLite --iEcm 8 --iLumi 4429 --inDir Processed_Inclusive/ --outDir ~/www/TStauStau/MCvsData/2012B/inclusive --outFile plotter_2012B_Inclusive.root --json tstaustau_samples_inclusive_2012B.json
# Exclusive
runPlotterFWLite --iEcm 8 --iLumi 4429 --inDir Processed_Exclusive/ --outDir ~/www/TStauStau/MCvsData/2012B/exclusive --outFile plotter_2012B_Exclusive.root --json tstaustau_samples_exclusive_2012B.json

# Only Run2012C
# Inclusive
runPlotterFWLite --iEcm 8 --iLumi 7147 --inDir Processed_Inclusive/ --outDir ~/www/TStauStau/MCvsData/2012C/inclusive --outFile plotter_2012C_Inclusive.root --json tstaustau_samples_inclusive_2012C.json
# Exclusive
runPlotterFWLite --iEcm 8 --iLumi 7147 --inDir Processed_Exclusive/ --outDir ~/www/TStauStau/MCvsData/2012C/exclusive --outFile plotter_2012C_Exclusive.root --json tstaustau_samples_exclusive_2012C.json

# Only Run2012C
# Inclusive
runPlotterFWLite --iEcm 8 --iLumi 7317 --inDir Processed_Inclusive/ --outDir ~/www/TStauStau/MCvsData/2012D/inclusive --outFile plotter_2012D_Inclusive.root --json tstaustau_samples_inclusive_2012D.json
# Exclusive
runPlotterFWLite --iEcm 8 --iLumi 7317 --inDir Processed_Exclusive/ --outDir ~/www/TStauStau/MCvsData/2012D/exclusive --outFile plotter_2012D_Exclusive.root --json tstaustau_samples_exclusive_2012D.json

# Runs 2012{B,C,D}
# Inclusive
runPlotterFWLite --iEcm 8 --iLumi 18893 --inDir Processed_Inclusive/ --outDir ~/www/TStauStau/MCvsData/2012BCD/inclusive --outFile plotter_2012BCD_Inclusive.root --json tstaustau_samples_inclusive_2012BCD.json
# Exclusive
runPlotterFWLite --iEcm 8 --iLumi 18893 --inDir Processed_Exclusive/ --outDir ~/www/TStauStau/MCvsData/2012BCD/exclusive --outFile plotter_2012BCD_Exclusive.root --json tstaustau_samples_exclusive_2012BCD.json
