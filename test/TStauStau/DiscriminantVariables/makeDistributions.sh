#!/bin/bash

#Pass0
#makeDistributions --json tstaustau_samples_inclusive_2012BCD.json --variables variables.json --inDir ~/local-area/Results/ --signalSelection "stauMass-neutralinoMass==100" --plotExt .png --plotExt .pdf --plotExt .root --baseSelection "selected" --outDir DeltaM100/Pass0/
#makeDistributions --json tstaustau_samples_inclusive_2012BCD.json --variables variables.json --inDir ~/local-area/Results/ --signalSelection "stauMass-neutralinoMass==80" --plotExt .png --plotExt .pdf --plotExt .root --baseSelection "selected" --outDir DeltaM80/Pass0/
#makeDistributions --json tstaustau_samples_inclusive_2012BCD.json --variables variables.json --inDir ~/local-area/Results/ --signalSelection "stauMass-neutralinoMass==40" --plotExt .png --plotExt .pdf --plotExt .root --baseSelection "selected" --outDir DeltaM40/Pass0/

#Pass1
#makeDistributions --json tstaustau_samples_inclusive_2012BCD.json --variables variables.json --inDir ~/local-area/Results/ --signalSelection "stauMass-neutralinoMass==100" --plotExt .png --plotExt .pdf --plotExt .root --baseSelection "selected&&MT2>90" --outDir DeltaM100/Pass1/
#makeDistributions --json tstaustau_samples_inclusive_2012BCD.json --variables variables.json --inDir ~/local-area/Results/ --signalSelection "stauMass-neutralinoMass==80" --plotExt .png --plotExt .pdf --plotExt .root --baseSelection "selected&&MT2>90" --outDir DeltaM80/Pass1/
#makeDistributions --json tstaustau_samples_inclusive_2012BCD.json --variables variables.json --inDir ~/local-area/Results/ --signalSelection "stauMass-neutralinoMass==40" --plotExt .png --plotExt .pdf --plotExt .root --baseSelection "selected&&met.Et()>160" --outDir DeltaM40/Pass1/


makeDistributions --json tstaustau_samples_inclusive_2012BCD.json --variables variables.json --inDir ~/local-area/Results/ --signalSelection "stauMass-neutralinoMass==100" --extraSignal "stauMass-neutralinoMass==40" --plotExt .png --plotExt .pdf --plotExt .root --baseSelection "selected" --outDir DeltaM100/Pass0/
