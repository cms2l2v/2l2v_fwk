#!/bin/bash

OUTDIR=/lustre/ncg.ingrid.pt/cmslocal/cbeiraod/Distributions
INDIR=/lustre/ncg.ingrid.pt/cmslocal/cbeiraod/NewObject

makeDistributions --json tstaustau_samples_full.json --variables variables.json --inDir $INDIR/ --outDir $OUTDIR/ --signalSelection "stauMass-neutralinoMass==50" --extraSignal "stauMass-neutralinoMass==200" --plotExt .png --plotExt .root --baseSelection "selected" --iLumi 19672 --pointVar "stauMass*1000+neutralinoMass" --sigXSec 0.1 --sigNInitEvents 10000 --unblind
