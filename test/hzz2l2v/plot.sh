#used to take additional arguments from the command line
arguments=''; for var in "$@"; do arguments=$arguments" "$var; done
if [ $# -ge 1 ]; then echo "Additional arguments will be considered: "$arguments ;fi

runPlotter --iEcm 13 --iLumi 20.791 --inDir $CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v/results/ --outDir $CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v/plots/ --outFile $CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v/plotter.root  --json $CMSSW_BASE/src/UserCode/llvv_fwk/test/hzz2l2v/samples.json --no2D $arguments
