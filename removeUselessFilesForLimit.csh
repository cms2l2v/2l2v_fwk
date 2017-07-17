#!/bin/csh

set list=`cat list_file_to_remove_for_limit`

foreach i ($list)
	echo $i
	rm $i
end

sed -i '/extractPhotonWeights_UsingBins.cc/d' bin/BuildFile.xml
sed -i '/fitTTbarCrossSection.cc/d' bin/BuildFile.xml
sed -i '/runFixedPlotter.cc/d' bin/BuildFile.xml
sed -i '/computeLeptonEfficency.cc/d' bin/BuildFile.xml
sed -i '/runPhotonZClosure.cc/d' bin/BuildFile.xml
sed -i '/runHZZ2l2vAnalysis.cc/d' bin/BuildFile.xml
sed -i '/extractFRWeights.cc/d' bin/BuildFile.xml
sed -i '/extractPhotonWeights.cc/d' bin/BuildFile.xml
sed -i '/runEfficencyMiniAODExample.cc/d' bin/BuildFile.xml
sed -i '/MELA/d' bin/BuildFile.xml
sed -i '/SVfitStandalone/d' bin/BuildFile.xml
