#/bin/bash
PDs=(
    BTag2012A BTag2012B BTag2012C BTag2012D
#    DoubleElectron2012A DoubleElectron2012B DoubleElectron2012C DoubleElectron2012D
#    DoubleMu2012A DoubleMu2012B DoubleMu2012C DoubleMu2012D
#    MuEG2012A MuEG2012B MuEG2012C MuEG2012D
)
#mbxsec=70300
#PUJSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/PileUp/pileup_latest.txt

for j in ${PDs[@]}; do 
    #crab -report -c Data8TeV_${j}
    #mv Data8TeV_${j}/res/missingLumiSummary.json Data8TeV_${j}.missing.json
    #mv Data8TeV_${j}/res/lumiSummary.json Data8TeV_${j}.json
    lumiCalc2.py overview -i Data8TeV_${j}.json > Data8TeV_${j}.lumi
    #pileupCalc.py -i Data8TeV_${j}.json --inputLumiJSON ${PUJSON} --calcMode observed --minBiasXsec ${mbxsec} --maxPileupBin 100 --numPileupBins 100  Data8TeV_${j}_targetpu_${mbxsec}.root &
done