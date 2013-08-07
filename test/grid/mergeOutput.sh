#!/bin/bash

crab_working_dir=${1}
nsplit=${2}
output_dir=${3}

crab -status -c ${crab_working_dir}
crab -get -USER.dontCheckSpaceLeft=1 -c ${crab_working_dir}
crab -report -c ${crab_working_dir}
mv ${crab_working_dir}/res/lumiSummary.json ${crab_working_dir}.json
mv ${crab_working_dir}/res/missingLumiSummary.json ${crab_working_dir}-missing.json
job_outputs=(`ls ${crab_working_dir}/res/*.root`)
totalOutputs=${#job_outputs[@]} 

if [ -z ${nsplit} ]; then
    echo "Setting nsplit=1"
    nsplit=1
fi

if [ -z ${output_dir} ]; then
    echo "Setting output to local"
    output_dir="./"
fi


let "step=${totalOutputs}/${nsplit}"
let "totalstep=${nsplit}-1"
echo "${totalOutputs} files will be hadded in groups of ${step}"
for i in `seq 0 ${totalstep}`; do
    let "startOutput=${i}*step"
    let "endOutput=(${i}+1)*step-1"

    if [[ ${i} -eq ${totalstep} ]]; then
	let "endOutput=${totalOutputs}-1"
    fi
    echo "${startOutput}-${endOutput}"

    outputList=""
    for j in `seq ${startOutput} ${endOutput}`; do
	newFile=${job_outputs[${j}]}
	if [ -e ${newFile} ]; then
	    outputList="${outputList} ${newFile}"
	fi
    done
    hadd -f -k /tmp/psilva/${crab_working_dir}_${i}.root ${outputList}
    cmsStage -f /tmp/psilva/${crab_working_dir}_${i}.root ${output_dir} 
    rm /tmp/psilva/${crab_working_dir}_${i}.root
done


