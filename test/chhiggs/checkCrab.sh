#!/bin/bash


for CRABDIR in ${1}/FARM/inputs/crab_*; do
    echo "Executing:"
    echo "crab status -d ${CRABDIR}"
    echo "----------------------------------------------------------------------------------------------"
    crab status -d ${CRABDIR}
done
exit 0