#!/bin/bash

if [[ -z "$1" ]]; then
  echo "You must specify a json file"
  exit
fi

FILE="$(readlink -m $1)"

TEMP="grep 'iLumi' $FILE | wc -l"
NRounds="$(eval $TEMP)"
CWD="$(pwd)"
CWD="$(readlink -m $CWD)"

echo "Found $NRounds rounds in file $1"
if [[ $NRounds == 0 ]]; then
  echo "There are no rounds to submit jobs for"
  exit
fi

if [[ -d ./Results ]]; then
  echo "A results directory already exists."
  echo "Do you want to delete the directory [y/n]: (the script will stop if no is chosen)"
  read answer
  if [[ $answer == 'y' ]]; then
    echo "Deleting ./Results/"
    rm -Rf ./Results
  else
    echo "Terminating script"
    exit
  fi
fi
mkdir ./Results

for (( round=0; round<$NRounds; round++ ))
do
  echo "Processing Round " $round

  TEMP="cp runRound.sh.templ ./Results/runRound$round.sh"
  eval $TEMP

  sed -i -e "s|#ROUND|$round|g" "./Results/runRound$round.sh"
  sed -i -e "s|#CWD|$CWD|g" "./Results/runRound$round.sh"
  sed -i -e "s|#FILE|$FILE|g" "./Results/runRound$round.sh"
done
