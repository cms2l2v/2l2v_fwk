#!/bin/bash

for stauM in {50..500..10}
do
  for neutM in {0..480..10}
  do
    echo "Stau: " $stauM "; Neutralino: " $neutM

    file=""
    printf -v file './Results/SignalPoint_%d%03d.txt' "$stauM" "$neutM"
    if [[ -f $file ]]; then
      combine -M Asymptotic --run=blind $file -n S$stauM-N$neutM
    fi
  done
done
