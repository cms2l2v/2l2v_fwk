jsonTools --in ../../data/all_samples.json --keepFromFile SUMMARY --out mergeEDMs.json
mergeEDMtuples.py -i  /storage/data/cms/store/user/quertenmont/14_03_20_2l2nu_EDMtuples -j mergeEDMs.json -o  $PWD -D True
