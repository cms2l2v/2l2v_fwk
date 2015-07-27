cp test/chhiggs/results_ttbar/*SingleMuon*.json > myjson.json
<edit myjson.json to remove the additional "}{"
python scripts/lumi/lcr2/lcr2.py -i myjson.json