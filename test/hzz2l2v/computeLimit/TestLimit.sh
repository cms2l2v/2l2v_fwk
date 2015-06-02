mkdir -p TEST;
cd TEST;
computeLimit --m 400 --blind --histo mt_shapes  --in $PWD/../../plotter.root --syst  --index 10 --json $PWD/../../samples.json  --shapeMin 0 --shapeMax 9999  --bins eq0jets,geq1jets,vbf --systpostfix _13TeV --rebin 8 --dropBckgBelow 0.0 ;
sh combineCards.sh;
combine -M Asymptotic -m 400 --run expected card_combined.dat > COMB.log;
cd ..
