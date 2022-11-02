nohup python postproc.py --input nanoAOD-2018MC.root --isMC 1 --era 2018 > nohup.log 2>&1 &

# activate ROOT env to run the script

nohup python skimTuple.py --input nanoAOD-2018MC_Skim.root --config ./2018trigger.json --selection DeepTau --output test.root --type mc --pu PileupHistogram-goldenJSON-13tev-2018-66000ub-99bins.root > nohup-skim.log 2>&1 &
