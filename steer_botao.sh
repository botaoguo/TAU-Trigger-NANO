nohup python postproc.py --input nanoAOD.root --isMC 1 --era 2017 > nohup.log 2>&1 &

# activate ROOT env to run the script

nohup python skimTuple.py --input nanoAOD-2018MC_Skim.root --config ./2018trigger.json --selection DeepTau --output test.root --type mc > nohup-skim.log 2>&1 &
