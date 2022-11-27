# event tuple for 2018
nohup python postproc.py --input /eos/cms/store/group/phys_tau/TauFW/nanoV10/Run2_2018/SingleMuon_Run2018B/nano_18.root --isMC 0 --era 2018 > nohup.log 2>&1 &


# MC sample
/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv2-106X_upgrade2018_realistic_v15_L1v1-v1/NANOAODSIM

nohup python postproc.py --input nanoMC-1.root --isMC 1 --era 2018 > nohup.log 2>&1 &
# run event tuple using MC

nohup python postproc.py --input nanoAOD-2018MC-2.root --isMC 1 --era 2018 > nohup.log 2>&1 &

# run event tuple using data

nohup python postproc.py --input nanoAOD-2018data-1.root --isMC 0 --era 2018 > nohup.log 2>&1 &

# activate ROOT env to run the script

##############
# skim tuple #
##############
nohup python skimTuple.py --input nanoAOD-2018MC_Skim.root --config ./2018trigger.json --selection DeepTau --output test.root --type mc --pu PileupHistogram-goldenJSON-13tev-2018-66000ub-99bins.root > nohup-skim.log 2>&1 &

# test for MC
nohup python skimTuple.py --input nanoAOD-2018MC-2_Skim.root --config ./2018trigger.json --selection DeepTau --output nanoAOD-2018MC-2_Skim-skimtuple.root --type mc --pu PileupHistogram-goldenJSON-13tev-2018-66000ub-99bins.root > nohup-skim-MC.log 2>&1 &
# test for data
nohup python skimTuple.py --input nanoAOD-2018data-1_Skim.root --config ./2018trigger.json --selection DeepTau --output nanoAOD-2018data-1_Skim-skimtuple.root --type data > nohup-skim-data.log 2>&1 &

##################
# create turn on #
##################
nohup python createTurnOn.py --input-data nanoAOD-2018data-2_Skim-skimtuple.root --input-dy-mc nanoAOD-2018MC-2_Skim-skimtuple.root --output TrunOn > nohup-createTurnOn.log 2>&1 &

###############
# fit turn on #
###############
nohup python fitTurnOn.py --input TurnOn.root --output fitTrunOn --decay-mode 'all,0,1' > nohup-fitTurnOn.log 2>&1 &
