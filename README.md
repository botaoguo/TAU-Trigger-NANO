# TAU-Trigger-NANO
```
python postproc.py --input nanoAOD.root --isMC 1 --era 2018
```

This tool is based on the CMS nanoAOD tool:
https://github.com/cms-nanoAOD/nanoAOD-tools

The input file nanoAOD.root can be found in DAS : 
/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv2-106X_upgrade2018_realistic_v15_L1v1-v1/NANOAODSIM

Then, it will output a event tuple which named nanoAOD_Skim.root

Activate ROOT env in conda base env to run the skimTuple.py
```
python skimTuple.py --input nanoAOD_Skim.root --config ./2018trigger.json --selection DeepTau --output test.root --type mc --pu PileupHistogram-goldenJSON-13tev-2018-66000ub-99bins.root
```

Run on CRAB
```
crab submit -c crab_cfg.py
```
