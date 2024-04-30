# Setup the framework
1. This tool is based on the CMS nanoAOD tool: https://github.com/cms-nanoAOD/nanoAOD-tools, release a CMSSW (take CMSSW_10_6_29 as exmaple),
```
cd $CMSSW_BASE/src
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
cd PhysicsTools/NanoAODTools
cmsenv
scram b
```
2. Then clone this framework using command below
```
git clone https://github.com/botaoguo/TAU-Trigger-NANO.git
```

# TAU-Trigger-NANO
1. postprocessing the nanoAOD file
```
cd TAU-Trigger-NANO
```
Then you can using two type of input parameter to run the processing, notice that file.txt should contains one root file each line, like **/your/path/your_input.root** in each line of the txt file

```
python postproc.py --input nanoAOD.root --isMC 1 --era 2022 --output ./
python postproc.py --inputFileList file.txt --isMC 1 --era 2022 --output ./
```
The ntuple you get would be "nanoAOD_Skim.root"


2. loop the tuple to match Trigger Object filterbit and some other cut
In this case, you should activate a env that contains ROOT package and make sure that **RDF**
, **numpy**, **matplotlib**, **sklearn** and other packages imported in the file **fitTurnOn.py** could be used.
run the command below to do trigger object match and get the "numerator and denominator" in efficiency measurements
```
python skimTuple_2022preEE.py --input nanoAOD_Skim.root --output sk_mc --selection DeepTau --type mc --pudata pileupfile_data.root --pumc pileupfile_mc.root
python skimTuple_2022preEE.py --input nanoAOD_data_Skim.root --output sk_data --selection DeepTau --type data
```
Then you will get two files, one for mc, another for data.

3. plot the TurnOn curves
run the command below
```
python createTurnOn_multi.py --input-data sk_data.root --input-mc sk_mc.root --output TurnOn
```

4. fit the TurnOn curves
run the command below
```
python fitTurnOn_multi.py --intput TurnOn.root --output fitTurnOn
```