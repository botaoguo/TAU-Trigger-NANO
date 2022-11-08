from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config

config = Configuration()

config.section_("General")
config.General.requestName = 'DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8'
config.General.workArea = "crab2018-test3"
config.General.transferLogs = False
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'crab_script.sh'
#config.JobType.numCores = 1
#config.JobType.maxMemoryMB = 2000
# hadd nano will not be needed once nano tools are in cmssw
config.JobType.inputFiles = ['crab_postproc.py', '../scripts/haddnano.py', 'summaryProducer.py', 'selectionFilter.py', 'tupleProducer.py', 'keep_and_drop.txt', 'output_branch.txt', 'triggerDescriptor.py', 'TriggerConfig.py', 'DAS_filesearch.py', '2018trigger.json', '2018filterName.txt', '2018filterHash.txt']
config.JobType.sendPythonFolder = True
config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv2-106X_upgrade2018_realistic_v15_L1v1-v1/NANOAODSIM'
#config.Data.inputDBS = 'phys03'
config.Data.inputDBS = 'global'
#config.Data.splitting = 'Automatic'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 1
#config.Data.totalUnits = -1 
config.Data.ignoreLocality = True
config.Data.allowNonValidInputDataset = True

#config.Data.outLFNDirBase = '/store/user/boguo/NanoPost'
config.Data.publication = False
config.Data.outputDatasetTag = 'DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8'
config.section_("Site")
config.Site.storageSite = "T2_IT_Pisa"
config.Site.whitelist = ["T2_IT_Pisa"]
#config.Site.storageSite = "T2_CH_CERN"
#config.section_("User")
#config.User.voGroup = 'dcms'
