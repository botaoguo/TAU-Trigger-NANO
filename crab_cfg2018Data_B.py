from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config

config = Configuration()

config.section_("General")
config.General.requestName = 'SingleMuon_Run2018UL_B'
config.General.workArea = "crab2018-0427"
config.General.transferLogs = False
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'crab_script2018Data.sh'
config.JobType.maxJobRuntimeMin = 2750
#config.JobType.numCores = 1
config.JobType.maxMemoryMB = 2500
# hadd nano will not be needed once nano tools are in cmssw
config.JobType.inputFiles = ['crab_postproc.py', '../scripts/haddnano.py', 'summaryProducer.py', 'selectionFilter.py', 'tupleProducer.py', 'keep_and_drop.txt', 'output_branch.txt', 'triggerDescriptor.py', 'TriggerConfig.py', '2018trigger.json']
config.JobType.sendPythonFolder = True
config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
config.Data.inputDataset = '/SingleMuon/Run2018B-UL2018_MiniAODv1_NanoAODv2-v2/NANOAOD'
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
config.Data.outputDatasetTag = 'SingleMuon_Run2018UL_B'
config.section_("Site")
config.Site.storageSite = "T2_IT_Pisa"
config.Site.whitelist = ["T2_IT_Pisa"]
#config.Site.storageSite = "T2_CH_CERN"
#config.section_("User")
#config.User.voGroup = 'dcms'
