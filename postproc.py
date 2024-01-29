# write a script for postproc
from summaryProducer import *
from selectionFilter import *
from tupleProducer import *
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from importlib import import_module
import os
import sys
import ROOT
import argparse

import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Config as cms

#sys.path.append('Common')
#from TriggerConfig import *
import TriggerConfig
#from PhysicsTools.NanoAODTools.postprocessing.tauana.Common.TriggerConfig import *

ROOT.PyConfig.IgnoreCommandLineOptions = True

def main():
    
    parser = argparse.ArgumentParser(description='Post Processing.')
    parser.add_argument('--input', required=False, type=str, help="NANO input")
    parser.add_argument('--inputFileList', required=False, type=str, help="NANO input file list")
    parser.add_argument('--output', required=True, type=str, help="eventTuple output")
    parser.add_argument('--isMC', required=True, type=int, help="judge if isMC")
    # parser.add_argument('--nanoVer', required=False, type=int, default=9, help="NanoAOD Version")
    parser.add_argument('--era', required=True, type=str, help="")
    args = parser.parse_args()
    print "args = ",args

    if (args.input is None) and (args.inputFileList is None):
        raise RuntimeError("Please check the input!")
    if (args.input is not None) and (args.inputFileList is not None):
        raise RuntimeError("Please check the input!")

    isMC = args.isMC
    output = args.output
    era = args.era
    # nanoVer = args.nanoVer
    if args.input:
        files = [ args.input ]
    if args.inputFileList:
        f = open(args.inputFileList, "r")
        files = f.read().splitlines()
    print(files)
    # exit(0)
    
    #trigFile = './2017trigger.json'
    #hltPaths, tagHltPaths = TriggerConfig.LoadAsVPSet(trigFile)

    #print("hltPaths: ", hltPaths)

    import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
    #process = cms.Process('NANO')
    #process.metSequence = cms.Sequence()
    #process.hltfilter = hlt.hltHighLevel.clone(
    #    TriggResultsTag = cms.InputTag("TriggerResult", "", "HLT"),
    #    HLTPaths = [p + '*' for p in tagHltPaths],
    #    andOr = cms.bool(True),
    #    throw = cms.bool(True)
    #)

    #process.patTrigger = cms.EDProducer("PATTriggerObjectStandAloneUnpacker",
    #    patTriggerObjectStandAlone = cms.InputTag("")
    #)
    
    #Modules = [tuple2017MC()]
    if isMC:
        jsoninput = None,
        if era == '2016':
            Modules = [summary2016MC(), selection2016MC(), tuple2016MC()]
        elif era == '2017':
            Modules = [summary2017MC(), selection2017MC(), tuple2017MC()]
        elif era == '2018':
            Modules = [summary2018MC(), selection2018MC(), tuple2018MC()]
        elif era == '2022':
            Modules = [summary2022MC(), selection2022MC(), tuple2022MC()]
        else:
            raise RuntimeError("Please check the right Year!")
        p = PostProcessor(output, files, "1", 
                        branchsel = "keep_and_drop.txt", 
                        modules= Modules, 
                        provenance=True,
                        outputbranchsel = "output_branch.txt"
        )

    else:
        if era == '2016':
            Modules = [summary2016data(), selection2016data(), tuple2016data()]
            lumisToProcess = cms.untracked.VLuminosityBlockRange( LumiList.LumiList(filename="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt").getCMSSWString().split(',') )
        elif era == '2017':
            Modules = [summary2017data(), selection2017data(), tuple2017data()]
            lumisToProcess = cms.untracked.VLuminosityBlockRange( LumiList.LumiList(filename="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt").getCMSSWString().split(',') )
        elif era == '2018':
            Modules = [summary2018data(), selection2018data(), tuple2018data()]
            lumisToProcess = cms.untracked.VLuminosityBlockRange( LumiList.LumiList(filename="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt").getCMSSWString().split(',') )
        elif era == '2022':
            Modules = [summary2022data(), selection2022data(), tuple2022data()]
            lumisToProcess = cms.untracked.VLuminosityBlockRange( LumiList.LumiList(filename="./Cert_Collisions2022_355100_362760_GoldenJSON.txt").getCMSSWString().split(',') )
        else:
            raise RuntimeError("Please check the right Year!")
    
        runsAndLumis_special = {}
        for l in lumisToProcess:
            if "-" in l:
                start, stop = l.split("-")
                rstart, lstart = start.split(":")
                rstop, lstop = stop.split(":")
            else:
                rstart, lstart = l.split(":")
                rstop, lstop = l.split(":")
            if rstart != rstop:
                raise Exception(
                    "Cannot convert '%s' to runs and lumis json format" % l)
            if rstart not in runsAndLumis_special:
                runsAndLumis_special[rstart] = []
            runsAndLumis_special[rstart].append([int(lstart), int(lstop)])
        
        jsoninput = runsAndLumis_special
        
        p = PostProcessor(output, files, "1", 
                        branchsel = "keep_and_drop.txt", 
                        modules= Modules, 
                        jsonInput=jsoninput,
                        provenance=True,
                        outputbranchsel = "output_branch.txt"
        )
    p.run()
    print("Done !")

if __name__ == '__main__':
    main()
