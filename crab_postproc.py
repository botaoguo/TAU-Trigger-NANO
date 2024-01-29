# write a script for postproc
from summaryProducer import *
from selectionFilter import *
from tupleProducer import *
# this takes care of converting the input files from CRAB
from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import inputFiles, runsAndLumis
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from importlib import import_module
import os
import sys
import ROOT
import argparse

import FWCore.ParameterSet.Config as cms

#sys.path.append('Common')
#from TriggerConfig import *
import TriggerConfig
#from PhysicsTools.NanoAODTools.postprocessing.tauana.Common.TriggerConfig import *

ROOT.PyConfig.IgnoreCommandLineOptions = True

def main():
    
    parser = argparse.ArgumentParser(description='Post Processing.')
    #parser.add_argument('--input', required=True, type=str, help="NANO input")
    #parser.add_argument('--output', required=True, type=str, help="eventTuple output")
    parser.add_argument('--isMC', required=True, type=int, help="judge if isMC")
    parser.add_argument('--era', required=True, type=str, help="")
    args = parser.parse_args()
    print "args = ",args

    isMC = args.isMC
    era = args.era
    #files = [ args.input ]
    
    #trigFile = './2017trigger.json'
    #hltPaths, tagHltPaths = TriggerConfig.LoadAsVPSet(trigFile)

    #print("hltPaths: ", hltPaths)

    #import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
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
    else:
        if era == '2016':
            Modules = [summary2016data(), selection2016data(), tuple2016data()]
        elif era == '2017':
            Modules = [summary2017data(), selection2017data(), tuple2017data()]
        elif era == '2018':
            Modules = [summary2018data(), selection2018data(), tuple2018data()]
        elif era == '2022':
            Modules = [summary2022data(), selection2022data(), tuple2022data()]
        else:
            raise RuntimeError("Please check the right Year!")
        

    p = PostProcessor(".", inputFiles(), "1", 
                      branchsel = "keep_and_drop.txt", 
                      modules= Modules, 
                      provenance=True,
                      fwkJobReport=True,
                      jsonInput=runsAndLumis(),
                      outputbranchsel = "output_branch.txt"
    )
    p.run()
    print("Done !")

if __name__ == '__main__':
    main()
