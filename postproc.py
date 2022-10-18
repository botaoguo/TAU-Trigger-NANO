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

import FWCore.ParameterSet.Config as cms

#sys.path.append('Common')
#from TriggerConfig import *
import TriggerConfig
#from PhysicsTools.NanoAODTools.postprocessing.tauana.Common.TriggerConfig import *

ROOT.PyConfig.IgnoreCommandLineOptions = True

def main():
    
    parser = argparse.ArgumentParser(description='Post Processing.')
    parser.add_argument('--input', required=True, type=str, help="NANO input")
    #parser.add_argument('--output', required=True, type=str, help="eventTuple output")
    parser.add_argument('--isMC', required=True, type=int, help="judge if isMC")
    parser.add_argument('--era', required=True, type=str, help="")
    args = parser.parse_args()
    print "args = ",args

    isMC = args.isMC
    era = args.era
    files = [ args.input ]
    
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
    Modules = [summary2018MC(), selection2018MC(), tuple2018MC()]

    p = PostProcessor(".", files, "1", 
                      branchsel = "keep_and_drop.txt", 
                      modules= Modules, 
                      provenance=True,
                      outputbranchsel = "output_branch.txt"
    )
    p.run()
    print("Done !")

if __name__ == '__main__':
    main()
