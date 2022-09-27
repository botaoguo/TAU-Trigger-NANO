# write a script for postproc
from summaryProducer import *
from selectionFilter import *
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from importlib import import_module
import os
import sys
import ROOT
import argparse

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
    

    p = PostProcessor(".", files, "1", "keep_and_drop.txt", modules=[summary2017MC(), selection2017MC()], provenance=True)
    p.run()
    print("Done !")

if __name__ == '__main__':
    main()
