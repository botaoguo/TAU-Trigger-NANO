# write a Module for summary producer
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import *
import ROOT as R

ROOT.PyConfig.IgnoreCommandLineOptions = True

class summaryProducer(Module):
    def __init__(self, isMC, era):
        #self.writeHistFile = True

        self.isMC = isMC
        self.era = era
        self.tot_evt = 0
        self.ctr_evt_processed = 0
        self.tot_evt_genweight = 0
        
        # cutflow hist
        #self.cutflow_hist = R.TH1F('cutflow','Cutflow',20,0,20)
        #self.cutflow_assigned = 0
        pass
    
    def beginJob(self, histFile=None, histDirName=None):
        #Module.beginJob(self, histFile, histDirName)
        #self.cutflow_hist = R.TH1F('cutflow','Cutflow',20,0,20)
        #self.addObject(self.cutflow_hist)
        #self.outtree.branch("Cutflow",   "B");
        pass
    
    def endJob(self):
        pass
    
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.tot_evt = inputTree.GetEntries()
        print("total event : {}".format(self.tot_evt))
        self.out = wrappedOutputTree

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        #self.cutflow_hist.Write()
        pass

    def fill_cut(self, cutnm):
        _ibin = self.cutflow_hist.GetXaxis().FindBin(cutnm)
        if _ibin == -1:
            self.cutflow_assigned += 1
            self.cutflow_hist.GetXaxis().SetBinLabel(self.cutflow_assigned, cutnm)
            self.cutflow_hist.SetBinContent( self.cutflow_assigned, 1 )
        else:
            self.cutflow_hist.SetBinContent( _ibin, self.cutflow_hist.GetBinContent(_ibin) + 1 )

    def analyze(self, event):
        # process event, return True (go to next module) or False (fail, go to next event)
        
        #electrons = Collection(event, "Electron")
        #muons = Collection(event, "Muon")
        #jets = Collection(event, "Jet")
        #taus = Collection(event, "Tau")
        generator = Object(event, "Generator")

        self.ctr_evt_processed += 1
        self.tot_evt_genweight += generator.weight
        #self.fill_cut('tot_evt')

        if self.ctr_evt_processed % 10000 ==0:
            print("Processed evt = {}".format(self.ctr_evt_processed))

        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

summary2016MC = lambda : summaryProducer(True,"2016")
summary2017MC = lambda : summaryProducer(True,"2017")
summary2018MC = lambda : summaryProducer(True,"2018")
summary2016data = lambda : summaryProducer(False,"2016")
summary2017data = lambda : summaryProducer(False,"2017")
summary2018data = lambda : summaryProducer(False,"2018")