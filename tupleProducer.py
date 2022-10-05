# write a Module for tuple producer
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import *
from selectionFilter import *
from summaryProducer import *
import ROOT as R

ROOT.PyConfig.IgnoreCommandLineOptions = True

import TriggerConfig

class tupleProducer(Module):
    def __init__(self, isMC, era):
        self.isMC = isMC
        self.era = era
        
        # producer hist
        self.producer_hist = R.TH1F('producer_selection','producer_selection',20,0,20)
        self.producer_assigned = 0
        pass
    
    def beginJob(self, histFile=None, histDirName=None):
        pass
    
    def endJob(self):
        pass
    
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.producer_hist.Write()
        pass

    def fill_cut(self, cutnm):
        _ibin = self.producer_hist.GetXaxis().FindBin(cutnm)
        if _ibin == -1:
            self.producer_assigned += 1
            self.producer_hist.GetXaxis().SetBinLabel(self.producer_assigned, cutnm)
            self.producer_hist.SetBinContent( self.producer_assigned, 1 )
        else:
            self.producer_hist.SetBinContent( _ibin, self.producer_hist.GetBinContent(_ibin) + 1 )

    def matchTriggerObject(self):
        pass

    def analyze(self, event):
        # process event, return True (go to next module) or False (fail, go to next event)
        #electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        #jets = Collection(event, "Jet")
        taus = Collection(event, "Tau")
        generator = Object(event, "Generator")
        HLT = Object(event, "HLT")

        signalMuon = None
        has_muon = False
        tag_trig_match = False
        has_tau = True
        btag_veto = True

        if self.isMC:
            pass
        # genEventWeight = ...
        
        # has muon
        if len(muons) >=1:
            has_muon = True
            #signalMuon = selectionFilter(Module).signalMuon
        
        # tag trigger match
        trigFile = './2017trigger.json'
        hltPaths, tagHltPaths = TriggerConfig.LoadAsVPSet(trigFile)

        #print("hltPaths: ", hltPaths)
        tag_trig_match = True

        # has tau
        for _t, _tau in enumerate(taus):
            _tau_v4 = _tau.p4()



        #if self.ctr_evt_processed % 10000 ==0:
        #    print("Processed evt = {}".format(self.ctr_evt_processed))

        self.fill_cut('total')
        if has_muon:
            self.fill_cut('has_muon')
            if tag_trig_match:
                self.fill_cut('tag_trig_match')
                if has_tau:
                    self.fill_cut('has_tau')
                    if btag_veto:
                        self.fill_cut('btag_veto')
                        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

tuple2016MC = lambda : tupleProducer(True,"2016")
tuple2017MC = lambda : tupleProducer(True,"2017")
tuple2018MC = lambda : tupleProducer(True,"2018")
tuple2016data = lambda : tupleProducer(False,"2016")
tuple2017data = lambda : tupleProducer(False,"2017")
tuple2018data = lambda : tupleProducer(False,"2018")