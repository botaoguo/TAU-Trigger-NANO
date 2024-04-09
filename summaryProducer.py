# write a Module for summary producer
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.output import OutputTree
from PhysicsTools.NanoAODTools.postprocessing.tools import *
import ROOT as R
from array import array
import triggerDescriptor as trig
ROOT.PyConfig.IgnoreCommandLineOptions = True

class summaryProducer(Module):
    def __init__(self, isMC, era):
        #self.writeHistFile = True

        self.isMC = isMC
        self.era = era
        self.ctr_evt_processed = 0
        #self.tot_evt_genweight = 0
        
        self.summary = R.TTree("summary", "summary")
        pass
    
    def beginJob(self, histFile=None, histDirName=None):
        # trigFile = './{}trigger.json'.format(self.era)
        # hltPaths, tagHltPaths = trig.LoadAsList(trigFile)
        #filtName_sel = self.getFilterName('./2018filterName.txt', False)
        #filtHash_sel = self.getFilterName('./2018filterHash.txt', True) 

        # index = R.std.vector('unsigned int')()
        # pattern = R.std.vector('string')()
        #filt_hash = R.std.vector('unsigned int')()
        #filt_name = R.std.vector('string')()
        
        #for _n in range( len(filtName_sel) ):
        #    filt_name.push_back( filtName_sel[_n] )
        #for _h in range( len(filtHash_sel) ):
            #print("filter_Hash : ", filtHash_sel[_h])
            #filt_hash.push_back( int(filtHash_sel[_h]) )
        # for _i, _hltp in enumerate(hltPaths):
        #     index.push_back(_i)
        #     pattern.push_back(_hltp[3])
        
        #self.summary.Branch("filter_name", filt_name)
        #self.summary.Branch("filter_hash", filt_hash)
        # self.summary.Branch("trigger_index", index)
        # self.summary.Branch("trigger_pattern", pattern)
        # self.summary.Fill()
        pass
    
    def endJob(self):
        #self.summary.Fill()
        pass
    
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.tot_evt = inputTree.GetEntries()
        print("total event : {}".format(self.tot_evt))
        self.out = wrappedOutputTree

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.summary.Write()
        pass

    def fill_cut(self, cutnm):
        _ibin = self.cutflow_hist.GetXaxis().FindBin(cutnm)
        if _ibin == -1:
            self.cutflow_assigned += 1
            self.cutflow_hist.GetXaxis().SetBinLabel(self.cutflow_assigned, cutnm)
            self.cutflow_hist.SetBinContent( self.cutflow_assigned, 1 )
        else:
            self.cutflow_hist.SetBinContent( _ibin, self.cutflow_hist.GetBinContent(_ibin) + 1 )

    def getFilterName(self, filt_nm_file, ishash):
        _filt_sel = []
        if ishash:
            for line in open(filt_nm_file, 'r'):
                rs = line.rstrip('\n')
                _filt_sel.append(rs)
        else:
            for line in open(filt_nm_file, 'r'):
                rs = line.rstrip('\n')
                _filt_sel.append(rs)
        return _filt_sel

    def analyze(self, event):
        # process event, return True (go to next module) or False (fail, go to next event)
        
        #electrons = Collection(event, "Electron")
        #muons = Collection(event, "Muon")
        #jets = Collection(event, "Jet")
        #taus = Collection(event, "Tau")
        generator = Object(event, "Generator")

        self.ctr_evt_processed += 1
        #self.tot_evt_genweight += generator.weight
        #self.fill_cut('tot_evt')

        if self.ctr_evt_processed % 1000 ==0:
            print("Processed evt = {}".format(self.ctr_evt_processed))
            #exit(0)

        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

summary2016MC = lambda : summaryProducer(True,"2016")
summary2017MC = lambda : summaryProducer(True,"2017")
summary2018MC = lambda : summaryProducer(True,"2018")
summary2022MC = lambda : summaryProducer(True,"2022")
summary2023MC = lambda : summaryProducer(True,"2023")
summary2016data = lambda : summaryProducer(False,"2016")
summary2017data = lambda : summaryProducer(False,"2017")
summary2018data = lambda : summaryProducer(False,"2018")
summary2022data = lambda : summaryProducer(False,"2022")
summary2023data = lambda : summaryProducer(False,"2023")