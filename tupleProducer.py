# write a Module for tuple producer
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import *
import selectionFilter
from summaryProducer import *
import ROOT as R
import math

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
        self.out.branch("muon_pt", "F")
        self.out.branch("muon_eta", "F")
        self.out.branch("muon_phi", "F")
        self.out.branch("muon_mass", "F")
        self.out.branch("muon_iso", "F")
        self.out.branch("muon_charge", "I")
        self.out.branch("muon_mt", "F")
        
        self.out.branch("muon_gen_pt", "F")
        self.out.branch("muon_gen_eta", "F")
        self.out.branch("muon_gen_phi", "F")
        self.out.branch("muon_gen_mass", "F")
        self.out.branch("muon_gen_match", "I")
        self.out.branch("muon_gen_charge", "I")
        
        self.out.branch("tau_pt", "F")
        self.out.branch("tau_eta", "F")
        self.out.branch("tau_phi", "F")
        self.out.branch("tau_mass", "F")
        self.out.branch("tau_charge", "I")
        
        self.out.branch("tau_decayMode", "I")
        self.out.branch("tau_dxy", "F")
        self.out.branch("tau_dz", "F")
        
        self.out.branch("tau_idDeepTau2017v2p1VSe", "I")
        self.out.branch("tau_idDeepTau2017v2p1VSmu", "I")
        self.out.branch("tau_idDeepTau2017v2p1VSjet", "I")

        self.out.branch("tau_rawDeepTau2017v2p1VSe", "F")
        self.out.branch("tau_rawDeepTau2017v2p1VSmu", "F")
        self.out.branch("tau_rawDeepTau2017v2p1VSjet", "F")
                
        self.out.branch("tau_gen_pt", "F")
        self.out.branch("tau_gen_eta", "F")
        self.out.branch("tau_gen_phi", "F")
        self.out.branch("tau_gen_mass", "F")
        self.out.branch("tau_gen_match", "I")
        self.out.branch("tau_gen_charge", "I")
                
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

    def deltaR2(self, l1_v4, l2_v4):
        dphi = l1_v4.Phi() - l2_v4.Phi()
        while dphi<-math.pi:
            dphi+=2*math.pi
        while dphi>math.pi:
            dphi-=2*math.pi        
        deta = l1_v4.Eta() - l2_v4.Eta()
        return (dphi**2 + deta**2)

    def lepton_gen_match(self, lep, genLeptons):
        result = None
        ele_pdgid = 11
        mu_pdgid = 13
        tau_pdgid = 15
        genMatches = dict()
        genMatches['ele'] = [ele_pdgid, False, 1]
        genMatches['tauele'] = [ele_pdgid, True, 2]
        genMatches['mu'] = [mu_pdgid, False, 3]
        genMatches['taumu'] = [mu_pdgid, True, 4]
        genMatches['tau'] = [tau_pdgid, 5]
        dR2_threshold = 0.2**2
        best_match_dr2 = dR2_threshold
        for _g, _genlep in enumerate(genLeptons):
            #print("genlep pdgid : ",_genlep.pdgId)
            dr2 = self.deltaR2(lep.p4(), _genlep.p4())
            if dr2 > best_match_dr2:
                continue
            best_match_dr2 = dr2
            result = _genlep
        #if ( (genLeptons.hasTauAnc) )
        if result is not None:
            result.match = 6 # 6 stand no match
            if ( (abs(_genlep.pdgId) == ele_pdgid) and (_genlep.p4().Pt() > 8) ):
                if (_genlep.hasTauAnc):
                    result.match = genMatches['tauele'][2]
                else:
                    result.match = genMatches['ele'][2]
            if ( (abs(_genlep.pdgId) == mu_pdgid) and (_genlep.p4().Pt() > 8) ):
                if (_genlep.hasTauAnc):
                    result.match = genMatches['taumu'][2]
                else:
                    result.match = genMatches['mu'][2]
            if ( (abs(_genlep.pdgId) == tau_pdgid) and (_genlep.p4().Pt() > 15) ):
                result.match = genMatches['tau'][1]
            #print("result match : {}".format(result.match))
        return result                

    def selectGenLeg(self, gen_vistaus):
        result = None
        vis_ptmax = 0
        for _g, _gen_vis_tau in enumerate(gen_vistaus):
            if _gen_vis_tau.p4().Pt() is not None:
                if _gen_vis_tau.p4().Pt() > vis_ptmax:
                    result = _gen_vis_tau
                    vis_ptmax = _gen_vis_tau.p4().Pt()
                result.match = 5
            else:
                result.match = 6
        return result
            

    def analyze(self, event):
        # process event, return True (go to next module) or False (fail, go to next event)
        #electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        jets = Collection(event, "Jet")
        met = Object(event, "MET")
        taus = Collection(event, "Tau")
        genparts = Collection(event, "GenPart")     # gen particles ?
        genleptons = Collection(event, "GenDressedLepton")    # gen leptons ?
        gen_vis_taus = Collection(event, "GenVisTau")
        generator = Object(event, "Generator")
        HLT = Object(event, "HLT")

        signalMu = None
        signalMu_v4 = None
        has_muon = False
        tag_trig_match = False
        has_tau = True
        btag_veto = False
        best_tau = dict()
        sele_tau = None        

        if self.isMC:
            pass
            #genWeight = ...
        
        # has muon
        signal_muon_sel = []
        sele_filt = selectionFilter.selectionFilter(True,"2017")
        signal_muon_sel = sele_filt.find_signal_muon(muons)
        signalMu_v4 = signal_muon_sel[0][0]     # muon_ref_p4
        signalMu = signal_muon_sel[0][1]
        if len(signal_muon_sel) >=1:
            gen_muon = self.lepton_gen_match(signalMu, genleptons)
            has_muon = True
            #print("signal muon : {}",signalMuon)
        
        # tag trigger match
        trigFile = './2017trigger.json'
        hltPaths, tagHltPaths = TriggerConfig.LoadAsVPSet(trigFile)

        #print("hltPaths: ", hltPaths)
        tag_trig_match = True

        # has tau
        # TODO
        for _t, _tau in enumerate(taus):
            _tau_v4 = _tau.p4()
            if (_tau_v4.Pt()) > 18 and (abs(_tau_v4.Eta()) < 2.3) and (self.deltaR2(signalMu_v4, _tau_v4) > 0.5*0.5):
                #pass_mva_sel = (_tau.rawMVAoldDM2017v2 > 0)
                pass_mva_sel = (_tau.idAntiMu > 0.5)
                pass_deep_sel = ( (_tau.rawDeepTau2017v2p1VSjet > 0) and (_tau.idDeepTau2017v2p1VSe > 0.5) and (_tau.idDeepTau2017v2p1VSmu > 0.5) )
                if (pass_mva_sel or pass_deep_sel) and ( ('pt' not in best_tau.items()) or (best_tau['pt'].p4().Pt() < _tau.p4().Pt()) ):
                    best_tau['pt'] = _tau
                if pass_mva_sel and ( ('mva' not in best_tau.items()) or (best_tau['mva'].idAntiMu < _tau.idAntiMu) ):
                    best_tau['mva'] = _tau
                if pass_deep_sel and ( ('deepTau' not in best_tau.items()) or (best_tau['deepTau'].rawDeepTau2017v2p1VSjet < _tau.rawDeepTau2017v2p1VSjet) ):
                    best_tau['deepTau'] = _tau
        #sele_gen_tau = None
        sele_gen_tau = self.selectGenLeg(gen_vis_taus)
        #has_selected_gen_tau = False
        #has_selected_gen_tau = (sele_gen_tau.match != 6)
        #selected_gen_tau_stored = False
        #for _b, _best_tau in enumerate(best_tau):
        #    if _best_tau not in taus:
        #        gen_tau = None
        #        gen_tau = self.lepton_gen_match(_best_tau, genleptons)
        #        has_gen_tau = False
        #        has_gen_tau = (gen_tau.match != 6)
        
        #if best_tau is not None:
        if 'pt' in best_tau:
            sele_tau = best_tau['pt']
        elif 'mva' in best_tau:
            sele_tau = best_tau['mva']
        elif 'deepTau' in best_tau:
            sele_tau = best_tau['deepTau']
        else:
            has_tau = False

        # btag veto
        if has_tau:
            btagThreshold = 0.9
            if btagThreshold > 0:
                for _jet in jets:
                    _j_v4 = _jet.p4()
                    if (self.deltaR2(signalMu_v4, _j_v4) > 0.5*0.5) and (self.deltaR2(sele_tau.p4(), _j_v4) > 0.5*0.5) and (_j_v4.Pt() > 20) and (abs(_j_v4.Eta()) < 2.4) and (_jet.btagDeepFlavB > btagThreshold):
                        btag_veto = True
                        break

        #if self.ctr_evt_processed % 10000 ==0:
        #    print("Processed evt = {}".format(self.ctr_evt_processed))

        self.fill_cut('total')
        if has_muon:
            self.fill_cut('has_muon')
            if tag_trig_match:
                self.fill_cut('tag_trig_match')
                if has_tau:
                    self.fill_cut('has_tau')
                    if not btag_veto:
                        self.fill_cut('btag_veto')
                        # fill muon and gen muon
                        if gen_muon is not None:
                            has_gen_muon = False
                            has_gen_muon = (gen_muon.match != 6)
                            self.out.fillBranch("muon_gen_match", gen_muon.match)
                            if has_gen_muon:
                                self.out.fillBranch("muon_gen_pt",gen_muon.p4().Pt())
                                self.out.fillBranch("muon_gen_eta",gen_muon.p4().Eta())
                                self.out.fillBranch("muon_gen_phi",gen_muon.p4().Phi())
                                self.out.fillBranch("muon_gen_mass",gen_muon.p4().M())
                                self.out.fillBranch("muon_gen_charge",(-gen_muon.pdgId) / abs(gen_muon.pdgId) )
                            else:
                                self.out.fillBranch("muon_gen_pt", -999.0)
                                self.out.fillBranch("muon_gen_eta", -999.0)
                                self.out.fillBranch("muon_gen_phi", -999.0)
                                self.out.fillBranch("muon_gen_mass", -999.0)
                                self.out.fillBranch("muon_gen_charge", -999)
                        if len(signal_muon_sel) >=1:
                            self.out.fillBranch("muon_pt",signalMu_v4.Pt())
                            self.out.fillBranch("muon_eta",signalMu_v4.Eta())
                            self.out.fillBranch("muon_phi",signalMu_v4.Phi())
                            self.out.fillBranch("muon_mass",signalMu_v4.M())
                            self.out.fillBranch("muon_iso",signalMu.pfRelIso04_all)
                            self.out.fillBranch("muon_charge",signalMu.charge)
                            self.out.fillBranch("muon_mt",sele_filt.calc_MT(signalMu_v4,sele_filt.calc_met(met)))
                        else:
                            self.out.fillBranch("muon_pt", -999.0)
                            self.out.fillBranch("muon_eta", -999.0)
                            self.out.fillBranch("muon_phi", -999.0)
                            self.out.fillBranch("muon_mass", -999.0)
                            self.out.fillBranch("muon_iso", -999.0)
                            self.out.fillBranch("muon_charge", -999)
                            self.out.fillBranch("muon_mt", -999.0)
                        # fill gen tau
                        has_gen_tau = False
                        has_gen_tau = ( (sele_gen_tau is not None) and (sele_gen_tau.match != 6) )
                        if has_gen_tau:
                            self.out.fillBranch("tau_gen_pt", sele_gen_tau.p4().Pt())
                            self.out.fillBranch("tau_gen_eta", sele_gen_tau.p4().Eta())
                            self.out.fillBranch("tau_gen_phi", sele_gen_tau.p4().Phi())
                            self.out.fillBranch("tau_gen_mass", sele_gen_tau.p4().M())
                            self.out.fillBranch("tau_gen_charge", sele_gen_tau.charge)
                            self.out.fillBranch("tau_gen_match", sele_gen_tau.match)
                        else:
                            self.out.fillBranch("tau_gen_pt", -999.0)
                            self.out.fillBranch("tau_gen_eta", -999.0)
                            self.out.fillBranch("tau_gen_phi", -999.0)
                            self.out.fillBranch("tau_gen_mass", -999.0)
                            self.out.fillBranch("tau_gen_charge", -999)
                            self.out.fillBranch("tau_gen_match", 6)
                        # fill tau
                        if sele_tau is not None:
                            self.out.fillBranch("tau_pt", sele_tau.p4().Pt())
                            self.out.fillBranch("tau_eta", sele_tau.p4().Eta())
                            self.out.fillBranch("tau_phi", sele_tau.p4().Phi())
                            self.out.fillBranch("tau_mass", sele_tau.p4().M())
                            self.out.fillBranch("tau_charge", sele_tau.charge)
                            self.out.fillBranch("tau_decayMode", sele_tau.decayMode)
                            self.out.fillBranch("tau_dxy", sele_tau.dxy)
                            self.out.fillBranch("tau_dz", sele_tau.dz)
                            self.out.fillBranch("tau_idDeepTau2017v2p1VSe", sele_tau.idDeepTau2017v2p1VSe)
                            self.out.fillBranch("tau_idDeepTau2017v2p1VSmu", sele_tau.idDeepTau2017v2p1VSmu)
                            self.out.fillBranch("tau_idDeepTau2017v2p1VSjet", sele_tau.idDeepTau2017v2p1VSjet)
                            self.out.fillBranch("tau_rawDeepTau2017v2p1VSe", sele_tau.rawDeepTau2017v2p1VSe)
                            self.out.fillBranch("tau_rawDeepTau2017v2p1VSmu", sele_tau.rawDeepTau2017v2p1VSmu)
                            self.out.fillBranch("tau_rawDeepTau2017v2p1VSjet", sele_tau.rawDeepTau2017v2p1VSjet)
                        else:
                            self.out.fillBranch("tau_pt", -999.0)
                            self.out.fillBranch("tau_eta", -999.0)
                            self.out.fillBranch("tau_phi", -999.0)
                            self.out.fillBranch("tau_mass", -999.0)
                            self.out.fillBranch("tau_charge", -999)
                            self.out.fillBranch("tau_decayMode", -999)
                            self.out.fillBranch("tau_dxy", -999.0)
                            self.out.fillBranch("tau_dz", -999.0)
                            self.out.fillBranch("tau_idDeepTau2017v2p1VSe", -999)
                            self.out.fillBranch("tau_idDeepTau2017v2p1VSmu", -999)
                            self.out.fillBranch("tau_idDeepTau2017v2p1VSjet", -999)
                            self.out.fillBranch("tau_rawDeepTau2017v2p1VSe", -999.0)
                            self.out.fillBranch("tau_rawDeepTau2017v2p1VSmu", -999.0)
                            self.out.fillBranch("tau_rawDeepTau2017v2p1VSjet", -999.0)
                        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

tuple2016MC = lambda : tupleProducer(True,"2016")
tuple2017MC = lambda : tupleProducer(True,"2017")
tuple2018MC = lambda : tupleProducer(True,"2018")
tuple2016data = lambda : tupleProducer(False,"2016")
tuple2017data = lambda : tupleProducer(False,"2017")
tuple2018data = lambda : tupleProducer(False,"2018")

