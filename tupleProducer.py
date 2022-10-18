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
        self.count = 0
        
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
                
        self.out.branch("tau_gen_vis_pt", "F")
        self.out.branch("tau_gen_vis_eta", "F")
        self.out.branch("tau_gen_vis_phi", "F")
        self.out.branch("tau_gen_vis_mass", "F")
        
        self.out.branch("tau_gen_rad_pt", "F")
        self.out.branch("tau_gen_rad_eta", "F")
        self.out.branch("tau_gen_rad_phi", "F")
        self.out.branch("tau_gen_rad_energy", "F")
        
        self.out.branch("tau_gen_n_gammas", "I")
        self.out.branch("tau_gen_n_gammas_rad", "I")
        
        self.out.branch("tau_gen_match", "I")
        self.out.branch("tau_gen_charge", "I")
        self.out.branch("tau_sel", "I")
        self.out.branch("vis_mass", "F")
                
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
    
    def findTerminalCopy(self, _gen, gen_part, gen_particles, first):
        part = gen_part
        _gen_num = _gen
        #isFirstCopy = ( ((part.statusFlags) >> 12) & 1 )
        #isLastCopy = ( ((part.statusFlags) >> 13) & 1 )
        while ( ( ((part.statusFlags) >> 13) & 1 ) == 0 ):
            nextCopyFound = False
            daughter_sel = []
            # loop for finding the daughter particle
            for _g, _gpart in enumerate(gen_particles):
                if _gen_num == _gpart.genPartIdxMother:
                    daughter_sel.append( (_g, _gpart) )
            # find the copy particle
            for _index, _dau in enumerate(daughter_sel):
                if gen_part.pdgId == daughter_sel[_index][1].pdgId:
                    part = daughter_sel[_index][1]
                    _gen_num = daughter_sel[_index][0]
                    nextCopyFound = True
                    break
            if nextCopyFound == False:
                raise RuntimeError("Unable to find the last particle!")
        return _gen_num, part

    def findLeptonGenMatch(self, _gen, gen_part, gen_particles):
        result = None
        # return match result
        ele_pdgid = 11
        mu_pdgid = 13
        tau_pdgid = 15
        lepton_pdgid = [11, 13, 15]
        genMatches = dict()
        genMatches[ele_pdgid, False] = [ele_pdgid, False, 1]    # ele
        genMatches[mu_pdgid, False] = [mu_pdgid, False, 2]      # muon
        genMatches[ele_pdgid, True] = [ele_pdgid, True, 3]      # ele from tau
        genMatches[mu_pdgid, True] = [mu_pdgid, True, 4]        # muon from tau
        genMatches[tau_pdgid, False] = [tau_pdgid, False, 5]
        genMatches[tau_pdgid, True] = [tau_pdgid, True, 5]
        pt_thresholds = dict()
        pt_thresholds[ele_pdgid] = [ele_pdgid, 8]
        pt_thresholds[mu_pdgid] = [mu_pdgid, 8]
        pt_thresholds[tau_pdgid] = [tau_pdgid, 15]
        isTauProduct = ( ((gen_part.statusFlags) >> 5) & 1 )
        isPrompt = (gen_part.statusFlags & 1)
        isFirstCopy = ( ((gen_part.statusFlags) >> 12) & 1 )
        #print("statusFlags: {}".format(isTauProduct))
        if not ( ( isPrompt or isTauProduct ) and isFirstCopy ):
            return None
        abs_pdg = abs(gen_part.pdgId)
        if abs_pdg not in pt_thresholds.keys():
            return None
        # find Terminal Copy for finding the last copy particle
        gen_part_lastcopy_idx, gen_part_lastcopy = self.findTerminalCopy(_gen, gen_part, gen_particles, False)
        # find stable particle for lastcopy particle
        decay_sel = [] 
        for _idx, _part in enumerate(gen_particles):
            if _part.status == 1:   # status == 1 stand for stable
                _stable_idx = _idx
                _stable_part = _part
                while _part.genPartIdxMother != -1:
                    if _part.genPartIdxMother == gen_part_lastcopy_idx:
                        decay_sel.append( (_stable_idx, _stable_part) )
                        break
                    _idx = _part.genPartIdxMother
                    _part = gen_particles[_idx]
        if len(decay_sel) >=1:
            for _i, _decay_part in enumerate(decay_sel):
                print("decay sel pdgid: {}".format(_decay_part[1].pdgId))
        # find stable part for particle exclude lastcopy
        rad_sel = []
        for _idx, _part in enumerate(gen_particles):
            if _part.status == 1:   # status == 1 stand for stable
                _stable_idx = _idx
                _stable_part = _part
                while _part.genPartIdxMother != -1:
                    if _part.genPartIdxMother != gen_part_lastcopy_idx:
                        if _part.genPartIdxMother == _gen:
                            rad_sel.append( (_stable_idx, _stable_part) )
                            break
                    else:  # mother_idx == last copy particle, not radiation
                        break
                    _idx = _part.genPartIdxMother
                    _part = gen_particles[_idx]
        # vis and vis_rad fill
        neutrinos_pdg = [12, 14, 16]
        light_lepton = [11, 13]
        gamma = 22
        neu_lep_gamma = [11, 12, 13, 14, 16, 22]
        # fill vis
        light_lepton_sel = []
        decay_gamma_sel = []
        hadron_sel = []
        vis_p4 = R.TLorentzVector()
        for _i, _decay_part in enumerate(decay_sel):
            if abs(_decay_part[1].pdgId) not in neutrinos_pdg:
                vis_p4 += _decay_part[1].p4()
            # pick light lepton if there exist
            if abs(_decay_part[1].pdgId) in light_lepton:
                light_lepton_sel.append(decay_sel[_i])
            # pick gamma if there exist
            if abs(_decay_part[1].pdgId) == gamma:
                decay_gamma_sel.append(decay_sel[_i])
            # pick charged hadron and neutral hadron if there exist
            if abs(_decay_part[1].pdgId) not in neu_lep_gamma:
                hadron_sel.append(decay_sel[_i])
                 
        if len(light_lepton_sel) >=2 :
            print("more than 2 light lepton!")
        # fill vis_rad
        rad_gamma_sel =[]
        vis_rad_p4 = R.TLorentzVector()
        for _i, _rad_part in enumerate(rad_sel):
            if abs(_rad_part[1].pdgId) not in neutrinos_pdg:
                vis_rad_p4 += _rad_part[1].p4()
            # pick rad gamma if there exist
            if abs(_rad_part[1].pdgId) == gamma:
                rad_gamma_sel.append(rad_sel[_i])
        # tot vis p4
        total_vis_p4 = vis_p4 + vis_rad_p4
        # match light lepton
        if ( (abs_pdg == tau_pdgid) and (len(light_lepton_sel) == 1) ):
            abs_lep_pdg = abs(light_lepton_sel[0][1].pdgId)
            pt_thr = pt_thresholds[abs_lep_pdg][1]
            print("pt thr : {}".format(pt_thr))
            if ( (light_lepton_sel[0][1].p4().Pt() > pt_thr) or ( total_vis_p4.Pt() < pt_thr ) ):
                return None
            match = genMatches[abs_lep_pdg, True][2]
            print("light lepton match : {}".format(genMatches[abs_lep_pdg, True][2]))
        else:
            if total_vis_p4.Pt() <= pt_thresholds[abs_pdg][1]:
                return None
            if abs_pdg not in lepton_pdgid:
                print("no match!")
                match = 6  # 6 stand for no match
            else:
                match = genMatches[abs_pdg, isTauProduct][2]
                print("match : {}".format(genMatches[abs_pdg, isTauProduct][2]))
        # fill result
        result = gen_part
        result.match = match
        result.gen_particle_firstCopy = gen_part
        result.gen_particle_lastCopy = gen_part_lastcopy
        result.visible_daughters = decay_sel
        result.visible_radiation = rad_sel
        result.visible_p4 = vis_p4
        result.visible_rad_p4 = vis_rad_p4
        result.n_hadrons = len(hadron_sel)
        result.n_gammas = len(decay_gamma_sel)
        result.n_gammas_rad = len(rad_gamma_sel)
        return result

    def collectGenLeptons(self, gen_particles):
        lt_gen_lepton_sel = []
        for _gen, _gen_part in enumerate(gen_particles):
            result = self.findLeptonGenMatch(_gen, _gen_part, gen_particles)
            if result is not None:
                lt_gen_lepton_sel.append(result)
        return lt_gen_lepton_sel
    
    def lepton_gen_match(self, lep, gen_leptons):
        result = None
        dR2_threshold = 0.2**2
        best_match_dr2 = dR2_threshold
        for _g, _genlep in enumerate(gen_leptons):
            total_vis_p4 = _genlep.visible_p4 + _genlep.visible_rad_p4
            dr2_vis = self.deltaR2(lep.p4(), _genlep.visible_p4)
            dr2_tot_vis = self.deltaR2(lep.p4(), total_vis_p4)
            dr2 = min(dr2_vis, dr2_tot_vis)
            if dr2 >= best_match_dr2:
                continue
            best_match_dr2 = dr2
            result = _genlep
        return result

    def collectTaus(self, signalMu_v4, taus, genleptons, deltaR2Thr):
        gen = 1
        pt = 2
        mva = 4
        deepTau = 8
        best_tau = dict()
        for _t, _tau in enumerate(taus):
            _tau_v4 = _tau.p4()
            if (_tau_v4.Pt()) > 18 and (abs(_tau_v4.Eta()) < 2.3) and (self.deltaR2(signalMu_v4, _tau_v4) > deltaR2Thr):
                pass_mva_sel = (_tau.rawMVAoldDM2017v2 > 0)    # and tau.tauID("againstMuonLoose3") > 0.5f ?
                pass_deep_sel = ( (_tau.rawDeepTau2017v2p1VSjet > 0) and (_tau.idDeepTau2017v2p1VSe > 0.5) and (_tau.idDeepTau2017v2p1VSmu > 0.5) )
                if (pass_mva_sel or pass_deep_sel) and ( (pt not in best_tau.items()) or (best_tau[pt].p4().Pt() < _tau.p4().Pt()) ):
                    best_tau[pt] = _tau
                if pass_mva_sel and ( (mva not in best_tau.items()) or (best_tau[mva].rawMVAoldDM2017v2 < _tau.rawMVAoldDM2017v2) ):
                    best_tau[mva] = _tau
                if pass_deep_sel and ( (deepTau not in best_tau.items()) or (best_tau[deepTau].rawDeepTau2017v2p1VSjet < _tau.rawDeepTau2017v2p1VSjet) ):
                    best_tau[deepTau] = _tau
        #
        selected_gen_tau = self.selectGenLeg(genleptons, True)
        if selected_gen_tau is None:
            has_selected_gen_tau = False
        else:
            has_selected_gen_tau = (selected_gen_tau.match != 6)
        # 
        selected_gen_tau_stored = False
        # dict loop
        selected_taus_selection = 0
        selected_taus = []
        for _i, _item in enumerate(best_tau.items()):
            reco_tau_id = _item[0]  # _item[0] means pt, mva, deepTau
            reco_tau = _item[1]     # _item[1] means best tau 's values
            if (selected_taus is None) or (reco_tau not in (selected_taus[_x][0] for _x in range( len(selected_taus) )) ):
                gen_tau = self.lepton_gen_match(reco_tau, genleptons)
                if gen_tau is None:
                    has_gen_tau = False
                else:
                    has_gen_tau = (gen_tau.match != 6)
                selected_taus.append( (reco_tau, gen_tau) )
                if ( has_selected_gen_tau and has_gen_tau and selected_gen_tau.gen_particle_firstCopy == gen_tau.gen_particle_firstCopy ):
                    selected_gen_tau_stored = True
                    selected_taus_selection |= gen
            selected_taus_selection |= reco_tau_id
        if (has_selected_gen_tau) and (not selected_gen_tau_stored):
            reco_tau = None
            for _t, _tau in enumerate(taus):
                gen_tau = self.lepton_gen_match(_tau, genleptons)
                if (gen_tau is not None) and (gen_tau.match !=6) and (gen_tau.gen_particle_firstCopy == selected_gen_tau.gen_particle_firstCopy):
                    reco_tau = _tau
                    break
            if (reco_tau is not None) and (selected_taus is not None) and (reco_tau in (selected_taus[_x][0] for _x in range( len(selected_taus) ))):
                raise RuntimeError("Inconsistency in CollectTaus algorithm.")
            selected_taus.append( (reco_tau, selected_gen_tau) )
            selected_taus_selection = gen
        return selected_taus, selected_taus_selection

    def selectGenLeg(self, genleptons, True):
        matches = 5  # 5 stand for tau
        leg = None
        for _l, lepton in enumerate(genleptons):
            if (lepton.match == 5) and ( (leg is None) or (leg.visible_p4.Pt() < lepton.visible_p4.Pt()) ):
                leg = lepton
        if leg is not None:
            print("leg visible pdgid,pt :{0}, {1}".format(leg.pdgId, leg.visible_p4.Pt()))
        else:
            print("leg is None, match==6")
        return leg

    def analyze(self, event):
        # process event, return True (go to next module) or False (fail, go to next event)
        #electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        jets = Collection(event, "Jet")
        met = Object(event, "MET")
        taus = Collection(event, "Tau")
        genparts = Collection(event, "GenPart")     # gen particles ?
        gendressedlep = Collection(event, "GenDressedLepton")    # gen leptons ?
        gen_vis_taus = Collection(event, "GenVisTau")
        generator = Object(event, "Generator")
        HLT = Object(event, "HLT")

        signalMu = None
        signalMu_v4 = None
        has_muon = False
        tag_trig_match = False
        has_tau = True
        btag_veto = False
        deltaR2Thr = 0.5*0.5
        sele_tau = None        

        if self.isMC:
            pass
            #genWeight = ...
                        
        # collect gen leptons
        genleptons = self.collectGenLeptons(genparts)
        self.count += 1
        print("**********")
        #if self.count == 300:
        #    exit(0)
        
        # has muon
        signal_muon_sel = []
        sele_filt = selectionFilter.selectionFilter(True,"2018")
        signal_muon_sel = sele_filt.find_signal_muon(muons)
        signalMu_v4 = signal_muon_sel[0][0]     # muon_ref_p4
        signalMu = signal_muon_sel[0][1]
        if len(signal_muon_sel) >=1:
            gen_muon = self.lepton_gen_match(signalMu, genleptons)
            if gen_muon is not None:
                print("gen_muon_match:{}".format(gen_muon.match))
            else:
                print("gen muon is None")
            has_muon = True
            #print("signal muon : {}",signalMuon)
        
        # tag trigger match
        # TODO
        trigFile = './2017trigger.json'
        hltPaths, tagHltPaths = TriggerConfig.LoadAsVPSet(trigFile)

        #print("hltPaths: ", hltPaths)
        tag_trig_match = True

        # sele tau as a list? (reco_tau, gen_tau, reco_tau_id)
        selected_taus, sel = self.collectTaus(signalMu_v4, taus, genleptons, deltaR2Thr)
        if len(selected_taus) <= 0:
            has_tau = False
            tau = None
            gen_tau = None
        else:
            # pick selected_taus[-1]
            tau = selected_taus[-1][0]       # reco tau
            gen_tau = selected_taus[-1][1]   # gen_tau
        if tau is None:
            has_reco_tau = False
        else:
            has_reco_tau = True
        if gen_tau is None:
            has_gen_tau = False
        else:
            has_gen_tau = (gen_tau.match != 6)
        # btag veto
        tau_ref_p4 = R.TLorentzVector()
        if has_reco_tau:
            tau_ref_p4 = tau.p4()
        elif has_gen_tau:
            tau_ref_p4 = gen_tau.visible_p4
        btagThreshold = 0.9
        if btagThreshold > 0:
            for _jet in jets:
                _j_v4 = _jet.p4()
                if (self.deltaR2(signalMu_v4, _j_v4) > 0.5*0.5) and (self.deltaR2(tau_ref_p4, _j_v4) > 0.5*0.5) and (_j_v4.Pt() > 20) and (abs(_j_v4.Eta()) < 2.4) and (_jet.btagDeepFlavB > btagThreshold):
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
                        #if (gen_muon is not None) and (gen_muon.match != 6):
                        if (gen_muon is not None):
                            self.out.fillBranch("muon_gen_match", gen_muon.match)
                            self.out.fillBranch("muon_gen_pt",gen_muon.p4().Pt())
                            self.out.fillBranch("muon_gen_eta",gen_muon.p4().Eta())
                            self.out.fillBranch("muon_gen_phi",gen_muon.p4().Phi())
                            self.out.fillBranch("muon_gen_mass",gen_muon.p4().M())
                            self.out.fillBranch("muon_gen_charge",(-gen_muon.pdgId) / abs(gen_muon.pdgId) )
                        else:
                            self.out.fillBranch("muon_gen_match", 6)
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
                        if has_gen_tau:
                            self.out.fillBranch("tau_gen_vis_pt", gen_tau.visible_p4.Pt())
                            self.out.fillBranch("tau_gen_vis_eta", gen_tau.visible_p4.Eta())
                            self.out.fillBranch("tau_gen_vis_phi", gen_tau.visible_p4.Phi())
                            self.out.fillBranch("tau_gen_vis_mass", gen_tau.visible_p4.M())
                            
                            self.out.fillBranch("tau_gen_rad_pt", gen_tau.visible_rad_p4.Pt())
                            self.out.fillBranch("tau_gen_rad_eta", gen_tau.visible_rad_p4.Eta())
                            self.out.fillBranch("tau_gen_rad_phi", gen_tau.visible_rad_p4.Phi())
                            self.out.fillBranch("tau_gen_rad_energy", gen_tau.visible_rad_p4.E())
                            
                            self.out.fillBranch("tau_gen_n_gammas", gen_tau.n_gammas)
                            self.out.fillBranch("tau_gen_n_gammas_rad", gen_tau.n_gammas_rad)
                            self.out.fillBranch("tau_gen_charge", (-gen_tau.pdgId) / abs(gen_tau.pdgId) )
                            self.out.fillBranch("tau_gen_match", gen_tau.match)
                        else:
                            self.out.fillBranch("tau_gen_vis_pt", -999.0)
                            self.out.fillBranch("tau_gen_vis_eta", -999.0)
                            self.out.fillBranch("tau_gen_vis_phi", -999.0)
                            self.out.fillBranch("tau_gen_vis_mass", -999.0)
                            
                            self.out.fillBranch("tau_gen_rad_pt", -999.0)
                            self.out.fillBranch("tau_gen_rad_eta", -999.0)
                            self.out.fillBranch("tau_gen_rad_phi", -999.0)
                            self.out.fillBranch("tau_gen_rad_energy", -999.0)
                            
                            self.out.fillBranch("tau_gen_n_gammas", -999)
                            self.out.fillBranch("tau_gen_n_gammas_rad", -999)
                            self.out.fillBranch("tau_gen_charge", -999)
                            self.out.fillBranch("tau_gen_match", 6)
                        # fill tau
                        self.out.fillBranch("tau_sel", sel)
                        if has_reco_tau:
                            self.out.fillBranch("tau_pt", tau.p4().Pt())
                            self.out.fillBranch("tau_eta", tau.p4().Eta())
                            self.out.fillBranch("tau_phi", tau.p4().Phi())
                            self.out.fillBranch("tau_mass", tau.p4().M())
                            self.out.fillBranch("tau_charge", tau.charge)
                            self.out.fillBranch("tau_decayMode", tau.decayMode)
                            self.out.fillBranch("tau_dxy", tau.dxy)
                            self.out.fillBranch("tau_dz", tau.dz)
                            self.out.fillBranch("tau_idDeepTau2017v2p1VSe", tau.idDeepTau2017v2p1VSe)
                            self.out.fillBranch("tau_idDeepTau2017v2p1VSmu", tau.idDeepTau2017v2p1VSmu)
                            self.out.fillBranch("tau_idDeepTau2017v2p1VSjet", tau.idDeepTau2017v2p1VSjet)
                            self.out.fillBranch("tau_rawDeepTau2017v2p1VSe", tau.rawDeepTau2017v2p1VSe)
                            self.out.fillBranch("tau_rawDeepTau2017v2p1VSmu", tau.rawDeepTau2017v2p1VSmu)
                            self.out.fillBranch("tau_rawDeepTau2017v2p1VSjet", tau.rawDeepTau2017v2p1VSjet)
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
                        # fill visible mass
                        self.out.fillBranch("vis_mass", (signalMu_v4 + tau_ref_p4).M())
                        return True
        return False

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

tuple2016MC = lambda : tupleProducer(True,"2016")
tuple2017MC = lambda : tupleProducer(True,"2017")
tuple2018MC = lambda : tupleProducer(True,"2018")
tuple2016data = lambda : tupleProducer(False,"2016")
tuple2017data = lambda : tupleProducer(False,"2017")
tuple2018data = lambda : tupleProducer(False,"2018")

