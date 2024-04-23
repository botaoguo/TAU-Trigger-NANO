import os, sys
import ROOT
import argparse

parser = argparse.ArgumentParser(description='Skim full tuple.')
parser.add_argument('--input', required=True, type=str, nargs='+', help="input files")
# parser.add_argument('--type', required=True, type=str, help="pnet or deeptau")
args = parser.parse_args()


path_prefix = '' if 'TAU-Trigger-NANO' in os.getcwd() else 'TAU-Trigger-NANO/'
sys.path.insert(0, path_prefix + 'Common/python')
from AnalysisTools import *

ROOT.gInterpreter.Declare('#include "{}interface/tau_ntupler.h"'.format(path_prefix))

input_vec = ListToStdVector(args.input)
df = ROOT.RDataFrame('Events', input_vec)

# Tag and Probe
# tag the signal muon
df = df.Define("sig_muon_idx","SelectMuon(Muon_pt, Muon_eta, Muon_phi, Muon_pfRelIso04_all, Muon_mediumId)")\
       .Filter("sig_muon_idx != -1")
df = df.Filter("HLT_IsoMu24==1")
df = (
    df.Define("sig_muon_pt", "Muon_pt[sig_muon_idx]")
    .Define("sig_muon_eta", "Muon_eta[sig_muon_idx]")
    .Define("sig_muon_phi", "Muon_phi[sig_muon_idx]")
    .Define("sig_muon_mass", "Muon_mass[sig_muon_idx]")
    .Define("sig_muon_iso", "Muon_pfRelIso04_all[sig_muon_idx]")
    .Define("sig_muon_mediumId", "Muon_mediumId[sig_muon_idx]")
)
# match the trigobj and muon with dR < 0.5 and objId==13
df = df.Define("match_sig_muon","Muon_match(TrigObj_id, TrigObj_eta, TrigObj_phi, sig_muon_eta, sig_muon_phi)")\
       .Filter("match_sig_muon == 1")

# probe tau
# tau_pt > 20, tau_eta < 2.3, VSe >= 2, VSmu >= 1, VSjet >= 5
df = df.Define("probe_tau_idx","SelectTau(Tau_pt, Tau_eta, Tau_idDeepTau2018v2p5VSe, Tau_idDeepTau2018v2p5VSmu, Tau_idDeepTau2018v2p5VSjet)")\
       .Filter("probe_tau_idx != -1")

# Find the second leading tau
# df = df.Define("subleading_tau_idx", "FindSubLeadingTau(Tau_pt, probe_tau_idx)")


# Create new columns with the kinematics of good leptons
df = (
    df.Define("leading_tau_pt", "Tau_pt[probe_tau_idx]")
    .Define("leading_tau_eta", "Tau_eta[probe_tau_idx]")
    .Define("leading_tau_phi", "Tau_phi[probe_tau_idx]")
    .Define("leading_tau_mass", "Tau_mass[probe_tau_idx]")
    .Define("leading_tau_decayMode", "Tau_decayMode[probe_tau_idx]")
    .Define("leading_tau_idDeepTauVSe", "Tau_idDeepTau2018v2p5VSe[probe_tau_idx]")
    .Define("leading_tau_idDeepTauVSmu", "Tau_idDeepTau2018v2p5VSmu[probe_tau_idx]")
    .Define("leading_tau_idDeepTauVSjet", "Tau_idDeepTau2018v2p5VSjet[probe_tau_idx]")
)
# match tau with trig dR < 0.5
df = df.Define("match_probe_tau","Tau_match(TrigObj_id, TrigObj_eta, TrigObj_phi, leading_tau_eta, leading_tau_phi)")\
       .Filter("match_probe_tau == 1")
df = df.Filter("leading_tau_decayMode != 5 && leading_tau_decayMode != 6")

# electron veto
df = df.Define("veto_ele", "VetoEle(Electron_pt, Electron_eta, Electron_mvaIso)")\
       .Filter("veto_ele == 1")

# sig_muon and probe_tau dR > 0.5
df = df.Define("muon_tau_dR", "deltaR(sig_muon_eta, leading_tau_eta, sig_muon_phi, leading_tau_phi)")\
       .Filter("muon_tau_dR > 0.5")
# mt betweeb muon and met should < 30 GeV
df = df.Define("muon_mt", "MT(sig_muon_pt, PuppiMET_pt, sig_muon_phi, PuppiMET_phi)")\
       .Filter("muon_mt > 30")
# vis_mass in [40,80] GeV
df = df.Define("vis_mass", "VisMass(sig_muon_pt, sig_muon_eta, sig_muon_phi, sig_muon_mass, leading_tau_pt, leading_tau_eta, leading_tau_phi, leading_tau_mass)")\
       .Filter("vis_mass >= 40 && vis_mass <= 80")

# df = (
#     df.Define("subleading_tau_pt", "Tau_pt[subleading_tau_idx]")
#     .Define("subleading_tau_eta", "Tau_eta[subleading_tau_idx]")
#     .Define("subleading_tau_phi", "Tau_phi[subleading_tau_idx]")
#     .Define("subleading_tau_mass", "Tau_mass[subleading_tau_idx]")
#     .Define("subleading_tau_decayMode", "Tau_decayMode[subleading_tau_idx]")
#     .Define("subleading_tau_idDeepTauVSe", "Tau_idDeepTau2018v2p5VSe[subleading_tau_idx]")
#     .Define("subleading_tau_idDeepTauVSmu", "Tau_idDeepTau2018v2p5VSmu[subleading_tau_idx]")
#     .Define("subleading_tau_idDeepTauVSjet", "Tau_idDeepTau2018v2p5VSjet[subleading_tau_idx]")
# )

# df = df.Define("has_vbfjets","SelectVBFJets(Jet_pt, Jet_eta, Jet_phi, Jet_mass, leading_tau_eta, leading_tau_phi)")
# df = df.Define("has_jet","SelectJet(Jet_pt, Jet_eta, Jet_phi, Jet_mass, leading_tau_eta, leading_tau_phi)")

df = df.Define("weight","1")

df_pnet = df.Define("pass_ditau","PassDiTauPNet(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi)")
# df_pnet = df_pnet.Define("pass_singletau","HLT_SinglePNetTauhPFJet130_Loose_L2NN_eta2p3==1")
# df_pnet = df_pnet.Define("pass_ditaujet","has_jet==1 && HLT_DoublePNetTauhPFJet26_L2NN_eta2p3_PFJet60==1")
# df_pnet = df_pnet.Define("pass_vbfsingletau","has_vbfjets==1 && HLT_VBF_DiPFJet45_Mjj650_PNetTauhPFJet45_L2NN_eta2p3==1")

df_deeptau = df.Define("pass_ditau","PassDiTauDeepTau(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi)")
# df_deeptau = df_deeptau.Define("pass_singletau","HLT_LooseDeepTauPFTauHPS180_L2NN_eta2p1==1")
# df_deeptau = df_deeptau.Define("pass_ditaujet","has_jet==1 && HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60==1")
# df_deeptau = df_deeptau.Define("pass_vbfsingletau","has_vbfjets==1 && HLT_VBF_DiPFJet45_Mjj650_MediumDeepTauPFTauHPS45_L2NN_eta2p1==1")

# Show the leading tau's pt and id for the first few events
skim_branches = [
    "sig_muon_idx","sig_muon_pt", "sig_muon_eta", "sig_muon_phi", "sig_muon_mass",
    "sig_muon_iso", "sig_muon_mediumId", "match_sig_muon",
    "veto_ele", "muon_tau_dR", "muon_mt", "vis_mass",
    
    "probe_tau_idx","leading_tau_idDeepTauVSjet","leading_tau_decayMode", 
    "leading_tau_pt","leading_tau_eta","leading_tau_phi","leading_tau_mass", "match_probe_tau",

    "HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Medium_L2NN_eta2p3_CrossL1",
    "HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS35_L2NN_eta2p1_CrossL1",
    "HLT_IsoMu24",
    "weight", "pass_ditau",
]

df_pnet.Snapshot("Events", "zvbftest/test_pnet.root", skim_branches)
df_deeptau.Snapshot("Events", "zvbftest/test_deeptau.root", skim_branches)