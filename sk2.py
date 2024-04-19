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

# Open the ROOT file containing the tree
# file = ROOT.TFile.Open("/eos/user/b/ballmond/Tau_Trigger/Run2024C/Tau_Run2024C-PromptReco-v1_NANOAOD_1200b/f861fcac-160b-4965-9a6f-80a8e06b8d1f.root")
# tree = file.Get("Events")

# Create an RDataFrame from the tree
# df = ROOT.RDataFrame(tree)
input_vec = ListToStdVector(args.input)
df = ROOT.RDataFrame('Events', input_vec)

# Filter events with exactly 2 good tau candidates
# df = df.Filter("nTau >= 2")
# tau_pt > 20, tau_eta < 2.3, VSe >= 2, VSmu >= 1, VSjet >= 5
df = df.Define("has_tau","SelectTau(Tau_pt, Tau_eta, Tau_idDeepTau2018v2p5VSe, Tau_idDeepTau2018v2p5VSmu, Tau_idDeepTau2018v2p5VSjet)")\
       .Filter("has_tau == 1")
# Find a single good lepton, otherwise return -1 as index
df = df.Define("leading_tau_idx", "FindLeadingTau(Tau_pt)")\
       .Filter("leading_tau_idx != -1")
# Find the second leading tau
df = df.Define("subleading_tau_idx", "FindSubLeadingTau(Tau_pt, leading_tau_idx)")
# df = df.Define("leading_tau_pt","nTau == 2")

# Create new columns with the kinematics of good leptons
df = (
    df.Define("leading_tau_pt", "Tau_pt[leading_tau_idx]")
    .Define("leading_tau_eta", "Tau_eta[leading_tau_idx]")
    .Define("leading_tau_phi", "Tau_phi[leading_tau_idx]")
    .Define("leading_tau_mass", "Tau_mass[leading_tau_idx]")
    .Define("leading_tau_decayMode", "Tau_decayMode[leading_tau_idx]")
    .Define("leading_tau_idDeepTauVSe", "Tau_idDeepTau2018v2p5VSe[leading_tau_idx]")
    .Define("leading_tau_idDeepTauVSmu", "Tau_idDeepTau2018v2p5VSmu[leading_tau_idx]")
    .Define("leading_tau_idDeepTauVSjet", "Tau_idDeepTau2018v2p5VSjet[leading_tau_idx]")
)
df = (
    df.Define("subleading_tau_pt", "Tau_pt[subleading_tau_idx]")
    .Define("subleading_tau_eta", "Tau_eta[subleading_tau_idx]")
    .Define("subleading_tau_phi", "Tau_phi[subleading_tau_idx]")
    .Define("subleading_tau_mass", "Tau_mass[subleading_tau_idx]")
    .Define("subleading_tau_decayMode", "Tau_decayMode[subleading_tau_idx]")
    .Define("subleading_tau_idDeepTauVSe", "Tau_idDeepTau2018v2p5VSe[subleading_tau_idx]")
    .Define("subleading_tau_idDeepTauVSmu", "Tau_idDeepTau2018v2p5VSmu[subleading_tau_idx]")
    .Define("subleading_tau_idDeepTauVSjet", "Tau_idDeepTau2018v2p5VSjet[subleading_tau_idx]")
)

# apply offline leading tau pt cut
df = df.Filter("leading_tau_pt>20 && abs(leading_tau_eta)<2.3")
df = df.Filter("leading_tau_idDeepTauVSe>=2 && leading_tau_idDeepTauVSmu >=1 && leading_tau_idDeepTauVSjet >=5")
df = df.Filter("leading_tau_decayMode != 5 && leading_tau_decayMode != 6")
# apply to the second leading tau
df = df.Filter("subleading_tau_pt>20 && abs(subleading_tau_eta)<2.3")
df = df.Filter("subleading_tau_idDeepTauVSe>=2 && subleading_tau_idDeepTauVSmu >=1 && subleading_tau_idDeepTauVSjet >=5")
df = df.Filter("subleading_tau_decayMode != 5 && subleading_tau_decayMode != 6")


df = df.Define("has_vbfjets","SelectVBFJets(Jet_pt, Jet_eta, Jet_phi, Jet_mass, leading_tau_eta, leading_tau_phi)")
df = df.Define("has_jet","SelectJet(Jet_pt, Jet_eta, Jet_phi, Jet_mass, leading_tau_eta, leading_tau_phi)")

df = df.Define("weight","1")

# df = df.Define("pass_ditau_testpnet", "HLT_DoublePNetTauhPFJet30_Medium_L2NN_eta2p3==1 && PassDiTauPNet(nTrigObj, TrigObj_l1pt, TrigObj_l1iso, TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi)")
# df = df.Define("pass_ditau_testdeeptau", "HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1==1 && PassDiTauDeepTau(nTrigObj, TrigObj_l1pt, TrigObj_l1iso, TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi)")

# df = df.Define("pass_ditaujet_testpnet", "HLT_DoublePNetTauhPFJet26_L2NN_eta2p3_PFJet60==1 && PassDiTaujetPNet(nTrigObj, TrigObj_l1pt, TrigObj_l1iso, TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi)")
# df = df.Define("pass_ditaujet_testdeeptau", "HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60==1 && PassDiTaujetDeepTau(nTrigObj, TrigObj_l1pt, TrigObj_l1iso, TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi)")

df_pnet = df.Define("pass_ditau","HLT_DoublePNetTauhPFJet30_Medium_L2NN_eta2p3==1")
df_pnet = df_pnet.Define("pass_singletau","HLT_SinglePNetTauhPFJet130_Loose_L2NN_eta2p3==1")
df_pnet = df_pnet.Define("pass_ditaujet","has_jet==1 && HLT_DoublePNetTauhPFJet26_L2NN_eta2p3_PFJet60==1")
df_pnet = df_pnet.Define("pass_vbfsingletau","has_vbfjets==1 && HLT_VBF_DiPFJet45_Mjj650_PNetTauhPFJet45_L2NN_eta2p3==1")

df_deeptau = df.Define("pass_ditau","HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1==1")
df_deeptau = df_deeptau.Define("pass_singletau","HLT_LooseDeepTauPFTauHPS180_L2NN_eta2p1==1")
df_deeptau = df_deeptau.Define("pass_ditaujet","has_jet==1 && HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60==1")
df_deeptau = df_deeptau.Define("pass_vbfsingletau","has_vbfjets==1 && HLT_VBF_DiPFJet45_Mjj650_MediumDeepTauPFTauHPS45_L2NN_eta2p1==1")

# Show the leading tau's pt and id for the first few events
skim_branches = [
    "Tau_pt","nTau",
    "leading_tau_idx","leading_tau_idDeepTauVSjet","leading_tau_decayMode", 
    "leading_tau_pt","leading_tau_eta","leading_tau_phi","leading_tau_mass",
    "subleading_tau_idx","subleading_tau_idDeepTauVSjet","subleading_tau_decayMode", 
    "subleading_tau_pt","subleading_tau_eta","subleading_tau_phi","subleading_tau_mass",
    "has_tau", "has_vbfjets", "has_jet",
    
    # "HLT_DoublePNetTauhPFJet30_Medium_L2NN_eta2p3", # PNet, ditau
    # "HLT_SinglePNetTauhPFJet130_Loose_L2NN_eta2p3", # PNet, singletau
    # "HLT_DoublePNetTauhPFJet26_L2NN_eta2p3_PFJet60", # PNet, ditau+jet
    # "HLT_VBF_DiPFJet45_Mjj650_PNetTauhPFJet45_L2NN_eta2p3", # PNet, VBF singletau
    
    # "HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1", # DeepTau, ditau
    # "HLT_LooseDeepTauPFTauHPS180_L2NN_eta2p1", # DeepTau, singletau
    # "HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60", # DeepTau, ditau+jet
    # "HLT_VBF_DiPFJet45_Mjj550_MediumDeepTauPFTauHPS45_L2NN_eta2p1", # DeepTau, VBF singletau ?TODO
    # "HLT_VBF_DiPFJet45_Mjj650_MediumDeepTauPFTauHPS45_L2NN_eta2p1", # replace vbf singletau
    "pass_ditau", "pass_ditaujet", "pass_singletau", "pass_vbfsingletau",
    "weight",
    # "pass_ditau_testpnet", "pass_ditau_testdeeptau",
    # "pass_ditaujet_testpnet", "pass_ditaujet_testdeeptau",
]

df_pnet.Snapshot("Events", "zvbftest/test_pnet.root", skim_branches)
df_deeptau.Snapshot("Events", "zvbftest/test_deeptau.root", skim_branches)