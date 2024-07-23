import os, sys
import ROOT
import argparse
import time

parser = argparse.ArgumentParser(description='Skim full tuple.')
parser.add_argument('--input', required=False, type=str, nargs='+', help="input files")
parser.add_argument('--inputlist', required=False, type=str, help="input files")
parser.add_argument('--output', required=True, type=str, help="output file's dir")
parser.add_argument('--version', required=True, type=str, help="")
args = parser.parse_args()

path_prefix = '' if 'TAU-Trigger-NANO' in os.getcwd() else 'TAU-Trigger-NANO/'
sys.path.insert(0, path_prefix + 'Common/python')
from AnalysisTools import *

ROOT.gInterpreter.Declare('#include "{}interface/tau_ntupler.h"'.format(path_prefix))

# Record the start time
start_time = time.time()
if args.inputlist is not None:
    f = open(args.inputlist, "r")
    files = f.read().splitlines()
elif args.input is not None:
    files = args.input
input_vec = ListToStdVector(files)

df = ROOT.RDataFrame('Events', input_vec)


# filter for 382730 and 382811
df = df.Filter("( run==382730 && (luminosityBlock >= 1 && luminosityBlock <= 52) ) || ( run==382811 && (luminosityBlock >= 1 && luminosityBlock <= 396) )")

# tau_pt > 40, tau_eta < 2.3, VSe >= 1, VSmu >= 1, VSjet >= 1
# tag the signal muon
df = df.Define("sig_muon_idx","SelectMuon(Muon_pt, Muon_eta, Muon_phi, Muon_pfRelIso04_all, Muon_mediumId)")\
       .Filter("sig_muon_idx != -1")
df = df.Filter("HLT_IsoMu24==1")

df = (
    df.Define("sig_muon_pt", "Muon_pt[sig_muon_idx]")
    .Define("sig_muon_eta", "Muon_eta[sig_muon_idx]")
    .Define("sig_muon_phi", "Muon_phi[sig_muon_idx]")
    .Define("sig_muon_mass", "Muon_mass[sig_muon_idx]")
    .Define("sig_muon_charge", "Muon_charge[sig_muon_idx]")
    .Define("sig_muon_iso", "Muon_pfRelIso04_all[sig_muon_idx]")
    .Define("sig_muon_mediumId", "Muon_mediumId[sig_muon_idx]")
)
# match the trigobj and muon with dR < 0.5 and objId==13
df = df.Define("match_sig_muon","Muon_match(TrigObj_id, TrigObj_eta, TrigObj_phi, sig_muon_eta, sig_muon_phi)")\
       .Filter("match_sig_muon == 1")

# probe tau
df = df.Define("probe_tau_idx","Check_SelectTau(Tau_pt, Tau_eta, Tau_idDeepTau2018v2p5VSe, Tau_idDeepTau2018v2p5VSmu, Tau_idDeepTau2018v2p5VSjet)")\
       .Filter("probe_tau_idx != -1")

df = (
    df.Define("leading_tau_pt", "Tau_pt[probe_tau_idx]")
    .Define("leading_tau_eta", "Tau_eta[probe_tau_idx]")
    .Define("leading_tau_phi", "Tau_phi[probe_tau_idx]")
    .Define("leading_tau_mass", "Tau_mass[probe_tau_idx]")
    .Define("leading_tau_charge", "Tau_charge[probe_tau_idx]")
    .Define("leading_tau_decayMode", "Tau_decayMode[probe_tau_idx]")
    .Define("leading_tau_idDeepTauVSe", "Tau_idDeepTau2018v2p5VSe[probe_tau_idx]")
    .Define("leading_tau_idDeepTauVSmu", "Tau_idDeepTau2018v2p5VSmu[probe_tau_idx]")
    .Define("leading_tau_idDeepTauVSjet", "Tau_idDeepTau2018v2p5VSjet[probe_tau_idx]")
)
df = df.Define("weight","1")

# ditau monitoring
# ditau
df = df.Define("pass_ditau_pnet_medium","HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Medium_L2NN_eta2p3_CrossL1==1 && PassDiTauPNet(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_l1pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi, 1)")
df = df.Define("pass_ditau_pnet_tight","HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Tight_L2NN_eta2p3_CrossL1==1 && PassDiTauPNet(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_l1pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi, 2)")
df = df.Define("pass_ditau_deeptau","HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS35_L2NN_eta2p1_CrossL1==1 && PassDiTauDeepTau(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_l1pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi)")


skim_branches = [
    "probe_tau_idx", "sig_muon_idx", "match_sig_muon",
    
    "sig_muon_pt", "sig_muon_eta", "sig_muon_phi", "sig_muon_mass",
    "sig_muon_charge", "sig_muon_iso", "sig_muon_mediumId",
    
    "leading_tau_idDeepTauVSjet", "leading_tau_idDeepTauVSe", "leading_tau_idDeepTauVSmu",
    "leading_tau_decayMode", "leading_tau_charge",
    "leading_tau_pt","leading_tau_eta","leading_tau_phi","leading_tau_mass",
    
    # ditau monitoring
    "HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Medium_L2NN_eta2p3_CrossL1",
    "HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Tight_L2NN_eta2p3_CrossL1",
    "HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS35_L2NN_eta2p1_CrossL1",
    
    "weight", "luminosityBlock", "run",
    "pass_ditau_deeptau", "pass_ditau_pnet_medium", "pass_ditau_pnet_tight",
]

output_file = args.output + "/" + "skimtuple_" + args.version + ".root"
# output_file = args.output + "/" + "skimtuple.root"
df.Snapshot("Events", output_file, skim_branches)

# Record the end time
end_time = time.time()
print("Totally cost {} seconds...".format(end_time - start_time))