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


# probe tau

# tau_pt > 40, tau_eta < 2.3, VSe >= 1, VSmu >= 1, VSjet >= 1
df = df.Define("probe_tau_idx","Check_SelectTau(Tau_pt, Tau_eta, Tau_idDeepTau2018v2p5VSe, Tau_idDeepTau2018v2p5VSmu, Tau_idDeepTau2018v2p5VSjet)")\
       .Filter("probe_tau_idx != -1")

df = df.Define("probe_secondtau_idx","Check_SelectSecondTau(Tau_pt, Tau_eta, Tau_idDeepTau2018v2p5VSe, Tau_idDeepTau2018v2p5VSmu, Tau_idDeepTau2018v2p5VSjet, probe_tau_idx)")\
       .Filter("probe_secondtau_idx != -1")

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
df = (
    df.Define("subleading_tau_pt", "Tau_pt[probe_secondtau_idx]")
    .Define("subleading_tau_eta", "Tau_eta[probe_secondtau_idx]")
    .Define("subleading_tau_phi", "Tau_phi[probe_secondtau_idx]")
    .Define("subleading_tau_mass", "Tau_mass[probe_secondtau_idx]")
    .Define("subleading_tau_charge", "Tau_charge[probe_secondtau_idx]")
    .Define("subleading_tau_decayMode", "Tau_decayMode[probe_secondtau_idx]")
    .Define("subleading_tau_idDeepTauVSe", "Tau_idDeepTau2018v2p5VSe[probe_secondtau_idx]")
    .Define("subleading_tau_idDeepTauVSmu", "Tau_idDeepTau2018v2p5VSmu[probe_secondtau_idx]")
    .Define("subleading_tau_idDeepTauVSjet", "Tau_idDeepTau2018v2p5VSjet[probe_secondtau_idx]")
)
# df = df.Define("has_2taus","probe_tau_idx != probe_secondtau_idx")
df = df.Define("weight","1")

# ditau
# deeptau wp==3
df = df.Define("trig_deeptau1_idx","MatchFirstTau(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_l1pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi, 3)")
df = df.Define("trig_deeptau2_idx","MatchSecondTau(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_l1pt, TrigObj_eta, TrigObj_phi, subleading_tau_pt, subleading_tau_eta, subleading_tau_phi, 3, trig_deeptau1_idx)")
df = df.Define("pass_ditau_deeptau_medium","HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1==1 && trig_deeptau1_idx!=-1 && trig_deeptau2_idx!=-1")
# pnet wp==4
df = df.Define("trig_pnettau1_idx","MatchFirstTau(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_l1pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi, 4)")
df = df.Define("trig_pnettau2_idx","MatchSecondTau(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_l1pt, TrigObj_eta, TrigObj_phi, subleading_tau_pt, subleading_tau_eta, subleading_tau_phi, 4, trig_pnettau1_idx)")
df = df.Define("pass_ditau_pnet_medium", "HLT_DoublePNetTauhPFJet30_Medium_L2NN_eta2p3==1 && trig_pnettau1_idx!=-1 && trig_pnettau2_idx!=-1")

skim_branches = [
    "probe_tau_idx", "probe_secondtau_idx",
    "leading_tau_idDeepTauVSjet", "leading_tau_idDeepTauVSe", "leading_tau_idDeepTauVSmu",
    "leading_tau_decayMode", "leading_tau_charge",
    "leading_tau_pt","leading_tau_eta","leading_tau_phi","leading_tau_mass",
    "subleading_tau_idDeepTauVSjet", "subleading_tau_idDeepTauVSe", "subleading_tau_idDeepTauVSmu",
    "subleading_tau_decayMode", "subleading_tau_charge",
    "subleading_tau_pt","subleading_tau_eta","subleading_tau_phi","subleading_tau_mass",
    
    # ditau trigger 
    "HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1",
    "HLT_DoublePNetTauhPFJet30_Medium_L2NN_eta2p3",
    
    "trig_deeptau1_idx", "trig_deeptau2_idx",
    "trig_pnettau1_idx", "trig_pnettau2_idx",
    "weight", "luminosityBlock", "run",
    "pass_ditau_deeptau_medium", "pass_ditau_pnet_medium",
]

output_file = args.output + "/" + "skimtuple_" + args.version + ".root"
# output_file = args.output + "/" + "skimtuple.root"
df.Snapshot("Events", output_file, skim_branches)

# Record the end time
end_time = time.time()
print("Totally cost {} seconds...".format(end_time - start_time))