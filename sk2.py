import os, sys
import ROOT
import argparse
import time

parser = argparse.ArgumentParser(description='Skim full tuple.')
parser.add_argument('--input', required=False, type=str, nargs='+', help="input files")
parser.add_argument('--inputlist', required=False, type=str, help="input files")
parser.add_argument('--type', required=False, type=str, default="data", help="mc or data")
parser.add_argument('--output', required=True, type=str, help="output file's dir")
parser.add_argument('--version', required=True, type=str, help="")
args = parser.parse_args()

path_prefix = '' if 'TAU-Trigger-NANO' in os.getcwd() else 'TAU-Trigger-NANO/'
sys.path.insert(0, path_prefix + 'Common/python')
from AnalysisTools import *

ROOT.gInterpreter.Declare('#include "tau_ntupler.h"')
# ROOT.gInterpreter.Declare('#include "{}interface/tau_ntupler.h"'.format(path_prefix))
# try:
#     ROOT.gInterpreter.Declare('#include "tau_ntupler.h"')
# except:
#     ROOT.gInterpreter.Declare('#include "{}interface/tau_ntupler.h"'.format(path_prefix))

# Record the start time
start_time = time.time()
if args.inputlist is not None:
    f = open(args.inputlist, "r")
    files = f.read().splitlines()
elif args.input is not None:
    files = args.input
input_vec = ListToStdVector(files)

df = ROOT.RDataFrame('Events', input_vec)

# JSON Filter
if args.type == "mc":
    pass
elif args.type == "data":
    df = df.Filter("jsonFilterlambda(run, luminosityBlock)")

# apply met filter
df = df.Filter("Flag_METFilters==1")

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
    .Define("sig_muon_charge", "Muon_charge[sig_muon_idx]")
    .Define("sig_muon_iso", "Muon_pfRelIso04_all[sig_muon_idx]")
    .Define("sig_muon_mediumId", "Muon_mediumId[sig_muon_idx]")
)
# match the trigobj and muon with dR < 0.5 and objId==13
df = df.Define("match_sig_muon","Muon_match(TrigObj_id, TrigObj_eta, TrigObj_phi, sig_muon_eta, sig_muon_phi)")\
       .Filter("match_sig_muon == 1")

# probe tau
# tau_pt > 20, tau_eta < 2.3, VSe >= 2, VSmu >= 4, VSjet >= 5
df = df.Define("probe_tau_idx","SelectTau(Tau_pt, Tau_eta, Tau_idDeepTau2018v2p5VSe, Tau_idDeepTau2018v2p5VSmu, Tau_idDeepTau2018v2p5VSjet)")\
       .Filter("probe_tau_idx != -1")

df = (
    df.Define("leading_tau_pt", "Tau_pt[probe_tau_idx]")
    .Define("leading_tau_eta", "Tau_eta[probe_tau_idx]")
    .Define("leading_tau_phi", "Tau_phi[probe_tau_idx]")
    .Define("leading_tau_mass", "Tau_mass[probe_tau_idx]")
    .Define("leading_tau_charge", "Tau_charge[probe_tau_idx]")
    .Define("leading_tau_decayMode", "Tau_decayMode[probe_tau_idx]")
    .Define("leading_tau_decayModePNet", "Tau_decayModePNet[probe_tau_idx]")
    .Define("leading_tau_idDeepTauVSe", "Tau_idDeepTau2018v2p5VSe[probe_tau_idx]")
    .Define("leading_tau_idDeepTauVSmu", "Tau_idDeepTau2018v2p5VSmu[probe_tau_idx]")
    .Define("leading_tau_idDeepTauVSjet", "Tau_idDeepTau2018v2p5VSjet[probe_tau_idx]")
)

# tau decay mode != 5 && != 6
df = df.Filter("leading_tau_decayMode != 5 && leading_tau_decayMode != 6")

# electron veto
df = df.Define("veto_ele", "VetoEle(Electron_pt, Electron_eta, Electron_mvaIso)")\
       .Filter("veto_ele == 1")

# second muon veto apply
df = df.Define("veto_SecMu","VetoSecondMu(Muon_pt, Muon_eta, Muon_looseId, sig_muon_idx)")\
       .Filter("veto_SecMu == 1")

# sig_muon and probe_tau dR > 0.5
df = df.Define("muon_tau_dR", "deltaR(sig_muon_eta, leading_tau_eta, sig_muon_phi, leading_tau_phi)")\
       .Filter("muon_tau_dR > 0.5")

# mt betweeb muon and met should < 30 GeV
df = df.Define("muon_mt", "MT(sig_muon_pt, PuppiMET_pt, sig_muon_phi, PuppiMET_phi)")\
       .Filter("muon_mt < 30")

# vis_mass in [40,80] GeV
df = df.Define("vis_mass", "VisMass(sig_muon_pt, sig_muon_eta, sig_muon_phi, sig_muon_mass, leading_tau_pt, leading_tau_eta, leading_tau_phi, leading_tau_mass)")\
       .Filter("vis_mass >= 40 && vis_mass <= 80")

# muon_charge + tau_charge == 0
# df = df.Filter("sig_muon_charge + leading_tau_charge == 0")

# apply to the mc numerator
# gen match with muon and tau
if args.type == "mc":
    df = df.Define("genmuon_idx","MuonGenMatch(GenPart_pdgId, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, sig_muon_eta, sig_muon_phi)")
    df = df.Define("gentau_idx","TauGenMatch(GenPart_pdgId, GenPart_statusFlags, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, leading_tau_eta, leading_tau_phi)")

# df = df.Define("weight","1")

skim_branches = [
    "sig_muon_idx","sig_muon_pt", "sig_muon_eta", "sig_muon_phi", "sig_muon_mass", "sig_muon_charge",
    "sig_muon_iso", "sig_muon_mediumId", "match_sig_muon",
    "veto_ele", "veto_SecMu", "muon_tau_dR", "muon_mt", "vis_mass",
    "probe_tau_idx","leading_tau_idDeepTauVSjet","leading_tau_decayMode", "leading_tau_charge",
    "leading_tau_pt","leading_tau_eta","leading_tau_phi","leading_tau_mass", "leading_tau_decayModePNet",
    "leading_tau_idDeepTauVSe", "leading_tau_idDeepTauVSmu",
    # ditau monitoring 
    "HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Medium_L2NN_eta2p3_CrossL1", # ditau pnet bit: 1,4,23
    "HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Tight_L2NN_eta2p3_CrossL1", # 2,4,23
    # "HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS35_L2NN_eta2p1_CrossL1", # ditau deeptau bit: 3,23
    # mutau monitoring
    "HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Loose_eta2p3_CrossL1", # mutau pnet bit: 0,4,13
    "HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Medium_eta2p3_CrossL1", # 1,4,13
    "HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Tight_eta2p3_CrossL1", # 2,4,13
    # "HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1", # mutau deeptau bit: 3,13
    # etau monitoring, same as mutau
    # etau pnet bit: loose 0,4,27
    # medium1,4,27
    # tight 2,4,27
    # etau deeptau bit: 3,27
    # ditaujet monitoring
    "HLT_IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1", # ditaujet pnet bit: 4,20
    # "HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_CrossL1", # ditaujet deeptau bit: 3,20
    # singletau
    "HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Loose_L2NN_eta2p3_CrossL1", # singletau pnet bit: 0,4,26
    "HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Medium_L2NN_eta2p3_CrossL1", # 1,4,26
    "HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Tight_L2NN_eta2p3_CrossL1", # 2,4,26
    # "HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS180_eta2p1", # singletau deeptau bit: 3,26
    # vbfsingletau
    "HLT_IsoMu24_eta2p1_PNetTauhPFJet45_L2NN_eta2p3_CrossL1", # vbfsingletau pnet bit: 4,30
    # "HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS45_L2NN_eta2p1_CrossL1", # vbfsingletau deeptau bit: 3,30
    # vbfditau
    "HLT_IsoMu24_eta2p1_PNetTauhPFJet20_eta2p2_SingleL1", # vbfditau pnet bit: 4,25
    # "HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS20_eta2p1_SingleL1", # vbfditau deeptau bit: 3,25
    
    "HLT_IsoMu24",
    # "weight", 
    "luminosityBlock", "run", "nTau",
    "nJet", "PV_npvs", "PV_npvsGood",
    
    # TrigObj related vars
    "TrigObj_eta", "TrigObj_filterBits", "TrigObj_id", "TrigObj_l1charge", "TrigObj_l1iso",
    "TrigObj_l1pt", "TrigObj_l1pt_2", "TrigObj_l2pt", "TrigObj_phi", "TrigObj_pt", "nTrigObj",
    
    # L1 related vars
    "L1_Mu18er2p1_Tau24er2p1", "L1_Mu18er2p1_Tau26er2p1",
    "L1_Mu22er2p1_IsoTau32er2p1", "L1_Mu22er2p1_IsoTau34er2p1",
    "L1_Mu22er2p1_IsoTau30er2p1", "L1_Mu22er2p1_Tau70er2p1",
    "L1_Mu18er2p1_Tau26er2p1_Jet55", "L1_Mu18er2p1_Tau26er2p1_Jet70",
    "L1_Mu22er2p1_IsoTau40er2p1", "L1_SingleMu22",
]
if args.type == "mc":
    # skim_branches.append("HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS35_L2NN_eta2p1_CrossL1_pHLT")
    # skim_branches.append("HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1_pHLT")
    # skim_branches.append("HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_CrossL1_pHLT")
    # skim_branches.append("HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS180_eta2p1_pHLT")
    # skim_branches.append("HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS45_L2NN_eta2p1_CrossL1_pHLT")
    # skim_branches.append("HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS20_eta2p1_SingleL1_pHLT")
    skim_branches.append("genmuon_idx")
    skim_branches.append("gentau_idx")
    
    

output_file = args.output + "/" + "skimtuple_" + args.version + ".root"
# output_file = args.output + "/" + "skimtuple.root"
df.Snapshot("Events", output_file, skim_branches)

# # Save the cutflow histogram to the output file
# output_root_file = ROOT.TFile(output_file, "UPDATE")  # Open in UPDATE mode to add objects
# cutflow.Write()  # Write the cutflow histogram
# output_root_file.Close()

# Record the end time
end_time = time.time()
print("Totally cost {} seconds...".format(end_time - start_time))
