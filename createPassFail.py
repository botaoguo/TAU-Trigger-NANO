import os, sys
import ROOT
import argparse
import time

parser = argparse.ArgumentParser(description='Skim full tuple.')
parser.add_argument('--input', required=False, type=str, nargs='+', help="input files")
parser.add_argument('--inputlist', required=False, type=str, help="input files")
parser.add_argument('--type', required=False, type=str, default="data", help="mc or data")
parser.add_argument('--output', required=True, type=str, help="output file's dir")
# parser.add_argument('--version', required=True, type=str, help="")
args = parser.parse_args()

path_prefix = '' if 'TAU-Trigger-NANO' in os.getcwd() else 'TAU-Trigger-NANO/'
sys.path.insert(0, path_prefix + 'Common/python')
# from AnalysisTools import *
ROOT.gInterpreter.Declare('#include "tau_ntupler.h"')

def ListToStdVector(l, elem_type='string'):
    v = ROOT.std.vector(elem_type)()
    for x in l:
        if elem_type in ['Int_t', 'UInt_t']:
            x = int(x)
        v.push_back(x)
    return v

# Record the start time
start_time = time.time()
if args.inputlist is not None:
    f = open(args.inputlist, "r")
    files = f.read().splitlines()
elif args.input is not None:
    files = args.input
input_vec = ListToStdVector(files)

df = ROOT.RDataFrame('Events', input_vec)

# mutau 
# L1 eff
df = df.Define("pass_mutau_L1", "L1_Mu18er2p1_Tau24er2p1==1 || L1_Mu18er2p1_Tau26er2p1==1")
df = df.Define("pass_ditau_L1", "L1_Mu22er2p1_IsoTau32er2p1==1 || L1_Mu22er2p1_IsoTau34er2p1==1 || L1_Mu22er2p1_Tau70er2p1==1")

# mutau
df = df.Define("pass_mutau_pnet_loose","HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Loose_eta2p3_CrossL1==1 && PassMuTauPNet(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi, 0)")
df = df.Define("pass_mutau_pnet_medium","HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Medium_eta2p3_CrossL1==1 && PassMuTauPNet(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi, 1)")
df = df.Define("pass_mutau_pnet_tight","HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Tight_eta2p3_CrossL1==1 && PassMuTauPNet(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi, 2)")
# etau
df = df.Define("pass_etau_pnet_loose_withmutaubit","HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Loose_eta2p3_CrossL1==1 && PassETauPNet_withMutaubit(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_l1iso, TrigObj_l1pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi, 0)")
df = df.Define("pass_etau_pnet_medium_withmutaubit","HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Medium_eta2p3_CrossL1==1 && PassETauPNet_withMutaubit(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_l1iso, TrigObj_l1pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi, 1)")
df = df.Define("pass_etau_pnet_tight_withmutaubit","HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Tight_eta2p3_CrossL1==1 && PassETauPNet_withMutaubit(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_l1iso, TrigObj_l1pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi, 2)")
# ditau
df = df.Define("pass_ditau_pnet_medium","HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Medium_L2NN_eta2p3_CrossL1==1 && PassDiTauPNet(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_l1pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi, 1)")
df = df.Define("pass_ditau_pnet_tight","HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Tight_L2NN_eta2p3_CrossL1==1 && PassDiTauPNet(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_l1pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi, 2)")
# ditaujet
df = df.Define("pass_ditaujet_pnet","HLT_IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1==1 && PassDiTaujetPNet(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_l1iso, TrigObj_l1pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi)")
# singletau no filter bit check
df = df.Define("pass_singletau_pnet_loose_nofilterbit","HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Loose_L2NN_eta2p3_CrossL1==1 && PassSingleTauPNet_nofilterbit(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_l1pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi, 0)")
df = df.Define("pass_singletau_pnet_medium_nofilterbit","HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Medium_L2NN_eta2p3_CrossL1==1 && PassSingleTauPNet_nofilterbit(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_l1pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi, 1)")
df = df.Define("pass_singletau_pnet_tight_nofilterbit","HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Tight_L2NN_eta2p3_CrossL1==1 && PassSingleTauPNet_nofilterbit(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_l1pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi, 2)")
# vbfsingletau check
df = df.Define("pass_vbfsingletau_pnet_nofilterbit","HLT_IsoMu24_eta2p1_PNetTauhPFJet45_L2NN_eta2p3_CrossL1==1 && PassVBFSingleTauPNet_nofilterbit(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_l1iso, TrigObj_l1pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi)")

if args.type == "mc":
    # df = df.Define("pass_mutau_deeptau","HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1_pHLT==1 && PassMuTauDeepTau(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi)")
    # df = df.Define("pass_etau_deeptau_withmutaubit","HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1_pHLT==1 && PassETauDeepTau_withMutaubit(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_l1iso, TrigObj_l1pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi)")
    # df = df.Define("pass_ditau_deeptau","HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS35_L2NN_eta2p1_CrossL1_pHLT==1 && PassDiTauDeepTau(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_l1pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi)")
    # df = df.Define("pass_ditaujet_deeptau","HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_CrossL1_pHLT==1 && PassDiTaujetDeepTau(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_l1iso, TrigObj_l1pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi)")
    # df = df.Define("pass_singletau_deeptau_nofilterbit","HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS180_eta2p1_pHLT==1 && PassSingleTauDeepTau_nofilterbit(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_l1pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi)")
    # df = df.Define("pass_vbfsingletau_deeptau_nofilterbit","HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS45_L2NN_eta2p1_CrossL1_pHLT==1 && PassVBFSingleTauDeepTau_nofilterbit(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_l1iso, TrigObj_l1pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi)")
    df = df.Filter("sig_muon_charge + leading_tau_charge == 0")
    df = df.Filter("genmuon_idx!=-1 && gentau_idx!=-1")
    # need to add pileup reweighting
    df = df.Define("weight","1")
else:
    # df = df.Define("pass_mutau_deeptau","HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1==1 && PassMuTauDeepTau(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi)")
    # df = df.Define("pass_etau_deeptau_withmutaubit","HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1==1 && PassETauDeepTau_withMutaubit(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_l1iso, TrigObj_l1pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi)")
    # df = df.Define("pass_ditau_deeptau","HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS35_L2NN_eta2p1_CrossL1==1 && PassDiTauDeepTau(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_l1pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi)")
    # df = df.Define("pass_ditaujet_deeptau","HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_CrossL1==1 && PassDiTaujetDeepTau(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_l1iso, TrigObj_l1pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi)")
    # df = df.Define("pass_singletau_deeptau_nofilterbit","HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS180_eta2p1==1 && PassSingleTauDeepTau_nofilterbit(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_l1pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi)")
    # df = df.Define("pass_vbfsingletau_deeptau_nofilterbit","HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS45_L2NN_eta2p1_CrossL1==1 && PassVBFSingleTauDeepTau_nofilterbit(TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_l1iso, TrigObj_l1pt, TrigObj_eta, TrigObj_phi, leading_tau_pt, leading_tau_eta, leading_tau_phi)")
    df = df.Define("weight", "sig_muon_charge != leading_tau_charge ? 1. : -1.")

output_file = args.output
# + "/" + "skimtuple_2025C.root"
# output_file = args.output + "/" + "skimtuple.root"
df.Snapshot("Events", output_file)

# # Save the cutflow histogram to the output file
# output_root_file = ROOT.TFile(output_file, "UPDATE")  # Open in UPDATE mode to add objects
# cutflow.Write()  # Write the cutflow histogram
# output_root_file.Close()

# Record the end time
end_time = time.time()
print("Totally cost {} seconds...".format(end_time - start_time))
