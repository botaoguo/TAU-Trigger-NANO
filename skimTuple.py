#!/usr/bin/env python

import argparse
from array import array
import math
import numpy as np
import os
import re
import sys
import ROOT

parser = argparse.ArgumentParser(description='Skim full tuple.')
parser.add_argument('--input', required=True, type=str, nargs='+', help="input files")
parser.add_argument('--config', required=True, type=str, help="config with triggers description")
parser.add_argument('--selection', required=True, type=str, help="tau selection")
parser.add_argument('--output', required=True, type=str, help="output file")
parser.add_argument('--type', required=True, type=str, help="data or mc")
parser.add_argument('--pu', required=False, type=str, default=None,
                    help="file with the pileup profile for the data taking period")
args = parser.parse_args()

path_prefix = '' if 'TAU-Trigger-NANO' in os.getcwd() else 'TAU-Trigger-NANO/'
sys.path.insert(0, path_prefix + 'Common/python')
from AnalysisTypes import *
from AnalysisTools import *
import TriggerConfig
ROOT.ROOT.EnableImplicitMT(4)
ROOT.gROOT.SetBatch(True)
ROOT.gInterpreter.Declare('#include "{}interface/PyInterface.h"'.format(path_prefix))
ROOT.gInterpreter.Declare('#include "{}interface/picoNtupler.h"'.format(path_prefix))

if args.type not in ['data', 'mc']:
    raise RuntimeError("Invalid sample type")

input_vec = ListToStdVector(args.input)
if args.type == 'mc':
    if args.pu is None:
        raise RuntimeError("Pileup file should be provided for mc.")
    data_pu_file = ROOT.TFile(args.pu, 'READ')
    data_pu = data_pu_file.Get('pileup')
    df_all = ROOT.RDataFrame('Events', input_vec)
    mc_pu = df_all.Histo1D(ROOT.RDF.TH1DModel(data_pu), 'npu')
    ROOT.PileUpWeightProvider.Initialize(data_pu, mc_pu.GetPtr())


selection_id = ParseEnum(TauSelection, args.selection)
df = ROOT.RDataFrame('Events', input_vec)
df = df.Filter('''
               (tau_sel & {}) != 0 && muon_pt > 27 && muon_iso < 0.1 && muon_mt < 30
               && tau_pt > 20 && abs(tau_eta) < 2.1 && tau_decayMode != 5 && tau_decayMode != 6
               && vis_mass > 40 && vis_mass < 80
               '''.format(selection_id))
if selection_id == TauSelection.DeepTau:
    df = df.Filter('( (tau_idDeepTau2017v2p1VSmu << 2) & (1 << {})) != 0'.format(DiscriminatorWP.Tight))
if args.type == 'mc':
    df = df.Filter('tau_charge + muon_charge == 0 && tau_gen_match == 5')
    df = df.Define('weight', "PileUpWeightProvider::GetDefault().GetWeight(npu) * 1.0")
else:
    df = df.Define('weight', "muon_charge != tau_charge ? 1. : -1.")

skimmed_branches = [
    'tau_pt', 'tau_eta', 'tau_phi', 'tau_mass', 'tau_charge', 'tau_decayMode', 'weight', 'tau_idDeepTau2017v2p1VSjet',
# 2018 MC
#    'HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1', 
#    'HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1', 
#    'HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1', 
#    'HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1', 
#    'HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1',
#    'HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1',
  
# 2017 data
     'HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1',
#     'HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1',
     'HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1',
     'HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1',
     'HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1'
]

# PassMuTauFilter = "PassMuTauFilter(nTrigObj, TrigObj_id, TrigObj_filterBits, \
#                    TrigObj_pt, TrigObj_eta, TrigObj_phi, \
#                    tau_pt, tau_eta, tau_phi)"
# PassDiTauFilter = "PassDiTauFilter(nTrigObj, TrigObj_id, TrigObj_filterBits, \
#                    TrigObj_pt, TrigObj_eta, TrigObj_phi, \
#                    tau_pt, tau_eta, tau_phi)"
# PassElTauFilter = "PassElTauFilter(nTrigObj, TrigObj_id, TrigObj_filterBits, \
#                    TrigObj_pt, TrigObj_eta, TrigObj_phi, TrigObj_l1pt, TrigObj_l1iso,\
#                    tau_pt, tau_eta, tau_phi)"

deltaRThr = 0.5
#for channel_name, channel_id in channels.items():
#    pass_branch = 'pass_' + channel_name
#    df = df.Define(pass_branch, '''TriggerMatchProvider::GetDefault().Pass({}, run, tau_eta, tau_phi, hlt_accept, {},
#                   hltObj_types, hltObj_pt, hltObj_eta, hltObj_phi, hltObj_hasPathName, filter_hltObj, filter_hash
#                   )'''.format(channel_id, deltaRThr))
#    skimmed_branches.append(pass_branch)

# for 2018MC
#df = df.Define("pass_ditau", "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1 == 1")
#df = df.Define("pass_etau", "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1 == 1 && TrigObj_l1pt[0] > 26 && TrigObj_l1iso[0] > 0")
#df = df.Define("pass_mutau", "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1 == 1")

# for 2017Data
df = df.Define("pass_ditau", "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1 == 1")
df = df.Define("pass_etau", "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1 == 1 && TrigObj_l1pt[muon_trig_obj_idx] > 26 && TrigObj_l1iso[muon_trig_obj_idx] > 0 && tau_trig_obj_idx >= 0 && TrigObj_pt[tau_trig_obj_idx] > 30")
df = df.Define("pass_mutau", "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1 == 1")

skimmed_branches.append("pass_ditau")
skimmed_branches.append("pass_etau")
skimmed_branches.append("pass_mutau")

df.Snapshot('Events', args.output, ListToStdVector(skimmed_branches))
