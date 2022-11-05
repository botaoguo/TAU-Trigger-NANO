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

if args.type not in ['data', 'mc']:
    raise RuntimeError("Invalid sample type")

input_vec = ListToStdVector(args.input)
print("CHECK POINT 1")
if args.type == 'mc':
    if args.pu is None:
        raise RuntimeError("Pileup file should be provided for mc.")
    data_pu_file = ROOT.TFile(args.pu, 'READ')
    print("CHECK POINT 2")
    data_pu = data_pu_file.Get('pileup')
    df_all = ROOT.RDataFrame('Events', input_vec)
    print("CHECK POINT 3")
    mc_pu = df_all.Histo1D(ROOT.RDF.TH1DModel(data_pu), 'npu')
    ROOT.PileUpWeightProvider.Initialize(data_pu, mc_pu.GetPtr())
print("CHECK POINT 4")
trig_descriptors, channel_triggers = TriggerConfig.Load(args.config)
trigger_dict, filter_dict = TriggerConfig.LoadTriggerDictionary(input_vec)
print("CHECK POINT 5")
triggerMatch = ROOT.TriggerMatchProvider.Initialize()
print("CHECK POINT 6")
channels = {}
for channel_name, channel_trig_descs in channel_triggers.items():
    channel_id = ParseEnum(Channel, channel_name)
    channels[channel_name] = channel_id
    for desc in channel_trig_descs:
        if 'sample_types' in desc and args.type not in desc['sample_types']: continue
        if desc['leg_types'][-1] != 'tau': continue
        print("CHECK POINT 7")
        match_desc = ROOT.TriggerMatchProvider.MatchDescriptor()
        print("CHECK POINT 8")
        pattern = '^{}.*'.format(desc['name'])
        hlt_paths = TriggerConfig.GetMatchedTriggers(trigger_dict, pattern)
        match_desc.match_mask = int(TriggerConfig.GetMatchMask(hlt_paths))
        print("CHECK POINT 9")
        filter_names = desc['filters'][-1]
        #match_desc.filter_hashes = ListToStdVector([ filter_dict[f] for f in filter_names ], elem_type='UInt_t')
        print("CHECK POINT 10")
        if 'min_run' in desc and args.type == 'data':
            match_desc.min_run = desc['min_run']
        if 'max_run' in desc and args.type == 'data':
            match_desc.max_run = desc['max_run']
        sel_name = 'selection_' + channel_name
        if sel_name in desc:
            if 'hltObj_pt' in desc[sel_name]:
                match_desc.hltObj_pt = desc[sel_name]['hltObj_pt']
            if 'l1Tau_pt' in desc[sel_name]:
                match_desc.l1Tau_pt = desc[sel_name]['l1Tau_pt']
            if 'l1Tau_hwIso' in desc[sel_name]:
                match_desc.l1Tau_hwIso = desc[sel_name]['l1Tau_hwIso']
        triggerMatch.Add(channel_id, match_desc)
print("CHECK POINT 11")
selection_id = ParseEnum(TauSelection, args.selection)
df = ROOT.RDataFrame('Events', input_vec)
df = df.Filter('''
               (tau_sel & {}) != 0 && muon_pt > 27 && muon_iso < 0.1 && muon_mt < 30
               && tau_pt > 20 && abs(tau_eta) < 2.1 && tau_decayMode != 5 && tau_decayMode != 6
               && vis_mass > 40 && vis_mass < 80
               '''.format(selection_id))
print("CHECK POINT 12")
if selection_id == TauSelection.DeepTau:
    df = df.Filter('( (tau_idDeepTau2017v2p1VSmu << 2) & (1 << {})) != 0'.format(DiscriminatorWP.Tight))
if args.type == 'mc':
    df = df.Filter('tau_charge + muon_charge == 0 && tau_gen_match == 5')
    df = df.Define('weight', "PileUpWeightProvider::GetDefault().GetWeight(npu) * 1.0")
else:
    df = df.Define('weight', "muon_charge != tau_charge ? 1. : -1.")

skimmed_branches = [
    'tau_pt', 'tau_eta', 'tau_phi', 'tau_mass', 'tau_charge', 'tau_decayMode', 'weight',
    'tau_idDeepTau2017v2p1VSjet'
]

deltaRThr = 0.5
for channel_name, channel_id in channels.items():
    pass_branch = 'pass_' + channel_name
    df = df.Define(pass_branch, '''TriggerMatchProvider::GetDefault().Pass({}, run, tau_eta, tau_phi, hlt_accept, {},
                   hltObj_types, hltObj_pt, hltObj_eta, hltObj_phi, hltObj_hasPathName, filter_hltObj, filter_hash
                   )'''.format(channel_id, deltaRThr))
    skimmed_branches.append(pass_branch)
print("CHECK POINT 13")
df.Snapshot('Events', args.output, ListToStdVector(skimmed_branches))
print("CHECK POINT 14")
