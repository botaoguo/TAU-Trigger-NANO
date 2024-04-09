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
# parser.add_argument('--config', required=True, type=str, help="config with triggers description")
parser.add_argument('--selection', required=True, type=str, help="tau selection")
parser.add_argument('--output', required=True, type=str, help="output file")
parser.add_argument('--type', required=True, type=str, help="data or mc")
parser.add_argument('--pudata', required=False, type=str, default=None,
                    help="file with the pileup profile for the data taking period")
parser.add_argument('--pumc', required=False, type=str, default=None,
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
    if args.pudata is None or args.pumc is None:
        raise RuntimeError("Pileup file should be provided for mc.")
    data_pu_file = ROOT.TFile(args.pudata, 'READ')
    data_pu = data_pu_file.Get('pileup')
    df_all = ROOT.RDataFrame('Events', input_vec)
    mc_pu_file = ROOT.TFile(args.pumc, 'READ')
    mc_pu = mc_pu_file.Get('pileup')
    # mc_pu = df_all.Histo1D(ROOT.RDF.TH1DModel(data_pu), 'npu')
    # ROOT.PileUpWeightProvider.Initialize(data_pu, mc_pu.GetPtr())
    ROOT.PileUpWeightProvider.Initialize(data_pu, mc_pu)


selection_id = ParseEnum(TauSelection, args.selection)
df = ROOT.RDataFrame('Events', input_vec)
df = df.Filter('''
               (tau_sel & {}) != 0 && muon_pt > 24 && muon_iso < 0.1 && muon_mt < 30
               && tau_pt > 20 && abs(tau_eta) < 2.1 && tau_decayMode != 5 && tau_decayMode != 6
               && vis_mass > 40 && vis_mass < 80
               '''.format(selection_id))
if selection_id == TauSelection.DeepTau:
    # df = df.Filter('tau_idDeepTau2017v2p1VSmu >= {}'.format(DiscriminatorWP.Medium))
    # df = df.Filter('( tau_idDeepTau2017v2p1VSmu  & (1 << {})) != 0'.format(DiscriminatorWP.Medium))
    # df = df.Filter('tau_idDeepTau2018v2p5VSmu >= {}'.format(DiscriminatorWP.Medium))
    # df = df.Filter('( tau_idDeepTau2018v2p5VSmu  & (1 << {})) != 0'.format(DiscriminatorWP.Medium))
    df = df.Filter('( tau_idDeepTau2018v2p5VSmu  & 4) != 0')
if args.type == 'mc':
    df = df.Filter('tau_charge + muon_charge == 0 && tau_gen_match == 5')
    df = df.Define('weight', "PileUpWeightProvider::GetDefault().GetWeight(npu) * 1.0")
else:
    df = df.Define('weight', "muon_charge != tau_charge ? 1. : -1.")

skimmed_branches = [
    'tau_pt', 'tau_eta', 'tau_phi', 'tau_mass', 'tau_charge', 'tau_decayMode', 'weight', 'tau_idDeepTau2017v2p1VSjet', 'tau_idDeepTau2018v2p5VSjet',
    'TrigObj_l1pt', 'TrigObj_l1iso', 'nTrigObj'
    # use monitoring path, as TnP won't work in HLT path
]

df = df.Define("pass_mutau", "PassMuTauTrig2022(nTrigObj, TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, tau_pt, tau_eta, tau_phi)")
df = df.Define("pass_etau", "PassEleTauTrig2022(nTrigObj, TrigObj_l1pt, TrigObj_l1iso, TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, tau_pt, tau_eta, tau_phi)")

# ditau -> drop !bit18 cut, change l1pt>32 with l1pt>=32
df = df.Define("pass_ditau", "PassDiTauTrig2022(nTrigObj, TrigObj_l1pt, TrigObj_l1iso, TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, tau_pt, tau_eta, tau_phi)")
# 
df = df.Define("pass_ditaujet", "PassDiTauJetTrig2022(nTrigObj, TrigObj_l1pt, TrigObj_l1iso, TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, tau_pt, tau_eta, tau_phi)")

skimmed_branches.append("pass_ditau")
skimmed_branches.append("pass_etau")
skimmed_branches.append("pass_mutau")
skimmed_branches.append("pass_ditaujet")

df.Snapshot('Events', args.output, ListToStdVector(skimmed_branches))
print("Check Point")