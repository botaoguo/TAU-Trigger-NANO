#!/usr/bin/env python

import argparse
from array import array
import math
import numpy as np
import os
import re
import sys
import ROOT

import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Create turn on curves.')
parser.add_argument('--input', required=True, type=str, nargs='+', help="input files")
parser.add_argument('--output', required=True, type=str, help="output file prefix")
parser.add_argument('--types-trigger', required=False, type=str, default='pnet_loose', help="trigger type")
# parser.add_argument('--types-pnet', required=False, type=str, default='pnet_loose,pnet_medium,pnet_tight', help="pnet type")
# parser.add_argument('--types-deeptau', required=False, type=str, default='deeptau', help="deeptau type")
parser.add_argument('--channels', required=False, type=str, default='mutau', help="channels to process")
parser.add_argument('--decay-modes', required=False, type=str, default='all,0,1,10,11,1011', help="decay modes to process")
parser.add_argument('--working-points', required=False, type=str,
                    default='Medium',
                    help="working points to process")
parser.add_argument('--var', required=False, type=str, default='leading_tau_pt', help="plot var")
parser.add_argument('--vbfditau_ptcut', required=False, type=int, default=30, help="VBF ditau pt cut on eta and phi")
parser.add_argument('--mutau_ptcut', required=False, type=int, default=27, help="mutau cut on eta and phi, PV_npvsGood")
args = parser.parse_args()

path_prefix = '' if 'TAU-Trigger-NANO' in os.getcwd() else 'TAU-Trigger-NANO/'
sys.path.insert(0, path_prefix + 'Common/python')
from AnalysisTypes import *
from AnalysisTools import *
import RootPlotting
ROOT.ROOT.EnableImplicitMT(4)
ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2()
RootPlotting.ApplyDefaultGlobalStyle()

def CreateBins(max_pt, for_fitting):
        #bins = np.arange(20, 100, step=10)
        
        bins = np.arange(20, 40, step=4)
        bins = np.append(bins, np.arange(40, 60, step=5))
        bins = np.append(bins, np.arange(60, 100, step=10))
        high_pt_bins = [ 100, 150, 200, 300, 400, 500, 650, 800, 1000 ]
        
        # bins = np.arange(20, 100, step=8)
        # high_pt_bins = [ 100, 200, 300, 400, 500, 650, 800, 1000 ]
        
        n = 0
        while n < len(high_pt_bins) and high_pt_bins[n] < max_pt:
            n += 1
        use_logx = max_pt > 200
        return np.append(bins, high_pt_bins[0:n+1]), use_logx

ptcut_dict = {
    "mutau": args.mutau_ptcut, # pnet 27, deeptau 27
    "etau": 30,
    "ditau": 35, # pnet 30, deeptau 35,
    # "ditau": 40,
    "ditaujet": 30, # pnet 26, deeptau 30
    "vbfditau": args.vbfditau_ptcut, # pnet 20, deeptau 20, extra cut 25, 30
    "vbfsingletau": 45, # pnet 45, deeptau 45
}    

class TurnOnData:
    def __init__(self):
        self.hist_total = None
        self.hist_passed = None
        self.eff = None
        self.pass_yield = 0
        self.total_yield = 0

def CreateHistograms(input_file, channels, decay_modes, discr_name, working_points, hist_models, label, var,
                     output_file, era):
    df = ROOT.RDataFrame('Events', input_file)
    turnOn_data = {}
    dm_labels = {}

    for dm in decay_modes:
        if dm == 'all':
            dm_labels[dm] = ''
            df_dm = df
        elif dm == '1011':
            dm_labels[dm] = '_dm{}'.format(dm)
            df_dm = df.Filter('leading_tau_decayMode == 10 | leading_tau_decayMode == 11')
        else:
            dm_labels[dm] = '_dm{}'.format(dm)
            df_dm = df.Filter('leading_tau_decayMode == {}'.format(dm))
        turnOn_data[dm] = {}
        for wp in working_points:
            wp_bit = ParseEnum(DiscriminatorWP, wp)
            # df_wp = df_dm.Filter('({} & (1 << {})) != 0'.format(discr_name, wp_bit))
            df_wp = df_dm.Filter('{0} >= {1}'.format(discr_name, wp_bit))
            turnOn_data[dm][wp] = {}
            for channel in channels:
                turnOn_data[dm][wp][channel] = {}
                # add tau pt cut for any variable not tau pt (eta, phi, PV, era sum)
                if "pt" not in var:
                    df_wp = df_wp.Filter("leading_tau_pt>{}".format(ptcut_dict[channel]))
                if var == "nTau":
                    df_wp = df_wp.Filter('abs(leading_tau_eta) < 0.2')
                    df_wp = df_wp.Filter('abs(leading_tau_phi) < 0.3')
                # df_wp = df_wp.Filter('run <= 382262')
                # df_wp = df_wp.Filter('abs(sig_muon_eta) < 1.1')
                df_ch = df_wp.Filter('pass_{}_{} > 0.5'.format(channel, label))
                for model_name, hist_model in hist_models.items():
                    turn_on = TurnOnData()
                    turn_on.hist_total = df_wp.Histo1D(hist_model, var, 'weight')
                    turn_on.hist_passed = df_ch.Histo1D(hist_model, var, 'weight')
                    turn_on.total_yield = turn_on.hist_total.Integral()
                    turn_on.pass_yield = turn_on.hist_passed.Integral()
                    # for i in range(turn_on.hist_passed.GetNbinsX() + 1):
                    #     print("pass event: {}".format(turn_on.hist_passed.GetBinContent(i)))
                    turnOn_data[dm][wp][channel][model_name] = turn_on

    for dm in decay_modes:
        for wp in working_points:
            for channel in channels:
                for model_name, hist_model in hist_models.items():
                    turn_on = turnOn_data[dm][wp][channel][model_name]
                    name_pattern = '{}_{}_{}{}_{}_{{}}_{}'.format(label, channel, wp, dm_labels[dm], model_name, era)
                    turn_on.name_pattern = name_pattern
                    if 'fit' in model_name:
                        passed, total, eff = AutoRebinAndEfficiency(turn_on.hist_passed.GetPtr(),
                                                                    turn_on.hist_total.GetPtr(), bin_scan_pairs)
                    else:
                        passed, total = turn_on.hist_passed.GetPtr(), turn_on.hist_total.GetPtr()
                        try:
                            FixEfficiencyBins(passed, total)
                        except:
                            continue
                        turn_on.eff = ROOT.TEfficiency(passed, total)
                        eff = turn_on.eff
                    output_file.WriteTObject(total, name_pattern.format('total'), 'Overwrite')
                    output_file.WriteTObject(passed, name_pattern.format('passed'), 'Overwrite')
                    output_file.WriteTObject(eff, name_pattern.format('eff'), 'Overwrite')

    return turnOn_data

out_dir = os.path.dirname(args.output)
os.system("mkdir -p {}".format(out_dir))
output_file = ROOT.TFile(args.output + '.root', 'RECREATE')

input_files = args.input
print("input file: {}".format(input_files))

# labels = args.types_pnet.split(',')
# put deeptau in the last postion, take as denominator
# labels += args.types_deeptau.split(',')
labels = args.types_trigger.split(',')

labels = labels*len(input_files)
colors_list = [ ROOT.kRed, ROOT.kGreen, ROOT.kViolet, ROOT.kBlue, ROOT.kOrange, ROOT.kCyan, ROOT.kGray, ROOT.kBlack]
# colors_list = [ ROOT.kRed, ROOT.kBlack, ROOT.kRed, ROOT.kBlack, ROOT.kRed, ROOT.kBlack]
colors = colors_list[0:len(labels)-1] + [colors_list[-1]]
marker = [20,21,22] # 20 circle, 21 square, 22 triangle
n_inputs = len(labels)
var = args.var
decay_modes = args.decay_modes.split(',')
channels = args.channels.split(',')
working_points = args.working_points.split(',')
if "pt" in var:
    bins, use_logx = CreateBins(200, False)
    title, x_title = '#tau p_{T}', '#tau p_{T} (GeV)'
elif "era" in var:
    bins = np.linspace(0, 1, 1) # single bin for a single value
    # D, E, F, G, H, I
    # add label map 0-D, 1-E, 2-F, 3-G, 4-H, 5-I
    title, x_title = 'Efficiency per Era', 'Era'
elif "eta" in var:
    bins = np.linspace(-2.3, 2.3, 12)
    # bins = np.arange(-2.3, 2.3, step=0.2)
    # bins = np.append(bins, [2.3])
    use_logx = False
    title, x_title = '#tau #eta', '#tau #eta'
elif "phi" in var:
    bins = np.linspace(-3.2, 3.2, 12)
    # bins = np.arange(-3.2, 3.2, step=0.2)
    # bins = np.append(bins, [3.2])
    use_logx = False
    title, x_title = '#tau #phi', '#tau #phi'
elif var == "nTau":
    bins = np.linspace(0, 10, 10)
    use_logx = False
    title, x_title = 'nTau', 'nTau'
elif "PV_npvs" in var:
    bins = np.linspace(0, 60, 7)
    bins = np.append(bins, [80])
    # bins = np.arange(-3.2, 3.2, step=0.2)
    # bins = np.append(bins, [3.2])
    use_logx = False
    if "Good" in var:
        title, x_title = 'PV_npvsGood', 'PV_npvsGood'
    else:
        title, x_title = 'PV_npvs', 'PV_npvs'

hist_models = {
    'plot': ROOT.RDF.TH1DModel(var, var, len(bins) - 1, array('d', bins)),
}
turnOn_data = [None] * n_inputs
# using label as deeptau and pnet etc.
# turnOn_data[0] -> pnet_loose
# turnOn_data[1] -> pnet_medium
# turnOn_data[2] -> pnet_tight
# turnOn_data[3] -> deeptau
for input_id in range(n_inputs):
    fileprefix = os.path.basename(input_files[int(input_id)])
    fileprefix = os.path.splitext(fileprefix)[0]
    era = fileprefix.strip("skimtuple_")
    print("Creating {0} histograms for {1}...".format(labels[input_id],era))
    turnOn_data[input_id] = CreateHistograms(input_files[int(input_id)], channels, decay_modes, 'leading_tau_idDeepTauVSjet', # tau_idDeepTau2017v2p1VSjet,
                                             working_points, hist_models, labels[input_id], var, output_file, era)
# DM, WP, Trigger
print(f"n_inputs: {n_inputs}")

canvas = RootPlotting.CreateCanvas()
n_plots = len(decay_modes) * len(channels) * len(working_points)
plot_id = 0
print(f"channels: {channels}")
print(f"working_points: {working_points}")
print(f"decay_modes: {decay_modes}")
for channel in channels:
    for wp in working_points:
        fig, ax = plt.subplots()
        for dm in decay_modes:
            print()
            print(f"dm = {dm}")
            if dm == 'all':
                dm_label = ''
                dm_plain_label = ''
            else:
                dm_label = ' DM={}'.format(dm)
                dm_plain_label = '_dm{}'.format(dm)
                continue
            ratio_graph = [None] * n_inputs
            ref_hist = hist_models['plot'].GetHistogram()
            ratio_ref_hist = ref_hist.Clone()
            turnOns = [None] * n_inputs
            curves = [None] * n_inputs
            pass_yield = [None] * n_inputs
            total_yield = [None] * n_inputs
            ratio_yield = [None] * n_inputs
            ratio_error = [None] * n_inputs # error by hand
            # turnOn_data[0] -> pnet_loose
            # turnOn_data[1] -> pnet_medium
            # turnOn_data[2] -> pnet_tight
            # turnOn_data[3] -> deeptau
            for input_id in range(n_inputs):
                turnOns[input_id] = turnOn_data[input_id][dm][wp][channel]['plot']
                curves[input_id] = turnOns[input_id].eff
                pass_yield[input_id] = turnOns[input_id].pass_yield
                total_yield[input_id] = turnOns[input_id].total_yield
                ratio_yield[input_id] = pass_yield[input_id]/total_yield[input_id]
                ratio_error[input_id] = ratio_yield[input_id]*np.sqrt(1/pass_yield[input_id] + 1/total_yield[input_id])
                # the error bars on a ratio plot of a histogram A divided by a histogram B are:
                # (A/B) * sqrt[ (errA / A)^2 + (errB / B)^2 ]
                # for data, error = sqrt [N] , where A and B are simply N events in a bin
                # so if numerator and denominaotr are both data
                # ratio error = (A/B) * sqrt [ (1/A) + (1/B) ]
            print("ratio_yield")
            print(ratio_yield)
            print("ratio_error")
            print(ratio_error)
            # baseline plot
            dummyX = np.arange(6)
            ax.errorbar(dummyX, ratio_yield, xerr=0.5, yerr=ratio_error, # dummy yerr for now
                        color="Black", marker="o", markersize=4, label=f"DM = {dm}") # linestyle="none"
            ax.set_ylim(0, 1.2)
            ax.set_ylabel("L1 + HLT Efficiency", fontsize=16)
            ax.legend()
            # manually overwriting tick labels with corresponding eras
            use_era_labels = True
            if use_era_labels:
              ax.set_xlabel("Eras", fontsize=16) # eras
              eras_label = ["D", "E", "F", "G", "H", "I"] 
              ax.set_xticks(dummyX)
              ax.set_xticklabels(labels=eras_label, ha="center", fontsize=10)
            else: #use_date_labels
              ax.set_xticklabels(labels=[]) # turn labels off if using dates
              #ax.set_xlabel("Date", fontsize=16)
              ax.text(0.95, -0.10, "Date", transform=ax.transAxes, fontsize=16, ha="center", va="center") # instead of set xlabel
              date_txt1, date_txt2, date_txt3, date_txt4 = "26/06", "11/07", "26/07", "10/08"
              ax.text(0.16, -0.05, date_txt1, transform=ax.transAxes, fontsize=10, ha="center", va="center") # between CD
              ax.text(0.42, -0.05, date_txt2, transform=ax.transAxes, fontsize=10, ha="center", va="center") # between EF
              ax.text(0.71, -0.05, date_txt3, transform=ax.transAxes, fontsize=10, ha="center", va="center") # between GI
              ax.text(1.00, -0.05, date_txt4, transform=ax.transAxes, fontsize=10, ha="center", va="center") # after H

            # on-figure caption
            trigger_text = r"$\mu\tau_h$ Trigger : Tau Leg Performance"
            ax.text(0.01, 0.95, trigger_text, transform=ax.transAxes, fontsize=9)
            offline_text1  = r"$p_T^{\tau_h}$ (Offline) > 30 GeV"
            offline_text2  = "Offline Medium WP Tau ID Applied"
            ax.text(0.5, 0.40, offline_text1, transform=ax.transAxes, fontsize=12, ha="center", va="center")
            ax.text(0.5, 0.35, offline_text2, transform=ax.transAxes, fontsize=12, ha="center", va="center")
            # styling for DPS ("make it pretty")
            ax.set_title(r"2024 (13.6 TeV), 110 $fb^{-1}$", loc='right', y=0.99)
            ax.text(0.01, 1.02, "CMS", transform=ax.transAxes, fontsize=16, weight='bold')
            ax.text(0.12, 1.02, "Preliminary", transform=ax.transAxes, fontsize=16, style='italic')
            ax.text(0.12, 1.02, "Preliminary", transform=ax.transAxes, fontsize=16, style='italic')
            ax.minorticks_on()
            ax.tick_params(which="both", top=True, bottom=True, right=True, direction="in")
            ax.grid(axis='x', linestyle="dotted")
            ax.grid(axis='y', linestyle="dotted")

        plt.savefig("test" + ".png", format="png", dpi=150)
        #plt.savefig("test" + ".eps", format="eps")
        # plt.savefig(f"test_{channel}_{wp}_{dm}" + ".pdf", format="pdf")
        plt.show()



