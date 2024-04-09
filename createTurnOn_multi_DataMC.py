#!/usr/bin/env python

import argparse
from array import array
import math
import numpy as np
import os
import re
import sys
import ROOT

parser = argparse.ArgumentParser(description='Create turn on curves.')
parser.add_argument('--input-data', required=True, type=str, help="skimmed data input")
parser.add_argument('--input-dy-mc', required=True, type=str, help="skimmed DY MC input")
parser.add_argument('--output', required=True, type=str, help="output file prefix")
parser.add_argument('--channels', required=False, type=str, default='etau,mutau,ditau,ditaujet', help="channels to process")
# parser.add_argument('--channels', required=False, type=str, default='ditau,ditau_withptiso_nobitcut,ditau_withptiso_bit1,ditau_withptiso_bit1_bit17,ditau_withptiso_bit1_bit17_0bit18', help="channels to process")
parser.add_argument('--decay-modes', required=False, type=str, default='all,0,1,10,11,1011', help="decay modes to process")
parser.add_argument('--working-points', required=False, type=str,
                    default='VVVLoose,VVLoose,VLoose,Loose,Medium,Tight,VTight,VVTight',
                    help="working points to process")
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

bin_scans = {
    2:  [ 0.01 ],
    5:  [ 0.01, 0.05 ],
    10: [ 0.05, 0.1 ],
    20: [ 0.1, 0.2 ],
    #50: [ 0.2, 0.4 ],
    100: [ 0.2 ],
}
bin_scan_pairs = []
for max_bin_delta_pt, max_rel_err_vec in bin_scans.items():
    for max_rel_err in max_rel_err_vec:
        bin_scan_pairs.append([max_bin_delta_pt, max_rel_err])
# debug for DM 11
bin_scans_mc = {
    # 2:  [ 0.01 ],
    # 5:  [ 0.01, 0.05 ],
    # 10: [ 0.05, 0.1 ],
    # 20: [ 0.1, 0.2 ],
    #50: [ 0.2, 0.4 ],
    100: [ 0.08 ],
}
bin_scan_pairs_mc = []

# mc 100: [0.08,0.1]
for max_bin_delta_pt_mc, max_rel_err_vec_mc in bin_scans_mc.items():
    for max_rel_err_mc in max_rel_err_vec_mc:
        bin_scan_pairs_mc.append([max_bin_delta_pt_mc, max_rel_err_mc])

def CreateBins(max_pt, for_fitting):
    if for_fitting:
        step=1
        return np.arange(20, 1000+step, step=step), False
        #high_pt_bins = np.arange(100, 501, step=5)
    else:
        #bins = np.arange(20, 100, step=10)
        bins = np.arange(20, 40, step=4)
        bins = np.append(bins, np.arange(40, 60, step=5))
        bins = np.append(bins, np.arange(60, 100, step=10))
        high_pt_bins = [ 100, 150, 200, 300, 400, 500, 650, 800, 1000 ]
        n = 0
        while n < len(high_pt_bins) and high_pt_bins[n] < max_pt:
            n += 1
        use_logx = max_pt > 200
        return np.append(bins, high_pt_bins[0:n+1]), use_logx

class TurnOnData:
    def __init__(self):
        self.hist_total = None
        self.hist_passed = None
        self.eff = None

def CreateHistograms(input_file, input_file_mc, channels, decay_modes, discr_name, working_points, hist_models, label, label_mc, var,
                     output_file):
    df = ROOT.RDataFrame('Events', input_file)
    df_mc = ROOT.RDataFrame('Events', input_file_mc)
    turnOn_data = {}
    turnOn_mc = {}
    dm_labels = {}

    for dm in decay_modes:
        if dm == 'all':
            dm_labels[dm] = ''
            df_dm = df
            df_dm_mc = df_mc
        elif dm == '1011':
            dm_labels[dm] = '_dm{}'.format(dm)
            df_dm = df.Filter('tau_decayMode == 10 | tau_decayMode == 11')
            df_dm_mc = df_mc.Filter('tau_decayMode == 10 | tau_decayMode == 11')
        else:
            dm_labels[dm] = '_dm{}'.format(dm)
            df_dm = df.Filter('tau_decayMode == {}'.format(dm))
            df_dm_mc = df_mc.Filter('tau_decayMode == {}'.format(dm))
        turnOn_data[dm] = {}
        turnOn_mc[dm] = {}
        for wp in working_points:
            wp_bit = ParseEnum(DiscriminatorWP, wp)
            # df_wp = df_dm.Filter('({} & (1 << {})) != 0'.format(discr_name, wp_bit))
            df_wp = df_dm.Filter('{0} >= {1}'.format(discr_name, wp_bit))
            df_wp_mc = df_dm_mc.Filter('{0} >= {1}'.format(discr_name, wp_bit))
            turnOn_data[dm][wp] = {}
            turnOn_mc[dm][wp] = {}
            for channel in channels:
                turnOn_data[dm][wp][channel] = {}
                turnOn_mc[dm][wp][channel] = {}
                df_ch = df_wp.Filter('pass_{} > 0.5'.format(channel))
                df_ch_mc = df_wp_mc.Filter('pass_{} > 0.5'.format(channel))
                for model_name, hist_model in hist_models.items():
                    turn_on = TurnOnData() # class TurnOnData
                    turn_on_mc = TurnOnData() # class TurnOnData
                    
                    turn_on.hist_total = df_wp.Histo1D(hist_model, var, 'weight')
                    turn_on_mc.hist_total = df_wp_mc.Histo1D(hist_model, var, 'weight')
                    
                    turn_on.hist_passed = df_ch.Histo1D(hist_model, var, 'weight')
                    turn_on_mc.hist_passed = df_ch_mc.Histo1D(hist_model, var, 'weight')
                    
                    turnOn_data[dm][wp][channel][model_name] = turn_on
                    turnOn_mc[dm][wp][channel][model_name] = turn_on_mc

    for dm in decay_modes:
        for wp in working_points:
            for channel in channels:
                for model_name, hist_model in hist_models.items():
                    
                    turn_on = turnOn_data[dm][wp][channel][model_name]
                    turn_on_mc = turnOn_mc[dm][wp][channel][model_name]
                    
                    name_pattern = '{}_{}_{}{}_{}_{{}}'.format(label, channel, wp, dm_labels[dm], model_name)
                    name_pattern_mc = '{}_{}_{}{}_{}_{{}}'.format(label_mc, channel, wp, dm_labels[dm], model_name)
                    turn_on.name_pattern = name_pattern
                    turn_on.name_pattern_mc = name_pattern_mc
                    if 'fit' in model_name:
                        passed, total, eff, passed_mc, total_mc, eff_mc = AutoRebinAndEfficiencyforDataMCboth(turn_on.hist_passed.GetPtr(),
                                                                                                      turn_on.hist_total.GetPtr(), bin_scan_pairs, turn_on_mc.hist_passed.GetPtr(), turn_on_mc.hist_total.GetPtr())
                    else:
                        passed, total = turn_on.hist_passed.GetPtr(), turn_on.hist_total.GetPtr()
                        passed_mc, total_mc = turn_on_mc.hist_passed.GetPtr(), turn_on_mc.hist_total.GetPtr()
                        # # new add by botao
                        # if (passed.Integral() < 0) or (total.Integral() < 0):
                        #     continue
                        # # end add
                        try:
                            FixEfficiencyBins(passed, total)
                            FixEfficiencyBins(passed_mc, total_mc)
                        except:
                            continue
                        turn_on.eff = ROOT.TEfficiency(passed, total)
                        turn_on_mc.eff = ROOT.TEfficiency(passed_mc, total_mc)
                        eff = turn_on.eff
                        eff_mc = turn_on_mc.eff
                    output_file.WriteTObject(total, name_pattern.format('total'), 'Overwrite')
                    output_file.WriteTObject(passed, name_pattern.format('passed'), 'Overwrite')
                    output_file.WriteTObject(eff, name_pattern.format('eff'), 'Overwrite')
                    
                    output_file.WriteTObject(total_mc, name_pattern_mc.format('total'), 'Overwrite')
                    output_file.WriteTObject(passed_mc, name_pattern_mc.format('passed'), 'Overwrite')
                    output_file.WriteTObject(eff_mc, name_pattern_mc.format('eff'), 'Overwrite')
                    # print(name_pattern)
                    # print('hist_total {}'.format(turn_on.hist_total.GetPtr().GetNbinsX()))
                    # for n in range(turn_on.hist_total.GetPtr().GetNbinsX() + 1):
                    #     print('\t{} {} +/- {} {} +/- {}'.format(turn_on.hist_total.GetPtr().GetBinLowEdge(n+1), turn_on.hist_total.GetPtr().GetBinContent(n+1), turn_on.hist_total.GetPtr().GetBinError(n+1), turn_on.hist_passed.GetPtr().GetBinContent(n+1), turn_on.hist_passed.GetPtr().GetBinError(n+1)))
                    # print('hist_passed {}'.format(turn_on.hist_passed.GetPtr().GetNbinsX()))
                    # for n in range(turn_on.hist_passed.GetPtr().GetNbinsX() + 1):
                    #     print('\t{} {}'.format(turn_on.hist_passed.GetPtr().GetBinLowEdge(n+1), ))
                    #if 'fit' not in model_name:
                    #    turn_on.eff = ROOT.TEfficiency(turn_on.hist_passed.GetPtr(), turn_on.hist_total.GetPtr())
                    #    output_file.WriteTObject(turn_on.eff, name_pattern.format('passed'), 'Overwrite')

    return turnOn_data,turnOn_mc

output_file = ROOT.TFile(args.output + '.root', 'RECREATE')
input_files = [ args.input_data, args.input_dy_mc ]
n_inputs = len(input_files)
labels = [ 'data', 'mc' ]
var = 'tau_pt'
title, x_title = '#tau p_{T}', '#tau p_{T} (GeV)'
decay_modes = args.decay_modes.split(',')
channels = args.channels.split(',')
working_points = args.working_points.split(',')
bins, use_logx = CreateBins(200, False)
bins_fit, _ = CreateBins(200, True)
hist_models = {
    'plot': ROOT.RDF.TH1DModel(var, var, len(bins) - 1, array('d', bins)),
    'fit': ROOT.RDF.TH1DModel(var, var, len(bins_fit) - 1, array('d', bins_fit))
}
turnOn_data = [None] * n_inputs
# input id = 0 -> data
# input id = 1 -> mc
# for input_id in range(n_inputs):
#     print("Creating {} histograms...".format(labels[input_id]))
#     turnOn_data[input_id] = CreateHistograms(input_files[input_id], channels, decay_modes, 'tau_idDeepTau2018v2p5VSjet', # tau_idDeepTau2017v2p1VSjet,
#                                              working_points, hist_models, labels[input_id], var, output_file)
turnOn_data[0], turnOn_data[1] = CreateHistograms(input_files[0], input_files[1], channels, decay_modes, 'tau_idDeepTau2018v2p5VSjet',
                                             working_points, hist_models, labels[0], labels[1], var, output_file)

colors = [ ROOT.kRed, ROOT.kBlack ]
canvas = RootPlotting.CreateCanvas()

n_plots = len(decay_modes) * len(channels) * len(working_points)
plot_id = 0
for channel in channels:
    for wp in working_points:
        for dm in decay_modes:
            if dm == 'all':
                dm_label = ''
                dm_plain_label = ''
            else:
                dm_label = ' DM={}'.format(dm)
                dm_plain_label = '_dm{}'.format(dm)
            ratio_graph = None
            ref_hist = hist_models['plot'].GetHistogram()
            ratio_ref_hist = ref_hist.Clone()
            turnOns = [None] * n_inputs
            curves = [None] * n_inputs
            for input_id in range(n_inputs):
                turnOns[input_id] = turnOn_data[input_id][dm][wp][channel]['plot']
                curves[input_id] = turnOns[input_id].eff
            y_min, y_max = (0, 1)
            y_title = 'Efficiency'
            title = '{} {}{}'.format(channel, wp, dm_label)
            plain_title = '{}_{}{}'.format(channel, wp, dm_plain_label)
            main_pad, ratio_pad, title_controls = RootPlotting.CreateTwoPadLayout(canvas, ref_hist, ratio_ref_hist,
                                                                                  log_x=use_logx, title=title)
            RootPlotting.ApplyAxisSetup(ref_hist, ratio_ref_hist, x_title=x_title, y_title=y_title,
                                        ratio_y_title='Ratio', y_range=(y_min, y_max * 1.1), max_ratio=1.5)
            legend = RootPlotting.CreateLegend(pos=(0.78, 0.28), size=(0.2, 0.15))
            for input_id in range(n_inputs):
                curve = curves[input_id]
                try:
                    curve.Draw('SAME')
                except:
                    continue
                RootPlotting.ApplyDefaultLineStyle(curve, colors[input_id])
                legend.AddEntry(curve, labels[input_id], 'PLE')

                if input_id < n_inputs - 1:
                    ratio_graph = RootPlotting.CreateEfficiencyRatioGraph(turnOns[input_id].hist_passed,
                                                                          turnOns[input_id].hist_total,
                                                                          turnOns[-1].hist_passed,
                                                                          turnOns[-1].hist_total)
                    if ratio_graph:
                        output_file.WriteTObject(ratio_graph, 'ratio_{}'.format(plain_title), 'Overwrite')
                        ratio_pad.cd()
                        ratio_color = colors[input_id] if n_inputs > 2 else ROOT.kBlack
                        RootPlotting.ApplyDefaultLineStyle(ratio_graph, ratio_color)
                        ratio_graph.Draw("0PE SAME")
                        main_pad.cd()
            legend.Draw()

            canvas.Update()
            output_file.WriteTObject(canvas, 'canvas_{}'.format(plain_title), 'Overwrite')
            RootPlotting.PrintAndClear(canvas, args.output + '.pdf', plain_title, plot_id, n_plots,
                                       [ main_pad, ratio_pad ])
            plot_id += 1
output_file.Close()
