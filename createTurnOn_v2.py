#!/usr/bin/env python

import argparse
from array import array
import math
import numpy as np
import os
import re
import sys
import ROOT
import cmsstyle as CMS

parser = argparse.ArgumentParser(description='Create turn on curves.')
parser.add_argument('--input', required=True, type=str, nargs='+', help="input files")
parser.add_argument('--output', required=True, type=str, help="output file prefix")
# parser.add_argument('--types-pnet', required=False, type=str, default='pnet_loose,pnet_medium,pnet_tight', help="pnet type")
# parser.add_argument('--types-deeptau', required=False, type=str, default=None, help="deeptau type")
parser.add_argument('--types-pnet', required=False, type=str, default='pnet_loose,pnet_medium,pnet_tight', help="pnet type")
parser.add_argument('--types-deeptau', required=False, type=str, default='deeptau', help="deeptau type")
parser.add_argument('--channels', required=False, type=str, default='mutau', help="channels to process")
parser.add_argument('--decay-modes', required=False, type=str, default='all,0,1,10,11,1011', help="decay modes to process")
parser.add_argument('--working-points', required=False, type=str,
                    default='Medium',
                    help="working points to process")
parser.add_argument('--var', required=False, type=str, default='leading_tau_pt', help="plot var")
parser.add_argument('--vbfditau_ptcut', required=False, type=int, default=30, help="VBF ditau pt cut on eta and phi")
parser.add_argument('--mutau_ptcut', required=False, type=int, default=30, help="mutau cut on eta and phi, PV_npvsGood")
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

# def CreateBins(max_pt, for_fitting):
#         #bins = np.arange(20, 100, step=10)
        
#         bins = np.arange(20, 40, step=4)
#         bins = np.append(bins, np.arange(40, 60, step=5))
#         bins = np.append(bins, np.arange(60, 100, step=10))
#         high_pt_bins = [ 100, 150, 200, 300, 400, 500, 650, 800, 1000 ]
        
#         # bins = np.arange(20, 100, step=8)
#         # high_pt_bins = [ 100, 200, 300, 400, 500, 650, 800, 1000 ]
        
#         n = 0
#         while n < len(high_pt_bins) and high_pt_bins[n] < max_pt:
#             n += 1
#         use_logx = max_pt > 200
#         return np.append(bins, high_pt_bins[0:n+1]), use_logx

def CreateBins(max_pt, for_fitting):
    if for_fitting:
        step=1
        return np.arange(20, 1000+step, step=step), False
        #high_pt_bins = np.arange(100, 501, step=5)
    else:
        #bins = np.arange(20, 100, step=10)
        # 10,20,24,28,32,36,40
        # bins = np.arange(20, 40, step=4)
        bins = [10,20,24,28,32,36,40]
        high_pt_bins = [ 40, 50, 70, 150, 200, 300, 400, 500, 650, 800, 1000 ]
        n = 0
        while n < len(high_pt_bins) and high_pt_bins[n] < max_pt:
            n += 1
        use_logx = max_pt > 200
        return np.append(bins, high_pt_bins[0:n+1]), use_logx

ptcut_dict = {
    # etau 30, ditau 35, mutau 30
    "mutau": 30, # pnet 27, deeptau 27
    "etau": 35,
    "ditau": 50, # pnet 30, deeptau 35,
    "ditaujet": 50, # pnet 26, deeptau 30
    "vbfditau": args.vbfditau_ptcut, # pnet 20, deeptau 20, extra cut 25, 30
    "vbfsingletau": 50, # pnet 45, deeptau 45
}
HLTptcut_dict = {
    "mutau": 27,
    "etau": 30,
    "ditau": 35,
    "ditaujet": 30,
    "vbfditau": 20,
    "vbfsingletau": 45,
    "singletau" : 180,
}
latex_dict = {
    "mutau" : "#mu#tau_{h}",
    "etau": "e#tau_{h}",
    "ditau": "Di-#tau_{h}",
    "ditaujet": "Di-#tau_{h}+Jet",
    "vbfditau": "VBF Di-#tau_{h}",
    "vbfsingletau": "VBF Single-#tau_{h}",
    "singletau": "Single-#tau_{h}",
}
labels_dict = {
    "pnet" : "ParticleNet",
    "pnet_nofilterbit" : "ParticleNet",
    "pnet_loose" : "ParticleNet",
    "pnet_loose_nofilterbit" : "ParticleNet",
    "pnet_loose_withmutaubit" : "ParticleNet",
    "pnet_medium" : "ParticleNet",
    "deeptau" : "DeepTau",
    "deeptau_nofilterbit" : "DeepTau",
    "deeptau_withmutaubit" : "DeepTau",
}

class TurnOnData:
    def __init__(self):
        self.hist_total = None
        self.hist_passed = None
        self.eff = None
        self.pass_yield = 0
        self.total_yield = 0

def CreateHistograms(input_file, channels, decay_modes, discr_name, working_points, hist_models, label, var,
                     output_file):
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
                # add tau pt cut for eta and phi
                if "pt" not in var:
                    df_wp = df_wp.Filter("leading_tau_pt>{}".format(ptcut_dict[channel]))
                # df_wp = df_wp.Filter('abs(sig_muon_eta) < 1.1')
                
                # if add L1 before, it would be (mu24 && L1 && HLT) / (mu24 && L1)
                # if "pnet" in label or "deeptau" in label:
                #     df_wp = df_wp.Filter("pass_{}_L1==1".format(channel))
                
                df_ch = df_wp.Filter('pass_{}_{} > 0.5 && pass_mutau_L1==1'.format(channel, label)) # df_wp is the total denominator
                # df_ch = df_wp.Filter('pass_{}_{} > 0.5'.format(channel, label)) # df_wp is the total denominator
                
                # if add L1 after, it would be (mu24 && HLT) / (mu24 && L1)
                # if label == "pnet" or label == "deeptau":
                #     df_wp = df_wp.Filter("pass_{}_L1==1".format(channel))
                
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
                    name_pattern = '{}_{}_{}{}_{}_{{}}'.format(label, channel, wp, dm_labels[dm], model_name)
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
labels = args.types_pnet.split(',')
# put deeptau in the last postion, take as denominator
# labels += args.types_deeptau.split(',')
# labels = args.types_deeptau.split(',')
colors_list = [ ROOT.kRed, ROOT.kGreen, ROOT.kViolet, ROOT.kBlue, ROOT.kCyan, ROOT.kOrange, ROOT.kGray, ROOT.kBlack]
colors = colors_list[0:len(labels)-1] + [colors_list[-1]]
# colors = [ROOT.kRed]
n_inputs = len(labels)
var = args.var
decay_modes = args.decay_modes.split(',')
channels = args.channels.split(',')
working_points = args.working_points.split(',')
if "pt" in var:
    if channels[0] == "singletau":
        bins = [50,90,120,150,200,250,500]
        use_logx = False
    else:
        bins, use_logx = CreateBins(150, False)
    # bins, use_logx = CreateBins(200, False)
    title, x_title = 'Offline #tau_{h} p_{T} [GeV]', 'Offline #tau_{h} p_{T} [GeV]'
elif "eta" in var:
    bins = np.linspace(-2.3, 2.3, 12)
    # bins = np.arange(-2.3, 2.3, step=0.2)
    # bins = np.append(bins, [2.3])
    use_logx = False
    title, x_title = 'Offline #tau_{h} #eta', 'Offline #tau_{h} #eta'
elif "phi" in var:
    bins = np.linspace(-3.2, 3.2, 12)
    # bins = np.arange(-3.2, 3.2, step=0.2)
    # bins = np.append(bins, [3.2])
    use_logx = False
    title, x_title = 'Offline #tau_{h} #phi', 'Offline #tau_{h} #phi'
elif "PV_npvs" in var:
    bins = np.linspace(0, 60, 7)
    bins = np.append(bins, [80])
    # bins = np.arange(-3.2, 3.2, step=0.2)
    # bins = np.append(bins, [3.2])
    use_logx = False
    if "Good" in var:
        title, x_title = 'Number of Offline Reconstructed Primary Vertices', 'Number of Offline Reconstructed Primary Vertices'
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
    print("Creating {} histograms...".format(labels[input_id]))
    turnOn_data[input_id] = CreateHistograms(input_files, channels, decay_modes, 'leading_tau_idDeepTauVSjet', # tau_idDeepTau2017v2p1VSjet,
                                             working_points, hist_models, labels[input_id], var, output_file)


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
            ratio_graph = [None] * n_inputs
            ref_hist = hist_models['plot'].GetHistogram()
            ratio_ref_hist = ref_hist.Clone()
            turnOns = [None] * n_inputs
            curves = [None] * n_inputs
            pass_yield = [None] * n_inputs
            total_yield = [None] * n_inputs
            # turnOn_data[0] -> pnet_loose
            # turnOn_data[1] -> pnet_medium
            # turnOn_data[2] -> pnet_tight
            # turnOn_data[3] -> deeptau
            for input_id in range(n_inputs):
                turnOns[input_id] = turnOn_data[input_id][dm][wp][channel]['plot']
                curves[input_id] = turnOns[input_id].eff
                pass_yield[input_id] = turnOns[input_id].pass_yield
                total_yield[input_id] = turnOns[input_id].total_yield
            y_min, y_max = (0, 1)
            y_title = 'L1+HLT Efficiency'
            # y_title = 'HLT Efficiency'
            # y_title = 'L1 Efficiency'
            title = '{} {}'.format(channel, dm_label) # title = '{} {}{}'.format(channel, wp, dm_label)
            # mutau_Tight
            plain_title = '{}_{}{}'.format(channel, wp, dm_plain_label)
            main_pad, ratio_pad, title_controls = RootPlotting.CreateTwoPadLayout(canvas, ref_hist, ratio_ref_hist,
                                                                                  log_x=use_logx, title=title)
            RootPlotting.ApplyAxisSetup(ref_hist, ratio_ref_hist, x_title=x_title, y_title=y_title,
                                        ratio_y_title='Ratio', y_range=(y_min, y_max * 1.1), max_ratio=1.5)
            # legend = RootPlotting.CreateLegend(pos=(0.25, 0.07), size=(0.2, 0.05*n_inputs))
            legend = RootPlotting.CreateLegend(text_size=0.05,pos=(0.43, 0.07), size=(0.2, 0.05*n_inputs))
            if "pt" not in var:
                txt = "p^{#tau}_{T}(Offline) " + "> {} GeV".format(ptcut_dict[channel])
                text = RootPlotting.DrawLabel(txt,text_size=0.05,pos=(0.55,0.37))
            else:
                txt = "p^{#tau}_{T}(HLT) " + "> {} GeV".format(HLTptcut_dict[channel])
                text2 = RootPlotting.DrawLabel(txt,text_size=0.05,pos=(0.55,0.37))
            text3 = RootPlotting.DrawLabel("Offline Medium WP Tau ID Applied",text_size=0.05,pos=(0.55,0.3))
            text4 = RootPlotting.DrawLabel("{} Trigger: Tau Performance".format(latex_dict[channel]),text_size=0.035,pos=(0.35,0.93))
            canvas.cd()
            
            cms_text = RootPlotting.DrawLabel("CMS",text_size=0.035,font=61,pos=(0.2,0.951))
            prelim_text = RootPlotting.DrawLabel("Preliminary",text_size=0.035,font=52,pos=(0.32,0.95))
            lumi_text = RootPlotting.DrawLabel("XXX fb^{-1}, 2025 (13.6 TeV)",text_size=0.035,font=42,pos=(0.75,0.95))
            canvas.Update()
            main_pad.cd()

            ratio_pad.cd()
            xmin = ratio_ref_hist.GetXaxis().GetXmin()
            xmax = ratio_ref_hist.GetXaxis().GetXmax()
            line = ROOT.TLine(xmin, 1.0, xmax, 1.0)
            line.SetLineColor(ROOT.kGray + 2)
            # line.SetLineStyle(2)
            line.SetLineWidth(3)
            line.Draw("same")
            main_pad.cd()
            
            for input_id in range(n_inputs):
                curve = curves[input_id]
                text_yield = pass_yield[input_id]
                try:
                    curve.Draw('SAME')
                except:
                    continue
                RootPlotting.ApplyDefaultLineStyle(curve, colors[input_id])
                legend.AddEntry(curve, labels[input_id] + '({})'.format(int(text_yield)), 'PLE')
                # legend.AddEntry(curve, labels_dict[labels[input_id]], 'PLE')

                if input_id < n_inputs - 1:
                    # turnOns[0], [1], [2] as numerator
                    # turnOns[-1] -> stands deeptau, take as denominator
                    ratio_graph[input_id] = RootPlotting.CreateEfficiencyRatioGraph(turnOns[input_id].hist_passed,
                                                                          turnOns[input_id].hist_total,
                                                                          turnOns[-1].hist_passed,
                                                                          turnOns[-1].hist_total)
                    if ratio_graph[input_id]:
                        output_file.WriteTObject(ratio_graph[input_id], 'ratio_{0}_{1}'.format(labels[input_id],plain_title), 'Overwrite')
                        ratio_pad.cd()
                        ratio_color = colors[input_id] if n_inputs > 2 else ROOT.kBlack
                        RootPlotting.ApplyDefaultLineStyle(ratio_graph[input_id], ratio_color)
                        ratio_graph[input_id].Draw("0PE SAME")
                        main_pad.cd()
            legend.Draw()

            canvas.Update()
            output_file.WriteTObject(canvas, 'canvas_{}'.format(plain_title), 'Overwrite')
            RootPlotting.PrintAndClear(canvas, args.output + '.pdf', plain_title, plot_id, n_plots, dm_plain_label,
                                       [ main_pad, ratio_pad ])
            plot_id += 1
output_file.Close()
