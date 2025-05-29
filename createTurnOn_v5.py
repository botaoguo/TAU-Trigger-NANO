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
parser.add_argument('--mutau_ptcut', required=False, type=int, default=27, help="mutau cut on eta and phi, PV_npvsGood")
parser.add_argument('--etau_ptcut', required=False, type=int, default=30, help="etau cut on eta and phi, PV_npvsGood")
parser.add_argument('--ditau_ptcut', required=False, type=int, default=35, help="ditaujet cut on eta and phi, PV_npvsGood")
parser.add_argument('--ditaujet_ptcut', required=False, type=int, default=30, help="ditaujet cut on eta and phi, PV_npvsGood")
parser.add_argument('--vbfsingletau_ptcut', required=False, type=int, default=45, help="VBF ditau pt cut on eta and phi")
parser.add_argument('--vbfditau_ptcut', required=False, type=int, default=30, help="VBF ditau pt cut on eta and phi")

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
    "ditau": 30,
    "ditaujet": 26,
    "vbfditau": 20,
    "vbfsingletau": 45,
    "singletau" : 130,
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
                # add tau pt cut for eta and phi
                if "pt" not in var:
                    df_wp = df_wp.Filter("leading_tau_pt>{}".format(ptcut_dict[channel]))
                # OSminusSS, OSplusSS, OSonly
                if "OSonly" in era:
                    df_wp = df_wp.Filter("sig_muon_charge != leading_tau_charge")
                # df_wp = df_wp.Filter('run <= 382262')
                # df_wp = df_wp.Filter('abs(sig_muon_eta) < 1.1')
                if "mc" in input_file:
                    df_ch = df_wp.Filter('pass_{}_{} > 0.5 && genmuon_idx!=-1 && gentau_idx!=-1'.format(channel, label))
                else:
                    df_ch = df_wp.Filter('pass_{}_{} > 0.5'.format(channel, label))                    
                for model_name, hist_model in hist_models.items():
                    turn_on = TurnOnData()
                    if "OSplusSS" in era:
                        df_wp = df_wp.Define("abs_weight", "abs(weight)")
                        df_ch = df_ch.Define("abs_weight", "abs(weight)")
                        turn_on.hist_total = df_wp.Histo1D(hist_model, var, 'abs_weight')
                        turn_on.hist_passed = df_ch.Histo1D(hist_model, var, 'abs_weight')
                    else:
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
    title, x_title = 'Offline #eta^{#tau}', 'Offline #eta^{#tau}'
elif "phi" in var:
    bins = np.linspace(-3.2, 3.2, 12)
    # bins = np.arange(-3.2, 3.2, step=0.2)
    # bins = np.append(bins, [3.2])
    use_logx = False
    title, x_title = 'Offline #phi^{#tau}', 'Offline #phi^{#tau}'
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
    fileprefix = os.path.basename(input_files[int(input_id)])
    fileprefix = os.path.splitext(fileprefix)[0]
    era = fileprefix.strip("skimtuple_")
    print("Creating {0} histograms for {1}...".format(labels[input_id],era))
    turnOn_data[input_id] = CreateHistograms(input_files[int(input_id)], channels, decay_modes, 'leading_tau_idDeepTauVSjet', # tau_idDeepTau2017v2p1VSjet,
                                             working_points, hist_models, labels[input_id], var, output_file, era)


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
            title = '{} {} {} {}'.format(channel, wp, labels[0], dm_label)
            # title = '{} {}{}'.format(channel, wp, dm_label)
            # title = '{} {} run <= 382262'.format(channel, dm_label)
            # mutau_Tight
            plain_title = '{}_{}{}'.format(channel, wp, dm_plain_label)
            title_controls = RootPlotting.CreateNoPadLayout(canvas, ref_hist, ratio_ref_hist,
                                                                                  log_x=use_logx, title=title, x_title=x_title, y_title=y_title)
            # RootPlotting.ApplyAxisSetup2(ref_hist, ratio_ref_hist, x_title=x_title, y_title=y_title,
            #                             ratio_y_title='Ratio', y_range=(y_min, y_max * 1.1), max_ratio=1.5)
            legend = RootPlotting.CreateLegend(text_size=0.035, pos=(0.37, 0.27), size=(0.2, 0.05*n_inputs))

            if "pt" not in var:
                txt = "p^{#tau}_{T}(Offline) " + "> {} GeV".format(ptcut_dict[channel])
                text = RootPlotting.DrawLabel(txt,text_size=0.035,pos=(0.55,0.45))
            else:
                txt = "p^{#tau}_{T}(HLT) " + "> {} GeV".format(HLTptcut_dict[channel])
                text2 = RootPlotting.DrawLabel(txt,text_size=0.035,pos=(0.55,0.45))
            text3 = RootPlotting.DrawLabel("Offline Medium WP Tau ID Applied",text_size=0.035,pos=(0.55,0.4))
            if channels[0] == "vbfsingletau":
                text4 = RootPlotting.DrawLabel("{} Trigger: Tau Performance".format(latex_dict[channel]),text_size=0.032,pos=(0.4,0.86))
            else:
                text4 = RootPlotting.DrawLabel("{} Trigger: Tau Performance".format(latex_dict[channel]),text_size=0.032,pos=(0.35,0.86))
            canvas.cd()
            
            cms_text = RootPlotting.DrawLabel("CMS",text_size=0.035,font=61,pos=(0.14,0.921))
            prelim_text = RootPlotting.DrawLabel("Preliminary",text_size=0.035,font=52,pos=(0.26,0.92))
            lumi_text = RootPlotting.DrawLabel("109 fb^{-1}, 2024 (13.6 TeV)",text_size=0.035,font=42,pos=(0.72,0.92))
            canvas.Update()

            for input_id in range(n_inputs):
                curve = curves[input_id]
                text_yield = pass_yield[input_id]
                # define the era
                fileprefix = os.path.basename(input_files[int(input_id)])
                fileprefix = os.path.splitext(fileprefix)[0]
                era = fileprefix.strip("skimtuple_2024_")
                try:
                    curve.Draw('SAME')
                except:
                    continue
                RootPlotting.ApplyDefaultLineStyle(curve, colors[input_id])
                # RootPlotting.ApplyDefaultLineStyleMarker(curve, colors[input_id], marker[int(input_id)])
                # legend.AddEntry(curve, labels[input_id] + " {} ".format(era) + "({})".format(int(text_yield)), 'PLE')
                legend.AddEntry(curve, labels_dict[labels[input_id]] + " " + era, 'PLE')

                if input_id < n_inputs - 1:
                    # turnOns[0], [1], [2] as numerator
                    # turnOns[-1] -> stands deeptau, take as denominator
                    ratio_graph[input_id] = RootPlotting.CreateEfficiencyRatioGraph(turnOns[input_id].hist_passed,
                                                                          turnOns[input_id].hist_total,
                                                                          turnOns[-1].hist_passed,
                                                                          turnOns[-1].hist_total)
            legend.Draw()

            canvas.Update()
            output_file.WriteTObject(canvas, 'canvas_{}'.format(plain_title), 'Overwrite')
            RootPlotting.PrintAndClear(canvas, args.output + '.pdf', plain_title, plot_id, n_plots, dm_plain_label,
                                       [ ])
            plot_id += 1
output_file.Close()
