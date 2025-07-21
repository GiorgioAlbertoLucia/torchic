import numpy as np
from ROOT import TH2F, RooRealVar, TFile, TGraphErrors, TGraph
from torchic import HistLoadInfo
from torchic import histogram
from deprecated.torchic.core.roofitter import Roofitter

def purity_proton_TPC(output_file: TFile):

    h2_nsigmaTPC_info = HistLoadInfo('data/purity_AnalysisResults_PbPb.root', 'lithium4analysis/QA/h2NsigmaProtonTPC_preselection')
    h2_nsigmaTPC = histogram.load_hist(h2_nsigmaTPC_info)
    xmin = 0.2
    xmax = 0.8
    ymin = -6.
    ymax = 6.

    x = RooRealVar('x', 'x', ymin, ymax)
    roofitter = Roofitter(x, ['exp_offset', 'crystal_ball', 'pol1']) # fit with a gaussian and an exponential bkg
    print('params: \n', roofitter.pdf_params)
    roofitter.init_param('crystal_ball_1_mean', 0., -2.5, 2.5)
    roofitter.init_param('crystal_ball_1_sigma', 0., 1e-3, 1e3)
    #roofitter.init_param('exp_offset_0_sigma', 0., 1e-3, 1e3)
    roofitter.init_param('exp_0_alpha', 0.5, 0., 10.)
    roofitter.init_param('exp_0_offset', -0.5, -10., 0.)

    pts = []
    signal = []
    bkg = []
    tot = []
    purity = []

    out_dir = output_file.mkdir('Pr_TPC')

    for pt_bin in range(h2_nsigmaTPC.GetXaxis().FindBin(xmin), h2_nsigmaTPC.GetXaxis().FindBin(xmax)):

        h_nsigmaTPC = h2_nsigmaTPC.ProjectionY(f'h_nsigmaTPC_{pt_bin}', pt_bin, pt_bin, 'e')
        h_nsigmaTPC.Rebin(4)
        pt = h2_nsigmaTPC.GetXaxis().GetBinCenter(pt_bin)
        funcs_to_fit = ['pol1_2', 'crystal_ball_1'] if pt < 0.5 else ['exp_0', 'crystal_ball_1']

        fractions = roofitter.fit(h_nsigmaTPC, ymin, ymax, funcs_to_fit=funcs_to_fit)

        integrals = roofitter.functions_integral(-2, 2)
        
        pts.append(np.abs(pt))
        signal.append(integrals['crystal_ball_1'])
        ibkg = integrals['pol1_2'] if pt < 0.5 else integrals['exp_0']
        bkg.append(ibkg)
        itot = integrals['crystal_ball_1'] + ibkg
        tot.append(itot)
        purity.append(integrals['crystal_ball_1'] / itot)

        roofitter.plot(out_dir, canvas_name=f'cNSigmaTPC_pt{np.abs(pt)}', funcs_to_plot=funcs_to_fit)

    pt_bin_width = h2_nsigmaTPC.GetXaxis().GetBinWidth(1)/2.
    graph_signal = TGraphErrors(len(pts), np.array(pts, dtype=np.float32), np.array(signal, dtype=np.float32), np.array([pt_bin_width]*len(pts), dtype=np.float32), np.array([0.]*len(pts), dtype=np.float32))
    graph_bkg = TGraphErrors(len(pts), np.array(pts, dtype=np.float32), np.array(bkg, dtype=np.float32), np.array([pt_bin_width]*len(pts), dtype=np.float32), np.array([0.]*len(pts), dtype=np.float32))
    graph_tot = TGraphErrors(len(pts), np.array(pts, dtype=np.float32), np.array(tot, dtype=np.float32), np.array([pt_bin_width]*len(pts), dtype=np.float32), np.array([0.]*len(pts), dtype=np.float32))
    graph_purity = TGraphErrors(len(pts), np.array(pts, dtype=np.float32), np.array(purity, dtype=np.float32), np.array([pt_bin_width]*len(pts), dtype=np.float32), np.array([0.]*len(pts), dtype=np.float32))

    out_dir.cd()
    graph_signal.Write('g_signal')
    graph_bkg.Write('g_bkg')
    graph_tot.Write('g_tot')
    graph_purity.Write('g_purity')
    h2_nsigmaTPC.Write()


def purity_proton_TOF(output_file: TFile):
    
    h2_nsigmaTOF_info = HistLoadInfo('data/data_visual_selectionsPr_purity.root', 'TOF/NSigmaTOFvsPtPr')
    h2_nsigmaTOF = histogram.load_hist(h2_nsigmaTOF_info)
    xmin = -2
    xmax = -0.7
    ymin = -8.
    ymax = 8.

    x = RooRealVar('x', 'x', ymin, ymax)
    roofitter = Roofitter(x, ['exp_offset', 'gaus', 'pol1']) # fit with a gaussian and an exponential bkg
    print('params: \n', roofitter.pdf_params)
    roofitter.init_param('gaus_1_mean', 0., -2.5, 2.5)
    roofitter.init_param('gaus_1_sigma', 0., 1e-3, 1e3)
    #roofitter.init_param('exp_offset_0_sigma', 0., 1e-3, 1e3)
    roofitter.init_param('exp_0_alpha', 0.5, 0., 10.)
    roofitter.init_param('exp_0_offset', -0.5, -10., 0.)

    pts = []
    signal = []
    bkg = []    
    tot = []
    purity = []

    out_dir = output_file.mkdir('Pr_TOF')

    for pt_bin in range(h2_nsigmaTOF.GetXaxis().FindBin(xmin), h2_nsigmaTOF.GetXaxis().FindBin(xmax)):

        h_nsigmaTOF = h2_nsigmaTOF.ProjectionY(f'h_nsigmaTOF_{pt_bin}', pt_bin, pt_bin, 'e')
        h_nsigmaTOF.Rebin(4)
        pt = h2_nsigmaTOF.GetXaxis().GetBinCenter(pt_bin)
        funcs_to_fit = ['exp_0', 'gaus_1'] if pt < -1.5 else ['pol1_2', 'gaus_1']

        fractions = roofitter.fit(h_nsigmaTOF, ymin, ymax, funcs_to_fit=funcs_to_fit)

        integrals = roofitter.functions_integral(-2, 2)
        
        pts.append(np.abs(pt))
        signal.append(integrals['gaus_1'])
        ibkg = integrals['exp_0'] if pt < -1.5 else integrals['pol1_2']
        bkg.append(ibkg)
        itot = integrals['gaus_1'] + ibkg
        tot.append(itot)
        purity.append(integrals['gaus_1'] / itot)

        roofitter.plot(out_dir, canvas_name=f'cNSigmaTOF_pt{np.abs(pt)}', funcs_to_plot=funcs_to_fit)

    pt_bin_width = h2_nsigmaTOF.GetXaxis().GetBinWidth(1)/2.
    graph_signal = TGraphErrors(len(pts), np.array(pts, dtype=np.float32), np.array(signal, dtype=np.float32), np.array([pt_bin_width]*len(pts), dtype=np.float32), np.array([0.]*len(pts), dtype=np.float32))
    graph_bkg = TGraphErrors(len(pts), np.array(pts, dtype=np.float32), np.array(bkg, dtype=np.float32), np.array([pt_bin_width]*len(pts), dtype=np.float32), np.array([0.]*len(pts), dtype=np.float32))
    graph_tot = TGraphErrors(len(pts), np.array(pts, dtype=np.float32), np.array(tot, dtype=np.float32), np.array([pt_bin_width]*len(pts), dtype=np.float32), np.array([0.]*len(pts), dtype=np.float32))
    graph_purity = TGraphErrors(len(pts), np.array(pts, dtype=np.float32), np.array(purity, dtype=np.float32), np.array([pt_bin_width]*len(pts), dtype=np.float32), np.array([0.]*len(pts), dtype=np.float32))

    out_dir.cd()
    graph_signal.Write('g_signal')
    graph_bkg.Write('g_bkg')
    graph_tot.Write('g_tot')
    graph_purity.Write('g_purity')
    h2_nsigmaTOF.Write()

def purity_he3_TPC(output_file: TFile):

    h2_nsigmaTPC_info = HistLoadInfo('data/purity_AnalysisResults_PbPb.root', 'lithium4analysis/QA/h2NsigmaHe3TPC_preselection')
    h2_nsigmaTPC = histogram.load_hist(h2_nsigmaTPC_info)
    xmin = 0.8
    xmax = 2.
    ymin = -5.
    ymax = 5.

    x = RooRealVar('x', 'x', ymin, ymax)
    roofitter = Roofitter(x, ['exp_offset', 'gaus', 'gaus']) # fit with a gaussian and an exponential bkg
    print('params: \n', roofitter.pdf_params)
    roofitter.init_param('gaus_1_mean', 0.2, 0., 1.5)
    roofitter.init_param('gaus_1_sigma', 0.2, 1e-3, 0.8)
    roofitter.init_param('exp_0_alpha', 0.5, 0., 10.)
    roofitter.init_param('exp_0_offset', -0.5, -10., 0.)
    roofitter.init_param('gaus_2_mean', -4., -3.2, -7.)
    roofitter.init_param('gaus_2_sigma', 0.2, 1e-3, 0.5)
    #roofitter.init_param('gaus_2_tau', -0.5, -10., 0.)

    pts = []
    signal = []
    bkg = []
    tot = []
    purity = []

    out_dir = output_file.mkdir('He_TPC')

    for pt_bin in range(h2_nsigmaTPC.GetXaxis().FindBin(xmin), h2_nsigmaTPC.GetXaxis().FindBin(xmax)):

        h_nsigmaTPC = h2_nsigmaTPC.ProjectionY(f'h_nsigmaTPC_{pt_bin}', pt_bin, pt_bin, 'e')
        h_nsigmaTPC.Rebin(4)
        pt = h2_nsigmaTPC.GetXaxis().GetBinCenter(pt_bin)
        funcs_to_fit = ['gaus_2', 'gaus_1'] if pt < 1.4 else ['exp_0', 'gaus_1']

        mean = histogram.get_mean(h_nsigmaTPC, -6., -2.)
        rms = histogram.get_rms(h_nsigmaTPC, -6., -2.)
        roofitter.init_param('gaus_2_mean', mean, mean-2*rms, mean+2*rms)
        roofitter.init_param('gaus_2_sigma', rms, 1e-3, 0.5)

        fractions = roofitter.fit(h_nsigmaTPC, ymin, ymax, funcs_to_fit=funcs_to_fit)

        integrals = roofitter.functions_integral(-2, 2)
        
        pts.append(np.abs(pt))
        signal.append(integrals['gaus_1'])
        ibkg = integrals['gaus_2'] if pt < 1.4 else integrals['exp_0']
        bkg.append(ibkg)
        itot = integrals['gaus_1'] + ibkg
        tot.append(itot)
        purity.append(integrals['gaus_1'] / itot)

        roofitter.plot(out_dir, canvas_name=f'cNSigmaTPC_pt{np.abs(pt)}', funcs_to_plot=funcs_to_fit)

    pt_bin_width = h2_nsigmaTPC.GetXaxis().GetBinWidth(1)/2.
    graph_signal = TGraphErrors(len(pts), np.array(pts, dtype=np.float32), np.array(signal, dtype=np.float32), np.array([pt_bin_width]*len(pts), dtype=np.float32), np.array([0.]*len(pts), dtype=np.float32))
    graph_bkg = TGraphErrors(len(pts), np.array(pts, dtype=np.float32), np.array(bkg, dtype=np.float32), np.array([pt_bin_width]*len(pts), dtype=np.float32), np.array([0.]*len(pts), dtype=np.float32))
    graph_tot = TGraphErrors(len(pts), np.array(pts, dtype=np.float32), np.array(tot, dtype=np.float32), np.array([pt_bin_width]*len(pts), dtype=np.float32), np.array([0.]*len(pts), dtype=np.float32))
    graph_purity = TGraphErrors(len(pts), np.array(pts, dtype=np.float32), np.array(purity, dtype=np.float32), np.array([pt_bin_width]*len(pts), dtype=np.float32), np.array([0.]*len(pts), dtype=np.float32))

    out_dir.cd()
    graph_signal.Write('g_signal')
    graph_bkg.Write('g_bkg')
    graph_tot.Write('g_tot')
    graph_purity.Write('g_purity')
    h2_nsigmaTPC.Write()

if __name__ == '__main__':
    
    output_file = TFile('data/purity.root', 'recreate')
    purity_proton_TPC(output_file)
    purity_proton_TOF(output_file)
    purity_he3_TPC(output_file)
    output_file.Close()