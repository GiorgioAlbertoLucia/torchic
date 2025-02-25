import os
import numpy as np
from math import erf
import pandas as pd
from copy import deepcopy
from ROOT import TH2F, TGraphErrors, TDirectory, TF1, gInterpreter, TCanvas, TPaveText
from ROOT import RooDataHist, RooArgList, RooRealVar, RooFit

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BETHEBLOCH_DIR = os.path.join(CURRENT_DIR, 'BetheBloch.hh')
gInterpreter.ProcessLine(f'#include "{BETHEBLOCH_DIR}"')
from ROOT import BetheBloch

from torchic.core.roofitter import Roofitter
from torchic.utils.terminal_colors import TerminalColors as tc

DEFAULT_BETHEBLOCH_PARS = { # params for TPC He3 pp
                            'kp1': -241.490, 
                            'kp2': 0.374245,
                            'kp3': 1.397847,
                            'kp4': 1.078250,
                            'kp5': 2.048336
                          }

def py_BetheBloch(betagamma, kp1, kp2, kp3, kp4, kp5):
    '''
        Python implementation of the Bethe-Bloch formula.
    '''
    beta = betagamma / np.sqrt(1 + betagamma**2)
    aa = beta**kp4
    bb = (1/betagamma)**kp5
    bb = np.log(bb + kp3)
    return (kp2 - aa - bb) * kp1 / aa

def cluster_size_parametrisation(betagamma, kp1, kp2, kp3, charge, kp4):
    '''
        Python implementation of a simil Bethe-Bloch formula: kp1 / betagamma**kp2 + kp3
    '''
    return (kp1 / betagamma**kp2 + kp3) * charge ** kp4

def cluster_size_resolution(betagamma, rp0, rp1, rp2):
    '''
        Python implementation of the resolution function.
    '''
    return rp0 * erf((betagamma - rp1) / rp2)
np_cluster_size_resolution = np.vectorize(cluster_size_resolution)

def bethe_bloch_calibration(x: RooRealVar, h2: TH2F, output_file: TDirectory, fit_func, fit_params, **kwargs) -> dict:
    '''
        Perform a Bethe-Bloch calibration on a 2D histogram.
        The histogram is sliced along the x-axis and fitted with a Gaussian.
        The mean and sigma of the Gaussian are stored in a TGraphErrors.
        The bin error is calculated as the bin width.
        The mean error is calculated as mean * expected resolution.
        The histogram, curve and TGraphErrors are stored in the output file.

        IMPORTANT REQUIREMENT: if not provided, the names of the parameters of the gaussian must be 'mean', 'sigma'

        Parameters:
        - h2: TH2F
            The 2D histogram to be calibrated.
        - output_file: TDirectory
            The output file where the TGraphErrors will be stored.
        - **kwargs:
            Additional arguments to be passed to the fit_by_slices function.
            -> first_bin_fit_by_slices: int
                First bin to be fitted by slices.
            -> last_bin_fit_by_slices: int
                Last bin to be fitted by slices.
            -> mean_label: str
                The name of the mean parameter of the gaussian.
            -> sigma_label: str
                The name of the sigma parameter of the gaussian.
            -> f'{signal_func_name}_sigma_err': str
                The name of the error of the sigma parameter of the gaussian.
    '''

    fit_results = pd.DataFrame()
    first_fit_bin = kwargs.get('first_bin_fit_by_slices')
    for xbin in range(kwargs.get('first_bin_fit_by_slices'), kwargs.get('last_bin_fit_by_slices')+1):
        h = h2.ProjectionY(f'h_{xbin}', xbin, xbin, 'e')
        bin_center = h2.GetXaxis().GetBinCenter(xbin)
            
        datahist =  RooDataHist('data', 'data', RooArgList(x), h)

        fit_ranges = kwargs.get('fit_ranges', None)
        if fit_ranges is not None:
            fit_range = fit_ranges[xbin-first_fit_bin]
            fit_func.fitTo(datahist, RooFit.Save(), RooFit.Range(fit_range[0], fit_range[1]))
        else:
            fit_func.fitTo(datahist, RooFit.Save())

        frame = x.frame()
        datahist.plotOn(frame)
        fit_func.plotOn(frame)

        ## put fit parameters on the plot
        text = TPaveText(0.7, 0.7, 0.9, 0.9, 'NDC')
        for name, param in fit_params.items():
            text.AddText(f'{name} = {param.getVal():.4f} +/- {param.getError():.4f}')
        frame.addObject(text)
        frame.Draw()
        output_file.cd()
        frame.Write(f'frame_{xbin}')

        bin_fit_results = pd.DataFrame.from_dict({key: [value.getVal()] for key, value in fit_params.items()})
        bin_fit_results_error = pd.DataFrame.from_dict({key+'_err': [value.getError()] for key, value in fit_params.items()})
        bin_fit_results['bin_center'] = bin_center
        bin_fit_results['integral'] = fit_func.createIntegral(RooArgList(x)).getVal()
        bin_fit_results['unnorm_integral'] = fit_func.createIntegral(RooArgList(x)).getVal() * h.Integral()
        bin_fit_results = pd.concat([bin_fit_results, bin_fit_results_error], axis=1)
        fit_results = pd.concat([fit_results, bin_fit_results], ignore_index=True)

    bin_error = (fit_results['bin_center'][1] - fit_results['bin_center'][0])/2.
    fit_results['bin_error'] = bin_error

    signal_func_name = kwargs.get('signal_func_name', 'gaus')
    fit_results['mean_err'] = fit_results[f'{signal_func_name}_sigma'] / np.sqrt(fit_results['unnorm_integral'])
    fit_results['res'] = fit_results[f'{signal_func_name}_sigma'] / fit_results[f'{signal_func_name}_mean']
    fit_results['res_err'] = np.sqrt((fit_results[f'{signal_func_name}_sigma_err']/fit_results[f'{signal_func_name}_mean'])**2 + (fit_results[f'{signal_func_name}_sigma']*fit_results['mean_err']/fit_results[f'{signal_func_name}_mean']**2)**2)

    graph_mean = TGraphErrors(len(fit_results), np.array(fit_results['bin_center']), np.array(fit_results[f'{signal_func_name}_mean']), np.array(fit_results['bin_error']), np.array(fit_results['mean_err']))
    graph_res = TGraphErrors(len(fit_results), np.array(fit_results['bin_center']), np.array(fit_results['res']), np.array(fit_results['bin_error']),  np.array(fit_results['res_err']))
    xmin = h2.GetXaxis().GetBinLowEdge(kwargs.get('first_bin_fit_by_slices'))
    xmax = h2.GetXaxis().GetBinUpEdge(kwargs.get('last_bin_fit_by_slices'))

    bethe_bloch_func = TF1('bethe_bloch_func', BetheBloch, xmin, xmax, 5)
    bethe_bloch_pars = kwargs.get('bethe_bloch_pars', deepcopy(DEFAULT_BETHEBLOCH_PARS))
    bethe_bloch_func.SetParNames(*bethe_bloch_pars.keys())
    bethe_bloch_func.SetParameters(*bethe_bloch_pars.values())
    graph_mean.Fit(bethe_bloch_func, 'RMS+')

    const_fit = TF1('const_fit', '[0]', xmin, xmax)
    const_fit.SetParameter(0, 0.09)
    graph_res.Fit(const_fit, 'RMS+')

    print(tc.GREEN+'[INFO]:'+tc.RESET+'-------- BETHE BLOCH PARAMETRISATION --------')
    for ipar, par in bethe_bloch_pars.items():
        bethe_bloch_pars[ipar] = bethe_bloch_func.GetParameter(ipar)
        print(tc.GREEN+'[INFO]:'+tc.RED+f'{ipar}:'+tc.RESET, bethe_bloch_func.GetParameter(ipar))
    print(tc.GREEN+'[INFO]:'+tc.RED+f'\tchi2 / NDF:'+tc.RESET, bethe_bloch_func.GetChisquare(), '/', bethe_bloch_func.GetNDF())

    canvas = TCanvas('canvas', '')
    canvas.cd()
    h2.Draw('colz')
    bethe_bloch_func.Draw('same')
    
    output_file.cd()
    graph_mean.Write('g_FitBySlices')
    graph_res.Write('g_ResBySlices')
    h2.Write('h2_dEdx')
    bethe_bloch_func.Write('f_BetheBlochCurve')
    canvas.Write('c_dEdxAndBetheBloch')

    return bethe_bloch_pars

def cluster_size_calibration(clsize: RooRealVar, h2: TH2F, output_file: TDirectory, fit_func, fit_params, fit_mc: bool = False, **kwargs) -> dict:
    '''
        Perform a calibration fit on a 2D histogram.
        The histogram is sliced along the x-axis and fitted with a double Gaussian.
        The variables of the fit are stored in a TGraphErrors.
        The bin error is calculated as the bin width.
        The mean error is calculated as sigma / sqrt(n_entries).

        IMPORTANT REQUIREMENT: if not provided, the names of the parameters of the gaussian must be 'mean1', 'sigma1'

        Parameters:
        - h2: TH2F
            The 2D histogram to be calibrated.
        - output_file: TDirectory
            The output file where the TGraphErrors will be stored.
        - fitter: Roofitter
            The fitter object to be used.
        - charge: float
            The charge of the particle.
        - fit_charge: bool
            If True, the charge exponent will be fit, while the other parameters will be fixed.
            If False, the charge exponent will be fixed, while the other parameters will be fit.
        - fit_mc: bool
            If True, only the constant paramter will be left free.
        - **kwargs:
            Additional arguments to be passed to the fit_by_slices function.
            -> first_bin_fit_by_slices: int
                First bin to be fitted by slices.
            -> last_bin_fit_by_slices: int
                Last bin to be fitted by slices.
            -> mean_label: str
                The name of the mean parameter of the gaussian.
            -> sigma_label: str
                The name of the sigma parameter of the gaussian.
            -> sigma_err_label: str
                The name of the error of the sigma parameter of the gaussian.
            -> simil_bethe_bloch_pars: dict
                The parameters of the simil Bethe-Bloch function.
                The parameters are kp1, kp2, kp3, charge, kp4.
    '''

    fit_results = pd.DataFrame()
    first_fit_bin = kwargs.get('first_bin_fit_by_slices')
    for xbin in range(kwargs.get('first_bin_fit_by_slices'), kwargs.get('last_bin_fit_by_slices')+1):
        h = h2.ProjectionY(f'h_{xbin}', xbin, xbin, 'e')
        bin_center = h2.GetXaxis().GetBinCenter(xbin)
            
        datahist =  RooDataHist('data', 'data', RooArgList(clsize), h)

        fit_ranges = kwargs.get('fit_ranges', None)
        if fit_ranges is not None:
            fit_range = fit_ranges[xbin-first_fit_bin]
            fit_func.fitTo(datahist, RooFit.Save(), RooFit.Range(fit_range[0], fit_range[1]))
        else:
            fit_func.fitTo(datahist, RooFit.Save())

        frame = clsize.frame()
        datahist.plotOn(frame)
        fit_func.plotOn(frame)

        ## put fit parameters on the plot
        text = TPaveText(0.7, 0.7, 0.9, 0.9, 'NDC')
        for name, param in fit_params.items():
            text.AddText(f'{name} = {param.getVal():.4f} +/- {param.getError():.4f}')
        frame.addObject(text)
        frame.Draw()
        output_file.cd()
        frame.Write(f'frame_{xbin}')

        bin_fit_results = pd.DataFrame.from_dict({key: [value.getVal()] for key, value in fit_params.items()})
        bin_fit_results_error = pd.DataFrame.from_dict({key+'_err': [value.getError()] for key, value in fit_params.items()})
        bin_fit_results['bin_center'] = bin_center
        bin_fit_results['integral'] = fit_func.createIntegral(RooArgList(clsize)).getVal()
        bin_fit_results['unnorm_integral'] = fit_func.createIntegral(RooArgList(clsize)).getVal() * h.Integral()
        bin_fit_results = pd.concat([bin_fit_results, bin_fit_results_error], axis=1)
        fit_results = pd.concat([fit_results, bin_fit_results], ignore_index=True)

    bin_error = (fit_results['bin_center'][1] - fit_results['bin_center'][0])/2.
    fit_results['bin_error'] = bin_error

    signal_func_name = kwargs.get('signal_func_name', 'gaus')
    fit_results['mean_err'] = fit_results[f'{signal_func_name}_sigma'] / np.sqrt(fit_results['unnorm_integral'])
    fit_results['res'] = fit_results[f'{signal_func_name}_sigma'] / fit_results[f'{signal_func_name}_mean']
    fit_results['res_err'] = np.sqrt((fit_results[f'{signal_func_name}_sigma_err']/fit_results[f'{signal_func_name}_mean'])**2 + (fit_results[f'{signal_func_name}_sigma']*fit_results['mean_err']/fit_results[f'{signal_func_name}_mean']**2)**2)

    graph_mean = TGraphErrors(len(fit_results), np.array(fit_results['bin_center']), np.array(fit_results[f'{signal_func_name}_mean']), np.array(fit_results['bin_error']), np.array(fit_results['mean_err']))
    graph_res = TGraphErrors(len(fit_results), np.array(fit_results['bin_center']), np.array(fit_results['res']), np.array(fit_results['bin_error']),  np.array(fit_results['res_err']))
    graph_int = TGraphErrors(len(fit_results), np.array(fit_results['bin_center']), np.array(fit_results['integral']), np.array(fit_results['bin_error']), np.zeros(len(fit_results['bin_center'])))
    #graph_tau = TGraphErrors(len(fit_results), np.array(fit_results['bin_center']), np.array(fit_results['tau1']), np.array(fit_results['bin_error']), np.array(fit_results['tau1_err']))
    xmin = h2.GetXaxis().GetBinLowEdge(kwargs.get('first_bin_fit_by_slices'))
    xmax = h2.GetXaxis().GetBinUpEdge(kwargs.get('last_bin_fit_by_slices'))
    
    simil_bethe_bloch_func = TF1('simil_bethe_bloch_func', '[0]/x^[1] + [2]', xmin, xmax)
    DEFAULT_PARAMS = {'kp1': 2.6, 'kp2': 2., 'kp3': 5.5}
    simil_bethe_bloch_pars = kwargs.get('simil_bethe_bloch_pars', deepcopy(DEFAULT_PARAMS))
    simil_bethe_bloch_func.SetParameters(*simil_bethe_bloch_pars.values())

    if fit_mc:
        simil_bethe_bloch_func.FixParameter(0, simil_bethe_bloch_pars['kp1'])
        simil_bethe_bloch_func.FixParameter(1, simil_bethe_bloch_pars['kp2'])
        simil_bethe_bloch_func.SetParameter(2, simil_bethe_bloch_pars['kp3'])
    else:
        simil_bethe_bloch_func.SetParLimits(0, 0., 10.)
        simil_bethe_bloch_func.SetParLimits(1, 0., 10.)

    simil_bethe_bloch_func.SetParNames(*simil_bethe_bloch_pars.keys())
    graph_mean.Fit(simil_bethe_bloch_func, 'RMS+')
    
    resolution_fit = TF1('resolution_fit', '[0]*ROOT::Math::erf((x - [1])/[2])', xmin, xmax)
    resolution_fit.SetParameter(0, 0.24)
    resolution_fit.SetParameter(1, -0.32)
    resolution_fit.SetParameter(2, 1.53)
    graph_res.Fit(resolution_fit, 'RMS+')
    resolution_params = {'rp0': resolution_fit.GetParameter(0),
                         'rp1': resolution_fit.GetParameter(1),
                         'rp2': resolution_fit.GetParameter(2)
                        }

    print(tc.GREEN+'[INFO]:'+tc.RESET+'-------- BETHE BLOCH PARAMETRISATION --------')
    for ipar, par in enumerate(simil_bethe_bloch_pars.keys()):
        simil_bethe_bloch_pars[par] = simil_bethe_bloch_func.GetParameter(ipar)
        print(tc.GREEN+'[INFO]:'+tc.RED+f'{par}:'+tc.RESET, simil_bethe_bloch_pars[par])
    print(tc.GREEN+'[INFO]:'+tc.RED+f'\tchi2 / NDF:'+tc.RESET, simil_bethe_bloch_func.GetChisquare(), '/', simil_bethe_bloch_func.GetNDF())
    
    output_file.cd()
    graph_mean.SetTitle('; ; #mu_{1} [a.u.]')
    graph_mean.Write('g_MeanBySlices')
    graph_res.SetTitle('; ; #sigma_{1}/#mu_{1}')
    graph_res.Write('g_ResBySlices')
    graph_int.Write('g_IntegralBySlices')

    return simil_bethe_bloch_pars, resolution_params