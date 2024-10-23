import os
import numpy as np
from copy import deepcopy
from ROOT import TH2F, TGraphErrors, TDirectory, TF1, gInterpreter, TCanvas

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BETHEBLOCH_DIR = os.path.join(CURRENT_DIR, 'BetheBloch.hh')
gInterpreter.ProcessLine(f'#include "{BETHEBLOCH_DIR}"')
from ROOT import BetheBloch

from torchic.core.fitter import Fitter, fit_by_slices
from torchic.utils.terminal_colors import TerminalColors as tc

DEFAULT_BETHEBLOCH_PARS = { # params for TPC He3 pp
                            'kp1': -241.490, 
                            'kp2': 0.374245,
                            'kp3': 1.397847,
                            'kp4': 1.078250,
                            'kp5': 2.048336
                          }
def pyBetheBloch(betagamma, kp1, kp2, kp3, kp4, kp5):
    beta = betagamma / np.sqrt(1 + betagamma**2)
    aa = beta**kp4
    bb = (1/betagamma)**kp5
    bb = np.log(bb + kp3)
    return (kp2 - aa - bb) * kp1 / aa

def pySimilBetheBloch(betagamma, kp1, kp2):
    return kp1 / betagamma**2 + kp2

def bethe_bloch_calibration(h2: TH2F, output_file: TDirectory, **kwargs) -> dict:
    '''
        Perform a Bethe-Bloch calibration on a 2D histogram.
        The histogram is sliced along the x-axis and fitted with a Gaussian.
        The mean and sigma of the Gaussian are stored in a TGraphErrors.
        The bin error is calculated as the bin width.
        The mean error is calculated as mean * expected resolution.

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
    '''

    fitter = Fitter('gaus')
    kwargs['fit_options'] = 'QL+' # fit option for the slice fit
    fit_results = fit_by_slices(h2, fitter, init_mode='gaus', **kwargs)
    bin_error = (fit_results['bin_center'][1] - fit_results['bin_center'][0])/2.
    fit_results['bin_error'] = bin_error
    fit_results['mean_err'] = fit_results['mean'] * 0.09

    graph_mean = TGraphErrors(len(fit_results), np.array(fit_results['bin_center']), np.array(fit_results['mean']), np.array(fit_results['bin_error']), np.array(fit_results['mean_err']))
    xmin = h2.GetXaxis().GetBinLowEdge(kwargs.get('first_bin_fit_by_slices'))
    xmax = h2.GetXaxis().GetBinUpEdge(kwargs.get('last_bin_fit_by_slices'))
    bethe_bloch_func = TF1('bethe_bloch_func', BetheBloch, xmin, xmax, 5)
    bethe_bloch_pars = kwargs.get('bethe_bloch_pars', deepcopy(DEFAULT_BETHEBLOCH_PARS))
    bethe_bloch_func.SetParNames(*bethe_bloch_pars.keys())
    bethe_bloch_func.SetParameters(*bethe_bloch_pars.values())
    graph_mean.Fit(bethe_bloch_func, 'RMS+')

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
    h2.Write('h2_dEdx')
    bethe_bloch_func.Write('f_BetheBlochCurve')
    canvas.Write('c_dEdxAndBetheBloch')

    return bethe_bloch_pars

def cluster_size_calibration(h2: TH2F, output_file: TDirectory, **kwargs) -> dict:
    '''
        Perform a calibration fit on a 2D histogram.
        The histogram is sliced along the x-axis and fitted with a double Gaussian.
        The variables of the fit are stored in a TGraphErrors.
        The bin error is calculated as the bin width.
        The mean error is calculated as sigma / sqrt(n_entries).

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
    '''

    #fitter = Fitter('[norm0]*exp(-0.5*(x-[mean0])*(x-[mean0])/([sigma0]*[sigma0])) + [norm1]*exp(-0.5*(x-[mean1])*(x-[mean1])/([sigma1]*[sigma1]))')
    fitter = Fitter('x < ([mean0] + [tau0]*[sigma0]) ? [norm0]*exp( -0.5*(x-[mean0])*(x-[mean0])/([sigma0]*[sigma0]) ) : [norm0]*exp(-(x-[mean0]-0.5*[sigma0]*[tau0])*[tau0]/[sigma0]) \
                    + [norm1]*exp(-0.5*(x-[mean1])*(x-[mean1])/([sigma1]*[sigma1]))')
    #                # + x < ([mean1] + [tau1]*[sigma1]) ? [norm1]*exp( -0.5*(x-[mean1])*(x-[mean1])/([sigma1]*[sigma1]) ) : [norm1]*exp(-(x-[mean1]-0.5*[sigma1]*[tau1])*[tau1]/[sigma1])')
    fitter.set_param_range('norm0', 2e3, 1e5)
    fitter.set_param_range('norm1', 100., 1e5)
    fitter.set_param_range('mean0', 2.5, 1.5, 3.5)
    fitter.set_param_range('mean1', 7., 5.5, 10.)
    fitter.set_param_range('sigma0', 0.5, 0.2, 1.5)
    fitter.set_param_range('sigma1', 0.5, 0.2, 2.)
    fitter.set_param_range('tau0', 0., 4.)
    #fitter.set_param_range('tau1', 0., 2.)
    kwargs['fit_options'] = 'QLR+' # fit option for the slice fit
    #kwargs['init_mode'] = 'multi_gaus'
    kwargs['n_components'] = 2
    if kwargs.get('save_slices', False):
        th1_dir = output_file.mkdir('Slices')
        kwargs['output_dir'] = th1_dir
    fit_results = fit_by_slices(h2, fitter, **kwargs)

    bin_error = (fit_results['bin_center'][1] - fit_results['bin_center'][0])/2.
    fit_results['bin_error'] = bin_error
    gaussian_integral = lambda norm, sigma: norm * np.sqrt(2*np.pi) * sigma
    fit_results['integral1'] = gaussian_integral(fit_results['norm1'], fit_results['sigma1'])
    fit_results['mean_err1'] = fit_results['sigma1'] / np.sqrt(fit_results['integral1'])
    fit_results['res1'] = fit_results['sigma1'] / fit_results['mean1']
    fit_results['res1_err'] = np.sqrt((fit_results['sigma1_err']/fit_results['mean1'])**2 + (fit_results['sigma1']*fit_results['mean_err1']/fit_results['mean1']**2)**2)

    graph_mean = TGraphErrors(len(fit_results), np.array(fit_results['bin_center']), np.array(fit_results['mean1']), np.array(fit_results['bin_error']), np.array(fit_results['mean_err1']))
    graph_sigma = TGraphErrors(len(fit_results), np.array(fit_results['bin_center']), np.array(fit_results['sigma1']), np.array(fit_results['bin_error']), np.array(fit_results['sigma1_err']))
    graph_res = TGraphErrors(len(fit_results), np.array(fit_results['bin_center']), np.array(fit_results['res1']), np.array(fit_results['bin_error']),  np.array(fit_results['res1_err']))
    #graph_tau = TGraphErrors(len(fit_results), np.array(fit_results['bin_center']), np.array(fit_results['tau1']), np.array(fit_results['bin_error']), np.array(fit_results['tau1_err']))
    xmin = h2.GetXaxis().GetBinLowEdge(kwargs.get('first_bin_fit_by_slices'))
    xmax = h2.GetXaxis().GetBinUpEdge(kwargs.get('last_bin_fit_by_slices'))
    
    simil_bethe_bloch_func = TF1('simil_bethe_bloch_func', '[kp1]/x^2 + [kp2]', xmin, xmax)
    simil_bethe_bloch_pars = kwargs.get('simil_bethe_bloch_pars', deepcopy(DEFAULT_BETHEBLOCH_PARS))
    simil_bethe_bloch_pars = {k: v for k, v in list(simil_bethe_bloch_pars.items())[:2]}
    simil_bethe_bloch_func.SetParNames(*simil_bethe_bloch_pars.keys())
    #simil_bethe_bloch_func.SetParameters(*simil_bethe_bloch_pars.values())
    graph_mean.Fit(simil_bethe_bloch_func, 'RMS+')
    
    resolution_fit = TF1('resolution_fit', '[0]', xmin, xmax)
    resolution_fit.SetParameter(0, 0.1)
    graph_res.Fit(resolution_fit, 'RMS+')
    resolution = resolution_fit.GetParameter(0)

    print(tc.GREEN+'[INFO]:'+tc.RESET+'-------- BETHE BLOCH PARAMETRISATION --------')
    for ipar, par in simil_bethe_bloch_pars.items():
        simil_bethe_bloch_pars[ipar] = simil_bethe_bloch_func.GetParameter(ipar)
        print(tc.GREEN+'[INFO]:'+tc.RED+f'{ipar}:'+tc.RESET, simil_bethe_bloch_func.GetParameter(ipar))
    print(tc.GREEN+'[INFO]:'+tc.RED+f'\tchi2 / NDF:'+tc.RESET, simil_bethe_bloch_func.GetChisquare(), '/', simil_bethe_bloch_func.GetNDF())
    
    output_file.cd()
    graph_mean.SetTitle('; ; #mu_{1} [a.u.]')
    graph_mean.Write('g_MeanBySlices')
    graph_res.SetTitle('; ; #sigma_{1}/#mu_{1}')
    graph_res.Write('g_ResBySlices')

    return simil_bethe_bloch_pars, resolution