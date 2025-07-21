import numpy as np
import pandas as pd
from ROOT import TFile, TF1, TCanvas, TGraphErrors
from ROOT import RooRealVar

from torchic.core.histogram import AxisSpec, HistLoadInfo, load_hist
from deprecated.torchic.core.roofitter import Roofitter

def purity_lambda(output_file: TFile):
    '''
        The following function performs a ITS calibration (of protons) with a simil-Bethe-Bloch formula.
        The calibration is performed on the <ITS Cluster size> x <cos#lambda> vs. #beta*#gamma plot and the data is retrieved from a tree (O2Physics table).
        The data is preselected with several cuts.
    '''

    input_file = 'data/AnalysisResultsLambdaOmega.root'    
    input_hist = 'lf-tree-creator-cluster-studies/LFTreeCreator/massLambda'
    hist_load_info = HistLoadInfo(input_file, input_hist)
    h2Lambda = load_hist(hist_load_info)

    x = RooRealVar('x', 'x', 1.10, 1.16)
    fitter = Roofitter(x, ['crystal_ball', 'comp_exp'])
    fitter.init_param('crystal_ball_0_mean', 1.115, 1.100, 1.130)
    fitter.init_param('crystal_ball_0_sigma', 0.002, 0.0001, 0.1)
    fitter.init_param('comp_exp_1_alpha', -1, -10, 10)

    dir_lambda = output_file.mkdir('Lambda')

    fit_results = pd.DataFrame()
    
    for ipbin in range(h2Lambda.GetNbinsX()):
        pt = h2Lambda.GetXaxis().GetBinCenter(ipbin)
        if np.abs(pt) < 0.7:
            continue
        hLambda = h2Lambda.ProjectionY(f'hLambda_{pt}', ipbin, ipbin)

        fractions = fitter.fit(hLambda, 1.10, 1.16)
        fit_fractions = fitter.fit_fractions
        dir_lambda.cd()
        fitter.plot(dir_lambda, canvas_name=f'cLambda_{pt}', funcs_to_plot=['crystal_ball_0', 'comp_exp_1'])

        bin_fit_results = pd.DataFrame.from_dict({key: [value] for key, value in fitter.fit_results.items()})
        bin_fit_results['bin_center'] = pt
        bin_fit_results['unnorm_integral'] = bin_fit_results['integral'] * hLambda.Integral()
        bin_fit_results['purity'] = fit_fractions['crystal_ball_0']
        fit_results = pd.concat([fit_results, bin_fit_results], ignore_index=True)
    
    bin_error = (fit_results['bin_center'][1] - fit_results['bin_center'][0]) / 2
    fit_results['bin_error'] = bin_error
    fit_results['mean_error'] = fit_results['crystal_ball_0_sigma'] / fit_results['unnorm_integral']
    fit_results['resolution'] = fit_results['crystal_ball_0_sigma'] / fit_results['crystal_ball_0_mean']
    fit_results['resolution_error'] = fit_results['resolution'] * np.sqrt((fit_results['mean_error'] / fit_results['crystal_ball_0_mean'])**2 + (fit_results['crystal_ball_0_sigma'] / fit_results['crystal_ball_0_sigma'])**2)
        
    graph_mean = TGraphErrors(len(fit_results), np.array(fit_results['bin_center']), np.array(fit_results[f'crystal_ball_0_mean']), np.array(fit_results['bin_error']), np.array(fit_results['mean_error']))
    graph_res = TGraphErrors(len(fit_results), np.array(fit_results['bin_center']), np.array(fit_results['resolution']), np.array(fit_results['bin_error']),  np.array(fit_results['resolution_error']))
    graph_purity = TGraphErrors(len(fit_results), np.array(fit_results['bin_center']), np.array(fit_results['purity']), np.array(fit_results['bin_error']), np.zeros(len(fit_results['bin_center'])))

    dir_lambda.cd()
    graph_mean.Write('gMeanFromFit')
    graph_res.Write('gResolution')
    graph_purity.Write('gPurity')
    h2Lambda.Write()

def purity_omega(output_file: TFile):
    '''
        The following function performs a ITS calibration (of protons) with a simil-Bethe-Bloch formula.
        The calibration is performed on the <ITS Cluster size> x <cos#lambda> vs. #beta*#gamma plot and the data is retrieved from a tree (O2Physics table).
        The data is preselected with several cuts.
    '''

    input_file = 'data/AnalysisResultsLambdaOmega.root'    
    input_hist = 'lf-tree-creator-cluster-studies/LFTreeCreator/massOmega'
    hist_load_info = HistLoadInfo(input_file, input_hist)
    h2Omega = load_hist(hist_load_info)
    h2Omega.RebinX(2)

    x = RooRealVar('x', 'x', 1.62, 1.72)
    fitter = Roofitter(x, ['crystal_ball', 'pol0'])
    fitter.init_param('crystal_ball_0_mean', 1.672, 1.660, 1.684)
    fitter.init_param('crystal_ball_0_sigma', 0.005, 0.0001, 0.1)

    dir_omega = output_file.mkdir('Omega')

    fit_results = pd.DataFrame()
    
    for ipbin in range(h2Omega.GetNbinsX()):
        pt = h2Omega.GetXaxis().GetBinCenter(ipbin)
        if np.abs(pt) < 2:
            continue
        hOmega = h2Omega.ProjectionY(f'hOmega_{pt}', ipbin, ipbin)

        fractions = fitter.fit(hOmega, hOmega.GetXaxis().GetXmin(), hOmega.GetXaxis().GetXmax())
        fit_fractions = fitter.fit_fractions
        dir_omega.cd()
        fitter.plot(dir_omega, canvas_name=f'cOmega_{pt}', funcs_to_plot=['crystal_ball_0', 'pol0_1'])

        bin_fit_results = pd.DataFrame.from_dict({key: [value] for key, value in fitter.fit_results.items()})
        bin_fit_results['bin_center'] = pt
        bin_fit_results['unnorm_integral'] = bin_fit_results['integral'] * hOmega.Integral()
        bin_fit_results['purity'] = fit_fractions['crystal_ball_0']
        fit_results = pd.concat([fit_results, bin_fit_results], ignore_index=True)
    
    bin_error = (fit_results['bin_center'][1] - fit_results['bin_center'][0]) / 2
    fit_results['bin_error'] = bin_error
    fit_results['mean_error'] = fit_results['crystal_ball_0_sigma'] / fit_results['unnorm_integral']
    fit_results['resolution'] = fit_results['crystal_ball_0_sigma'] / fit_results['crystal_ball_0_mean']
    fit_results['resolution_error'] = fit_results['resolution'] * np.sqrt((fit_results['mean_error'] / fit_results['crystal_ball_0_mean'])**2 + (fit_results['crystal_ball_0_sigma'] / fit_results['crystal_ball_0_sigma'])**2)
        
    graph_mean = TGraphErrors(len(fit_results), np.array(fit_results['bin_center']), np.array(fit_results[f'crystal_ball_0_mean']), np.array(fit_results['bin_error']), np.array(fit_results['mean_error']))
    graph_res = TGraphErrors(len(fit_results), np.array(fit_results['bin_center']), np.array(fit_results['resolution']), np.array(fit_results['bin_error']),  np.array(fit_results['resolution_error']))
    graph_purity = TGraphErrors(len(fit_results), np.array(fit_results['bin_center']), np.array(fit_results['purity']), np.array(fit_results['bin_error']), np.zeros(len(fit_results['bin_center'])))

    dir_omega.cd()
    graph_mean.Write('gMeanFromFit')
    graph_res.Write('gResolution')
    graph_purity.Write('gPurity')
    h2Omega.Write()


if __name__ == '__main__':

    output_file_path = 'data/purity_lambda_omega.root'
    output_file = TFile(output_file_path, 'recreate')

    purity_lambda(output_file)
    purity_omega(output_file)

    output_file.Close()
