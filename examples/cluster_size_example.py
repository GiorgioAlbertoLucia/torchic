import numpy as np
from ROOT import TFile, TF1, TCanvas
from ROOT import RooRealVar

from torchic.physics.calibration import cluster_size_parametrisation, cluster_size_calibration
from torchic.physics.ITS import average_cluster_size
from torchic.core.dataset import Dataset
from torchic.core.histogram import AxisSpec
from torchic.core.roofitter import Roofitter

def cluster_size_calibration_example(output_file: TFile):
    '''
        The following function performs a ITS calibration (of protons) with a simil-Bethe-Bloch formula.
        The calibration is performed on the <ITS Cluster size> x <cos#lambda> vs. #beta*#gamma plot and the data is retrieved from a tree (O2Physics table).
        The data is preselected with several cuts.
    '''

    input_file = 'data/LHC23_PbPb_pass4_long_same.root'    
                                            # required tree from O2Physics table that contains the following branches: 
                                            # fEtaPr, fChi2TPCPr, fInnerParamTPCPr, fItsClusterSizePr
    folder_name = 'DF*'                     # name of the folder(s) containing the tree (O2Physics table)
    tree_name = 'O2lithium4table'             # name of the tree (O2Physics table) in the input file

    dataset = Dataset(input_file, folder_name=folder_name, tree_name=tree_name, columns=['fPtPr', 'fEtaPr', 'fChi2TPCPr', 'fInnerParamTPCPr', 'fItsClusterSizePr', 'fMassTOFPr', 'fNSigmaTPCPr'])

    dataset.query('fPtPr < 0.0', inplace=True)
    dataset.query('abs(fEtaPr) < 0.8', inplace=True)
    dataset.query('fChi2TPCPr > 0.5', inplace=True)

    resolution_params = {
            'kp0': 1.22204e-02,
            'kp1': 7.48467e-01,
        }
    expTOFmassPr = 0.9487
    dataset['fExpTOFMassPr'] = expTOFmassPr
    dataset['fSigmaTOFMassPr'] = (resolution_params['kp0'] * np.exp(resolution_params['kp1'] * np.abs(dataset['fPtPr']))) * expTOFmassPr
    dataset['fNSigmaTOFPr'] = (dataset['fMassTOFPr'] - dataset['fExpTOFMassPr']) / dataset['fSigmaTOFMassPr']
    dataset.query('(fPtPr < 0.8) or (-1 < fNSigmaTOFPr < 1)', inplace=True)
    dataset.query('(fPtPr > 0.8) or (-5 < fNSigmaTOFPr < 5)', inplace=True)

    dataset['fITSAvgClSize'], __ = average_cluster_size(dataset['fItsClusterSizePr'])
    dataset['fITSClSizeCosLam'] = dataset['fITSAvgClSize'] / np.cosh(dataset['fEtaPr'])

    dataset['fBetaGamma'] = dataset['fInnerParamTPCPr'] / 0.93827

    cluster_size_dir = output_file.mkdir('ClusterSizeCalibration')

    axis_spec_p = AxisSpec(100, 0., 5, 'h2_nSigma', ';#it{p}_{TPC} (GeV/#it{c}); n#sigma_{TPC}')
    axis_spec_bg = AxisSpec(40, 0., 4, 'h2_dEdxBetaGamma', ';#beta#gamma; Expected dE/dx_{TPC}')
    axis_spec_nsigma_its = AxisSpec(100, -5., 5, 'h2_nSigmaITS', ';#it{p}_{TPC} (GeV/#it{c}); n#sigma_{ITS}')
    axis_spec_cl = AxisSpec(62, -0.5, 15.5, 'h2_ItsClusterSize', ';#beta#gamma; #LT ITS cluster size #GT #times #LT cos#lambda #GT')
    axis_spec_cl_exp = AxisSpec(62, -0.5, 15.5, 'h2_ExpItsClusterSize', ';#beta#gamma; #LT ITS cluster size #GT #times #LT cos#lambda #GT')

    # cluster size calibration
    h2_clsize = dataset.build_hist('fBetaGamma', 'fITSClSizeCosLam', axis_spec_bg, axis_spec_cl)
    
    cl_options = {'first_bin_fit_by_slices': h2_clsize.GetXaxis().FindBin(0.5),
                  'last_bin_fit_by_slices': h2_clsize.GetXaxis().FindBin(2.7),
                  'output_dir': cluster_size_dir,
                  'fit_range': [1., 9.5],
                  'simil_bethe_bloch_pars': {'kp1': 2.781,
                                             'kp2': 1.159,
                                             'kp3': 5.116,
                                            },
                  'mean_label': 'gaus_0_mean',
                  'sigma_label': 'gaus_0_sigma',
                  'sigma_err_label': 'gaus_0_sigma_err',
                  'integral_label': 'gaus_0_integral',
    }

    x = RooRealVar('x', 'x', 1., 9.5)
    fitter = Roofitter(x, ['gaus'])
    fitter.init_param('gaus_0_mean', 5., 2.5, 6.)
    fitter.init_param('gaus_0_sigma', 0.5, 0.2, 2.)
    cluster_size_pars, its_resolution = cluster_size_calibration(h2_clsize, cluster_size_dir, fitter, **cl_options)

    simil_bethe_bloch_func = TF1('f_SimilBetheBloch', '[kp1]/x^[kp2] + [kp3]', 0.3, 4.) # function used in cluster_size_calibration
    for ipar, par in cluster_size_pars.items():
        simil_bethe_bloch_func.SetParameter(ipar, par)
    c_clsize = TCanvas('c_ClSizeAndCurve', 'canvas', 800, 600)
    h2_clsize.Draw('colz')
    simil_bethe_bloch_func.Draw('same')

    # draw results
    cluster_size_dir.cd()
    h2_clsize.Write()
    simil_bethe_bloch_func.Write()
    c_clsize.Write()

    dataset['fExpClSizeCosLam'] = cluster_size_parametrisation(dataset['fBetaGamma'], *cluster_size_pars.values())
    dataset['fNSigmaITS'] = (dataset['fITSClSizeCosLam'] - dataset['fExpClSizeCosLam']) / (its_resolution*dataset['fExpClSizeCosLam'])
    h2_clsize_exp = dataset.build_hist('fBetaGamma', 'fExpClSizeCosLam', axis_spec_bg, axis_spec_cl_exp)
    h2_nsigma_its = dataset.build_hist('fInnerParamTPCPr', 'fNSigmaITS', axis_spec_p, axis_spec_nsigma_its)
    
    # draw results
    cluster_size_dir.cd()
    h2_clsize_exp.Write()
    h2_nsigma_its.Write()

def cluster_size_cut_example(output_file: TFile):
    '''
        Demonstration of the cluster size cut implementation and its effect on the data.
    '''

    input_file = 'data/LHC23_PbPb_pass4_long_same.root'    
                                            # required tree from O2Physics table that contains the following branches: 
                                            # fEtaPr, fChi2TPCPr, fInnerParamTPCPr, fItsClusterSizePr
    folder_name = 'DF*'                     # name of the folder(s) containing the tree (O2Physics table)
    tree_name = 'O2lithium4table'             # name of the tree (O2Physics table) in the input file

    dataset = Dataset(input_file, folder_name=folder_name, tree_name=tree_name, columns=['fPtPr', 'fEtaPr', 'fChi2TPCPr', 'fInnerParamTPCPr', 'fItsClusterSizePr', 'fMassTOFPr', 'fNSigmaTPCPr'])

    dataset.query('fPtPr < 0.0', inplace=True)
    dataset.query('abs(fEtaPr) < 0.8', inplace=True)
    dataset.query('fChi2TPCPr > 0.5', inplace=True)
    dataset.query('abs(fNSigmaTPCPr) < 2', inplace=True)

    resolution_params = {
            'kp0': 1.22204e-02,
            'kp1': 7.48467e-01,
        }
    expTOFmassPr = 0.9487
    dataset['fExpTOFMassPr'] = expTOFmassPr
    dataset['fSigmaTOFMassPr'] = (resolution_params['kp0'] * np.exp(resolution_params['kp1'] * np.abs(dataset['fPtPr']))) * expTOFmassPr
    dataset['fNSigmaTOFPr'] = (dataset['fMassTOFPr'] - dataset['fExpTOFMassPr']) / dataset['fSigmaTOFMassPr']
    dataset.query('(fPtPr < 0.8) or (-1 < fNSigmaTOFPr < 1)', inplace=True)
    dataset.query('(fPtPr > 0.8) or (-5 < fNSigmaTOFPr < 5)', inplace=True)

    dataset['fITSAvgClSize'], __ = average_cluster_size(dataset['fItsClusterSizePr'])
    dataset['fITSClSizeCosLam'] = dataset['fITSAvgClSize'] / np.cosh(dataset['fEtaPr'])

    dataset['fBetaGamma'] = dataset['fInnerParamTPCPr'] / 0.93827

    cluster_size_dir = output_file.mkdir('ClusterSizeCut')

    axis_spec_p = AxisSpec(100, 0., 5, 'h2_nSigma', ';#it{p}_{TPC} (GeV/#it{c}); n#sigma_{TPC}')
    axis_spec_bg = AxisSpec(80, 0., 4, 'h2_dEdxBetaGamma', ';#beta#gamma; Expected dE/dx_{TPC}')
    axis_spec_nsigma_its = AxisSpec(100, -5., 5, 'h2_nSigmaITS', ';#it{p}_{TPC} (GeV/#it{c}); n#sigma_{ITS}')
    axis_spec_cl = AxisSpec(120, 0., 15., 'h2_ItsClusterSize', ';#beta#gamma; #LT ITS cluster size #GT #times #LT cos#lambda #GT')
    axis_spec_cl_exp = AxisSpec(62, -0.5, 15.5, 'h2_ExpItsClusterSize', ';#beta#gamma; #LT ITS cluster size #GT #times #LT cos#lambda #GT')

    # cluster size calibration
    h2_clsize_before_cut = dataset.build_hist('fBetaGamma', 'fITSClSizeCosLam', axis_spec_bg, axis_spec_cl)

    cluster_size_pars = {
        'kp1': 0.903,
        'kp2': 2.014,
        'kp3': 2.440,
    }

    simil_bethe_bloch_func = TF1('f_SimilBetheBloch', '[kp1]/x^[kp2] + [kp3]', 0.3, 4.) # function used in cluster_size_calibration
    for ipar, par in cluster_size_pars.items():
        simil_bethe_bloch_func.SetParameter(ipar, par)
    c_clsize = TCanvas('c_ClSizeAndCurve', 'canvas', 800, 600)
    h2_clsize_before_cut.Draw('colz')
    simil_bethe_bloch_func.Draw('same')

    # draw results
    cluster_size_dir.cd()
    h2_clsize_before_cut.Write('h2_clsize_before_cut')
    simil_bethe_bloch_func.Write()
    c_clsize.Write()

    its_resolution = 0.2
    dataset['fExpClSizeCosLam'] = cluster_size_parametrisation(dataset['fBetaGamma'], *cluster_size_pars.values())
    dataset['fNSigmaITS'] = (dataset['fITSClSizeCosLam'] - dataset['fExpClSizeCosLam']) / (its_resolution*dataset['fExpClSizeCosLam'])
    h2_clsize_exp = dataset.build_hist('fBetaGamma', 'fExpClSizeCosLam', axis_spec_bg, axis_spec_cl_exp)
    h2_nsigma_its = dataset.build_hist('fInnerParamTPCPr', 'fNSigmaITS', axis_spec_p, axis_spec_nsigma_its)

    dataset.query('fNSigmaITS > -1.5', inplace=True)
    h2_clsize_after_cut = dataset.build_hist('fBetaGamma', 'fITSClSizeCosLam', axis_spec_bg, axis_spec_cl)
    
    # draw results
    cluster_size_dir.cd()
    h2_clsize_exp.Write()
    h2_nsigma_its.Write()
    h2_clsize_after_cut.Write('h2_clsize_after_cut')    

if __name__ == '__main__':
    output_file_path = 'data/cluster_size_calibration_example.root'
    output_file = TFile(output_file_path, 'recreate')
    
    cluster_size_calibration_example(output_file)
    cluster_size_cut_example(output_file)

    output_file.Close()
