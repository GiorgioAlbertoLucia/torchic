import numpy as np
from math import erf
from ROOT import TFile, TF1, TCanvas
from ROOT import RooRealVar

from deprecated.torchic.physics.calibration import cluster_size_parametrisation, np_cluster_size_resolution, cluster_size_calibration
from torchic.physics.ITS import average_cluster_size
from torchic.core.dataset import Dataset
from torchic.core.histogram import AxisSpec
from deprecated.torchic.core.roofitter import Roofitter

PR_MASS = 0.93827   # proton mass in GeV/c^2

def cluster_size_calibration_example_proton(output_file: TFile):
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

    dataset = Dataset.from_root(input_file, folder_name=folder_name, tree_name=tree_name, columns=['fPtPr', 'fEtaPr', 'fChi2TPCPr', 'fInnerParamTPCPr', 'fItsClusterSizePr', 'fMassTOFPr', 'fNSigmaTPCPr'])

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

    dataset['fBetaGamma'] = dataset['fInnerParamTPCPr'] / PR_MASS

    cluster_size_dir = output_file.mkdir('ClusterSizeCalibration')

    axis_spec_p = AxisSpec(100, 0., 5, 'h2_nSigma', ';#it{p}_{TPC} (GeV/#it{c}); n#sigma_{TPC}')
    axis_spec_bg = AxisSpec(40, 0., 4, 'h2_dEdxBetaGamma', ';#beta#gamma; Expected dE/dx_{TPC}')
    axis_spec_nsigma_its = AxisSpec(100, -5., 5, 'h2_nSigmaITS', ';#it{p}_{TPC} (GeV/#it{c}); n#sigma_{ITS}')
    axis_spec_cl = AxisSpec(62, -0.5, 15.5, 'h2_ItsClusterSize', ';#beta#gamma; #LT ITS cluster size #GT #times #LT cos#lambda #GT')
    axis_spec_cl_exp = AxisSpec(62, -0.5, 15.5, 'h2_ExpItsClusterSize', ';#beta#gamma; #LT ITS cluster size #GT #times #LT cos#lambda #GT')

    # cluster size calibration
    h2_clsize = dataset.build_th2('fBetaGamma', 'fITSClSizeCosLam', axis_spec_bg, axis_spec_cl)
    
    cl_options = {'first_bin_fit_by_slices': h2_clsize.GetXaxis().FindBin(0.5),
                  'last_bin_fit_by_slices': h2_clsize.GetXaxis().FindBin(2.7),
                  'output_dir': cluster_size_dir,
                  'fit_range': [1., 9.5],
                  'simil_bethe_bloch_pars': {'kp1': 1.055,
                                             'kp2': 1.570,
                                             'kp3': 2.317,
                                             'charge': 1.,
                                             'kp4': 1.,
                                            },
                  'signal_func_name': 'exp_mod_gaus_0',
    }

    x = RooRealVar('x', 'x', 1., 9.5)
    fitter = Roofitter(x, ['exp_mod_gaus'])
    fitter.init_param('exp_mod_gaus_0_mean', 5., 2.5, 6.)
    fitter.init_param('exp_mod_gaus_0_sigma', 0.5, 0.2, 2.)
    cluster_size_pars, its_resolution_pars = cluster_size_calibration(h2_clsize, cluster_size_dir, fitter, **cl_options)

    simil_bethe_bloch_func = TF1('f_SimilBetheBloch', '([0]/x^[1] + [2]) * [3]^[4]', 0.3, 4.) # function used in cluster_size_calibration
    for ipar, par in enumerate(cluster_size_pars.values()):
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
    dataset['fResolutionITS'] = np_cluster_size_resolution(dataset['fBetaGamma'], *its_resolution_pars.values())
    dataset['fNSigmaITS'] = (dataset['fITSClSizeCosLam'] - dataset['fExpClSizeCosLam']) / (dataset['fResolutionITS']*dataset['fExpClSizeCosLam'])
    h2_clsize_exp = dataset.build_th2('fBetaGamma', 'fExpClSizeCosLam', axis_spec_bg, axis_spec_cl_exp)
    h2_nsigma_its = dataset.build_th2('fInnerParamTPCPr', 'fNSigmaITS', axis_spec_p, axis_spec_nsigma_its)
    
    # draw results
    cluster_size_dir.cd()
    h2_clsize_exp.Write()
    h2_nsigma_its.Write()

def cluster_size_calibration_example_proton_pid(output_file: TFile):
    '''
        The following function performs a ITS calibration (of protons) with a simil-Bethe-Bloch formula.
        The calibration is performed on the <ITS Cluster size> x <cos#lambda> vs. #beta*#gamma plot and the data is retrieved from a tree (O2Physics table).
        The data is preselected with several cuts.
    '''

    input_file = 'data/LHC22o_pass7_minBias_small.root'    
                                            # required tree from O2Physics table that contains the following branches: 
                                            # fEtaPr, fChi2TPCPr, fInnerParamTPCPr, fItsClusterSizePr
    folder_name = 'DF*'                     # name of the folder(s) containing the tree (O2Physics table)
    tree_name = 'O2clsttable'             # name of the tree (O2Physics table) in the input file

    dataset = Dataset.from_root(input_file, folder_name=folder_name, tree_name=tree_name, columns=['fP', 'fEta', 'fPhi', 'fItsClusterSize', 'fPartID'])

    print(f'\n\n{dataset['fPartID'].unique()=}\n\n')
    dataset.query('fPartID == 4', inplace=True)
    dataset.query('abs(fEta) < 0.8', inplace=True)

    dataset['fITSAvgClSize'], __ = average_cluster_size(dataset['fItsClusterSize'])
    dataset['fITSClSizeCosLam'] = dataset['fITSAvgClSize'] / np.cosh(dataset['fEta'])

    
    dataset['fBetaGamma'] = np.abs(dataset['fP']) / PR_MASS

    cluster_size_dir = output_file.mkdir('ClusterSizeCalibrationPid')

    axis_spec_p = AxisSpec(100, 0., 5, 'h2_nSigma', ';#it{p}_{TPC} (GeV/#it{c}); n#sigma_{TPC}')
    axis_spec_bg = AxisSpec(40, 0., 4, 'h2_dEdxBetaGamma', ';#beta#gamma; Expected dE/dx_{TPC}')
    axis_spec_nsigma_its = AxisSpec(100, -5., 5, 'h2_nSigmaITS', ';#it{p}_{TPC} (GeV/#it{c}); n#sigma_{ITS}')
    axis_spec_cl = AxisSpec(62, -0.5, 15.5, 'h2_ItsClusterSize', ';#beta#gamma; #LT ITS cluster size #GT #times #LT cos#lambda #GT')
    axis_spec_cl_exp = AxisSpec(62, -0.5, 15.5, 'h2_ExpItsClusterSize', ';#beta#gamma; #LT ITS cluster size #GT #times #LT cos#lambda #GT')

    # cluster size calibration
    h2_clsize = dataset.build_th2('fBetaGamma', 'fITSClSizeCosLam', axis_spec_bg, axis_spec_cl)
    
    cl_options = {'first_bin_fit_by_slices': h2_clsize.GetXaxis().FindBin(0.5),
                  'last_bin_fit_by_slices': h2_clsize.GetXaxis().FindBin(2.7),
                  'output_dir': cluster_size_dir,
                  'fit_range': [1., 9.5],
                  'simil_bethe_bloch_pars': {'kp1': 0.903,
                                             'kp2': 2.014,
                                             'kp3': 2.440,
                                             'charge': 1.,
                                             'kp4': 1.,
                                            },
                  'signal_func_name': 'crystal_ball_0',
                  #'signal_func_name': 'exp_mod_gaus_0',
    }

    x = RooRealVar('x', 'x', 0., 15.)
    #fitter = Roofitter(x, ['exp_mod_gaus'])
    #fitter.init_param('exp_mod_gaus_0_mean', 5., 2., 6.)
    #fitter.init_param('exp_mod_gaus_0_sigma', 0.5, 0., 10.)
    #fitter.init_param('exp_mod_gaus_0_tau', 1., -10., 10.)
    fitter = Roofitter(x)
    # signal
    fitter.build_pdf('crystal_ball', double_sided=False)
    fitter.init_param('crystal_ball_0_mean', 5., 0., 15.)
    fitter.init_param('crystal_ball_0_sigma', 0.5, 0., 2.)
    fitter.init_param('crystal_ball_0_alpha', 1., -10., 0.)
    fitter.init_param('crystal_ball_0_n', 1., 0., 30.)
    # bkg
    fitter.build_pdf('exp', exp_offset=False)
    fitter.init_param('exp_1_alpha', -0.1, -1., 1.)
    cluster_size_pars, its_resolution_pars = cluster_size_calibration(h2_clsize, cluster_size_dir, fitter, **cl_options)

    simil_bethe_bloch_func = TF1('f_SimilBetheBloch', '([0]/x^[1] + [2]) * [3]^[4]', 0.3, 4.) # function used in cluster_size_calibration
    for ipar, par in enumerate(cluster_size_pars.values()):
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
    dataset['fResolutionITS'] = np_cluster_size_resolution(dataset['fBetaGamma'], *its_resolution_pars.values())
    dataset['fNSigmaITS'] = (dataset['fITSClSizeCosLam'] - dataset['fExpClSizeCosLam']) / (dataset['fResolutionITS']*dataset['fExpClSizeCosLam'])
    h2_clsize_exp = dataset.build_th2('fBetaGamma', 'fExpClSizeCosLam', axis_spec_bg, axis_spec_cl_exp)
    h2_nsigma_its = dataset.build_th2('fP', 'fNSigmaITS', axis_spec_p, axis_spec_nsigma_its)
    
    # draw results
    cluster_size_dir.cd()
    h2_clsize_exp.Write()
    h2_nsigma_its.Write()

def cluster_size_calibration_example_proton_pid_mc(output_file: TFile):
    '''
        The following function performs a ITS calibration (of protons) with a simil-Bethe-Bloch formula.
        The calibration is performed on the <ITS Cluster size> x <cos#lambda> vs. #beta*#gamma plot and the data is retrieved from a tree (O2Physics table).
        The data is preselected with several cuts.
    '''

    input_file = 'data/MC_LHC24f3_small_pid.root'    
                                            # required tree from O2Physics table that contains the following branches: 
                                            # fEtaPr, fChi2TPCPr, fInnerParamTPCPr, fItsClusterSizePr
    folder_name = 'DF*'                     # name of the folder(s) containing the tree (O2Physics table)
    tree_name = 'O2clsttablemcext'             # name of the tree (O2Physics table) in the input file

    dataset = Dataset.from_root(input_file, folder_name=folder_name, tree_name=tree_name)

    print(f'\n\n{dataset.columns=}')
    print(f'{dataset['fPartID'].unique()=}\n\n')
    dataset.query('fPartID == 4', inplace=True)
    dataset.query('abs(fEta) < 0.8', inplace=True)

    dataset['fITSAvgClSize'], __ = average_cluster_size(dataset['fItsClusterSize'])
    dataset['fITSClSizeCosLam'] = dataset['fITSAvgClSize'] / np.cosh(dataset['fEta'])

    dataset['fBetaGamma'] = np.abs(dataset['fP']) / PR_MASS

    cluster_size_dir = output_file.mkdir('ClusterSizeCalibrationPidMC')

    axis_spec_p = AxisSpec(100, 0., 5, 'h2_nSigma', ';#it{p}_{TPC} (GeV/#it{c}); n#sigma_{TPC}')
    axis_spec_bg = AxisSpec(40, 0., 4, 'h2_dEdxBetaGamma', ';#beta#gamma; Expected dE/dx_{TPC}')
    axis_spec_nsigma_its = AxisSpec(100, -5., 5, 'h2_nSigmaITS', ';#it{p}_{TPC} (GeV/#it{c}); n#sigma_{ITS}')
    axis_spec_cl = AxisSpec(62, -0.5, 15.5, 'h2_ItsClusterSize', ';#beta#gamma; #LT ITS cluster size #GT #times #LT cos#lambda #GT')
    axis_spec_cl_exp = AxisSpec(62, -0.5, 15.5, 'h2_ExpItsClusterSize', ';#beta#gamma; #LT ITS cluster size #GT #times #LT cos#lambda #GT')

    # cluster size calibration
    h2_clsize = dataset.build_th2('fBetaGamma', 'fITSClSizeCosLam', axis_spec_bg, axis_spec_cl)
    
    cl_options = {'first_bin_fit_by_slices': h2_clsize.GetXaxis().FindBin(0.5),
                  'last_bin_fit_by_slices': h2_clsize.GetXaxis().FindBin(2.7),
                  'output_dir': cluster_size_dir,
                  'fit_range': [1., 9.5],
                  'simil_bethe_bloch_pars': {'kp1': 0.903,
                                             'kp2': 2.014,
                                             'kp3': 2.440,
                                             'charge': 1.,
                                             'kp4': 1.,
                                            },
                  'signal_func_name': 'exp_mod_gaus_0',
    }

    x = RooRealVar('x', 'x', 1., 9.5)
    fitter = Roofitter(x, ['exp_mod_gaus'])
    fitter.init_param('exp_mod_gaus_0_mean', 5., 2., 6.)
    fitter.init_param('exp_mod_gaus_0_sigma', 0.5, 0., 10.)
    fitter.init_param('exp_mod_gaus_0_tau', 1., -10., 10.)    
    cluster_size_pars, its_resolution_pars = cluster_size_calibration(h2_clsize, cluster_size_dir, fitter, **cl_options)

    simil_bethe_bloch_func = TF1('f_SimilBetheBloch', '([0]/x^[1] + [2]) * [3]^[4]', 0.3, 4.) # function used in cluster_size_calibration
    for ipar, par in enumerate(cluster_size_pars.values()):
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
    dataset['fResolutionITS'] = np_cluster_size_resolution(dataset['fBetaGamma'], *its_resolution_pars.values())
    dataset['fNSigmaITS'] = (dataset['fITSClSizeCosLam'] - dataset['fExpClSizeCosLam']) / (dataset['fResolutionITS']*dataset['fExpClSizeCosLam'])
    h2_clsize_exp = dataset.build_th2('fBetaGamma', 'fExpClSizeCosLam', axis_spec_bg, axis_spec_cl_exp)
    h2_nsigma_its = dataset.build_th2('fP', 'fNSigmaITS', axis_spec_p, axis_spec_nsigma_its)
    
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

    dataset = Dataset.from_root(input_file, folder_name=folder_name, tree_name=tree_name, columns=['fPtPr', 'fEtaPr', 'fChi2TPCPr', 'fInnerParamTPCPr', 'fItsClusterSizePr', 'fMassTOFPr', 'fNSigmaTPCPr'])

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

    dataset['fBetaGamma'] = dataset['fInnerParamTPCPr'] / PR_MASS

    cluster_size_dir = output_file.mkdir('ClusterSizeCut')

    axis_spec_p = AxisSpec(100, 0., 5, 'h2_nSigma', ';#it{p}_{TPC} (GeV/#it{c}); n#sigma_{TPC}')
    axis_spec_bg = AxisSpec(80, 0., 4, 'h2_dEdxBetaGamma', ';#beta#gamma; Expected dE/dx_{TPC}')
    axis_spec_nsigma_its = AxisSpec(100, -5., 5, 'h2_nSigmaITS', ';#it{p}_{TPC} (GeV/#it{c}); n#sigma_{ITS}')
    axis_spec_cl = AxisSpec(120, 0., 15., 'h2_ItsClusterSize', ';#beta#gamma; #LT ITS cluster size #GT #times #LT cos#lambda #GT')
    axis_spec_cl_exp = AxisSpec(62, -0.5, 15.5, 'h2_ExpItsClusterSize', ';#beta#gamma; #LT ITS cluster size #GT #times #LT cos#lambda #GT')

    # cluster size calibration
    h2_clsize_before_cut = dataset.build_th2('fBetaGamma', 'fITSClSizeCosLam', axis_spec_bg, axis_spec_cl)

    cluster_size_pars = {
        'kp1': 0.903,
        'kp2': 2.014,
        'kp3': 2.440,
        'charge': 1.,
        'kp4': 1.,
    }

    simil_bethe_bloch_func = TF1('f_SimilBetheBloch', '([0]/x^[1] + [2]) * [3]^[4]', 0.3, 4.) # function used in cluster_size_calibration
    for ipar, par in enumerate(cluster_size_pars.values()):
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
    h2_clsize_exp = dataset.build_th2('fBetaGamma', 'fExpClSizeCosLam', axis_spec_bg, axis_spec_cl_exp)
    h2_nsigma_its = dataset.build_th2('fInnerParamTPCPr', 'fNSigmaITS', axis_spec_p, axis_spec_nsigma_its)

    dataset.query('fNSigmaITS > -1.5', inplace=True)
    h2_clsize_after_cut = dataset.build_th2('fBetaGamma', 'fITSClSizeCosLam', axis_spec_bg, axis_spec_cl)
    
    # draw results
    cluster_size_dir.cd()
    h2_clsize_exp.Write()
    h2_nsigma_its.Write()
    h2_clsize_after_cut.Write('h2_clsize_after_cut')  

def cluster_size_calibration_example_he(output_file: TFile):
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

    dataset = Dataset.from_root(input_file, folder_name=folder_name, tree_name=tree_name, columns=['fPtHe3', 'fEtaHe3', 'fChi2TPCHe3', 'fInnerParamTPCHe3', 'fItsClusterSizeHe3', 'fNSigmaTPCHe3', 'fDCAxyHe3', 'fDCAzHe3'])

    # selection cuts from llbariogl/NucleiFlow
    dataset.query('abs(fEtaHe3) < 0.8', inplace=True)
    dataset.query('fPtHe3 < 0', inplace=True)    
    dataset.query('abs(fDCAxyHe3) < 0.1', inplace=True)
    dataset.query('abs(fDCAzHe3) < 1.', inplace=True)
    dataset.query('abs(fNSigmaTPCHe3) < 2.', inplace=True)

    dataset['fITSAvgClSize'], __ = average_cluster_size(dataset['fItsClusterSizeHe3'])
    dataset['fITSClSizeCosLam'] = dataset['fITSAvgClSize'] / np.cosh(dataset['fEtaHe3'])

    dataset['fInnerParamTPCHe3'] = dataset['fInnerParamTPCHe3'] * 2.
    dataset['fBetaGamma'] = np.abs(dataset['fInnerParamTPCHe3']) / 2.8092

    cluster_size_dir = output_file.mkdir('ClusterSizeCalibrationHe')

    axis_spec_p = AxisSpec(100, 0., 5, 'h2_nSigma', ';#it{p}_{TPC} (GeV/#it{c}); n#sigma_{TPC}')
    axis_spec_bg = AxisSpec(40, 0., 4, 'h2_dEdxBetaGamma', ';#beta#gamma; Expected dE/dx_{TPC}')
    axis_spec_nsigma_its = AxisSpec(100, -5., 5, 'h2_nSigmaITS', ';#it{p}_{TPC} (GeV/#it{c}); n#sigma_{ITS}')
    axis_spec_cl = AxisSpec(62, -0.5, 15.5, 'h2_ItsClusterSize', ';#beta#gamma; #LT ITS cluster size #GT #times #LT cos#lambda #GT')
    axis_spec_cl_exp = AxisSpec(62, -0.5, 15.5, 'h2_ExpItsClusterSize', ';#beta#gamma; #LT ITS cluster size #GT #times #LT cos#lambda #GT')

    # cluster size calibration
    h2_clsize = dataset.build_th2('fBetaGamma', 'fITSClSizeCosLam', axis_spec_bg, axis_spec_cl)
    
    cl_options = {'first_bin_fit_by_slices': h2_clsize.GetXaxis().FindBin(0.5),
                  'last_bin_fit_by_slices': h2_clsize.GetXaxis().FindBin(2.7),
                  'output_dir': cluster_size_dir,
                  'signal_range': [4.5, 14],
                  'simil_bethe_bloch_pars': {'kp1': 0.903,
                                             'kp2': 2.014,
                                             'kp3': 2.440,
                                             'charge': 2.,
                                             'kp4': 2.3,
                                            },
                  'signal_func_name': 'exp_mod_gaus_0',
    }

    ITS_cluster_size = RooRealVar('x', 'x', 0., 15.)
    fitter = Roofitter(ITS_cluster_size, ['exp_mod_gaus'])
    fitter.init_param('exp_mod_gaus_0_mean', 8., 5., 10.)
    fitter.init_param('exp_mod_gaus_0_sigma', 0.5, 0.2, 2.)
    fitter.init_param('exp_mod_gaus_0_tau', 1, -10, 10.)
    cluster_size_pars, its_resolution_pars = cluster_size_calibration(h2_clsize, cluster_size_dir, fitter, charge=2., fit_charge=True, **cl_options)

    simil_bethe_bloch_func = TF1('f_SimilBetheBloch', '([0]/x^[1] + [2]) * [3]^[4]', 0.3, 4.) # function used in cluster_size_calibration
    for ipar, par in enumerate(cluster_size_pars.values()):
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
    dataset['fResolutionITS'] = np_cluster_size_resolution(dataset['fBetaGamma'], *its_resolution_pars.values())
    dataset['fNSigmaITS'] = (dataset['fITSClSizeCosLam'] - dataset['fExpClSizeCosLam']) / (dataset['fResolutionITS']*dataset['fExpClSizeCosLam'])
    h2_clsize_exp = dataset.build_th2('fBetaGamma', 'fExpClSizeCosLam', axis_spec_bg, axis_spec_cl_exp)
    h2_nsigma_its = dataset.build_th2('fInnerParamTPCHe3', 'fNSigmaITS', axis_spec_p, axis_spec_nsigma_its)
    
    # draw results
    cluster_size_dir.cd()
    h2_clsize_exp.Write()
    h2_nsigma_its.Write()  

def cluster_size_calibration_example_he_free(output_file: TFile):
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

    dataset = Dataset.from_root(input_file, folder_name=folder_name, tree_name=tree_name, columns=['fPtHe3', 'fEtaHe3', 'fPhiHe3', 'fChi2TPCHe3', 'fInnerParamTPCHe3', 'fItsClusterSizeHe3', 'fNSigmaTPCHe3', 'fDCAxyHe3', 'fDCAzHe3'])

    # selection cuts from llbariogl/NucleiFlow
    dataset.query('abs(fEtaHe3) < 0.8', inplace=True)
    dataset.query('fPtHe3 < 0', inplace=True)    
    dataset.query('abs(fDCAxyHe3) < 0.1', inplace=True)
    dataset.query('abs(fDCAzHe3) < 1.', inplace=True)
    dataset.query('abs(fNSigmaTPCHe3) < 2.', inplace=True)

    dataset.eval('fPxHe3 = fPtHe3 * cos(fPhiHe3)', inplace=True)
    dataset.eval('fPyHe3 = fPtHe3 * sin(fPhiHe3)', inplace=True)
    dataset.eval('fPzHe3 = fPtHe3 * sinh(fEtaHe3)', inplace=True)
    dataset.eval('fPHe3 = fPtHe3 * cosh(fEtaHe3)', inplace=True)

    dataset['fITSAvgClSize'], __ = average_cluster_size(dataset['fItsClusterSizeHe3'])
    dataset['fITSClSizeCosLam'] = dataset['fITSAvgClSize'] / np.cosh(dataset['fEtaHe3'])

    dataset['fBetaGamma'] = np.abs(dataset['fPHe3']) / 2.8092

    cluster_size_dir = output_file.mkdir('ClusterSizeCalibrationHeFree')

    axis_spec_p = AxisSpec(100, 0., 5, 'h2_nSigma', ';#it{p} (GeV/#it{c}); n#sigma_{TPC}')
    axis_spec_bg = AxisSpec(40, 0., 4, 'h2_dEdxBetaGamma', ';#beta#gamma; Expected dE/dx_{TPC}')
    axis_spec_nsigma_its = AxisSpec(100, -5., 5, 'h2_nSigmaITS', ';#it{p} (GeV/#it{c}); n#sigma_{ITS}')
    axis_spec_cl = AxisSpec(62, -0.5, 15.5, 'h2_ItsClusterSize', ';#beta#gamma; #LT ITS cluster size #GT #times #LT cos#lambda #GT')
    axis_spec_cl_exp = AxisSpec(62, -0.5, 15.5, 'h2_ExpItsClusterSize', ';#beta#gamma; #LT ITS cluster size #GT #times #LT cos#lambda #GT')

    # cluster size calibration
    h2_clsize = dataset.build_th2('fBetaGamma', 'fITSClSizeCosLam', axis_spec_bg, axis_spec_cl)
    
    cl_options = {'first_bin_fit_by_slices': h2_clsize.GetXaxis().FindBin(0.5),
                  'last_bin_fit_by_slices': h2_clsize.GetXaxis().FindBin(2.7),
                  'output_dir': cluster_size_dir,
                  'signal_range': [4.5, 14],
                  'simil_bethe_bloch_pars': {'kp1': 1.225,
                                             'kp2': 1.694,
                                             'kp3': 2.089,
                                             'charge': 1.,
                                             'kp4': 1.,
                                            },
                  'signal_func_name': 'exp_mod_gaus_0',
    }

    ITS_cluster_size = RooRealVar('x', 'x', 0., 15.)
    fitter = Roofitter(ITS_cluster_size, ['exp_mod_gaus'])
    fitter.init_param('exp_mod_gaus_0_mean', 8., 5., 10.)
    fitter.init_param('exp_mod_gaus_0_sigma', 0.5, 0.2, 2.)
    fitter.init_param('exp_mod_gaus_0_tau', 1, -10, 10.)
    cluster_size_pars, its_resolution_pars = cluster_size_calibration(h2_clsize, cluster_size_dir, fitter, charge=1., fit_charge=False, **cl_options)

    simil_bethe_bloch_func = TF1('f_SimilBetheBloch', '([0]/x^[1] + [2]) * [3]^[4]', 0.3, 4.) # function used in cluster_size_calibration
    for ipar, par in enumerate(cluster_size_pars.values()):
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
    dataset['fResolutionITS'] = np_cluster_size_resolution(dataset['fBetaGamma'], *its_resolution_pars.values())
    dataset['fNSigmaITS'] = (dataset['fITSClSizeCosLam'] - dataset['fExpClSizeCosLam']) / (dataset['fResolutionITS']*dataset['fExpClSizeCosLam'])
    h2_clsize_exp = dataset.build_th2('fBetaGamma', 'fExpClSizeCosLam', axis_spec_bg, axis_spec_cl_exp)
    h2_nsigma_its = dataset.build_th2('fPHe3', 'fNSigmaITS', axis_spec_p, axis_spec_nsigma_its)
    
    # draw results
    cluster_size_dir.cd()
    h2_clsize_exp.Write()
    h2_nsigma_its.Write()  

def cluster_size_cut_example_he(output_file: TFile):
    '''
        Demonstration of the cluster size cut implementation and its effect on the data.
    '''

    input_file = 'data/LHC23_pass4_skimmed_he_pid.root'    
                                            # required tree from O2Physics table that contains the following branches: 
                                            # fEtaPr, fChi2TPCPr, fInnerParamTPCPr, fItsClusterSizePr
    folder_name = 'DF*'                     # name of the folder(s) containing the tree (O2Physics table)
    tree_name = 'O2clsttableextra'             # name of the tree (O2Physics table) in the input file

    dataset = Dataset.from_root(input_file, folder_name=folder_name, tree_name=tree_name, columns=['fP', 'fEta', 'fPhi', 'fItsClusterSize', 'fPartID', 'fTpcNSigma'])

    # selection cuts from llbariogl/NucleiFlow
    dataset.query('abs(fEta) < 0.8', inplace=True)
    dataset.query('fP < 0', inplace=True)    
    
    dataset['fITSAvgClSize'], __ = average_cluster_size(dataset['fItsClusterSize'])
    dataset['fITSClSizeCosLam'] = dataset['fITSAvgClSize'] / np.cosh(dataset['fEta'])

    dataset['fP'] = dataset['fP'] * 2.
    dataset['fBetaGamma'] = np.abs(dataset['fP']) / 2.8092

    cluster_size_dir = output_file.mkdir('ClusterSizeCutHe')

    axis_spec_p = AxisSpec(100, 0., 5, 'h2_nSigma', ';#it{p} (GeV/#it{c}); n#sigma_{TPC}')
    axis_spec_bg = AxisSpec(40, 0., 4, 'h2_dEdxBetaGamma', ';#beta#gamma; Expected dE/dx_{TPC}')
    axis_spec_nsigma_its = AxisSpec(100, -5., 5, 'h2_nSigmaITS', ';#it{p} (GeV/#it{c}); n#sigma_{ITS}')
    axis_spec_nsigma_tpc = AxisSpec(100, -5., 5, 'h_nSigmaTPC', '; n#sigma_{TPC}')
    axis_spec_cl = AxisSpec(62, -0.5, 15.5, 'h2_ItsClusterSize', ';#beta#gamma; #LT ITS cluster size #GT #times #LT cos#lambda #GT')
    axis_spec_cl_exp = AxisSpec(62, -0.5, 15.5, 'h2_ExpItsClusterSize', ';#beta#gamma; #LT ITS cluster size #GT #times #LT cos#lambda #GT')

    # cluster size calibration
    h2_clsize_before_cut = dataset.build_th2('fBetaGamma', 'fITSClSizeCosLam', axis_spec_bg, axis_spec_cl)
    h_nsigma_tpc_before_cut = dataset.build_th2('fTpcNSigma', axis_spec_nsigma_tpc)
    h_nsigma_tpc_before_cut.SetName('h_nSigmaTPC_before_cut')

    cluster_size_pars = {
        'kp1': 2.8752,
        'kp2': 1.1246,
        'kp3': 5.0259,
        'charge': 1.,
        'kp4': 1.,
    }
    its_resolution_pars = {
        'rp0': 0.2431,
        'rp1': -0.3293,
        'rp2': 1.533,
    }

    simil_bethe_bloch_func = TF1('f_SimilBetheBloch', '([0]/x^[1] + [2]) * [3]^[4]', 0.3, 4.) # function used in cluster_size_calibration
    for ipar, par in enumerate(cluster_size_pars.values()):
        simil_bethe_bloch_func.SetParameter(ipar, par)
    c_clsize = TCanvas('c_ClSizeAndCurve', 'canvas', 800, 600)
    h2_clsize_before_cut.Draw('colz')
    simil_bethe_bloch_func.Draw('same')

    # draw results
    cluster_size_dir.cd()
    h2_clsize_before_cut.Write('h2_clsize_before_cut')
    simil_bethe_bloch_func.Write()
    h_nsigma_tpc_before_cut.Write()
    c_clsize.Write()

    dataset['fITSResolution'] = np_cluster_size_resolution(dataset['fBetaGamma'], *its_resolution_pars.values())
    dataset['fExpClSizeCosLam'] = cluster_size_parametrisation(dataset['fBetaGamma'], *cluster_size_pars.values())
    dataset['fNSigmaITS'] = (dataset['fITSClSizeCosLam'] - dataset['fExpClSizeCosLam']) / (dataset['fITSResolution']*dataset['fExpClSizeCosLam'])
    h2_clsize_exp = dataset.build_th2('fBetaGamma', 'fExpClSizeCosLam', axis_spec_bg, axis_spec_cl_exp)
    h2_nsigma_its = dataset.build_th2('fP', 'fNSigmaITS', axis_spec_p, axis_spec_nsigma_its)

    dataset.query('fNSigmaITS > -1.5', inplace=True)
    h2_clsize_after_cut = dataset.build_th2('fBetaGamma', 'fITSClSizeCosLam', axis_spec_bg, axis_spec_cl)
    h_nsigma_tpc_after_cut = dataset.build_th2('fTpcNSigma', axis_spec_nsigma_tpc)
    
    # draw results
    cluster_size_dir.cd()
    h2_clsize_exp.Write()
    h2_nsigma_its.Write()
    h2_clsize_after_cut.Write('h2_clsize_after_cut')  
    h_nsigma_tpc_after_cut.Write('h_nSigmaTPC_after_cut')

if __name__ == '__main__':
    output_file_path = 'data/cluster_size_calibration_example_newest.root'
    output_file = TFile(output_file_path, 'recreate')
    
    #cluster_size_calibration_example_proton(output_file)
    #cluster_size_cut_example(output_file)
    cluster_size_calibration_example_proton_pid(output_file)
    cluster_size_calibration_example_proton_pid_mc(output_file)
    #cluster_size_calibration_example_he_free(output_file)
    #cluster_size_cut_example_he(output_file)
    #cluster_size_calibration_example_he_free(output_file)

    output_file.Close()
