import numpy as np
from ROOT import TFile, TF1, TCanvas, RooRealVar

from deprecated.torchic.physics.calibration import bethe_bloch_calibration, py_BetheBloch, cluster_size_parametrisation, cluster_size_calibration
from torchic.physics.ITS import average_cluster_size
from torchic import Dataset, AxisSpec
from deprecated.torchic.core.roofitter import Roofitter

def getSign(flags):
    if flags & 256:
        return 1
    else:
        return -1


# vectorised version of getSign
getSign_vectorised = np.vectorize(getSign)

def bethe_bloch_calibration_example():
    '''
        The following function performs a TPC calibration (of He3) with the Bethe-Bloch formula.
        The calibration is performed on the dE/dx vs. beta*gamma plot and the data is retrieved from a tree (O2Physics table).
        The data is preselected with several cuts, notably a ITS cluster size cut.
    '''

    input_file = 'data/AO2D_merged.root'    # required tree from O2Physics table that contains the following branches: 
                                            # fPt, fTPCsignal, fTPCInnerParam, fDCAxy, fDCAz, fTPCnCls, fZvertex, fITSclusterSizes, fEta
    folder_name = 'DF*'                     # name of the folder(s) containing the tree (O2Physics table)
    tree_name = 'O2nucleitable'             # name of the tree (O2Physics table) in the input file

    dataset = Dataset.from_root(input_file, folder_name=folder_name, tree_name=tree_name, columns=['fPt', 'fTPCsignal', 'fTPCInnerParam', 'fDCAxy', 'fDCAz', 'fTPCnCls', 'fZvertex', 'fITSclusterSizes', 'fEta', 'fFlags'])
    dataset['fSign'] = getSign_vectorised(dataset['fFlags'])

    # selection cuts from llbariogl/NucleiFlow
    dataset.query('abs(fEta) < 0.8', inplace=True)
    dataset.query('fSign < 0', inplace=True)    
    dataset.query('abs(fDCAxy) < 0.1', inplace=True)
    dataset.query('abs(fDCAz) < 1.', inplace=True)
    dataset.query('fTPCnCls > 99', inplace=True)
    dataset.query('abs(fZvertex) < 10', inplace=True)

    dataset['fITSAvgClSize'], __ = average_cluster_size(dataset['fITSclusterSizes'])
    dataset['fITSClSizeCosLam'] = dataset['fITSAvgClSize'] / np.cosh(dataset['fEta'])

    dataset['fTPCInnerParam'] = dataset['fTPCInnerParam'] * 2.
    dataset['fBetaGamma'] = dataset['fTPCInnerParam'] / 2.80923

    output_file_path = 'data/bethe_bloch_calibration_example.root'
    output_file = TFile(output_file_path, 'recreate')
    bethe_bloch_dir = output_file.mkdir('BetheBlochCalibration')
    bethe_bloch_slices_dir = bethe_bloch_dir.mkdir('Slices')
    cluster_size_dir = output_file.mkdir('ClusterSizeCalibration')
    cluster_size_slices_dir = cluster_size_dir.mkdir('Slices')

    axis_spec_p = AxisSpec(100, 0., 5, 'h2_nSigma', ';#it{p}_{TPC} (GeV/#it{c}); n#sigma_{TPC}')
    axis_spec_bg = AxisSpec(40, 0., 4, 'h2_dEdxBetaGamma', ';#beta#gamma; Expected dE/dx_{TPC}')
    axis_spec_nsigma_tpc = AxisSpec(100, -5., 5, 'h2_nSigmaTPC', ';#it{p}_{TPC} (GeV/#it{c}); n#sigma_{TPC}')
    axis_spec_nsigma_its = AxisSpec(100, -5., 5, 'h2_nSigmaITS', ';#it{p}_{TPC} (GeV/#it{c}); n#sigma_{ITS}')
    axis_spec_dEdx = AxisSpec(300, 0., 3000., 'h2_dEdxBetaGamma', ';#beta#gamma; Expected dE/dx_{TPC}')
    axis_spec_cl = AxisSpec(62, -0.5, 15.5, 'h2_ItsClusterSize', ';#beta#gamma; #LT ITS cluster size #GT #times #LT cos#lambda #GT')
    axis_spec_cl_exp = AxisSpec(62, -0.5, 15.5, 'h2_ExpItsClusterSize', ';#beta#gamma; #LT ITS cluster size #GT #times #LT cos#lambda #GT')

    # cluster size calibration
    h2_clsize = dataset.build_hist('fBetaGamma', 'fITSClSizeCosLam', axis_spec_bg, axis_spec_cl)
    
    cl_options = {'first_bin_fit_by_slices': h2_clsize.GetXaxis().FindBin(0.6),
                  'last_bin_fit_by_slices': h2_clsize.GetXaxis().FindBin(2.8),
                  'output_dir': cluster_size_slices_dir,
                  'signal_func_name': 'exp_mod_gaus_0',
    }
    ITS_cluster_size = RooRealVar('x', 'x', 1., 11.)
    fitter = Roofitter(ITS_cluster_size, ['exp_mod_gaus', 'gaus'])
    fitter.init_param('exp_mod_gaus_0_mean', 8., 5., 10.)
    fitter.init_param('exp_mod_gaus_0_sigma', 0.5, 0.2, 2.)
    fitter.init_param('gaus_1_mean', 2.5, 1.5, 3.5)
    fitter.init_param('gaus_1_sigma', 0.5, 0.2, 1.5)
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
    h2_nsigma_its = dataset.build_hist('fTPCInnerParam', 'fNSigmaITS', axis_spec_p, axis_spec_nsigma_its)
    h2_dEdx_no_clsize_cut = dataset.build_hist('fBetaGamma', 'fTPCsignal', axis_spec_bg, axis_spec_dEdx)
    h2_dEdx_no_clsize_cut.SetName('h2_dEdxBetaGammaNoClSizeCut')
    dataset.query('fNSigmaITS > 0', inplace=True)
    
    # energy loss calibration
    h2_dEdx = dataset.build_hist('fBetaGamma', 'fTPCsignal', axis_spec_bg, axis_spec_dEdx)
    
    bb_options = {'first_bin_fit_by_slices': h2_dEdx.GetXaxis().FindBin(0.55),
                  'last_bin_fit_by_slices': h2_dEdx.GetXaxis().FindBin(4.),
                  'output_dir': bethe_bloch_slices_dir,
                  'signal_func_name': 'gaus_0',
                  'signal_range': (0., 3000.),
                 }
    dEdx = RooRealVar('x', 'x', 0., 3000.)
    fitter = Roofitter(dEdx, ['gaus'])
    fitter.init_param('gaus_0_mean', 500., 200., 800.)
    fitter.init_param('gaus_0_sigma', 50, 0., 500.)

    bethe_bloch_pars = bethe_bloch_calibration(h2_dEdx, bethe_bloch_dir, fitter, **bb_options)

    dataset['fExpSignalTPC'] = py_BetheBloch(dataset['fBetaGamma'], *bethe_bloch_pars.values())
    dataset['fNSigmaTPC'] = (dataset['fTPCsignal'] - dataset['fExpSignalTPC']) / (0.09*dataset['fExpSignalTPC'])
    h2_bg = dataset.build_hist('fBetaGamma', 'fExpSignalTPC', axis_spec_bg, axis_spec_dEdx)
    h2_nsigma_tpc = dataset.build_hist('fTPCInnerParam', 'fNSigmaTPC', axis_spec_p, axis_spec_nsigma_tpc)

    # draw results
    cluster_size_dir.cd()
    h2_clsize_exp.Write()
    h2_nsigma_its.Write()

    bethe_bloch_dir.cd()
    h2_dEdx_no_clsize_cut.Write()
    h2_nsigma_tpc.Write()
    h2_bg.Write()
    
    output_file.Close()
    

if __name__ == '__main__':
    bethe_bloch_calibration_example()

