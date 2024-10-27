import os 
import numpy as np
from ROOT import TFile, TF1, TCanvas

from torchic.physics.calibration import bethe_bloch_calibration, pyBetheBloch, pySimilBetheBloch, cluster_size_calibration
from torchic.physics.ITS import average_cluster_size
from torchic.core.dataset import Dataset
from torchic.core.histogram import AxisSpec

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

    dataset = Dataset(input_file, folder_name=folder_name, tree_name=tree_name, columns=['fPt', 'fTPCsignal', 'fTPCInnerParam', 'fDCAxy', 'fDCAz', 'fTPCnCls', 'fZvertex', 'fITSclusterSizes', 'fEta'])
    
    # selection cuts from llbariogl/NucleiFlow
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
    cluster_size_dir = output_file.mkdir('ClusterSizeCalibration')

    axis_spec_p = AxisSpec(100, 0., 5, 'h2_nSigma', ';#it{p}_{TPC} (GeV/#it{c}); n#sigma_{TPC}')
    axis_spec_bg = AxisSpec(40, 0., 4, 'h2_dEdxBetaGamma', ';#beta#gamma; Expected dE/dx_{TPC}')
    axis_spec_nsigma_tpc = AxisSpec(100, -5., 5, 'h2_nSigmaTPC', ';#it{p}_{TPC} (GeV/#it{c}); n#sigma_{TPC}')
    axis_spec_nsigma_its = AxisSpec(100, -5., 5, 'h2_nSigmaITS', ';#it{p}_{TPC} (GeV/#it{c}); n#sigma_{ITS}')
    axis_spec_dEdx = AxisSpec(300, 0., 3000., 'h2_dEdxBetaGamma', ';#beta#gamma; Expected dE/dx_{TPC}')
    axis_spec_cl = AxisSpec(62, -0.5, 15.5, 'h2_ItsClusterSize', ';#beta#gamma; #LT ITS cluster size #GT #times #LT cos#lambda #GT')
    axis_spec_cl_exp = AxisSpec(62, -0.5, 15.5, 'h2_ExpItsClusterSize', ';#beta#gamma; #LT ITS cluster size #GT #times #LT cos#lambda #GT')

    # cluster size calibration
    h2_clsize = dataset.build_hist('fBetaGamma', 'fITSClSizeCosLam', axis_spec_bg, axis_spec_cl)
    
    cl_options = {'first_bin_fit_by_slices': h2_clsize.GetXaxis().FindBin(0.5),
                  'last_bin_fit_by_slices': h2_clsize.GetXaxis().FindBin(3.2),
                  'save_slices': True,
                  'fit_range': [1., 11.],
    }
    cluster_size_pars, its_resolution = cluster_size_calibration(h2_clsize, cluster_size_dir, **cl_options)

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

    dataset['fExpClSizeCosLam'] = pySimilBetheBloch(dataset['fBetaGamma'], *cluster_size_pars.values())
    dataset['fNSigmaITS'] = (dataset['fITSClSizeCosLam'] - dataset['fExpClSizeCosLam']) / (its_resolution*dataset['fExpClSizeCosLam'])
    h2_clsize_exp = dataset.build_hist('fBetaGamma', 'fExpClSizeCosLam', axis_spec_bg, axis_spec_cl_exp)
    h2_nsigma_its = dataset.build_hist('fTPCInnerParam', 'fNSigmaITS', axis_spec_p, axis_spec_nsigma_its)
    h2_dEdx_no_clsize_cut = dataset.build_hist('fBetaGamma', 'fTPCsignal', axis_spec_bg, axis_spec_dEdx)
    h2_dEdx_no_clsize_cut.SetName('h2_dEdxBetaGammaNoClSizeCut')
    dataset.query('fNSigmaITS > -1', inplace=True)
    
    # energy loss calibration
    h2_dEdx = dataset.build_hist('fBetaGamma', 'fTPCsignal', axis_spec_bg, axis_spec_dEdx)
    
    bb_options = {'first_bin_fit_by_slices': h2_dEdx.GetXaxis().FindBin(0.55),
               'last_bin_fit_by_slices': h2_dEdx.GetXaxis().FindBin(4.)
    }
    bethe_bloch_pars = bethe_bloch_calibration(h2_dEdx, bethe_bloch_dir, **bb_options)

    dataset['fExpSignalTPC'] = pyBetheBloch(dataset['fBetaGamma'], *bethe_bloch_pars.values())
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

