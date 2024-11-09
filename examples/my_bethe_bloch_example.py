import os 
import numpy as np
from ROOT import TFile, TF1, gInterpreter, TCanvas
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BETHEBLOCH_DIR = os.path.join(CURRENT_DIR, '..', 'torchic', 'physics', 'BetheBloch.hh')
gInterpreter.ProcessLine(f'#include "{BETHEBLOCH_DIR}"')
from ROOT import BetheBloch

from torchic.physics.calibration import bethe_bloch_calibration, py_BetheBloch, cluster_size_parametrisation, DEFAULT_BETHEBLOCH_PARS, cluster_size_calibration
from torchic.physics.ITS import average_cluster_size
from torchic.core.dataset import Dataset
from torchic.core.histogram import load_hist, HistLoadInfo, AxisSpec

def bethe_bloch_calibration_example():

    hist_info = HistLoadInfo('data/data_visual_selectionsHe3.root', 'TPC/dEdXvsBetaGammaHe3')
    h2 = load_hist(hist_info)
    h2.RebinX(2)

    output_file_path = 'data/bethe_bloch_calibration.root'
    output_file = TFile(output_file_path, 'recreate')
    bethe_bloch_dir = output_file.mkdir('BetheBlochCalibration')

    options = {'first_bin_fit_by_slices': h2.GetXaxis().FindBin(0.3),
               'last_bin_fit_by_slices': h2.GetXaxis().FindBin(2.)
    }
    bethe_bloch_pars = bethe_bloch_calibration(h2, bethe_bloch_dir, **options)

    # draw nsigma distribution
    dataset = Dataset('data/LHC23zzh_pass4_small_mixing.root', tree_name='O2lithium4table', folder_name='DF*', columns=['fSignalTPCHe3', 'fInnerParamTPCHe3'])
    dataset['fInnerParamTPCHe3'] = dataset['fInnerParamTPCHe3'] * 2.
    dataset['fBetaGammaHe3'] = dataset['fInnerParamTPCHe3'] / 2.80923
    dataset['fExpSignalTPCHe3'] = py_BetheBloch(dataset['fBetaGammaHe3'], *bethe_bloch_pars.values())
    dataset['fNSigmaTPCHe3'] = (dataset['fSignalTPCHe3'] - dataset['fExpSignalTPCHe3']) / (0.09*dataset['fExpSignalTPCHe3'])

    #dataset.query('fBetaGammaHe3 > 0.5 & fBetaGammaHe3 < 1.5', inplace=True)
    #print(dataset.head(10))

    axis_spec_p = AxisSpec(100, 0., 5, 'h2_nSigma', ';#it{p}_{TPC} (GeV/#it{c}); n#sigma_{TPC}')
    axis_spec_bg = AxisSpec(80, 0., 4, 'h2_dEdxBetaGamma', ';#beta#gamma; Expected dE/dx_{TPC}')
    axis_spec_nsigma = AxisSpec(100, -5., 5, 'h2_nSigma', ';#it{p}_{TPC} (GeV/#it{c}); n#sigma_{TPC}')
    axis_spec_dEdx = AxisSpec(750, 0., 3000., 'h2_dEdxBetaGamma', ';#beta#gamma; Expected dE/dx_{TPC}')
    h2_nsigma = dataset.build_hist('fInnerParamTPCHe3', 'fNSigmaTPCHe3', axis_spec_p, axis_spec_nsigma)
    h2_bg = dataset.build_hist('fBetaGammaHe3', 'fExpSignalTPCHe3', axis_spec_bg, axis_spec_dEdx)

    bethe_bloch_func = TF1('bethe_bloch_func', BetheBloch, 0.3, 2., 5)
    bethe_bloch_func.SetParameters(*DEFAULT_BETHEBLOCH_PARS.values())
    bethe_bloch_func.SetParNames(*DEFAULT_BETHEBLOCH_PARS.keys())
    c_old_params_BB = TCanvas('c_OldParamsBetheBloch', 'canvas', 800, 600)
    h2_bg.Draw('colz')
    bethe_bloch_func.Draw('same')

    bethe_bloch_dir.cd()
    h2_nsigma.Write()
    h2_bg.Write()
    c_old_params_BB.Write()
    
    output_file.Close()

def bethe_bloch_calibration_lbariogl_example():

    dataset = Dataset('data/AO2D_merged.root', tree_name='O2nucleitable', folder_name='DF*', columns=['fPt', 'fTPCsignal', 'fTPCInnerParam', 'fDCAxy', 'fDCAz', 'fTPCnCls', 'fZvertex', 'fITSclusterSizes', 'fEta'])
    
    # selection cuts from llbariogl/NucleiFlow
    dataset.query('abs(fDCAxy) < 0.1', inplace=True)
    dataset.query('abs(fDCAz) < 1.', inplace=True)
    dataset.query('fTPCnCls > 99', inplace=True)
    dataset.query('abs(fZvertex) < 10', inplace=True)

    dataset['fITSAvgClSize'], __ = average_cluster_size(dataset['fITSclusterSizes'])
    dataset['fITSClSizeCosLam'] = dataset['fITSAvgClSize'] / np.cosh(dataset['fEta'])
    #dataset.query('(fPt < 2.4 and fITSClSizeCosLam > 6) or (2.4 < fPt < 2.8 and fITSClSizeCosLam > 5.5) or (2.8 < fPt < 100 and fITSClSizeCosLam > 5)', inplace=True)

    dataset['fTPCInnerParam'] = dataset['fTPCInnerParam'] * 2.
    dataset['fBetaGamma'] = dataset['fTPCInnerParam'] / 2.80923

    output_file_path = 'data/bethe_bloch_calibration_lbariogl.root'
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

    simil_bethe_bloch_func = TF1('f_SimilBetheBloch', '[kp1]/x^2 + [kp2]', 0.3, 4.)
    for ipar, par in cluster_size_pars.items():
        simil_bethe_bloch_func.SetParameter(ipar, par)
    c_clsize = TCanvas('c_ClSizeAndCurve', 'canvas', 800, 600)
    h2_clsize.Draw('colz')
    simil_bethe_bloch_func.Draw('same')

    cluster_size_dir.cd()
    h2_clsize.Write()
    simil_bethe_bloch_func.Write()
    c_clsize.Write()

    dataset['fExpClSizeCosLam'] = cluster_size_parametrisation(dataset['fBetaGamma'], *cluster_size_pars.values())
    dataset['fNSigmaITS'] = (dataset['fITSClSizeCosLam'] - dataset['fExpClSizeCosLam']) / (its_resolution*dataset['fExpClSizeCosLam'])
    h2_clsize_exp = dataset.build_hist('fBetaGamma', 'fExpClSizeCosLam', axis_spec_bg, axis_spec_cl_exp)
    h2_nsigma_its = dataset.build_hist('fTPCInnerParam', 'fNSigmaITS', axis_spec_p, axis_spec_nsigma_its)
    h2_dEdx_old = dataset.build_hist('fBetaGamma', 'fTPCsignal', axis_spec_bg, axis_spec_dEdx)
    h2_dEdx_old.SetName('h2_dEdxBetaGammaOld')
    dataset.query('fNSigmaITS > -1', inplace=True)
    
    # energy loss calibration
    h2_dEdx = dataset.build_hist('fBetaGamma', 'fTPCsignal', axis_spec_bg, axis_spec_dEdx)
    
    bb_options = {'first_bin_fit_by_slices': h2_dEdx.GetXaxis().FindBin(0.55),
               'last_bin_fit_by_slices': h2_dEdx.GetXaxis().FindBin(4.)
    }
    bethe_bloch_pars = bethe_bloch_calibration(h2_dEdx, bethe_bloch_dir, **bb_options)

    dataset['fExpSignalTPC'] = py_BetheBloch(dataset['fBetaGamma'], *bethe_bloch_pars.values())
    dataset['fNSigmaTPC'] = (dataset['fTPCsignal'] - dataset['fExpSignalTPC']) / (0.09*dataset['fExpSignalTPC'])
    h2_bg = dataset.build_hist('fBetaGamma', 'fExpSignalTPC', axis_spec_bg, axis_spec_dEdx)
    h2_nsigma_tpc = dataset.build_hist('fTPCInnerParam', 'fNSigmaTPC', axis_spec_p, axis_spec_nsigma_tpc)

    bethe_bloch_func_old = TF1('bethe_bloch_func_old', BetheBloch, 0.3, 4., 5)
    bethe_bloch_func_old.SetParameters(*DEFAULT_BETHEBLOCH_PARS.values())
    bethe_bloch_func_old.SetParNames(*DEFAULT_BETHEBLOCH_PARS.keys())
    bethe_bloch_func_new = TF1('bethe_bloch_func_new', BetheBloch, 0.3, 4., 5)
    bethe_bloch_func_new.SetParameters(*bethe_bloch_pars.values())
    bethe_bloch_func_new.SetParNames(*bethe_bloch_pars.keys())
    c_old_params_BB = TCanvas('c_OldParamsBetheBloch', 'canvas', 800, 600)
    h2_dEdx.Draw('colz')
    bethe_bloch_func_old.Draw('same')
    bethe_bloch_func_new.SetLineColor(3)
    bethe_bloch_func_new.Draw('same')

    cluster_size_dir.cd()
    h2_clsize_exp.Write()
    h2_nsigma_its.Write()

    bethe_bloch_dir.cd()
    h2_dEdx_old.Write()
    h2_nsigma_tpc.Write()
    h2_bg.Write()
    c_old_params_BB.Write()
    
    output_file.Close()
    

if __name__ == '__main__':
    bethe_bloch_calibration_lbariogl_example()

