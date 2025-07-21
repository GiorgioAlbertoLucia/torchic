import numpy as np
from ROOT import TDirectory, TFile

from torchic.core.dataset import Dataset
from torchic.core.histogram import AxisSpec
from torchic.physics.ITS import average_cluster_size
from deprecated.torchic.physics.calibration import bethe_bloch_calibration

def compute_efficiency(data: Dataset, output_file: TDirectory, **kwargs):

    HE3_MASS = 2.80923
    PROTON_MASS = 0.938272

    rec_data = None

    if 'particle' in kwargs and kwargs['particle'] == 'Li':

        data['fSignedPtPr'] = data['fPtPr']
        data['fSignedPtHe3'] = data['fPtHe3']
        data['fPtPr'] = np.abs(data['fPtPr'])   
        data['fPtHe3'] = np.abs(data['fPtHe3'])
        data['fSignHe3'] = data['fSignedPtHe3'] / data['fPtHe3']
        
        data['fPxHe3'] = data['fPtHe3'] * np.cos(data['fPhiHe3'])
        data['fPxPr'] = data['fPtPr'] * np.cos(data['fPhiPr'])
        data['fPyHe3'] = data['fPtHe3'] * np.sin(data['fPhiHe3'])
        data['fPyPr'] = data['fPtPr'] * np.sin(data['fPhiPr'])
        data['fPzHe3'] = data['fPtHe3'] * np.sinh(data['fEtaHe3'])
        data['fPzPr'] = data['fPtPr'] * np.sinh(data['fEtaPr'])
        data['fPHe3'] = data['fPtHe3'] * np.cosh(data['fEtaHe3'])
        data['fPPr'] = data['fPtPr'] * np.cosh(data['fEtaPr'])
        data['fEHe3'] = np.sqrt(data['fPHe3']**2 + HE3_MASS**2)
        data['fEPr'] = np.sqrt(data['fPPr']**2 + PROTON_MASS**2)
        data['fDeltaPhi'] = data['fPhiHe3'] - data['fPhiPr']

        # invariant mass 
        data['fPt'] = np.sqrt(data['fPtHe3']**2 + data['fPtPr']**2 + 2*data['fPtHe3']*data['fPtPr']*np.cos(data['fDeltaPhi'])) * data['fSignHe3']
        data['fPtMC'] = data['fSignedPtMC']
        data['fMassInvLi'] = np.sqrt( (data['fEHe3'] + data['fEPr'])**2 - ((data['fPxHe3'] + data['fPxPr'])**2 + (data['fPyHe3'] + data['fPyPr'])**2 + (data['fPzHe3'] + data['fPzPr'])**2) )

        rec_data = data.query('fPt < 99', inplace=False)
    else:
        rec_data = data.query('fPt > -990', inplace=False)
    if 'particle' in kwargs:
        if kwargs['particle'] == 'He':
            rec_data.query('fInnerParamTPC > 0.8', inplace=True)
            rec_data.query('-2. < fNSigmaTPC < 2.', inplace=True)
    gen_data = data

    if kwargs.get('with_selections', False) and 'particle' in kwargs:
        if kwargs['particle'] == 'He':
            pass
        if kwargs['particle'] == 'Pr':
            rec_data.query('-3. < fNSigmaTPC < 3.', inplace=True)
            #rec_data.query('fNSigmaTOF < -990. or -3 < fNSigmaTOF < 3', inplace=True)

    # do it in two steps for matter and antimatter
    conditions_dict = {
                       'antimatter': 'fPtMC < 0',
                       'matter': 'fPtMC > 0', 
                       }

    for key, condition in conditions_dict.items():
        tmp_rec_data = rec_data.query(condition, inplace=False)
        tmp_gen_data = gen_data.query(condition, inplace=False)

        tmp_rec_data['fPt'] = np.abs(tmp_rec_data['fPt'])
        tmp_gen_data['fPtMC'] = np.abs(tmp_gen_data['fPtMC'])

        axis_spec_pt_rec = AxisSpec(50, 0, 10, 'pt_rec', '#it{p}_{T}^{rec} (GeV/#it{c})')
        axis_spec_pt_gen = AxisSpec(50, 0, 10, 'pt_gen', '#it{p}_{T}^{gen} (GeV/#it{c})')

        h_rec = tmp_rec_data.build_hist('fPt', axis_spec_pt_rec)
        h_gen = tmp_gen_data.build_hist('fPtMC', axis_spec_pt_gen)
        h_efficiency = h_rec.Clone('hEff')
        h_efficiency.Reset()

        for ibin in range(1, h_rec.GetNbinsX() + 1):
            rec = h_rec.GetBinContent(ibin)
            gen = h_gen.GetBinContent(ibin)
            eff, eff_err = 0, 0
            if gen > 0:
                eff = rec / gen
                eff_err = np.sqrt(eff * (1 - eff) / gen) if eff < 1 else 0
            h_efficiency.SetBinContent(ibin, eff) 
            h_efficiency.SetBinError(ibin, eff_err)

        output_file.cd()
        suffix = '_selections' if kwargs.get('with_selections', False) else ''
        h_efficiency.Write(f'hEfficiency_{key}{suffix}')

        del tmp_rec_data, tmp_gen_data, h_rec, h_gen, h_efficiency

def tpc_calibration_proton(data: Dataset, output_file: TDirectory, particle: str):

    axis_spec_dedx = AxisSpec(200, 0, 2000, 'h2_dEdx', ';#beta#gamma;dE/dx (a.u.)')
    axis_spec_bg = AxisSpec(50, 0, 5, '', '#beta#gamma')

    rec_data = data.query('fPt > -990', inplace=False)

    PROTON_MASS = 0.938272
    HE3_MASS = 2.80923
    if particle == 'Pr': 
        particle_mass = PROTON_MASS
    elif particle == 'He':
        particle_mass = HE3_MASS
        rec_data['fInnerParamTPC'] = rec_data['fInnerParamTPC'] * 2
    else:
        raise ValueError('Invalid particle type')
    
    rec_data['fBetaGamma'] = rec_data['fInnerParamTPC'] / particle_mass

    h2_dEdx = rec_data.build_hist('fBetaGamma', 'fSignalTPC', axis_spec_bg, axis_spec_dedx)
    if particle == 'Pr':    
        bb_options = {'first_bin_fit_by_slices': h2_dEdx.GetXaxis().FindBin(0.2),
               'last_bin_fit_by_slices': h2_dEdx.GetXaxis().FindBin(4.)
    }
    elif particle == 'He':
        bb_options = {'first_bin_fit_by_slices': h2_dEdx.GetXaxis().FindBin(0.5),
               'last_bin_fit_by_slices': h2_dEdx.GetXaxis().FindBin(1.8)
    }
    
    bethe_bloch_pars = bethe_bloch_calibration(h2_dEdx, output_file, **bb_options)

    return bethe_bloch_pars

def data_visual(data: Dataset, output_file: TDirectory):

    axis_spec_pt_rec = AxisSpec(100, -10, 10, 'h_PtRec', '#it{p}_{T}^{rec} (GeV/#it{c})')
    axis_spec_pt_gen = AxisSpec(100, -10, 10, 'h_PtGen', '#it{p}_{T}^{gen} (GeV/#it{c})')
    axis_spec_its = AxisSpec(62, -0.5, 15.5, 'h2_ItsClusterSize', ';#beta#gamma; #LT ITS cluster size #GT #times #LT cos#lambda #GT')
    axis_spec_nsigma_tpc = AxisSpec(100, -5, 5, 'h2_NSigmaTPC', ';#it{p}_{T}^{rec} (GeV/#it{c});n#sigma_{TPC}')

    rec_data = data.query('fPt > -990', inplace=False)
    rec_data['fItsAvgClSize'], __ = average_cluster_size(rec_data['fItsClusterSize'])
    rec_data['fItsClSizeCosLam'] = rec_data['fItsAvgClSize'] / np.cosh(rec_data['fEta'])
    
    h_pt_gen = data.build_hist('fPtMC', axis_spec_pt_gen)
    h_pt_rec = rec_data.build_hist('fPt', axis_spec_pt_rec)
    h2_its_pt = rec_data.build_hist('fPt', 'fItsClSizeCosLam', axis_spec_pt_rec, axis_spec_its)
    h2_nsigma_tpc = rec_data.build_hist('fPt', 'fNSigmaTPC', axis_spec_pt_rec, axis_spec_nsigma_tpc)

    output_file.cd()
    h_pt_gen.Write()
    h_pt_rec.Write()
    h2_its_pt.Write()
    h2_nsigma_tpc.Write()
    

if __name__ == '__main__':

    input_files = {'He': ['/home/galucia/antiLithium4/task/MCWorkflowFindables/output/MC_efficiency_he.root'],
                   'Pr': ['/home/galucia/antiLithium4/task/MCWorkflowFindables/output/MC_efficiency_pr.root'],
                   'Li': ['/home/galucia/antiLithium4/task/MCWorkflowFindables/output/MC_efficiency_li.root']
                   }

    output_file_path = '/home/galucia/antiLithium4/analysis/output/MC/mc.root'
    output_file = TFile(output_file_path, 'recreate')

    for part in ['He', 'Pr', 'Li']:

        tree_name = 'O2lithium4tablemc' if part == 'Li' else 'O2lithium4findmc'
        folder_name = 'DF*'
        data = Dataset(input_files[part], tree_name=tree_name, folder_name=folder_name)

        tdir = output_file.mkdir(part)
        if part == 'Pr' or part == 'He':
            tpc_calibration_proton(data, tdir, particle=part)

        compute_efficiency(data, tdir, particle=part)
        compute_efficiency(data, tdir, particle=part, with_selections=True)
        
        if part != 'Li':
            data_visual(data, tdir)
    
    output_file.Close()