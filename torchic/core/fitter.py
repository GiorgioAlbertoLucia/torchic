'''

'''

import re
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from ROOT import TF1, TH1F, TH2F, TCanvas

from torchic.utils.overload import overload, signature

class Fitter:

    _preloaded_exprs = {
        'gaus': '[norm]*exp(-0.5*(x-[mean])*(x-[mean])/([sigma]*[sigma]))'
    }

    def __init__(self, fit_expr) -> None:

        self._fit_expr = fit_expr
        preloaded_expr = self._check_preloaded_expr(fit_expr)
        if preloaded_expr:
            self._fit_expr = preloaded_expr

        self._param_dict = {}
        self._update_param_dict()
        self._data = None

    @property
    def data(self) -> TH1F:
        return self._data
    
    @property
    def fit_expr(self) -> str:
        return self._fit_expr
    
    @property
    def param_dict(self) -> dict:
        return self._param_dict

    def load_data(self, data) -> None:  
        if 'Series' in str(type(data)):
            raise NotImplementedError
        elif 'TH1' in str(type(data)):
            self._data = data

    def _check_preloaded_expr(self, expr) -> str | None:
        if expr in self._preloaded_exprs:
            return self._preloaded_exprs[expr]
        return None

    def _update_param_dict(self) -> None:
        param_list = re.findall(r'\[([^\]]+)\]', self._fit_expr)        
        self._param_dict = {param: 0. for param in param_list}

    def initialise_params(self, params) -> None:
        for param, value in params.items():
            self._param_dict[param] = value
    
    def set_param(self, param, value) -> None:
        if param not in self._param_dict:
            raise ValueError(f'Parameter {param} not found in the fit expression.')
        self._param_dict[param] = value

    @overload
    @signature('str', 'float', 'float')
    def set_param_range(self, param, min_val, max_val) -> None:
        if param not in self._param_dict:
            raise ValueError(f'Parameter {param} not found in the fit expression.')
        self._param_dict[param + '_min'] = min_val
        self._param_dict[param + '_max'] = max_val
    
    @set_param_range.overload
    @signature('str', 'float', 'float', 'float')
    def set_param_range(self, param, value, min_val, max_val) -> None:
        if param not in self._param_dict:
            raise ValueError(f'Parameter {param} not found in the fit expression.')
        self._param_dict[param] = value
        self._param_dict[param + '_min'] = min_val
        self._param_dict[param + '_max'] = max_val


    def auto_initialise(self, init_mode, **kwargs) -> None:
        if init_mode == 'gaus':
            self.auto_gaus_initialise()
        elif init_mode == 'multi_gaus':
            self.auto_multi_gaus_initialise(kwargs.get('n_components', 2))
        else:
            raise NotImplementedError
    
    def auto_gaus_initialise(self) -> None:
        '''
            Automatically initialise the parameters for a Gaussian fit.
            Requires the data to be loaded and that the parameters in the fit expression are named 'mean', 'sigma' and 'norm'.
        '''
        if not self._data:
            raise ValueError('No data loaded.')
        
        mean = self._data.GetMean()
        sigma = self._data.GetRMS()
        norm = self._data.GetBinContent(self._data.FindBin(mean))

        self._param_dict['mean'] = mean
        self._param_dict['sigma'] = sigma
        self._param_dict['norm'] = norm

    def auto_multi_gaus_initialise(self, n_components=2) -> None:
        '''
            Automatically initialise the parameters for a gaussian fit.
            Expects the function with lower mean to be the first in the list.
        '''

        data_points = []
        for ibin in range(1, self._data.GetNbinsX()+1):
            data_points.extend([self._data.GetBinCenter(ibin)] * int(self._data.GetBinContent(ibin)))

        data_points = np.array(data_points)
        if len(data_points) <= 0:
            print('No data points to fit')
            return

        kmeans = KMeans(n_clusters=n_components, init='k-means++', n_init='auto').fit(data_points.reshape(-1, 1))
        centers = kmeans.cluster_centers_
        labels = kmeans.labels_

        covariances = []
        for icomp in range(n_components):
            comp_data = data_points[np.where(np.array(labels)==icomp)[0]]
            covariances.append(np.cov(comp_data.T))

        weights = np.array([np.sum(labels==icomp) for icomp in range(n_components)])
        max_bin_contents = [np.max(data_points[labels == icomp]) for icomp in range(n_components)]

        # Sort centers and get the sorted indices
        sorted_indices = np.argsort(centers.flatten()) 

        # Reorder centers, covariances, and weights based on the sorted indices
        centers = centers[sorted_indices]
        covariances = [covariances[i] for i in sorted_indices]
        weights = weights[sorted_indices]
        max_bin_contents = [np.max(data_points[labels == icomp]) for icomp in range(n_components)]

        for icomp, label in enumerate(set(labels)):
            mean = centers[label][0]
            std = np.sqrt(covariances[label])
            peak_height = max_bin_contents[label]
            norm = peak_height * np.sqrt(2*np.pi) * std
            self._param_dict[f'mean{icomp}'] = mean
            self._param_dict[f'sigma{icomp}'] = std
            self._param_dict[f'norm{icomp}'] = norm

    def fit(self, **kwargs) -> dict:

        if not self._data:
            raise ValueError('No data loaded.')

        fit_func = TF1('fit_func', self._fit_expr)
        for param, value in self._param_dict.items():
            if param.endswith('_min') or param.endswith('_max'):
                continue
            fit_func.SetParameter(param, value)
            if str(f'{param}_min') in self._param_dict and str(f'{param}_max') in self._param_dict:
                fit_func.SetParLimits(fit_func.GetParNumber(param), self._param_dict[f'{param}_min'], self._param_dict[f'{param}_max'])
        if kwargs.get('fit_range'):
            fit_func.SetRange(*kwargs['fit_range'])
        fit_status = self._data.Fit(fit_func, kwargs.get('fit_options', ''))
        
        fit_result = {}
        for param, value in self._param_dict.items():
            if param.endswith('_min') or param.endswith('_max'):
                continue
            fit_result[param] = fit_func.GetParameter(param)
            fit_result[param + '_err'] = fit_func.GetParError(fit_func.GetParNumber(param))
        fit_result['chi2'] = fit_func.GetChisquare()
        fit_result['ndf'] = fit_func.GetNDF()
        fit_result['successful'] = True #fit_status.IsValid()
        
        return fit_result

# Standalone functions

def fit_TH1(h1: TH1F, fitter: Fitter, **kwargs) -> dict:
    '''
        Fit a TH1F histogram using a Fitter object.
        Returns a dictionary with the fit results.
    '''

    fitter.load_data(h1)
    if kwargs.get('init_mode', None):
        fitter.auto_initialise(kwargs['init_mode'], n_components=kwargs.get('n_components', 2))
    fit_result = fitter.fit(**kwargs)

    if fit_result['successful']:
        return fit_result
    return None

def fit_by_slices(h2: TH2F, fitter: Fitter, **kwargs) -> pd.DataFrame:
    '''
        Fit a TH2F histogram by slicing it along the x-axis and fitting each slice.
        Returns a pandas DataFrame with the fit results

        Parameters:
        - h2: TH2F
            The 2D histogram to be fitted by slices.
        - fitter: Fitter
            The Fitter object to be used for the fits.
        - **kwargs:
            Additional arguments to be passed to the fit_TH1 function.
            -> first_bin_fit_by_slices: int
                First bin to be fitted by slices.
            -> last_bin_fit_by_slices: int
                Last bin to be fitted by slices.
            -> init_mode: str
                Argument of Fitter.init_mode()
            -> fit_range: tuple
                Range for the slice fit
            -> fit_options: str
                Options for the slice fit
    '''

    fit_results = pd.DataFrame()
    for ibin in range(kwargs.get('first_bin_fit_by_slices', 1), kwargs.get('last_bin_fit_by_slices', h2.GetNbinsX() + 1)):
        h1 = h2.ProjectionY(f'proj_{ibin}', ibin, ibin)
        bin_fit_result = fit_TH1(h1, fitter, **kwargs)
        if bin_fit_result:
            if kwargs.get('output_dir', None):
                output_dir = kwargs['output_dir']
                output_dir.cd()
                h1.Write(f'h1_{ibin}')
            bin_fit_result['integral'] = h1.Integral(1, h1.GetNbinsX())
            bin_fit_result['bin_center'] = h2.GetXaxis().GetBinCenter(ibin)
            df_bin_fit_result = pd.DataFrame.from_dict({key: [val] for key, val in bin_fit_result.items()})
            fit_results = pd.concat([fit_results, df_bin_fit_result], ignore_index=True)
    return fit_results