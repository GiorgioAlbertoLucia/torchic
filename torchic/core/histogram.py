'''
    Various utility functions for creating histograms with ROOT
'''

from dataclasses import dataclass
from ROOT import TH1F, TH2F, TFile
import boost_histogram as bh
from torchic.utils.overload import overload, signature

import numpy as np

@dataclass
class AxisSpec:

    nbins: int
    xmin: float
    xmax: float
    name: str = ''
    title: str = ''

    @classmethod
    def from_dict(cls, d: dict):
        return cls(d['nbins'], d['xmin'], d['xmax'], d['name'], d['title'])
    
@dataclass
class HistLoadInfo:
    hist_file_path: str
    hist_name: str

def build_TH1(data, axis_spec_x: AxisSpec) -> TH1F:
    '''
        Build a histogram with one axis

        Args:
            data (pd.Series): The data to be histogrammed
            axis_spec_x (AxisSpec): The specification for the x-axis

        Returns:
            TH1F: The histogram
    '''

    hist = TH1F(axis_spec_x.name, axis_spec_x.title, axis_spec_x.nbins, axis_spec_x.xmin, axis_spec_x.xmax)
    for x in data:
        hist.Fill(x)
    return hist

def build_TH2(data_x, data_y, axis_spec_x: AxisSpec, axis_spec_y: AxisSpec) -> TH2F:
    '''
        Build a histogram with two axes

        Args:
            data_x (pd.Series): The data to be histogrammed on the x-axis
            data_y (pd.Series): The data to be histogrammed on the y-axis
            axis_spec_x (AxisSpec): The specification for the x-axis
            axis_spec_y (AxisSpec): The specification for the y-axis

        Returns:
            TH1F: The histogram
    '''

    hist = TH2F(axis_spec_y.name, axis_spec_y.title, axis_spec_x.nbins, axis_spec_x.xmin, axis_spec_x.xmax, axis_spec_y.nbins, axis_spec_y.xmin, axis_spec_y.xmax)
    for x, y in zip(data_x, data_y):
        hist.Fill(x, y)
    return hist

def fill_TH1(data, hist: TH1F):
    '''
        Fill a histogram with data

        Args:
            data (pd.Series): The data to fill the histogram with
            hist (TH1F): The histogram to fill
    '''
    for x in data:
        hist.Fill(x)
    
def fill_TH2(data_x, data_y, hist: TH2F):
    '''
        Fill a 2D histogram with data

        Args:
            data_x (pd.Series): The data to fill the x-axis of the histogram with
            data_y (pd.Series): The data to fill the y-axis of the histogram with
            hist (TH2F): The histogram to fill
    '''
    for x, y in zip(data_x, data_y):
        hist.Fill(x, y)

def build_boost1(data, axis_spec_x: AxisSpec) -> bh.Histogram:
    '''
        Build a histogram with one axis

        Args:
            data (pd.Series): The data to be histogrammed
            axis_spec_x (AxisSpec): The specification for the x-axis

        Returns:
            TH1F: The histogram
    '''

    hist = bh.Histogram(
        bh.axis.Regular(axis_spec_x.nbins, axis_spec_x.xmin, axis_spec_x.xmax,
                        metadata=axis_spec_x.title)
    )
    hist.fill(data)
    return hist

def build_boost2(data_x, data_y, axis_spec_x: AxisSpec, axis_spec_y: AxisSpec) -> bh.Histogram:

    hist = bh.Histogram(
        bh.axis.Regular(axis_spec_x.nbins, axis_spec_x.xmin, axis_spec_x.xmax,
                        metadata=axis_spec_x.title),
        bh.axis.Regular(axis_spec_y.nbins, axis_spec_y.xmin, axis_spec_y.xmax,
                        metadata=axis_spec_y.title)
    )
    hist.fill(data_x, data_y)
    return hist

@overload
@signature('HistLoadInfo')
def load_hist(hist_load_info: HistLoadInfo):
    '''
        Load a histogram from a ROOT file

        Args:
            hist_load_info (HistLoadInfo): The information needed to load the histogram

        Returns:
            TH1F: The histogram
    '''

    hist_file = TFile(hist_load_info.hist_file_path, 'READ')
    hist = hist_file.Get(hist_load_info.hist_name)
    hist.SetDirectory(0)
    hist_file.Close()
    return hist

@overload
@signature(str, str)
def load_hist(hist_file_path: str, hist_name: str) -> TH1F:
    '''
        Load a histogram from a ROOT file

        Args:
            hist_load_info (HistLoadInfo): The information needed to load the histogram

        Returns:
            TH1F: The histogram
    '''

    hist_file = TFile(hist_file_path, 'READ')
    hist = hist_file.Get(hist_name)
    hist.SetDirectory(0)
    hist_file.Close()
    return hist

def build_efficiency(hist_tot: TH1F, hist_sel: TH1F, name: str = None, xtitle: str = None, ytitle: str = 'Efficiency') -> TH1F:
    '''
        Compute the efficiency of a selection

        Args:
            hist_tot, hist_sel (TH1F): The total and selected histograms (denominator, numerator)
            name (str): The name of the efficiency plot
            xtitle (str): The x-axis title
            ytitle (str): The y-axis title

        Returns:
            TH1F: The efficiency histogram
    '''
    if name is None:
        name = hist_sel.GetName() + "_eff"
    if xtitle is None:
        xtitle = hist_sel.GetXaxis().GetTitle()
    hist_eff = TH1F(name, f'{name}; f{xtitle} ; f{ytitle}', hist_tot.GetNbinsX(), hist_tot.GetXaxis().GetXmin(), hist_tot.GetXaxis().GetXmax())
    for xbin in range(1, hist_tot.GetNbinsX()+1):
            if hist_tot.GetBinContent(xbin) > 0:
                eff = hist_sel.GetBinContent(xbin)/hist_tot.GetBinContent(xbin)
                if eff <= 1:
                    eff_err = np.sqrt(eff * (1 - eff) / hist_tot.GetBinContent(xbin))
                    hist_eff.SetBinError(xbin, eff_err)
                else:
                    hist_eff.SetBinError(xbin, 0)
                hist_eff.SetBinContent(xbin, eff)
    return hist_eff

def normalize_hist(hist: TH1F, low_edge: float = None, high_edge: float = None, option: str = '') -> None:
    '''
        Return normalized histogram

        Args:
            hist (TH1F): The histogram to normalize

        Returns:
            TH1F: The efficiency histogram
    '''
    if low_edge is None or high_edge is None:
        low_edge = hist.GetXaxis().GetXmin()
        high_edge = hist.GetXaxis().GetXmax()
    integral = hist.Integral(hist.FindBin(low_edge), hist.FindBin(high_edge), option)
    if integral > 0:
        hist.Scale(1./integral, option)
        
def project_hist(hist2D: TH2F, min_int: float, max_int: float, name: str = None, var_to_project: str = 'Y') -> TH1F:
    '''
        Return the y or x projection of a 2D histogram
        
        Args:
            hist2D (TH2F): The 2D histogram to project
            name (str): The name of the projection
            min_int (float): The minimum value of the projection
            max_int (float): The maximum value of the projection
            var_to_project (str): The variable to project
            
        Returns:
            TH1F: The projected histogram
    '''
    if name is None:
        name = hist2D.GetName() + f'_proj_{var_to_project}'
    if var_to_project == 'Y':
        hist = hist2D.ProjectionY(name, hist2D.GetXaxis().FindBin(min_int), hist2D.GetXaxis().FindBin(max_int), 'e')
    elif var_to_project == 'X':
        hist = hist2D.ProjectionX(name, hist2D.GetYaxis().FindBin(min_int), hist2D.GetYaxis().FindBin(max_int), 'e')
    else:
        raise ValueError('var_to_project must be either X or Y')
    return hist
    
def scale_hist_axis(old_hist: TH1F, scale_factor: float, **kwargs) -> TH1F | None:
    '''
        Return a histogram with scaled axis
        
        Args:
            old_hist (TH1F): The histogram to scale
            scale_factor (float): The factor by which the axis is divided
        Kwargs:
            xmin, xmax (float): The low and high edges in which the histogram is scaled
            inplace (bool): If True, the scaling is done in place
            name (str): The name of the new histogram
            title (str): The title of the new histogram
            nbins (int): The number of bins of the new histogram
            xtitle (str): The x-axis title
            ytitle (str): The y-axis title
            
        Returns:
            TH1F: The scaled histogram
    '''
    
    nbins = kwargs.get('nbins', old_hist.GetNbinsX())
    xmin = kwargs.get('xmin', old_hist.GetXaxis().GetXmin())
    xmax = kwargs.get('xmax', old_hist.GetXaxis().GetXmax())
    new_hist = TH1F(kwargs.get('name', f'{old_hist.GetName()}_scaled'), kwargs.get('title', old_hist.GetTitle()), nbins, xmin, xmax)
    if 'xtitle' in kwargs:
        new_hist.GetXaxis().SetTitle(kwargs.get('xtitle'))
    new_hist.GetYaxis().SetTitle(kwargs.get('ytitle', 'Counts (a.u.)'))

    for ibin in range(old_hist.FindBin(xmin),old_hist.FindBin(xmax)):
        bin_center = (old_hist.GetXaxis().GetBinCenter(ibin)) / scale_factor
        bin_content = old_hist.GetBinContent(ibin)
        bin_index=new_hist.FindBin(bin_center)
        new_hist.SetBinContent(bin_index, bin_content)
        new_hist.SetBinError(bin_index, old_hist.GetBinError(ibin))
    
    if kwargs.get('inplace', False):
        old_hist.Clone(new_hist)
        del new_hist
        old_hist.SetName(kwargs.get('name', f'{old_hist.GetName()}_scaled'))
    else:
        return new_hist
    
def get_mean(hist: TH1F, low_edge: float = None, high_edge: float = None) -> float:
    '''
        Return the mean of a histogram in a given range
        
        Args:
            hist (TH1F): The histogram to calculate the mean
            low_edge (float): The lower edge of the histogram
            high_edge (float): The upper edge of the histogram
            
        Returns:
            float: The mean of the histogram
    '''
    if low_edge is None:
        low_edge = hist.GetXaxis().GetXmin()
    if high_edge is None:
        high_edge = hist.GetXaxis().GetXmax()
    total_content = 0
    total_value = 0
    for ibin in range(hist.FindBin(low_edge), hist.FindBin(high_edge)):
        total_content += hist.GetBinContent(ibin)
        total_value += hist.GetBinContent(ibin) * hist.GetBinCenter(ibin)
    if total_content > 0:
        return total_value / total_content
    else:
        return 0
    
def get_rms(hist: TH1F, low_edge: float = None, high_edge: float = None) -> float:
    '''
        Return the RMS of a histogram in a given range
        
        Args:
            hist (TH1F): The histogram to calculate the RMS
            low_edge (float): The lower edge of the histogram
            high_edge (float): The upper edge of the histogram
            
        Returns:
            float: The RMS of the histogram
    '''
    if low_edge is None:
        low_edge = hist.GetXaxis().GetXmin()
    if high_edge is None:
        high_edge = hist.GetXaxis().GetXmax()
    total_content = 0
    total_value = 0
    total_value2 = 0
    for ibin in range(hist.FindBin(low_edge), hist.FindBin(high_edge)):
        total_content += hist.GetBinContent(ibin)
        total_value += hist.GetBinContent(ibin) * hist.GetBinCenter(ibin)
        total_value2 += hist.GetBinContent(ibin) * hist.GetBinCenter(ibin) * hist.GetBinCenter(ibin)
    if total_content > 0:
        return np.sqrt(total_value2 / total_content - (total_value / total_content) ** 2)
    else:
        return 0
    