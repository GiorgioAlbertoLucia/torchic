'''
    Various utility functions for creating histograms with ROOT
'''

from dataclasses import dataclass
from ROOT import TH1F, TH2F, TFile
from torchic.utils.overload import overload, signature

import numpy as np
import pandas as pd

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

def build_efficiency(hist_tot: TH1F, hist_sel: TH1F, name: str = None, xtitle: str = None, ytitle: str = "Efficiency") -> TH1F:
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