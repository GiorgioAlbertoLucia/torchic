'''
    Various utility functions for creating histograms with ROOT
'''

from dataclasses import dataclass
from ROOT import TH1F, TH2F, TFile
from torchic.utils.overload import overload, signature

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

def getEfficiency(histTot, histSel, name:str, xtitle="p_{T} (GeV/c)", ytitle="Efficiency"):
    '''
        Compute the efficiency of a selection

        Args:
            histTot, histSel (TH1F): The total and selected histograms
            name (str): The name of the efficiency plot
            xtitle (str): The x-axis title
            ytitle (str): The y-axis title

        Returns:
            TH1F: The efficiency histogram
    '''
    hEff = TH1F(name, f'{name}; f{xtitle} ; f{ytitle}', histTot.GetNbinsX(), histTot.GetXaxis().GetXmin(), histTot.GetXaxis().GetXmax())
    for xbin in range(1, histTot.GetNbinsX()+1):
            if histTot.GetBinContent(xbin) > 0:
                eff = histSel.GetBinContent(xbin)/histTot.GetBinContent(xbin)
                if(eff<=1):
                    effErr = np.sqrt(eff*(1-eff)/histTot.GetBinContent(xbin))
                    hEff.SetBinError(xbin, effErr)
                else:
                    hEff.SetBinError(xbin, 0)
                hEff.SetBinContent(xbin, eff)
    return hEff

def getNormalizeHist(hist):
    '''
        Return normalized histogram

        Args:
            hist (TH1F): The histogram to normalize

        Returns:
            TH1F: The efficiency histogram
    '''
    integral=hist.Integral()
    if integral>0:
        hist.Scale(1./integral)
    return hist