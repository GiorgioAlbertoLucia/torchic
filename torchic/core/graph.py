from functools import singledispatch
import pandas as pd
from ROOT import TGraphErrors, TFile

from torchic.core.histogram import HistLoadInfo

@singledispatch
def create_graph(arg, *args, **kwargs):
    raise NotImplementedError(f"Unsupported type: {type(arg)}")

@create_graph.register
def _(df: pd.DataFrame, x: str, y: str, ex, ey, name:str='', title:str='') -> TGraphErrors:
        '''
            Create a TGraphErrors from the input DataFrame

            Parameters
            ----------
            x (str): x-axis variable
            y (str): y-axis variable
            ex (str): x-axis error
            ey (str): y-axis error
        '''

        if len(df) == 0:
            return TGraphErrors()
        graph = TGraphErrors(len(df[x]))
        for irow, row in df.iterrows():
            graph.SetPoint(irow, row[x], row[y])
            xerr = row[ex] if ex != 0 else 0.
            yerr = row[ey] if ey != 0 else 0.
            graph.SetPointError(irow, xerr, yerr)
        
        graph.SetName(name)
        graph.SetTitle(title)

        return graph

@create_graph.register
def _(xs: list, ys: list, exs: list, eys: list, name:str='', title:str='') -> TGraphErrors:
        '''
            Create a TGraphErrors from the input DataFrame

            Parameters
            ----------
            x (str): x-axis variable
            y (str): y-axis variable
            ex (str): x-axis error
            ey (str): y-axis error
        '''

        if len(xs) == 0:
            return TGraphErrors()
        graph = TGraphErrors(len(xs))
        for idx in range(len(xs)):
            graph.SetPoint(idx, xs[idx], ys[idx])
            xerr = exs[idx] if exs[idx] != 0 else 0.
            yerr = eys[idx] if eys[idx] != 0 else 0.
            graph.SetPointError(idx, xerr, yerr)
        
        graph.SetName(name)
        graph.SetTitle(title)

        return graph

@singledispatch
def load_graph(arg, *args, **kwargs):
    raise NotImplementedError(f"Unsupported type: {type(arg)}")

@load_graph.register
def _(graph_load_info: HistLoadInfo):
    '''
        Load a graph from a ROOT file

        Args:
            graph_load_info (HistLoadInfo): The information needed to load the graph

        Returns:
            TH1F: The graph
    '''

    graph_file = TFile(graph_load_info.graph_file_path, 'READ')
    hist = graph_file.Get(graph_load_info.graph_name)
    graph_file.Close()
    return hist

@load_graph.register
def _(graph_file_path: str, graph_name: str):
    '''
        Load a graph from a ROOT file

        Args:
            graph_load_info (HistLoadInfo): The information needed to load the graph

        Returns:
            TH1F: The graph
    '''

    graph_file = TFile(graph_file_path, 'READ')
    hist = graph_file.Get(graph_name)
    graph_file.Close()
    return hist
