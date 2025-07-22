from  torchic.utils.overload import overload, signature
from ROOT import TGraphErrors

@overload
@signature('DataFrame', str, str, str, str, str, str)
def create_graph(df, x: str, y: str, ex, ey, name:str='', title:str='') -> TGraphErrors:
        '''
            Create a TGraphErrors from the input DataFrame

            Parameters
            ----------
            x (str): x-axis variable
            y (str): y-axis variable
            ex (str): x-axis error
            ey (str): y-axis error
        '''

        # eliminate None values on x, y
        #df = df.filter(df[x].is_not_null())
        #df = df.filter(df[y].is_not_null())

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

@overload
@signature('list', 'list', 'list', 'list', str, str)
def create_graph(xs, ys, exs, eys, name:str='', title:str='') -> TGraphErrors:
        '''
            Create a TGraphErrors from the input DataFrame

            Parameters
            ----------
            x (str): x-axis variable
            y (str): y-axis variable
            ex (str): x-axis error
            ey (str): y-axis error
        '''

        # eliminate None values on x, y
        #df = df.filter(df[x].is_not_null())
        #df = df.filter(df[y].is_not_null())

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