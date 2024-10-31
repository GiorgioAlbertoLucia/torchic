from ROOT import TH1F, TCanvas, TDirectory
from ROOT import RooRealVar, RooGaussian, RooAddPdf, RooGenericPdf, RooArgList, RooDataHist

DEFAULT_COLORS = [
    797,    # kOrange-3
    418,    # kGreen+2
    632,    # kRed+2
    430,    # kCyan-2
]
N_COLORS = len(DEFAULT_COLORS)

class Roofitter:
    '''
        Class to fit a RooFit model to data. Multiple functions can be combined.
    '''
    def __init__(self, x: RooRealVar, pdfs):
        self._x = x
        self._data_hist = None
        
        self._pdf_counter = 0 # Counter to keep track of the number of pdfs to assign them a unique name
        self._pdfs = {}
        self._pdf_params = {}

        for pdf in pdfs:
            self.build_pdf(pdf)

        self._model = None

    def build_pdf(self, pdf, args = None, return_function: bool = False, **kwargs):
        '''
            Add a pdf to the list of pdfs to be combined
        '''
        returned_function = None
        if pdf ==  'gaus':
            returned_function = self._build_gaus(return_function=return_function)   
        elif pdf == 'exp':
            returned_function = self._build_exp(return_function=return_function, exp_offset=kwargs.get('exp_offset', False))
        elif pdf == 'exp_offset':
            returned_function = self._build_exp(return_function=return_function, exp_offset=True)
        else:
            raise ValueError(f'pdf {pdf} not recognized')
        
        if return_function:
            return returned_function

    def _build_gaus(self, x: RooRealVar = None, return_function: bool = False):

        if x is None:
            x = self._x
        mean = RooRealVar(f'mean_{self._pdf_counter}', f'mean_{self._pdf_counter}', 0, -10, 10)
        sigma = RooRealVar(f'sigma_{self._pdf_counter}', f'sigma_{self._pdf_counter}', 1, 0.001, 10)
        self._pdf_params[f'gaus_{self._pdf_counter}_mean'] = mean
        self._pdf_params[f'gaus_{self._pdf_counter}_sigma'] = sigma
        gaus = RooGaussian(f'gaus_{self._pdf_counter}', f'gaus_{self._pdf_counter}', x, self._pdf_params[f'gaus_{self._pdf_counter}_mean'], self._pdf_params[f'gaus_{self._pdf_counter}_sigma'])
        self._pdfs[f'gaus_{self._pdf_counter}'] = gaus
        self._pdf_counter += 1

        if return_function:
            return gaus, mean, sigma
        else:
            return None
        
    def _build_exp(self, x: RooRealVar = None, return_function: bool = False, exp_offset: bool = False):
        
        alpha = RooRealVar(f'alpha_{self._pdf_counter}', f'alpha_{self._pdf_counter}', -0.5, -10, 0)
        offset = None
        exp = RooGenericPdf(f'exp_{self._pdf_counter}', f'exp_{self._pdf_counter}', f'exp(-alpha_{self._pdf_counter}*x)', RooArgList(self._x, alpha))
        self._pdf_params[f'exp_{self._pdf_counter}_alpha'] = alpha
        self._pdfs[f'exp_{self._pdf_counter}'] = exp
        if exp_offset:
            offset = RooRealVar(f'offset_{self._pdf_counter}', f'offset_{self._pdf_counter}', 1, -100, 100)
            exp_offset = RooGenericPdf(f'exp_{self._pdf_counter}', f'exp_{self._pdf_counter}', f'exp(-alpha_{self._pdf_counter}*(x + offset_{self._pdf_counter}))', RooArgList(self._x, alpha, offset))
            self._pdf_params[f'exp_{self._pdf_counter}_offset'] = offset
            self._pdfs[f'exp_{self._pdf_counter}'] = exp_offset
        self._pdf_counter += 1

        if return_function:
            return exp, alpha, offset
        else:
            return None
        
    def fit(self, hist: TH1F, xmin: float = None, xmax: float = None, **kwargs) -> list:
        '''
            Fit the pdf to the data
        '''
        if xmin is not None and xmax is not None:
            self._x.setRange('fit_range', xmin, xmax)
        
        if 'funcs_to_fit' in kwargs:
            funcs_to_fit = kwargs['funcs_to_fit']
        else:
            funcs_to_fit = list(self._pdfs.keys())

        fractions = [RooRealVar(f'fraction_{func}', f'fraction_{func}', 0.5, 0, 1) for func in funcs_to_fit]
        self._model = RooAddPdf('model', kwargs.get('title', 'model'), RooArgList(*[self._pdfs[func] for func in funcs_to_fit]), RooArgList(*fractions))
        
        self._data_hist = RooDataHist('data_hist', 'data_hist', RooArgList(self._x), hist)
        self._model.fitTo(self._data_hist, PrintLevel=kwargs.get('fit_print_level', -1))

        return fractions

    def plot(self, output_file: TDirectory, **kwargs) -> None:

        canvas = TCanvas(kwargs.get('canvas_name', 'canvas'), 'canvas', 800, 600)
        frame = self._x.frame()
        self._data_hist.plotOn(frame)
        self._model.plotOn(frame)
        self._model.paramOn(frame)
        for icomp, component in enumerate(self._pdfs.values()):
            self._model.plotOn(frame, Components={component}, LineColor={DEFAULT_COLORS[icomp%N_COLORS]}, LineStyle='--')
        frame.GetXaxis().SetTitle(kwargs.get('xtitle', ''))
        frame.Draw('same')

        output_file.cd()
        canvas.Write()
