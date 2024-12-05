import os
from ROOT import TH1F, TCanvas, TDirectory, gInterpreter
from ROOT import RooRealVar, RooGaussian, RooCrystalBall, RooAddPdf, RooGenericPdf, RooArgList, RooDataHist, RooArgSet

from torchic.core.histogram import get_mean, get_rms

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOCUSTOMPDFS_DIR = os.path.join(CURRENT_DIR, 'RooCustomPdfs')
gInterpreter.ProcessLine(f'#include "{ROOCUSTOMPDFS_DIR}/RooGausExp.cxx"')
from ROOT import RooGausExp

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
        Available functions are:
            - 'gaus': Gaussian
            - 'exp_mod_gaus': Exponential modified Gaussian
            - 'exp': Exponential
            - 'exp_offset': Exponential with offset
            - 'comp_exp': Complementary exponential (i.e. 1 - exp(-alpha*x))
            - 'crystal_ball': Crystal Ball
            - 'polN': Polynomial of order N
    '''
    def __init__(self, x: RooRealVar, pdfs):
        self._x = x
        self._data_hist = None
        
        self._pdf_counter = 0 # Counter to keep track of the number of pdfs to assign them a unique name
        self._pdfs = {}
        self._pdf_params = {}
        self._fit_results = {}
        self._fit_fractions = {}

        for pdf in pdfs:
            self.build_pdf(pdf)

        self._model = None

    @property
    def pdf_params(self):
        return self._pdf_params

    @property
    def fit_results(self):
        return self._fit_results

    @property
    def fit_fractions(self):
        return self._fit_fractions

    @property
    def pdfs(self):
        return self._pdfs

    def init_param(self, name: str, value: float, min: float = None, max: float = None) -> None:
        '''
            Initialise the value of a RooRealVar parameter
        '''
        self._pdf_params[name].setVal(value)
        if min is not None and max is not None:
            self._pdf_params[name].setRange(min, max)   

    def build_pdf(self, pdf, args = None, return_function: bool = False, **kwargs):
        '''
            Add a pdf to the list of pdfs to be combined
        '''
        returned_function = None
        if pdf == 'gaus':
            returned_function = self._build_gaus(return_function=return_function)   
        elif pdf == 'exp_mod_gaus':
            returned_function = self._build_exp_mod_gaus(return_function=return_function)
        elif pdf == 'exp':
            returned_function = self._build_exp(return_function=return_function, exp_offset=kwargs.get('exp_offset', False))
        elif pdf == 'exp_offset':
            returned_function = self._build_exp(return_function=return_function, exp_offset=True)
        elif pdf == 'comp_exp':
            returned_function = self._build_comp_exp(return_function=return_function)
        elif pdf == 'crystal_ball':
            returned_function = self._build_crystal_ball(return_function=return_function)
        elif 'pol' in pdf:
            returned_function = self._build_polynomial(int(pdf.split('pol')[1]), return_function=return_function)
        else:
            raise ValueError(f'pdf {pdf} not recognized')
        
        if return_function:
            return returned_function

    def _build_gaus(self, x: RooRealVar = None, return_function: bool = False):

        if x is None:
            x = self._x

        self._pdf_params[f'gaus_{self._pdf_counter}_mean'] = RooRealVar(f'mean_{self._pdf_counter}', f'mean_{self._pdf_counter}', 0, -10, 10)
        self._pdf_params[f'gaus_{self._pdf_counter}_sigma'] = RooRealVar(f'sigma_{self._pdf_counter}', f'sigma_{self._pdf_counter}', 1, 0.001, 10)
        gaus = RooGaussian(f'gaus_{self._pdf_counter}', f'gaus_{self._pdf_counter}', x, self._pdf_params[f'gaus_{self._pdf_counter}_mean'], self._pdf_params[f'gaus_{self._pdf_counter}_sigma'])
        self._pdfs[f'gaus_{self._pdf_counter}'] = gaus
        self._pdf_counter += 1

        if return_function:
            return gaus, self._pdf_params[f'gaus_{self._pdf_counter}_mean'], self._pdf_params[f'gaus_{self._pdf_counter}_sigma']
        else:
            return None
    
    def _build_exp_mod_gaus(self, x: RooRealVar = None, return_function: bool = False) -> tuple | None:
        if x is None:
            x = self._x

        self._pdf_params[f'exp_mod_gaus_{self._pdf_counter}_mean'] = RooRealVar(f'mean_{self._pdf_counter}', f'mean_{self._pdf_counter}', 0, -10, 10)
        self._pdf_params[f'exp_mod_gaus_{self._pdf_counter}_sigma'] = RooRealVar(f'sigma_{self._pdf_counter}', f'sigma_{self._pdf_counter}', 1, 0.001, 10)
        self._pdf_params[f'exp_mod_gaus_{self._pdf_counter}_tau'] = RooRealVar(f'tau_{self._pdf_counter}', f'tau_{self._pdf_counter}', -0.5, -10, 0)
        exp_mod_gaus = RooGausExp(f'exp_mod_gaus_{self._pdf_counter}', f'exp_mod_gaus_{self._pdf_counter}',
                                    x, self._pdf_params[f'exp_mod_gaus_{self._pdf_counter}_mean'], 
                                    self._pdf_params[f'exp_mod_gaus_{self._pdf_counter}_sigma'], self._pdf_params[f'exp_mod_gaus_{self._pdf_counter}_tau'])
        self._pdfs[f'exp_mod_gaus_{self._pdf_counter}'] = exp_mod_gaus
        self._pdf_counter += 1

        if return_function:
            return exp_mod_gaus, self._pdf_params[f'exp_mod_gaus_{self._pdf_counter}_mean'], self._pdf_params[f'exp_mod_gaus_{self._pdf_counter}_sigma'], self._pdf_params[f'exp_mod_gaus_{self._pdf_counter}_tau']
        else:
            return None

    def _build_exp(self, x: RooRealVar = None, return_function: bool = False, exp_offset: bool = False) -> tuple | None:
        
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
        
    def _build_comp_exp(self, x: RooRealVar = None, return_function: bool = False) -> tuple | None:

        alpha = RooRealVar(f'alpha_{self._pdf_counter}', f'alpha_{self._pdf_counter}', -0.5, -10, 0)
        offset = None
        exp = RooGenericPdf(f'comp_exp_{self._pdf_counter}', f'comp_exp_{self._pdf_counter}', f'1 - exp(-alpha_{self._pdf_counter}*x)', RooArgList(self._x, alpha))
        self._pdf_params[f'comp_exp_{self._pdf_counter}_alpha'] = alpha
        self._pdfs[f'comp_exp_{self._pdf_counter}'] = exp
        #if exp_offset:
        #    offset = RooRealVar(f'offset_{self._pdf_counter}', f'offset_{self._pdf_counter}', 1, -100, 100)
        #    exp_offset = RooGenericPdf(f'exp_{self._pdf_counter}', f'exp_{self._pdf_counter}', f'1 - exp(-alpha_{self._pdf_counter}*(x + offset_{self._pdf_counter}))', RooArgList(self._x, alpha, offset))
        #    self._pdf_params[f'exp_{self._pdf_counter}_offset'] = offset
        #    self._pdfs[f'exp_{self._pdf_counter}'] = exp_offset
        self._pdf_counter += 1

        if return_function:
            return exp, alpha, offset
        else:
            return None
    
    def _build_crystal_ball(self, x: RooRealVar = None, return_function: bool = False) -> tuple | None:
        if x is None:
            x = self._x

        self._pdf_params[f'crystal_ball_{self._pdf_counter}_mean'] = RooRealVar(f'mean_{self._pdf_counter}', f'mean_{self._pdf_counter}', 0, -10, 10)
        self._pdf_params[f'crystal_ball_{self._pdf_counter}_sigma'] = RooRealVar(f'sigma_{self._pdf_counter}', f'sigma_{self._pdf_counter}', 1, 0.001, 10)
        self._pdf_params[f'crystal_ball_{self._pdf_counter}_alphaL'] = RooRealVar(f'alphaL_{self._pdf_counter}', f'alphaL_{self._pdf_counter}', 1, 0, 10)
        self._pdf_params[f'crystal_ball_{self._pdf_counter}_nL'] = RooRealVar(f'nL_{self._pdf_counter}', f'nL_{self._pdf_counter}', 1, 0, 10)
        self._pdf_params[f'crystal_ball_{self._pdf_counter}_alphaR'] = RooRealVar(f'alphaR_{self._pdf_counter}', f'alphaR_{self._pdf_counter}', 1, 0, 10)
        self._pdf_params[f'crystal_ball_{self._pdf_counter}_nR'] = RooRealVar(f'nR_{self._pdf_counter}', f'nR_{self._pdf_counter}', 1, 0, 10)

        crystal_ball = RooCrystalBall(f'crystal_ball_{self._pdf_counter}', f'crystal_ball_{self._pdf_counter}', x, 
                                      self._pdf_params[f'crystal_ball_{self._pdf_counter}_mean'], self._pdf_params[f'crystal_ball_{self._pdf_counter}_sigma'], 
                                      self._pdf_params[f'crystal_ball_{self._pdf_counter}_alphaL'], self._pdf_params[f'crystal_ball_{self._pdf_counter}_nL'], 
                                      self._pdf_params[f'crystal_ball_{self._pdf_counter}_alphaR'], self._pdf_params[f'crystal_ball_{self._pdf_counter}_nR'])
        self._pdfs[f'crystal_ball_{self._pdf_counter}'] = crystal_ball
        self._pdf_counter += 1

        if return_function:
            return crystal_ball, self._pdf_params[f'crystal_ball_{self._pdf_counter}_mean'], self._pdf_params[f'crystal_ball_{self._pdf_counter}_sigma'], self._pdf_params[f'crystal_ball_{self._pdf_counter}_alphaL'], self._pdf_params[f'crystal_ball_{self._pdf_counter}_nL'], self._pdf_params[f'crystal_ball_{self._pdf_counter}_alphaR'], self._pdf_params[f'crystal_ball_{self._pdf_counter}_nR']
        else:
            return None

    def _build_polynomial(self, order: int, x: RooRealVar = None, return_function: bool = False) -> tuple | None:
        if x is None:
            x = self._x

        for i in range(order+1):
            self._pdf_params[f'pol{order}_{self._pdf_counter}_coeff{i}'] = RooRealVar(f'coeff{i}_{self._pdf_counter}', f'coeff{i}_{self._pdf_counter}', 0, -10, 10)

        polynomial = RooGenericPdf(f'pol{order}_{self._pdf_counter}', f'pol{order}_{self._pdf_counter}', 
                                   '+'.join([f'coeff{i}_{self._pdf_counter}*pow(x, {i})' for i in range(order+1)]), 
                                   RooArgList(x, *[self._pdf_params[f'pol{order}_{self._pdf_counter}_coeff{i}'] for i in range(order+1)]))
        self._pdfs[f'pol{order}_{self._pdf_counter}'] = polynomial
        self._pdf_counter += 1

        if return_function:
            return polynomial, *[self._pdf_params[f'pol_{order}_{self._pdf_counter}_coeff{i}'] for i in range(order+1)]
        else:
            return None

    def init_gaus(self, hist: TH1F, func_name: str, xmin: float = None, xmax: float = None) -> None:
        '''
            Initialise the parameters of a Gaussian function from a histogram
        '''
        mean = get_mean(hist, xmin, xmax)
        sigma = get_rms(hist, xmin, xmax)
        self.init_param(f'{func_name}_mean', mean, mean - 3*sigma, mean + 3*sigma)
        self.init_param(f'{func_name}_sigma', sigma, 0.1, 3*sigma)

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

        fractions = [RooRealVar(f'fraction_{func}', f'fraction_{func}', 0.5, 0, 1) for func in funcs_to_fit[:-1]]
        self._model = RooAddPdf('model', kwargs.get('title', 'model'), RooArgList(*[self._pdfs[func] for func in funcs_to_fit]), RooArgList(*fractions))
        self._model.fixCoefNormalization(RooArgSet(self._x))
        
        self._data_hist = RooDataHist('data_hist', 'data_hist', RooArgList(self._x), hist)
        self._model.fitTo(self._data_hist, PrintLevel=kwargs.get('fit_print_level', -1))
        fractions.append(RooRealVar(f'fraction_{funcs_to_fit[-1]}', f'fraction_{funcs_to_fit[-1]}', 1 - sum([frac.getVal() for frac in fractions]), 0, 1))
        self._fit_fractions = {func_to_fit: fractions[ifrac].getVal() for ifrac, func_to_fit in enumerate(funcs_to_fit)}

        for parname, par in self._pdf_params.items():
            self._fit_results[parname] = par.getVal()
            self._fit_results[parname + '_err'] = par.getError()
        chi2 = self._model.createChi2(self._data_hist)
        self._fit_results['integral'] = self._model.createIntegral(RooArgList(self._x)).getVal()
        self._fit_results['chi2'] = chi2.getVal()
        self._fit_results['ndf'] = hist.GetNbinsX() - len(self._pdf_params)

        return fractions

    def plot(self, output_file: TDirectory, **kwargs) -> None:

        if 'funcs_to_plot' in kwargs:
            funcs_to_plot = kwargs['funcs_to_plot']
        else:
            funcs_to_plot = list(self._pdfs.keys())

        canvas = TCanvas(kwargs.get('canvas_name', 'canvas'), 'canvas', 800, 600)
        frame = self._x.frame()
        if not self._data_hist:
            print('No data histogram to plot')
        else: 
            print('Plotting data histogram')
        self._data_hist.plotOn(frame)
        if not self._model:
            print('No model to plot')
        else:
            print('Plotting model')
        self._model.plotOn(frame)
        self._model.paramOn(frame)
        for icomp, component in enumerate(funcs_to_plot):
            self._model.plotOn(frame, Components={self._pdfs[component]}, LineColor={DEFAULT_COLORS[icomp%N_COLORS]}, LineStyle={'--'})
        frame.GetXaxis().SetTitle(kwargs.get('xtitle', ''))
        frame.Draw('same')

        output_file.cd()
        canvas.Write()

    def functions_integral(self, xmin: float, xmax: float) -> float:
        '''
            Compute the integral of the functions in the model in a given range. 
            Useful for computing the signal and background integrals or purity.
        '''

        self._x.setRange('integral_range', xmin, xmax)
        integrals = {}

        pdf_list = self._model.pdfList()
        for pdf in pdf_list:
            norm_integral = pdf.createIntegral(self._x, self._x, 'integral_range').getVal()
            exp_events = self._model.expectedEvents(RooArgSet(self._x)) 
            integral = norm_integral * self._fit_fractions[pdf.GetName()]#* exp_events
            integrals[pdf.GetName()] = integral
        return integrals
