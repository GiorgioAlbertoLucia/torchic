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
    def __init__(self, x: RooRealVar, pdfs:list = []) -> None:
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
            Initialise the value of a RooRealVar parameter.
            Default names of the parameters are in the format 'pdfname_counter_paramname'
            Parameters associated with functions are:
            * GAUS: mean, sigma
            * EXP_MOD_GAUS: mean, sigma, tau
            * EXP: alpha, offset (if exp_offset is True)
            * COMP_EXP (1 - exp): alpha 
            * CRYSTAL_BALL: mean, sigma, alphaL, nL, alphaR, nR (if double_sided is True)
            * CRYSTAL_BALL: mean, sigma, alpha, n (if double_sided is False)
            * POLN: coeff{i} for i in range(N+1)
        '''
        self._pdf_params[name].setVal(value)
        if min is not None and max is not None:
            self._pdf_params[name].setRange(min, max)   

    def build_pdf(self, pdf, args = None, **kwargs):
        '''
            Add a pdf to the list of pdfs to be combined

            Args:
                pdf (str): Name of the pdf to build
                args (list): List of arguments to pass to the pdf function (NOT IMPLEMENTED)
                kwargs (dict): Dictionary of keyword arguments to pass to the pdf function
                Accepted kwargs:
                * exp_offset (bool): If True, the exponential function will have an offset (default: False)
                * double_sided (bool): If True, the crystal ball function will be double sided (default: True)
        '''
        returned_function = None
        if pdf == 'gaus':
            returned_function = self._build_gaus()   
        elif pdf == 'exp_mod_gaus':
            returned_function = self._build_exp_mod_gaus()
        elif pdf == 'exp':
            returned_function = self._build_exp(exp_offset=kwargs.get('exp_offset', False))
        elif pdf == 'exp_offset':
            returned_function = self._build_exp(exp_offset=True)
        elif pdf == 'comp_exp':
            returned_function = self._build_comp_exp()
        elif pdf == 'crystal_ball':
            returned_function = self._build_crystal_ball(double_sided=kwargs.get('double_sided', True))
        elif 'pol' in pdf:
            returned_function = self._build_polynomial(int(pdf.split('pol')[1]))
        else:
            raise ValueError(f'pdf {pdf} not recognized')
        
        return returned_function

    def _build_gaus(self, x: RooRealVar = None):

        if x is None:
            x = self._x
        pdf_counter = self._pdf_counter
        self._pdf_counter += 1

        self._pdf_params[f'gaus_{pdf_counter}_mean'] = RooRealVar(f'mean_{pdf_counter}', f'mean_{pdf_counter}', 0, -10, 10)
        self._pdf_params[f'gaus_{pdf_counter}_sigma'] = RooRealVar(f'sigma_{pdf_counter}', f'sigma_{pdf_counter}', 1, 0.001, 10)
        gaus = RooGaussian(f'gaus_{pdf_counter}', f'gaus_{pdf_counter}', x, self._pdf_params[f'gaus_{pdf_counter}_mean'], self._pdf_params[f'gaus_{pdf_counter}_sigma'])
        self._pdfs[f'gaus_{pdf_counter}'] = gaus

        return gaus, self._pdf_params[f'gaus_{pdf_counter}_mean'], self._pdf_params[f'gaus_{pdf_counter}_sigma']
        
    def _build_exp_mod_gaus(self, x: RooRealVar = None) -> tuple | None:
        if x is None:
            x = self._x
        pdf_counter = self._pdf_counter
        self._pdf_counter += 1

        self._pdf_params[f'exp_mod_gaus_{pdf_counter}_mean'] = RooRealVar(f'mean_{pdf_counter}', f'mean_{pdf_counter}', 0, -10, 10)
        self._pdf_params[f'exp_mod_gaus_{pdf_counter}_sigma'] = RooRealVar(f'sigma_{pdf_counter}', f'sigma_{pdf_counter}', 1, 0.001, 10)
        self._pdf_params[f'exp_mod_gaus_{pdf_counter}_tau'] = RooRealVar(f'tau_{pdf_counter}', f'tau_{pdf_counter}', -0.5, -10, 0)
        exp_mod_gaus = RooGausExp(f'exp_mod_gaus_{pdf_counter}', f'exp_mod_gaus_{pdf_counter}',
                                    x, self._pdf_params[f'exp_mod_gaus_{pdf_counter}_mean'], 
                                    self._pdf_params[f'exp_mod_gaus_{pdf_counter}_sigma'], self._pdf_params[f'exp_mod_gaus_{pdf_counter}_tau'])
        self._pdfs[f'exp_mod_gaus_{pdf_counter}'] = exp_mod_gaus

        return exp_mod_gaus, self._pdf_params[f'exp_mod_gaus_{pdf_counter}_mean'], self._pdf_params[f'exp_mod_gaus_{pdf_counter}_sigma'], self._pdf_params[f'exp_mod_gaus_{pdf_counter}_tau']
        
    def _build_exp(self, x: RooRealVar = None, exp_offset: bool = False) -> tuple | None:
        
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

        return exp, alpha, offset
        
    def _build_comp_exp(self, x: RooRealVar = None) -> tuple | None:

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

        return exp, alpha, offset
        
    def _build_crystal_ball(self, x: RooRealVar = None, double_sided:bool=True) -> tuple | None:
        if x is None:
            x = self._x
            pdf_counter = self._pdf_counter
            self._pdf_counter += 1

        if double_sided:
            self._pdf_params[f'crystal_ball_{pdf_counter}_mean'] = RooRealVar(f'mean_{pdf_counter}', f'mean_{pdf_counter}', 0, -10, 10)
            self._pdf_params[f'crystal_ball_{pdf_counter}_sigma'] = RooRealVar(f'sigma_{pdf_counter}', f'sigma_{pdf_counter}', 1, 0.001, 10)
            self._pdf_params[f'crystal_ball_{pdf_counter}_alphaL'] = RooRealVar(f'alphaL_{pdf_counter}', f'alphaL_{pdf_counter}', 1, 0, 10)
            self._pdf_params[f'crystal_ball_{pdf_counter}_nL'] = RooRealVar(f'nL_{pdf_counter}', f'nL_{pdf_counter}', 1, 0, 10)
            self._pdf_params[f'crystal_ball_{pdf_counter}_alphaR'] = RooRealVar(f'alphaR_{pdf_counter}', f'alphaR_{pdf_counter}', 1, 0, 10)
            self._pdf_params[f'crystal_ball_{pdf_counter}_nR'] = RooRealVar(f'nR_{pdf_counter}', f'nR_{pdf_counter}', 1, 0, 10)

            crystal_ball = RooCrystalBall(f'crystal_ball_{pdf_counter}', f'crystal_ball_{pdf_counter}', x, 
                                          self._pdf_params[f'crystal_ball_{pdf_counter}_mean'], self._pdf_params[f'crystal_ball_{pdf_counter}_sigma'], 
                                          self._pdf_params[f'crystal_ball_{pdf_counter}_alphaL'], self._pdf_params[f'crystal_ball_{pdf_counter}_nL'], 
                                          self._pdf_params[f'crystal_ball_{pdf_counter}_alphaR'], self._pdf_params[f'crystal_ball_{pdf_counter}_nR'])
            self._pdfs[f'crystal_ball_{pdf_counter}'] = crystal_ball
            return crystal_ball, self._pdf_params[f'crystal_ball_{pdf_counter}_mean'], self._pdf_params[f'crystal_ball_{pdf_counter}_sigma'], self._pdf_params[f'crystal_ball_{pdf_counter}_alphaL'], self._pdf_params[f'crystal_ball_{pdf_counter}_nL'], self._pdf_params[f'crystal_ball_{pdf_counter}_alphaR'], self._pdf_params[f'crystal_ball_{pdf_counter}_nR']
        
        else: 
            self._pdf_params[f'crystal_ball_{pdf_counter}_mean'] = RooRealVar(f'mean_{pdf_counter}', f'mean_{pdf_counter}', 0, -10, 10)
            self._pdf_params[f'crystal_ball_{pdf_counter}_sigma'] = RooRealVar(f'sigma_{pdf_counter}', f'sigma_{pdf_counter}', 1, 0.001, 10)
            self._pdf_params[f'crystal_ball_{pdf_counter}_alpha'] = RooRealVar(f'alpha_{pdf_counter}', f'alpha_{pdf_counter}', 1, 0, 10)
            self._pdf_params[f'crystal_ball_{pdf_counter}_n'] = RooRealVar(f'n_{pdf_counter}', f'n_{pdf_counter}', 1, 0, 10)
            
            crystal_ball = RooCrystalBall(f'crystal_ball_{pdf_counter}', f'crystal_ball_{pdf_counter}', x, 
                                          self._pdf_params[f'crystal_ball_{pdf_counter}_mean'], self._pdf_params[f'crystal_ball_{pdf_counter}_sigma'], 
                                          self._pdf_params[f'crystal_ball_{pdf_counter}_alpha'], self._pdf_params[f'crystal_ball_{pdf_counter}_n'], 
                                          doubleSided=double_sided)
            self._pdfs[f'crystal_ball_{pdf_counter}'] = crystal_ball
            return crystal_ball, self._pdf_params[f'crystal_ball_{pdf_counter}_mean'], self._pdf_params[f'crystal_ball_{pdf_counter}_sigma'], self._pdf_params[f'crystal_ball_{pdf_counter}_alpha'], self._pdf_params[f'crystal_ball_{pdf_counter}_n']
        
    def _build_polynomial(self, order: int, x: RooRealVar = None) -> tuple | None:
        if x is None:
            x = self._x
        pdf_counter = self._pdf_counter
        self._pdf_counter += 1

        for i in range(order+1):
            self._pdf_params[f'pol{order}_{pdf_counter}_coeff{i}'] = RooRealVar(f'coeff{i}_{pdf_counter}', f'coeff{i}_{pdf_counter}', 0, -10, 10)

        polynomial = RooGenericPdf(f'pol{order}_{pdf_counter}', f'pol{order}_{pdf_counter}', 
                                   '+'.join([f'coeff{i}_{pdf_counter}*pow(x, {i})' for i in range(order+1)]), 
                                   RooArgList(x, *[self._pdf_params[f'pol{order}_{pdf_counter}_coeff{i}'] for i in range(order+1)]))
        self._pdfs[f'pol{order}_{pdf_counter}'] = polynomial

        return polynomial, *[self._pdf_params[f'pol{order}_{pdf_counter}_coeff{i}'] for i in range(order+1)]
        
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
