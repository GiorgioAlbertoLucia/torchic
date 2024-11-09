import pandas as pd
from copy import deepcopy
from ROOT import TH1F, TH2F, TCanvas, TDirectory
from ROOT import RooRealVar, RooGaussian, RooCrystalBall, RooAddPdf, RooGenericPdf, RooArgList, RooDataHist, RooGExpModel, RooArgSet

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
        exp_mod_gaus = RooGExpModel(f'exp_mod_gaus_{self._pdf_counter}', f'exp_mod_gaus_{self._pdf_counter}',
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
        self._data_hist.plotOn(frame)
        self._model.plotOn(frame)
        self._model.paramOn(frame)
        for icomp, component in enumerate(funcs_to_plot):
            self._model.plotOn(frame, Components={self._pdfs[component]}, LineColor={DEFAULT_COLORS[icomp%N_COLORS]}, LineStyle='--')
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

# Standalone functions

def fit_by_slices_roofit(h2: TH2F, fitter: Roofitter, **kwargs) -> pd.DataFrame:
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
            -> output_dir: TDirectory
                Output directory to save the histograms
            -> xmin: float
                Minimum x-value to be used in the fit
            -> xmax: float
                Maximum x-value to be used in the fit
            -> funcs_to_fit: list
                List of functions to be fitted
    '''

    fit_results = pd.DataFrame()
    for ibin in range(kwargs.get('first_bin_fit_by_slices', 1), kwargs.get('last_bin_fit_by_slices', h2.GetNbinsX() + 1)):
        h1 = h2.ProjectionY(f'proj_{ibin}', ibin, ibin)
        fitter.fit(h1, **kwargs)

        if kwargs.get('output_dir', None):
            output_dir = kwargs['output_dir']
            fitter.plot(output_dir, canvas_name=f'c_{h1.GetName()}')
        bin_fit_results = deepcopy(fitter.fit_results)
        bin_fit_results['integral'] = h1.Integral(1, h1.GetNbinsX())
        bin_fit_results['bin_center'] = h2.GetXaxis().GetBinCenter(ibin)

        df_bin_fit_result = pd.DataFrame.from_dict({key: [value] for key, value in bin_fit_results.items()})
        fit_results = pd.concat([fit_results, df_bin_fit_result], ignore_index=True)
    return fit_results
