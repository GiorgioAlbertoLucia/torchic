from ROOT import RooRealVar, RooCrystalBall, RooGaussian, RooAddPdf, RooChebychev, RooGenericPdf

from torchic.roopdf import load_fit_modules
load_fit_modules()

def init_crystal_ball(var: RooRealVar, **kwargs):

    name = kwargs.get('name', 'pdf')
    pdf_pars = {
        'mean': kwargs.get('mean', RooRealVar(f'{name}_mean', 'mean', 0, '')),
        'sigma': kwargs.get('sigma', RooRealVar(f'{name}_sigma', 'sigma', 1., '')),
        'aL': kwargs.get('aL', RooRealVar(f'{name}_aL', 'aL', 0.7, 30.)),
        'nL': kwargs.get('nL', RooRealVar(f'{name}_nL', 'nL', 0.3, 30.)),
        'aR': kwargs.get('aR', RooRealVar(f'{name}_aR', 'aR', 0.7, 30.)),
        'nR': kwargs.get('nR', RooRealVar(f'{name}_nR', 'nR', 0.3, 30.)),
    }
    pdf = RooCrystalBall(name, name, var, *pdf_pars.values())
    if kwargs.get('double_sided', False):
        pdf = RooCrystalBall(name, name, var, *pdf_pars[['mean', 'sigma', 'aL', 'nL']].values(),
                             doubleSided=True)
    return pdf, pdf_pars

def init_gaus_exp(var: RooRealVar, **kwargs):

    from torchic.roopdf import RooGausExp
    name = kwargs.get('name', 'pdf')
    pdf_pars = {
        'mean': kwargs.get('mean', RooRealVar(f'{name}_mean', 'mean', 0, '')),
        'sigma': kwargs.get('sigma', RooRealVar(f'{name}_sigma', 'sigma', 1., '')),
        'rlife': kwargs.get('rlife', RooRealVar(f'{name}_rlife', 'rlife', 2., 0., 10.)),
    }
    pdf = RooGausExp(name, name, var, *pdf_pars.values())
    return pdf, pdf_pars

def init_gaus(var: RooRealVar, **kwargs):

    name = kwargs.get('name', 'pdf')
    pdf_pars = {
        'mean': kwargs.get('mean', RooRealVar(f'{name}_mean', 'mean', 0, '')),
        'sigma': kwargs.get('sigma', RooRealVar(f'{name}_sigma', 'sigma', 1., '')),
    }
    pdf = RooGaussian(name, name, var, *pdf_pars.values())
    return pdf, pdf_pars

def init_double_gaus(var: RooRealVar, **kwargs):
    
    name = kwargs.get('name', 'pdf')
    pdf_pars = {
        'mean': kwargs.get('mean', RooRealVar(f'{name}_mean', 'mean', 0, '')),
        'sigma1': kwargs.get('sigma1', RooRealVar(f'{name}_sigma1', 'sigma1', 1., '')),
        'sigma2': kwargs.get('sigma2', RooRealVar(f'{name}_sigma2', 'sigma2', 1., '')),
        'fraction': kwargs.get('fraction', RooRealVar(f'{name}_fraction', 'fraction', 0.5, 0., 1.)),
    }
    gaus1 = RooGaussian(f'{name}_gaus1', f'{name}_gaus1', var, pdf_pars['mean'], pdf_pars['sigma1'])
    gaus2 = RooGaussian(f'{name}_gaus2', f'{name}_gaus2', var, pdf_pars['mean'], pdf_pars['sigma2'])
    pdf = RooAddPdf(name, name, [gaus1, gaus2], [pdf_pars['fraction']])
    pdf_pars['gaus1'] = gaus1
    pdf_pars['gaus2'] = gaus2
    return pdf, pdf_pars

def init_roochebychev_n(var: RooRealVar, n: int, **kwargs):
    name = kwargs.get('name', 'pdf')
    pdf_pars = {}
    for i in range(n+1):
        pdf_pars[f'p{i}'] = kwargs.get(f'p{i}', RooRealVar(f'{name}_p{i}', f'p{i}', 0., -5., 5.))
    pdf = RooChebychev(name, name, var, list(pdf_pars.values()))
    return pdf, pdf_pars

def init_polynomial_n(var: RooRealVar, n: int, **kwargs):
    name = kwargs.get('name', 'pdf')
    if n < 0:
        raise ValueError(f"Polynomial degree n must be non-negative, got {n}.")
    pdf_pars = {'p0': kwargs.get('p0', RooRealVar(f'{name}_p0', 'p0', 1., -1.e3, 1.e3))}
    pdf_string = f'{name}_p0'
    for i in range(1, n+1):
        pdf_pars[f'p{i}'] = kwargs.get(f'p{i}', RooRealVar(f'{name}_p{i}', f'p{i}', 1., -1.e3, 1.e3))
        pdf_string += f' + {name}_p{i} * TMath::Power({var.GetName()}, {i})'
    print(f"Polynomial PDF string: {pdf_string}")
    pdf = RooGenericPdf(name, name, pdf_string, [var, *pdf_pars.values()])
    return pdf, pdf_pars
    
def init_roopdf(pdf: str, var: RooRealVar, **kwargs):

    if pdf == 'crystal_ball':       return init_crystal_ball(var, **kwargs)
    elif pdf == 'gaus_exp':         return init_gaus_exp(var, **kwargs)
    elif pdf == 'gaus':             return init_gaus(var, **kwargs)
    elif pdf == 'double_gaus':      return init_double_gaus(var, **kwargs)
    elif pdf == 'chebychev':        return init_roochebychev_n(var, **kwargs)
    elif pdf == 'polynomial':       return init_polynomial_n(var, **kwargs)
    else:   raise ValueError(f'Unknown function: {pdf}. Supported functions are "crystal_ball", "gaus_exp", "gaus", "double_gaus".')
