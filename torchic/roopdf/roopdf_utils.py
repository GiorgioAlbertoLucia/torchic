from ROOT import RooRealVar, RooCrystalBall, RooGaussian

from torchic.roopdf import RooGausExp

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

def init_roopdf(pdf: str, var: RooRealVar, **kwargs):

    if pdf == 'crystal_ball':   return init_crystal_ball(var, **kwargs)
    elif pdf == 'gaus_exp':     return init_gaus_exp(var, **kwargs)
    elif pdf == 'gaus':         return init_gaus(var, **kwargs)
    else:   raise ValueError(f'Unknown function: {function}. Supported functions are "crystal_ball", "gaus_exp", "gaus".')
