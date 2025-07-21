import numpy as np
from math import erf

DEFAULT_BETHEBLOCH_PARS = { # params for TPC He3 pp
                            'kp1': -241.490, 
                            'kp2': 0.374245,
                            'kp3': 1.397847,
                            'kp4': 1.078250,
                            'kp5': 2.048336
                          }

def py_BetheBloch(betagamma, kp1, kp2, kp3, kp4, kp5):
    '''
        Python implementation of the Bethe-Bloch formula.
    '''
    beta = betagamma / np.sqrt(1 + betagamma**2)
    aa = beta**kp4
    bb = (1/betagamma)**kp5
    bb = np.log(bb + kp3)
    return (kp2 - aa - bb) * kp1 / aa

def cluster_size_parametrisation(betagamma, kp1, kp2, kp3, charge, kp4):
    '''
        Python implementation of a simil Bethe-Bloch formula: kp1 / betagamma**kp2 + kp3
    '''
    return (kp1 / betagamma**kp2 + kp3) * charge ** kp4

def cluster_size_resolution(betagamma, rp0, rp1, rp2):
    '''
        Python implementation of the resolution function.
    '''
    return rp0 * erf((betagamma - rp1) / rp2)
np_cluster_size_resolution = np.vectorize(cluster_size_resolution)
