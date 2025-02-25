''''
    Functions and classes tailored to manage the data collected by ALICE ITS
'''

import numpy as np
import pandas as pd
from scipy.special import erf

N_ITS_LAYERS = 7
PID_ITS_PARAMETERS = {
    '''
        Parameters for the PID of ITS
        [kp1, kp2, kp3, res1, res2, res3]

        avg_cluster_size = kp1 / bg^kp2 + kp3
        sigma_cluster_size = avg_cluster_size * (res1 * erf((bg - res2) / res3))
    '''
    'Pr': [1.18941, 1.53792, 1.69961, 1.94669e-01, -2.08616e-01, 1.30753],
    'He': [2.35117, 1.80347, 5.14355, 8.74371e-02, -1.82804, 5.06449e-01],
}

def unpack_cluster_sizes(cluster_sizes, layer) -> list:
    '''
        Unpack the cluster size from the data
    '''
    return (cluster_sizes >> layer*4) & 0b1111

def average_cluster_size(cluster_sizes: pd.Series) -> tuple:
    '''
        Compute the average cluster size. A truncated mean will be used to avoid the presence of outliers.
    '''
    
    np_cluster_sizes = cluster_sizes.to_numpy()
    avg_cluster_size = np.zeros(len(np_cluster_sizes))
    max_cluster_size = 0
    n_hits = np.zeros(len(np_cluster_sizes))
    for ilayer in range(N_ITS_LAYERS):
        cluster_size_layer = np.right_shift(np_cluster_sizes, 4*ilayer) & 0b1111
        avg_cluster_size += cluster_size_layer
        n_hits += (cluster_size_layer > 0).astype(int)
        max_cluster_size = np.maximum(max_cluster_size, cluster_size_layer)
    
    avg_cluster_size = (avg_cluster_size - max_cluster_size) / (n_hits - 1)
    # avg_cluster_size /= n_hits

    return avg_cluster_size, n_hits

def expected_cluster_size(beta_gamma: pd.Series, pid_parameters: tuple = None, particle: str = None) -> pd.Series:
    '''
        Compute the expected cluster size for a given particle
    '''

    if pid_parameters is None and particle is None:
        raise ValueError('Either pid_parameters or particle must be provided')
    if pid_parameters is None and particle is not None:
        pid_parameters = PID_ITS_PARAMETERS[particle]
    
    kp1, kp2, kp3, res1, res2, res3 = pid_parameters
    avg_cluster_size = kp1 / beta_gamma**kp2 + kp3
    return avg_cluster_size

def sigma_its(beta_gamma: pd.Series, pid_parameters: tuple = None, particle: str = None) -> pd.Series:
    '''
        Compute the number of sigmas of the cluster size for a given particle
    '''

    if pid_parameters is None and particle is None:
        raise ValueError('Either pid_parameters or particle must be provided')
    if pid_parameters is None and particle is not None:
        pid_parameters = PID_ITS_PARAMETERS[particle]
    
    kp1, kp2, kp3, res1, res2, res3 = pid_parameters
    avg_cluster_size = kp1 / beta_gamma**kp2 + kp3
    sigma_cluster_size = avg_cluster_size * (res1 * erf((beta_gamma - res2) / res3))
    return sigma_cluster_size
