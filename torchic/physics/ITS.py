''''
    Functions and classes tailored to manage the data collected by ALICE ITS
'''

import numpy as np
import pandas as pd

N_ITS_LAYERS = 7

def unpack_cluster_sizes(cluster_sizes, layer) -> list:
    '''
        Unpack the cluster size from the data
    '''
    return (cluster_sizes >> layer*4) & 0b1111

def average_cluster_size(cluster_sizes: pd.Series) -> tuple:
    '''
        Compute the average cluster size
    '''
    
    np_cluster_sizes = cluster_sizes.to_numpy()
    avg_cluster_size = np.zeros(len(np_cluster_sizes))
    n_hits = np.zeros(len(np_cluster_sizes))
    for ilayer in range(N_ITS_LAYERS):
        avg_cluster_size += np.right_shift(np_cluster_sizes, 4*ilayer) & 0b1111    
        n_hits += np.right_shift(np_cluster_sizes, 4*ilayer) & 0b1111 > 0
    avg_cluster_size /= n_hits

    return avg_cluster_size, n_hits