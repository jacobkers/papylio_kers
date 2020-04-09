# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 00:41:20 2020

@author: iason
"""

if __name__ == '__main__':
    import sys
    from pathlib import Path
    p = Path(__file__).parents[2]  # include two levels up directory where the trace_analysis package lies
    sys.path.insert(0, str(p))

import numpy as np
import pandas as pd
from trace_analysis import Experiment

import matplotlib.pyplot as plt
import seaborn as sns

def plot_fret_histogram(molecules, Emin=0.05, Emax=0.8, Irmin=40, Iroff=0,
                        Igoff=0, **hist_kwargs):
    plt.style.use('seaborn-dark')
    plt.style.use('seaborn-colorblind')
    data = np.concatenate([molecule.E(Irmin, Iroff, Igoff) for molecule in molecules])
    data = data[data>Emin]
    data = data[data<Emax]
    centers, bins, _ = plt.hist(data, **hist_kwargs)



def get_fret_per_offtime(molecule, trace_type='red', Emin=0.05, Emax=0.8,
                        Irmin=40, Iroff=0, Igoff=0):
    '''
    Function that gives the FRET trace with FRET calculcated only for the
    selected offtime interval. The FRET value outside this interval is put to
    zero

    '''

    if mol.steps is None:
        return





def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


if __name__ == '__main__':

    mainPath = r'C:/Users/iason/Desktop/traceanalysis/trace_analysis/traces/test_data/DNA04'

    exp = Experiment(mainPath)
    file = exp.files[0]
    filename = mainPath + '/hel12_dwells_data.xlsx'
    dwell_data =  pd.read_excel(filename, index_col=[0, 1], dtype={'kon' :np.str})
    fret_hist = plot_fret_histogram(file.molecules, bins=100)

