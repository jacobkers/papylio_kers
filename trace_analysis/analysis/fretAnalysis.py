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



def get_fret_per_offtime(molecule, trace_type='red', Emin=0.05, Emax=1,
                        Irmin=40, Iroff=0, Igoff=0):
    '''
    Function that gives the FRET trace with FRET calculcated only for the
    selected offtime interval. The FRET value outside this interval is put to
    zero
    '''
    if mol.steps is None:
        return
    time_axis = mol.file.time
    indices = np.empty(mol.steps.time.values.size, dtype=int)
    for i, t in enumerate(mol.steps.time.values):
        print(t)
        idx = find_nearest(time_axis, t)
        indices[i] = idx
    fret = molecule.E(Irmin, Iroff, Igoff)
    idx_start = indices[::2]
    idx_stop = indices[1::2]
    mask = np.zeros(fret.size, dtype=int)
    for i in range(idx_start.size):
        mask[idx_start[i]+1:idx_stop[i]-1] = 1

    fret_masked = np.multiply(fret, mask)
    return fret, mask


def find_nearest(array, value):
    # array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return int(idx)


if __name__ == '__main__':

    mainPath = r'C:/Users/iason/Desktop/traceanalysis/trace_analysis/traces/test_data/DNA04'

    exp = Experiment(mainPath)
    file = exp.files[0]
    filename = mainPath + '/hel12_dwells_data.xlsx'
    dwell_data =  pd.read_excel(filename, index_col=[0, 1], dtype={'kon' :np.str})
    # fret_hist = plot_fret_histogram(file.molecules, bins=100)
    mol = file.molecules[0]
    fret, mask = get_fret_per_offtime(mol)

