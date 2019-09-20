# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 23:19:30 2019

@author: iason
"""
import os
import numpy as np
import sys
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="dark")
sns.set_color_codes()

import traceAnalysisCode as analysis



def plot_dwells(dwells_dataframe, dwelltype='dwelltime', save=True):
    data = dwells_dataframe
    dwells = data[dwelltype].values
    dwells = dwells[~np.isnan(dwells)]

    dwells_l = data[dwelltype][data.side == 'l'].values
    dwells_m = data[dwelltype][data.side == 'm'].values
    dwells_r = data[dwelltype][data.side == 'r'].values
    dwells_r = dwells_r[~np.isnan(dwells_r)]

    tau_all = np.average(dwells)

    tau_l = np.average(dwells_l)
    tau_m = np.average(dwells_m)
    tau_r = np.average(dwells_r)


    t = file.time

    plt.figure(figsize=(10,5))
    colors = ['b', 'y', 'g', 'r']
    print('{:.1f}, {:.1f}, {:.1f}, {:.1f})'.format(tau_all, tau_l, tau_m, tau_r))
    for tau, d, lab, c in zip([tau_all, tau_l, tau_m, tau_r],
                           [dwells, dwells_l, dwells_m, dwells_r],
                           ['all', 'l', 'm', 'r'], colors):

        values, bins = np.histogram(d, bins=50, density=True)
        centers = (bins[1:] + bins[:-1]) / 2.0
        plt.semilogy(centers, values, '.', c=c, label=fr'$\tau_{lab} = ${tau:.1f}' )
        tau = np.average(d)
        exp = 1/tau*np.exp(-t/tau)
        plt.plot(t, exp, c=c)
        plt.xlim((0, 150))

    plt.legend( prop={'size': 16})
    plt.xlabel('time (s)')
    plt.ylabel('Prob.')
    plt.title(f'{dwelltype} histogram')





os.chdir(os.path.dirname(os.path.abspath(__file__)))
#mainPath = './traces'
mainPath = './traces'
exp = analysis.Experiment(mainPath, 0.1)
file = exp.files[0]
filename = './hel2_dwells_data.xlsx'
data = pd.read_excel(filename, index_col=[0,1], dtype={'kon':np.str})

dwelltype = 'dwelltime'
plot_dwells(data, dwelltype=dwelltype)
plt.savefig(f'{file.name}_{dwelltype}_hist.png', facecolor=None, dpi=300)