# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 23:19:30 2019

@author: iason

modified by Pim Sun Sept 22 2019
"""
#####################################################################
#Code to plot histograms of dwelltimes from all files in directory: 
#count vs dwelltime and log(Prob.) vs dwelltime
#####################################################################
import os
import numpy as np
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="dark")
sns.set_color_codes()

import traceAnalysisCode as analysis


def plot_dwells(file, dwells_array, dwelltype='offtime',nbins=20, save=True):
    data = dwells_array
    dwells = data[dwelltype].values
    dwells = dwells[~np.isnan(dwells)]
    ndwells=len(dwells)

    dwells_l = data[dwelltype][data.side == 'l'].values
    dwells_m = data[dwelltype][data.side == 'm'].values
    dwells_r = data[dwelltype][data.side == 'r'].values
    dwells_r = dwells_r[~np.isnan(dwells_r)]

    tau_all = np.average(dwells)

    tau_l = np.average(dwells_l)
    tau_m = np.average(dwells_m)
    tau_r = np.average(dwells_r)


    t = file.time
        
    plt.figure(1, facecolor='white')
    plt.hist(dwells,bins=nbins)
    plt.title(f'Histogram {dwelltype} bins:{nbins} dwells:{ndwells}')
    plt.xlabel(f'{dwelltype} (sec)')
    plt.ylabel('count')
    plt.savefig(f'{file.name}_histogram.png', facecolor='white', dpi=300)

    plt.figure(num=2,figsize=(10,5), facecolor='white')
    colors = ['b', 'y', 'g', 'r']
    print('tau_all={:.1f}, tau_l={:.1f}, tau_m={:.1f}, tau_r={:.1f}'.format(tau_all, tau_l, tau_m, tau_r))
    for tau, d, lab, c in zip([tau_all, tau_l, tau_m, tau_r],
                           [dwells, dwells_l, dwells_m, dwells_r],
                           ['all', 'l', 'm', 'r'], colors):

        values, bins = np.histogram(d, bins=nbins, density=True)
        centers = (bins[1:] + bins[:-1]) / 2.0
        plt.semilogy(centers, values, '.', c=c, label=fr'$\tau_{lab} = ${tau:.1f}' )
        tau = np.average(d)
        exp = 1/tau*np.exp(-t/tau)
        plt.plot(t, exp, c=c)
        #plt.xlim((0, 150))

    plt.legend(prop={'size': 16})
    plt.xlabel(f'{dwelltype} (sec)')
    plt.ylabel('log(Prob.)')
    plt.title(f'log(Prob.) vs {dwelltype} (s) bins:{nbins} dwells:{ndwells}')

if __name__ == '__main__':
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    mainPath = './traces'
    dwelltype = 'offtime'
    nbins=100

    exp = analysis.Experiment(mainPath)
    file = exp.files[0]
    filename = './'+file.name+'_dwells_data.xlsx'
    print(filename)
    data=pd.read_excel(filename, index_col=[0,1], dtype={'kon':np.str})

    if len(exp.files)>1:
        for file in exp.files[0:]:
            filename = './'+file.name+'_dwells_data.xlsx'
            print(filename)
            data2=(pd.read_excel(filename, index_col=[0,1], dtype={'kon':np.str}))
            data=data.append(data2,ignore_index=True)

    plot_dwells(file, data, dwelltype, nbins)
    if len(exp.files)>1:
        plt.savefig(f'{len(exp.files)}files_{dwelltype}_loghist.png', facecolor='white', dpi=300)
    else:
        plt.savefig(f'{file.name}_{dwelltype}_loghist.png', facecolor='white', dpi=300)