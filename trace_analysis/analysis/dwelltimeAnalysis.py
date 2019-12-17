# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 11:32:58 2019

@author: ikatechis
"""

import os
import numpy as np
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="ticks")
sns.set_color_codes()

#from trace_analysis import Experiment

def analyze(dwells_data, dist, configuration):
    conf = configuration
    d = apply_config_to_data(dwells_data, dist, conf)

    figures = []
    for key in d.keys():

        if d[key].isempty:  # check if the dataframe is empty
            continue

        dwells = d[key].loc[dist].values
        figure = plot(dwells, trace=key, binsize=conf['binsize'],
                      scale=conf['scale'], style=conf['PlotType'],
                      fit_data=None)
        figures.append(figure)

    return d, figures

def plot(dwells, trace='red', binsize='auto', scale='log',
         style='dots', fit_data=None):
    pass



def apply_config_to_data(dwells_data, dist, config):
    d = dwells_data
    # Select the requested sides
    side_list = ['l'*bool(config['side']['left']),
               'm'*bool(config['side']['middle']),
               'r'*bool(config['side']['right'])]

    if dist == 'offtime':
        d = d[d.side.isin(side_list)]
    if dist == 'ontime':
        d = d[d.onside.isin(side_list)]
    # apply min, max conditions
    if config['max'] in ['Max', 'max']:
        d = d[d[dist] > float(config['min'])]
    else:
        d = d[d[dist] > float(config['min']) & d[dist] < float(config['max'])]

    data = {}

    for key in config['trace'].keys():
        if config['trace'][key]:
            data[key] = d[d['trace'] == key]
        else:
            pass

    return data




def get_average_dwell(dwells):
    #  correct for dwell exceeding the measurement time minus 5 sec
    Tmax = dwells.max() - 5
    Ntot = dwells.size
    Ncut = dwells[dwells > Tmax].size
    avrg_dwell = np.average(dwells[dwells < Tmax])
    avrg_dwell = avrg_dwell + Ncut*Tmax/Ntot
    return avrg_dwell

if __name__ == '__main__':
    filename = 'F:/Google Drive/PhD/Programming - Data Analysis/traceanalysis/traces/'
    filename += 'hel0_dwells_data.xlsx'

    data = pd.read_excel(filename, index_col=[0, 1], dtype={'kon' :np.str})

    config = {'trace': {'red': True, 'green': False, 'total': False, 'FRET': True},
         'side': {'left': True, 'middle': True, 'right': True},
         'min': '0', 'max': 'Max',
         'scale': 'Normal',
         'PlotType': 'Dots',
         'binsize': 'Auto',
         'FitBool': False,
         'TmaxBool': False,
         'BootBool': False,
         'model': '1Exp',
         'Nfits': '10',
         'BootRepeats': '100'}

    result = analyze(data, 'offtime', config)

    #  select only the ones that don't exceed the total measurement time minus 10 sec
#    dwells_in = dwells[dwells < dwells.max() - 10]
    # avrg_dwell = get_average_dwell(dwells)

    # values, bins = np.histogram(dwells, bins=nbins, density=True)
    # centers = (bins[1:] + bins[:-1]) / 2.0
    # if not plot:
    #     return values, centers

    # if plot:
    #     plt.figure(dwelltype)
    #     line = plt.plot(centers, values, '.',
    #                     label=extra_label+fr'$\tau = ${avrg_dwell:.1f} s' )[0]

    #     plt.xlabel('time (s)')
    #     plt.ylabel('Prob.')
    #     plt.title(f'{dwelltype} histogram: nbins={nbins} N={dwells.size}')
    #     plt.legend(prop={'size': 16})
    #     # plot a 1exp ML 'fit' for the average dwelltime
    #     t = np.arange(0, dwells.max(), 0.1)
    #     exp = 1/avrg_dwell*np.exp(-t/avrg_dwell)
    #     if log:
    #         plt.semilogy(t, exp, color=line.get_color())
    #     else:
    #         plt.plot(t, exp, color=line.get_color())


    #     return line