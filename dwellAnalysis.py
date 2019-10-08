# -*- coding: utf-8 -*-
"""
Created on Tue May 28 12:05:13 2019

@author: iason
"""
import time as timetime
import numpy as np
import sys
#import seaborn as sns
#sns.set(style="dark")
#sns.set_color_codes()

sys.path.append('..')
import pandas as pd
import os
#import matplotlib.pyplot as plt
import traceAnalysisCode as analysis


def analyze_dwelltimes(exp_file, save=True, filename=None):
    exp_file.importExcel()  # this should not be needed normally
    max_time = exp_file.time[-1]
    for i, mol in enumerate(exp_file.molecules):
        if mol.steps is None:
            continue
        times = mol.steps.time.sort_values().values
        try:
            times1 = times.reshape((int(times.size/2), 2))
        except ValueError:
            print(times)
        dwells = np.diff(times1, axis=1).flatten()

        labels = []
        for i, d in enumerate(dwells):
            lab = 'm'
            if times[0] == 0 and i == 0:  # first loop
                lab = 'l'
            if max_time - times[-1] < 0.1  and i == len(dwells) - 1:  # last loop
                lab = 'r'

            labels.append(lab)
        dwells = pd.DataFrame({'offtime': dwells, 'side': labels} )


        # Calculate the on times
        ontimes = []
        labels = []
        if times[0] != 0:  # append the left kon if it exists
            ontimes.append(times[0])
            labels.append('l')


        #for i, t in zip(range(1,times.size, 2), times[2:-1:2]):
        for i in range(2,times.size, 2):
            ontimes.append(times[i] - times[i-1])
            labels.append('m')

        if max_time - times[-1] > 0.1:  # append the right kon if it exists
            ontimes.append(max_time - times[-1])
            labels.append('r')

        ontimes = pd.DataFrame({'ontime': ontimes,
                            'onside': labels,
                            'order': np.arange(0, len(ontimes))} )

        mol.steps = pd.concat([mol.steps, dwells, ontimes], axis=1)

    if save:
        data = exp_file.savetoExcel(filename=exp_file.name+'_dwells_data.xlsx')
    else:
        data = exp_file.savetoExcel(save=False)
    return data

if __name__ == '__main__':

    start = timetime.time()
#    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    mainPath = 'G:/SM-data/20190923_dcas9_DNA10-11-20/#4.11_streptavidin_1nM_dcas9-crRNA-Cy5_8nM_DNA10-Cy3_G_movies_photobleach_0.3exp.time_coloc'
    #mainPath = './iasonas/cas9_simulations/DNA11-20_30exposure'
    exp = analysis.Experiment(mainPath)
    file = exp.files[1]
    file.exposure_time = 0.3
    data = analyze_dwelltimes(file, save=True)
#    for file in exp.files:
##    file = exp.files[1]
#
#        df = file.importExcel()
#        #
#        #
#        #
#        data = analyze_dwelltimes(file, save=True)
#        #
#
#
    print(f'Analysis time: {timetime.time() - start} sec')




