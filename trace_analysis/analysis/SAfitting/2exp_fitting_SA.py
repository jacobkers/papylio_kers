# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 18:59:11 2019

@author: pimam
"""
import numpy as np
import matplotlib.pylab as plt
import os
import SA_simple as SA


def P(dwells, params):
    P1, tau1, tau2 = params
    return P1/tau1*np.exp(-dwells/tau1)+(1-P1)/tau2*np.exp(-dwells/tau2)


def LogLikeLihood(xdata, params, model):
    Pi = model(xdata, params)
    LLike = np.sum(-np.log(Pi))
    return LLike


def update_temp(T, alpha):
    '''
    Exponential cooling scheme
    :param T: current temperature.
    :param alpha: Cooling rate.
    :return: new temperature
    '''
    T *= alpha
    return T


def simmulated_annealing(data, objective_function, model, x_initial, lwrbnd, uprbnd, Tstart=200., Tfinal=0.001, delta=0.05, alpha=0.9):

    T = Tstart
    step = 0
    x = x_initial
    while T > Tfinal:
        step += 1
        if (step % 100 == 0):
            T = update_temp(T, alpha)
        print(x)
        x = np.log(x)
        x_trial = np.zeros(len(x))
        # x_trial = x + np.random.uniform(-delta, delta, size=len(x))
        for i in range(0, len(x)):
            x_trial[i] = np.random.uniform(np.max([x[i] - delta, lwrbnd[i]]),
                                           np.min([x[i] + delta, uprbnd[i]]))
        x_trial = np.exp(x_trial)
        x = np.exp(x)
        x = Metropolis(objective_function, model, x, x_trial, T, data)
    return x


def Metropolis(f, model, x, x_trial, T, data):
    # Metropolis Algorithm to decide if you accept the trial solution.
    Vnew = f(data, x_trial, model)
    Vold = f(data, x, model)
    if (np.random.uniform() < np.exp(-(Vnew - Vold) / T)):
        x = x_trial
    return x


if __name__ == '__main__':
    try:
        os.remove('./init_monitor.txt')
        os.remove('./monitor.txt')
    except IOError:
        print('monitor file not present')

    dwells_rec = np.load('./data/2exp_N=1000_rep=1_tau1=10_tau2=100_a=0.5.npy')
    dwells_cut = []
    Nfits = 5
    max_dwells = dwells_rec.max()
    avg_dwells = np.average(dwells_rec)
    x_initial = [0.5, avg_dwells, avg_dwells]
    lwrbnd = [0, 0, 0]
    uprbnd = [1, max_dwells, max_dwells]

    fitdata = simmulated_annealing(data=dwells_rec, objective_function=LogLikeLihood, model=P, x_initial=x_initial, lwrbnd=lwrbnd, uprbnd=uprbnd)
#    fitparams = SA.sim_anneal_fit(dwells_rec, dwells_cut, Earray, p_start, objective_function=LogLikeLihood, model=P, Tstart=100, Tfinal=0.001, delta=1, alpha=0.85)
    print("fit found: ", str(fitdata))
    fitparams = [fitdata]
    for i in range(1, Nfits):
        fitdata = simmulated_annealing(data=dwells_rec, objective_function=LogLikeLihood, model=P, x_initial=x_initial, lwrbnd=lwrbnd, uprbnd=uprbnd)
        print("fit found: ", str(fitdata))
        fitparams = np.concatenate((fitparams, [fitdata]), axis=0)

    plt.figure()
    values, bins = np.histogram(dwells_rec, bins=20, density=True)
    centers = (bins[1:] + bins[:-1]) / 2.0
    plt.plot(centers, values, '.', label='All dwells')

    timearray = np.linspace(0, max_dwells, num=200)
    for i in [ 1, 2, 3, 4]:#range(0, np.size(fitparams, 0)):
        fit = P(timearray, fitparams[i])
        plt.plot(timearray, fit, label='fit'+str(i))
    plt.legend()
