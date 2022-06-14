import os
import numpy as np
import xarray as xr
import sys
import pandas as pd
import matplotlib.pyplot as plt
import scipy.optimize
import seaborn as sns

def dwell_frames_from_classification(classification):
    # This assumes continuous monitoring with a specific cycle time.
    single_true_array = np.ones((classification.shape[0],1)).astype(bool)
    is_state_transition = np.hstack([single_true_array, classification[:,:-1] != classification[:,1:], single_true_array])
    state_transition_molecule, state_transition_frame = np.where(is_state_transition)
    is_endpoint = state_transition_frame == classification.shape[1]
    dwell_states = classification[state_transition_molecule[~is_endpoint], state_transition_frame[~is_endpoint]]
    dwell_frames = np.diff(state_transition_frame)[~is_endpoint[:-1]]
    dwell_molecules = state_transition_molecule[~is_endpoint]
    return dwell_molecules, dwell_states, dwell_frames

def determine_dwell_means(traces_flattened, dwell_frames):
    mean_trace = np.mean(traces_flattened)
    values_cumsum = np.concatenate([[0], np.cumsum(traces_flattened-mean_trace)])
    oneD_indices = np.concatenate([[0],dwell_frames.cumsum()])
    dwell_means = np.diff(values_cumsum[oneD_indices]) / dwell_frames + mean_trace
    return dwell_means

def set_states(dwell_molecules, dwell_states, at_trace_edges=True, around_negative_states=True, to_state=-2):
    states_to_set = np.zeros(len(dwell_states)).astype(bool)

    switched_molecule = np.diff(dwell_molecules).astype(bool)
    start_and_end_trace = np.concatenate([[True], switched_molecule]) | np.concatenate([switched_molecule, [True]])

    negative_states = dwell_states < 0

    if at_trace_edges:
        states_to_set |= start_and_end_trace

    if around_negative_states:
        negative_state_neighbours = np.concatenate([[False], negative_states[:-1]]) | \
                                    np.concatenate([negative_states[1:], [False]])
        states_to_set |= negative_state_neighbours & ~start_and_end_trace

    states_to_set[negative_states] = False

    dwell_states[states_to_set] = to_state
    return dwell_states


def dwell_times_from_classification(classification, traces=None, cycle_time=None, inactivate_start_and_end_states=True):
    if isinstance(classification, xr.DataArray):
        classification = classification.values
    dwell_molecules, dwell_states, dwell_frames = dwell_frames_from_classification(classification)
    if inactivate_start_and_end_states:
        dwell_states = set_states(dwell_molecules, dwell_states)

    ds = xr.Dataset({'molecule': ('dwell', dwell_molecules),
                     'state': ('dwell', dwell_states),
                     'frame_count': ('dwell', dwell_frames)})

    if cycle_time is not None:
        dwell_times = dwell_frames * cycle_time
        ds['duration'] = xr.DataArray(dwell_times, dims=['dwell'])

    if traces is not None:
        if isinstance(traces, xr.DataArray):
            traces = traces.values
        dwell_means = determine_dwell_means(traces.flatten(), dwell_frames)
        ds['mean'] = xr.DataArray(dwell_means, dims=['dwell'])

    return ds


def single_decaying_exponential(t, A, tau):
    return A * np.exp(-t/tau)

def analyze_dwells(dwells, fit_function, cycle_time=1, plot=False, axes=None):
    states = np.unique(dwells.state)
    positive_states = states[states>=0]

    bins=50

    if plot and axes is None:
        fig, axes = plt.subplots(1,len(positive_states), figsize=(len(positive_states)*3, 2), tight_layout=True, sharey=True)

    fit_values = {}
    for i, state in enumerate(positive_states):
        dwells_with_state = dwells.sel(dwell=dwells.state==state)
        c, t_edges = np.histogram(dwells_with_state.duration, bins=bins+1, range=[-cycle_time/2, (bins+1/2)*cycle_time])
        t = (t_edges[:-1]+t_edges[1:])/2
        popt, pcov = scipy.optimize.curve_fit(fit_function, t[1:], c[1:])

        if plot:
            axes[i].bar(t, c, width=cycle_time)
            axes[i].plot(t, fit_function(t, *popt), c='r')
            axes[i].set_title(str(np.round(popt[-1],3))+ ',' + str(np.round(dwells_with_state['mean'].mean().item(),3)))

        fit_values[state] = popt

    return fit_values



if __name__ == '__main__':
    classification = np.array([-1,-1,2,2,2,2,1,1,2,2,2,1,1,1,-1,-1,0,0,0,0,0,0,0,0,0,0,2,2,0,0,0,0,0])
    classification = np.repeat(classification[None,:], 5, axis=0)
    traces = np.random.random(classification.shape)

    dwell_times_from_classification(classification, traces=traces, cycle_time=0.1)

#
#
# classification = np.array([-1,-1,2,2,2,2,1,1,2,2,2,1,1,1,2,2,0,0,0,0,0,0,0,0,0,0,-1,-1,0,0,0,0,0])
# trace = np.random.random(classification.shape)
#
#
# def dwell_times_from_classification(classification, cycle_time=1):
#     # This assumes continuous monitoring with a specific cycle time.
#     is_state_transition = np.concatenate([[True], classification[:-1] != classification[1:], [True]])
#     state_transition_indices = np.where(is_state_transition)[0]
#     dwell_states = classification[state_transition_indices[:-1]]
#     dwell_frames = np.diff(state_transition_indices)
#     dwell_times = dwell_frames * cycle_time
#
#     mean_trace = np.mean(trace)
#     values_cumsum = np.concatenate([[0], np.cumsum(trace-mean_trace)])
#     dwell_means = np.diff(values_cumsum[state_transition_indices]) / dwell_frames + mean_trace
#
#     return dwell_states, dwell_times, dwell_means
#
#
# dwell_times_from_classification(classification)














# if __name__ == '__main__':
#     import SAfitting
#     import common_PDF
# else:
#     from trace_analysis.analysis import SAfitting
#     from trace_analysis.analysis import common_PDF
# # import SAfitting
# sns.set(style="ticks")
# sns.set_color_codes()
#
#
# def analyze(dwells_data, dataset_name, dist, configuration):
#     conf = configuration
#     # find the Tmax until which data is selected
#     d = apply_config_to_data(dwells_data, dist, conf)
#     figures = []
#     fit_data = []
#     keys_with_data = []  # keys refer to 'red', 'green', 'total', 'FRET'
#     for key in d.keys():
#         if d[key].empty:  # check if the dataframe is empty
#             print(f'{dist} dataFrame for {key} is empty')
#             continue
#         dwells = d[key].loc[:,dist].values
#         dwells = dwells[dwells>0]
#         print(np.size(dwells), 'dwells selected')
#         if conf['FitBool']:
#             fit_res = fit(dwells, model=conf['model'], dataset_name=dataset_name,
#                           Nfits=int(conf['Nfits']),
#                           include_over_Tmax=conf['TmaxBool'],
#                           bootstrap=conf['BootBool'],
#                           boot_repeats=int(conf['BootRepeats']))
#             fit_data.append(fit_res)
#         else:
#             fit_res = None
#         print(f'plotting {key} {dist}')
#         figure = plot(dwells, dataset_name, dist, trace=key, binsize=conf['binsize'],
#                       scale=conf['scale'], style=conf['PlotType'],
#                       fit_result=fit_res)
#         figures.append(figure)
#         keys_with_data.append(key)
#
#     if fit_data != []:
#         fit_data = pd.concat(fit_data, axis=1, keys=keys_with_data)
#     return d, figures, fit_data
#
# def fit(dwells, model='1Exp', dataset_name='Dwells', Nfits=1,
#         include_over_Tmax=True, bootstrap=True, boot_repeats=100):
#
#     if model == '1Exp+2Exp':
#         fit_result = []
#         for model in ['1Exp', '2Exp']:
#             result, boots = SAfitting.fit(dwells, model, dataset_name, Nfits,
#                                    include_over_Tmax, bootstrap, boot_repeats)
#             fit_result.append(result)
#         fit_result = pd.concat(fit_result, axis=1, ignore_index=True)
#         return fit_result
#
#     fit_result, boots = SAfitting.fit(dwells, model, dataset_name, Nfits, include_over_Tmax,
#                                   bootstrap, boot_repeats)
#     # print(fit_result)
#     return fit_result
#
#
# def plot(dwells, name, dist='offtime', trace='red', binsize='auto', scale='log',
#          style='dots', color='from_trace', fit_result=None):
#
#     if fit_result is not None:
#         if fit_result.Ncut[0] > 0:
#             Tcut = dwells.max() - 5  # 5 sec is kind of arbitrary here
#             dwells = dwells[dwells < Tcut]
#
#     try:
#         bsize = float(binsize)
#         if scale == 'Log-Log':
#             bin_edges = 10**(np.arange(np.log10(min(dwells)), np.log10(max(dwells)) + bsize, bsize))
#         else:
#             bin_edges = np.arange(min(dwells), max(dwells) + bsize, bsize)
#     except ValueError:
#         if binsize == 'Auto':
#             binsize = 'auto'
#         bin_edges = binsize
#     values, bins = np.histogram(dwells, bins=bin_edges, density=True)
#
#     # Determine position of bins
#     if scale == 'Log-Log':
#         centers = (bins[1:] * bins[:-1])**0.5  # geometric average of bin edges
#     else:
#         centers = (bins[1:] + bins[:-1]) / 2.0
#
#     # combine bins until they contain at least one data point (for y-log plots)
#     if scale in ['Log', 'Log-Log']:
#         izeros = np.where(values == 0)[0]
#         print('izeros', izeros)
#         j = 0
#         while j < len(izeros):
#             i = j
#             j += 1
#             while j < len(izeros) and izeros[j] - izeros[j-1] == 1:
#                 j += 1
#             # print('jstart ', izeros[i])
#             # print('jend ', izeros[i]+(j-i))
#             # print('values ', values[izeros[i]:(izeros[i]+j-i+1)])
#             # print('mean value', np.sum(values[izeros[i]:(izeros[i]+j-i+1)])/(j-i+1))
#             values[izeros[i]:(izeros[i]+j-i+1)] = np.sum(values[izeros[i]:(izeros[i]+j-i+1)])/(j-i+1)
#
#     fig = plt.figure(f'Histogram {trace} {dist}s {name}', figsize=(4, 3), dpi=200)
#
#     if color == 'from_trace':
#         if dist == 'offtime':
#             color = 'r'*(trace == 'red') + 'g'*(trace == 'green') + \
#                     'b'*(trace == 'FRET') + 'sandybrown'*(trace == 'total')
#         if dist == 'ontime':
#             color = 'firebrick'*(trace == 'red') + 'olive'*(trace == 'green') + \
#                     'darkviolet'*(trace == ' FRET') + 'saddlebrown'*(trace == 'total')
#     label = f'{dist} pdf, N={dwells.size}'
#     if style == 'dots':
#         plt.plot(centers, values, '.', color=color, label=label)
#     if style == 'bars':
#         plt.bar(centers, values, color=color, label=label,
#                 width=(bins[1] - bins[0]))
#     if style == 'line':
#         plt.plot(centers, values, '-', lw=2, color=color, label=label)
#
#     if fit_result is not None:
#         if fit_result.model[0] == '1Exp':
#             tau = fit_result.value[0]
#             error = fit_result.error[0]
#             Ncut = fit_result.Ncut[0]
#             print(f'plotting 1Exp fit')
#             time, fit = common_PDF.Exp1(tau,
#                                         Tmax=centers[-1]+(bins[1]-bins[0])/2)
#             label = f'\n tau={tau:.1f}'
#             if error != 0:
#                 label += f'$\pm$ {error:.1f}'
# #            plt.plot(time, fit, color='r', label=f'1expFit, Ncut={int(Ncut)} \n {label}')
#
#         elif fit_result.model[0] == '2Exp':
#             p, errp = fit_result.value[0], fit_result.error[0]
#             tau1, err1 = fit_result.value[1], fit_result.error[1]
#             tau2, err2 = fit_result.value[2], fit_result.error[2]
#             Ncut = fit_result.Ncut[0]
#             print(fit_result)
#             print(f'errors: ', errp, err1, err2)
#             time, fit = common_PDF.Exp2(p, tau1, tau2, Tmax=centers[-1])
#             label = f'\n p={p:.2f}, tau1={tau1:.1f}, tau2={int(tau2)}'
# #            plt.plot(time, fit, color='r', label=f'2expFit, Ncut={int(Ncut)} \n {label}')
#
#         elif fit_result.model[0] == '3Exp':
#             p1, errp1 = fit_result.value[0], fit_result.error[0]
#             p2, errp2 = fit_result.value[1], fit_result.error[1]
#             tau1, err1 = fit_result.value[2], fit_result.error[2]
#             tau2, err2 = fit_result.value[3], fit_result.error[3]
#             tau3, err3 = fit_result.value[4], fit_result.error[4]
#             Ncut = fit_result.Ncut[0]
#             print(fit_result)
#             print(f'errors: ', errp1, errp2, err1, err2, err3)
#             time, fit = common_PDF.Exp3(p1, p2, tau1, tau2, tau3,
#                                         Tmax=centers[-1])
#             label = f'\n p1={p1:.2f}, p2={p2:.2f}, tau1={tau1:.1f}, tau2={int(tau2)}, tau3={int(tau3)}'
# #            plt.plot(time, fit, color='r', label=f'3expFit, Ncut={int(Ncut)} \n {label}')
#
#         if fit_result.Ncut[0] > 0:
#             label = f', Ncut={int(Ncut)}' + label
#
#         plt.plot(time, fit, color='k', label=f'{fit_result.model[0]}fit{label}')
#
#     if scale in ['Log', 'Log-Log']:
#         plt.yscale('log')
#
#     if scale == 'Log-Log':
#         plt.xscale('log')
#
#     plt.legend()
#     plt.ylabel('Probability')
#     plt.xlabel(f'{dist} (s)')
#     # plt.locator_params(axis='y', nbins=3)
#     plt.tight_layout()
#     plt.show()
#     return fig
#
# def apply_config_to_data(dwells_data, dist, config):
#     d = dwells_data
#     # Select the requested sides
#     side_list = ['l'*bool(config['side']['left']),
#                'm'*bool(config['side']['middle']),
#                'r'*bool(config['side']['right'])]
#
#     if dist == 'offtime':
#         d = d[d.side.isin(side_list)]
#     if dist == 'ontime':
#         d = d[d.onside.isin(side_list)]
#     # apply min, max conditions
#     if config['max'] in ['Max', 'max']:
#         d = d[d[dist] > float(config['min'])]
#     else:
#         d = d[d[dist] > float(config['min'])]
#         d = d[d[dist] < float(config['max'])]
#
#     data = {}
#
#     for key in config['trace'].keys():
#         if config['trace'][key]:
#             data[key] = d[d['trace'] == key]
#         else:
#             pass
#
#     return data
#
#
# if __name__ == '__main__':
#     filename = 'C:/Users/iason/Desktop/traceanalysis/trace_analysis/traces/'
#     filename += 'hel0_dwells_data.xlsx'
#
#     data = pd.read_excel(filename, index_col=[0, 1], dtype={'kon' :np.str})
#     print(data.shape)
#     config = {'trace': {'red': True, 'green': False, 'total': False, 'FRET': False},
#          'side': {'left': True, 'middle': True, 'right': True},
#          'min': '0', 'max': 'max',
#          'scale': 'Normal',
#          'PlotType': 'dots',
#          'binsize': 'auto',
#          'FitBool': True,
#          'TmaxBool': False,
#          'BootBool': False,
#          'model': '2Exp',
#          'Nfits': '1',
#          'BootRepeats': '5'}
#
#     result = analyze(data, 'test', 'offtime', config)
