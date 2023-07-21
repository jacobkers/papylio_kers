import os
import numpy as np
import xarray as xr
import sys
import pandas as pd
import matplotlib.pyplot as plt
import scipy.optimize
import inspect
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
    states_to_set = np.zeros(len(dwell_states), dtype=bool)

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
        molecule_coords = {n: c.values for n, c in classification.coords.items() if c.dims[0] == 'molecule' and len(c.dims) == 1}
        classification = classification.values
    else:
        molecule_coords = None
    dwell_molecules, dwell_states, dwell_frames = dwell_frames_from_classification(classification)
    if inactivate_start_and_end_states:
        dwell_states = set_states(dwell_molecules, dwell_states, to_state=-128)
        # Probably better to indicate for each dwell whether it is at a trace edge or around a negative state

    ds = xr.Dataset()
    if molecule_coords is not None:
        for n, c in molecule_coords.items():
            ds[n] = ('dwell', c[dwell_molecules])
    else:
        ds['molecule'] = ('dwell', dwell_molecules)

    ds['state'] = ('dwell', dwell_states)
    ds['frame_count'] = ('dwell', dwell_frames)

    if cycle_time is not None:
        dwell_times = dwell_frames * cycle_time
        ds['duration'] =('dwell', dwell_times)

    if traces is not None:
        # if isinstance(traces, xr.Dataset):
        #     for name, da in traces.data_vars.items():
        #         dwell_means = determine_dwell_means(da.values.flatten(), dwell_frames)
        #         ds['mean'] = xr.DataArray(dwell_means, dims=['dwell'])
        name = ''
        if isinstance(traces, xr.DataArray):
            if traces.name is not None:
                name = '_' + traces.name
            traces = traces.values

        dwell_means = determine_dwell_means(traces.flatten(), dwell_frames)
        ds['mean'+name] = ('dwell', dwell_means)

    return ds


# def single_decaying_exponential(t, A, tau):
#     return A * np.exp(-t/tau)
# single_decaying_exponential.bounds = ((0,0),(np.inf,np.inf))
# def p0(t, y):
#     return y.max(), t.mean()
# single_decaying_exponential.p0 = p0

def single_decaying_exponential(t, A, k):
    return A * np.exp(-k*t)
single_decaying_exponential.bounds = ((0,0),(np.inf,np.inf))
def p0(t, y):
    return y.max(), 1/t.mean()
single_decaying_exponential.p0 = p0

def analyze_dwells(dwells, fit_function=single_decaying_exponential, cycle_time=1, plot=False,
                   axes=None, state_names={0: 'Low FRET state', 1: 'High FRET state'}, logy=False):
    # states = np.unique(dwells.state)
    states = np.array(list(state_names.keys()))
    positive_states = states[states>=0]

    bins=50

    if plot and axes is None:
        fig, axes = plt.subplots(1,len(positive_states), figsize=(len(positive_states)*3, 2), layout='constrained', sharey=True)
    else:
        axes = None

    fit_parameters = list(inspect.signature(fit_function).parameters)[1:]
    fit_values = xr.Dataset(coords={'state': positive_states, 'parameter': fit_parameters})
    fit_values['optimal_value'] = xr.DataArray(np.nan, dims=('state', 'parameter'), coords={'state': positive_states, 'parameter': fit_parameters})
    fit_values['error'] = xr.DataArray(np.nan, dims=('state', 'parameter'), coords={'state': positive_states, 'parameter': fit_parameters})
    fit_values['covariance'] = xr.DataArray(np.nan, dims=('state', 'parameter','parameter'),
                                       coords={'state': positive_states, 'parameter': fit_parameters, 'parameter': fit_parameters})

    fit_values.attrs['fit_function'] = fit_function.__name__

    for i, state in enumerate(positive_states):
        dwells_with_state = dwells.sel(dwell=dwells.state==state)
        c, t_edges = np.histogram(dwells_with_state.duration, bins=bins+1, range=[-cycle_time/2, (bins+1/2)*cycle_time])
        t = (t_edges[:-1]+t_edges[1:])/2
        try:
            popt, pcov = scipy.optimize.curve_fit(fit_function, t[2:], c[2:], p0=fit_function.p0(t[2:], c[2:]),
                                                  bounds=fit_function.bounds, absolute_sigma=True)
            perr = np.sqrt(np.diag(pcov))
            fit_values['optimal_value'][dict(state=state)] = popt
            fit_values['covariance'][dict(state=state)] = pcov
            fit_values['error'][dict(state=state)] = perr

            # fit_values[state] = {fit_parameter: {'value': value, 'error': error} for fit_parameter, value, error in
            #                      zip(fit_parameters, popt, perr)}

        except RuntimeError:
            popt = pcov = perr = None

        if plot:
            # TODO: make this a separate function
            axes[i].bar(t, c, width=cycle_time)
            # axes[i].set_title(+ ',' + str())
            if popt is not None:
                axes[i].plot(t, fit_function(t, *popt), c='r', label=r'$count = Ae^{-\frac{t}{\tau}}$')
                text_string = f'A={popt[0]:.5}\nτ={popt[-1]:.5}\n1/τ={1/popt[-1]:.5}' #\nmean={dwells_with_state["mean"].mean().item():.5}'
                text_string = f'A={popt[0]:.0f}±{perr[0]:.0f}\nk={popt[-1]:.2f}±{perr[-1]:.2f}\nτ={1 / popt[-1]:.1f}'  # \nmean={dwells_with_state["mean"].mean().item():.5}'
            else:
                text_string = 'No fit found'

            text_string = state_names[state] + '\n' + text_string

            # axes[i].text(0.95, 0.95, text_string, transform=axes[i].transAxes, fontsize=14,
            #     verticalalignment='top')# , bbox=props)
            axes[i].annotate(text_string, (0.55, 0.9), xycoords='axes fraction', fontsize=9,
                        verticalalignment='top')


            if axes[i].get_subplotspec().is_last_row():
                axes[i].set_xlabel('Dwell time (s)')
            axes[i].set_xlim([-cycle_time/2,30*cycle_time])
            if axes[i].get_subplotspec().is_first_col():
                # axes[i].set_ylim(0,c[1:].max()*1.1)
                axes[i].set_ylabel('Dwell count')
                if logy:
                    axes[i].set_yscale('log')
                    axes[i].set_ylim(0.5, axes[i].get_ylim()[1])
            # if axes[i].get_subplotspec().is_last_col():
            #     axes[i].legend()

    return fit_values, axes



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
