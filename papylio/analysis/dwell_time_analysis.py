import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import numbers

import papylio

class ExponentialDistribution:
    def __init__(self, number_of_exponentials, P_bounds=(-np.inf,np.inf), k_bounds=(1e-9,np.inf), truncation=None,
                 sampling_interval=None):
        self.number_of_exponentials = number_of_exponentials
        self.P_bounds = P_bounds
        self.k_bounds = k_bounds
        # self.bounds = np.array([(0, np.inf)] * (2 * self.number_of_exponentials - 1))
        self.bounds = np.array([P_bounds] * (self.number_of_exponentials - 1) + [k_bounds] * self.number_of_exponentials)
        self.truncation = truncation
        self.sampling_interval = sampling_interval

    def __call__(self, t, *parameters):
        return self.pdf(t, *parameters)

    @property
    def parameter_names(self):
        return ([f'P{i}' for i in range(self.number_of_exponentials-1)] +
                [f'k{i}' for i in range(self.number_of_exponentials)])

    @property
    def parameter_names_full(self):
        return ([f'P{i}' for i in range(self.number_of_exponentials)] +
                [f'k{i}' for i in range(self.number_of_exponentials)])

    # def normalize_P(self, P, k):
    #     # P = np.abs(P)
    #     P = np.array(P)
    #     k = np.array(k)
    #     # return P / np.sum(P * (1/k))
    #     return P

    def P_and_k_from_parameters(self, parameters):
        parameters = np.array(parameters)
        P = parameters[0:(self.number_of_exponentials-1)]
        P = np.hstack([P, [1-np.sum(P)]])
        k = parameters[(self.number_of_exponentials - 1):]
        return P, k

    def pdf(self, t, *parameters):
        # Parameters given as P1, P2, k1, k2, k3
        t = np.array(t)
        pdf = np.zeros_like(t).astype(float)
        P, k = self.P_and_k_from_parameters(parameters)
        for i in range(self.number_of_exponentials):
            pdf += P[i] * k[i] * np.exp(-k[i] * t)

        if self.truncation is not None:
            truncated_probability = (self.cdf_untruncated(self.truncation[1], *parameters) -
                                     self.cdf_untruncated(self.truncation[0], *parameters) + 1e-10)
            pdf = pdf / truncated_probability
            pdf[(t < self.truncation[0]) | (t > self.truncation[1])] = 0

        return pdf

    def cdf_untruncated(self, t, *parameters):
        t = np.array(t)
        cdf = np.ones_like(t).astype(float)
        P, k = self.P_and_k_from_parameters(parameters)
        for i in range(self.number_of_exponentials):
            cdf -= P[i] * np.exp(-k[i] * t)
        return cdf

    def cdf(self, t, *parameters):
        cdf = self.cdf_untruncated(t, *parameters)

        if self.truncation is not None:
            cdf_at_truncation = self.cdf_untruncated(self.truncation, *parameters)
            cdf = (cdf - cdf_at_truncation[0]) / (cdf_at_truncation[1] - cdf_at_truncation[0])
            cdf[(t < self.truncation[0])] = 0
            cdf[(t > self.truncation[1])] = 1

        return cdf

    def pdf_binned(self, t, *parameters):
        # Parameters given as P1, P2, k1, k2, k3
        # result = np.zeros_like(t).astype(float)
        # P, k = self.P_and_k_from_parameters(parameters)
        # for i in range(self.number_of_exponentials):
        #     result += 2 * np.sinh(k[i] * self.bin_width/2) * P[i] * np.exp(-k[i] * t)
        # result /= self.bin_width

        pdf_binned = self.cdf(t+self.bin_width/2, *parameters) - self.cdf(t-self.bin_width/2, *parameters)
        pdf_binned /= self.bin_width
        return pdf_binned

    def likelihood(self, parameters, t):
        return np.prod(self.pdf(t, *parameters))

    def loglikelihood(self, parameters, t):
        if self.truncation is not None:
            t = t[(t > self.truncation[0]) | (t < self.truncation[1])]

        probability_density = self.pdf(t, *parameters) + 1e-10
        # probability_density = self.pdf_binned(t, *parameters) + 1e-10

        loglikelihood = np.sum(np.log(probability_density))

        return loglikelihood

    def negative_loglikelihood(self, parameters, t):
        return -self.loglikelihood(parameters, t)

    def loglikelihood_binned(self, parameters, bin_centers, counts):
        if self.truncation is not None:
            bin_centers = bin_centers[(bin_centers > self.truncation[0]) | (bin_centers < self.truncation[1])]

        probability_density = self.pdf_binned(bin_centers, *parameters) + 1e-10
        loglikelihood = np.sum(counts * np.log(probability_density * self.bin_width))

        return loglikelihood

    def negative_loglikelihood_binned(self, parameters, bin_centers, counts):
        return -self.loglikelihood_binned(parameters, bin_centers, counts)

    def mle(self, *args, **kwargs):
        return self.maximum_likelihood_estimation(self, *args, **kwargs)

    def maximum_likelihood_estimation(self, dwell_times, scipy_optimize_method='minimize', **kwargs):
        scipy_optimize_kwargs = dict(bounds = self.bounds)

        if scipy_optimize_method in ['minimize', 'basinhopping', 'dual_annealing', 'differential_evolution']:
            scipy_optimize_kwargs['x0'] = self.parameter_guess(dwell_times)

        # def constraint(parameters):
        #     P, k = self.P_and_k_from_parameters(parameters)
        #     return np.sum(P * k)
        # if scipy_optimization_method in ['minimize']:
        #     scipy_optimize_kwargs['constraints'] = scipy.optimize.NonlinearConstraint(constraint, 0, np.inf)
        # elif scipy_optimization_method == 'shgo':
        #     scipy_optimize_kwargs['constraints'] = {'type': 'ineq', 'fun': constraint} #TODO: change in scipy v1.11.0

        scipy_optimize_kwargs.update(kwargs)

        optimal_parameters = getattr(scipy.optimize, scipy_optimize_method)(self.negative_loglikelihood,
                                     args = (dwell_times,), **scipy_optimize_kwargs).x
        BIC = self.BIC(dwell_times, optimal_parameters)

        dwell_analysis = self.parameters_to_dataset(optimal_parameters, BIC=BIC)
        dwell_analysis.attrs['fit_method'] =  'maximum_likelihood_estimation'
        dwell_analysis.attrs['scipy_optimize_method'] = scipy_optimize_method

        for key, item in scipy_optimize_kwargs.items():
            if isinstance(item, np.ndarray):
                item = item.tolist()
            dwell_analysis.attrs['scipy_optimize_'+key] = item
        # dwell_analysis.attrs['scipy_optimize_kwargs'] = scipy_optimize_kwargs

        return dwell_analysis

    def hist_fit(self, *args, **kwargs):
        return self.histogram_fit(*args, **kwargs)

    def histogram_fit(self, dwell_times, bins='auto_discrete', remove_first_bins=None, **kwargs):
        if bins == 'auto_discrete':
            bin_edges = auto_bin_size_for_discrete_data(dwell_times)
        else:
            bin_edges = np.histogram_bin_edges(dwell_times, bins=bins)

        if remove_first_bins is not None:
            bins = bins[remove_first_bins:]

        self.bin_width = bin_edges[1] - bin_edges[0]
        weights = (1 / len(dwell_times) / self.bin_width,) * len(dwell_times)  # Assuming evenly spaced bins
        counts, bin_edges = np.histogram(dwell_times, bins=bin_edges, density=False, weights=weights)

        scipy_curve_fit_kwargs = dict(p0 = self.parameter_guess(dwell_times),
                                      bounds = self.bounds.T, absolute_sigma=True)
        scipy_curve_fit_kwargs.update(kwargs)

        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        optimal_parameters, parameter_covariances = scipy.optimize.curve_fit(self.pdf_binned, bin_centers, counts,
                                                                            **scipy_curve_fit_kwargs)
        parameter_errors = np.sqrt(np.diag(parameter_covariances))
        BIC = self.BIC_histogram(bin_centers, counts*len(dwell_times)*self.bin_width, optimal_parameters)

        dwell_analysis = self.parameters_to_dataset(optimal_parameters, parameter_errors=parameter_errors, BIC=BIC)
        dwell_analysis.attrs['fit_method'] =  'histogram_fit'
        dwell_analysis.attrs['bins'] = bins
        dwell_analysis.attrs['remove_first_bins'] = remove_first_bins

        for key, item in scipy_curve_fit_kwargs.items():
            if isinstance(item, np.ndarray):
                item = item.tolist()
            dwell_analysis.attrs['scipy_curve_fit_'+key] = item

        return dwell_analysis
        #
        # fig, ax = plt.subplots()
        # ax.hist(t, bins=bins, density=True)
        # t_fit = np.linspace(0, bin_centers.max(), 1000)
        # t_hist = bin_centers
        # y_fit = self.pdf_binned(t_fit, *optimal_parameters)
        # y_hist = counts
        # ax.plot(t_fit, y_fit)
        # ax.plot(t_hist, y_hist)

    def cdf_fit(self, dwell_times, **kwargs):
        t, ecdf = empirical_cdf(dwell_times, self.sampling_interval)

        scipy_curve_fit_kwargs = dict(p0 = self.parameter_guess(dwell_times),
                                      bounds = self.bounds.T, absolute_sigma=True)
        scipy_curve_fit_kwargs.update(kwargs)

        optimal_parameters, parameter_covariances = scipy.optimize.curve_fit(self.cdf, t, ecdf,
                                                                             **scipy_curve_fit_kwargs)
        parameter_errors = np.sqrt(np.diag(parameter_covariances))
        BIC = self.BIC(dwell_times, optimal_parameters) #TODO: change or check BIC calculation

        dwell_analysis = self.parameters_to_dataset(optimal_parameters, parameter_errors=parameter_errors, BIC=BIC)
        dwell_analysis.attrs['fit_method'] =  'cdf_fit'

        for key, item in scipy_curve_fit_kwargs.items():
            if isinstance(item, np.ndarray):
                item = item.tolist()
            dwell_analysis.attrs['scipy_curve_fit_'+key] = item

        return dwell_analysis

    def parameter_guess(self, dwell_times):
        guess_P = [1 / (i + 1) for i in range(1, self.number_of_exponentials)]
        guess_k = [1 / (dwell_times.mean() * (i + 1)) for i in range(self.number_of_exponentials)]
        parameters = guess_P + guess_k
        return parameters

    def parameters_full(self, parameters):
        parameters = list(parameters)
        P = parameters[0:self.number_of_exponentials-1]
        parameters.insert(self.number_of_exponentials - 1, 1 - sum(P))
        return parameters

    def BIC(self, dwell_times, optimal_parameters):
        return len(optimal_parameters) * np.log(len(dwell_times)) - 2 * self.loglikelihood(optimal_parameters, dwell_times)

    def BIC_histogram(self, bin_centers, counts, optimal_parameters):
        k = len(optimal_parameters)
        N = np.sum(counts)
        BIC = k * np.log(N) - 2 * self.loglikelihood_binned(optimal_parameters, bin_centers, counts)
        return BIC

    # def parameters_to_dataframe(self, parameters, bic=True, dwell_times=None):
    #     dwell_analysis = {}
    #     dwell_analysis['Exponential'] = [self.number_of_exponentials] * self.number_of_exponentials #+
    #     dwell_analysis['Component'] = np.arange(self.number_of_exponentials)+1# [''] * (e-1)
    #     dwell_analysis['P'] = np.hstack([parameters[0:self.number_of_exponentials-1] , np.nan])
    #     dwell_analysis['P'][-1] = 1-np.nansum(dwell_analysis['P'])
    #     dwell_analysis['k'] = parameters[self.number_of_exponentials - 1:]
    #     # if parameter_errors is not None:
    #     #     dwell_analysis['P error'] = np.hstack([parameter_errors[0:e-1] ,np.nan])
    #     #     dwell_analysis['k error'] = parameter_errors[e - 1:]
    #     if bic and dwell_times is not None:
    #         dwell_analysis['BIC']  = [self.BIC(dwell_times, *parameters)] + [np.nan] * (self.number_of_exponentials-1)
    #
    #     return pd.DataFrame(dwell_analysis).set_index('Exponential')
    #
    # def dataframe_to_parameters(self, dataframe):
    #     parameters = dataframe.loc[self.number_of_exponentials][['P', 'k']].values.T.flatten()
    #     # return parameters[~np.isnan(parameters)]
    #     parameters = list(parameters)
    #     parameters.pop(self.number_of_exponentials-1)
    #     return parameters

    def parameters_to_dataset(self, parameters, parameter_errors=None, BIC=None):
        # coords = dict(parameter=pd.MultiIndex.from_product((['P','k'], np.arange(self.number_of_exponentials)), names=['name', 'component']))
        # coords = dict(parameter=('parameter', np.repeat(['P','k'],2)),
        #               component=('parameter', list(range(self.number_of_exponentials))*2))
        # parameters_da = xr.DataArray(self.parameters_full(parameters), dims=('parameter'),
        #                              coords=coords,
        #                              name='parameters').expand_dims('fit')
        #
        # bic = xr.DataArray([self.BIC(dwell_times, *parameters)], dims=('fit'))
        # dwell_analysis = xr.Dataset(dict(parameter=parameters, bic=bic))
        data_vars = {} # dict(P=P, k=k, BIC=BIC, fit_function=fit_function, fit_method=fit_method, number_of_components=number_of_components)
        P = self.parameters_full(parameters)[0:self.number_of_exponentials]
        k = self.parameters_full(parameters)[self.number_of_exponentials:]

        data_vars['P'] = xr.DataArray(P, dims=('component')).expand_dims('fit')
        data_vars['P'].attrs['units'] = ''

        data_vars['k'] = xr.DataArray(k, dims=('component')).expand_dims('fit')
        data_vars['k'].attrs['units'] = 's⁻¹'

        if parameter_errors is not None:
            P_error = self.parameters_full(parameter_errors)[0:self.number_of_exponentials]
            k_error = self.parameters_full(parameter_errors)[self.number_of_exponentials:]

            data_vars['P_error'] = xr.DataArray(P_error, dims=('component')).expand_dims('fit')
            data_vars['P_error'].attrs['units'] = ''

            data_vars['k_error'] = xr.DataArray(k_error, dims=('component')).expand_dims('fit')
            data_vars['k_error'].attrs['units'] = 's⁻¹'

        if BIC is not None:
            data_vars['BIC'] = xr.DataArray([BIC], dims=('fit'))
            data_vars['BIC'].attrs['units'] = ''

        data_vars['number_of_components'] = xr.DataArray([self.number_of_exponentials], dims='fit')

        coords = dict(component=np.arange(self.number_of_exponentials))
        dwell_analysis = xr.Dataset(data_vars=data_vars, coords=coords)
        dwell_analysis.attrs['papylio_version'] = papylio.__version__
        dwell_analysis.attrs['fit_function'] = 'exponential'
        dwell_analysis.attrs['truncation'] = self.truncation
        dwell_analysis.attrs['P_bounds'] = self.P_bounds
        dwell_analysis.attrs['k_bounds'] = self.k_bounds
        return dwell_analysis

    def dataset_to_parameters(self, dataset):
        # dataset = dataset.sel(fit=number_of_components == self.number_of_exponentials) # Could be used to assure the correct dataset is passed.
        parameters = dataset[['P', 'k']].to_array().values.flatten()
        parameters = parameters[~np.isnan(parameters)]
        parameters = list(parameters)
        parameters.pop(self.number_of_exponentials-1)
        return parameters

def auto_bin_size_for_discrete_data(dwell_times):
    dwell_times.sort()
    d = np.diff(dwell_times)
    bin_width_min = d[d > 0][1]

    Q1 = np.percentile(dwell_times, 25)
    Q3 = np.percentile(dwell_times, 75)
    IQR = Q3 - Q1
    bin_width = 2 * IQR / len(dwell_times) ** (1 / 3)

    bin_width = np.ceil(bin_width / bin_width_min) * bin_width_min

    # plot_range = (bin_width / 2, np.percentile(dwell_times, 99))
    # plot_range = (np.min(dwell_times) / 2, np.max(dwell_times))
    plot_range = (np.min(dwell_times) / 2, np.max(dwell_times))

    bin_edges = np.arange(plot_range[0], plot_range[1], bin_width)

    return bin_edges

def plot_dwell_time_histogram(dwell_times, bins='auto_discrete', range=None, ax=None, **hist_kwargs):
    if ax is None:
        fig, ax = plt.subplots()

    if bins == 'auto_discrete':
        bins = auto_bin_size_for_discrete_data(dwell_times)
        bins = bins[(bins>range[0])&(bins<range[1])]
    else:
        bins = np.histogram_bin_edges(dwell_times, bins=bins, range=range)

    bin_width = bins[1] - bins[0]
    weights = (1/len(dwell_times)/ bin_width,) * len(dwell_times) # Assuming evenly spaced bins
    counts, bin_edges, _ = ax.hist(dwell_times, bins=bins, range=range, weights=weights, density=False, **hist_kwargs)
    # relative_counts = counts/(len(dwell_times) * np.diff(bin_edges))
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    ax.set_ylim(0,counts.max()*1.03)
    return counts, bin_centers

def empirical_cdf(dwell_times, sampling_interval):
    # t = np.hstack([0, np.sort(dwell_times)])
    # empirical_cdf = np.arange(0, len(t)) / (len(t) - 1)

    t, counts  = np.unique(dwell_times, return_counts=True)
    t = t + sampling_interval / 2
    empirical_cdf = counts.cumsum() / counts.sum()
    # t = np.hstack([[0,], t])
    # empirical_cdf = np.hstack([[0, ], empirical_cdf])

    return t, empirical_cdf

def plot_empirical_cdf(dwell_times, sampling_interval, ax=None, **plot_kwargs):
    if ax is None:
        fig, ax = plt.subplots()

    t, ecdf = empirical_cdf(dwell_times, sampling_interval)

    ax.plot(t, ecdf, **plot_kwargs)

def fit_dwell_times(dwell_times, method='maximum_likelihood_estimation', number_of_exponentials=[1,2],
                    P_bounds=(-np.inf,np.inf), k_bounds=(0,np.inf), sampling_interval=None, fit_dwell_times_kwargs={}):
    if isinstance(number_of_exponentials, numbers.Number):
        number_of_exponentials = [number_of_exponentials]

    if sampling_interval is None:
        sampling_interval = np.min(dwell_times)

    dwell_analysis = []
    for i in number_of_exponentials:
        distribution = ExponentialDistribution(i, P_bounds, k_bounds, truncation=(sampling_interval/2, np.inf),
                                               sampling_interval=sampling_interval)
        # self.bounds[self.bounds[:, 1] > bounds_max, 1] = bounds_max
        dwell_analysis_i = getattr(distribution, method)(dwell_times, **fit_dwell_times_kwargs)
        # dwell_analysis.append(distribution.parameters_to_dataframe(optimal_parameters, BIC=True, dwell_times=dwell_times))
        dwell_analysis.append(dwell_analysis_i)

    # dwell_analysis[e] = {'parameters': optimal_parameters,
    #                   'parameter_errors': parameter_errors,
    #                   'BIC': distribution.BIC(dwell_times, *optimal_parameters) }

    # return pd.concat(dwell_analysis)
    dwell_analysis = xr.concat(dwell_analysis, dim='fit', combine_attrs="drop_conflicts")
    dwell_analysis.attrs['sampling_interval'] = sampling_interval
    return dwell_analysis

def plot_dwell_analysis_state(dwell_analysis, dwell_times, plot_type='pdf_binned', plot_range=None, bins='auto_discrete', log=False, ax=None):
    if ax is None:
        fig, ax = plt.subplots()

    if plot_range is None:
        plot_range = (0, np.max(dwell_times)) #np.quantile(dwell_times, 0.99))

    sampling_interval = dwell_analysis.sampling_interval

    if plot_type in ['pdf', 'pdf_binned']:
        counts, bin_centers = plot_dwell_time_histogram(dwell_times, bins=bins, range=plot_range, log=log, ax=ax)
        bin_width = bin_centers[1] - bin_centers[0]
    elif plot_type == 'cdf':
        plot_empirical_cdf(dwell_times, sampling_interval, ax=ax)
        ax.set_xlim(*plot_range)
        bin_width = None
        if log:
            ax.set_yscale('log')

    dwell_analysis_formatted = dwell_analysis[['P','k','BIC','number_of_components']].to_dataframe().dropna(subset=['P','k']) #dwell_analysis.copy()
    dwell_analysis_formatted['P'] = [f'{float(x):.2f}' if pd.notna(x) else '' for x in dwell_analysis_formatted['P']]
    dwell_analysis_formatted['k'] = [f'{float(x):.4f} /s' if pd.notna(x) else '' for x in dwell_analysis_formatted['k']]
    dwell_analysis_formatted['BIC'] = [f'{float(x):.0f}' if pd.notna(x) else '' for x in dwell_analysis_formatted['BIC']]
    dwell_analysis_formatted.loc[dwell_analysis_formatted.index.get_level_values('component') > 0, 'BIC'] = ''

    # dwell_analysis_formatted = dwell_analysis_formatted.reset_index()
    # exponentials = dwell_analysis_formatted['Exponential'].values
    # dwell_analysis_formatted.loc[
    #     np.hstack([np.where(exponentials == e)[0][1:] for e in np.unique(exponentials)]), 'Exponential'] = ''
    # dwell_analysis_formatted['Exponential'] = dwell_analysis_formatted['Exponential'].replace([1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    #                                                ['Single', 'Double', 'Triple', 'Quadruple', 'Quintuple',
    #                                                 'Sextuple', 'Septuple', 'Octuple', 'Nonuple',
    #                                                 'Decuple'])

    dwell_analysis_formatted = dwell_analysis_formatted.reset_index()
    number_of_components = dwell_analysis_formatted['number_of_components'].values
    dwell_analysis_formatted['Exponential'] = dwell_analysis_formatted['number_of_components'].replace([1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                                                   ['Single', 'Double', 'Triple', 'Quadruple', 'Quintuple',
                                                    'Sextuple', 'Septuple', 'Octuple', 'Nonuple',
                                                    'Decuple'])
    dwell_analysis_formatted.loc[dwell_analysis_formatted.component > 0, 'Exponential'] = ''
    dwell_analysis_formatted = dwell_analysis_formatted.set_index('number_of_components')[['Exponential','P','k','BIC']]

    labels = dwell_analysis_formatted.to_string(header=True, index=False, justify='left').split('\n')
    header = labels.pop(0)
    labels = np.array(labels)

    ax.plot([], [], color='k', linestyle=None, lw=0, marker=None, label=header)

    # for i, number_of_exponentials in enumerate(dwell_analysis.index.unique().values):
    #     distribution = ExponentialDistribution(number_of_exponentials)
    #
    #     label = '\n'.join(labels[dwell_analysis.index == number_of_exponentials])
    #
    #     ax.plot(t, distribution.pdf(t, *distribution.dataframe_to_parameters(dwell_analysis.loc[[number_of_exponentials]])),
    #             label=label)

    for i, number_of_exponentials in enumerate(np.unique(number_of_components)):
        distribution = ExponentialDistribution(number_of_exponentials, truncation=(sampling_interval/2, np.inf),
                                               sampling_interval=sampling_interval)
        distribution.bin_width = bin_width
        label = '\n'.join(labels[dwell_analysis_formatted.reset_index().number_of_components == number_of_exponentials])

        parameters = distribution.dataset_to_parameters(dwell_analysis.sel(fit=dwell_analysis.number_of_components==number_of_exponentials))

        if plot_type == 'pdf_binned':
            t = np.vstack([bin_centers - bin_width / 2, bin_centers + bin_width / 2]).T.flatten()
            y = getattr(distribution, plot_type)(bin_centers, *parameters)
            y = np.repeat(y, 2)
        else:
            t = np.linspace(plot_range[0], plot_range[1], 1000)
            y = getattr(distribution, plot_type)(t, *parameters)

        ax.plot(t, y, label=label)

    legend = ax.legend(prop={"family": "monospace"}, labelcolor='linecolor', frameon=False)
    ax.set_xlabel('Dwell time (s)')
    if plot_type == 'pdf':
        ax.set_ylabel('Normalized counts')
    elif plot_type == 'cdf':
        ax.set_ylabel('Normalized cumulative counts')

    return ax.figure, ax

def analyze_dwells(dwells, method='maximum_likelihood_estimation', number_of_exponentials=[1,2,3], state_names=None,
                   P_bounds=(-np.inf,np.inf), k_bounds=(1e-9,np.inf), sampling_interval=None, fit_dwell_times_kwargs={}):
    # number_of_exponentials can be given per state as {0: [1,2,3], 1: [1,2]}
    if state_names is None:
        states = np.unique(dwells.state)
        states = states[states >= 0]
        # state_names = {state: '' for state in states}
    else:
        states = np.array(list(state_names.keys()))

    if not isinstance(number_of_exponentials, dict):
        number_of_exponentials = {state: number_of_exponentials for state in states}

    # number_of_exponentials_max = np.concatenate(list(number_of_exponentials.values())).max()
    #
    # # fit_parameters = list(inspect.signature(fit_function).parameters)[1:]
    # fit_parameters = ['P', 'k']
    # fit_exponentials = np.arange(number_of_exponentials_max)+1
    # coords = {'state': states, 'fit': ,'exponential': fit_exponentials, 'parameter': fit_parameters}
    # dwell_analysis = xr.Dataset(coords=coords)
    # dwell_analysis['optimal_value'] = xr.DataArray(np.nan, dims=coords.keys(), coords=coords)
    # # dwell_analysis['error'] = xr.DataArray(np.nan, dims=('state', 'parameter'), coords={'state': positive_states, 'parameter': fit_parameters})
    # # dwell_analysis['covariance'] = xr.DataArray(np.nan, dims=('state', 'parameter','parameter'),
    # #                                    coords={'state': positive_states, 'parameter': fit_parameters, 'parameter': fit_parameters})
    #
    # # dwell_analysis.attrs['fit_function'] = fit_function.__name__
    # dwell_analysis['BIC'] = xr.DataArray(np.nan, dims=('state'), coords=coords)
    # dwell_analysis.attrs['version'] = papylio.__version__
    # dwell_analysis.attrs['method'] = method

    # dwell_analysis = {}
    # for i, state in enumerate(states):
    #     dwells_with_state = dwells.sel(dwell=dwells.state==state)
    #
    #     dwell_times = dwells_with_state.duration.values
    #     dwell_analysis_state = fit_dwell_times(dwell_times, method=method, number_of_exponentials=number_of_exponentials[state])
    #     dwell_analysis[(state, state_names[state])] = dwell_analysis_state
    #
    # dwell_analysis = pd.concat(dwell_analysis, names=('State', 'State name', 'Exponential'))
    #
    # return dwell_analysis


    dwell_analysis = []
    for i, state in enumerate(states):
        dwells_with_state = dwells.sel(dwell=dwells.state==state)

        dwell_times = dwells_with_state.duration.values
        dwell_analysis_state = fit_dwell_times(dwell_times, method=method, number_of_exponentials=number_of_exponentials[state],
                                               P_bounds=P_bounds, k_bounds=k_bounds, sampling_interval=sampling_interval,
                                               fit_dwell_times_kwargs=fit_dwell_times_kwargs)
        dwell_analysis.append(dwell_analysis_state.expand_dims({'state': [state]}))

    dwell_analysis = xr.concat(dwell_analysis, dim='state')
    if state_names is not None:
        dwell_analysis = dwell_analysis.assign_coords(state_name=('state', list(state_names.values())))

    return dwell_analysis


# [['P','k']].to_dataframe().dropna()
def plot_dwell_analysis(dwell_analysis, dwells, plot_type='pdf_binned', plot_range=None, axes=None, bins='auto_discrete',
                        log=False, sharey=True):

    # states = dwell_analysis.index.get_level_values('State').unique().values
    states = dwell_analysis.state.values

    if bins is None or isinstance(bins, numbers.Number) or isinstance(bins, str):
        bins = [bins] * len(states)

    if isinstance(plot_type, str):
        plot_type = [plot_type] * len(states)

    if plot_range is None or isinstance(plot_range[0], numbers.Number):
        plot_range = [plot_range] * len(states)

    if axes is None:
        fig, axes = plt.subplots(1,len(states), figsize=(len(states)*4.5, 4), layout='constrained', sharey=sharey)

    from collections.abc import Iterable
    if isinstance(axes, Iterable):
        axes = [axes]

    for i, (state, dwell_analysis_state) in enumerate(dwell_analysis.groupby('state')):

        dwell_times = dwells.sel(dwell=dwells.state==state).duration.values
        plot_dwell_analysis_state(dwell_analysis_state, dwell_times, plot_type=plot_type[i], plot_range=plot_range[i], bins=bins[i], log=log, ax=axes[i])
        if 'state_name' in dwell_analysis_state.data_vars:
            axes[i].set_title(dwell_analysis_state.state_name[0])
        else:
            axes[i].set_title(state)

        if i > 0:
            axes[i].set_ylabel('')

    return axes[0].figure, axes

#
# def plot_dwell_analysis(dwell_analysis, plot_range, axes=None, sharey=True):
#     states = dwell_analysis.index.get_level_values('State').unique().values
#     if axes is None:
#         fig, axes = plt.subplots(1,len(states), figsize=(len(states)*4, 4), layout='constrained', sharey=sharey)
#
#     for i, (state, fit_result_state) in enumerate(dwell_analysis.groupby('State')):
#         state_name = fit_result_state.index.get_level_values('State name')[0]
#         fit_result_state = fit_result_state.droplevel(level='State').droplevel(level='State name')
#         plot_dwell_analysis_state(fit_result_state, plot_range, ax=axes[i])
#         if i == 0:
#             axes[i].set_title(state_name)
#
#     return fig, axes



# if plot_range is None:
#     plot_range = (0, np.max(dwell_times))

 # plot_dwell_time_histogram(dwell_times, bins=bins, log=log, ax=ax)
