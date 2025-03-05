import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import pandas as pd

class ExponentialDistribution:
    def __init__(self, number_of_exponentials):
        self.number_of_exponentials = number_of_exponentials
        self.bounds = np.array([(0, np.inf)] * (2 * self.number_of_exponentials - 1))

    def __call__(self, t, *parameters):
        return self.pdf(t, *parameters)

    @property
    def parameter_names(self):
        return ([f'P{i}' for i in range(self.number_of_exponentials-1)] +
                [f'k{i}' for i in range(self.number_of_exponentials-1)])

    def pdf(self, t, *parameters):
        # Parameters given as P1, P2, k1, k2, k3
        result = np.zeros_like(t).astype(float)
        P = parameters[0:(self.number_of_exponentials-1)]
        P += (1-np.sum(P),)
        k = parameters[(self.number_of_exponentials-1):]
        for i in range(self.number_of_exponentials):
            result += P[i] * k[i] * np.exp(-k[i] * t)
        return result

    def cdf(self, t, *parameters):
        result = np.ones_like(t).astype(float)
        P = parameters[0:(self.number_of_exponentials-1)]
        P += (1-np.sum(P),)
        k = parameters[(self.number_of_exponentials-1):]
        for i in range(self.number_of_exponentials):
            result -= P[i] * np.exp(-k[i] * t)
        return result

    def likelihood(self, parameters, t):
        return np.prod(self.pdf(t, *parameters))

    def loglikelihood(self, parameters, t):
        return np.sum(np.log(self.pdf(t, *parameters)))

    def negative_loglikelihood(self, parameters, t):
        return -self.loglikelihood(parameters, t)

    def maximum_likelihood_estimation(self, t):
        optimal_parameters = scipy.optimize.minimize(fun=self.negative_loglikelihood, x0=self.parameter_guess(t),
                                                     args=t, bounds=self.bounds).x
        return optimal_parameters, None

    def histogram_fitting(self, t, bins='auto'):
        counts, bin_edges = np.histogram(t, bins=bins, density=True)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        optimal_parameters, parameter_covariances = scipy.optimize.curve_fit(self.pdf, bin_centers, counts,
                                                                            p0=self.parameter_guess(t),
                                                                            bounds=self.bounds.T, absolute_sigma=True)
        parameter_errors = np.sqrt(np.diag(parameter_covariances))
        return optimal_parameters, parameter_errors

    def cdf_fitting(self, t):
        t = np.sort(t)
        empirical_cdf = np.arange(1, len(t)+1) / len(t)
        optimal_parameters, parameter_covariances = scipy.optimize.curve_fit(self.cdf, t, empirical_cdf,
                                                                            p0=self.parameter_guess(t),
                                                                            bounds=self.bounds.T, absolute_sigma=True)
        parameter_errors = np.sqrt(np.diag(parameter_covariances))
        return optimal_parameters, parameter_errors

    def parameter_guess(self, t):
        guess_P = [1 / (i + 1) for i in range(1, self.number_of_exponentials)]
        guess_k = [1 / (t.mean() * (i + 1)) for i in range(self.number_of_exponentials)]
        parameters = guess_P + guess_k
        return parameters

    def BIC(self, t, *optimal_parameters):
        return len(optimal_parameters) * np.log(len(t)) - 2 * self.loglikelihood(optimal_parameters, t)

    def parameters_to_dataframe(self, parameters, bic=True, dwell_times=None):
        fit_results = {}
        fit_results['Exponential'] = [self.number_of_exponentials] * self.number_of_exponentials #+ [''] * (e-1)
        fit_results['P'] = np.hstack([parameters[0:self.number_of_exponentials-1] , np.nan])
        fit_results['P'][-1] = 1-np.nansum(fit_results['P'])
        fit_results['k'] = parameters[self.number_of_exponentials - 1:]
        # if parameter_errors is not None:
        #     fit_results['P error'] = np.hstack([parameter_errors[0:e-1] ,np.nan])
        #     fit_results['k error'] = parameter_errors[e - 1:]
        if bic and dwell_times is not None:
            fit_results['BIC']  = [self.BIC(dwell_times, *parameters)] + [np.nan] * (self.number_of_exponentials-1)

        return pd.DataFrame(fit_results).set_index('Exponential')

    def dataframe_to_parameters(self, dataframe):
        parameters = dataframe.loc[self.number_of_exponentials][['P', 'k']].values.T.flatten()
        # return parameters[~np.isnan(parameters)]
        parameters = list(parameters)
        parameters.pop(self.number_of_exponentials-1)
        return parameters


def plot_dwell_time_histogram(dwelltimes, bins='auto', log=False, ax=None):
    if ax is None:
        fig, ax = plt.subplots()
    counts, bin_edges, _ = ax.hist(dwelltimes, bins=bins, log=log, density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    return counts, bin_centers


def fit_dwell_times(dwell_times, method='maximum_likelihood_estimation', number_of_exponentials=[1,2]):
    if isinstance(number_of_exponentials, np.number):
        number_of_exponentials = [number_of_exponentials]

    fit_results = []
    for e in number_of_exponentials:
        distribution = ExponentialDistribution(e)
        optimal_parameters, parameter_errors = getattr(distribution, method)(dwell_times)
        fit_results.append(distribution.parameters_to_dataframe(optimal_parameters, bic=True, dwell_times=dwell_times))

    # fit_results[e] = {'parameters': optimal_parameters,
    #                   'parameter_errors': parameter_errors,
    #                   'BIC': distribution.BIC(dwell_times, *optimal_parameters) }

    return pd.concat(fit_results)

def plot_fit_results(dwell_times, fit_results, bins='auto', plot_range=None, log=False, ax=None):
    if plot_range is None:
        plot_range = (0, np.max(dwell_times))

    if ax is None:
        fig, ax = plt.subplots()

    # if semi_log:
    #     ax.set_yscale('log')

    t = np.linspace(plot_range[0], plot_range[1], 1000)

    plot_dwell_time_histogram(dwell_times, bins=bins, log=log, ax=ax)


    fit_results_formatted = fit_results.copy()
    fit_results_formatted['P'] = [f'{float(x):.2f}' if pd.notna(x) else '' for x in fit_results_formatted['P']]
    fit_results_formatted['k'] = [f'{float(x):.4f} /s' if pd.notna(x) else '' for x in fit_results_formatted['k']]
    fit_results_formatted['BIC'] = [f'{float(x):.0f}' if pd.notna(x) else '' for x in fit_results_formatted['BIC']]


    fit_results_formatted = fit_results_formatted.reset_index()
    exponentials = fit_results_formatted['Exponential'].values
    fit_results_formatted.loc[
        np.hstack([np.where(exponentials == e)[0][1:] for e in np.unique(exponentials)]), 'Exponential'] = ''
    fit_results_formatted['Exponential'] = fit_results_formatted['Exponential'].replace([1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                                                   ['Single', 'Double', 'Triple', 'Quadruple', 'Quintuple',
                                                    'Sextuple', 'Septuple', 'Octuple', 'Nonuple',
                                                    'Decuple'])

    labels = fit_results_formatted.to_string(header=True, index=False, justify='left').split('\n')
    header = labels.pop(0)
    labels = np.array(labels)

    ax.plot([], [], color='k', linestyle=None, lw=0, marker=None, label=header)

    for i, number_of_exponentials in enumerate(fit_results.index.unique().values):
        distribution = ExponentialDistribution(number_of_exponentials)

        label = '\n'.join(labels[fit_results.index == number_of_exponentials])

        ax.plot(t, distribution.pdf(t, *distribution.dataframe_to_parameters(fit_results.loc[[number_of_exponentials]])),
                label=label)

    legend = ax.legend(prop={"family": "monospace"}, labelcolor='linecolor', frameon=False)
    ax.set_xlabel('Dwell time (s)')
    ax.set_ylabel('Normalized counts')

def analyze_dwells(dwells, method='maximum_likelihood_estimation', number_of_exponentials=[1,2,3], state_names=None,
                   plot=False, axes=None, log=False, sharey=True):
    if state_names is None:
        states = np.unique(dwells.state)
        states = states[states >= 0]
        state_names = {state: state for state in states}
    else:
        states = np.array(list(state_names.keys()))

    if plot and axes is None:
        fig, axes = plt.subplots(1,len(states), figsize=(len(states)*4, 4), layout='constrained', sharey=sharey)

    # fit_parameters = list(inspect.signature(fit_function).parameters)[1:]
    # fit_values = xr.Dataset(coords={'state': positive_states, 'parameter': fit_parameters})
    # fit_values['optimal_value'] = xr.DataArray(np.nan, dims=('state', 'parameter'), coords={'state': positive_states, 'parameter': fit_parameters})
    # fit_values['error'] = xr.DataArray(np.nan, dims=('state', 'parameter'), coords={'state': positive_states, 'parameter': fit_parameters})
    # fit_values['covariance'] = xr.DataArray(np.nan, dims=('state', 'parameter','parameter'),
    #                                    coords={'state': positive_states, 'parameter': fit_parameters, 'parameter': fit_parameters})

    # fit_values.attrs['fit_function'] = fit_function.__name__

    for i, state in enumerate(states):
        dwells_with_state = dwells.sel(dwell=dwells.state==state)

        dwell_times = dwells_with_state.duration.values
        fit_results = fit_dwell_times(dwell_times, method=method, number_of_exponentials=number_of_exponentials)

        if plot:
            plot_fit_results(dwell_times, fit_results, log=log, ax=axes[i])
            axes[i].set_title(state_names[state])

    return fit_results, axes