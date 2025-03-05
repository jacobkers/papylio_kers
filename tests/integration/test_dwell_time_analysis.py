import pytest
import tifffile
import numpy as np


from papylio.analysis.dwelltime_analysis import ExponentialDistribution, fit_dwell_times, plot_fit_results, analyze_dwells
from tests.integration.test_dwell_time_extraction import dwells

@pytest.fixture
def dwell_times_single_exponential():
    return np.random.exponential(1, 1000)

@pytest.fixture
def dwell_times_double_exponential():
    return np.hstack([np.random.exponential(1 / 0.3, 1000), np.random.exponential(1 / 0.01, 1000)])

def test_exponential_distribution():
    for i in range(3):
        ExponentialDistribution(i)

def test_maximum_likelihood_estimation(dwell_times_double_exponential):
    optimal_parameters = ExponentialDistribution(2).maximum_likelihood_estimation(dwell_times_double_exponential)

def test_histogram_fitting(dwell_times_double_exponential):
    optimal_parameters, _ = ExponentialDistribution(2).histogram_fitting(dwell_times_double_exponential)

def test_cdf_fitting(dwell_times_double_exponential):
    optimal_parameters, _ = ExponentialDistribution(2).cdf_fitting(dwell_times_double_exponential)

def test_fit_dwell_times(dwell_times_double_exponential):
    fit_results = fit_dwell_times(dwell_times_double_exponential, 'maximum_likelihood_estimation', number_of_exponentials=[1,2,3])

def test_plot_fit_results(dwell_times_double_exponential):
    fit_results = fit_dwell_times(dwell_times_double_exponential, 'maximum_likelihood_estimation',
                                  number_of_exponentials=[1,2,3])
    plot_fit_results(dwell_times_double_exponential, fit_results, log=True)

def test_analyze_dwells(dwells):
    fit_results, axis = analyze_dwells(dwells, plot=True, axes=None, state_names=None, log=False, sharey=False)