import xarray as xr
import numpy as np
from itertools import accumulate, groupby
from hmmlearn import hmm
from tqdm import tqdm
import trace_analysis as ta


# exp = ta.Experiment(r'D:\20200918 - Test data\Single-molecule data small')
# exp = ta.Experiment(r'P:\SURFdrive\Promotie\Data\Test data')

# ds = xr.open_mfdataset(file_paths, concat_dim='molecule', combine='nested')
# classification_background = ~(ds.intensity < 5000).any('channel')

# file = files_green_laser[0]
# for file in files:
#
#     traces = file.dataset.FRET #ds.intensity.sel(channel=0, drop=True)

# n_components=2, covariance_type="full", n_iter=100

def hmm_traces(traces, initial_boolean_classification=None, **kwargs):
    if initial_boolean_classification is None:
        initial_boolean_classification = [None]*len(traces)

    classification = -xr.ones_like(traces).astype(int)
    classification.name = 'classification'
    hmm_transition_matrix = xr.DataArray(np.ones((len(traces.molecule), 2, 2))*np.nan, coords={},
                                         dims=['molecule', 'from_state', 'to_state'], name='hmm_transition_matrix')
    hmm_means = xr.DataArray(np.ones((len(traces.molecule), 2))*np.nan, coords={}, dims=['molecule','state'], name='hmm_means')
    hmm_log_probs = xr.DataArray(np.ones(len(traces.molecule))*np.nan, coords={}, dims='molecule', name='hmm_log_probabilities')

    for molecule in tqdm(traces.molecule.values):
        hmm_transition_matrix[molecule],  hmm_means[molecule], classification[molecule], hmm_log_probs[molecule] = \
            hmm_trace(traces[molecule], initial_boolean_classification=initial_boolean_classification[molecule], **kwargs)

    return xr.merge([hmm_transition_matrix, hmm_means, hmm_log_probs, classification])


    # hmm_transition_matrix.name = 'hmm_transition_matrix'
    # hmm_transition_matrix.to_netcdf(file.relativeFilePath.with_suffix('.nc'), engine='h5netcdf', mode='a')
    # hmm_means.name = 'hmm_means'
    # hmm_means.to_netcdf(file.relativeFilePath.with_suffix('.nc'), engine='h5netcdf', mode='a')
    # classification.name = 'classification'
    # classification.to_netcdf(file.relativeFilePath.with_suffix('.nc'), engine='h5netcdf', mode='a')

    # nds = xr.Dataset({'hmm_transition_matrix': hmm_transition_matrix, 'hmm_means': hmm_means, 'classification': classification})
    # nds.load()
    # nds.to_netcdf(file.relativeFilePath.with_suffix('.nc'), engine='h5netcdf', mode='a')

def hmm_traces(traces, **kwargs):
    hmm_transition_matrix, hmm_means, classification, hmm_log_probs =\
        xr.apply_ufunc(
            hmm_trace,
            traces,
            kwargs=kwargs,
            input_core_dims=[['frame']],
            output_core_dims=[['from_state','to_state'],['state'],['frame'],[]],
            # dask="parallelized",
            output_dtypes=[np.ndarray, np.ndarray, np.ndarray, float],
            dask_gufunc_kwargs=dict(output_sizes={"from_state": 2, "to_state": 2, 'state': 2}),
            vectorize=True
        )
    hmm_transition_matrix.name = 'hmm_transition_matrix'
    hmm_means.name = 'hmm_means'
    classification.name = 'classification'
    hmm_log_probs.name = 'hmm_log_probabilities'
    return xr.merge([hmm_transition_matrix, hmm_means, hmm_log_probs, classification])


def hmm_trace(trace, initial_boolean_classification=None, **kwargs):
    # TODO: Test handling of initial_boolean_classification
    # print(trace)
    # FRET = np.atleast_2d(traces.sel(molecule=molecule_index)).T
    trace_ = np.array(trace).reshape(-1, 1)
    # print(trace_)
    if initial_boolean_classification is not None:
        # initial_boolean_classification = np.ones_like(trace).astype(bool)
        trace_ = trace_[initial_boolean_classification]
        segment_lengths = [sum(1 for _ in group) for value, group in groupby(initial_boolean_classification) if value]
    else:
        segment_lengths = None
    model = hmm.GaussianHMM(**kwargs)
    model.fit(trace_, segment_lengths)
    try:
        log_prob, classification_HMM = model.decode(trace_)
        #classification_HMM = xr.DataArray(classification_HMM, coords=trace.coords)
        if initial_boolean_classification is not None:
            classification_HMM[~initial_boolean_classification] = -1
    except ValueError:
        print('!!!')
        log_prob = np.nan
        classification_HMM = -np.ones_like(trace).astype(int)

    return model.transmat_, model.means_.squeeze(), classification_HMM, log_prob



# p_stay = np.diagonal(hmm_transition_matrix, axis1=1, axis2=2)
#
# dt = 0.1 #s
# tau = -dt/np.log(p_stay)


# for file_path in file_paths:
#     with xr.open_dataset(file_path) as ds:
#         # ds = ds.reset_index('molecule', drop=True)
#         # # ds = ds.reset_index('dim_1', drop=True)
#         # ds = ds.reset_index('frame', drop=True)
#         # ds = ds.reset_index('from_state', drop=True)
#         # ds = ds.reset_index('to_state', drop=True)
#         ds.load()
#     ds['channel'] = [0,1]
#
#     ds.to_netcdf(file_path, engine='h5netcdf', mode='w')


