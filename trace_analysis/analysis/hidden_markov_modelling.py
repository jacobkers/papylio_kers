import xarray as xr
import numpy as np
from itertools import accumulate, groupby
from hmmlearn import hmm
import trace_analysis as ta


exp = ta.Experiment(r'D:\20200918 - Test data\Single-molecule data small')
# exp = ta.Experiment(r'P:\SURFdrive\Promotie\Data\Test data')

file_paths = [p for p in exp.nc_file_paths if '561' in str(p)]
# ds = xr.open_mfdataset(file_paths, concat_dim='molecule', combine='nested')
files = [file for file in exp.files if '561 ' in str(file.relativeFilePath)]

# file=exp.files[0]
for file in files:
    with xr.open_dataset(file.absoluteFilePath.with_suffix('.nc')) as ds:
        traces = ds.FRET #ds.intensity.sel(channel=0, drop=True)
        traces_sel = traces.sel(molecule=ds.sequence_name.str.contains('HJ'))
        traces_sel.load()
        # ds_sel = ds.sel(molecule=ds.sequence_name == 'HJ7_G')

        # classification_background = ~(ds.intensity < 5000).any('channel')


        classification_background = xr.ones_like(traces).astype(bool)

        classification = -xr.ones_like(traces).astype(int)
        hmm_transition_matrix = xr.DataArray(np.ones((len(traces.molecule), 2, 2))*np.nan, coords={}, dims=['molecule', 'from_state', 'to_state'])
        hmm_means = xr.DataArray(np.ones((len(traces.molecule), 2))*np.nan, coords={}, dims=['molecule','state'])
        logprobs = xr.DataArray(np.ones(len(traces.molecule))*np.nan, coords={}, dims='molecule')

        for molecule_index in traces_sel.molecule_in_file:
            cbg = classification_background.sel(molecule=molecule_index).values

            FRET = np.atleast_2d(traces.sel(molecule=molecule_index)).T
            segment_lengths = [sum(1 for _ in g) for v, g in groupby(cbg) if v]

            model = hmm.GaussianHMM(n_components=2, covariance_type="full", n_iter=100)
            model.fit(FRET[cbg], segment_lengths)
            logprob, classification_HMM = model.decode(FRET)
            classification_HMM[~cbg] = -1

            hmm_transition_matrix[traces.molecule == molecule_index] = model.transmat_
            hmm_means[traces.molecule == molecule_index] = model.means_.T
            classification[traces.molecule == molecule_index] = classification_HMM

        # hmm_transition_matrix.name = 'hmm_transition_matrix'
        # hmm_transition_matrix.to_netcdf(file.relativeFilePath.with_suffix('.nc'), engine='h5netcdf', mode='a')
        # hmm_means.name = 'hmm_means'
        # hmm_means.to_netcdf(file.relativeFilePath.with_suffix('.nc'), engine='h5netcdf', mode='a')
        # classification.name = 'classification'
        # classification.to_netcdf(file.relativeFilePath.with_suffix('.nc'), engine='h5netcdf', mode='a')

        nds = xr.Dataset({'hmm_transition_matrix': hmm_transition_matrix, 'hmm_means': hmm_means, 'classification': classification})
        nds.load()
    nds.to_netcdf(file.relativeFilePath.with_suffix('.nc'), engine='h5netcdf', mode='a')

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


