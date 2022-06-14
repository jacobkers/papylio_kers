


# def pre_classification(ds):
#     classification_da = ds['classification'].copy()
#     classification_da.name = 'pre_classification'
#     for molecule_index in ds.molecule:
#         molecule = ds.sel(molecule=molecule_index)
#         test = molecule.intensity.to_pandas()
#         corr = test.iloc[0, :].rolling(10, center=True, min_periods=1).corr(test.iloc[1, :])
#         # corr[0:9] = corr[9]
#         # corr2 = moving_average(corr, 15)
#         corr2 = corr.rolling(15, center=True, min_periods=1).mean()
#         selection = (corr2 < -0.60).rolling(15, center=True, min_periods=1).max().astype(int)
#         # corr2 = scipy.signal.medfilt(corr, 15)
#         selection2 = ds['intensity_total'].sel(molecule=molecule_index) < 55000
#
#         classification_da[dict(molecule=molecule_index)] = selection*selection2-1
#     return classification_da



def rolling_correlation(traces, rolling_dim='frame', correlation_dim='channel', window=10):
    windows = traces.rolling(dim={rolling_dim: window}, center=True, min_periods=1).construct(window_dim='section', stride=1, keep_attrs=None)

    mean_windows = windows.mean('section')
    windows_minus_mean = windows-mean_windows

    a = windows_minus_mean.prod(correlation_dim, skipna=False).sum('section')
    b = (windows_minus_mean**2).sum('section').prod(correlation_dim)**(1/2)
    p = a/b

    return p

def classify_correlation(traces, rolling_dim='frame', correlation_dim='channel', window=10, rolling_mean_window=10, threshold=0.75):
    rc = rolling_correlation(traces, rolling_dim=rolling_dim, correlation_dim=correlation_dim, window=window)
    rcm = rc.rolling(dim={rolling_dim: rolling_mean_window}, center=True, min_periods=1).mean()
    classification = (rcm > threshold).astype(int).rolling(dim={rolling_dim: rolling_mean_window}, center=True, min_periods=1).max()
    classification.name = 'classification'
    return classification


def classify_anticorrelation(traces, rolling_dim='frame', correlation_dim='channel', window=10, rolling_mean_window=10, threshold=-0.75):
    rc = rolling_correlation(traces, rolling_dim=rolling_dim, correlation_dim=correlation_dim, window=window)
    rcm = rc.rolling(dim={rolling_dim: rolling_mean_window}, center=True, min_periods=1).mean() # To smooth out variations
    classification = (rcm < threshold).astype(int).rolling(dim={rolling_dim: rolling_mean_window}, center=True, min_periods=1).max() # To widen the window
    classification.name = 'classification'
    return classification