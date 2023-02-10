import numpy as np
import xarray as xr


def trace_classification_threshold(traces, threshold):
    classification_lower = classification_upper = True
    if threshold[0] is not None:
        classification_lower = np.vstack([(trace > threshold[0]) for trace in traces.values])
    if threshold[1] is not None:
        classification_upper = np.vstack([(trace < threshold[1]) for trace in traces.values])

    classification = xr.DataArray((classification_upper & classification_lower),
                                  dims=('molecule', 'frame'), name='classification')
    return classification


def trace_selection_threshold(traces, threshold):
    classification = trace_classification_threshold(traces, threshold)
    return classification.all(dim='frame')



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