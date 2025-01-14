import numpy as np #scientific computing with Python
import matplotlib.pyplot as plt

import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable

# from trace_analysis.molecule import Molecule

def histogram(da, axis=None, **hist_kwargs):
    if axis is None:
        figure, axis = plt.subplots()
    else:
        figure = axis.figure

    if 'channel' in da.dims:
        das = [da.sel(channel=channel) for channel in da.channel]
    else:
        das = [da]

    if len(das) > 1:
        hist_kwargs['histtype'] = 'step'

    for da in das:
        da.plot.hist(ax=axis, **hist_kwargs)
    axis.set_ylabel('Count')
    axis.set_title('')

    # if save:
    #     fig.savefig(self.absoluteFilePath.with_name(f'{self.name}_{parameter}_histogram').with_suffix('.png'))

    return figure, axis



def histogram_FRET(data, axis, **kwargs):
    axis.hist(data, range=(0, 1), **kwargs)
    axis.set_xlim((0, 1))
    axis.set_xlabel('FRET')
    axis.set_ylabel('Count')

def fit_hist(data, axis):
    hist, bin_edges = np.histogram(data, 100, range=(0, 1))
    bin_centers = (bin_edges[0:-1] + bin_edges[1:]) / 2

    # plt.plot(bin_centers,hist)

    from scipy.signal import butter
    from scipy.signal import filtfilt
    b, a = butter(2, 0.2, 'low')
    output_signal = filtfilt(b, a, hist)
    plt.plot(bin_centers, output_signal)

    from scipy.signal import find_peaks
    peaks, properties = find_peaks(output_signal, prominence=5, width=7)  # prominence=1
    plt.plot(bin_centers[peaks], hist[peaks], "x")

    def func(x, a, b, c, d, e, f):
        return a * np.exp(-(x - b) ** 2 / (2 * c ** 2)) + d * np.exp(-(x - e) ** 2 / (2 * f ** 2))

    from scipy.optimize import curve_fit
    popt, pcov = curve_fit(func, bin_centers, hist, method='trf',
                           p0=[hist[peaks[0]], bin_centers[peaks[0]], 0.1, hist[peaks[1]], bin_centers[peaks[1]], 0.1],
                           bounds=(0, [np.inf, 1, 1, np.inf, 1, 1]))

    axis.plot(bin_centers, func(bin_centers, *popt))
    # plt.plot(bin_centers,func(bin_centers, 10000,0.18,0.1,5000,0.5,0.2))

# uniqueFileNames = list(set([re.search('hel[0-9]*',fileName).group() for fileName in fileNames]))


def show_image_3d(image, figure=None):
    if not figure:
        figure = plt.figure()

    from matplotlib import cm
    axis = figure.gca(projection='3d')
    X = np.arange(image.shape[1])
    Y = np.arange(image.shape[0])
    X, Y = np.meshgrid(X, Y)
    axis.plot_surface(X, Y, image, cmap=cm.coolwarm,
                      linewidth=0, antialiased=False)

