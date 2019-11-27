import numpy as np #scientific computing with Python
import matplotlib.pyplot as plt
from trace_analysis.coordinate_transformations import transform
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable

def histogram(input, axis, makeFit=False):
    if not input: return None
    if not axis: axis = plt.gca()
    #    if not isinstance(input,list): input = [input]
    #
    #    molecules = list()
    #
    #    for i in input:
    #        if isinstance(i, Molecule):
    #            molecules.append(i)
    #        else:
    #            molecules.append(i.molecules)
    molecules = input

    # data = np.concatenate([molecule.intensity[0,:] for molecule in molecules])
    # axis.hist(data,100)
    data = np.concatenate([molecule.E() for molecule in molecules])
    axis.hist(data, 100, range=(0, 1))

    if makeFit:
        fit_hist(data, axis)


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

def scatter_coordinates(pointsets):
    for pointset in pointsets:
        plt.scatter(pointset[:,0], pointset[:,1])

def show_point_connections(pointset1,pointset2):
    for coordinate1, coordinate2 in zip(pointset1, pointset2):
        plt.plot([coordinate1[0],coordinate2[0]],[coordinate1[1],coordinate2[1]], color='r')


def plot_match(match, write_path, name, unit = 'um', ax1pos = None, ax2pos = None):

    def MiSeq_pixels_to_um(pixels):
        return 958 / 2800 * (pixels - 1000) / 10

    source = match.transform_source_to_destination
    destination = match.destination
    transformationMatrix = match.transformation

    # ps2 = 5.2 / 30 * ps2

    image_size_in_source = np.array([[1024,0],[2048,2048]])
    image_size_in_destination = transform(image_size_in_source, transformationMatrix)

    if unit == 'um':
        source = MiSeq_pixels_to_um(source)
        destination = MiSeq_pixels_to_um(destination)
        image_size_in_destination = MiSeq_pixels_to_um(image_size_in_destination)

    # fig = plt.figure(figsize = (8,4))
    # ax1, ax2 = fig.subplots(1, 2)

    fig, ax1 = plt.subplots(figsize = (8,4))
    fig.subplots_adjust(0.05,0.05,0.95,0.95)

    divider = make_axes_locatable(ax1)
    ax2 = divider.append_axes('right', size='60%', pad=0.5)
    #
    # from mpl_toolkits.axes_grid1 import ImageGrid
    #
    # fig = plt.figure(figsize = (8,4))
    # grid = ImageGrid(fig, 111, nrows_ncols=(1,2), axes_pad=0.1, add_all=True, label_mode='L')
    #
    # ax1 = grid[0]
    # ax2 = grid[1]

    ax1.scatter(source[:,0],source[:,1],c='#40A535',marker = 'x')
    ax1.scatter(destination[:,0],destination[:,1], marker = '.', facecolors = 'k', edgecolors='k')
    ax1.set_facecolor('white')
    ax1.set_aspect('equal')


    ax1.set_xlim([0, 30000])
    ax1.set_ylim([0, 30000])
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')

    ax1.ticklabel_format(axis='both', style='sci', scilimits=(5,6), useOffset=None, useLocale=None, useMathText=True)

    if unit == 'um':
        start = 0
        end = 1001
        stepsize = 250

        ax1.set_xlim([0, 1000])
        ax1.set_ylim([0, 1000])
        ax1.set_xlabel('x (\u03BCm)')
        ax1.set_ylabel('y (\u03BCm)')
        ax1.xaxis.set_ticks(np.arange(start, end, stepsize))
        ax1.yaxis.set_ticks(np.arange(start, end, stepsize))

        ax1.ticklabel_format(axis='both', style='sci', scilimits=(0,4), useOffset=None, useLocale=None, useMathText=True)


    p = patches.Rectangle(
        #(left, bottom), width, height,
        (image_size_in_destination[0,0],image_size_in_destination[0,1]),
        image_size_in_destination[1,0]-image_size_in_destination[0,0],
        image_size_in_destination[1,1]-image_size_in_destination[0,1],
        linewidth=2,edgecolor='#40A535',facecolor='none', fill='false'
        )

    ax1.add_patch(p)


    #plt.tight_layout()

    # fig.savefig(write_path.joinpath(name + '.pdf'), bbox_inches='tight')
    # fig.savefig(write_path.joinpath(name + '.png'), bbox_inches='tight', dpi=1000)


    ax2.scatter(source[:,0],source[:,1],c='#40A535',marker = 'x')
    ax2.scatter(destination[:,0],destination[:,1], marker = '.', facecolors = 'k', edgecolors='k')
    ax2.set_facecolor('white')
    ax2.set_aspect('equal')

    ax2.set_xlim([image_size_in_destination[1,0],image_size_in_destination[0,0]])
    ax2.set_ylim([image_size_in_destination[0,1], image_size_in_destination[1,1]])

    ax2.set_xlabel('x')
    ax2.set_ylabel('y')

    ax2.ticklabel_format(axis='both', style='sci', scilimits=(5,6), useOffset=None, useLocale=None, useMathText=True)

    if unit == 'um':
        ax2.set_xlim([image_size_in_destination[1,0],image_size_in_destination[0,0]])
        ax2.set_ylim([image_size_in_destination[0,1],image_size_in_destination[1,1]])
        ax2.set_xlabel('x (\u03BCm)')
        ax2.set_ylabel('y (\u03BCm)')

        ax2.ticklabel_format(axis='both', style='sci', scilimits=(0,4), useOffset=None, useLocale=None, useMathText=True)

    if ax1pos is not None:
        ax1.set_position(ax1pos)
        ax2.set_position(ax2pos)

    for spine in ax2.spines.values():
        spine.set_edgecolor('#40A535')

    #plt.tight_layout()

    # fig.savefig('zoomInScatter.pdf', bbox_inches='tight', transparent=True)
    # fig.savefig('zoomInScatter.png', bbox_inches='tight')



    from matplotlib.patches import ConnectionPatch
    con = ConnectionPatch(xyA=(image_size_in_destination[0,0],image_size_in_destination[0,1]),
                          xyB=(image_size_in_destination[0,0],image_size_in_destination[0,1]),
                          coordsA="data", coordsB="data",
                          axesA=ax1, axesB=ax2, color="k")
    ax1.add_artist(con)

    con = ConnectionPatch(xyA=(image_size_in_destination[0,0],image_size_in_destination[1,1]),
                          xyB=(image_size_in_destination[0,0],image_size_in_destination[1,1]),
                          coordsA="data", coordsB="data",
                          axesA=ax1, axesB=ax2, color="k")
    ax1.add_artist(con)

    con = ConnectionPatch(xyA=(image_size_in_destination[1,0],image_size_in_destination[0,1]),
                          xyB=(image_size_in_destination[1,0],image_size_in_destination[0,1]),
                          coordsA="data", coordsB="data",
                          axesA=ax1, axesB=ax2, color="k")
    ax1.add_artist(con)

    con = ConnectionPatch(xyA=(image_size_in_destination[1,0],image_size_in_destination[1,1]),
                          xyB=(image_size_in_destination[1,0],image_size_in_destination[1,1]),
                          coordsA="data", coordsB="data",
                          axesA=ax1, axesB=ax2, color="k")
    ax1.add_artist(con)

    for artist in ax1.artists:
        artist.set_linestyle((0,(5,5)))
        artist.set_linewidth(0.5)
        artist.set_edgecolor('grey')

    plt.show()

    n = name.replace('\\', '_')
    fig.savefig(write_path.joinpath(n + '.pdf'), bbox_inches='tight')
    fig.savefig(write_path.joinpath(n + '.png'), bbox_inches='tight', dpi=1000)

    return ax1, ax2