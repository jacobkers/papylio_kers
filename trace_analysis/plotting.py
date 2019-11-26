import numpy as np #scientific computing with Python
import matplotlib.pyplot as plt
from trace_analysis.coordinate_transformations import transform
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection

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


def plot_match(ps1, ps3, transformationMatrix, write_path, name, ax1pos = None, ax2pos = None):
    ps1 = 958 / 2800 * (ps1 - 1000) / 10
    # ps2 = 5.2 / 30 * ps2
    ps3 = 958 / 2800 * (ps3 - 1000) / 10

    #%%
    imageSize2 = np.array([[1024,0],[2048,2048]])
    imageSize2r = transform(imageSize2, reflection=0) # This reflects with repect to axis 0.
    # transformationMatrix = np.array(
    #         [[ 1.17434823e-02,  5.15661570e+00,  2.30346061e+04],
    #          [-5.15661570e+00,  1.17434823e-02,  1.45601679e+04],
    #          [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])

    imageSize3 = transform(imageSize2r, transformationMatrix)
    imageSize3 = 958/2800*(imageSize3-1000)/10


    fig = plt.figure(figsize = (10,4))
    ax1, ax2 = fig.subplots(1, 2)



    ax1.scatter(ps3[:,0],ps3[:,1],c='#40A535',marker = 'x')
    ax1.scatter(ps1[:,0],ps1[:,1], marker = '.', facecolors = 'k', edgecolors='k')

    ax1.set_facecolor('white')

    ax1.set_aspect('equal')
    #ax.set_title('Mapped')
    # ax.set_xlim([0, 30000])
    # ax.set_ylim([0, 30000])
    #
    # ax.set_xlabel('x')
    # ax.set_ylabel('y')
    #
    # ax.ticklabel_format(axis='both', style='sci', scilimits=(5,6), useOffset=None, useLocale=None, useMathText=True)

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
        (imageSize3[0,0],imageSize3[0,1]),imageSize3[1,0]-imageSize3[0,0],imageSize3[1,1]-imageSize3[0,1],
        linewidth=2,edgecolor='#40A535',facecolor='none', fill='false'
        )

    ax1.add_patch(p)


    #plt.tight_layout()

    # fig.savefig(write_path.joinpath(name + '.pdf'), bbox_inches='tight')
    # fig.savefig(write_path.joinpath(name + '.png'), bbox_inches='tight', dpi=1000)


    ax2.scatter(ps3[:,0],ps3[:,1],c='#40A535',marker = 'x')
    ax2.scatter(ps1[:,0],ps1[:,1], marker = '.', facecolors = 'k', edgecolors='k')

    ax2.set_facecolor('white')

    ax2.set_aspect('equal')
    #ax.set_title('Mapped')
    # ax.set_xlim([imageSize3[0,0], imageSize3[1,0]])
    # ax.set_ylim([imageSize3[0,1], imageSize3[1,1]])
    #
    # ax.set_xlabel('x')
    # ax.set_ylabel('y')
    #
    # ax.ticklabel_format(axis='both', style='sci', scilimits=(5,6), useOffset=None, useLocale=None, useMathText=True)

    ax2.set_xlim([imageSize3[1,0],imageSize3[0,0]])
    ax2.set_ylim([imageSize3[0,1],imageSize3[1,1]])
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
    con = ConnectionPatch(xyA=(imageSize3[0,0],imageSize3[0,1]), xyB=(imageSize3[0,0],imageSize3[0,1]),
                          coordsA="data", coordsB="data",
                          axesA=ax1, axesB=ax2, color="k")
    ax1.add_artist(con)

    con = ConnectionPatch(xyA=(imageSize3[0,0],imageSize3[1,1]), xyB=(imageSize3[0,0],imageSize3[1,1]),
                          coordsA="data", coordsB="data",
                          axesA=ax1, axesB=ax2, color="k")
    ax1.add_artist(con)

    con = ConnectionPatch(xyA=(imageSize3[1,0],imageSize3[0,1]), xyB=(imageSize3[1,0],imageSize3[0,1]),
                          coordsA="data", coordsB="data",
                          axesA=ax1, axesB=ax2, color="k")
    ax1.add_artist(con)

    con = ConnectionPatch(xyA=(imageSize3[1,0],imageSize3[1,1]), xyB=(imageSize3[1,0],imageSize3[1,1]),
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