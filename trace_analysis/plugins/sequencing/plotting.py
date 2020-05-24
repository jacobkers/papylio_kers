import numpy as np
import matplotlib.pyplot as plt
from trace_analysis.coordinate_transformations import transform
import matplotlib.patches as patches
from matplotlib.patches import ConnectionPatch
from mpl_toolkits.axes_grid1 import make_axes_locatable

def plot_sequencing_match(match, write_path, name, unit = 'um', ax1pos = None, ax2pos = None):

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