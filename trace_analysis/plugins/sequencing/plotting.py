import numpy as np
import matplotlib.pyplot as plt
from trace_analysis.coordinate_transformations import transform
import matplotlib.patches as patches
from matplotlib.patches import ConnectionPatch
from mpl_toolkits.axes_grid1 import make_axes_locatable

def plot_sequencing_match(match, write_path, title, filename, unit = 'um', MiSeq_pixels_to_um = None, Fluo_pixels_to_um = None, save=True):
    source = match.source
    source_in_destination = match.transform_coordinates(source)
    destination = match.destination
    destination_in_source = match.transform_coordinates(destination, inverse=True)

    source_vertices = match.source_vertices
    destination_vertices = match.transform_coordinates(match.source_vertices)

    if unit == 'um':
        source = Fluo_pixels_to_um(source)
        source_in_destination = MiSeq_pixels_to_um(source_in_destination)
        destination = MiSeq_pixels_to_um(destination)
        destination_in_source = Fluo_pixels_to_um(destination_in_source)

        source_vertices = Fluo_pixels_to_um(source_vertices)
        destination_vertices = MiSeq_pixels_to_um(destination_vertices)

    # fig = plt.figure(figsize = (8,4))
    # ax1, ax2 = fig.subplots(1, 2)

    fig, ax1 = plt.subplots(figsize = (8,4))
    fig.subplots_adjust(0.05,0.05,0.95,0.93)

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

    ax1.scatter(source_in_destination[:,0],source_in_destination[:,1],c='#40A535',marker = 'x')
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

    p = patches.Polygon(destination_vertices,
                        linewidth=2, edgecolor='#40A535', facecolor='none', fill='false'
                        )

    ax1.add_patch(p)


    #plt.tight_layout()

    # fig.savefig(write_path.joinpath(name + '.pdf'), bbox_inches='tight')
    # fig.savefig(write_path.joinpath(name + '.png'), bbox_inches='tight', dpi=1000)


    ax2.scatter(source[:,0],source[:,1],c='#40A535',marker = 'x')
    ax2.scatter(destination_in_source[:,0],destination_in_source[:,1], marker = '.', facecolors = 'k', edgecolors='k')
    ax2.set_facecolor('white')
    ax2.set_aspect('equal')

    image_size = np.array([np.min(source_vertices, axis=0),np.max(source_vertices, axis=0)])

    ax2.set_xlim([image_size[0,0], image_size[1,0]])
    ax2.set_ylim([image_size[0,1], image_size[1,1]])

    ax2.set_xlabel('x')
    ax2.set_ylabel('y')

    ax2.ticklabel_format(axis='both', style='sci', scilimits=(5,6), useOffset=None, useLocale=None, useMathText=True)

    if unit == 'um':
        # ax2.set_xlim([image_size[0, 0], image_size[1, 0]])
        # ax2.set_ylim([image_size[0, 1], image_size[1, 1]])
        ax2.set_xlabel('x (\u03BCm)')
        ax2.set_ylabel('y (\u03BCm)')

        ax2.ticklabel_format(axis='both', style='sci', scilimits=(0,4), useOffset=None, useLocale=None, useMathText=True)

    # if ax1pos is not None:
    #     ax1.set_position(ax1pos)
    #     ax2.set_position(ax2pos)

    for spine in ax2.spines.values():
        spine.set_edgecolor('#40A535')

    def connect_vertices_in_axis(verticesA, verticesB, axisA, axisB, **kwargs):
        for vertexA, vertexB in zip(verticesA, verticesB):
            con = ConnectionPatch(xyA=vertexA, xyB=vertexB, coordsA="data", coordsB="data",
                                  axesA=axisA, axesB=axisB, **kwargs)
            axisB.add_artist(con)

    connect_vertices_in_axis(source_vertices, destination_vertices, ax2, ax1)

    for artist in ax1.artists:
        artist.set_linestyle((0,(5,5)))
        artist.set_linewidth(0.5)
        artist.set_edgecolor('grey')

    fig.suptitle(title, fontsize='medium')

    plt.show()

    if save:
        n = filename.replace('\\', '_')
        fig.savefig(write_path.joinpath(n + '.pdf'), bbox_inches='tight')
        fig.savefig(write_path.joinpath(n + '.png'), bbox_inches='tight', dpi=1000)

    return ax1, ax2