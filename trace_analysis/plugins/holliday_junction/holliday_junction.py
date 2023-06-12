import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import numpy as np
import itertools
import tqdm
from matplotlib import patches


def basepair_count_per_position(basepair_count, save_path):
    title = (basepair_count.name.replace('basepair_count','basepair_count_per_position'))

    basepair_count_stacked = basepair_count.stack(basepair=('base_0', 'base_1'))
    basepair_count_stacked['basepair'] = xr.concat([basepair_count_stacked.base_0, basepair_count_stacked.base_1],
                                                   dim='dummy').str.join(dim='dummy').values

    basepaired_position = xr.concat(
        [basepair_count_stacked.sel(position_0=p0, position_1=p1) for p0, p1 in [[7, 0], [1, 2], [3, 4], [5, 6]]],
        dim='position')

    cm = plt.colormaps['tab20c']
    figure, axis = plt.subplots(figsize=(10,3), layout='constrained')
    axis.cla()
    x = np.arange(len(basepaired_position.position)) * (len(basepaired_position.basepair) + 4)
    labels = []
    x_ticks = []
    for i, basepair in enumerate(basepaired_position.basepair.values):
        axis.bar(x + i, basepaired_position.sel(basepair=basepair), width=1, label=basepair, color=cm.colors[i])
        labels += [basepair[0] + '\n' + basepair[1]] * len(x)
        x_ticks += (x+i).tolist()
    # ax.legend()
    axis.set_xticks(x_ticks, labels)
    axis.set_ylabel('Count')
    axis.set_title(title)

    axis.annotate('Basepair 8-1', xy=(8,-40), xycoords=('data','axes points'), ha='center')
    axis.annotate('Basepair 2-3', xy=(28,-40), xycoords=('data','axes points'), ha='center')
    axis.annotate('Basepair 4-5', xy=(48,-40), xycoords=('data','axes points'), ha='center')
    axis.annotate('Basepair 6-7', xy=(68,-40), xycoords=('data','axes points'), ha='center')
    savefile_path = save_path / title
    figure.savefig(savefile_path.with_suffix('.png'))
    figure.savefig(savefile_path.with_suffix('.pdf'))

def basepaired_sequence_subset_count(sequence_subset_count, save_path):
    title = (sequence_subset_count.name.replace('sequence_count', 'basepaired_sequence_count'))
    variable = sequence_subset_count.attrs['variable']

    sequence_subset_count2 = sequence_subset_count.reindex(**{variable: all_basepaired_subsets(), 'fill_value': 0})

    # bases = ['A', 'T', 'C', 'G']
    basepairs = ['AT', 'TA', 'CG', 'GC']
    # possible_sequence_subsets = [''.join(ss) for ss in
    #                              itertools.product(bases, basepairs, bases, bases, basepairs, bases)]
    # possible_sequence_subset_count = sequence_subset_count.reindex(sequence_subset=possible_sequence_subsets,
    #                                                                fill_value=0)

    sequence_subset_count_bp = xr.DataArray(0, dims=('bp0', 'bp1', 'bp2', 'bp3'),
                                            coords={'bp0': basepairs, 'bp1': basepairs, 'bp2': basepairs,
                                                    'bp3': basepairs}).astype(sequence_subset_count2.dtype)
    for bpc in tqdm.tqdm(itertools.product(basepairs, basepairs, basepairs, basepairs)):
        sequence_subset_count_bp.loc[dict(bp0=bpc[0], bp1=bpc[1], bp2=bpc[2], bp3=bpc[3])] = \
            sequence_subset_count2.sel(**{variable: (''.join(bpc)[1:] + ''.join(bpc)[0])})

    figure, axis = plt.subplots(figsize=(7.2, 6.5))#, layout='tight')
    # axis = axes[0]
    axis.cla()
    data = sequence_subset_count_bp.stack(bp02=('bp0', 'bp2')).stack(bp13=('bp1', 'bp3')).T
    image = axis.imshow(data.values, vmin=0, vmax=sequence_subset_count_bp.mean() * 2, cmap='coolwarm')
    # axis.images[0].set_data(data.values)

    for x, bp in zip(np.arange(0, 16, 4) + 1.5, basepairs):
        axis.annotate(bp[::-1], xy=(x, 16), xycoords=('data', 'data'), ha='center', va='center')
    for y, bp in zip(np.arange(0, 16, 4) + 1.5, basepairs):
        axis.annotate(bp, xy=(-1, y), xycoords=('data', 'data'), ha='center', va='center')
    for x, bp in zip(np.arange(0, 16), basepairs * 4):
        axis.annotate(bp, xy=(x, -1), xycoords=('data', 'data'), ha='center', va='center')
    for y, bp in zip(np.arange(0, 16), basepairs * 4):
        axis.annotate(bp, xy=(16, y), xycoords=('data', 'data'), ha='center', va='center')

    axis.set_xticks(np.arange(-0.5, 16, 4)[1:-1], minor=False)
    # axis.set_xticks(np.arange(0,16), basepairs*4, minor=True)
    axis.set_yticks(np.arange(-0.5, 16, 4)[1:-1], minor=False)
    # axis.set_yticks(np.arange(0,16), basepairs*4, minor=True)
    axis.tick_params(which="major", top=False, labeltop=False, right=False, labelright=False,
                     bottom=False, labelbottom=False, left=False, labelleft=False)
    axis.grid(which='major', color="w", linestyle='-', linewidth=2)
    for spine in axis.spines.values():
        spine.set_visible(False)
    axis.annotate('Basepair 1-8', xy=(7.5, 17), xycoords=('data', 'data'), ha='center', va='center')
    axis.annotate('Basepair 2-3', xy=(-2, 7.5), xycoords=('data', 'data'), ha='center', va='center',
                  rotation='vertical')
    axis.annotate('Basepair 4-5', xy=(7.5, -2), xycoords=('data', 'data'), ha='center', va='center')
    axis.annotate('Basepair 6-7', xy=(17, 7.5), xycoords=('data', 'data'), ha='center', va='center',
                  rotation='vertical')

    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(axis)
    cax = divider.append_axes("right", "4%", pad="15%")
    figure.colorbar(image, cax=cax)

    # figure.colorbar(image, aspect=30)
    # cax = figure.axes[-1]
    cax.set_ylabel('Count')
    cax.axes.ticklabel_format(scilimits=(0,0))
    for spine in cax.spines.values():
        spine.set_visible(False)

    fontsize = 8
    title2 = title + ' - scaled_mean_basepaired'
    figure.suptitle(title2, fontsize=fontsize)
    savefile_path = save_path / title2
    figure.savefig(savefile_path.with_suffix('.png'))
    figure.savefig(savefile_path.with_suffix('.pdf'))

    image.set_clim(sequence_subset_count_bp.min(), np.percentile(sequence_subset_count.values, 20))
    title2 = title + ' - scaled_bottom_20%_of_all'
    figure.suptitle(title2, fontsize=fontsize)
    savefile_path = save_path / title2
    figure.savefig(savefile_path.with_suffix('.png'))
    figure.savefig(savefile_path.with_suffix('.pdf'))

    image.set_clim(sequence_subset_count_bp.min(), sequence_subset_count.mean().item() * 2)
    title2 = title + ' - scaled_mean_of_all'
    figure.suptitle(title2, fontsize=fontsize)
    savefile_path = save_path / title2
    figure.savefig(savefile_path.with_suffix('.png'))
    figure.savefig(savefile_path.with_suffix('.pdf'))

    image.set_clim(sequence_subset_count_bp.min(), sequence_subset_count_bp.max())
    title2 = title + ' - scaled_min_max_basepaired'
    figure.suptitle(title2, fontsize=fontsize)
    savefile_path = save_path / title2
    figure.savefig(savefile_path.with_suffix('.png'))
    figure.savefig(savefile_path.with_suffix('.pdf'))




# def plot_basepaired_holliday_junction(data, save_path, **imshow_kwargs):
#
#     data = data.reindex(**{'sequence_subset': all_basepaired_subsets(), 'fill_value': np.nan})
#
#     # bases = ['A', 'T', 'C', 'G']
#     basepairs = ['AT', 'TA', 'CG', 'GC']
#
#     data_bp = xr.DataArray(0, dims=('bp0', 'bp1', 'bp2', 'bp3'),
#                            coords={'bp0': basepairs, 'bp1': basepairs, 'bp2': basepairs,
#                                    'bp3': basepairs}).astype(data.dtype)
#     for bpc in tqdm.tqdm(itertools.product(basepairs, basepairs, basepairs, basepairs)):
#         data_bp.loc[dict(bp0=bpc[0], bp1=bpc[1], bp2=bpc[2], bp3=bpc[3])] = \
#             data.sel(**{'sequence_subset': (''.join(bpc)[1:] + ''.join(bpc)[0])})
#
#     figure, axis = plt.subplots(figsize=(7.9, 6.5), layout='constrained')
#     # axis = axes[0]
#     axis.cla()
#     data_bp_stacked = data_bp.stack(bp02=('bp0', 'bp2')).stack(bp13=('bp1', 'bp3')).T
#     image = axis.imshow(data_bp_stacked.values, cmap='coolwarm', **imshow_kwargs)
#     # axis.images[0].set_data(data.values)
#
#     for x, bp in zip(np.arange(0, 16, 4) + 1.5, basepairs):
#         axis.annotate(bp[::-1], xy=(x, 16), xycoords=('data', 'data'), ha='center', va='center')
#     for y, bp in zip(np.arange(0, 16, 4) + 1.5, basepairs):
#         axis.annotate(bp, xy=(-1, y), xycoords=('data', 'data'), ha='center', va='center')
#     for x, bp in zip(np.arange(0, 16), basepairs * 4):
#         axis.annotate(bp, xy=(x, -1), xycoords=('data', 'data'), ha='center', va='center')
#     for y, bp in zip(np.arange(0, 16), basepairs * 4):
#         axis.annotate(bp, xy=(16, y), xycoords=('data', 'data'), ha='center', va='center')
#
#     axis.set_xticks(np.arange(-0.5, 16, 4)[1:-1], minor=False)
#     # axis.set_xticks(np.arange(0,16), basepairs*4, minor=True)
#     axis.set_yticks(np.arange(-0.5, 16, 4)[1:-1], minor=False)
#     # axis.set_yticks(np.arange(0,16), basepairs*4, minor=True)
#     axis.tick_params(which="major", top=False, labeltop=False, right=False, labelright=False,
#                      bottom=False, labelbottom=False, left=False, labelleft=False)
#     axis.grid(which='major', color="w", linestyle='-', linewidth=2)
#     for spine in axis.spines.values():
#         spine.set_visible(False)
#     axis.annotate('Basepair 1-8', xy=(7.5, 17), xycoords=('data', 'data'), ha='center', va='center')
#     axis.annotate('Basepair 2-3', xy=(-2, 7.5), xycoords=('data', 'data'), ha='center', va='center',
#                   rotation='vertical')
#     axis.annotate('Basepair 4-5', xy=(7.5, -2), xycoords=('data', 'data'), ha='center', va='center')
#     axis.annotate('Basepair 6-7', xy=(17, 7.5), xycoords=('data', 'data'), ha='center', va='center',
#                   rotation='vertical')
#
#     from mpl_toolkits.axes_grid1 import make_axes_locatable
#     divider = make_axes_locatable(axis)
#     cax = divider.append_axes("right", "4%", pad="15%")
#     figure.colorbar(image, cax=cax)
#
#     # figure.colorbar(image, aspect=30)
#     # cax = figure.axes[-1]
#     cax.set_ylabel(data.name)
#     cax.axes.ticklabel_format(scilimits=(0,0))
#     for spine in cax.spines.values():
#         spine.set_visible(False)
#
#     title = 'basepaired_HJ_' + data.name
#
#     fontsize = 8
#     # title2 = title + ' - scaled_mean_basepaired'
#     figure.suptitle(title, fontsize=fontsize)
#     savefile_path = save_path / title
#     figure.savefig(savefile_path.with_suffix('.png'))
#     figure.savefig(savefile_path.with_suffix('.pdf'))
#
#     # image.set_clim(sequence_subset_count_bp.min(), np.percentile(sequence_subset_count.values, 20))
#     # title2 = title + ' - scaled_bottom_20%_of_all'
#     # figure.suptitle(title2, fontsize=fontsize)
#     # savefile_path = save_path / title2
#     # figure.savefig(savefile_path.with_suffix('.png'))
#     # figure.savefig(savefile_path.with_suffix('.pdf'))
#     #
#     # image.set_clim(sequence_subset_count_bp.min(), sequence_subset_count.mean().item() * 2)
#     # title2 = title + ' - scaled_mean_of_all'
#     # figure.suptitle(title2, fontsize=fontsize)
#     # savefile_path = save_path / title2
#     # figure.savefig(savefile_path.with_suffix('.png'))
#     # figure.savefig(savefile_path.with_suffix('.pdf'))
#     #
#     # image.set_clim(sequence_subset_count_bp.min(), sequence_subset_count_bp.max())
#     # title2 = title + ' - scaled_min_max_basepaired'
#     # figure.suptitle(title2, fontsize=fontsize)
#     # savefile_path = save_path / title2
#     # figure.savefig(savefile_path.with_suffix('.png'))
#     # figure.savefig(savefile_path.with_suffix('.pdf'))
#
#


def plot_basepaired_holliday_junction(data, size=1, name=None, s2max=None, geometry='square',
                                      axis_facecolor="lightgrey", save_path=None, vmin=None, vmax=None):
    if name is None:
        name = data.name

    if data.ndim == 1:
        data = data.expand_dims('loop_dim', 0)

    data = xr.Dataset(dict(data=data))
    if np.isscalar(size):
        data['size'] = xr.DataArray(size, coords=data.data.coords)
    elif isinstance(size, (list, tuple)):
        data['size'] = xr.DataArray(np.full_like(data.data, np.array(size)[None,:].T, dtype='object'), dims=data.data.dims)
    elif isinstance(size, xr.DataArray):
        if size.ndim == 1:
            size = size.expand_dims({data.data.dims[0]: data.data.shape[0]}, 0)
        data['size'] = size

    if isinstance(geometry, str):
        data['geometry'] = xr.DataArray(geometry, coords=data.data.coords)
    elif isinstance(geometry, (list, tuple)):
        data['geometry'] = xr.DataArray(np.full_like(data.data, np.array(geometry)[None,:].T, dtype='object'), dims=data.data.dims)
    elif isinstance(geometry, xr.DataArray):
        if geometry.ndim == 1:
            geometry = geometry.expand_dims({data.data.dims[0]: data.data.shape[0]}, 0)
        data['geometry'] = geometry

    data = data.reindex(**{'sequence_subset': all_basepaired_subsets(), 'fill_value': np.nan})

    # bases = ['A', 'T', 'C', 'G']
    basepairs = ['AT', 'TA', 'CG', 'GC']

    bp0 = np.char.add(data.sequence_subset.str[7:8].values, data.sequence_subset.str[0:1].values)
    bp1 = data.sequence_subset.str[1:3].values
    bp2 = data.sequence_subset.str[3:5].values
    bp3 = data.sequence_subset.str[5:7].values
    data_bp = data.assign_coords(bp=('sequence_subset',pd.MultiIndex.from_arrays([bp0, bp1, bp2, bp3], names=['bp0','bp1','bp2','bp3'])))
    data_bp = data_bp.swap_dims(sequence_subset='bp').unstack('bp')
    data_bp = data_bp.sel(bp0=basepairs, bp1=basepairs, bp2=basepairs, bp3=basepairs)


    # data_bp = xr.DataArray(0, dims=('bp0', 'bp1', 'bp2', 'bp3'),
    #                        coords={'bp0': basepairs, 'bp1': basepairs, 'bp2': basepairs,
    #                                'bp3': basepairs}).astype(data.dtype)
    # for bpc in tqdm.tqdm(itertools.product(basepairs, basepairs, basepairs, basepairs)):
    #     data_bp.loc[dict(bp0=bpc[0], bp1=bpc[1], bp2=bpc[2], bp3=bpc[3])] = \
    #         data.sel(**{'sequence_subset': (''.join(bpc)[1:] + ''.join(bpc)[0])})

    figure, axis = plt.subplots(figsize=(7.9, 6.5), layout='constrained')
    # axis = axes[0]
    axis.cla()
    data_bp_stacked = data_bp.stack(bp13=('bp1', 'bp3')).stack(bp02=('bp0', 'bp2'))
    # image = axis.imshow(data_bp_stacked.values, cmap='coolwarm', **imshow_kwargs)

    data_bp_stacked2 = data_bp_stacked.copy()
    data_bp_stacked2 = data_bp_stacked2.assign_coords(dict(bp02i=('bp02', np.arange(16)), bp13i=('bp13', np.arange(16))))\
        .unstack(('bp13', 'bp02')).stack(bp0123=('bp0', 'bp1', 'bp2', 'bp3'))
    axis.cla()
    # axis.scatter(x=data_bp_stacked2.bp02i, y=data_bp_stacked2.bp13i, c=data_bp_stacked2.data.values, cmap='coolwarm',
    #              vmin=0, vmax=1, s=data_bp_stacked2.size_data.values, marker='s')

    x = data_bp_stacked2.bp02i.values
    y = data_bp_stacked2.bp13i.values
    from matplotlib.colors import Normalize
    if vmin is None:
        vmin = 0
    if vmax is None:
        vmax = 1
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = 'coolwarm'
    from matplotlib.cm import ScalarMappable
    scalar_mappable = ScalarMappable(norm=norm, cmap=cmap)
    # c = plt.get_cmap('coolwarm')(norm(data_bp_stacked2.data.values))

    for i in range(data_bp_stacked2.data.shape[0]):
        c = scalar_mappable.to_rgba(data_bp_stacked2.data[i].values)

        s = np.sqrt(data_bp_stacked2.size[i].values)
        if s2max is None:
            s2max = np.nanmax(s)**2
        s = np.minimum(s/np.sqrt(s2max),1) #s/s.max()

        g = data_bp_stacked2.geometry[i].values

        for xi, yi, ci, si, gi in zip(x,y,c,s,g):
            shape_kwargs = dict(facecolor=ci, linewidth=None)
            if gi == 'square':
                rect = plt.Rectangle([xi - si / 2, yi - si / 2], si, si, **shape_kwargs)
            elif gi == 'circle':
                rect = plt.Circle([xi, yi], radius=si/2, **shape_kwargs)
            elif gi == 'diamond':
                si2 = si/np.sqrt(2)
                rect = plt.Rectangle([xi, yi - si / 2], si2, si2, angle=45, **shape_kwargs)
            elif gi == 'triangle_top_left':
                rect = plt.Polygon([[xi - si / 2, yi - si / 2],
                                    [xi - si / 2, yi + si / 2],
                                    [xi + si / 2, yi - si / 2]], closed=True, **shape_kwargs)
            elif gi == 'triangle_bottom_right':
                rect = plt.Polygon([[xi + si / 2, yi + si / 2],
                                    [xi - si / 2, yi + si / 2],
                                    [xi + si / 2, yi - si / 2]], closed=True, **shape_kwargs)
            elif gi == 'semi_circle_top_left':
                rect = patches.Wedge([xi, yi], r=si / 2, theta1=-225, theta2=-45, **shape_kwargs)
            elif gi == 'semi_circle_bottom_right':
                rect = patches.Wedge([xi, yi], r=si / 2, theta1=-45, theta2=-225, **shape_kwargs)

            axis.add_patch(rect)
    axis.set_aspect(1)
    half_shape_size = np.nanmax(s)/2
    axis.set_xlim(x.min()-half_shape_size, x.max()+half_shape_size)
    axis.set_ylim(y.min()-half_shape_size, y.max()+half_shape_size)
    # axis.autoscale_view()

    # axis.images[0].set_data(data.values)

    for x, bp in zip(np.arange(0, 16, 4) + 1.5, basepairs):
        axis.annotate(bp[::-1], xy=(x, 16), xycoords=('data', 'data'), ha='center', va='center')
    for y, bp in zip(np.arange(0, 16, 4) + 1.5, basepairs):
        axis.annotate(bp, xy=(-1, y), xycoords=('data', 'data'), ha='center', va='center')
    for x, bp in zip(np.arange(0, 16), basepairs * 4):
        axis.annotate(bp, xy=(x, -1), xycoords=('data', 'data'), ha='center', va='center')
    for y, bp in zip(np.arange(0, 16), basepairs * 4):
        axis.annotate(bp, xy=(16, y), xycoords=('data', 'data'), ha='center', va='center')

    axis.set_xticks(np.arange(-0.5, 16, 4)[1:-1], minor=False)
    axis.set_yticks(np.arange(-0.5, 16, 4)[1:-1], minor=False)
    axis.tick_params(which="major", top=False, labeltop=False, right=False, labelright=False,
                     bottom=False, labelbottom=False, left=False, labelleft=False)
    axis.grid(which='major', color="w", linestyle='-', linewidth=2)
    # axis.grid(which='major', color="k", linestyle='-', linewidth=2)

    axis.set_xticks(np.arange(-0.5,16.5), minor=True)
    axis.set_yticks(np.arange(-0.5,16.5), minor=True)
    axis.tick_params(which="minor", top=False, labeltop=False, right=False, labelright=False,
                     bottom=False, labelbottom=False, left=False, labelleft=False)
    axis.grid(which='minor', color="k", linestyle='-', linewidth=0)

    for spine in axis.spines.values():
        spine.set_visible(False)
        spine.set_lw(2)

    axis.annotate('Basepair 1-8', xy=(7.5, 17), xycoords=('data', 'data'), ha='center', va='center')
    axis.annotate('Basepair 2-3', xy=(-2, 7.5), xycoords=('data', 'data'), ha='center', va='center',
                  rotation='vertical')
    axis.annotate('Basepair 4-5', xy=(7.5, -2), xycoords=('data', 'data'), ha='center', va='center')
    axis.annotate('Basepair 6-7', xy=(17, 7.5), xycoords=('data', 'data'), ha='center', va='center',
                  rotation='vertical')

    axis.set_facecolor(axis_facecolor)
    axis.invert_yaxis()

    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(axis)
    cax = divider.append_axes("right", "4%", pad="15%")
    figure.colorbar(scalar_mappable, cax=cax)

    # figure.colorbar(image, aspect=30)
    # cax = figure.axes[-1]
    if 'unit' in data.data.attrs and data.data.attrs['unit'] !='':
        unit_string = f' ({data.data.attrs["unit"]})'
    else:
        unit_string = ''
    cax.set_ylabel(name + unit_string)
    cax.axes.ticklabel_format(scilimits=(0,0))
    for spine in cax.spines.values():
        spine.set_visible(False)

    title = 'basepaired_HJ_' + name

    fontsize = 8
    # title2 = title + ' - scaled_mean_basepaired'
    figure.suptitle(title, fontsize=fontsize)
    if save_path is not None:
        savefile_path = save_path / title
        figure.savefig(savefile_path.with_suffix('.png'))
        figure.savefig(savefile_path.with_suffix('.pdf'))


def format_sequence_subset(sequence_subset):
    ss = sequence_subset
    return np.array(
    [f'  ||  ',
    f'  {ss[3]}{ss[4]}  ',
    f'-{ss[2]}  {ss[5]}-',
    f'-{ss[1]}  {ss[6]}-',
    f'  {ss[0]}{ss[7]}  ',
    f'  ||  ']
    )

def format_sequence_subsets(sequence_subsets):
    spacer = np.array(['  ']*6)
    sequence_subsets_formatted = np.array([format_sequence_subset(ss) for ss in sequence_subsets])
    spacers = np.array([spacer]*len(sequence_subsets))

    final = char_add([np.char.add(*s) for s in zip(sequence_subsets_formatted, spacers)])
    final = '\n'.join(final)
    return final

# From itertools recipes
def batched(iterable, n):
    "Batch data into tuples of length n. The last batch may be shorter."
    # batched('ABCDEFG', 3) --> ABC DEF G
    if n < 1:
        raise ValueError('n must be at least one')
    it = iter(iterable)
    while batch := tuple(itertools.islice(it, n)):
        yield batch

def print_sequence_subsets(sequence_subsets, rows=None):
    if isinstance(sequence_subsets, str):
        sequence_subsets = [sequence_subsets]
    if rows is not None:
        sequence_subsets = batched(sequence_subsets, rows)
    else:
        sequence_subsets = [sequence_subsets]
    for sequence_subsets_row in sequence_subsets:
        print(format_sequence_subsets(sequence_subsets_row), '\n\n')


def char_add(strings):
    for i, s in enumerate(strings):
        if i == 0:
            final = s
        else:
            final = np.char.add(final, s)
    return final
