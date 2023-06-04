import itertools
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

from trace_analysis.plugins.holliday_junction.holliday_junction import basepaired_subsets

def sequence_subset_with_comparable_structure(sequence_subset):
    # sequence_subset = sequence_subset[7] + sequence_subset[:7]
    sequence_subsets = []
    for i in range(4):
        sequence_subsets.append(sequence_subset[i*2:] + sequence_subset[:i*2])
    for ss in sequence_subsets.copy():
        change_base_within_basepair = str.maketrans('CG', 'GC')
        sequence_subsets.append(ss.translate(change_base_within_basepair))
    for ss in sequence_subsets.copy():
        change_base_within_basepair = str.maketrans('AT', 'TA')
        sequence_subsets.append(ss.translate(change_base_within_basepair))
    for ss in sequence_subsets.copy():
        switch_base_pairs = str.maketrans('ATCG', 'GCAT')
        sequence_subsets.append(ss.translate(switch_base_pairs))
    return set(sequence_subsets)


def sequence_subset_structure_category():
    subset_category_dict = {}
    sequence_subsets = set(basepaired_subsets())

    structure_index = 0
    while len(sequence_subsets) > 0:
        ss_with_structure = sequence_subset_with_comparable_structure(next(iter(sequence_subsets)))
        for ss in ss_with_structure:
            subset_category_dict[ss] = structure_index
        sequence_subsets = sequence_subsets.difference(ss_with_structure)
        structure_index += 1
    return subset_category_dict

def unique_subset_structures():
    unique_subset_category_dict = {}
    for key, value in sequence_subset_structure_category().items():
        if value not in unique_subset_category_dict.values():
            unique_subset_category_dict[key] = value
    return unique_subset_category_dict

def basepairing(sequences):
    s = sequences.astype('U').view('U1').reshape(-1, 8)
    s1 = xr.DataArray(s, coords=[('sequence', sequences), ('nucleotide_1', range(8))])
    s2 = xr.DataArray(s, coords=[('sequence', sequences), ('nucleotide_2', range(8))])
    s11, s22 = xr.broadcast(s1, s2)
    basepairs = xr.concat([s11,s22], dim='b').str.join(dim='b')
    basepairing = (basepairs == 'AT').astype(int) + (basepairs == 'TA').astype(int) + (basepairs == 'CG').astype(int) + (basepairs == 'GC').astype(int)
    return basepairing


def plot_basepairing(basepairing, title='', ax=None, save_path=None, max_linewidth=10):
    data = basepairing.mean('sequence')

    n = 8
    phi = 9/8*np.pi
    angle = np.linspace(0,2*np.pi,n+1)[:n]
    coordinates = np.vstack([np.sin(angle+phi), np.cos(angle+phi)]).T

    if ax is None:
        fig, ax = plt.subplots(tight_layout=True)
    else:
        fig = ax.figure
    ax.set_aspect(1)
    ax.axes.set_axis_off()
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-1.1, 1.1)
    for i,j in itertools.combinations(range(n),2):
        # ax.annotate("",
        #             xy=coordinates[i], xycoords='data',
        #             xytext=coordinates[j], textcoords='data',
        #             arrowprops=dict(arrowstyle="-"),
        #             )
        ax.plot(*coordinates[[i,j]].T, c='k', linewidth=data[i,j]*max_linewidth)
    ax.scatter(*coordinates.T, s=400, c='k')
    for i, c in enumerate(coordinates):
        ax.annotate(i+1, c, ha='center', va='center', c='white')
        # ax.text(*c, i + 1, size=15, ha='center', va='center', c='white',
        #         bbox=dict(boxstyle="circle,pad=0.3", fc="cyan", ec="b", lw=2))
    ax.set_title(title)

    if save_path is not None:
        fig.savefig(save_path.joinpath(f'Basepairing - '+ title +'.png'))


def plot_basepairing_individual(sequences, name='', rows=5, columns=10, titles=None, save_path=None, max_linewidth=5):
    fig, axes = plt.subplots(rows, columns, figsize=(columns*2.5, rows*2.5))#, tight_layout=True)
    fig.subplots_adjust(hspace=0.3, wspace=0.3, left=0.02, right=0.98, bottom=0.02, top=0.95)
    axes = axes.flatten()

    for ax in axes:
        ax.set_axis_off()

    for i, (sequence, ax) in enumerate(zip(sequences, axes)):
        if titles is not None:
            title = f'{titles[i]} - {sequence}'
        else:
            title = f'{sequence}'
        plot_basepairing(basepairing(np.array([sequence])), title=title, ax=ax, max_linewidth=max_linewidth)
    if save_path is not None:
        fig.savefig(save_path.joinpath(f'Basepairing - {name} - individual.png'))


