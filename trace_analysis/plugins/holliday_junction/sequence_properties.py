import re
import numpy as np
import xarray as xr


def datatype_conversion(dtype=list):
    def datatype_conversion_decorator(func):
        def wrapper(sequences, *args, **kwargs):
            if isinstance(sequences, xr.DataArray):
                coords = sequences.coords
                sequences = dtype(sequences.values)
                return xr.DataArray(func(sequences, *args, **kwargs), coords=coords)
            elif isinstance(sequences, np.ndarray):
                sequences = dtype(sequences)
                return np.array(func(sequences, *args, **kwargs))
            else:
                sequences = dtype(sequences)
                return list(func(sequences, *args, **kwargs))
        return wrapper
    return datatype_conversion_decorator


@datatype_conversion()
def fraction_GC(sequences):
   return [len(re.findall('G|C', sequence))/len(sequence) for sequence in sequences]

@datatype_conversion()
def number_of_neighboring_bases(sequences, base_type):
    if base_type == 'purine':
        search_string = 'AG|GA|AA|GG'
    elif base_type == 'pyrimidine':
        search_string = 'TC|CT|TT|CC'
    elif len(base_type) == 1:
        search_string = base_type + base_type
    return [len(re.findall(search_string, sequence + sequence[0])) for sequence in sequences]

@datatype_conversion()
def number_of_bases(sequences, base_type, positions='all'):
    if not positions=='all':
        sequences_new = []
        sequences_new.append([])
    if base_type == 'purine':
        search_string = 'AG|GA|AA|GG'
    elif base_type == 'pyrimidine':
        search_string = 'TC|CT|TT|CC'
    elif len(base_type) == 1:
        search_string = base_type + base_type
    return [len(re.findall(search_string, sequence + sequence[0])) for sequence in sequences]

@datatype_conversion(np.array)
def base_count(sequences, positions='all', bases='all'):
    # TODO: make it possible to
    if positions == 'all':
        positions = np.arange(len(sequences[0]))
    else:
        positions = np.array(positions)

    if bases == 'all':
        bases = np.array(['A', 'T', 'C', 'G'])
    elif bases == 'purines':
        bases = np.array(['A', 'G'])
    elif bases == 'pyrimidines':
        bases = np.array(['T', 'C'])
    else:
        bases = np.array(bases)

    sequence_array = np.array(sequences).astype('U').view('U1').reshape(-1, len(sequences[0]))
    return (sequence_array[:, positions, None] == bases[None, None, :]).any(axis=2).sum(axis=1)

@datatype_conversion(np.array)
def base_combination_count(sequences, position_pairs, bases='purines'):
    base_combination_count = \
        np.array([(base_count(sequences, positions=position_pair, bases=bases)==2).astype(int) for position_pair in position_pairs])
    return base_combination_count.sum(axis=0)
