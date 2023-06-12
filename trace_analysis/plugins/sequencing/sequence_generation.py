import numpy as np

def generate_sequences(base_composition):
    """

    Parameters
    ----------
    base_composition : list of str
        A list with the length of the sequence, where each entry indicates the possible bases at that position.
        e.g. ['N', 'AC', 'GC', 'ACTG']

    Returns
    -------
    sequences : numpy.array
        Array of sequences
    """

    base_composition_2 = []
    for bases in base_composition:
        if bases == 'N':
            bases = 'ACTG'
        base_composition_2.append(list(bases))

    return np.array(np.meshgrid(*base_composition_2)).T.reshape(-1, len(base_composition_2)).view(f'U{len(base_composition_2)}').squeeze()