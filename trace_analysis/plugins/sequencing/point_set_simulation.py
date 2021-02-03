import numpy as np

def simulate_point_set(number_of_points, bounds = [[0,1],[0,1]]):
    bounds = np.array(bounds)
    dimensions = len(bounds)
    unit_coordinates = np.random.rand(number_of_points, dimensions)
    coordinates = unit_coordinates * (bounds[:,1]-bounds[:,0]) + bounds[:,0]

    return coordinates


def random_selection_from_point_set(coordinates, fraction):
    """ Obtain a random part of a point set
    
    Parameters
    ----------
    coordinates : Nx2 numpy.ndarray
        Coordinates to sample from
    fraction: float
        Fraction of points that should be selected

    Returns
    -------
    Nx2 numpy.ndarray
        Selected coordinates
    """

    random_generator = np.random.default_rng()
    size = round(fraction*len(coordinates))
    return random_generator.choice(coordinates, size, replace=False, axis=0, shuffle=False)

