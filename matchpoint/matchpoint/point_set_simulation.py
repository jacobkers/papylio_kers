import numpy as np
from skimage.transform import AffineTransform, PolynomialTransform

def simulate_point_set(number_of_points, bounds=([0, 0], [1, 1]), random_generator=None):
    """

    Parameters
    ----------
    number_of_points : int
        Number of points to simulate
    bounds : 2x2 numpy.ndarray or list
        Bounds for coordinate simulation; structured like coordinates,
        i.e. columns are x and y dimensions, rows are minimum and maximum values.

    Returns
    -------
    Nx2 numpy.ndarray
        Coordinates with error applied

    """
    if random_generator is None:
        random_generator = np.random.default_rng()

    bounds = np.array(bounds)
    dimensions = bounds.shape[1]
    unit_coordinates = random_generator.random((number_of_points, dimensions))
    coordinates = unit_coordinates * (bounds[1] - bounds[0]) + bounds[0]
    return coordinates


def random_selection_from_point_set(coordinates, fraction, random_generator=None):
    """ Obtain a random part of a point set

    Parameters
    ----------
    coordinates : Nx2 numpy.ndarray
        Coordinates; N is the number of points
    fraction : float
        Fraction of points that should be selected

    Returns
    -------
    Nx2 numpy.ndarray
        Subset of coordinates
    """
    if random_generator is None:
        random_generator = np.random.default_rng()
    size = round(fraction * len(coordinates))
    return random_generator.choice(coordinates, size, replace=False, axis=0, shuffle=False)


def add_uncertainty_to_point_set(coordinates, sigma=1, random_generator=None):
    """ Add random errors to the coordinates of a point set.
    For each point an error is randomly chosen from within a circle with radius maximum_error.

    Parameters
    ----------
    coordinates : Nx2 numpy.ndarray
        Coordinates; N is the number of points
    maximum_error : float
        Maximum random error applied to each point.

    Returns
    -------
    Nx2 numpy.ndarray
        Coordinates with error applied
    """

    # random_generator = np.random.default_rng()
    # error_magnitudes = random_generator.random(len(coordinates)) * maximum_error
    # error_angles = random_generator.random(len(coordinates)) * 2 * np.pi
    # errors = error_magnitudes[:, np.newaxis] * np.column_stack([np.cos(error_angles), np.sin(error_angles)])

    if random_generator is None:
        random_generator = np.random.default_rng()
    errors = random_generator.normal(0, sigma, size=coordinates.shape)
    # errors = np.random.normal(0, sigma, size=coordinates.shape)

    return coordinates + errors


def crop_point_set(coordinates, bounds):
    """ Crop point set

    Parameters
    ----------
    coordinates : Nx2 numpy.ndarray
        Coordinates; N is the number of points
    bounds : 2x2 numpy.ndarray or list
        Bounds used for cropping; structured like coordinates,
        i.e. columns are x and y dimensions, rows are minimum and maximum values.

    Returns
    -------
    Nx2 numpy.ndarray
        Cropped coordinates

    """
    bounds = np.array(bounds)
    crop_selection = np.all(np.hstack([coordinates > bounds[0], coordinates < bounds[1]]), axis=1)
    return coordinates[crop_selection]


def simulate_mapping_test_point_set(number_of_points, transformation, bounds=([0, 0], [1, 1]),
                                    crop_bounds=(None, None),
                                    fraction_missing=(0,0),
                                    error_sigma=(0,0), shuffle=True, seed=45896):
    """Simulate test point set for mapping

    A source point set is randomly generated between the given source_bounds. To obtain a corresponding destination
    point set the source is cropped using source_crop_bounds and transformed using transformation.
    Errors can be introduced by adding moving points randomly with respect to their center
    or by removing points from source or destination point sets.

    Note:   For now the random position errors are drawn from a uniform distribution within the maximum_error
            around each point. Possibly a normal distribution would be more appropriate.

    Parameters
    ----------
    number_of_source_points : int
        Number of coordinates in the source point set
    transformation : skimage.transform.AffineTransform or skimage.transform.PolynomialTransform
        Transformation applied to obtain the destination point set from the source point set.
    source_bounds
        Bounds of the source point set, structured like coordinates,
        i.e. columns are x and y dimensions, rows are minimum and maximum values.
    source_crop_bounds
        Bounds used for cropping of the source point set to obtain the destination point set;
        structured like coordinates, i.e. columns are x and y dimensions, rows are minimum and maximum values.
    fraction_missing_source : float
        Fraction of points that is deleted only in the source.
    fraction_missing_destination
        Fraction of points that is deleted only in the destination.
    maximum_error_source : float
        Maximum position error introduced in the source point set,
        i.e. the maximum shift from the original position in any direction.
    maximum_error_destination : float
        Maximum position error introduced in the destination point set,
        i.e. the maximum shift from the original position in any direction.
    shuffle : bool
        If true then the points in the destination will be shuffled.
        If false the order of the destination points in the source and destination will be identical.
        Default is 'True'.

    Returns
    -------
    source : Nx2 numpy.ndarray
        Coordinates of the source point set
    destination : Nx2 numpy.ndarray
        Coordinates of the destination point set

    """
    random_generator = np.random.default_rng(seed)

    number_of_point_sets = 2

    complete_point_set = simulate_point_set(number_of_points, bounds, random_generator)
    point_sets = []
    for i in range(number_of_point_sets):
        if crop_bounds[i] is not None:
            crop_bound = np.array(crop_bounds[i])
            point_set = crop_point_set(complete_point_set, crop_bound)
        else:
            point_set = complete_point_set.copy()

        point_set = random_selection_from_point_set(point_set, 1 - fraction_missing[i], random_generator)
        point_set = add_uncertainty_to_point_set(point_set, error_sigma[i], random_generator)

        # We could also transform both source and destination
        if i > 0:
            point_set = transformation(point_set)

            if shuffle:
                random_generator.shuffle(point_set, axis=0)

        point_sets.append(point_set)

    return point_sets


if __name__ == "__main__":
    from matchpoint.core import MatchPoint

    number_of_points = 100
    transformation = AffineTransform(scale=[0.75, 0.75], rotation=4 / 360 * 2 * np.pi, translation=[100, 0])
    bounds = [[0, 0], [100, 200]]
    crop_bounds = (None, [[25, 50], [100, 175]])
    fraction_missing = (0, 0)
    error_sigma = (0, 0)
    shuffle = True

    destination, source = simulate_mapping_test_point_set(number_of_points, transformation.inverse,
                                                          bounds, crop_bounds, fraction_missing, error_sigma, shuffle)

    m = MatchPoint(source, destination)
    m.transformation = transformation
    m.show_mapping_transformation(show_source=True)
