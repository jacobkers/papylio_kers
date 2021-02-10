import numpy as np
from skimage.transform import AffineTransform, PolynomialTransform
from trace_analysis.mapping.mapping import Mapping2


def simulate_point_set(number_of_points, bounds=([0, 0], [1, 1])):
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
    bounds = np.array(bounds)
    dimensions = bounds.shape[1]
    unit_coordinates = np.random.rand(number_of_points, dimensions)
    coordinates = unit_coordinates * (bounds[1] - bounds[0]) + bounds[0]
    return coordinates


def random_selection_from_point_set(coordinates, fraction):
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

    random_generator = np.random.default_rng()
    size = round(fraction * len(coordinates))
    return random_generator.choice(coordinates, size, replace=False, axis=0, shuffle=False)


def add_uncertainty_to_point_set(coordinates, maximum_error):
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

    random_generator = np.random.default_rng()
    error_magnitudes = random_generator.random(len(coordinates)) * maximum_error
    error_angles = random_generator.random(len(coordinates)) * 2 * np.pi
    errors = error_magnitudes[:, np.newaxis] * np.column_stack([np.cos(error_angles), np.sin(error_angles)])
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


def simulate_mapping_test_point_set(number_of_source_points, transformation, source_bounds=([0, 0], [1, 1]),
                                    source_crop_bounds=None,
                                    fraction_missing_source=0, fraction_missing_destination=0,
                                    maximum_error_source=0, maximum_error_destination=0, shuffle=True):
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
    source = simulate_point_set(number_of_source_points, source_bounds)
    if source_crop_bounds:
        source_crop_bounds = np.array(source_crop_bounds)
        destination = crop_point_set(source, source_crop_bounds)
    else:
        destination = source.copy()

    destination = transformation(destination)

    source = random_selection_from_point_set(source, 1 - fraction_missing_source)
    destination = random_selection_from_point_set(destination, 1 - fraction_missing_destination)

    # np.max(bounds[:, 1] - bounds[:, 0]) * fraction
    source = add_uncertainty_to_point_set(source, maximum_error_source)
    destination = add_uncertainty_to_point_set(destination, maximum_error_destination)

    if shuffle:
        random_generator = np.random.default_rng()
        random_generator.shuffle(destination, axis=0)

    return source, destination


if __name__ == "__main__":
    number_of_source_points = 100
    transformation = AffineTransform(scale=[0.75, 0.75], rotation=4 / 360 * 2 * np.pi, translation=[100, 0])
    source_bounds = [[0, 0], [100, 200]]
    source_crop_bounds = ([25, 50], [100, 175])
    fraction_missing_source = 0
    fraction_missing_destination = 0
    maximum_error_source = 0
    maximum_error_destination = 0
    shuffle = True

    destination, source = simulate_mapping_test_point_set(number_of_source_points, transformation.inverse,
                                                          source_bounds, source_crop_bounds,
                                                          fraction_missing_source, fraction_missing_destination,
                                                          maximum_error_source, maximum_error_destination, shuffle)
    m = Mapping2(source, destination)
    m.transformation = transformation
    m.show_mapping_transformation(show_source=True)
