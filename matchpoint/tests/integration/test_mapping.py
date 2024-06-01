import pytest
import tifffile
import pytest
import numpy as np
from skimage.transform import SimilarityTransform, AffineTransform
from matchpoint import MatchPoint


@pytest.fixture
def mapping():
    translation = np.array([10, -10])
    rotation = 1 / 360 * 2 * np.pi
    scale = [0.98, 0.98]
    transformation = SimilarityTransform(translation=translation, rotation=rotation, scale=scale)
    mapping = MatchPoint.simulate(number_of_points=200, transformation=transformation,
                                bounds=([0, 0], [256, 512]), crop_bounds=((50, 200), None), fraction_missing=(0.1, 0.1),
                                error_sigma=(0.5, 0.5), shuffle=True, seed=10252)
    return mapping


def test_determine_matched_pairs(mapping):
    mapping.transformation = mapping.transformation_correct
    mapping.determine_matched_pairs(distance_threshold=5)
    assert mapping.number_of_matched_points == 23

def test_vertices(mapping):
    assert (mapping.get_source_vertices() == mapping.source_vertices).all()
    assert (mapping.get_destination_vertices() == mapping.destination_vertices).all()

    assert (mapping.get_source_vertices(crop=True) == mapping.source_cropped_vertices).all()
    assert (mapping.get_destination_vertices(crop=True) == mapping.destination_cropped_vertices).all()


def test_iterative_closest_point():
    mapping = MatchPoint.simulate()
    mapping.transformation = AffineTransform(translation=(256,0))
    mapping.iterative_closest_point(5)
    assert mapping.transformation_is_similar_to_correct_transformation(translation_error=1, rotation_error=0.001, scale_error=0.01)


def test_iterative_closest_point_polynomial():
    mapping = MatchPoint.simulate()
    mapping.transformation = AffineTransform(translation=(256,0))
    mapping.transformation_type = 'polynomial'
    mapping.iterative_closest_point(5)
    #TODO add comparison for polynomial transforms, or based on point distances

@pytest.fixture()
def cross_correlation_mapping():
    transformation = SimilarityTransform(translation=[50, 50],  rotation=1 / 360 * 2 * np.pi, scale=[1, 1])
    mapping = MatchPoint.simulate(number_of_points=200, transformation=transformation,
                 bounds=[[0, 0], [256, 512]], crop_bounds=([[50,50], [150,150]], None), fraction_missing=(0.1, 0.1),
                 error_sigma=(0.5, 0.5), shuffle=True, seed=10532, show_correct=True)
    return mapping


def test_cross_correlation_space(cross_correlation_mapping):
    mapping = cross_correlation_mapping
    mapping.cross_correlation(space='source')
    assert mapping.transformation_is_similar_to_correct_transformation(translation_error=5, rotation_error=0.02, scale_error=0.03)
    mapping.cross_correlation(space='destination')
    assert mapping.transformation_is_similar_to_correct_transformation(translation_error=5, rotation_error=0.02, scale_error=0.03)


def test_cross_correlation_normalized(cross_correlation_mapping):
    mapping = cross_correlation_mapping
    mapping.cross_correlation(normalize=True, subtract_background=False)
    assert mapping.transformation_is_similar_to_correct_transformation(translation_error=5, rotation_error=0.02, scale_error=0.03)


def test_cross_correlation_subtract_background(cross_correlation_mapping):
    mapping = cross_correlation_mapping
    mapping.cross_correlation(subtract_background='minimum_filter')
    assert mapping.transformation_is_similar_to_correct_transformation(translation_error=5, rotation_error=0.02, scale_error=0.03)
    mapping.cross_correlation(subtract_background='median_filter')
    assert mapping.transformation_is_similar_to_correct_transformation(translation_error=5, rotation_error=0.02, scale_error=0.03)
    mapping.cross_correlation(subtract_background='expected_signal')
    assert mapping.transformation_is_similar_to_correct_transformation(translation_error=5, rotation_error=0.02,
                                                                       scale_error=0.03)
    mapping.cross_correlation(subtract_background='expected_signal_rough')
    assert mapping.transformation_is_similar_to_correct_transformation(translation_error=5, rotation_error=0.02,
                                                                       scale_error=0.03)
    mapping.cross_correlation(subtract_background=True)
    assert mapping.transformation_is_similar_to_correct_transformation(translation_error=5, rotation_error=0.02,
                                                                       scale_error=0.03)


def test_geometric_hash_table():
    translation = np.array([256, 10])
    rotation = 125 / 360 * 2 * np.pi
    scale = np.array([10, 10])
    transformation = SimilarityTransform(translation=translation, rotation=rotation, scale=scale)
    mapping = MatchPoint.simulate(number_of_points=200, transformation=transformation,
                                bounds=([0, 0], [256, 512]), crop_bounds=((50, 200), None), fraction_missing=(0.1, 0.1),
                                error_sigma=(0.5, 0.5), shuffle=True, seed=10252)
    mapping.geometric_hashing(method='one_by_one', tuple_size=4, maximum_distance_source=100, maximum_distance_destination=1000,
                              alpha=0.9, sigma=10, K_threshold=10e9, hash_table_distance_threshold=0.01,
                              magnification_range=None, rotation_range=None)

    assert mapping.transformation_is_similar_to_correct_transformation(translation_error=50, rotation_error=0.01,
                                                                       scale_error=0.2)

    mapping.transformation = transformation # mapping.transformation_correct
    mapping.geometric_hashing(method='abundant_transformations', tuple_size=4, maximum_distance_source=100, maximum_distance_destination=1000,
                              hash_table_distance_threshold=0.01,
                              parameters=['translation', 'rotation', 'scale']
                              )

    assert mapping.transformation_is_similar_to_correct_transformation(translation_error=10, rotation_error=0.01, scale_error=0.01)


def test_kernel_correlations():
    translation = np.array([10, -10])
    rotation = 1 / 360 * 2 * np.pi
    scale = [0.98, 0.98]
    transformation = AffineTransform(translation=translation, rotation=rotation, scale=scale)
    mapping = MatchPoint.simulate(number_of_points=10000, transformation=transformation,
                                bounds=([0, 0], [256, 512]), crop_bounds=(None, None), fraction_missing=(0.95, 0.6),
                                error_sigma=(0.5, 0.5), shuffle=True, seed=10252)

    minimization_bounds = ((0.97, 1.02), (-0.05, 0.05), (-20, 20), (-20, 20))
    mapping.kernel_correlation(minimization_bounds, sigma=1, crop=False, plot=False,
                               strategy='best1bin', maxiter=1000, popsize=50, tol=0.01, mutation=0.25, recombination=0.7,
                               seed=None, callback=None, disp=False, polish=True, init='sobol', atol=0,
                               updating='immediate', workers=1, constraints=())

    assert mapping.transformation_is_similar_to_correct_transformation(translation_error=1, rotation_error=0.001, scale_error=0.001)


def test_save(shared_datadir):
    mapping = MatchPoint.simulate()
    mapping.save(shared_datadir / 'mapping' / 'test_mapping.nc')
    # mapping = MatchPoint.load(r'C:\Users\ivoseverins\Scan 123 - HJ general - Tile 1101.mapping')
    # mapping.save(r'C:\Users\ivoseverins\test.nc.mapping', filetype='nc')
    return mapping


def test_load(shared_datadir):
    mapping_saved = test_save(shared_datadir)
    mapping_loaded = MatchPoint.load(shared_datadir / 'mapping' / 'test_mapping.nc')
    assert mapping_loaded == mapping_saved
