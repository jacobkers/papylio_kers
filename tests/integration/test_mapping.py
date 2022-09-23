import pytest
import tifffile
import numpy as np
from skimage.transform import SimilarityTransform
from trace_analysis.mapping.mapping import Mapping2
from trace_analysis.mapping.geometric_hashing import GeometricHashTable


def test_geometric_hash_table():
    translation = np.array([256, 10])
    rotation = 125/360*2*np.pi
    scale = np.array([10, 10])
    transformation = SimilarityTransform(translation=translation, rotation=rotation, scale=scale)
    mapping = Mapping2.simulate(number_of_points=200, transformation=transformation,
                 bounds=([0, 0], [256, 512]), crop_bounds=((50, 200), None), fraction_missing=(0.1, 0.1),
                 error_sigma=(0.5, 0.5), shuffle=True, seed=10252)

    hash_table = GeometricHashTable([mapping.destination], source_vertices=mapping.source_vertices,
                 number_of_source_bases=20, number_of_destination_bases='all',
                 tuple_size=4, maximum_distance_source=100, maximum_distance_destination=1000)

    found_transformation = hash_table.query(mapping.source, distance=15, alpha=0.9, sigma=10, K_threshold=10e9, hash_table_distance_threshold=0.01,
                     magnification_range=None, rotation_range=None)

    assert (np.abs(found_transformation.translation - transformation.translation) < 50).all()
    assert np.abs(found_transformation.rotation - transformation.rotation) < 0.01
    assert (np.abs(found_transformation.scale - transformation.scale) < 0.2).all()

    found_transformation = hash_table.query_tuple_transformations([mapping.source], hash_table_distance_threshold=0.01,
                                                                  parameters=['translation', 'rotation', 'scale'])
    assert (np.abs(found_transformation.translation - transformation.translation) < 10).all()
    assert np.abs(found_transformation.rotation - transformation.rotation) < 0.01
    assert (np.abs(found_transformation.scale - transformation.scale) < 0.1).all()

def test_geometric_hash_table():
    translation = np.array([256, 10])
    rotation = 125 / 360 * 2 * np.pi
    scale = np.array([10, 10])
    transformation = SimilarityTransform(translation=translation, rotation=rotation, scale=scale)
    mapping = Mapping2.simulate(number_of_points=200, transformation=transformation,
                                bounds=([0, 0], [256, 512]), crop_bounds=((50, 200), None), fraction_missing=(0.1, 0.1),
                                error_sigma=(0.5, 0.5), shuffle=True, seed=10252)
    mapping.geometric_hashing(method='test_one_by_one', tuple_size=4, maximum_distance_source=100, maximum_distance_destination=1000,
                              distance=15, alpha=0.9, sigma=10, K_threshold=10e9, hash_table_distance_threshold=0.01,
                              magnification_range=None, rotation_range=None)

    assert (np.abs(mapping.transformation.translation - mapping.transformation_correct.translation) < 50).all()
    assert np.abs(mapping.transformation.rotation - mapping.transformation_correct.rotation) < 0.01
    assert (np.abs(mapping.transformation.scale - mapping.transformation_correct.scale) < 0.2).all()

    mapping.transformation = mapping.transformation_correct
    mapping.geometric_hashing(method='abundant_transformations', tuple_size=4, maximum_distance_source=100, maximum_distance_destination=1000,
                              hash_table_distance_threshold=0.01,
                              parameters=['translation', 'rotation', 'scale']
                              )

    assert (np.abs(mapping.transformation.translation - mapping.transformation_correct.translation) < 10).all()
    assert np.abs(mapping.transformation.rotation - mapping.transformation_correct.rotation) < 0.01
    assert (np.abs(mapping.transformation.scale - mapping.transformation_correct.scale) < 0.1).all()


