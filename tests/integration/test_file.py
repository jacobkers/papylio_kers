import pytest
import tifffile
import numpy as np


@pytest.fixture
def experiment(shared_datadir):
    from trace_analysis import Experiment
    return Experiment(shared_datadir / 'BN_TIRF')

@pytest.fixture
def file(experiment):
    return experiment.files[1]

@pytest.fixture
def experiment_output(shared_datadir):
    from trace_analysis import Experiment
    return Experiment(shared_datadir / 'BN_TIRF_output_test_file')

@pytest.fixture
def file_output(experiment_output):
    return experiment_output.files[0]


def test_projection_image(file, shared_datadir):
    image_newly_made = file.projection_image
    image_from_original_file = tifffile.imread(shared_datadir / 'BN_TIRF_output_test_file' / 'TIRF 561 0001_ave_f0-10_i0.tif')
    assert (image_newly_made == image_from_original_file).all()
    image_loaded = file.projection_image
    assert (image_loaded == image_from_original_file).all()

def test_average_image(file, shared_datadir):
    image_newly_made = file.average_image
    image_from_original_file = tifffile.imread(shared_datadir / 'BN_TIRF_output_test_file' / 'TIRF 561 0001_ave_f0-10_i0.tif')
    assert (image_newly_made == image_from_original_file).all()
    image_loaded = file.average_image
    assert (image_loaded == image_from_original_file).all()

def test_maximum_projection_image(file, shared_datadir):
    image_newly_made = file.maximum_projection_image
    image_from_original_file = tifffile.imread(shared_datadir / 'BN_TIRF_output_test_file' / 'TIRF 561 0001_max_f0-10_i0.tif')
    assert (image_newly_made == image_from_original_file).all()
    image_loaded = file.maximum_projection_image
    assert (image_loaded == image_from_original_file).all()

def test_find_molecules(file):
    file.find_coordinates()

def test_extract_traces(file):
    file.find_coordinates()
    file.extract_traces()

def test_property_coordinates(file_output):
    file_output.coordinates
