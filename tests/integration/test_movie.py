import pytest
import tifffile
import numpy as np
from trace_analysis.movie.movie import Movie


@pytest.fixture
def movie(shared_datadir):
    movie = Movie(shared_datadir / 'BN_TIRF' / 'TIRF 561 0001.tif')
    movie.rot90 = 1
    return movie


@pytest.fixture
def experiment(shared_datadir):
    from trace_analysis import Experiment
    return Experiment(shared_datadir / 'BN_TIRF')


def test_movie_name(movie):
    assert movie.name == 'TIRF 561 0001'


def test_make_projection_image(movie, shared_datadir):
    image_from_method = movie.make_projection_image(projection_type='average', frame_range=(0, 20), illumination=None, write=True,
                                                    return_image=True, flatten_channels=True)
    assert (shared_datadir / 'BN_TIRF' / 'TIRF 561 0001_ave_f0-20_i0.tif').is_file()
    image_from_file = tifffile.imread(shared_datadir / 'BN_TIRF' / 'TIRF 561 0001_ave_f0-20_i0.tif')
    assert (image_from_file == image_from_method).all()
    image_from_original_file = tifffile.imread(shared_datadir / 'BN_TIRF_result' / 'TIRF 561 0001_ave_f0-20_i0.tif')
    assert (image_from_original_file == image_from_method).all()
    raw_images = tifffile.imread(shared_datadir / 'BN_TIRF' / 'TIRF 561 0001.tif', key=range(0, 20))
    raw_images = np.rot90(raw_images, axes=(1,2))
    assert ((image_from_file - raw_images.mean(axis=0)) < 1e-4).all()  # Not sure 1e-4 is accurate enough


def test_determine_temporal_background_correction(experiment, shared_datadir):
    movie = experiment.files[0].movie
    method = 'fit_background_peak'
    background_correction = movie.determine_temporal_background_correction(method)
    assert (shared_datadir / 'BN_TIRF' / 'TIRF 561 0001_corrections.nc').is_file()