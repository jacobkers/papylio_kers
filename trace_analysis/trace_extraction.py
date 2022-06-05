import PIL.ImageFilter
import numpy as np
import xarray as xr
from tqdm import tqdm
import dask_image.ndfilters
import scipy.ndimage

def make_gaussian_mask(size, offsets, sigma=1.291):
    # TODO: Explain calculation in docstring
    # It is to keep the photon number the same after applying the mask.
    # If there is a PSF of N photons, which is nothing but a 2D Gauss function with given sigma and amplitude,
    # the sum of the pixel is N. The idea is that the pixel sum should be the same after applying the mask.
    # The normalization factor is calculated to compensate the amplitude of 2D Gaussian after applying the mask.
    # The normalization factor should be different for different PSF size (i.e. different magnification or setup).
    # So N = sum(mask * (psf_single_photon*N)), and so sum(mask*psf_single_photon)
    # Both the mask and the psf are 2d Gaussians
    import xarray as xr
    roi = xr.DataArray(np.mgrid[0:size,0:size]-size//2, dims=('dimension','y','x'), coords={'dimension': ['y', 'x']})
    masks = np.exp(-((roi - offsets) ** 2).sum('dimension') / sigma**2 / 2).transpose('molecule','channel','y','x')
    psfs_single_photon = masks/masks.sum(dim=('x', 'y'))
    norm_factors = (masks*psfs_single_photon).sum(dim=('x','y'))
    masks = masks/norm_factors
    return masks

def extract_traces(movie, coordinates, background=None, mask_size=1.291, neighbourhood_size=11, correct_illumination=False):
    # go through all images, extract donor and acceptor signal
    coordinates['dimension'] = coordinates.dimension.astype('U')
    with movie:
        movie.read_header()

        intensity = xr.DataArray(np.empty((len(coordinates.molecule), len(coordinates.channel), movie.number_of_frames)),
                                 dims=['molecule', 'channel', 'frame'],
                                 coords=coordinates.drop('dimension').coords, name='intensity')

        # channel_offsets = xr.DataArray(np.vstack([channel.origin for channel in movie.channels]),
        #                                dims=('channel', 'dimension'),
        #                                coords={'channel': [channel.index for channel in movie.channels],
        #                                        'dimension': ['x', 'y']}) # TODO: Move to Movie
        # coordinates = coordinates - channel_offsets

        if background is None:
            background = xr.DataArray(dims=['molecule','channel'], coords={'molecule': coordinates.molecule, 'channel': coordinates.channel})

        offsets = coordinates % 1
        twoD_gaussians = make_gaussian_mask(size=neighbourhood_size, offsets=offsets, sigma=mask_size)

        coordinates_floored = (coordinates // 1).astype(int)

        roi_indices_general = xr.DataArray(np.mgrid[:neighbourhood_size, :neighbourhood_size] - neighbourhood_size // 2,
                                           dims=('dimension', 'y', 'x'),
                                           coords={'dimension': ['y', 'x']})#.transpose()

        roi_indices = coordinates_floored + roi_indices_general

        if correct_illumination:
            illumination_correction = IlluminationCorrection(movie.number_of_frames,
                                                             filter_function=scipy.ndimage.minimum_filter,
                                                             size=15, mode='wrap')

        oneD_indices = (roi_indices.sel(dimension='y')*movie.width+roi_indices.sel(dimension='x')).stack(peak=('molecule','channel')).stack(i=('y','x'))
        for frame_index in tqdm(range(movie.number_of_frames), desc=movie.name, leave=True):  # self.number_of_frames also works for pm, len(self.movie_file_object.filelist) not
            # print(frame_number)
            # if frame_number % 13 == 0:
            #     sys.stdout.write(f'\r   Frame {frame_number} of {movie.number_of_frames}')

            image = movie.read_frame(frame_index).astype('uint16')
            frame = xr.DataArray(image, dims=('y','x'))

            # TODO: Proper background subtraction

            if correct_illumination:
                illumination_correction.add_frame(frame_index, frame)
                # TODO: Determine how illumination correction is dependent on background

            if 'illumination' in background.dims:
                # TODO: Make this work properly and rename illumination in Movie
                # Do background subtraction on entire frame instead???
                frame_background = background.sel(illumination=movie.illumination.sel(frame=frame_index))
                # frame_background = background[frame_number % number_illumination]
            else:
                frame_background = background

            #intensity[:, :, frame_index] = extract_intensity_from_frame(frame, background, roi_indices, twoD_gaussians)
            intensity[:, :, frame_index] = extract_intensity_from_frame(frame, frame_background, oneD_indices, twoD_gaussians)

        # sys.stdout.write(f'\r   Frame {frame_number+1} of {movie.number_of_frames}\n')
        dataset = intensity.to_dataset()

        if correct_illumination:
            dataset['illumination_correction'] = illumination_correction.illumination_correction

    return dataset

# def extract_intensity_from_frame(frame, background, roi_indices, twoD_gaussians):
#     intensities = frame.sel(x=roi_indices.sel(dimension='x'), y=roi_indices.sel(dimension='y'))
#     intensities = intensities - background
#     weighted_intensities = intensities * twoD_gaussians
#     intensity_in_frame = weighted_intensities.sum(dim=('x', 'y'))
#     return intensity_in_frame

# A ufunc is probably better here
# def extract_intensity_from_frame(frame, background, roi_indices, twoD_gaussians):
#     intensities = frame.values[roi_indices.values[:,:,1,:,:], roi_indices.values[:,:,0,:,:]]
#     intensities = intensities - background.values[:,:,None,None]
#     weighted_intensities = intensities * twoD_gaussians.values
#     intensity_in_frame = weighted_intensities.sum(axis=(2,3))
#     return intensity_in_frame

def extract_intensity_from_frame(frame, background, oneD_indices, twoD_gaussians):  # extract traces
    intensities = frame.values.take(oneD_indices.values).reshape(twoD_gaussians.shape)
    intensities = intensities - background.values[:,:,None,None]
    weighted_intensities = intensities * twoD_gaussians.values
    intensity_in_frame = weighted_intensities.sum(axis=(2,3))
    return intensity_in_frame


class IlluminationCorrection:
    def __init__(self, number_of_frames, filter_function=scipy.ndimage.minimum_filter, **kwargs):
        self.filter_function = filter_function
        self.filter_kwargs = kwargs
        self._illumination_correction = np.empty(number_of_frames)

    def add_frame(self, index, frame):
        filtered_frame = self.filter_function(np.array(frame), **self.filter_kwargs)
        self._illumination_correction[index] = filtered_frame.sum()

    @property
    def illumination_correction(self):
        correction = self._illumination_correction.max() / self._illumination_correction
        return xr.DataArray(correction, dims=('frame',), name='illumination_correction')



# def illumination_intensity_from_frames(frames=None, filter_neighbourhood_size=15):
#     if frames is None:
#         frames = self.movie.read_frames_raw()
#
#     if not frames.chunks:
#         filtered_images = minimum_filter(frames, size=(1, 1, filter_neighbourhood_size, filter_neighbourhood_size))
#     else:
#         filtered_images = dask_image.ndfilters.minimum_filter(frames.data, size=(1, 1, filter_neighbourhood_size, filter_neighbourhood_size))
#
#     filtered_images = xr.DataArray(filtered_images, coords=frames.coords, name='illumination_correction')
#     illumination_intensity = filtered_images.sum(dim=('x','y'))
#     illumination_correction = (illumination_intensity.max(dim='frame') / illumination_intensity).T
#     # illumination_correction = illumination_correction.reset_index('frame', drop=True)
#     illumination_correction.to_netcdf(self.absoluteFilePath.with_suffix('.nc'), engine='h5netcdf', mode='a')