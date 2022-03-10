import numpy as np
import xarray as xr
from tqdm import tqdm


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


def extract_traces(movie, coordinates, background=None, mask_size=1.291, neighbourhood_size=11):
    # go through all images, extract donor and acceptor signal

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


        oneD_indices = (roi_indices.sel(dimension='y')*movie.width+roi_indices.sel(dimension='x')).stack(peak=('molecule','channel')).stack(i=('y','x'))
        for frame_number in tqdm(range(movie.number_of_frames), desc=movie.name, leave=True):  # self.number_of_frames also works for pm, len(self.movie_file_object.filelist) not
            # print(frame_number)
            # if frame_number % 13 == 0:
            #     sys.stdout.write(f'\r   Frame {frame_number} of {movie.number_of_frames}')

            image = movie.read_frame(frame_number)
            frame = xr.DataArray(image, dims=('y','x'))

            #intensity[:, :, frame_number] = extract_intensity_from_frame(frame, background, roi_indices, twoD_gaussians)
            intensity[:, :, frame_number] = extract_intensity_from_frame(frame, background, oneD_indices, twoD_gaussians)

        # sys.stdout.write(f'\r   Frame {frame_number+1} of {movie.number_of_frames}\n')

    return intensity

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