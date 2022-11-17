import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve
from skimage.transform import AffineTransform



def make_gaussian_mask(size, center=None, offset=(0, 0), sigma=1.291):
    # TODO: Change this to a gaussian only
    # This is actually not what we are looking for here, we just need a gaussian

    x = np.arange(0, size, 1, float)
    y = x[:, np.newaxis]

    if center is None: center = [size // 2, size //2]

    mask = np.exp(-((x - center[0] - offset[0]) ** 2 + (y - center[0] - offset[1]) ** 2) / sigma**2 / 2)
    psf_single_photon = mask/np.sum(mask)
    norm_factor = np.sum(np.multiply(mask, psf_single_photon))
    mask = np.divide(mask, norm_factor)

    return mask

def coordinates_to_image(coordinates, gaussian_width=7, divider=5):
    gauss = make_gaussian_mask(gaussian_width)

    min_x, min_y = coordinates.min(axis=0)


    transformation = AffineTransform(translation=[-min_x, -min_y]) + AffineTransform(scale=1/divider)
    coordinates = transformation(coordinates)

    max_x, max_y = coordinates.max(axis=0)

    image_width = int(np.ceil(max_x)) + 1
    image_height = int(np.ceil(max_y)) + 1

    image = np.zeros((image_height, image_width))
    indices = coordinates.round().astype(int)
    image[indices[:,1], indices[:,0]] = 1

    image_with_gaussians = fftconvolve(image, gauss)

    # def image_to_original_coordinates(image_coordinates):
    #     return image_coordinates+[[min_x, min_y]]

    return image_with_gaussians, transformation

def cross_correlate(source, destination, gaussian_width=7, divider=5, plot=False):
    pseudo_image_source, transfomation_source = coordinates_to_image(source, gaussian_width=gaussian_width, divider=divider) #/ 5)
    pseudo_image_destination, transfomation_destination = coordinates_to_image(destination, gaussian_width=gaussian_width, divider=divider) #/ 5)

    from scipy.signal import correlate, correlation_lags

    correlation = correlate(pseudo_image_destination, pseudo_image_source, mode='full')

    if plot is not False:
        if type(plot) is bool:
            axes = []
            for i in range(3):
                figure, axis = plt.subplots()
                axes.append(axis)
        else:
            axes = plot

        bounds_source = transfomation_source.inverse(np.array([[0, 0], pseudo_image_source.shape[::-1]])).T
        axes[0].imshow(pseudo_image_source, origin='lower', extent=bounds_source.flatten())
        bounds_destination = transfomation_destination.inverse(np.array([[0, 0], pseudo_image_destination.shape[::-1]])).T
        axes[1].imshow(pseudo_image_destination, origin='lower', extent=bounds_destination.flatten())

        axes[2].imshow(correlation, origin='lower')
        # plt.show()

    def correlation_coordinates_to_translation_coordinates(correlation_peak_coordinates):
        # return back_conversion_destination(correlation_peak_coordinates - np.array(pseudo_image_source.shape)[::-1])
        transfomation_correlation = AffineTransform(translation=correlation_peak_coordinates - (np.array(pseudo_image_source.shape)[::-1]-1))
        transfomation_destination_inverse = AffineTransform(transfomation_destination._inv_matrix)
        return transfomation_source + transfomation_correlation + transfomation_destination_inverse

    return correlation, correlation_coordinates_to_translation_coordinates