import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve
from skimage.transform import AffineTransform
from scipy.signal import correlate


def gaussian_kernel(size, center=None, sigma=1.291):
    x = np.arange(0, size, 1, float)+0.5
    y = x[:, np.newaxis]

    if center is None:
        center = [size / 2, size / 2]

    return np.exp(-((x - center[0]) ** 2 + (y - center[1]) ** 2) / (2 * sigma ** 2))


def coordinates_to_image(coordinates, kernel_size=7, gaussian_sigma=1, divider=5):
    gauss = gaussian_kernel(kernel_size, sigma=gaussian_sigma)

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


def cross_correlate(source, destination, kernel_size=7, gaussian_sigma=1, divider=5, subtract_background=True, plot=False, axes=None):
    pseudo_image_source, transformation_source = \
        coordinates_to_image(source, kernel_size=kernel_size, gaussian_sigma=gaussian_sigma, divider=divider)
    pseudo_image_destination, transformation_destination = \
        coordinates_to_image(destination, kernel_size=kernel_size, gaussian_sigma=gaussian_sigma, divider=divider)

    correlation_raw = correlate(pseudo_image_destination, pseudo_image_source, mode='full')

    if subtract_background:
        import scipy.ndimage.filters as filters
        correlation = correlation_raw - filters.minimum_filter(correlation_raw, 2 * kernel_size)
    else:
        correlation = correlation_raw
    # np.min(correlation.shape) / 200)

    if plot:
        if axes is None:
            axes = []
            for i in range(4):
                figure, axis = plt.subplots()
                axes.append(axis)

        bounds_source = transformation_source.inverse(np.array([[0, 0], pseudo_image_source.shape[::-1]])).T
        axes[0].imshow(pseudo_image_source, origin='lower', extent=bounds_source.flatten())
        bounds_destination = transformation_destination.inverse(np.array([[0, 0], pseudo_image_destination.shape[::-1]])).T
        axes[1].imshow(pseudo_image_destination, origin='lower', extent=bounds_destination.flatten())
        bounds_correlation = np.array([[0, 0], np.array(correlation.shape[::-1])*divider]).T
        bounds_correlation -= np.array([pseudo_image_source.shape[::-1]]).T*divider
        axes[2].imshow(correlation_raw, origin='lower', extent=bounds_correlation.flatten())
        if len(axes) > 3:
            axes[3].imshow(correlation, origin='lower', extent=bounds_correlation.flatten())

    def correlation_coordinates_to_translation_coordinates(correlation_peak_coordinates):
        # return back_conversion_destination(correlation_peak_coordinates - np.array(pseudo_image_source.shape)[::-1])
        transformation_correlation = AffineTransform(translation=correlation_peak_coordinates - (np.array(pseudo_image_source.shape)[::-1]-1))
        transformation_destination_inverse = AffineTransform(transformation_destination._inv_matrix)
        return transformation_source + transformation_correlation + transformation_destination_inverse

    return correlation, correlation_coordinates_to_translation_coordinates
