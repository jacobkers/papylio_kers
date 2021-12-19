import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve

from trace_analysis.trace_extraction import make_gaussian_mask


def coordinates_to_image(coordinates, gaussian_width=7):
    gauss = make_gaussian_mask(gaussian_width)

    min_x, min_y = coordinates.min(axis=0)
    max_x, max_y = coordinates.max(axis=0)

    coordinates = coordinates-[[min_x, min_y]]

    image_width = int(np.ceil(max_x)-np.floor(min_x))+1
    image_height = int(np.ceil(max_y)-np.floor(min_y))+1

    image = np.zeros((image_height, image_width))
    indices = coordinates.round().astype(int)
    image[indices[:,1], indices[:,0]] = 1

    image_with_gaussians = fftconvolve(image, gauss)

    def image_to_original_coordinates(image_coordinates):
        return image_coordinates+[[min_x, min_y]]

    return image_with_gaussians, image_to_original_coordinates

def cross_correlate(source, destination):
    pseudo_image_source, back_conversion_source = coordinates_to_image(source) #/ 5)
    pseudo_image_destination, back_conversion_destination = coordinates_to_image(destination) #/ 5)

    plt.figure()
    plt.imshow(pseudo_image_source, origin='lower')
    plt.figure()
    plt.imshow(pseudo_image_destination, origin='lower')

    from scipy.signal import correlate, correlation_lags

    correlation = correlate(pseudo_image_destination, pseudo_image_source, mode='full')
    plt.figure()
    plt.imshow(correlation)
    plt.show()

    def correlation_coordinates_to_translation_coordinates(correlation_peak_coordinates):
        return back_conversion_destination(correlation_peak_coordinates - np.array(pseudo_image_source.shape)[::-1])

    return correlation, correlation_coordinates_to_translation_coordinates