import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve
from skimage.transform import AffineTransform

from trace_analysis.trace_extraction import make_gaussian_mask


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

    if plot:
        plt.figure()
        plt.imshow(pseudo_image_source, origin='lower')
        plt.figure()
        plt.imshow(pseudo_image_destination, origin='lower')

        plt.figure()
        plt.imshow(correlation)
        plt.show()

    def correlation_coordinates_to_translation_coordinates(correlation_peak_coordinates):
        # return back_conversion_destination(correlation_peak_coordinates - np.array(pseudo_image_source.shape)[::-1])
        transfomation_correlation = AffineTransform(translation=correlation_peak_coordinates - (np.array(pseudo_image_source.shape)[::-1]-1))
        transfomation_destination_inverse = AffineTransform(transfomation_destination._inv_matrix)
        return transfomation_source + transfomation_correlation + transfomation_destination_inverse

    return correlation, correlation_coordinates_to_translation_coordinates