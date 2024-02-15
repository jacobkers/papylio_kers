import numpy as np
import matplotlib.pyplot as plt
import tqdm
from scipy.signal import fftconvolve
from skimage.transform import AffineTransform
from scipy.signal import correlate


def gaussian_kernel(size, center=None, sigma=1.291):
    x = np.arange(0, size, 1, float)+0.5
    y = x[:, np.newaxis]

    if center is None:
        center = [size / 2, size / 2]

    return np.exp(-((x - center[0]) ** 2 + (y - center[1]) ** 2) / (2 * sigma ** 2))

def correlate_normalized(in1, in2, normalize=True, zero_outer_pixels=3):
    # As adapted from equation 2 in Fast Normalized Cross-Correlation by J. P. Lewis
    shape_full = np.array(in1.shape)+np.array(in2.shape)-1
    out_full = np.zeros(shape_full)

    for v in tqdm.tqdm(range(shape_full[0])):
        for u in range(shape_full[1]):
            in1_overlap = in1[max(in1.shape[0]-v-1,0):min(shape_full[0]-v,in1.shape[0]),
                              max(in1.shape[1]-u-1,0):min(shape_full[1]-u,in1.shape[1])]
            in2_overlap = in2[max(0, v+1-in1.shape[0]):min(v+1, in2.shape[0]),
                              max(0, u+1-in1.shape[1]):min(u+1, in2.shape[1])]

            if normalize:
                in1_overlap = in1_overlap - in1_overlap.mean()
                in2_overlap = in2_overlap - in2_overlap.mean()
                nominator = np.sum(in1_overlap * in2_overlap)
                denominator = np.sqrt(np.sum(in1_overlap**2) * np.sum(in2_overlap**2))
                out_full[v, u] = nominator / denominator

            else:
                out_full[v,u] = np.sum(in1_overlap*in2_overlap)

    out_full = out_full[::-1, ::-1]
    out_full = np.nan_to_num(out_full)

    if zero_outer_pixels:
        out_full[:zero_outer_pixels, :] = 0
        out_full[-zero_outer_pixels:, :] = 0
        out_full[:, :zero_outer_pixels] = 0
        out_full[:, -zero_outer_pixels:] = 0

    return out_full


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

def cross_correlate(source, destination, kernel_size=7, gaussian_sigma=1, divider=5, subtract_background=True,
                    normalize=False, plot=False, axes=None):
    pseudo_image_source, transformation_source = \
        coordinates_to_image(source, kernel_size=kernel_size, gaussian_sigma=gaussian_sigma, divider=divider)
    pseudo_image_destination, transformation_destination = \
        coordinates_to_image(destination, kernel_size=kernel_size, gaussian_sigma=gaussian_sigma, divider=divider)

    if normalize:
        correlation_raw = correlate_normalized(pseudo_image_destination, pseudo_image_source)
    else:
        correlation_raw = correlate(pseudo_image_destination, pseudo_image_source, mode='full')

    if subtract_background == 'minimum_filter':
        import scipy.ndimage.filters as filters
        correlation = correlation_raw - filters.minimum_filter(correlation_raw, 2 * kernel_size)
    elif subtract_background == 'median_filter':
        import scipy.ndimage.filters as filters
        correlation = correlation_raw - filters.median_filter(correlation_raw, 2 * kernel_size)
    elif (subtract_background == 'expected_signal') or (subtract_background == True):
        ones_source = np.ones_like(pseudo_image_source)
        ones_destination = np.ones_like(pseudo_image_destination)
        normalize_source = correlate(pseudo_image_source, ones_destination, mode='full')
        normalize_destination = correlate(ones_source, pseudo_image_destination, mode='full')
        correlate_ones = correlate(ones_source, ones_destination, mode='full')
        correlation = correlation_raw - (normalize_source*normalize_destination/correlate_ones)
    elif subtract_background == 'expected_signal_rough':
        ones_source = np.ones_like(pseudo_image_source)
        ones_destination = np.ones_like(pseudo_image_destination)
        correlate_ones = correlate(ones_source, ones_destination, mode='full')
        correlation = correlation_raw - correlate_ones * pseudo_image_source.mean() * pseudo_image_destination.mean()
    else:
        correlation = correlation_raw

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
