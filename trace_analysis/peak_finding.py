import numpy as np
import matplotlib.pyplot as plt
import cv2
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import math

from trace_analysis.image_adapt.analyze_label import analyze # note analyze label is differently from the approach in pick spots

def find_peaks(image=None, method='AKAZE', **kwargs):
    if method == 'AKAZE':
        coordinates = analyze(image)[2]
    elif method == 'absolute-threshold':
        coordinates = find_peaks_absolute_threshold(image, **kwargs)
    elif method == 'adaptive-threshold':
        coordinates = find_peaks_adaptive_threshold(image, **kwargs)
    elif method == 'local-maximum':
        coordinates = find_peaks_local_maximum(image, **kwargs)
    return coordinates

def find_peaks_absolute_threshold(image, threshold = None, minimum_area = 5, maximum_area = 15):
    if threshold is None: threshold = (np.max(image) + np.min(image)) / 2
    image_thresholded = image > threshold
    coordinates = coordinates_from_contours(image_thresholded, minimum_area, maximum_area)
    return coordinates

def find_peaks_adaptive_threshold(image, minimum_area = 5, maximum_area = 15):
    # This may be needed if we go from a 16-bit image to an 8 bit image
    # if bounds is None:
    #     lower_bound = np.min(image)
    #     upper_bound = np.percentile(image.flatten(), 99.999)
    # else:
    #     lower_bound = bounds[0]
    #     upper_bound = bounds[1]
    #
    # # Change threshold and image to 8-bit, as cv2 can only analyse 8-bit images
    # # threshold = ((threshold - lower_bound) / (upper_bound - lower_bound) * 255)
    # # threshold = np.clip(threshold, 0, 255).astype('uint8')
    # image = ((image - lower_bound) / (upper_bound - lower_bound) * 255)
    # image = np.clip(image, 0, 255).astype('uint8')

    image_thresholded = cv2.adaptiveThreshold(image.astype(np.uint8),
                                              maxValue=1,
                                              adaptiveMethod=cv2.ADAPTIVE_THRESH_GAUSSIAN_C,
                                              thresholdType=cv2.THRESH_BINARY,
                                              blockSize=9,
                                              C=0)
    coordinates = coordinates_from_contours(image_thresholded, minimum_area, maximum_area)

    return coordinates

def find_peaks_local_maximum(image, threshold = 25, threshold_max = math.inf, filter_neighbourhood_size = 10):
    image_max = filters.maximum_filter(image, filter_neighbourhood_size)
    maxima = (image == image_max)
    image_min = filters.minimum_filter(image, filter_neighbourhood_size)
    # Probably I need to make the neighbourhood_size of the minimum filter larger.

    diff_threshold = ((image_max - image_min) > threshold)
    diff_threshold_max = ((image_max - image_min) < threshold_max)
    diff=diff_threshold*diff_threshold_max
    maxima[diff == 0] = 0

    labeled, num_objects = ndimage.label(maxima)
    if num_objects > 0:
        coordinates = np.fliplr(np.array(ndimage.center_of_mass(image, labeled, range(1, num_objects + 1))))
    else:
        coordinates = np.array([])
        print('No peaks found')

    return coordinates

def coordinates_from_contours(image_thresholded, minimum_area=5, maximum_area=15):
    contours, hierarchy = cv2.findContours(image_thresholded.astype(np.uint8),
                                           cv2.RETR_EXTERNAL,
                                           cv2.CHAIN_APPROX_NONE)
    x = []
    y = []

    # colorImg = cv2.cvtColor(image, cv2.COLOR_GRAY2BGR)

    coordinates = []

    for c in contours:
        # Calculate moments for each contour
        M = cv2.moments(c)

        # Calculate x, y coordinate of center

        if M["m00"] != 0:
            cX = int(M["m10"] / M["m00"])
            cY = int(M["m01"] / M["m00"])

            area = M["m00"]  # Is the same as cv2.contourArea(c) # Is the same as M["m00"]

            if (area > minimum_area) and (area < maximum_area):
                x = np.append(x, cX)
                y = np.append(y, cY)
                coordinates.append(np.array([cX, cY]))
        else:
            cX, cY = 0, 0

        # cv2.circle(colorImg, (cX, cY), 8, (0, 0, 255), thickness=1)

    return np.array(coordinates)




if __name__ == '__main__':
    from trace_analysis import Experiment
    exp = Experiment(r'D:\ivoseverins\SURFdrive\Promotie\Code\Python\traceAnalysis\twoColourExampleData\20141017 - Holliday junction - Copy\HJC-50pM')
    file = exp.files[2]
    movie = file.movie
    image = movie.make_average_image(write=False)

    coordinates = find_peaks(image=image, method='adaptive-threshold', minimum_area=5, maximum_area=15)

    plt.imshow(image)
    plt.scatter(coordinates[:,0],coordinates[:,1],color='r')

