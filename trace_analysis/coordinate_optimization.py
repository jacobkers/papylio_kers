import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def coordinates_within_margin(coordinates,  image = None, bounds = None, margin=10):
    if coordinates.size == 0: return np.array([])

    if image is not None:
        bounds = np.array([[0,image.shape[1]], [0,image.shape[0]]])

    criteria = np.array([(coordinates[:, 0] > (bounds[0,0] + margin)),
                         (coordinates[:, 0] < (bounds[0,1] - margin)),
                         (coordinates[:, 1] > (bounds[1,0] + margin)),
                         (coordinates[:, 1] < (bounds[1,1] - margin))
                         ])

    return coordinates[criteria.all(axis=0)]

def circle(r):
    d = 2*r + 1
    rx, ry = d/2, d/2
    x, y = np.indices((d, d))
    return (np.abs(np.hypot(r - x, r - y)-r) < 0.5).astype(int)

def coordinates_without_intensity_at_radius(coordinates, image, radius, cutoff, fraction_of_peak_max = 0.25):
    if cutoff == 'image_median': cutoff = np.median(image)
    circle_matrix = circle(radius)
    new_coordinates = []

    coordinates = coordinates_within_margin(coordinates, image=image, margin=radius+1)

    for i, coordinate in enumerate(coordinates):
        # Could use the coordinates_within_margin for this [IS 01-11-2019]
        #if np.all(coordinate > radius+1) and \
        #        np.all(coordinate < (np.array(image.shape)-radius-1)):

        cropped_peak = crop(image, coordinate, radius*2+1)
        if np.all((cropped_peak * circle_matrix) <
                  (cutoff + fraction_of_peak_max * np.max(cropped_peak))):
            new_coordinates.append(coordinate)

    return np.array(new_coordinates)

def crop(image, center, width):
    center = np.round(center).astype(int)
    return image[(center[1]-width//2):(center[1]+width//2+1),(center[0]-width//2):(center[0]+width//2+1)]

def twoD_gaussian(M, offset, amplitude, x0, y0, sigma_x, sigma_y):
    x, y = M
    return offset + amplitude * np.exp(- ((x-x0)/(2*sigma_x))**2 - ((y-y0)/(2*sigma_y))**2)

def fit_twoD_gaussian(Z):
    height, width = Z.shape
    x, y = np.arange(width)-width//2, np.arange(height)-height//2
    X, Y = np.meshgrid(x, y)
    xdata = np.vstack((X.ravel(), Y.ravel()))

    p0 = [20,20,0,0,1,1]
    popt, pcov = curve_fit(twoD_gaussian, xdata, Z.ravel(), p0) #input: function, xdata, ydata,p0
      
    # The offset can potentially be used for background subtraction
    return popt

def coordinates_after_gaussian_fit(coordinates, image, gaussian_width = 9):
    new_coordinates = []
    coordinates = coordinates_within_margin(coordinates, image=image, margin=gaussian_width//2+1)

    for i, coordinate in enumerate(coordinates):
        # Could use the coordinates_within_margin for this [IS 01-11-2019]
        #if np.all(coordinate > gaussian_width//2+1) and \
        #        np.all(coordinate < np.array(image.shape)-gaussian_width//2-1):
        cropped_peak = crop(image, coordinate, gaussian_width)
        try:
            coefficients = fit_twoD_gaussian(cropped_peak)
            #new_coordinates.append(coordinate + coefficients[2:4])
            tmp=coordinate + coefficients[2:4]
            if np.sum(np.abs(coefficients[2:4]))<gaussian_width*2:
                new_coordinates.append(tmp)
            # else: #MD do nothing, you don't want to include fits with a center far outside the cropped image
        except RuntimeError:
            pass
       #     # MD: leave this print out? 
       #    print('RuntimeError: No 2dGaussian fit possible')
    return np.array(new_coordinates)

if __name__ == '__main__':
    from trace_analysis import Experiment
    from trace_analysis.peak_finding import find_peaks

    exp = Experiment(r'D:\ivoseverins\SURFdrive\Promotie\Code\Python\traceAnalysis\twoColourExampleData\20141017 - Holliday junction - Copy\HJC-50pM')
    file = exp.files[2]
    movie = file.movie
    image = movie.make_average_image(write=False)

    coordinates = find_peaks(image=image, method='adaptive-threshold', minimum_area=5, maximum_area=15)

    plt.imshow(image)
    plt.scatter(coordinates[:, 0], coordinates[:, 1], color='b')

    coordinates = coordinates_within_margin(coordinates, bounds = np.array([[0,255], [0,511]]), margin=10)
    plt.scatter(coordinates[:, 0], coordinates[:, 1], color='g')

    coordinates = coordinates_after_gaussian_fit(coordinates, image)
    plt.scatter(coordinates[:, 0], coordinates[:, 1], color='y')

    coordinates = coordinates_without_intensity_at_radius(coordinates,
                                                          image,
                                                          radius=4,
                                                          cutoff=np.median(image),
                                                          fraction_of_peak_max=0.35) # was 0.25 in IDL code
    plt.scatter(coordinates[:, 0], coordinates[:, 1], color='r')
