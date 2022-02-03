import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from sklearn.neighbors import NearestNeighbors
from skimage.transform import AffineTransform, PolynomialTransform
from scipy.optimize import minimize

from trace_analysis.plotting import scatter_coordinates, show_point_connections
from trace_analysis.mapping.polywarp import PolywarpTransform

def ComputeKDE(P, resolution, min_val, max_val):
    grids = np.round((max_val-min_val)/resolution).astype(int)+20
    KDE = np.zeros(grids)

    start = min_val - 10*resolution*np.ones(2)
    within_range = np.all([P[:, 0] >= min_val[0] - 6 * resolution, P[:, 0] <= max_val[0] + 6 * resolution,
                           P[:, 1] >= min_val[1] - 6 * resolution, P[:, 1] <= max_val[1] + 6 * resolution], axis=0)



    for point in P[within_range]:
        center = np.round((point-start)/resolution).astype(int)
        x_range = np.arange(center[0] - 3, center[0] + 3 + 1)
        y_range = np.arange(center[1] - 3, center[1] + 3 + 1)
        x_val = start[0] + x_range * resolution - point[0]
        y_val = start[1] + y_range * resolution - point[1]
        kernel_x = np.exp(- x_val * x_val / (resolution * resolution))
        kernel_x = kernel_x / np.sum(kernel_x)
        kernel_y = np.exp(- y_val * y_val / (resolution * resolution))
        kernel_y = kernel_y / np.sum(kernel_y)
        KDE[x_range[0]:x_range[-1]+1,y_range[0]:y_range[-1]+1] = KDE[x_range[0]:x_range[-1]+1,y_range[0]:y_range[-1]+1] + np.outer(kernel_x,kernel_y)

    nm = np.sqrt(np.sum(KDE**2))
    KDE = KDE/nm

    return KDE


def ComputeKC(param, Model, Scene, resolution, display_it, SceneKDE, min_val, max_val):
    tr = AffineTransform(matrix=np.hstack([param, [0,0,1]]).reshape(3,3))
    PT = tr(Model)
    MKDE = ComputeKDE(PT, resolution, min_val, max_val)
    KCVal = -np.sum(MKDE*SceneKDE)

    if callable(display_it):
        kwargs = {'iteration': 1, 'error':KCVal, 'X': PT, 'Y': Scene}
        callback(**kwargs)
    elif display_it==True:
        plt.figure()
        plt.scatter(*PT.T)
        plt.scatter(*Scene.T)
        plt.title(f'KC value: {-KCVal}')

    return KCVal

def KCReg(M, S, h, display=False, motion='affine'):
    min_val = S.min(axis=0)
    max_val = S.max(axis=0)
    SceneKDE = ComputeKDE(S, h, min_val, max_val)

    if callable(display):
        kwargs = {'iteration': 1, 'error':0, 'X': M, 'Y': S}
        callback(**kwargs)
    elif display == True:
        plt.figure()
        plt.scatter(*M.T)
        plt.scatter(*S.T)
        plt.title('initial setup')

    if motion=='affine':
        initial_transformation = AffineTransform()
        res = minimize(ComputeKC, initial_transformation.params.flatten()[0:6], args=(M, S, h, display, SceneKDE, min_val, max_val), tol=1e-6, options={'maxiter': 100})
    else:
        raise ValueError

    return AffineTransform(matrix=np.hstack([res.x, [0,0,1]]).reshape(3,3))

#
# from scipy.spatial import cKDtree
# from scipy.spatial import distance_matrix
# def ComputeKC(transformation_parameters, source, destination, plot=False, axis=None):
#     transformation = AffineTransform(matrix=np.hstack([transformation_parameters, [0,0,1]]).reshape(3,3))
#     source_transformed = transformation(source)
#     distances = distance_matrix(source_transformed, destination)
#     KCVal = np.exp(-distances**2)
#
#     if plot:
#         if axis is None:
#             axis = plt.gca()
#         axis.cla()
#         axis.scatter(*source_transformed.T, color='green', label='Source')
#         axis.scatter(*destination.T, color='red', label='Destination')
#         axis.set_title(f'KC value: {-KCVal}')
#         axis.draw()
#         plt.pause(0.001)
#
#     return KCVal
#
#
#
# def KCReg(source, destination, plot=False):
#     # SceneKDE = ComputeKDE(S, h, min_val, max_val)
#
#     initial_transformation = AffineTransform()
#     res = minimize(ComputeKC, initial_transformation.params.flatten()[0:6], args=(source, destination, plot), tol=1e-6, options={'maxiter': 100})
#
#     return AffineTransform(matrix=np.hstack([res.x, [0,0,1]]).reshape(3,3))
#




if __name__ == '__main__':
    from trace_analysis.plugins.sequencing.point_set_simulation import simulate_mapping_test_point_set

    # Simulate soure and destination point sets
    number_of_source_points = 400
    transformation = AffineTransform(translation=[5,0], rotation=2/360*2*np.pi, scale=[0.98, 0.98])
    source_bounds = ([0, 0], [256, 512])
    source_crop_bounds = None
    fraction_missing_source = 0.8
    fraction_missing_destination = 0.6
    maximum_error_source = 4
    maximum_error_destination = 4
    shuffle = True

    source, destination = simulate_mapping_test_point_set(number_of_source_points, transformation,
                                                          source_bounds, source_crop_bounds,
                                                          fraction_missing_source, fraction_missing_destination,
                                                          maximum_error_source, maximum_error_destination, shuffle)

    from trace_analysis.mapping.mapping import Mapping2
    im = Mapping2(source, destination)
    im.transformation = transformation
    im.show_mapping_transformation()

    def visualize(iteration, error, X, Y, ax):
        plt.cla()
        ax.scatter(X[:, 0], X[:, 1], color='red', label='Target')
        ax.scatter(Y[:, 0], Y[:, 1], color='blue', label='Source')
        plt.text(0.87, 0.92, 'Iteration: {:d}\nQ: {:06.4f}'.format(
            iteration, error), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes,
                 fontsize='x-large')
        ax.legend(loc='upper left', fontsize='x-large')
        plt.draw()
        plt.pause(0.001)


    from functools import partial
    fig = plt.figure()
    fig.add_axes([0, 0, 1, 1])
    callback = partial(visualize, ax=fig.axes[0])

    found_transformation = KCReg(source, destination, 5, display=callback) #callback)


    # # Perform icp on the simulated point sets
    # max_iterations = 20
    # tolerance = 0.0000001
    # cutoff = None
    # cutoff_final = 10
    # initial_transformation = None
    # transform = AffineTransform
    # transform_final = PolywarpTransform
    #
    # transformation, transformation_inverse, distances, i = icp(source, destination, max_iterations, tolerance,
    #                                                            cutoff, cutoff_final, initial_transformation,
    #                                                            transform, transform_final, show_plot=True)
