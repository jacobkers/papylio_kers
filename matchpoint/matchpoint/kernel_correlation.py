import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from sklearn.neighbors import NearestNeighbors
from skimage.transform import AffineTransform, PolynomialTransform, EuclideanTransform, SimilarityTransform, ProjectiveTransform
from scipy.optimize import minimize


# def ComputeKDE(P, resolution, min_val, max_val):
#     grids = np.round((max_val-min_val)/resolution).astype(int)+20
#     KDE = np.zeros(grids)
#
#     start = min_val - 10*resolution*np.ones(2)
#     within_range = np.all([P[:, 0] >= min_val[0] - 6 * resolution, P[:, 0] <= max_val[0] + 6 * resolution,
#                            P[:, 1] >= min_val[1] - 6 * resolution, P[:, 1] <= max_val[1] + 6 * resolution], axis=0)
#
#
#
#     for point in P[within_range]:
#         center = np.round((point-start)/resolution).astype(int)
#         x_range = np.arange(center[0] - 3, center[0] + 3 + 1)
#         y_range = np.arange(center[1] - 3, center[1] + 3 + 1)
#         x_val = start[0] + x_range * resolution - point[0]
#         y_val = start[1] + y_range * resolution - point[1]
#         kernel_x = np.exp(- x_val * x_val / (resolution * resolution))
#         kernel_x = kernel_x / np.sum(kernel_x)
#         kernel_y = np.exp(- y_val * y_val / (resolution * resolution))
#         kernel_y = kernel_y / np.sum(kernel_y)
#         KDE[x_range[0]:x_range[-1]+1,y_range[0]:y_range[-1]+1] = KDE[x_range[0]:x_range[-1]+1,y_range[0]:y_range[-1]+1] + np.outer(kernel_x,kernel_y)
#
#     nm = np.sqrt(np.sum(KDE**2))
#     KDE = KDE/nm
#
#     return KDE
#
#
# def ComputeKC(param, Model, Scene, resolution, display_it, SceneKDE, min_val, max_val):
#     tr = AffineTransform(matrix=np.hstack([param, [0,0,1]]).reshape(3,3))
#     PT = tr(Model)
#     MKDE = ComputeKDE(PT, resolution, min_val, max_val)
#     KCVal = -np.sum(MKDE*SceneKDE)
#
#     if callable(display_it):
#         kwargs = {'iteration': 1, 'error':KCVal, 'X': PT, 'Y': Scene}
#         callback(**kwargs)
#     elif display_it==True:
#         plt.figure()
#         plt.scatter(*PT.T)
#         plt.scatter(*Scene.T)
#         plt.title(f'KC value: {-KCVal}')
#
#     return KCVal
#
# def KCReg(M, S, h, display=False, motion='affine'):
#     min_val = S.min(axis=0)
#     max_val = S.max(axis=0)
#     SceneKDE = ComputeKDE(S, h, min_val, max_val)
#
#     if callable(display):
#         kwargs = {'iteration': 1, 'error':0, 'X': M, 'Y': S}
#         callback(**kwargs)
#     elif display == True:
#         plt.figure()
#         plt.scatter(*M.T)
#         plt.scatter(*S.T)
#         plt.title('initial setup')
#
#     if motion=='affine':
#         initial_transformation = AffineTransform()
#         res = minimize(ComputeKC, initial_transformation.params.flatten()[0:6], args=(M, S, h, display, SceneKDE, min_val, max_val), tol=1e-6, options={'maxiter': 100})
#     else:
#         raise ValueError
#
#     return AffineTransform(matrix=np.hstack([res.x, [0,0,1]]).reshape(3,3))

#
import time
#from scipy.spatial import cKDTree
from scipy.spatial import distance_matrix, cKDTree

def parameter_to_transformation(transformation_parameters):
    if len(transformation_parameters) == 2:
        transformation = EuclideanTransform(rotation=0,
                                            translation=transformation_parameters[1:2])
    elif len(transformation_parameters) == 3:
        transformation = EuclideanTransform(rotation=transformation_parameters[0],
                                            translation=transformation_parameters[1:3])
    elif len(transformation_parameters) == 4:
        transformation = SimilarityTransform(scale=transformation_parameters[0], rotation=transformation_parameters[1],
                                             translation=transformation_parameters[2:4])
    elif len(transformation_parameters) == 5:
        transformation = AffineTransform(scale=transformation_parameters[0], rotation=transformation_parameters[1],
                                         shear=transformation_parameters[2], translation=transformation_parameters[3:5])
    elif len(transformation_parameters) == 6:
        transformation = AffineTransform(scale=transformation_parameters[0:2], rotation=transformation_parameters[2],
                                         shear=transformation_parameters[3], translation=transformation_parameters[4:6])
    return transformation


def compute_kernel_correlation(transformation, source, destination, sigma=1, plot=False, axis=None, per_point_pair=False):

    # transformation = AffineTransform(matrix=np.hstack([transformation_parameters, [0,0,1]]).reshape(3,3))
    # t = []
    # t.append(time.time())
    if not isinstance(transformation, ProjectiveTransform):
        transformation = parameter_to_transformation(transformation)

    # t.append(time.time())
    # print(t[-1]-t[-2])

    source_transformed = transformation(source)
    # t.append(time.time())
    # print(t[-1]-t[-2])

    # distances = distance_matrix(source_transformed, destination.data, threshold=1e9)
    # t.append(time.time())
    # print(t[-1]-t[-2])
    #
    # KCVal = -np.exp(-distances**2/(4*sigma**2)).sum()
    # t.append(time.time())
    # print(t[-1]-t[-2])
    # print(KCVal)

    #destination_tree = cKDTree(destination)
    source_transformed_tree = cKDTree(source_transformed)
    if isinstance(destination, cKDTree):
        destination_tree = destination
    else:
        destination_tree = cKDTree(destination)

    #distances2 = destination_tree.sparse_distance_matrix(source_tree, 3 * sigma, output_type='dok_matrix').values()
    #distances2 = np.fromiter(destination_tree.sparse_distance_matrix(source_tree, 3 * sigma, output_type='dict').values(),
                # dtype=float)
    if not per_point_pair:
        distances = destination_tree.sparse_distance_matrix(source_transformed_tree, 5 * sigma, output_type='ndarray')['v']
    else:
        distances = destination_tree.sparse_distance_matrix(source_transformed_tree, 5 * sigma, output_type='dok_matrix').power(-1).toarray()**-1
    # t.append(time.time())
    # print(t[-1]-t[-2])

    KCVal = -np.exp(-distances ** 2 / (4 * sigma ** 2))

    if not per_point_pair:
        KCVal = KCVal.sum()
    # t.append(time.time())
    # print(t[-1]-t[-2])

    # print(KCVal)


    if plot:
        if axis is None:
            axis = plt.gca()
        axis.cla()
        axis.scatter(*source_transformed.T, color='green', label='Source')
        axis.scatter(*destination.T, color='red', label='Destination')
        axis.set_title(f'KC value: {-KCVal}')
        #axis.draw()
        plt.pause(0.001)

    return KCVal

from scipy.optimize import shgo, dual_annealing, differential_evolution, basinhopping
def kernel_correlation(source, destination, bounds, sigma=1, plot=False, **kwargs):
    # if transform is None:
    #     initial_transformation = SimilarityTransform()
    # parameters = initial_transformation.params.flatten()[0:6]
    # parameters = np.hstack([initial_transformation.scale, initial_transformation.rotation, initial_transformation.translation])

    translation_to_origin = AffineTransform(translation=-destination.mean(axis=0))
    source = translation_to_origin(source)
    destination = translation_to_origin(destination)

    destination_cKDTree = cKDTree(destination)

    # res = minimize(compute_kernel_correlation, parameters, args=(source, destination_cKDTree, sigma, plot), tol=1e-6,
    #                bounds=bounds, method='L-BFGS-B',
    #                options={'disp': None, 'maxcor': 20, 'ftol': 2.220446049250313e-09, 'gtol': 1e-05, 'eps': 1e-08,
    #                         'maxfun': 15000, 'maxiter': 15000, 'iprint': - 1, 'maxls': 20, 'finite_diff_rel_step': None})
    # res = minimize(compute_kernel_correlation, parameters, args=(source, destination_cKDTree, sigma, plot), tol=1e-6,
    #                bounds=bounds, method='Powell',
    #                options={'xtol': 0.00001, 'ftol': 0.000001, 'maxiter': None, 'maxfev': None,
    #                         'disp': False, 'direc': None, 'return_all': False})

    # res = minimize(compute_kernel_correlation, parameters, args=(source, destination_cKDTree, sigma, plot), tol=1e-9,
    #                bounds=bounds, method='BFGS',
    #                options={'gtol': 1e-09, 'norm': np.inf, 'eps': 1.4901161193847656e-03, 'maxiter': None, 'disp': False,
    #                         'return_all': False, 'finite_diff_rel_step': None})

    # res = shgo(compute_kernel_correlation, bounds, args=(source, destination_cKDTree, sigma, plot), constraints=None, n=None, iters=1,
    #            callback=None, minimizer_kwargs={'method': 'Powell', 'options': {'ftol': 1e-9, 'gtol': 1e-5}},
    #             options={'minimize_every_iter':False}, sampling_method='simplicial')
    #
    # res = dual_annealing(compute_kernel_correlation, bounds, args=(source, destination_cKDTree, sigma, plot),
    #                      maxiter=50,
    #                      local_search_options={}, initial_temp=5230,#5230.0,
    #                      restart_temp_ratio=2e-05, visit=2.62, accept=- 5.0, maxfun=10000000.0, seed=None,
    #                      no_local_search=False, callback=None, x0=None)

    result = differential_evolution(compute_kernel_correlation, bounds, args=(source, destination_cKDTree, sigma, plot), **kwargs)

    # res = basinhopping(compute_kernel_correlation, parameters, niter=100, T=1.0, stepsize=0.5, minimizer_kwargs={'args': (source, destination_cKDTree, sigma, plot)}, take_step=None,
    #                    accept_test=None, callback=None, interval=50, disp=False, niter_success=None, seed=None)


    #print(result)
    # return AffineTransform(matrix=np.hstack([res.x, [0,0,1]]).reshape(3,3))

    transformation = parameter_to_transformation(result.x)

    final_transformation = AffineTransform(
        matrix=(translation_to_origin._inv_matrix @ transformation.params @ translation_to_origin.params))

    return final_transformation, result


    # return parameter_to_transformation(result.x), result





if __name__ == '__main__':
    from point_set_simulation import simulate_mapping_test_point_set

    # Simulate source and destination point sets
    number_of_points = 10000
    transformation = AffineTransform(translation=[10,-10], rotation=1/360*2*np.pi, scale=[0.98, 0.98])
    bounds = ([0, 0], [256, 512])
    crop_bounds = (None, None)
    fraction_missing_source = 0.95
    fraction_missing_destination = 0.6
    maximum_error_source = 0.5
    maximum_error_destination = 0.5
    shuffle = True

    destination, source = simulate_mapping_test_point_set(number_of_points, transformation.inverse,
                                                          bounds, crop_bounds,
                                                          (fraction_missing_source, fraction_missing_destination),
                                                          (maximum_error_source, maximum_error_destination), shuffle)

    from matchpoint.core import MatchPoint
    im = MatchPoint(source, destination)
    im.transformation = transformation
    im.show_mapping_transformation()

    # def visualize(iteration, error, X, Y, ax):
    #     plt.cla()
    #     ax.scatter(X[:, 0], X[:, 1], color='red', label='Target')
    #     ax.scatter(Y[:, 0], Y[:, 1], color='blue', label='Source')
    #     plt.text(0.87, 0.92, 'Iteration: {:d}\nQ: {:06.4f}'.format(
    #         iteration, error), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes,
    #              fontsize='x-large')
    #     ax.legend(loc='upper left', fontsize='x-large')
    #     plt.draw()
    #     plt.pause(0.001)
    #
    #
    # from functools import partial
    # fig = plt.figure()
    # fig.add_axes([0, 0, 1, 1])
    # callback = partial(visualize, ax=fig.axes[0])
    #
    # found_transformation = KCReg(source, destination, 5, display=False) #callback)
    #plt.figure()
    minimization_bounds = ((0.97,1.02),(-0.05,0.05),(-20,20),(-20,20))
    found_transformation, result = kernel_correlation(source, destination, minimization_bounds, 1, plot=False,
                                 strategy='best1bin', maxiter=1000, popsize=50, tol=0.01,
                                 mutation=0.25, recombination=0.7, seed=None, callback=None, disp=False,
                                 polish=True, init='sobol', atol=0, updating='immediate', workers=1,
                                 constraints=()
                                 )
    fm = MatchPoint(source, destination)
    fm.transformation = found_transformation
    fm.show_mapping_transformation()
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


