import numpy as np
import matplotlib.pyplot as plt

import itertools
from scipy.spatial import cKDTree
import random
from trace_analysis.mapping.point_set import crop_coordinates
from coordinate_transformations import translate, rotate, magnify, reflect, transform
# from trace_analysis.plotting import scatter_coordinates
from skimage.transform import AffineTransform, SimilarityTransform
import time


class GeometricHashTable:
    def __init__(self, destinations, source_vertices=None, initial_source_transformation=AffineTransform(),
                 number_of_source_bases=20, number_of_destination_bases='all',
                 tuple_size=4, maximum_distance_source=None, maximum_distance_destination=None):
        self.initial_source_transformation = initial_source_transformation

        self.tuple_size = tuple_size
        self.maximum_distance_source = maximum_distance_source
        self.maximum_distance_destination = maximum_distance_destination

        # self.mode = mode
        # if mode == 'translation': self.hashTableRange = [-10000, 10000]
        # else: self.hashTableRange = [-1,1]
        # self.nBins = nBins

        #
        # self.number_of_source_bases = number_of_source_bases
        # self.number_of_destination_bases = number_of_destination_bases

        # self._hashTable = None
        # self._matches = None
        #
        # self.rotationRange = rotationRange
        # self.magnificationRange = magnificationRange
        #

        self.destinations = destinations
        self.destination_KDTrees = [cKDTree(destination) for destination in destinations]
        self.source_vertices = source_vertices

        self.number_of_hashtable_entries_per_destination = []
        self.number_of_hashtable_entries_per_basis = []
        self.number_of_bases_per_destination = []

        self.create_hashtable()

    def create_hashtable(self):
        destination_hash_data = geometric_hash(self.destinations, self.maximum_distance_destination, self.tuple_size)
        self.destination_KDTrees, self.destination_tuple_sets, self.destination_hash_table_KDTree, \
        self.destination_transformation_matrices = destination_hash_data

    def query(self, source, distance=15, alpha=0.9, sigma=10, K_threshold=10e9, hash_table_distance_threshold=0.01,
              magnification_range=None, rotation_range=None):

        return find_match_after_hashing(source, self.maximum_distance_source, self.tuple_size, self.source_vertices,
                                        self.destination_KDTrees, self.destination_tuple_sets,
                                        self.destination_hash_table_KDTree,
                                        hash_table_distance_threshold, alpha, sigma, K_threshold,
                                        magnification_range, rotation_range)

    def query_tuple_transformations(self, sources, hash_table_distance_threshold=0.01, parameters=['rotation', 'scale'],
                                    plot=False, method='dbscan', **method_kwargs):
        # np.vstack(sources)

        source_hash_data = geometric_hash(sources, self.maximum_distance_source, self.tuple_size)
        _, _, source_hash_table_KDTree, source_transformation_matrices = source_hash_data

        return compare_tuple_transformations(source_hash_table_KDTree, source_transformation_matrices,
                                              self.destination_hash_table_KDTree, self.destination_transformation_matrices,
                                              hash_table_distance_threshold, parameters, plot=plot, method=method, **method_kwargs)


def geometric_hash(point_sets, maximum_distance=100, tuple_size=4):
    # TODO: Add minimum_distance and implement
    # TODO: Make invariant to mirroring
    # TODO: Make usable with multiple point-sets in a single hash table
    # TODO: Implement names of point_sets, possibly through a dictionary and adding a attribute to each KDtree
    start_time = time.time()

    if type(point_sets) is not list:
        point_sets = [point_sets]

    # TODO: do this only if point_sets are not cKDTrees already, pass cKDtree instead of point_set ?
    point_set_KDTrees = [cKDTree(point_set) for point_set in point_sets]

    hash_tables = []
    point_tuple_sets = []
    for point_set_KDTree in point_set_KDTrees:
        point_tuples = generate_point_tuples(point_set_KDTree, maximum_distance, tuple_size)
        hash_table, transformation_matrices = geometric_hash_table(point_set_KDTree, point_tuples)

        hash_tables.append(hash_table)
        point_tuple_sets.append(point_tuples)

    hash_table_KDTree = cKDTree(np.vstack(hash_tables))

    # print("--- %s seconds ---" % (time.time() - start_time))

    return point_set_KDTrees, point_tuple_sets, hash_table_KDTree, transformation_matrices


def generate_point_tuples(point_set_KDTree, maximum_distance, tuple_size):
    # TODO: If maximum distance is None, calculate maximum distance (i.e. containing all points, or giving a good density
    if maximum_distance is None:
        # Smallest enclosing circle would be better: https://rosettacode.org/wiki/Smallest_enclosing_circle_problem
        maximum_distance = np.sqrt(((point_set_KDTree.maxes-point_set_KDTree.mins)**2).sum())
    pairs = list(point_set_KDTree.query_pairs(maximum_distance))

    point_tuples = []

    # TODO: improve speed / different method
    for pair in pairs:
        pair_coordinates = point_set_KDTree.data[list(pair)]
        center = (pair_coordinates[0] + pair_coordinates[1]) / 2
        distance = np.linalg.norm(pair_coordinates[0] - pair_coordinates[1])
        internal_points = point_set_KDTree.query_ball_point(center, distance / 2 * 0.99)
        for internal_point_combination in itertools.combinations(internal_points, tuple_size-2):
            point_tuples.append(pair + tuple(internal_point_combination))
            # yield pair + tuple(internal_point_combination)

    return point_tuples


def geometric_hash_table(point_set_KDTree, point_tuples):
    pt = np.array(point_tuples)
    d = point_set_KDTree.data[pt]
    d0 = np.array([[0,0],[1,1]])


    T = np.repeat(np.expand_dims(np.diag([1.,1.,1.],k=0), axis=0), len(d), axis=0)
    T[:, :2, 2] = d0[0] - d[:, 0, :]


    diffs_d = d[:,1,:] - d[:,0,:]
    diffs_d0 = d0[1] - d0[0]
    m_d = np.linalg.norm(diffs_d, axis=1)
    m_d0 = np.linalg.norm(diffs_d0)

    unit_diffs_d = np.divide(diffs_d, m_d[:, None])
    unit_diffs_d0 = diffs_d0 / m_d0

    D = np.repeat(np.expand_dims(np.diag([1., 1., 1.], k=0), axis=0), len(d), axis=0)
    D[:,0,0] = D[:,1,1] = 1 / m_d * m_d0

    dot_product = np.dot(unit_diffs_d, unit_diffs_d0) # for mapping to (0,0),(1,1) could be replaced by sum to increase speed
    cross_product = np.cross(unit_diffs_d, unit_diffs_d0) # for mapping to (0,0),(1,1) could be replaced by subtraction to increase speed
    R = np.repeat(np.expand_dims(np.diag([1., 1., 1.], k=0), axis=0), len(d), axis=0)
    R[:,0,0] = R[:,1,1] = dot_product
    R[:,1,0] = cross_product
    R[:,0,1] = -cross_product

    transformation_matrices = R@D@T

    # transformation_matrices @ np.dstack([d, np.ones((d.shape[0], d.shape[1], 1))]).swapaxes(1,2)
    hash_code = transformation_matrices[:,:2,:] @ np.dstack([d[:,len(d0):,:], np.ones((d.shape[0], d.shape[1]-len(d0), 1))]).swapaxes(1, 2)

    # Break similarity of pair
    break_similarity = hash_code.swapaxes(1,2)[:, :, 0].sum(axis=1) > ((pt.shape[1] - d0.shape[0]) / 2)
    # this is not yet suited for using more than two basis points
    switch_basis_points_matrix = np.array([[-1,0,1],[0,-1,1],[0,0,1]])
    transformation_matrices[break_similarity] = switch_basis_points_matrix[None,:,:]@transformation_matrices[break_similarity]
    hash_code[break_similarity] = transformation_matrices[break_similarity, :2, :] @ np.dstack(
        [d[break_similarity, len(d0):, :], np.ones((d[break_similarity].shape[0], d.shape[1] - len(d0), 1))]).swapaxes(1, 2)

    hash_code = hash_code.swapaxes(1, 2)

    # Break similarity of internal points
    s = hash_code[:, :, 0].argsort()
    hash_code = np.take_along_axis(hash_code, s[:, :, None], axis=1)

    return hash_code.reshape(len(hash_code), -1), transformation_matrices


def find_match_after_hashing(source, maximum_distance_source, tuple_size, source_vertices,
                             destination_KDTrees, destination_tuple_sets, destination_hash_table_KDTree,
                             hash_table_distance_threshold=0.01,
                             alpha=0.1, sigma=10, K_threshold=10e9,
                             scaling_range=None, rotation_range=None):

    if type(destination_KDTrees) is not list:
        destination_KDTrees = [destination_KDTrees]

    source_KDTree = cKDTree(source)
    source_tuples = generate_point_tuples(source_KDTree, maximum_distance_source, tuple_size)

    hash_table_distances_checked = 0
    tuples_checked = 0
    # for source_tuple_index in np.arange(len(source_tuples)):
    source_hash_table, source_transformation_matrices = geometric_hash_table(source_KDTree, source_tuples)
    for source_tuple, source_hash_code in zip(source_tuples, source_hash_table):
        # TODO: Get all destination tuple indices within a range
        distance, destination_tuple_index = destination_hash_table_KDTree.query(source_hash_code)

        hash_table_distances_checked += 1
        # if hash_table_distances_checked > 500:
        #     return

        if distance < hash_table_distance_threshold:
            # Find the destination tuple by determining in which destination pointset the destination tuple index is located
            # and by determining what the index is within that destination pointset.
            cumulative_tuples_per_destination = np.cumsum([0]+[len(point_tuples) for point_tuples in destination_tuple_sets])
            destination_index = np.where((cumulative_tuples_per_destination[:-1] <= destination_tuple_index) &
                                         (cumulative_tuples_per_destination[1:] > destination_tuple_index))[0][0]
            tuple_index_in_destination_set = destination_tuple_index - cumulative_tuples_per_destination[destination_index]
            destination_tuple = destination_tuple_sets[destination_index][tuple_index_in_destination_set]

            #Or list(itertools.chain.from_iterable(destination_tuple_sets))[destination_tuple_index]

            tuples_checked += 1

            destination_KDTree = destination_KDTrees[destination_index]
            found_transformation = tuple_match(source, destination_KDTree, source_tuple, destination_tuple)

            if not transformation_in_parameter_range(found_transformation, rotation=rotation_range, scale=scaling_range):
                continue

            source_indices_without_tuple = [i for i in range(len(source)) if i not in source_tuple]
            source_without_tuple = source[source_indices_without_tuple]

            if test_transformation(source_without_tuple, destination_KDTree, found_transformation, source_vertices,
                                alpha=alpha, sigma=sigma, K_threshold=K_threshold):
                return found_transformation

            #     match = Mapping2(source=source, destination=destination_KDTree.data, method='Geometric hashing',
            #                      transformation_type='linear', initial_transformation=None)
            #     match.transformation = found_transformation
            #     match.destination_index = destination_index
            #
            #     match.hash_table_distance = distance
            #     match.hash_table_distances_checked = hash_table_distances_checked
            #     match.tuples_checked = tuples_checked
            #     return match


def polygon_area(vertices):
    x = vertices[:,0]
    y = vertices[:,1]
    return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))


def tuple_match(source, destination_KDTree, source_tuple, destination_tuple):
    source_coordinate_tuple = source[list(source_tuple)]
    destination_coordinate_tuple = destination_KDTree.data[list(destination_tuple)]

    transformation = estimate_transformation(source_coordinate_tuple[:2], destination_coordinate_tuple[:2])

    # transformation = AffineTransform()
    # transformation.estimate(source_coordinate_tuple[:2], destination_coordinate_tuple[:2])

    # transformation = AffineTransform(transformation_matrix)
    return transformation


def estimate_transformation(start_points, end_points, tr=None, di=None, ro=None):
    start_points = np.atleast_2d(start_points)
    end_points = np.atleast_2d(end_points)
    if len(start_points) == 1 & len(end_points) == 1:
        tr = True
        ro = False
        di = False

    elif len(start_points) == 2 & len(end_points) == 2:
        if tr is None: tr = True
        if di is None: di = True
        if ro is None: ro = True

    transformation_matrix = np.identity(3)

    if tr:
        translation_matrix = translate(end_points[0] - start_points[0])
        transformation_matrix = translation_matrix @ transformation_matrix

    if di or ro:
        diffs = np.array([start_points[0] - start_points[1], end_points[0] - end_points[1]])
        diff_lengths = np.linalg.norm(diffs, axis=1, keepdims=True)
        unit_diffs = diffs / diff_lengths

        if di:
            dilation_matrix = magnify(diff_lengths[1] / diff_lengths[0], end_points[0])
            transformation_matrix = dilation_matrix @ transformation_matrix

        if ro:
            angle = -np.arctan2(np.linalg.det(unit_diffs), np.dot(unit_diffs[0], unit_diffs[1]))
            # angle = np.arccos(np.dot(diffs[0]/endLength,diffs[1]/startLength))
            rotation_matrix = rotate(angle, end_points[0])
            transformation_matrix = rotation_matrix @ transformation_matrix

    return SimilarityTransform(transformation_matrix)


def apply_transformation(point_set, transformation_matrix):
    # ~3x Faster than Affinetransform()(point_set)
    point_set = np.append(point_set, np.ones((point_set.shape[0], 1)), axis=1)
    transformed_point_set = (transformation_matrix @ point_set.T)[0:2, :].T
    return transformed_point_set


def transformation_in_parameter_range(transformation, **kwargs):
    for parameter, parameter_range in kwargs.items():
        if parameter_range is not None:
            parameter_value = getattr(transformation, parameter)
            if not (parameter_range[0] < np.array(parameter_value) < parameter_range[1]):
                return False
    return True


def test_transformation(source, destination, found_transformation, source_vertices, alpha=0.9, sigma=10, K_threshold=10e2):
    if not type(destination) is cKDTree:
        destination_KDTree = cKDTree(destination)
    else:
        destination_KDTree = destination

    # TODO: Make crop both source_transformed and destination based on convex hull (and thus make this independen of source_vertices)
    source_vertices_transformed = found_transformation(source_vertices)
    destination_cropped = crop_coordinates(destination_KDTree.data, source_vertices_transformed)

    source_transformed_area = polygon_area(source_vertices_transformed)

    pDB = 1 / source_transformed_area

    K=1

    # TODO: Remove basis points?
    for coordinate in found_transformation(source):
        distance, index = destination_KDTree.query(coordinate)

        # 2d Gaussian
        pDF = alpha / source_transformed_area + \
              (1 - alpha) / (2 * np.pi * sigma ** 2) * \
              np.exp(-(distance ** 2) / (2 * sigma ** 2)) / len(destination_cropped)

        K = K * pDF/pDB
        if K > K_threshold:
            return True

    return False


def compare_tuple_transformations(source_hash_table_KDTree, source_transformation_matrices, destination_hash_table_KDTree,
                                  destination_transformation_matrices, hash_table_distance_threshold=0.01,
                                  parameters=['rotation', 'scale'], plot=False, method='dbscan', **method_kwargs):
    tuple_matches = source_hash_table_KDTree.query_ball_tree(destination_hash_table_KDTree, hash_table_distance_threshold)

    # TODO: make this matrix multiplication
    transformation_matrices = []
    for source_index, destination_indices in enumerate(tuple_matches):
        source_transformation_matrix = source_transformation_matrices[source_index]
        # source_transformation_matrix = np.linalg.inv(source_transformation_matrix)
        for destination_index in destination_indices:
            destination_transformation_matrix = destination_transformation_matrices[destination_index]
            destination_transformation_matrix_inverse = np.linalg.inv(destination_transformation_matrix)
            transformation_matrices.append(destination_transformation_matrix_inverse @ source_transformation_matrix)
    transformation_matrices = np.stack(transformation_matrices)

    # plt.figure()
    # plt.hist(transformation_matrices[:, 0, 2], 100)

    transformations = [AffineTransform(transformation_matrix) for transformation_matrix in transformation_matrices]


    parameter_values = [np.array([getattr(t, parameter) for t in transformations]).T for parameter in parameters]

    sample = np.vstack(list(parameter_values)).T
    # sample = transformation_matrices[:, :2, :].reshape(-1, 6)

    if method == 'histogram_max':
        h, edges = np.histogramdd(sample, **method_kwargs)

        bin_centers = [(e[:-1]+e[1:])/2 for e in edges]
        # found_transformation = np.array([bc[h_index] for bc, h_index in zip(bin_centers, np.where(h==h.max()))]).reshape(2,3)
        # t = AffineTransform(np.vstack([found_transformation, [0, 0, 1]]))
        hist_max_index = np.where(h == h.max())
        if len(hist_max_index[0]) > 1:
            raise RuntimeError('No optimal transformation found')
        found_values = [bc[h_index][0] for bc, h_index in zip(bin_centers, hist_max_index)]


    elif method == 'dbscan':
        import sklearn.cluster

        sample_normalized = (sample - sample.min(axis=0)) / (sample.max(axis=0) - sample.min(axis=0))
        clustering = sklearn.cluster.DBSCAN(eps=0.01, min_samples=10, **method_kwargs).fit(sample_normalized)
        found_values = list(np.median(sample[clustering.labels_ == 0], axis=0))
        if clustering.labels_.max() > 0:
            raise RuntimeError('No optimal transformation found')
        if plot:
            plt.figure()
            plt.scatter(*sample[:, :2].T, c=clustering.labels_)
            plt.show()
    else:
        raise ValueError(f'Unknown method {method}')

    parameter_dict = {
        parameter: found_values.pop(0) if parameter == 'rotation' else [found_values.pop(0), found_values.pop(0)]
        for parameter in parameters}



    # found_transformation = AffineTransform(**parameter_dict)

    # Hc = H.copy()
    # Hc[i]=0
    # print(np.max(H)/np.max(Hc)>2)

    # x, y = bin_centers_x, bin_centers_y
    # X, Y = np.meshgrid(x, y)
    # xdata = np.vstack((X.ravel(), Y.ravel()))

    # def twoD_gaussian(M, offset, amplitude, x0, y0, sigma_x, sigma_y):
    #     x, y = M
    #     return offset + amplitude * np.exp(- ((x - x0) / (2 * sigma_x)) ** 2 - ((y - y0) / (2 * sigma_y)) ** 2)
    #
    # from scipy.optimize import curve_fit
    # p0 = [20,20,0,0,1,1]
    # popt, pcov = curve_fit(twoD_gaussian, xdata, H.ravel())#, p0) #input: function, xdata, ydata,p0
    #
    # coeff, var_matrix = curve_fit(gauss.mult_gaussFun_Fit, (bin_centers_x, bin_centers_y), H, p0=p0)

    # plt.figure()
    # plt.scatter(ms, rs)
    #
    # plt.figure()
    # plt.hist(ms, 100)
    #
    # plt.figure()
    # plt.hist(rs, 100)

    # return parameter_dict
    return AffineTransform(**parameter_dict)


