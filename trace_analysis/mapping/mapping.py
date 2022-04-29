# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 13:55:53 2019

@author: ivoseverins
"""
import json

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.collections import PatchCollection
from pathlib import Path
import yaml
import skimage.transform
from skimage.transform import AffineTransform, PolynomialTransform, SimilarityTransform
# import matplotlib.path as pth
from shapely.geometry import Polygon, MultiPoint, LineString, Point
from tqdm import tqdm
from scipy.spatial import distance_matrix, cKDTree

from icp import icp, nearest_neighbor_pair, nearest_neighbour_match, direct_match
from polywarp import PolywarpTransform
from polynomial import PolynomialTransform
from point_set_simulation import simulate_mapping_test_point_set
from kernel_correlation import kernel_correlation, compute_kernel_correlation
from cross_correlation import cross_correlate

class Mapping2:
    """Mapping class to find, improve, store and use the mapping between a source point set and a destination point set

    Attributes
    ----------

    source_name : str
        Name of the source point set
    source : Nx2 numpy.ndarray
        Coordinates of the source point set
    destination_name : str
        Name of the destination point set
    destination : Nx2 numpy.ndarray
        Coordinates of the destination point set
    method : str
        Method for finding the transformation between source and destination. Options: 'direct', 'nearest_neighbour'
        and 'iterative_closest_point'
    transformation_type : str
        Type of transformation used, linear or polynomial
    initial_transformation : AffineTransform or PolynomialTransform
        Initial transformation to perform as a starting point for the mapping algorithms
    transformation : skimage.transform.AffineTransform or skimage.transform.PolynomialTransform
        Transformation from source point set to destination point set
    transformation_inverse : skimage.transform.AffineTransform or skimage.transform.PolynomialTransform
        Inverse transformation, i.e. from destination point set to source point set

    """

    @classmethod
    def simulate(cls, number_of_points=200, transformation=None,
                 bounds=([0, 0], [256, 512]), crop_bounds=(None, None), fraction_missing=(0.1, 0.1),
                 error_sigma=(0.5, 0.5), shuffle=True):

        if transformation is None:
            transformation = SimilarityTransform(translation=[256, 10], rotation=1/360*2*np.pi, scale=[0.98, 0.98])


        source, destination = simulate_mapping_test_point_set(number_of_points, transformation,
                                                              bounds, crop_bounds, fraction_missing,
                                                              error_sigma,
                                                              shuffle)

        mapping = cls(source, destination)
        mapping.transformation_correct = transformation

        def show_correct_mapping_transformation(self, *args, **kwargs):
            transformation_temp = self.transformation
            self.transformation = self.transformation_correct
            self.show_mapping_transformation(*args, **kwargs)
            self.transformation = transformation_temp

        mapping.show_correct_mapping_transformation = show_correct_mapping_transformation.__get__(mapping)
        mapping.show_correct_mapping_transformation()

        return mapping

    @classmethod
    def load(cls, filepath):
        # TODO: Move the load part in __init__ to here
        return cls(load=filepath)


    transformation_types = {'linear': AffineTransform,
                            'nonlinear': PolywarpTransform,
                            'polynomial': PolynomialTransform,
                            'affine': AffineTransform,
                            'similarity': SimilarityTransform}

    def __init__(self, source=None, destination=None, method=None,
                 transformation_type='affine', initial_transformation=None, transformation=None, transformation_inverse=None,
                 source_name='source', destination_name='destination',
                 source_unit=None, destination_unit=None, destination_distance_threshold=0,
                 name=None, load=None):
        """Set passed object attributes

        Parameters
        ----------
        source : Nx2 numpy.ndarray
            Source coordinates
        destination : Nx2 numpy.ndarray
            Destination coordinates
        method : str
            Method for finding the transformation between source and destination. Options: 'direct', 'nearest_neighbour'
            and 'iterative_closest_point'
        transformation_type : str
            Type of transformation used, linear or polynomial
        initial_transformation : skimage.transform._geometric.GeometricTransform or dict
            Initial transformation used as starting point by the perform_mapping method.
            Can be given as dictionary specifying translation, scale, rotation and shear.
        load : str or pathlib.Path
            Path to file that has to be loaded
        """

        self.name = name
        self.label = ''
        self.source_name = source_name
        self.source = np.array(source) #source=donor=left side image
        self.source_unit = source_unit
        # self.source_distance_threshold = source_distance_threshold
        self.destination_distance_threshold = destination_distance_threshold
        self._source_vertices = None
        self.destination_name = destination_name
        self.destination_unit = destination_unit
        self.destination = np.array(destination) #destination=acceptor=right side image
        self._destination_vertices = None
        self.method = method
        self.matched_pairs = np.empty((0,2), dtype=int)
        self.save_path = None

        self.transformation_type = transformation_type

        if type(initial_transformation) is dict:
            initial_transformation = AffineTransform(**initial_transformation)
        self.initial_transformation = initial_transformation

        if transformation is None:
            self.transformation = self.transform()
            self.transformation_inverse = self.transform()
        else:
            if type(transformation) is dict:
                transformation = self.transform(**transformation)
            self.transformation = transformation

            if transformation_inverse is None:
                self.transformation_inverse = self.transform(matrix=self.transformation._inv_matrix)
            else:
                if type(transformation_inverse) is dict:
                    transformation_inverse = self.transform(**transformation_inverse)
                self.transformation_inverse = transformation_inverse


        # if (source is not None) and (destination is not None):
        #     if self.method is None: self.method = 'icp'
        #     self.perform_mapping()

        if load:
            filepath = Path(load)
            if filepath.suffix in ['.yml','.yaml','.json','.mapping']:
                if filepath.suffix in ['.json', '.mapping'] and filepath.open('r').read(1) == '{':
                    with filepath.open('r') as json_file:
                        attributes = json.load(json_file)
                elif filepath.suffix in ['.yml', '.yaml', '.mapping']:
                    with filepath.open('r') as yml_file:
                        attributes = yaml.load(yml_file, Loader=yaml.CSafeLoader)

                for key, value in attributes.items():
                    if type(value) == list:
                        value = np.array(value)
                    try:
                        setattr(self, key, value)
                    except AttributeError:
                        pass
                self.transform = self.transformation_types[self.transformation_type]
                self.transformation = self.transform(self.transformation)
                self.transformation_inverse = self.transform(self.transformation_inverse)

            self.name = filepath.with_suffix('').name
            self.save_path = filepath.parent

    # Function to make attributes from transformation available from the Mapping2 class
    def __getattr__(self, item):
        if ('transformation' in self.__dict__) and hasattr(self.transformation, item):
            return getattr(self.transformation, item)
        else:
            super().__getattribute__(item)

    # @property
    # def translation(self):
    #     return self.transformation[0:2,2]
    #
    # @property
    # def magnification(self):
    #     return np.linalg.norm(self.transformation[:,0:2],axis=0)
    #
    # @property
    # def rotation(self):
    #     # rotation_matrix = self.transformation[0:2, 0:2]/self.magnification[0]
    #     # return np.arctan2(rotation_matrix[0,1]-rotation_matrix[1,0],rotation_matrix[0,0]+rotation_matrix[1,1])/(2*np.pi)*360
    #     return np.arctan2(self.transformation[0, 1], self.transformation[0, 0]) / (2 * np.pi) * 360
    #
    # @property
    # def reflection(self):
    #     return np.array([np.sign(self.transformation[0, 0]), np.sign(self.transformation[1, 1])])

    @property
    def file_path(self):
        return self.save_path.joinpath(self.name).with_suffix('.mapping')

    @property
    def source_to_destination(self):
        """Nx2 numpy.ndarray : Source coordinates transformed to the destination axis"""

        return self.transform_coordinates(self.source)

    @property
    def destination_to_source(self):
        """Nx2 numpy.ndarray : Source coordinates transformed to the destination axis"""

        return self.transform_coordinates(self.destination, inverse=True)

    @property
    def source_vertices(self):
        if self._source_vertices is None:
            return determine_vertices(self.source, self.source_distance_threshold)
        else:
            return self._source_vertices

    @source_vertices.setter
    def source_vertices(self, vertices):
        self._source_vertices = vertices

    @property
    def destination_vertices(self):
        if self._source_vertices is None:
            return determine_vertices(self.destination, self.destination_distance_threshold)
        else:
            return self._destination_vertices

    @destination_vertices.setter
    def destination_vertices(self, vertices):
        self._destination_vertices = vertices

    @property
    def source_cropped_vertices(self): # or crop_vertices_in_source
        return overlap_vertices(self.source_vertices, self.transform_coordinates(self.destination_vertices, inverse=True))

    @property
    def source_cropped(self):
        return crop_coordinates(self.source, self.source_cropped_vertices)

    @property
    def destination_cropped_vertices(self): # or crop_vertices_in_destination
        return overlap_vertices(self.transform_coordinates(self.source_vertices), self.destination_vertices)

    @property
    def destination_cropped(self):
        return crop_coordinates(self.destination, self.destination_cropped_vertices)

    def get_source(self, crop=False, space='source'):
        if crop in ['destination', False]:
            source = self.source
        elif crop in ['source', True]:
            source = self.source_cropped

        if space in ['destination', self.destination_name]:
            source = self.transform_coordinates(source)

        return source

    def get_destination(self, crop=False, space='destination'):
        if crop in ['source', False]:
            destination = self.destination
        elif crop in ['destination', True]:
            destination = self.destination_cropped

        if space in ['source', self.source_name]:
            destination = self.transform_coordinates(destination, inverse=True)

        return destination

    def get_source_vertices(self, crop=False, space='source'):
        if crop in ['destination', False]:
            source_vertices = self.source_vertices
        elif crop in ['source', True]:
            source_vertices = self.source_cropped_vertices

        if space in ['destination', self.destination_name]:
            source_vertices = self.transformation(source_vertices)

        return source_vertices

    def get_destination_vertices(self, crop=False, space='source'):
        if crop in ['destination', False]:
            destination_vertices = self.destination_vertices
        elif crop in ['source', True]:
            destination_vertices = self.destination_cropped_vertices

        if space in ['source', self.source_name]:
            destination_vertices = self.transformation(destination_vertices, inverse=True)

        return destination_vertices

    @property
    def source_distance_threshold(self):
        return self.destination_distance_threshold / np.max(self.transformation.scale)

    def find_distance_threshold(self, method='single_match_optimization', **kwargs):
        if method == 'single_match_optimization':
            self.single_match_optimization(**kwargs)
        else:
            raise ValueError('Unknown method')

    def number_of_single_matches_for_radii(self, radii):
        distance_matrix_ = self.distance_matrix(crop=True)
        number_of_pairs = np.array([len(singly_matched_pairs_within_radius(distance_matrix_, r)) for r in radii])
        return number_of_pairs

    def single_match_optimization(self, maximum_radius=20, number_of_steps=100, plot=True):
        radii = np.linspace(0, maximum_radius, number_of_steps)
        number_of_pairs = self.number_of_single_matches_for_radii(radii)
        self.destination_distance_threshold = distance_threshold_from_number_of_matches(radii, number_of_pairs, plot=plot)

    def determine_matched_pairs(self, distance_threshold=None):
        #TODO: add crop
        if distance_threshold is None:
            distance_threshold = self.destination_distance_threshold

        distance_matrix_ = self.distance_matrix(crop=False)

        self.matched_pairs = singly_matched_pairs_within_radius(distance_matrix_, distance_threshold)

    def pair_coordinates(self, point_set='destination', space='destination'):
        if point_set == 'source':
            return self.get_source(space=space)[self.matched_pairs[:,0]]
        elif point_set == 'destination':
            return self.get_destination(space=space)[self.matched_pairs[:,1]]

    def pair_distances(self, space='destination', show=False, **kwargs):
        xy_distances = self.pair_coordinates(point_set='source', space=space) - \
                       self.pair_coordinates(point_set='destination', space=space)

        pair_distances = np.linalg.norm(xy_distances, axis=1)

        if show:
            figure, axis = plt.subplots()
            axis.hist(pair_distances, bins=100, **kwargs)
            axis.set_xlabel('Distance'+self.get_unit_label(space))
            axis.set_ylabel('Count')
            axis.set_title('Pair distances')

        return pair_distances

    @property
    def number_of_source_points(self):
        return self.source.shape[0]

    @property
    def number_of_destination_points(self):
        return self.destination.shape[0]

    @property
    def number_of_matched_points(self):
        """Number of matched points determined by finding the two-way nearest neigbours that are closer than a distance
        threshold.

        Parameters
        ----------
        distance_threshold : float
            Distance threshold for nearest neighbour match, i.e. only nearest neighbours with a distance smaller than
            the distance threshold are used.

        Returns
        -------
        int
            Number of matched points

        """
        # distances, source_indices, destination_indices = \
        #     nearest_neighbor_pair(self.source_to_destination, self.destination)
        #
        # return np.sum(distances < self.destination_distance_threshold)

        return self.matched_pairs.shape[0]

    def fraction_of_source_matched(self, crop=False):
        return self.number_of_matched_points / self.get_source(crop).shape[0]

    def fraction_of_destination_matched(self, crop=False):
        return self.number_of_matched_points / self.get_destination(crop).shape[0]

        # Possiblility to estimate area per point without source or destination vertices
        # from scipy.spatial import ConvexHull, convex_hull_plot_2d
        # hull = ConvexHull(points)
        # number_of_vertices = n = hull.vertices.shape[0]
        # number_of_points = hull.points.shape[0]
        # corrected_number_of_points_in_hull = number_of_points-number_of_vertices/2-1
        # area = hull.volume / corrected_number_of_points_in_hull * number_of_points
        # # Number of vertices = nv
        # # Sum of vertice angles = n*360
        # # Sum of inner vertice angles = (nv-2)*180
        # # Part of point area inside hull = (nv-2)*180/(nv*360)=(nv-2)/(2nv)
        # # Points inside the hull = (nv-2)/(2nv)*nv+np-nv=nv/2-1+np-nv=np-nv/2-1

        # Other method would be using voronoi diagrams, and calculating area of inner points

    @property
    def source_area(self):
        return area(self.source_vertices)

    @property
    def source_in_destination_area(self):
        return area(self.transformation(self.source_vertices))

    @property
    def source_cropped_area(self):
        return area(self.source_cropped_vertices)

    @property
    def destination_area(self):
        return area(self.destination_vertices)

    @property
    def destination_cropped_area(self):
        return area(self.destination_cropped_vertices)

    def get_source_area(self, crop=False, space='destination'):
        return area(self.get_source_vertices(crop=crop, space=space))

    def get_destination_area(self, crop=False, space='destination'):
        return area(self.get_destination_vertices(crop=crop, space=space))

    @property
    def transformation_type(self):
        return self._transformation_type

    @transformation_type.setter
    def transformation_type(self, value):
        self._transformation_type = value
        self.transform = self.transformation_types[self._transformation_type]

    # @property
    # def transform(self):
    #     """type : Transform class based on the set transformation_type"""
    #
    #     return self.transformation_types[self.transformation_type]


    # @property
    # def source_vertices(self):
    #     point_set = self.source
    #     from scipy.spatial import ConvexHull, convex_hull_plot_2d
    #     hull = ConvexHull(point_set)
    #     # number_of_vertices = n = hull.vertices.shape[0]
    #     # number_of_points = hull.points.shape[0]
    #     # corrected_number_of_points_in_hull = number_of_points-number_of_vertices/2-1
    #     # area = hull.volume / corrected_number_of_points_in_hull * number_of_points

    def perform_mapping(self, method=None, **kwargs):
        """Find transformation from source to destination points using one of the mapping methods

        The starting point for the mapping is the initial_transformation attribute.

        Parameters
        ----------
        method : str
            Mapping method, if not specified the object method is used.
        **kwargs
            Keyword arguments passed to mapping methods.
        """

        if method is None:
            method = self.method

        self.transformation = self.initial_transformation

        if method in ['icp', 'iterative_closest_point']: #icp should be default
            self.iterative_closest_point(**kwargs)
        elif method in ['direct']:
            self.direct_match()
        elif method in ['nn', 'nearest_neighbour']:
            self.nearest_neighbour_match(**kwargs)
        else:
            raise ValueError('Method not found')

        self.method = method

        self.show_mapping_transformation()

    def direct_match(self, transformation_type=None, **kwargs):
        """Find transformation from source to destination points by matching based on the point order

        Note: the number and the order of source points should be equal to the number and the order of destination points.

        Parameters
        ----------
        transformation_type : str
            Type of transformation used, either linear or polynomial can be chosen.
            If not specified the object transformation_type is used.
        **kwargs
            Keyword arguments passed to the direct match function.

        """
        if transformation_type is not None:
            self.transformation_type = transformation_type

        self.transformation, self.transformation_inverse, error = \
            direct_match(self.source, self.destination, transform=self.transform, return_inverse=True, **kwargs)

        print(f'Direct match\n'
              f'Mean-squared error: {error}')

    def nearest_neighbour_match(self, distance_threshold=None, transformation_type=None, **kwargs):
        """Find transformation from source to destination points by matching nearest neighbours

        Two-way nearest neighbours are detected, i.e. the source point should be the nearest neighbour of the
        destination point and vice versa. Only nearest neighbours closer than the distance threshold are used to find
        the transformation.

        Note
        ----
        The current transformation is used as starting point for the algorithm.

        Note
        ----
        The printed error is based on the points selected for matching.

        Parameters
        ----------
        distance_threshold : float
            Distance threshold for nearest neighbour match, i.e. only nearest neighbours with a distance smaller than
            the distance threshold are used.
        transformation_type : str
            Type of transformation used, either linear or polynomial can be chosen.
            If not specified the object transformation_type is used.
        **kwargs
            Keyword arguments passed to the nearest-neighbour match function.

        """
        if transformation_type is not None:
            self.transformation_type = transformation_type

        self.transformation, self.transformation_inverse, _, _, error = \
            nearest_neighbour_match(self.source, self.destination, transform=self.transform,
                                    initial_transformation=self.transformation, distance_threshold=distance_threshold,
                                    return_inverse=True, **kwargs)

        print(f'Nearest-neighbour match\n'
              f'Mean-squared error: {error}')

    def iterative_closest_point(self, distance_threshold=None, **kwargs):
        """Find transformation from source to destination points using an iterative closest point algorithm

        In the iterative closest point algorithm, the two-way nearest neigbhbours are found and these are used to
        find the most optimal transformation. Subsequently the source is transformed according to this
        transformation. This process is repeated until the changes detected are below a tolerance level.

        The iterative closest point algorithm can be used in situations when deviations between the two point sets
        are relatively small.

        Note
        ----
        The current transformation is used as starting point for the algorithm.

        Note
        ----
        The printed error is based on the points selected for matching.

        Parameters
        ----------
        distance_threshold : int or float
            Distance threshold applied to select nearest-neighbours in the final round of icp,
            i.e. nearest-neighbours with di.
        **kwargs
            Keyword arguments passed to icp.

        """

        self.transformation, self.transformation_inverse, error, number_of_iterations = \
            icp(self.source, self.destination, distance_threshold_final=distance_threshold,
                initial_transformation=self.transformation, transform_final=self.transform, **kwargs)

        print(f'Iterative closest point match\n'
              f'Mean-squared error: {error}\n'
              f'Number of iterations: {number_of_iterations}')

    def kernel_correlation(self, bounds=((0.97, 1.02), (-0.05, 0.05), (-10, 10), (-10, 10)), sigma=1, crop=False, **kwargs):
        transformation, result = kernel_correlation(self.get_source(crop, space='destination'), self.get_destination(crop),
                                                 bounds, sigma, plot=False, **kwargs)
        self.transformation = AffineTransform(matrix=(self.transformation + transformation).params)
        self.transformation_inverse = type(self.transformation)(self.transformation._inv_matrix)
        self.mapping_statistics = {'kernel_correlation_value': result.fun}

    def kernel_correlation_score(self, sigma=1, crop=False):
        return compute_kernel_correlation(self.transformation, self.get_source(crop), self.get_destination(crop), sigma=sigma)

    def cross_correlation(self, peak_detection='auto', gaussian_width=7, divider=5, crop=False, plot=False):

        if self.transformation is None:
            self.transformation = AffineTransform()
            self.transformation_inverse = AffineTransform()
        correlation, self.correlation_conversion_function = cross_correlate(self.get_source(crop, 'destination'), self.get_destination(crop),
                                                                            gaussian_width=gaussian_width, divider=divider, plot=plot)

        import scipy.ndimage.filters as filters
        corrected_correlation = correlation - filters.minimum_filter(correlation, 2*gaussian_width)#np.min(correlation.shape) / 200)
        if plot:
            plt.figure()
            plt.imshow(corrected_correlation)
            plt.show()

        if peak_detection == 'auto':
            correlation_peak_coordinates = np.array(np.where(corrected_correlation==corrected_correlation.max())).flatten()[::-1]
            plt.gca().scatter(correlation_peak_coordinates[0], correlation_peak_coordinates[1], marker='o', facecolors='none', edgecolors='r')
            self.set_correlation_peak_coordinates(correlation_peak_coordinates)

    def set_correlation_peak_coordinates(self, correlation_peak_coordinates):
        if not hasattr(self, 'correlation_conversion_function') and self.correlation_conversion_function is not None:
            raise RuntimeError('Run cross_correlation first')
        transformation = self.correlation_conversion_function(correlation_peak_coordinates) # is this the correct direction
        self.transformation = AffineTransform(matrix=(self.transformation + transformation).params)
        self.transformation_inverse = AffineTransform(matrix=self.transformation._inv_matrix)
        self.correlation_conversion_function = None

    def get_unit(self, space):
        return self.__getattribute__(f'{space}_unit')

    def get_unit_label(self, space):
        if self.get_unit(space) is not None:
            return f' ({self.get_unit(space)})'
        else:
            return ''

    def show_mapping_transformation(self, figure=None, show_source=False, show_destination=False, show_pairs=True,
                                    crop=False, inverse=False, source_colour='forestgreen', destination_colour='r',
                                    pair_colour='b', use_distance_threshold=True, save=False, save_path=None):
        """Show a point scatter of the source transformed to the destination points and the destination.

        Parameters
        ----------
        figure : matplotlib.figure.Figure
            If figure is passed then the scatter plot will be made in the currently active axis.
            If no figure is passed a new figure will be made.
        show_source : bool
            If True then the source points are also displayed

        """
        # Perhaps in case there is no transformation defined, we can show the source and destination, without transformed coordinates.
        if not figure:
            figure = plt.figure()

        source = self.get_source(crop)
        destination = self.get_destination(crop)

        if len(self.matched_pairs) == 0:
            show_pairs = False

        axis = figure.gca()

        if not inverse:
            all_transformed_coordinates = self.transform_coordinates(self.source)
            transformed_coordinates = self.transform_coordinates(source)
            transformed_coordinates_name = self.source_name
            transformed_coordinates_colour = source_colour
            distance_threshold = self.destination_distance_threshold
            show_destination = True
        else:
            all_transformed_coordinates = self.transform_coordinates(self.destination)
            transformed_coordinates = self.transform_coordinates(destination)
            transformed_coordinates_name = self.destination_name
            transformed_coordinates_colour = destination_colour
            distance_threshold = self.source_distance_threshold
            show_source = True

        if distance_threshold > 0 and use_distance_threshold:
            plot_circles(axis, transformed_coordinates, radius=distance_threshold, linewidth=1,
                         facecolor='none', edgecolor=transformed_coordinates_colour)
            if show_pairs:
                plot_circles(axis, all_transformed_coordinates[self.matched_pairs[:,0]],
                             radius=distance_threshold, linewidth=1,
                             facecolor='none', edgecolor=pair_colour)
        else:
            axis.scatter(*transformed_coordinates.T, facecolors='none', edgecolors=transformed_coordinates_colour,
                         linewidth=1, marker='o', label=f'{transformed_coordinates_name} transformed ({transformed_coordinates.shape[0]})')
            if show_pairs:
                axis.scatter(*all_transformed_coordinates[self.matched_pairs[:, 0]].T, facecolors='none',
                             edgecolors=pair_colour, linewidth=1, marker='o')

        if show_source:
            axis.scatter(*source.T, facecolors=source_colour, edgecolors='none', marker='.',
                         label=f'{self.source_name} ({source.shape[0]})')
            if show_pairs:
                axis.scatter(*self.source[self.matched_pairs[:,0]].T, facecolors=pair_colour, edgecolors='none', marker='.')

        if show_destination:
            axis.scatter(*destination.T, facecolors=destination_colour, edgecolors='none', marker='.',
                         label=f'{self.destination_name} ({destination.shape[0]})')
            if show_pairs:
                axis.scatter(*self.destination[self.matched_pairs[:,1]].T, facecolors=pair_colour, edgecolors='none', marker='.')

        axis.set_aspect('equal')

        if show_source and not show_destination:
            unit_label = self.get_unit_label('source')
        elif not show_source and show_destination:
            unit_label = self.get_unit_label('destination')
        elif show_source and show_destination:
            unit_label = ''

        axis.set_title(self.name)
        axis.set_xlabel('x'+unit_label)
        axis.set_ylabel('y'+unit_label)

        legend_dict = {label: handle for handle, label in zip(*axis.get_legend_handles_labels())}
        transformed_coordinates_marker = mlines.Line2D([], [], linewidth=0, markerfacecolor='none',
                                         markeredgecolor=transformed_coordinates_colour, marker='o')
        legend_dict[f'{transformed_coordinates_name} transformed ({transformed_coordinates.shape[0]})'] = transformed_coordinates_marker

        if show_pairs:
            pair_marker1 = mlines.Line2D([], [], linewidth=0, markerfacecolor='none',
                                         markeredgecolor=pair_colour, marker='o')
            pair_marker2 = mlines.Line2D([], [], linewidth=0, markerfacecolor=pair_colour,
                                         markeredgecolor='none', marker='.')
            legend_dict[f'matched pairs ({self.number_of_matched_points})'] = (pair_marker1, pair_marker2)

        axis.legend(legend_dict.values(),legend_dict.keys())

        if save:
            if save_path is None:
                save_path = self.save_path
            save_path = Path(save_path)
            figure.savefig(save_path.joinpath(self.name+'.png'), bbox_inches='tight', dpi=250)

    def get_transformation_direction(self, direction):
        """ Get inverse parameter based on direction

        Parameters
        ----------
        direction : str
            Another way of specifying the direction of the transformation, choose '<source_name>2<destination_name>' or
            '<destination_name>2<source_name>'

        Returns
        -------
        inverse : bool
            Specifier for the use of the forward or inverse transformation

        """

        if direction in [self.source_name + '2' + self.destination_name, 'source2destination']:
            inverse = False
        elif direction in [self.destination_name + '2' + self.source_name, 'destination2source']:
            inverse = True
        else:
            raise ValueError('Wrong direction')

        return inverse

    def calculate_inverse_transformation(self):
        self.transformation_inverse = type(self.transformation)(matrix=self.transformation._inv_matrix)

    def transform_coordinates(self, coordinates, inverse=False, direction=None):
        """Transform coordinates using the transformation

        Parameters
        ----------
        coordinates : Nx2 numpy.ndarray
            Coordinates to be transformed
        inverse : bool
            If True then the inverse transformation will be used (i.e. from destination to source)
        direction : str
            Another way of specifying the direction of the transformation, choose '<source_name>2<destination_name>' or
            '<destination_name>2<source_name>'

        Returns
        -------
        Nx2 numpy.ndarray
            Transformed coordinates

        """
        if direction is not None:
            inverse = self.get_transformation_direction(direction)

        if not inverse:
            current_transformation = self.transformation
        else:
            current_transformation = self.transformation_inverse

        return current_transformation(coordinates)

    def transform_image(self, image, inverse=False, direction=None):
        """Transform image using the transformation

            Parameters
            ----------
            image : NxM numpy.ndarray
                Image to be transformed
            inverse : bool
                If True then the inverse transformation will be used (i.e. from destination to source)
            direction : str
                Another way of specifying the direction of the transformation, choose '<source_name>2<destination_name>' or
                '<destination_name>2<source_name>'

            Returns
            -------
            NxM numpy.ndarray
                Transformed image

            """

        if direction is not None:
            inverse = self.get_transformation_direction(direction)

        # Note that this is the opposite direction of transformation
        if inverse:
            current_transformation = self.transformation
        else:
            current_transformation = self.transformation_inverse

        return skimage.transform.warp(image, current_transformation, preserve_range=True)

    def distance_matrix(self, crop=True, space='destination', max_distance=None, **kwargs):
        source = self.get_source(crop=crop, space=space)
        destination = self.get_destination(crop=crop, space=space)

        if max_distance is None:
            return distance_matrix(source, destination, **kwargs)
        else:
            source_tree = cKDTree(source)
            destination_tree = cKDTree(destination)
            return source_tree.sparse_distance_matrix(destination_tree, max_distance=max_distance, **kwargs).todense()


    def density_source(self, crop=False, space='source'):
        return self.get_source(crop).shape[0] / self.get_source_area(crop=crop, space=space)

    def density_destination(self, crop=False, space='destination'):
        return self.get_destination(crop).shape[0] / self.get_destination_area(crop=crop, space=space)

    def Ripleys_K(self, crop=True, space='destination'):
        from scipy.spatial import cKDTree
        # source_in_destination_kdtree = cKDTree(self.source_to_destination)
        # destination_kdtree = cKDTree(self.destination)


        # point_set_joint = np.vstack([point_set_1, point_set_2])
        # A = (point_set_joint.max(axis=0) - point_set_joint.min(axis=0)).prod()

        # source_vertices = self.get_source_vertices(crop=crop, space=space)
        # destination_vertices = self.get_destination_vertices(crop=crop, space=space)

        density_source = self.density_source(crop=crop, space=space)
        density_destination = self.density_destination(crop=crop, space=space)

        area_overlap = self.get_destination_area(crop=True, space=space)

        d = self.distance_matrix(crop=crop).flatten()
        d.sort()

        K = np.arange(len(d)) / (density_source * density_destination * area_overlap)

        return d, K

    def Ripleys_L_minus_d(self, crop=True, space='destination', plot=False):
        d, K = self.Ripleys_K(crop=crop, space=space)
        # K_random = np.pi * t ** 2
        L = np.sqrt(K / np.pi)

        if plot:
            fig, ax = plt.subplots()
            ax.plot(d, L - d)
            ax.set_xlabel('Distance')
            ax.set_ylabel('L-d')

        return d, L-d

    def Ripleys_L_minus_d_max(self, crop=True, space='destination'):
        d, L_minus_d = self.Ripleys_L_minus_d(crop=crop, space=space)
        max_L_minus_d = L_minus_d.max()
        d_at_max = d[np.where(max_L_minus_d == L_minus_d)][0]
        return d_at_max, max_L_minus_d

    def save(self, save_path=None, filetype='json'):
        """Save the current mapping in a file, so that it can be opened later.

        Parameters
        ----------
        filepath : str or pathlib.Path
            Path to file (including filename)
        filetype : str
            Choose classic to export a .coeff file for a linear transformation
            Choose yml to export all object attributes in a yml text file

        """

        if save_path is None and self.save_path is not None:
            save_path = self.save_path
        if not save_path.is_dir(): # save_path.suffix != '':
            self.name = save_path.with_suffix('').name
            save_path = save_path.parent

        filepath = Path(save_path).joinpath(self.name)

        if filetype == 'classic':
            if self.transformation_type == 'linear':
                coeff_filepath = filepath.with_suffix('.coeff')
                coefficients = self.transformation.params[[0, 0, 0, 1, 1, 1], [2, 0, 1, 2, 0, 1]]
                # np.savetxt(coeff_filepath, coefficients, fmt='%13.6g') # Same format used as in IDL code
                coefficients_inverse = self.transformation_inverse.params[[0, 0, 0, 1, 1, 1], [2, 0, 1, 2, 0, 1]]
                np.savetxt(coeff_filepath, np.concatenate((coefficients, coefficients_inverse)),
                           fmt='%13.6g')  # Same format used as in IDL code
            elif self.transformation_type == 'nonlinear':
                # saving kx,ky, still need to see how to read it in again
                map_filepath = filepath.with_suffix('.map')
                PandQ = self.transformation.params
                coefficients = np.concatenate((PandQ[0].flatten(), PandQ[1].flatten()), axis=None)
                # np.savetxt(map_filepath, coefficients, fmt='%13.6g') # Same format used as in IDL code
                PiandQi = self.transformation_inverse.params
                coefficients_inverse = np.concatenate((PiandQi[0].flatten(), PiandQi[1].flatten()), axis=None)
                np.savetxt(map_filepath, np.concatenate((coefficients, coefficients_inverse)),
                           fmt='%13.6g')  # Same format used as in IDL code
            else:
                raise TypeError('Mapping is not of type linear or nonlinear')

        elif filetype in ['yml', 'yaml', 'json']:
            attributes = self.__dict__.copy()
            for key in list(attributes.keys()):
                value = attributes[key]
                if type(value) in [str, int, float]:
                    continue
                elif isinstance(value, skimage.transform._geometric.GeometricTransform):
                    attributes[key] = value.params.tolist()
#TODO: solve issue with nonlinear mapping.transformation_type, which spits out a tuple of two arrays (4x4) instead np.shape(value.params.tolist())== (2, 15) 
                elif type(value).__module__ == np.__name__:
                    attributes[key] = value.tolist();
                else:
                    attributes.pop(key)

            if filetype in ['yml','yaml']:
                with filepath.with_suffix('.mapping').open('w') as yml_file:
                    yaml.dump(attributes, yml_file, sort_keys=False)
            elif filetype == 'json':
                with filepath.with_suffix('.mapping').open('w') as json_file:
                    json.dump(attributes, json_file, sort_keys=False)

        self.save_path = save_path

def overlap_vertices(vertices_A, vertices_B):
    polygon_A = Polygon(vertices_A)
    polygon_B = Polygon(vertices_B)
    if polygon_A.intersects(polygon_B):
        polygon_intersection = polygon_A.intersection(polygon_B)
    # return np.array(polygon_overlap.exterior.coords.xy).T[:-1]

        return np.array(polygon_intersection.boundary.coords)[:-1]
    else:
        return np.empty((0,2))

def area(vertices):
    if len(vertices) < 3:
        return 0
    else:
        return Polygon(vertices).area

def crop_coordinates_indices(coordinates, vertices):
    # return pth.Path(vertices).contains_points(coordinates)
    cropped_coordinates = crop_coordinates(coordinates, vertices)
    return np.array([c in cropped_coordinates for c in coordinates])

def crop_coordinates(coordinates, vertices):
    if len(vertices) > 0 and len(np.atleast_2d(coordinates)) > 0:
        # return np.atleast_2d(Polygon(vertices).intersection(MultiPoint(coordinates)))
        pointset_intersected = Polygon(vertices).intersection(MultiPoint(coordinates))
        if isinstance(pointset_intersected, MultiPoint):
            return np.atleast_2d(LineString(pointset_intersected.geoms).coords)
        elif isinstance(pointset_intersected, Point) and len(pointset_intersected.coords) > 0:
            return np.array(pointset_intersected.coords)
        else:
            return np.empty((0, 2))
        # Hopefully shapely 2.0 will be more consistent
    else:
        return np.empty((0, 2)) # np.atleast_2d([])

    # bounds.sort(axis=0)
    # selection = (coordinates[:, 0] > bounds[0, 0]) & (coordinates[:, 0] < bounds[1, 0]) & \
    #             (coordinates[:, 1] > bounds[0, 1]) & (coordinates[:, 1] < bounds[1, 1])
    # return coordinates[selection]


def determine_vertices(point_set, margin=0):
    return np.array(MultiPoint(point_set).convex_hull.buffer(margin, join_style=1).boundary.coords)[:-1]

import scipy.sparse

def singly_matched_pairs_within_radius(distance_matrix_, distance_threshold):
    matches = distance_matrix_ < distance_threshold
    sum_1 = matches.sum(axis=1) != 1
    sum_0 = matches.sum(axis=0) != 1
    matches[sum_1, :] = False
    matches[:, sum_0] = False
    return np.asarray(np.where(matches)).T

def distance_threshold_from_number_of_matches(radii, number_of_pairs, plot=True):
    # distance_threshold = np.sum(radii * number_of_pairs_summed) / np.sum(number_of_pairs_summed)
    distance_threshold = radii[np.where(number_of_pairs == number_of_pairs.max())][0]

    if plot:
        figure, axis = plt.subplots()
        axis.plot(radii, number_of_pairs)
        axis.axvline(distance_threshold)
        axis.set_xlabel('Radius')
        axis.set_ylabel('Count')

    return distance_threshold

def plot_circles(axis, coordinates, radius=6, **kwargs):
    circles = [plt.Circle((x, y), radius=radius) for x, y in coordinates]
    c = PatchCollection(circles, **kwargs)
    axis.add_collection(c)


if __name__ == "__main__":
    # Simulate source and destination point sets
    number_of_points = 10000
    transformation = AffineTransform(translation=[10, -10], rotation=1 / 360 * 2 * np.pi, scale=[0.98, 0.98])
    bounds = ([0, 0], [256, 512])
    crop_bounds = (None, None)
    fraction_missing = (0.95, 0.6)
    error_sigma = (0.5, 0.5)
    shuffle = True

    mapping = Mapping2.simulate(number_of_points, transformation, bounds, crop_bounds, fraction_missing,
                                error_sigma, shuffle)

    # mapping.transformation = AffineTransform(rotation=1/360*2*np.pi, scale=1.01, translation=[5,5])
    mapping.show_mapping_transformation()

    bounds = ((0.97, 1.1), (-0.05, 0.05), (-20, 20), (-20, 20))
    mapping.kernel_correlation(bounds, strategy='best1bin', maxiter=1000, popsize=50, tol=0.01,
                               mutation=0.25, recombination=0.7, seed=None, callback=None, disp=False,
                               polish=True, init='sobol', atol=0, updating='immediate', workers=1,
                               constraints=())
    mapping.show_mapping_transformation()

    mapping.find_distance_threshold()
    mapping.determine_matched_pairs()
    mapping.show_mapping_transformation()
