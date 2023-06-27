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

from tqdm import tqdm
from scipy.spatial import distance_matrix, cKDTree

from point_set import overlap_vertices, area, crop_coordinates, determine_vertices, vertices_with_margin
from icp import icp, nearest_neighbor_pair, nearest_neighbour_match, direct_match
from polywarp import PolywarpTransform
from polynomial import PolynomialTransform
from point_set_simulation import simulate_mapping_test_point_set
from kernel_correlation import kernel_correlation, compute_kernel_correlation
from cross_correlation import cross_correlate
from geometric_hashing import GeometricHashTable

from trace_analysis.decorators import return_none_when_executed_by_pycharm

class Mapping2:
    """Mapping class to find, improve, store and use a mapping between a source point set and a destination point set.

    Attributes
    ----------
    source : Nx2 numpy.ndarray
        Source coordinates.
    destination : Nx2 numpy.ndarray
        Destination coordinates.
    method : str
        Method for finding the transformation between source and destination.
    transformation_type : str
        Type of transformation used. Choose from: affine (default), similarity or polynomial.
    initial_transformation : skimage.transform._geometric.GeometricTransform or dict
        Initial transformation used as starting point by the perform_mapping method.
        Can be given as dictionary specifying translation, scale, rotation and shear.
    transformation  : skimage.transform._geometric.GeometricTransform or dict
        Transformation relating the source and destionaton point sets.
        Can be given as dictionary specifying translation, scale, rotation and shear.
    transformation_inverse : skimage.transform._geometric.GeometricTransform or dict
        Inverse transformation.
        Can be given as dictionary specifying translation, scale, rotation and shear.
    source_name : str
        Name of the source. Default is 'source'.
    destination_name : str
        Name of the source. Default is 'destination_name'.
    source_unit : str
        Unit of the source coordinates.
    destination_unit : str
        Unit of the destination coordinates.
    destination_distance_threshold : float
        Distance threshold in destination units for determining matched pairs between source and destination.
    name : str
        Name of the mapping, used as the default filename when saving the mapping object.
    label : str
        Additional attribute to store information.
    matched_pairs : Nx2 numpy.ndarray
        Array with the indices of source points and destination points that are matched.
    save_path : pathlib2.Path
        Path to the folder where the Mapping2 objects is or should be saved.


    """

    @classmethod
    def simulate(cls, number_of_points=200, transformation=None,
                 bounds=[[0, 0], [256, 512]], crop_bounds=(None, None), fraction_missing=(0.1, 0.1),
                 error_sigma=(0.5, 0.5), shuffle=True, seed=10532, show_correct=True):
        """
        Simulate a point set with a specific transformation and return as a Mapping2 object. It generates an original
        point set of which the source and destination point sets are subsets.

        Parameters
        ----------
        number_of_points : int
            Number of points in the original point set to simulate, the source and destination point sets are drawn from the original point set.
        transformation : skimage.transform._geometric.GeometricTransform
            Transformation to be used between the source and destination point sets. If None, a preset transformation
            will be used: SimilarityTransform(translation=[256, 10], rotation=1/360*2*np.pi, scale=[0.98, 0.98]).
        bounds : 2x2 numpy.ndarray or list
            Bounds of the original point set. Structured like coordinates,
            i.e. columns are x and y dimensions, rows are minimum and maximum values.
            [[x_min, y_min],[x_max, y_max]]
        crop_bounds : tuple of 2x2 numpy.ndarray or tuple of 2x2 list, optional
            Bounds to be used for cropping the overall dataset into the source and destination point set (before the transformation is applied).
        fraction_missing : tuple(float, float)
            Fraction of points that is removed from the original point set to obtain the source and destination point sets.
        error_sigma : tuple(float, float)
            Standard deviation of the Gaussian error applied to the original point set to obtain the source and destination point sets.
        shuffle : bool
            If True (default) then the points in the destination point set will be shuffled.
            If False the order of the destination points in the source and destination will be identical.
        seed : int
            Seed used for random number generator.
        show_correct : bool
            If True show the generated dataset with the correct transformation.

        Returns
        -------
        mapping : Mapping2
            Mapping2 object with the simulated source and destination point sets.
            While the "transformation" attribute is set with a unit transformation, mapping has an additional
            attribute "transformation_correct" containing the correct transformation. The returned mapping also has a
            "show_correct_mapping_transformation" method, to visualize the correct transformation.
        """


        if transformation is None:
            transformation = SimilarityTransform(translation=[256, 10], rotation=1/360*2*np.pi, scale=[0.98, 0.98])

        source, destination = simulate_mapping_test_point_set(number_of_points, transformation,
                                                              bounds, crop_bounds, fraction_missing,
                                                              error_sigma,
                                                              shuffle, seed=seed)

        class SimulatedMapping2(cls):
            def show_correct_mapping_transformation(self, *args, **kwargs):
                transformation_temp = self.transformation
                self.transformation = self.transformation_correct
                self.show_mapping_transformation(*args, **kwargs)
                self.transformation = transformation_temp

        mapping = SimulatedMapping2(source, destination)
        mapping.transformation_correct = transformation

        if show_correct:
            mapping.show_correct_mapping_transformation()

        return mapping

    @classmethod
    def load(cls, filepath):
        """
        Load a saved Mapping2 object

        Parameters
        ----------
        filepath : str or pathlib2.Path
            Filepath to the saved object.

        Returns
        -------
        mapping : Mapping2
            Saved Mapping2 object.
        """

        mapping = cls()
        filepath = Path(filepath)
        if filepath.suffix in ['.yml', '.yaml', '.json', '.mapping']:
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
                    setattr(mapping, key, value)
                except AttributeError:
                    pass

        elif filepath.suffix == '.nc':
            import xarray as xr
            ds = xr.load_dataset(filepath)
            for key, value in ds.attrs.items():
                setattr(mapping, key, value)
            for key, value in ds.items():
                setattr(mapping, key, value.values)

        mapping.transform = cls.transformation_types[mapping.transformation_type]
        mapping.transformation = mapping.transform(mapping.transformation)
        mapping.transformation_inverse = mapping.transform(mapping.transformation_inverse)
        if hasattr(mapping, 'transformation_correct'):
            mapping.transformation_correct = mapping.transform(mapping.transformation_correct)

        mapping.name = filepath.with_suffix('').name
        mapping.save_path = filepath.parent

        return mapping


    transformation_types = {'linear': AffineTransform,
                            'nonlinear': PolywarpTransform, # Change name to polywarp?
                            'polynomial': PolynomialTransform,
                            'affine': AffineTransform,
                            'similarity': SimilarityTransform}

    def __init__(self, source=None, destination=None, method=None,
                 transformation_type='affine', initial_transformation=None, transformation=None, transformation_inverse=None,
                 source_name='source', destination_name='destination',
                 source_unit=None, destination_unit=None, destination_distance_threshold=0,
                 name=None):
        """Initialize a Mapping2 object

        Parameters
        ----------
        source : Nx2 numpy.ndarray
            Source coordinates.
        destination : Nx2 numpy.ndarray
            Destination coordinates.
        method : str, optional
            Method for finding the transformation between source and destination.
        transformation_type : str
            Type of transformation used. Choose from: affine (default), similarity or polynomial.
        initial_transformation : skimage.transform._geometric.GeometricTransform or dict, optional
            Initial transformation used as starting point by the perform_mapping method.
            Can be given as dictionary specifying translation, scale, rotation and shear.
        transformation  : skimage.transform._geometric.GeometricTransform or dict, optional
            Transformation relating the source and destionaton point sets.
            Can be given as dictionary specifying translation, scale, rotation and shear.
        transformation_inverse : skimage.transform._geometric.GeometricTransform or dict, optional
            Inverse transformation.
            Can be given as dictionary specifying translation, scale, rotation and shear.
        source_name : str, optional
            Name of the source. Default is 'source'.
        destination_name : str, optional
            Name of the source. Default is 'destination_name'.
        source_unit : str, optional
            Unit of the source coordinates.
        destination_unit : str, optional
            Unit of the destination coordinates.
        destination_distance_threshold : float
            Distance threshold in destination units for determining matched pairs between source and destination.
        name : str
            Name of the mapping, used as the default filename when saving the mapping object.
        """

        self.name = name
        self.label = ''
        self.source_name = source_name
        self.source = np.array(source) #source=donor=left side image
        self.source_unit = source_unit
        self._source_distance_threshold = None
        self._destination_distance_threshold = destination_distance_threshold
        self._source_vertices = None
        self.destination_name = destination_name
        self.destination_unit = destination_unit
        self.destination = np.array(destination) #destination=acceptor=right side image
        self._destination_vertices = None
        self.method = method # Remove as input parameter and just set to None?
        self.matched_pairs = np.empty((0,2), dtype=int)
        self.save_path = None

        self.transformation_type = transformation_type

        if type(initial_transformation) is dict:
            initial_transformation = AffineTransform(**initial_transformation)
        self.initial_transformation = initial_transformation # Is this still necessary?

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


    # Function to make attributes from transformation available from the Mapping2 class
    def __getattr__(self, item):
        if ('transformation' in self.__dict__) and hasattr(self.transformation, item):
            return getattr(self.transformation, item)
        else:
            super().__getattribute__(item)

    def __eq__(self, other):
        return compare_objects(self, other) and compare_objects(other, self)

    # @property
    # def translation(self):
    #     return self.transformation[0:2,2]
    #
    # @property
    # @return_none_when_executed_by_pycharm
    # def magnification(self):
    #     return np.linalg.norm(self.transformation[:,0:2],axis=0)
    #
    # @property
    # @return_none_when_executed_by_pycharm
    # def rotation(self):
    #     # rotation_matrix = self.transformation[0:2, 0:2]/self.magnification[0]
    #     # return np.arctan2(rotation_matrix[0,1]-rotation_matrix[1,0],rotation_matrix[0,0]+rotation_matrix[1,1])/(2*np.pi)*360
    #     return np.arctan2(self.transformation[0, 1], self.transformation[0, 0]) / (2 * np.pi) * 360
    #
    # @property
    # def reflection(self):
    #     return np.array([np.sign(self.transformation[0, 0]), np.sign(self.transformation[1, 1])])

    @property
    @return_none_when_executed_by_pycharm
    def file_path(self):
        """
        pathlib2.Path : Path where the mapping file is (to be) saved.
        """
        # TODO: Make the extension dependent on the actual file
        return self.save_path.joinpath(self.name).with_suffix('.mapping')

    @property
    @return_none_when_executed_by_pycharm
    def source_to_destination(self):
        """Nx2 numpy.ndarray : Source coordinates transformed to destination space"""

        return self.transform_coordinates(self.source)

    @property
    @return_none_when_executed_by_pycharm
    def destination_to_source(self):
        """Nx2 numpy.ndarray : Destination coordinates transformed to source space"""

        return self.transform_coordinates(self.destination, inverse=True)

    @property
    @return_none_when_executed_by_pycharm
    def source_vertices(self):
        """Nx2 numpy.ndarray : Vertices of source

        If no source vertices are set, the vertices are determined from the convex hull.
        """
        if self._source_vertices is None:
            return determine_vertices(self.source)#, self.source_distance_threshold)
        else:
            return self._source_vertices

    @source_vertices.setter
    def source_vertices(self, vertices):
        self._source_vertices = vertices

    @property
    @return_none_when_executed_by_pycharm
    def destination_vertices(self):
        """Nx2 numpy.ndarray : Vertices of destination

        If no destination vertices are set, the vertices are determined from the convex hull.
        """
        if self._destination_vertices is None:
            return determine_vertices(self.destination)#, self.destination_distance_threshold)
        else:
            return self._destination_vertices

    @destination_vertices.setter
    def destination_vertices(self, vertices):
        self._destination_vertices = vertices

    @property
    @return_none_when_executed_by_pycharm
    def source_cropped_vertices(self): # or crop_vertices_in_source
        """Nx2 numpy.ndaarray : Vertices of the overlapping area (intersection) of the source and destination areas
        in source space"""
        return overlap_vertices(self.source_vertices, self.transform_coordinates(self.destination_vertices, inverse=True))

    @property
    @return_none_when_executed_by_pycharm
    def source_cropped(self):
        """Nx2 numpy.ndaarray : Source points that overlap with the destination in source space"""
        return crop_coordinates(self.source, self.source_cropped_vertices)

    @property
    @return_none_when_executed_by_pycharm
    def destination_cropped_vertices(self): # or crop_vertices_in_destination
        """Nx2 numpy.ndaarray : Vertices of the overlapping area (intersection) of the source and destination areas
        in destination space"""
        return overlap_vertices(self.transform_coordinates(self.source_vertices), self.destination_vertices)

    @property
    @return_none_when_executed_by_pycharm
    def destination_cropped(self):
        """Nx2 numpy.ndaarray : Destination points that overlap with the source in destination space"""
        return crop_coordinates(self.destination, self.destination_cropped_vertices)

    def get_source(self, crop=False, space='source', margin=None):
        """Getter for the source point set.

        Parameters
        ----------
        crop : bool or str, optional
            If True or 'source', the point set is cropped to the area of the destination.
        space : str, optional
            In which coordinate space to return the point set. Either 'source' or 'destination'.
        margin : float or int, optional
            Margin used for cropping.


        Returns
        -------
        Nx2 numpy.ndaarray
            Source point set
        """
        if crop in ['destination', False]:
            source = self.source
        elif crop in ['source', True]:
            # source = self.source_cropped
            source = crop_coordinates(self.source, self.get_source_vertices(crop=crop, margin=margin))

        if space in ['destination', self.destination_name]:
            source = self.transform_coordinates(source)

        return source

    def get_destination(self, crop=False, space='destination', margin=None):
        """Getter for the destination point set.

        Parameters
        ----------
        crop : bool or str, optional
            If True or 'destination', the point set is cropped to the area of the source.
        space : str, optional
            In which coordinate space to return the point set. Either 'source' or 'destination'.
        margin : float or int, optional
            Margin used for cropping.


        Returns
        -------
        Nx2 numpy.ndaarray
            Destination point set
        """
        if crop in ['source', False]:
            destination = self.destination
        elif crop in ['destination', True]:
            # destination = self.destination_cropped
            destination = crop_coordinates(self.destination, self.get_destination_vertices(crop=crop, margin=margin))

        if space in ['source', self.source_name]:
            destination = self.transform_coordinates(destination, inverse=True)

        return destination

    def get_source_vertices(self, crop=False, space='source', margin=None):
        """Getter for vertices of the source point set.

        Parameters
        ----------
        crop : bool or str, optional
            If True or 'source', the vertices of the overlapping area with the destination point set are given.
        space : str, optional
            In which coordinate space to return the vertices. Either 'source' or 'destination'.
        margin : float or int, optional
            Margin used for cropping.

        Returns
        -------
        Nx2 numpy.ndaarray
            Coordinates of the source vertices
        """
        if crop in ['destination', False]:
            source_vertices = self.source_vertices
        elif crop in ['source', True]:
            source_vertices = self.get_overlap_vertices(space='source')

        if space in ['destination', self.destination_name]:
            source_vertices = self.transform_coordinates(source_vertices)

        if margin is not None:
            source_vertices = vertices_with_margin(source_vertices, margin)

        return source_vertices

    def get_destination_vertices(self, crop=False, space='destination', margin=None):
        """Getter for vertices of the destination point set.

        Parameters
        ----------
        crop : bool or str, optional
            If True or 'destination', the vertices of the overlapping area with the source point set are given.
        space : str, optional
            In which coordinate space to return the vertices. Either 'source' or 'destination'.
        margin : float or int, optional
            Margin used for cropping.

        Returns
        -------
        Nx2 numpy.ndaarray
            Coordinates of the destination vertices
        """

        if crop in ['source', False]:
            destination_vertices = self.destination_vertices
        elif crop in ['destination', True]:
            destination_vertices = self.get_overlap_vertices(space='destination')

        if space in ['source', self.source_name]:
            destination_vertices = self.transform_coordinates(destination_vertices, inverse=True)

        if margin is not None:
            destination_vertices = vertices_with_margin(destination_vertices, margin)

        return destination_vertices

    def get_overlap_vertices(self, space='source'):
        """Getter for vertices of the overlap between the source and destination point sets.

        Parameters
        ----------
        space : str, optional
            In which coordinate space to return the vertices. Either 'source' or 'destination'.

        Returns
        -------
        Nx2 numpy.ndaarray
            Coordinates of the overlap vertices
        """
        return overlap_vertices(self.get_source_vertices(space=space), self.get_destination_vertices(space=space))

    @property
    @return_none_when_executed_by_pycharm
    def source_distance_threshold(self):
        """float : Distance threshold in source space.

        If destination_distance_threshold is set then the source_distance_threshold is derived using the scale of the
        transformation.
        """
        if self._source_distance_threshold is not None:
            return self._source_distance_threshold
        elif self._destination_distance_threshold is not None:
            return self._destination_distance_threshold / np.max(self.transformation.scale)
        else:
            raise ValueError('No distance threshold set')

    @source_distance_threshold.setter
    def source_distance_threshold(self, value):
        self._source_distance_threshold = value
        self._destination_distance_threshold = None

    @property
    @return_none_when_executed_by_pycharm
    def destination_distance_threshold(self):
        """float : Distance threshold in source space.

        If source_distance_threshold is set then the destination_distance_threshold is derived using the scale of the
        transformation.
        """
        if self._destination_distance_threshold is not None:
            return self._destination_distance_threshold
        elif self._source_distance_threshold is not None:
            return self._source_distance_threshold * np.max(self.transformation.scale)
        else:
            raise ValueError('No distance threshold set')

    @destination_distance_threshold.setter
    def destination_distance_threshold(self, value):
        self._source_distance_threshold = None
        self._destination_distance_threshold = value

    def find_distance_threshold(self, method='single_match_optimization', **kwargs):
        """Find distance optimal distance threshold and automatically sets the destination distance threshold.

        Parameters
        ----------
        method : str
            Method to use for finding the distance threshold. Choose from:
                -   'single_match_optimization': Optimizes the number of singly-matched pairs,
                    i.e. the source point having only a single destination point within the distance threshold and the
                    destination point having only a single source point within the distance threshold.
                    See the Mapping2.single_match_optimization method.
        kwargs
            Keyword arguments to pass to method.
        """
        if method == 'single_match_optimization':
            self.single_match_optimization(**kwargs)
        else:
            raise ValueError('Unknown method')

    def number_of_single_matches_for_radii(self, distance_thresholds):
        """The number of singly-matched pairs for the provided distance thresholds.

        Here singly-matched pairs are defined as the source point having only a single destination point within the
        distance threshold and the destination point having only a single source point within the distance threshold.

        Parameters
        ----------
        distance_thresholds : 1D numpy.ndaarray
            Distance thresholds to test.

        Returns
        -------
        1D numpy.ndaarray
            Number of singly-matched points for each distance threshold.
        """
        distance_matrix_ = self.distance_matrix(crop=True)
        number_of_pairs = np.array([len(singly_matched_pairs_within_radius(distance_matrix_, r)) for r in distance_thresholds])
        return number_of_pairs

    def single_match_optimization(self, maximum_radius=20, number_of_steps=100, plot=True):
        """ Find the distance threshold with the highest number of singly matched pairs.

        Test thresholds between 0 and maximum_radius.
        Destination_distance_threshold is automatically set to the found value.

        Parameters
        ----------
        maximum_radius : float
            Maximum distance threshold to test.
        number_of_steps : int
            Number of steps on the interval
        plot : bool
            If True, a histogram of the number of pairs per distance threshold will be plotted
        """
        distance_thresholds = np.linspace(0, maximum_radius, number_of_steps)
        number_of_pairs = self.number_of_single_matches_for_radii(distance_thresholds)
        self.destination_distance_threshold = distance_threshold_from_number_of_matches(distance_thresholds, number_of_pairs, plot=plot)

    def determine_matched_pairs(self, distance_threshold=None, point_set_name='all'):
        """Find pairs of source and destination points that closer than the distance threshold and are singly matched,
        i.e. each source point only has one destination point within the threshold and vice versa.

        Sets the matched_pairs attribute.

        Parameters
        ----------
        distance_threshold : int or float
            Distance threshold for match determination.
        point_set_name : str
            If 'source' then only the source point should have a single match in the destination point set.
            the destination point can have multiple matches in the source point set.
            If 'destination' then only the destination point should have a single match in the source point set,
            the source point can have multiple matches in the destination point set.
            If 'all' then both source and destination should have a single match.
        """

        #TODO: add crop
        if distance_threshold is None:
            distance_threshold = self.destination_distance_threshold

        distance_matrix_ = self.distance_matrix(crop=False)

        self.matched_pairs = singly_matched_pairs_within_radius(distance_matrix_, distance_threshold, point_set_name=point_set_name)

    def number_of_matches_for_source_and_destination(self, distance_threshold=None, matches_per_point=[0,1], crop=True):
        """Gives the number of source and destination points that respectively have a specified number of matches
        with the destination and source point sets.

        Parameters
        ----------
        distance_threshold : int, float, list or numpy.ndarray
            Distance threshold(s) to determine the number of matches for. If None (default) then the
            destination_distance_threshold attribute is used.
        matches_per_point : int or list of int
            Number of matches to report
        crop : bool or str, optional
            If True: the overlap between the source and destination point set is used.
            If 'source': only the source is cropped to the destination point set.
            If 'destination': only the destination is cropped to the source point set.
            If False: no cropping is applied.

        Returns
        -------
        xarray.DataArray
            Number of source or destination points with specific distance threshold

        """
        if distance_threshold is None:
            distance_threshold = self.destination_distance_threshold
        if not (isinstance(distance_threshold, list) or isinstance(distance_threshold, np.ndarray)):
            distance_threshold = [distance_threshold]

        import xarray as xr
        data = xr.DataArray(0, dims=('R', 'matches_per_point', 'reference'),
                            coords={'R': distance_threshold, 'matches_per_point': matches_per_point,
                                    'reference': ['source', 'destination']})

        margin = -np.max(distance_threshold) * 1.5

        distance_matrix_sd = self.distance_matrix(crop=crop, margin=(margin, 0))
        distance_matrix_ds = self.distance_matrix(crop=crop, margin=(0, margin))
        for m in matches_per_point:
            for i, R in enumerate(distance_threshold):
                data[i, m, 0] = np.array(number_of_matches_within_radius(distance_matrix_sd, R, matches_per_point=m))[0]
                data[i, m, 1] = np.array(number_of_matches_within_radius(distance_matrix_ds, R, matches_per_point=m))[1]
        return data

    def pair_coordinates(self, point_set_name='destination', space='destination'):
        """ Get the coordinates of the paired source or destination points.

        Parameters
        ----------
        point_set_name : str
            Name from the point set to use the coordinates from. Either 'source' or 'destination' (default).
        space : str, optional
            In which coordinate space to return the coordinates. Either 'source' or 'destination' (default).


        Returns
        -------
        Nx2 numpy.ndaarray
            Coordinates of paired points
        """
        if point_set_name == 'source':
            return self.get_source(space=space)[self.matched_pairs[:,0]]
        elif point_set_name == 'destination':
            return self.get_destination(space=space)[self.matched_pairs[:,1]]

    def pair_distances(self, space='destination', show=False, **kwargs):
        """ Distances between the paired source and destination points.

        Parameters
        ----------
        space : str, optional
            In which coordinate space to return the coordinates. Either 'source' or 'destination' (default).
        show : bool
            If True, show a histogram of the distances.
        kwargs
            Keyword arguments passed to matplotlib.pyplot.histogram.

        Returns
        -------
        1D numpy.ndaarray
            Distances between paired points.
        """
        xy_distances = self.pair_coordinates(point_set_name='source', space=space) - \
                       self.pair_coordinates(point_set_name='destination', space=space)

        pair_distances = np.linalg.norm(xy_distances, axis=1)

        if show:
            figure, axis = plt.subplots()
            axis.hist(pair_distances, bins=100, **kwargs)
            axis.set_xlabel('Distance'+self.get_unit_label(space))
            axis.set_ylabel('Count')
            axis.set_title('Pair distances')

        return pair_distances

    @property
    @return_none_when_executed_by_pycharm
    def number_of_source_points(self):
        """int : Number of points in the source point set.
        """
        return self.source.shape[0]

    @property
    @return_none_when_executed_by_pycharm
    def number_of_destination_points(self):
        """int : Number of points in the destination point set.
        """
        return self.destination.shape[0]

    @property
    @return_none_when_executed_by_pycharm
    def number_of_matched_points(self):
        """int : Number of matched pairs between the source and destination point set
        """
        # distances, source_indices, destination_indices = \
        #     nearest_neighbor_pair(self.source_to_destination, self.destination)
        #
        # return np.sum(distances < self.destination_distance_threshold)

        return self.matched_pairs.shape[0]

    def fraction_of_source_matched(self, crop=False):
        """Fraction of the source points that is matched to a destination point

        Parameters
        ----------
        crop : bool or str, optional
            If True: the overlap between the source and destination point set is used.
            If 'source': only the source is cropped to the destination point set.
            If 'destination': only the destination is cropped to the source point set.
            If False: no cropping is applied.

        Returns
        -------
        float
            Fraction of source that is matched
        """
        return self.number_of_matched_points / self.get_source(crop).shape[0]

    def fraction_of_destination_matched(self, crop=False):
        """Fraction of the destination points that is matched to a source point

        Parameters
        ----------
        crop : bool or str, optional
            If True: the overlap between the source and destination point set is used.
            If 'source': only the source is cropped to the destination point set.
            If 'destination': only the destination is cropped to the source point set.
            If False: no cropping is applied.

        Returns
        -------
        float
            Fraction of destination that is matched
        """
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
    @return_none_when_executed_by_pycharm
    def source_area(self):
        """float : Area of the source point set."""
        return area(self.source_vertices)

    @property
    @return_none_when_executed_by_pycharm
    def source_in_destination_area(self):
        """float : Area of the source point set when transformed to destination space."""
        return area(self.transformation(self.source_vertices))

    @property
    @return_none_when_executed_by_pycharm
    def source_cropped_area(self):
        """float : Area of the source point set cropped to the destination point set."""
        return area(self.source_cropped_vertices)

    @property
    @return_none_when_executed_by_pycharm
    def destination_area(self):
        """float : Area of the destination point set."""
        return area(self.destination_vertices)

    @property
    @return_none_when_executed_by_pycharm
    def destination_cropped_area(self):
        """float : Area of the destination point set cropped to the source point set."""
        return area(self.destination_cropped_vertices)

    def get_source_area(self, crop=False, space='destination'):
        """Get the area of the source point set, i.e. within the source point vertices.

        Parameters
        ----------
        crop : bool or str, optional
            If True: the overlap between the source and destination point set is used.
            If 'source': only the source is cropped to the destination point set.
            If 'destination': only the destination is cropped to the source point set.
            If False: no cropping is applied.
        space : str, optional
            In which coordinate space to return the point set. Either 'source' or 'destination'.

        Returns
        -------
        float
            Source area
        """
        return area(self.get_source_vertices(crop=crop, space=space))

    def get_destination_area(self, crop=False, space='destination'):
        """Get the area of the destination point set, i.e. within the destination point vertices.

           Parameters
           ----------
           crop : bool or str, optional
               If True: the overlap between the source and destination point set is used.
               If 'source': only the source is cropped to the destination point set.
               If 'destination': only the destination is cropped to the source point set.
               If False: no cropping is applied.
           space : str, optional
               In which coordinate space to return the point set. Either 'source' or 'destination'.

           Returns
           -------
           float
               Destination area
           """
        return area(self.get_destination_vertices(crop=crop, space=space))

    @property
    @return_none_when_executed_by_pycharm
    def transformation_type(self):
        """str : Transformation type used for mapping.

        Currently supported transformation types are:
        - 'linear' or 'affine': affine transform using skimage.transform.AffineTransform
        - 'similarity': similarity transform using skimage.transform.SimilarityTransform
        - 'nonlinear': polynomial transform corresponding to the IDL polywarp transform
        - 'polynomial': polynomial transform using skimage.transform.PolynomialTransform

        Available transformation types are stored in the Mapping2.transformation_types class attribute.

        If set, then the transform attribute is set to the transformation class corresponding to the set value.
        """
        #TODO: derive the transformation type from self.transform instead of the other way around.
        return self._transformation_type

    @transformation_type.setter
    def transformation_type(self, value):
        self.transform = self.transformation_types[value]
        self._transformation_type = value

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
        kwargs
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

        The found transformation is stored in the object transformation parameter.

        Note
        ----
        The number and the order of source points should be equal to the number and the order of destination points.

        Parameters
        ----------
        transformation_type : str
            Type of transformation used, either linear or polynomial can be chosen.
            If not specified the object transformation_type is used.
        kwargs
            Keyword arguments passed to the direct match function.

        """
        if transformation_type is not None:
            self.transformation_type = transformation_type

        self.transformation, self.transformation_inverse, error = \
            direct_match(self.source, self.destination, transform=self.transform, return_inverse=True, **kwargs)

        print(f'Direct match\n'
              f'Mean-squared error: {error}')

    def nearest_neighbour_match(self, distance_threshold=None, transformation_type=None, **kwargs):
        """Find transformation from source to destination points by matching nearest neighbours.

        Two-way nearest neighbours are detected, i.e. the source point should be the nearest neighbour of the
        destination point and vice versa. Only nearest neighbours closer than the distance threshold are used to find
        the transformation.

        The found transformation is stored in the object transformation parameter.

        Note
        ----
        The current transformation is used as starting point for the algorithm.

        Note
        ----
        The printed error is based on the points selected for matching.

        Parameters
        ----------
        distance_threshold : float
            Distance threshold for nearest neighbour match in destination space. Only nearest neighbours with a
            distance smaller than the distance threshold are used.
        transformation_type : str
            Type of transformation used. For options see Mapping2.transformation_type.
            If not specified the object transformation_type is used.
        kwargs
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

        In the iterative closest point algorithm, the two-way nearest neighbours are found and these are used to
        find the most optimal transformation. Subsequently, the source is transformed according to this
        transformation. This process is repeated until the changes detected are below a tolerance level.

        The iterative closest point algorithm can be used in situations when deviations between the two point sets
        are relatively small.

        The found transformation is stored in the object's transformation attribute.

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
            i.e. nearest-neighbours with distances smaller than distance_threshold.
        kwargs
            Keyword arguments passed to the algorithm.
        """

        self.transformation, self.transformation_inverse, error, number_of_iterations = \
            icp(self.source, self.destination, distance_threshold_final=distance_threshold,
                initial_transformation=self.transformation, transform_final=self.transform, **kwargs)

        print(f'Iterative closest point match\n'
              f'Mean-squared error: {error}\n'
              f'Number of iterations: {number_of_iterations}')

    def kernel_correlation(self, bounds=((0.97, 1.02), (-0.05, 0.05), (-10, 10), (-10, 10)), sigma=1, crop=False,
                           plot=False, **kwargs):
        """Find transformation from source to destination points using an kernel correlation algorithm.

        The found transformation is stored in the object's transformation attribute.

        Note
        ----
        The found kernel correlation value is stored in the objects mapping statistics with key
        "kernel_correlation_value".

        Parameters
        ----------
        bounds : tuple of tuples
            The bounds for variation of the transformation parameters.
            The number of tuples passed determines the number of varied parameters.
            - Two tuples indicates variation of the x and y translation.
            - Three tuples indicate rotation, x-translation and y-translation.
            - Four tuples indicate scale, rotation, x-translation and y-translation.
            - Five tuples indicate scale, rotation, shear, x-translation and y-translation.
            - Six tuples indicate scale-x, scale-y, rotation, shear, x-translation and y-translation.
            The rotation parameter is given in radians.
        sigma : float
            The standard deviation of the Gaussian kernel used.
        crop : bool or str, optional
               If True: the overlap between the source and destination point set is used.
               If 'source': only the source is cropped to the destination point set.
               If 'destination': only the destination is cropped to the source point set.
               If False: no cropping is applied.
        plot : bool
            Whether to plot the mapping process.
        kwargs
            Keyword arguments passed to the algorithm.
        """
        transformation, result = kernel_correlation(self.get_source(crop, space='destination'), self.get_destination(crop),
                                                    bounds, sigma, plot=plot, **kwargs)
        self.transformation = AffineTransform(matrix=(self.transformation + transformation).params)
        self.transformation_inverse = type(self.transformation)(self.transformation._inv_matrix)
        self.mapping_statistics = {'kernel_correlation_value': result.fun}

    def kernel_correlation_score(self, sigma=1, crop=False, per_point_pair=False):
        """Determine the kernel correlation score for the current transformation.

        Parameters
        ----------
        sigma : float
            The standard deviation of the Gaussian kernel used.
        crop : bool or str, optional
            If True: the overlap between the source and destination point set is used.
            If 'source': only the source is cropped to the destination point set.
            If 'destination': only the destination is cropped to the source point set.
            If False: no cropping is applied.
        per_point_pair : bool
            If True the kernel correlation score is returned for each point pair in a 2D matrix.
            If False a summed kernel correlation score is returned.

        Returns
        -------
        float or NxM numpy.ndarray
            Kernel correlation score
        """
        return compute_kernel_correlation(self.transformation, self.get_source(crop), self.get_destination(crop),
                                          sigma=sigma, per_point_pair=per_point_pair)

    def cross_correlation(self, peak_detection='auto', kernel_size=7, gaussian_sigma=1, divider=5, crop=False, space='destination',
                          plot=False, axes=None):
        """Perform cross correlation on synthesized images of the two datasets.

        The current transformation is used as a starting point. Because cross-correlation only varies translation,
        other transformation parameters such rotation and scale should already be correct.
        The found transformation is stored in the object's transformation attribute.

        Parameters
        ----------
        peak_detection : str
            If 'auto': the cross-correlation peak is automatically determined by finding the most intense pixel in the
            cross-correlation image.
            If 'manual': cross-correlation image is created, but the peak has to be determined manually.
            The found pixel coordinates of the peak in the image can be set using the method 'set_correlation_peak_coordinates',
            the correct transformation will then be calculated.
        kernel_size : int
            The size in pixels of the gaussian kernel. The kernel will have shape (kernel_size x kernel_size).
        gaussian_sigma : float
            The standard deviation of the gaussian kernel. The same standard deviation is used for x and y.
        divider : int or float
            The number by which the point set coordinates are divided to obtain pixel indices.
            (Translating the lowest coordinates to the pixel (0,0) is done automatically.)
        crop : int or float, optional
            If True: the overlap between the source and destination point set is used.
            If 'source': only the source is cropped to the destination point set.
            If 'destination': only the destination is cropped to the source point set.
            If False: no cropping is applied.
        space : str, optional
            In which coordinate space to perform the correlation. Either 'source' or 'destination'.
        plot : bool
            If True: plots of the synthetic images and the cross-correlation image are shown.
        axes : list of matplotlib.axis.Axis
            Three axes to use for plotting.

        Note
        ----
        Automatic background subtraction of the cross-correlation image is performed by subtracting a minimum filtered
        version of the cross-correlation image. This reduces the difference in intensity between the edges and the center
        that naturally occur when doing cross-correlation, and it thus helps to find the correct peak.
        """
        #TODO: Make Gaussian without multiplying mask.
        #TODO: Add standard deviation as a parameter. (And change gaussian_width name to kernal_width?)
        #TODO: Use maximum_number_of_pixels instead of divider.
        #TODO: Add option to automatically save plot.

        if self.transformation is None:
            self.transformation = AffineTransform()
            self.transformation_inverse = AffineTransform()

        if plot and axes is None:
            figure, axes = plt.subplots(1, 4)

        correlation, self.correlation_conversion_function = cross_correlate(self.get_source(crop, space), self.get_destination(crop, space),
                                                                            kernel_size=kernel_size, gaussian_sigma=gaussian_sigma, divider=divider,
                                                                            subtract_background=True, plot=plot, axes=axes)

        self.correlation_space = space

        if peak_detection == 'auto':
            #TODO: Fit peak to gaussian to determine location with sub-pixel accuracy???
            correlation_peak_coordinates = np.array(np.where(correlation==correlation.max())).flatten()[::-1]+0.5
            if plot:
                bounds = np.array([axes[2].get_xlim(), axes[2].get_ylim()])
                # pixel_size = np.diff(bounds).flatten()/(np.array(axes[3].get_images()[0].get_size()[::-1])+1)
                origin = bounds[:,0]
                peak_coordinates_in_image = origin + correlation_peak_coordinates * divider
                axes[2].plot(*peak_coordinates_in_image.T, marker='o', markerfacecolor='none', markeredgecolor='r')
                axes[3].plot(*peak_coordinates_in_image.T, marker='o', markerfacecolor='none', markeredgecolor='r')
            self.set_correlation_peak_coordinates(correlation_peak_coordinates)
        elif peak_detection == 'manual':
            pass
        else:
            raise ValueError(f'{peak_detection} is not a valid value for peak_detection. Use either "auto" or "manual".')

    def set_correlation_peak_coordinates(self, correlation_peak_coordinates):
        """Manually set the peak coordinates after cross-cross correlation.

        Note
        ----
        This function can only be used after cross-correlation is run.

        Parameters
        ----------
        correlation_peak_coordinates : list, tuple or 1D numpy.ndarray
            Pixel coordinates of the manually found peak in the cross-correlation image.
        """
        if not hasattr(self, 'correlation_conversion_function') and self.correlation_conversion_function is not None:
            raise RuntimeError('Run cross_correlation first')
        transformation = self.correlation_conversion_function(correlation_peak_coordinates) # is this the correct direction
        if self.correlation_space == 'source':
            self.transformation = AffineTransform(matrix=(transformation + self.transformation).params)
        elif self.correlation_space == 'destination':
            self.transformation = AffineTransform(matrix=(self.transformation + transformation).params)
        else:
            raise ValueError('Unkown correlation_space, use either "source" or "destination"')
        self.transformation_inverse = AffineTransform(matrix=self.transformation._inv_matrix)
        self.correlation_conversion_function = None

    def geometric_hashing(self, method='one_by_one', tuple_size=4, maximum_distance_source=None,
                          maximum_distance_destination=None, **kwargs):
        """Perform mapping by geometric hashing.

        The found transformation is stored in the object's transformation attribute.

        Parameters
        ----------
        method : str
            The method to use:
            - 'one_by_one': tuples are generated, matched and tested sequentially.
            - 'abundant_transformations': All tuples are generated and matched. Transformation parameters are determined
               for matched tuples and the highest density cluster is determined.
        tuple_size : int
            Number of points per tuple.
        maximum_distance_source : float
            Maximum distance between the two outermost points in a tuple of the source point set.
        maximum_distance_destination : float
            Maximum distance between the two outermost points in a tuple of the destination point set.
        kwargs
            Keyword arguments to pass to the geometric hashing query function.
        """
        hash_table = GeometricHashTable([self.destination], source_vertices=self.source_vertices,
                                        initial_source_transformation=self.transformation,
                                        number_of_source_bases=20, number_of_destination_bases='all',
                                        tuple_size=tuple_size, maximum_distance_source=maximum_distance_source,
                                        maximum_distance_destination=maximum_distance_destination)

        if method == 'one_by_one':
            found_transformation = hash_table.query(self.source, **kwargs)
        elif method == 'abundant_transformations':
            found_transformation = hash_table.query_tuple_transformations([self.source], **kwargs)
        else:
            raise ValueError(f'Unknown method {method}')

        self.transformation = found_transformation
        self.transformation_inverse = type(found_transformation)(matrix=self.transformation._inv_matrix)


    def get_unit(self, space):
        """Get the unit of the source or destination point set.

        Parameters
        ----------
        space : str
            For which coordinate space to obtain the unit. Either 'source' or 'destination'.

        Returns
        -------
        str
            The unit of the point set.
        """
        return self.__getattribute__(f'{space}_unit')

    def get_unit_label(self, space):
        """Get the string to add the axis label during plotting.
        This is either the unit between parentheses, or if no unit is set it is an empty string.

        Parameters
        ----------
        space : str
            For which coordinate space to obtain the unit. Either 'source' or 'destination'.

        Returns
        -------
        str
            Unit label
        """
        if self.get_unit(space) is not None:
            return f' ({self.get_unit(space)})'
        else:
            return ''

    def show_source(self, **kwargs):
        """Plot the source point set.

        Parameters
        ----------
        kwargs
            Keyword arguments to pass to the show function. Keyword arguments starting with "show" should not be used.

        Returns
        -------
        figure : matplotlib.figure.Figure
            Figure of the plot
        axis : matplotlib.axes.Axes
            Axis of the plot
        """
        return self.show(show_source=True, show_destination=False, show_transformed_coordinates=False, **kwargs)

    def show_destination(self, **kwargs):
        """Plot the destination point set.

        Parameters
        ----------
        kwargs
            Keyword arguments to pass to the show function. Keyword arguments starting with "show" should not be used.

        Returns
        -------
        figure : matplotlib.figure.Figure
            Figure of the plot
        axis : matplotlib.axes.Axes
            Axis of the plot
        """
        return self.show(show_source=False, show_destination=True, show_transformed_coordinates=False, **kwargs)

    def show_mapping_transformation(self, *args, **kwargs):
        """Show the mapping. Identical to `Mapping2.show()`.
        """
        return self.show(*args, **kwargs)

    def show(self, axis=None, show_source=False, show_destination=False, show_transformed_coordinates=True,
             show_pairs=True, crop=False, inverse=False, source_colour='forestgreen', destination_colour='r',
             pair_colour='b', use_distance_threshold=False, save=False, save_path=None, show_legend=False,
             return_plot=False):
        #TODO: Always return plot?
        """Show the mapping.

        Show the source and destination points in a scatter plot.

        Parameters
        ----------
        axis : matplotlib.axis.Axis
            The axis to use for plotting. If no axis is given a new figure and axis are created.
        show_source : bool
            If True: the original source point set is shown.
        show_destination : bool
            If True: the original destination point set is shown.
        show_transformed_coordinates : bool
            If True: the transformed point set is shown. If 'inverse' is set to False, this is the transformed source
            point set. If 'inverse' is set to True, this is the transformed destination point set.
        show_pairs : bool
            If True: the paired points are indicated in a different colour.
        crop : bool or str
            If True: the overlap between the source and destination point set is used.
            If 'source': only the source is cropped to the destination point set.
            If 'destination': only the destination is cropped to the source point set.
            If False: no cropping is applied.
        inverse : bool
            Indicates the direction of the transformation to use.
            If False: the source is transformed to destination space.
            If True: the destination is transformed to the source space.
        source_colour : str
            Colour used for the source points.
        destination_colour : str
            Colour used for the destination points.
        pair_colour : str
            Colour used to indicate paired points.
        use_distance_threshold : bool
            If True: the open circles of the transformed point set will have the size of the source or destination
            distance threshold. Note: using this option slows down plotting so it is not recommended for large datasets.
        save : bool
            If True: the figure is saved in the save path.
        save_path : str or pathlib2.Path
            Folderpath to save the figure to. If no path is given, the object's save_path is used.
            For the name of the file, the object's name attribute is used.
        show_legend : bool
            If True: the legend is shown.
        return_plot : bool
            If True: the figure and axis are returned by the function.

        Returns
        -------
        figure : matplotlib.figure.Figure
            Figure of the plot
        axis : matplotlib.axes.Axes
            Axis of the plot
        """

        if axis is None:
            figure, axis = plt.subplots()
        else:
            figure = axis.figure

        source = self.get_source(crop)
        destination = self.get_destination(crop)

        if len(self.matched_pairs) == 0:
            show_pairs = False

        if show_transformed_coordinates:
            if not inverse:
                all_transformed_coordinates = self.transform_coordinates(self.source)
                transformed_coordinates = self.transform_coordinates(source)
                transformed_coordinates_name = self.source_name
                transformed_coordinates_colour = source_colour
                distance_threshold = self.destination_distance_threshold
                show_destination = True
            else:
                all_transformed_coordinates = self.transform_coordinates(self.destination, inverse=True)
                transformed_coordinates = self.transform_coordinates(destination, inverse=True)
                transformed_coordinates_name = self.destination_name
                transformed_coordinates_colour = destination_colour
                distance_threshold = self.source_distance_threshold
                show_source = True

            if distance_threshold > 0 and use_distance_threshold:
                plot_circles(axis, transformed_coordinates, radius=distance_threshold, linewidth=1,
                             facecolor='none', edgecolor=transformed_coordinates_colour)
                if show_pairs:
                    plot_circles(axis, all_transformed_coordinates[self.matched_pairs[:,int(inverse)]],
                                 radius=distance_threshold, linewidth=1,
                                 facecolor='none', edgecolor=pair_colour)
                # axis.plot(*transformed_coordinates.T, markerfacecolor='none', markeredgecolor=transformed_coordinates_colour,
                #           markeredgewidth=1, marker='o', linestyle='None', markersize=distance_threshold,
                #           label=f'{transformed_coordinates_name} transformed ({transformed_coordinates.shape[0]})')
                # if show_pairs:
                #     axis.plot(*all_transformed_coordinates[self.matched_pairs[:, int(inverse)]].T, markerfacecolor='none',
                #               markeredgecolor=pair_colour, markeredgewidth=1, marker='o', linestyle='None', markertransform=t,
                #               markersize=distance_threshold,)
            else:
                axis.plot(*transformed_coordinates.T, markerfacecolor='none', markeredgecolor=transformed_coordinates_colour,
                          markeredgewidth=1, marker='o', linestyle='None',
                          label=f'{transformed_coordinates_name} transformed ({transformed_coordinates.shape[0]})')
                if show_pairs:
                    axis.plot(*all_transformed_coordinates[self.matched_pairs[:, int(inverse)]].T, markerfacecolor='none',
                                 markeredgecolor=pair_colour, markeredgewidth=1, marker='o', linestyle='None')
        # else:
        #     show_source = True
        #     show_destination = True

        if show_source:
            axis.plot(*source.T, markerfacecolor=source_colour, markeredgecolor='none', marker='.', linestyle='None',
                         label=f'{self.source_name} ({source.shape[0]})')
            if show_pairs:
                axis.plot(*self.source[self.matched_pairs[:,0]].T, markerfacecolor=pair_colour, markeredgecolor='none',
                          marker='.', linestyle='None')

        if show_destination:
            axis.plot(*destination.T, markerfacecolor=destination_colour, markeredgecolor='none', marker='.',
                         linestyle='None', label=f'{self.destination_name} ({destination.shape[0]})')
            if show_pairs:
                axis.plot(*self.destination[self.matched_pairs[:,1]].T, markerfacecolor=pair_colour,
                          markeredgecolor='none', marker='.', linestyle='None')

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
        if show_transformed_coordinates:
            transformed_coordinates_marker = mlines.Line2D([], [], linewidth=0, markerfacecolor='none',
                                             markeredgecolor=transformed_coordinates_colour, marker='o')
            legend_dict[f'{transformed_coordinates_name} transformed ({transformed_coordinates.shape[0]})'] = transformed_coordinates_marker

        if show_pairs:
            pair_marker1 = mlines.Line2D([], [], linewidth=0, markerfacecolor='none',
                                         markeredgecolor=pair_colour, marker='o')
            pair_marker2 = mlines.Line2D([], [], linewidth=0, markerfacecolor=pair_colour,
                                         markeredgecolor='none', marker='.')
            legend_dict[f'matched pairs ({self.number_of_matched_points})'] = (pair_marker1, pair_marker2)
            # .get_legend().get_lines()[1].set_markerfacecolor('g')
        if show_legend:
            axis.legend(legend_dict.values(), legend_dict.keys(), loc='upper right')

        if save:
            if save_path is None:
                save_path = self.save_path
            save_path = Path(save_path)
            figure.savefig(save_path.joinpath(self.name+'.png'), bbox_inches='tight', dpi=250)

        if return_plot:
            return figure, axis

    def show_outline(self, inverse=False, source_colour='forestgreen', destination_colour='r', axis=None):
        """Show the outline or border of the point sets.

        Parameters
        ----------
        inverse : bool
            Indicates the direction of the transformation to use.
            If False: the source is transformed to destination space.
            If True: the destination is transformed to the source space.
        source_colour : str
            Colour used for the source points.
        destination_colour : str
            Colour used for the destination points.
        axis : matplotlib.axis.Axis
            The axis to use for plotting. If no axis is given a new figure and axis are created.

        Returns
        -------
        figure : matplotlib.figure.Figure
            Figure of the plot
        axis : matplotlib.axes.Axes
            Axis of the plot
        """
        if axis is None:
            figure, axis = plt.subplots()
        else:
            figure = axis.figure

        if inverse:
            space = 'source'
        else:
            space = 'destination'

        source_vertices = self.get_source_vertices(space=space)
        destination_vertices = self.get_destination_vertices(space=space)

        source_vertices = np.vstack([source_vertices, source_vertices[0]])
        destination_vertices = np.vstack([destination_vertices, destination_vertices[0]])

        axis.plot(*source_vertices.T, c=source_colour)
        axis.plot(*destination_vertices.T, c=destination_colour)

        axis.set_aspect('equal')

        return figure, axis

    def get_transformation_direction(self, direction):
        """Get inverse parameter based on direction

        Parameters
        ----------
        direction : str
           The direction of the transformation, either '<source_name>2<destination_name>' or
            '<destination_name>2<source_name>'. 'source2destination' and 'destination2source' can also be used.

        Returns
        -------
        inverse : bool
            Specifier for the use of the forward or inverse transformation.
        """

        if direction in [self.source_name + '2' + self.destination_name, 'source2destination']:
            inverse = False
        elif direction in [self.destination_name + '2' + self.source_name, 'destination2source']:
            inverse = True
        else:
            raise ValueError('Wrong direction')

        return inverse

    def calculate_inverse_transformation(self):
        #TODO: Make transformation a property and automatically calculate the inverse transformation.
        """Calculates the reverse transformation from the forward transformation.

        Sets the object's transformation_inverse attribute.

        Note
        ----
        This is only possible for linear transformations.
        """
        self.transformation_inverse = type(self.transformation)(matrix=self.transformation._inv_matrix)

    def transform_coordinates(self, coordinates, inverse=False, direction=None):
        """Transform coordinates using the object's transformation.

        Parameters
        ----------
        coordinates : Nx2 numpy.ndarray
            Coordinates to be transformed
        inverse : bool
            If True: the inverse transformation will be used (i.e. from destination to source)
        direction : str
            Another way of specifying the direction of the transformation, choose 'source2destination',
            'destination2source' or use the object's source and destination names instead of 'source' and 'destination'.

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
        """Transform an image using the objects transformation.

        Parameters
        ----------
        image : NxM numpy.ndarray
            Image to be transformed.
        inverse : bool
            If True: the inverse transformation will be used (i.e. from destination to source)
        direction : str
            Another way of specifying the direction of the transformation, choose 'source2destination',
            'destination2source' or use the object's source and destination names instead of 'source' and 'destination'.

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

    def distance_matrix(self, crop=True, space='destination', margin=None, max_distance=None, **kwargs):
        """Distance matrix for all point combinations of the source and destination.

        Parameters
        ----------
        crop : int or float, optional
            If True: the overlap between the source and destination point set is used.
            If 'source': only the source is cropped to the destination point set.
            If 'destination': only the destination is cropped to the source point set.
            If False: no cropping is applied.
        space : str, optional
            In which coordinate space to obtain the distances. Either 'source' or 'destination'.
        margin : float or tuple of floats
            Margins to exclude for the point sets.
            A tuple should be given in the format (margin_source, margin_destination).
            If a single margin is passed, then it will be used for both the source and destination.
        max_distance : float
            Maximum distance to include in the distance matrix.
            If given, then a sparse distance matrix will be returned.
        kwargs
            Keyword arguments to pass to the scipy.spatial.cKDTree.distance_matrix or
            scipy.spatial.cKDTree.sparse_distance_matrix function.

        Returns
        -------
        distance_matrix : MxN numpy.ndarray or scipy.sparse.spmatrix
            Distance matrix.
        """
        if margin is not None and len(margin) == 2:
            margin_source, margin_destination = margin
        else:
            margin_source = margin_destination = margin
        source = self.get_source(crop=crop, space=space, margin=margin_source)
        destination = self.get_destination(crop=crop, space=space, margin=margin_destination)

        if max_distance is None:
            return distance_matrix(source, destination, **kwargs)
        else:
            source_tree = cKDTree(source)
            destination_tree = cKDTree(destination)
            return source_tree.sparse_distance_matrix(destination_tree, max_distance=max_distance, **kwargs)

    def density_source(self, crop=False, space='source'):
        """Source density estimation from the number of points and the area contained within the convex hull.

        Parameters
        ----------
         crop : int or float, optional
            If True: the overlap between the source and destination point set is used.
            If 'source': only the source is cropped to the destination point set.
            If 'destination': only the destination is cropped to the source point set.
            If False: no cropping is applied.
        space : str, optional
            In which coordinate space to obtain the source density. Either 'source' or 'destination'.

        Returns
        -------
        density_source : float
            Density of the source dataset
        """
        return self.get_source(crop).shape[0] / self.get_source_area(crop=crop, space=space)

    def density_destination(self, crop=False, space='destination'):
        """Destination density estimation from the number of points and the area contained within the convex hull.

        Parameters
        ----------
         crop : int or float, optional
            If True: the overlap between the source and destination point set is used.
            If 'source': only the source is cropped to the destination point set.
            If 'destination': only the destination is cropped to the source point set.
            If False: no cropping is applied.
        space : str, optional
            In which coordinate space to obtain the destination density. Either 'source' or 'destination'.

        Returns
        -------
        density_source : float
            Density of the destination dataset
        """
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

    def save(self, save_path=None, filetype='nc'):
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
        save_path = Path(save_path)
        if not save_path.is_dir(): # save_path.suffix != '':
            self.name = save_path.with_suffix('').name
            save_path = save_path.parent
        if self.name is None:
            self.name = 'mapping'

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
                    attributes[key] = value.tolist()
                else:
                    attributes.pop(key)

            if filetype in ['yml','yaml']:
                with filepath.with_suffix('.mapping').open('w') as yml_file:
                    yaml.dump(attributes, yml_file, sort_keys=False)
            elif filetype == 'json':
                with filepath.with_suffix('.mapping').open('w') as json_file:
                    json.dump(attributes, json_file, sort_keys=False)

        elif filetype == 'nc':
            import xarray as xr
            ds = xr.Dataset()
            attributes = self.__dict__.copy()

            for key in list(attributes.keys()):
                value = attributes[key]
                if type(value) is str or np.issubdtype(type(value), np.number):
                    ds.attrs[key] = value
                elif isinstance(value, skimage.transform._geometric.GeometricTransform):
                    ds[key] = (('transformation_dim_0', 'transformation_dim_1'), value.params)
                elif type(value).__module__ == np.__name__:
                    if key in ['source', 'destination', 'matched_pairs']:
                        ds[key] = ((key + '_index', 'dimension'), value.reshape(-1,2))
                    else:
                        print(key)
                    # ds[key] = value

            ds.to_netcdf(filepath.with_suffix('.nc'), 'w')

        self.save_path = save_path

    def transformation_is_similar_to_correct_transformation(self, **kwargs):
        return is_similar_transformation(self.transformation, self.transformation_correct, **kwargs)


def compare_objects(object1, object2):
    for key1, value1 in object1.__dict__.items():
        if hasattr(object2, key1):
            if isinstance(value1, np.ndarray):
                is_equal = np.all(value1 == getattr(object2, key1))
            elif isinstance(value1, skimage.transform._geometric.GeometricTransform):
                is_equal = np.all(value1.params == getattr(object2, key1).params)
            else:
                is_equal = value1 == getattr(object2, key1)
        else:
            is_equal = False

        if not is_equal:
            return False
    return True


def is_similar_transformation(transformation1, transformation2, translation_error, rotation_error, scale_error):
    translation_check = (np.abs(transformation1.translation - transformation2.translation) < translation_error).all()
    rotation_check = np.abs(transformation1.rotation - transformation2.rotation) < rotation_error
    scale_check = (np.abs(np.array(transformation1.scale) - np.array(transformation2.scale)) < scale_error).all()
    return translation_check & rotation_check & scale_check

import scipy.sparse

#
# def singly_matched_pairs_within_radius(distance_matrix_, distance_threshold):
#     matches = distance_matrix_ < distance_threshold
#     sum_1 = matches.sum(axis=1) != 1
#     sum_0 = matches.sum(axis=0) != 1
#     matches[sum_1, :] = False
#     matches[:, sum_0] = False
#     return np.asarray(np.where(matches)).T


def matches_within_radius(distance_matrix_, distance_threshold, matches_per_point=1):
    matches = distance_matrix_ < distance_threshold
    matched_source_reference = matches.sum(axis=1) == matches_per_point
    matched_destination_reference = matches.sum(axis=0) == matches_per_point
    return matches, matched_source_reference, matched_destination_reference

def number_of_matches_within_radius(distance_matrix_, distance_threshold, matches_per_point=1):
    matches, matched_source_reference, matched_destination_reference = \
        matches_within_radius(distance_matrix_, distance_threshold, matches_per_point=matches_per_point)
    return matched_source_reference.sum(), matched_destination_reference.sum()

def singly_matched_pairs_within_radius(distance_matrix_, distance_threshold, point_set_name='all'):
    matches, matched_source, matched_destination = matches_within_radius(distance_matrix_, distance_threshold)
    if point_set_name in ['source', 'all']:
        matches[~matched_source, :] = False
    if point_set_name in ['destination', 'all']:
        matches[:, ~matched_destination] = False
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
