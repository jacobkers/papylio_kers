# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 13:55:53 2019

@author: ivoseverins
"""
import json

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import yaml
import skimage.transform
from skimage.transform import AffineTransform, PolynomialTransform
import matplotlib.path as pth
from shapely.geometry import Polygon, MultiPoint

from trace_analysis.mapping.icp import icp, nearest_neighbor_pair, nearest_neighbour_match, direct_match
from trace_analysis.mapping.polywarp import PolywarpTransform
from trace_analysis.mapping.polynomial import PolynomialTransform

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

    transformation_types = {'linear': AffineTransform,
                            'nonlinear': PolywarpTransform,
                            'polynomial': PolynomialTransform}

    def __init__(self, source=None, destination=None, method=None,
                 transformation_type=None, initial_transformation=None,
                 source_name='source', destination_name='destination',
                 source_unit=None, destination_unit=None, name=None,
                 load=None):
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
        self.source_name = source_name
        self.source = source #source=donor=left side image
        self.source_unit = source_unit
        self.source_vertices = None # np.array(MultiPoint(self.source).envelope.exterior.xy).T[:-1]
        self.destination_name = destination_name
        self.destination_unit = destination_unit
        self.destination = destination #destination=acceptor=right side image
        self.destination_vertices = None # np.array(MultiPoint(self.destination).envelope.exterior.xy).T[:-1]
        self.method = method
        self.transformation_type = transformation_type

        if type(initial_transformation) is dict:
            initial_transformation = AffineTransform(**initial_transformation)
        self.initial_transformation = initial_transformation

        self.transformation = None
        self.transformation_inverse = None

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
                self.transformation = self.transform(self.transformation)
                self.transformation_inverse = self.transform(self.transformation_inverse)

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
    def source_to_destination(self):
        """Nx2 numpy.ndarray : Source coordinates transformed to the destination axis"""

        return self.transform_coordinates(self.source)

    @property
    def transform(self):
        """type : Transform class based on the set transformation_type"""

        return self.transformation_types[self.transformation_type]

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

    def cross_correlation(self, peak_detection='auto'):
        from trace_analysis.mapping.cross_correlation import cross_correlate
        if self.transformation is None:
            self.transformation = AffineTransform()
            self.transformation_inverse = AffineTransform()
        correlation, self.correlation_conversion_function = cross_correlate(self.source_to_destination, self.destination)

        if peak_detection == 'auto':
            correlation_peak_coordinates = np.array(np.where(correlation==correlation.max())).flatten()[::-1]
            self.set_correlation_peak_coordinates(correlation_peak_coordinates)

    def set_correlation_peak_coordinates(self, correlation_peak_coordinates):
        if not hasattr(self, 'correlation_conversion_function') and self.correlation_conversion_function is not None:
            raise RuntimeError('Run cross_correlation first')
        transformation = self.correlation_conversion_function(correlation_peak_coordinates) # is this the correct direction
        self.transformation += transformation
        self.correlation_conversion_function = None

    def number_of_matched_points(self, distance_threshold):
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
        distances, source_indices, destination_indices = \
            nearest_neighbor_pair(self.source_to_destination, self.destination)

        return np.sum(distances < distance_threshold)

    def show_mapping_transformation(self, figure=None, show_source=False, show_destination=False, crop=None,
                                    inverse=False, source_colour='forestgreen', destination_colour='r', save_path=None):
        """Show a point scatter of the source transformed to the destination points and the destination.

        Parameters
        ----------
        figure : matplotlib.figure.Figure
            If figure is passed then the scatter plot will be made in the currently active axis.
            If no figure is passed a new figure will be made.
        show_source : bool
            If True then the source points are also displayed

        """
        if not figure:
            figure = plt.figure()

        if not crop:
            source = self.source
            destination = self.destination
        else:
            source = self.source_cropped
            destination = self.destination_cropped

        axis = figure.gca()

        if not inverse:
            destination_from_source = self.transform_coordinates(source)
            axis.scatter(destination_from_source[:, 0], destination_from_source[:, 1], facecolors='none',
                         edgecolors=source_colour, linewidth=1, marker='o', label=f'{self.source_name} transformed')
            show_destination = True
        else:
            source_from_destination = self.transform_coordinates(destination, inverse=True)
            axis.scatter(source_from_destination[:, 0], source_from_destination[:, 1], facecolors='none',
                         edgecolors=destination_colour, linewidth=1, marker='o', label=f'{self.destination_name} transformed')
            show_source = True

        if show_source:
            axis.scatter(source[:, 0], source[:, 1], facecolors=source_colour, edgecolors='none', marker='.',
                         label=self.source_name)
            # axis.scatter(self.source[:, 0], self.source[:, 1], facecolors='none', edgecolors='g', s=70,
            #            label=self.source_name)
        # axis.scatter(destination_from_source[:, 0], destination_from_source[:, 1], facecolors='none',
        #              edgecolors='forestgreen', linewidth=1, marker='o', label=f'{self.source_name} transformed')
        # axis.scatter(self.destination[:, 0], self.destination[:, 1], facecolors='r', edgecolors='none', marker='.',
        #              label=self.destination_name)

        # axis.scatter(destination_from_source[:, 0], destination_from_source[:, 1], facecolors='none',
        #              edgecolors='forestgreen', linewidth=1, marker='o', label=f'{self.source_name} transformed')
        # axis.scatter(self.destination[:, 0], self.destination[:, 1], facecolors='r', edgecolors='none', marker='.',
        #              label=self.destination_name)
        # axis.scatter(destination_from_source[:, 0], destination_from_source[:, 1], facecolors='none',
        #              edgecolors='g', linewidth=1, marker='o', label=f'{self.source_name} transformed', s=70)
        # axis.scatter(self.destination[:, 0], self.destination[:, 1], facecolors='none', edgecolors='r', marker='o',
        #              label=self.destination_name, alpha=.5, s=60)
        if show_destination:
            axis.scatter(destination[:, 0], destination[:, 1], facecolors=destination_colour, edgecolors='none', marker='.',
                         label=self.destination_name)

        axis.set_aspect('equal')

        if self.destination_unit is not None:
            unit_label = f' ({self.destination_unit})'
        else:
            unit_label = ''

        axis.set_title(self.name)
        axis.set_xlabel('x'+unit_label)
        axis.set_ylabel('y'+unit_label)
        axis.legend()

        if save_path is not None:
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

    @property
    def source_vertices_in_destination(self):
        if self.source_vertices is not None:
            return self.transform_coordinates(self.source_vertices)
        else:
            raise AttributeError('Source vertices not set')

    @property
    def destination_vertices_in_source(self):
        if self.destination_vertices is not None:
            return self.transform_coordinates(self.destination_vertices, inverse=True)
        else:
            raise AttributeError('Destination vertices not set')

    @property
    def source_cropped(self):
        crop_vertices_in_source = overlap_vertices(self.source_vertices, self.destination_vertices_in_source)
        return crop_coordinates(self.source, crop_vertices_in_source)

    @property
    def destination_cropped(self):
        crop_vertices_in_destination = overlap_vertices(self.source_vertices_in_destination, self.destination_vertices)
        return crop_coordinates(self.destination, crop_vertices_in_destination)

    def fraction_of_points_matched(self, distance_threshold):
        number_of_matched_points = self.number_of_matched_points(distance_threshold)
        fraction_source_matched = number_of_matched_points / self.source_cropped.shape[0]
        fraction_destination_matched = number_of_matched_points / self.destination_cropped.shape[0]

        return fraction_source_matched, fraction_destination_matched

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

    def save(self, filepath, filetype='json'):
        """Save the current mapping in a file, so that it can be opened later.

        Parameters
        ----------
        filepath : str or pathlib.Path
            Path to file (including filename)
        filetype : str
            Choose classic to export a .coeff file for a linear transformation
            Choose yml to export all object attributes in a yml text file

        """
        filepath = Path(filepath)
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


def overlap_vertices(vertices_A, vertices_B):
    polygon_A = Polygon(vertices_A)
    polygon_B = Polygon(vertices_B)
    polygon_overlap = polygon_A.intersection(polygon_B)
    return np.array(polygon_overlap.exterior.coords.xy).T[:-1]

def crop_coordinates_indices(coordinates, vertices):
    return pth.Path(vertices).contains_points(coordinates)

def crop_coordinates(coordinates, vertices):
    indices = crop_coordinates_indices(coordinates, vertices)
    return coordinates[indices]

    # bounds.sort(axis=0)
    # selection = (coordinates[:, 0] > bounds[0, 0]) & (coordinates[:, 0] < bounds[1, 0]) & \
    #             (coordinates[:, 1] > bounds[0, 1]) & (coordinates[:, 1] < bounds[1, 1])
    # return coordinates[selection]


if __name__ == "__main__":
    # Create test point set (with some missing points)
    from trace_analysis.plugins.sequencing.point_set_simulation import simulate_mapping_test_point_set

    # Simulate source and destination point sets
    number_of_source_points = 40
    transformation = AffineTransform(translation=[256, 25])#, rotation=5 / 360 * 2 * np.pi, scale=[0.98, 0.98])
    source_bounds = ([0,0], [256, 512])
    source_crop_bounds = [[50,150], [200,300]]
    fraction_missing_source = 0
    fraction_missing_destination = 0
    maximum_error_source = 0
    maximum_error_destination = 0
    shuffle = True

    source, destination = simulate_mapping_test_point_set(number_of_source_points, transformation,
                                                          source_bounds, source_crop_bounds,
                                                          fraction_missing_source, fraction_missing_destination,
                                                          maximum_error_source, maximum_error_destination, shuffle)

    # Make a mapping object, perform the mapping and show the transformation
    mapping = Mapping2(source, destination, transformation_type='linear')
    #mapping.method = 'icp'
    #mapping.perform_mapping()
    mapping.cross_correlation()
    # correlation_peak_coordinates = [244, 487]
    # mapping.set_correlation_peak_coordinates(correlation_peak_coordinates)
    mapping.show_mapping_transformation(show_source=True)

    # # Create a test image
    # image = np.zeros((512,512))
    # image[100:110,:]=1
    # image[:,100:110]=1
    #
    # # Transform the image
    # image_transformed = mapping.transform_image(image)
    #
    # # Show the original and transformed image
    # plt.figure()
    # plt.imshow(image)
    # plt.figure()
    # plt.imshow(image_transformed)
