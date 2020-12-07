# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 13:55:53 2019

@author: ivoseverins
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import yaml
import skimage.transform
from skimage.transform import AffineTransform, PolynomialTransform

from trace_analysis.mapping.icp import icp, nearest_neighbor_pair
from trace_analysis.image_adapt.polywarp import polywarp, polywarp_apply #required for nonlinear
from trace_analysis.mapping.polywarp import PolywarpTransform

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

    def __init__(self, source=None, destination=None, method=None,
                 transformation_type=None, initial_transformation=None, load=None):
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
        initial_transformation : AffineTransform or PolynomialTransform
            Initial transformation to perform as a starting point for the mapping algorithms
        load : str or pathlib.Path
            Path to file that has to be loaded
        """

        self.source_name = 'source'
        self.source = source #source=donor=left side image
        self.destination_name = 'destination'
        self.destination = destination #destination=acceptor=right side image
        self.method = method
        self.transformation_type = transformation_type
        self.initial_transformation = initial_transformation
        self.transformation = None
        self.transformation_inverse = None

        # if (source is not None) and (destination is not None):
        #     if self.method is None: self.method = 'icp'
        #     self.perform_mapping()

        if load:
            filepath = Path(load)
            if filepath.suffix in ['.yml','.yaml']:
                with filepath.open('r') as yml_file:
                    attributes = yaml.load(yml_file, Loader=yaml.SafeLoader)
                    for key, value in attributes.items():
                        if type(value) == list:
                            value = np.array(value)
                        setattr(self, key, value)
                #TODO: Following if else statement can be simplified
                if self.transformation_type == 'linear':
                    self.transformation = AffineTransform(self.transformation)
                    self.transformation_inverse = AffineTransform(self.transformation_inverse)
                #TODO: make the full class using PolywarpTransform before uncommenting this
                # elif self.transformation_type == 'nonlinear':
                #     transformation = PolywarpTransform()
                #     transformation.params = self.transformation
                #     transformation_inverse = PolywarpTransform()
                #     transformation_inverse.params = self.transformation_inverse
                #     self.transformation = transformation
                #     self.transformation_inverse = transformation_inverse
                elif self.transformation_type == 'polynomial':
                    transformation = PolynomialTransform()
                    transformation.params = self.transformation
                    transformation_inverse = PolynomialTransform()
                    transformation_inverse.params = self.transformation_inverse
                    self.transformation = transformation
                    self.transformation_inverse = transformation_inverse

    def __getattr__(self, item):
        if hasattr(self.transformation, item):
            return(getattr(self.transformation, item))
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

    def perform_mapping(self, method=None):
        """Find transformation from source to destination points using one of the mapping methods

        Parameters
        ----------
        method : str
            Mapping method, if not specified the object method is used.

        """

        if method is None:
            method = self.method

        if method in ['icp', 'iterative_closest_point']: #icp should be default
            self.iterative_closest_point()
        elif method in ['direct']:
            self.direct_match()
        elif method in ['nn', 'nearest_neighbour']:
            self.nearest_neighbour_match()
        else:
            raise ValueError('Method not found')

        self.method = method

    def direct_match(self, transformation_type=None):
        """Find transformation from source to destination points by matching based on the point order

        Note: the number and the order of source points should be equal to the number and the order of destination points.

        Parameters
        ----------
        transformation_type : str
            Type of transformation used, either linear or polynomial can be chosen.
            If not specified the object transformation_type is used.

        """
        if transformation_type is None:
            transformation_type = self.transformation_type

        self.transformation_type = transformation_type

        if transformation_type=='linear':
            source = np.hstack([self.source, np.ones((len(self.source), 1))])
            destination = np.hstack([self.destination, np.ones((len(self.destination), 1))])
            T, res, rank, s = np.linalg.lstsq(source, destination, rcond=None)
            self.transformation = AffineTransform(T.T)
            self.transformation_inverse = AffineTransform(self.transformation._inv_matrix)
        elif transformation_type=='nonlinear':
            self.transformation = PolywarpTransform()
            self.transformation.estimate(self.source, self.destination, order=3)
            self.transformation_inverse = PolywarpTransform()
            self.transformation_inverse.estimate(self.destination, self.source, order=3)
        elif transformation_type=='polynomial':
            self.transformation = PolynomialTransform()
            self.transformation.estimate(self.source, self.destination)
            self.transformation_inverse = PolynomialTransform()
            self.transformation_inverse.estimate(self.destination, self.source)

    def nearest_neighbour_match(self, distance_threshold=1, transformation_type=None):
        """Find transformation from source to destination points by matching nearest neighbours

        Two-way nearest neighbours are detectected, i.e. the source point should be the nearest neighbour of the
        destination point and vice versa. Only nearest neighbours closer than the distance threshold are used to find
        the transformation.

        Note
        ----
        The current transformation is first applied and then the nearest neighbour match is performed, basically to
        improve the current transformation.

        Parameters
        ----------
        distance_threshold : float
            Distance threshold for nearest neighbour match, i.e. nearest neighbours with a distance larger than the
            distance threshold are not used.
        transformation_type : str
            Type of transformation used, either linear or polynomial can be chosen.
            If not specified the object transformation_type is used.

        """
        if not transformation_type:
            transformation_type = self.transformation_type

        # TODO: Probably this can be partly merged with the icp function

        source = self.source
        destination_from_source = self.transform_coordinates(source)
        destination = self.destination
        distances, source_indices, destination_indices = nearest_neighbor_pair(destination_from_source,destination)

        source = np.hstack([source, np.ones((len(source), 1))])
        destination = np.hstack([destination, np.ones((len(destination), 1))])

        source_points_for_matching = source[source_indices[distances < distance_threshold]]
        destination_points_for_matching = destination[destination_indices[distances<distance_threshold]]

        if transformation_type == 'linear':
            T, res, rank, s = np.linalg.lstsq(source_points_for_matching, destination_points_for_matching, rcond=None)
            # transformation = T.T
            #
            # self.transformation = transformation @ self.transformation
            self.transformation = AffineTransform(T.T)
            self.transformation_inverse = AffineTransform(self.transformation._inv_matrix)

        elif transformation_type == 'nonlinear':
            self.transformation = PolywarpTransform()
            self.transformation.estimate(source_points_for_matching[:, 0:2], destination_points_for_matching[:, 0:2],
                                         order=3)
            self.transformation_inverse = PolywarpTransform()
            self.transformation_inverse.estimate(destination_points_for_matching[:, 0:2],
                                                 source_points_for_matching[:, 0:2], order=3)
        elif transformation_type == 'polynomial':
            self.transformation = PolynomialTransform()
            self.transformation.estimate(source_points_for_matching[:, 0:2], destination_points_for_matching[:, 0:2], order=3)
            self.transformation_inverse = PolynomialTransform()
            self.transformation_inverse.estimate(destination_points_for_matching[:, 0:2], source_points_for_matching[:, 0:2], order=3)

        self.transformation_type = transformation_type

        # new_destination_from_source = transform(destination_from_source[:,0:2], transformationMatrix=transformation)
        #
        # figure = plt.figure()
        # axis = figure.gca()
        # #axis.scatter(self.source[source_indices, 0], self.source[:, 1], c='g')
        # axis.scatter(self.destination[destination_indices, 0], self.destination[destination_indices, 1], c='r')
        # axis.scatter(destination_from_source[source_indices, 0], destination_from_source[source_indices, 1], c='g')
        # axis.scatter(new_destination_from_source[source_indices, 0], new_destination_from_source[source_indices, 1], c='b')

    def iterative_closest_point(self, **kwargs):
        """Find transformation from source to destination points using an iterative closest point algorithm

        In the iterative closest point algorithm, the two-way nearest neigbhours are found and these are used to
        find the most optimal transformation. Subsequently the source is transformed according to this
        transformation. This process is repeated until the changes detected are below a tolerance level.

        The iterative closest point algorithm can be used in situations when deviations between the two point sets
        are relatively small.

        Parameters
        ----------
        max_iterations : int
            Maximum number of iterations to perform
        tolerance : float
            If the mean distance between the two point sets is below this tolerance, then the convergence is assumed.
        use_cutoff : bool or int
            If an int is passed, then the value will be used as a distance threshold for nearest neighbour detection.
            If use_cutoff is set to True, then the median + std of the distances will be used as threshold.

        """

        # TODO: Let icp pass and return AffineTransformation or PolynomialTransformation
        self.transformation, distances, iterations, self.transformation_inverse = \
            icp(self.source, self.destination, initial_transformation=self.initial_transformation,
                transformation_type=self.transformation_type, **kwargs)
        if self.transformation_type == 'linear':
            self.transformation = AffineTransform(self.transformation)
            self.transformation_inverse = AffineTransform(self.transformation_inverse)
        elif self.transformation_type == 'nonlinear':
            # def polywarp_to_polynomial_transform(polywarp_output):
            #
            #     order = polywarp_output[0].shape[0]
            #
            #     polynomial_transform_params_list = []
            #
            #     for polywarp_matrix in polywarp_output:
            #         polynomial_transform_params_list.append(
            #             [polywarp_matrix[j - i, i] if ((j - i) < order) and (i < order) else 0 for j in
            #              range(order * 2 - 1) for
            #              i in range(j + 1)]
            #         )
            #
            #     polynomial_transform_params = np.vstack(polynomial_transform_params_list)
            #
            #     return PolynomialTransform(params=polynomial_transform_params)

            self.transformation = PolywarpTransform(self.transformation)
            self.transformation_inverse = PolywarpTransform(self.transformation_inverse)

    def number_of_matched_points(self, distance_threshold):
        """Number of matched points determined by finding the two-way nearest neigbours that are closer than a distance
        threshold.

        Parameters
        ----------
        distance_threshold : float
            Distance threshold for nearest neighbour match, i.e. nearest neighbours with a distance larger than the
            distance threshold are not counted.

        Returns
        -------
        int
            Number of matched points

        """
        distances, source_indices, destination_indices = nearest_neighbor_pair(self.transform_source_to_destination,
                                                                               self.destination)
        return np.sum(distances<distance_threshold)

    def show_mapping_transformation(self, figure=None, show_source=False):
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

        destination_from_source = self.transform_coordinates(self.source)

        axis = figure.gca()

        if show_source:
            axis.scatter(self.source[:, 0], self.source[:, 1], facecolors='forestgreen', edgecolors='none', marker='.')
        axis.scatter(self.destination[:, 0], self.destination[:, 1], facecolors='none', edgecolors='forestgreen', linewidth=1, marker='o')
        axis.scatter(destination_from_source[:, 0], destination_from_source[:, 1], facecolors='r', edgecolors='none', marker='.')

        axis.set_aspect('equal')
        axis.set_xlabel('x')
        axis.set_ylabel('y')

    def get_transformation_direction(self, direction=None):
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

        if direction is not None:
            if direction == self.source_name + '2' + self.destination_name:
                inverse = False
            elif direction == self.destination_name + '2' + self.source_name:
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

    def save(self, filepath, filetype='yaml'):
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

        elif filetype in ['yml', 'yaml']:
            attributes = self.__dict__.copy()
            for key, value in attributes.items():
                if type(value) in [str, int, float]:
                    continue
                if isinstance(value, skimage.transform._geometric.GeometricTransform):
                    value = value.params
                if type(value).__module__ == np.__name__:
                    attributes[key] = value.tolist()

            with filepath.with_suffix('.yml').open('w') as yml_file:
                yaml.dump(attributes, yml_file, sort_keys=False)


if __name__ == "__main__":
    # Create test point set (with some missing points)
    number_of_points = 40
    np.random.seed(32)
    source = np.random.rand(number_of_points, 2) * 1000
    transformation = AffineTransform(translation=[50,25], rotation=2/360*2*np.pi, scale=[1.1,1.1])
    destination = transformation(source)[5:]
    source = source[:35]

    # Make a mapping object, perform the mapping and show the transformation
    mapping = Mapping2(source, destination, transformation_type='linear')
    mapping.method = 'icp'
    mapping.perform_mapping()
    mapping.show_mapping_transformation(show_source=True)

    # Create a test image
    image = np.zeros((250,250))
    image[100:110,:]=1
    image[:,100:110]=1

    # Transform the image
    image_transformed = mapping.transform_image(image)

    # Show the original and transformed image
    plt.figure()
    plt.imshow(image)
    plt.figure()
    plt.imshow(image_transformed)
