# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 13:55:53 2019

@author: ivoseverins
"""

import numpy as np
import matplotlib.pyplot as plt

from trace_analysis.mapping.icp import icp, nearest_neighbor_pair
from trace_analysis.mapping.icp_nonrigid import icp_nonrigid
from trace_analysis.coordinate_transformations import transform
from trace_analysis.plotting import plot_match

class Mapping2:
    def __init__(self, source = None, destination = None, method = None,
                 transformation_type = 'linear', initial_translation = None):
        self.source = source
        self.destination = destination
        self.method = method
        self.transformation_type = transformation_type
        self._transformation = None
        self.transformation_inverse = None
        self.initial_translation = initial_translation

        # if (source is not None) and (destination is not None):
        #     if self.method is None: self.method = 'icp'
        #     self.perform_mapping()

    @property
    def translation(self):
        return self.transformation[0:2,2]

    @property
    def magnification(self):
        return np.linalg.norm(self.transformation[:,0:2],axis=0)

    @property
    def rotation(self):
        rotation_matrix = self.transformation[0:2, 0:2]/self.magnification[0]
        return np.arctan2(rotation_matrix[0,1]-rotation_matrix[1,0],rotation_matrix[0,0]+rotation_matrix[1,1])/(2*np.pi)*360

    @property
    def transform_source_to_destination(self):
        return self.transform_coordinates(self.source)

    @property
    def transformation(self):
        return self._transformation

    @transformation.setter
    def transformation(self, transformation):
        self._transformation = transformation
        self.transformation_inverse = np.linalg.inv(self.transformation)

    def perform_mapping(self):
        if self.method == 'icp':
            self.transformation, distances, iterations = icp(self.source, self.destination, initial_translation=self.initial_translation)
        elif self.method == 'icp-non-rigid':
            self.transformation, distances, iterations = icp_nonrigid(self.source, self.destination)
        # elif method == 'manual'         : mapping_manual(source, destination)
        # elif method == 'automatic'      : mapping_automatic(source, destination)
        else: print('Method not found')

    def nearest_neighbour_match(self, distance_threshold = 1):
        destination_from_source = self.transform_coordinates(self.source)
        destination = self.destination
        distances, source_indices, destination_indices = nearest_neighbor_pair(destination_from_source,destination)

        destination_from_source = np.hstack([destination_from_source, np.ones((len(destination_from_source), 1))])
        destination = np.hstack([destination, np.ones((len(destination), 1))])

        T, res, rank, s = np.linalg.lstsq(destination_from_source[source_indices[distances<distance_threshold]],
                                          destination[destination_indices[distances<distance_threshold]], rcond=None)
        transformation = T.T

        self.transformation = transformation @ self.transformation

        # new_destination_from_source = transform(destination_from_source[:,0:2], transformationMatrix=transformation)
        #
        # figure = plt.figure()
        # axis = figure.gca()
        # #axis.scatter(self.source[source_indices, 0], self.source[:, 1], c='g')
        # axis.scatter(self.destination[destination_indices, 0], self.destination[destination_indices, 1], c='r')
        # axis.scatter(destination_from_source[source_indices, 0], destination_from_source[source_indices, 1], c='g')
        # axis.scatter(new_destination_from_source[source_indices, 0], new_destination_from_source[source_indices, 1], c='b')

    def show_mapping_transformation(self, figure=None):
        if not figure: figure = plt.figure()

        destination_from_source = self.transform_coordinates(self.source)

        axis = figure.gca()

        axis.scatter(self.source[:, 0], self.source[:, 1], c='g')
        axis.scatter(self.destination[:, 0], self.destination[:, 1], c='r')
        axis.scatter(destination_from_source[:, 0], destination_from_source[:, 1], c='g')

    def transform_coordinates(self, coordinates, inverse = False):
        if self.transformation_type == 'linear':
            if inverse is False: return transform(coordinates, self.transformation)
            elif inverse is True: return transform(coordinates, self.transformation_inverse)
        else: print('Transformation not found')


# from trace_analysis.icp_nonrigid import icp_nonrigid
# from trace_analysis.image_adapt.polywarp import polywarp_apply
# kx, ky, distances, i = icp_nonrigid(donor,acceptor, tolerance=0.00000001, max_iterations=50)
# acceptor_calculated = polywarp_apply(kx, ky, donor)
