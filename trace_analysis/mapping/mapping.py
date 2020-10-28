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

from trace_analysis.mapping.icp import icp, nearest_neighbor_pair

from trace_analysis.coordinate_transformations import transform
from trace_analysis.image_adapt.polywarp import polywarp, polywarp_apply #required for nonlinear

class Mapping2:
    def __init__(self, source = None, destination = None, method = None,
                 transformation_type = None, initial_transformation = None, load = None):

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

    @property
    def translation(self):
        return self.transformation[0:2,2]

    @property
    def magnification(self):
        return np.linalg.norm(self.transformation[:,0:2],axis=0)

    @property
    def rotation(self):
        # rotation_matrix = self.transformation[0:2, 0:2]/self.magnification[0]
        # return np.arctan2(rotation_matrix[0,1]-rotation_matrix[1,0],rotation_matrix[0,0]+rotation_matrix[1,1])/(2*np.pi)*360
        return np.arctan2(self.transformation[0, 1], self.transformation[0, 0]) / (2 * np.pi) * 360

    @property
    def reflection(self):
        return np.array([np.sign(self.transformation[0, 0]), np.sign(self.transformation[1, 1])])

    @property
    def transform_source_to_destination(self): 
        return self.transform_coordinates(self.source)

    def calculate_inverse_transformation(self):
        if self.transformation_type == 'linear':
            self.transformation_inverse = np.linalg.inv(self.transformation)
        else:
            raise RuntimeError('Inverse transformation cannot be determined. Transformation type must be linear')

    def perform_mapping(self):
        print(self.transformation_type)
        if self.method == 'icp': #icp should be default
            self.transformation, distances, iterations, self.transformation_inverse = \
                icp(self.source, self.destination, initial_transformation=self.initial_transformation,
                    transformation_type=self.transformation_type)
        else: print('Method not found')

    def direct_match(self, transformation_type=None):
        if not transformation_type:
            transformation_type = self.transformation_type

        self.transformation_type = transformation_type

        if transformation_type=='linear':
            source = np.hstack([self.source, np.ones((len(self.source), 1))])
            destination = np.hstack([self.destination, np.ones((len(self.destination), 1))])
            T, res, rank, s = np.linalg.lstsq(source, destination, rcond=None)
            self.transformation = T.T
            self.calculate_inverse_transformation()
        elif transformation_type=='nonlinear':
            kx, ky = polywarp(self.destination, self.source)
            kx_inv, ky_inv = polywarp(self.source, self.destination)
            self.transformation = (kx, ky)
            self.transformation_inverse = (kx_inv, ky_inv)

    def nearest_neighbour_match(self, distance_threshold=1, transformation_type=None):
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

        if transformation_type=='linear':
            T, res, rank, s = np.linalg.lstsq(source_points_for_matching, destination_points_for_matching, rcond=None)
            # transformation = T.T
            #
            # self.transformation = transformation @ self.transformation
            self.transformation = T.T
            self.calculate_inverse_transformation()

        elif transformation_type=='nonlinear':
            kx, ky = polywarp(destination_points_for_matching[:, 0:2], source_points_for_matching[:, 0:2])
            kx_inv, ky_inv = polywarp(source_points_for_matching[:, 0:2], destination_points_for_matching[:, 0:2])
            self.transformation = (kx, ky)
            self.transformation_inverse = (kx_inv, ky_inv)

        self.transformation_type = transformation_type

        # new_destination_from_source = transform(destination_from_source[:,0:2], transformationMatrix=transformation)
        #
        # figure = plt.figure()
        # axis = figure.gca()
        # #axis.scatter(self.source[source_indices, 0], self.source[:, 1], c='g')
        # axis.scatter(self.destination[destination_indices, 0], self.destination[destination_indices, 1], c='r')
        # axis.scatter(destination_from_source[source_indices, 0], destination_from_source[source_indices, 1], c='g')
        # axis.scatter(new_destination_from_source[source_indices, 0], new_destination_from_source[source_indices, 1], c='b')

    def number_of_matched_points(self, distance_threshold):
        distances, source_indices, destination_indices = nearest_neighbor_pair(self.transform_source_to_destination,
                                                                               self.destination)
        return np.sum(distances<distance_threshold)

    def show_mapping_transformation(self, figure=None):
        if not figure: figure = plt.figure()

        destination_from_source = self.transform_coordinates(self.source)

        axis = figure.gca()

        axis.scatter(self.source[:, 0], self.source[:, 1], c='g')
        axis.scatter(self.destination[:, 0], self.destination[:, 1], c='r')
        axis.scatter(destination_from_source[:, 0], destination_from_source[:, 1], c='y')

    def transform_coordinates(self, coordinates, inverse=False, direction=None):
        if direction is not None:
            if direction == self.source_name + '2' + self.destination_name:
                inverse = False
            elif direction == self.destination_name + '2' + self.source_name:
                inverse = True
            else:
                raise ValueError('Wrong direction')

        if not inverse:
            current_transformation = self.transformation
        else:
            current_transformation = self.transformation_inverse

        if self.transformation_type == 'linear':
            # Maybe we should rename transform to linear_transform [IS 05-03-2020]
            # transform(pointSet, transformationMatrix=None, **kwargs):
            return transform(coordinates[:, :2], current_transformation)
        elif self.transformation_type == 'nonlinear':
            return polywarp_apply(current_transformation[0], current_transformation[1], coordinates)

    def save(self, filepath, filetype='yaml'):
        filepath = Path(filepath)
        if filetype == 'classic':
            if self.transformation_type == 'linear':
                coeff_filepath = filepath.with_suffix('.coeff')
                coefficients = self.mapping.transformation[[0, 0, 0, 1, 1, 1], [2, 0, 1, 2, 0, 1]]
                np.savetxt(coeff_filepath, coefficients, fmt='%13.6g') # Same format used as in IDL code
            else:
                raise TypeError('Mapping is not of type linear')
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


# from trace_analysis.icp_nonrigid import icp_nonrigid
# from trace_analysis.image_adapt.polywarp import polywarp_apply
# kx, ky, distances, i = icp_nonrigid(donor,acceptor, tolerance=0.00000001, max_iterations=50)
# acceptor_calculated = polywarp_apply(kx, ky, donor)

