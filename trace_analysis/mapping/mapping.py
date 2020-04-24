# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 13:55:53 2019

@author: ivoseverins
"""

import numpy as np
import matplotlib.pyplot as plt


from trace_analysis.mapping.icp import icp

from trace_analysis.coordinate_transformations import transform
from trace_analysis.image_adapt.polywarp import polywarp, polywarp_apply #required for nonlinear

class Mapping2:
    def __init__(self, source = None, destination = None, method = None,
                 transformation_type = 'linear', initial_translation = None):
        self.source_name = 'source'
        self.source = source #source=donor=left side image
        self.destination_name = 'destination'
        self.destination = destination #destination=acceptor=right side image
        self.method = method
        self.transformation_type = transformation_type
        self.initial_translation = initial_translation
        self.transformation = None
        self.transformation_inverse = None

        if (source is not None) and (destination is not None):
            if self.method is None: self.method = 'icp'
            self.perform_mapping()

#    @property
#    def transformation_inverse(self):
#        return np.linalg.inv(self.transformation)

    @property
    def transform_source_to_destination(self): 
        return self.transform_coordinates(self.source)

    def perform_mapping(self):
        print(self.transformation_type)
        if self.method == 'icp': #icp should be default
            self.transformation, distances, iterations, self.transformation_inverse = \
                icp(self.source, self.destination, initial_translation=self.initial_translation,
                    transformation_type=self.transformation_type)
        else: print('Method not found')

    def show_mapping_transformation(self, figure=None): 
        if not figure: figure = plt.figure()
        destination_from_source = self.transform_coordinates(self.source)

        axis = figure.gca()

        axis.scatter(self.source[:, 0], self.source[:, 1], c='g')
        axis.scatter(self.destination[:, 0], self.destination[:, 1], c='r')
        axis.scatter(destination_from_source[:, 0], destination_from_source[:, 1], c='y')

    def transform_coordinates(self, coordinates, direction='source2destination'):
        print(self.transformation, self.method)
        if self.method == 'icp':
#            return icp_apply_transform(coordinates, direction, self.transformation,self.transformation_inverse, self.transformation_type)
            coords_transformed= coordinates.copy()
        
            if self.transformation_type == 'linear': #$$$$$$$$$$$$$$$$$
                if direction == 'source2destination': return transform(coords_transformed[:,:2], self.transformation)
                elif direction == 'destination2source' : return transform(coords_transformed[:,:2], self.transformation_inverse)
        
            elif self.transformation_type == 'nonlinear':
                    if direction == 'source2destination' : return polywarp_apply(self.transformation[0],self.transformation[1],coords_transformed)
                    elif direction == 'destination2source' : return polywarp_apply(self.transformation_inverse[0],self.transformation_inverse[1],coords_transformed)
            
            return coords_transformed[:,:1]
        else: print('transform_coordinates only works for icp')

