# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 13:55:53 2019

@author: ivoseverins
"""

import matplotlib.pyplot as plt

from trace_analysis.mapping.icp import icp
from trace_analysis.mapping.icp_nonrigid import icp_nonrigid
from trace_analysis.coordinate_transformations import transform

class Mapping2:
    def __init__(self, source = None, destination = None, method = None,
                 transformation_type = 'linear'):
        self.source = source
        self.destination = destination
        self.method = method
        self.transformation_type = transformation_type
        self.transformation = None

        if (source is not None) and (destination is not None):
            if self.method is None: self.method = 'icp'
            self.perform_mapping()

    def perform_mapping(self):
        if self.method == 'icp':
            self.transformation, distances, iterations = icp(self.source, self.destination)
        elif self.method == 'icp-non-rigid':
            self.transformation, distances, iterations = icp_nonrigid(self.source, self.destination)
        # elif method == 'manual'         : mapping_manual(source, destination)
        # elif method == 'automatic'      : mapping_automatic(source, destination)
        else: print('Method not found')

    def show_mapping_transformation(self):
        destination_from_source = self.transform_coordinates(self.source)

        plt.scatter(self.source[:, 0], self.source[:, 1], c='b')
        plt.scatter(self.destination[:, 0], self.destination[:, 1], c='r')
        plt.scatter(destination_from_source[:, 0], destination_from_source[:, 1], c='g')

    def transform_coordinates(self, coordinates):
        if self.transformation_type == 'linear':
            return transform(coordinates, self.transformation)
        else: print('Transformation not found')



# from trace_analysis.icp_nonrigid import icp_nonrigid
# from trace_analysis.image_adapt.polywarp import polywarp_apply
# kx, ky, distances, i = icp_nonrigid(donor,acceptor, tolerance=0.00000001, max_iterations=50)
# acceptor_calculated = polywarp_apply(kx, ky, donor)