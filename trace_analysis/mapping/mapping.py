# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 13:55:53 2019

@author: ivoseverins
"""

import numpy as np
import matplotlib.pyplot as plt

from trace_analysis.mapping.icp import icp, icp_apply_transform
#from trace_analysis.coordinate_transformations import transform
#from trace_analysis.image_adapt.polywarp import polywarp, polywarp_apply #required for nonlinear

class Mapping2:
    def __init__(self, source = None, destination = None, method = None,
                 transformation_type = 'linear', dest2source_translation = None):
        self.source = source #source=donor=left side image
        self.destination = destination #destination=acceptor=right side image
        self.method = method
        self.transformation_type = transformation_type
        self.transformation = None
        self.dest2source_translation = dest2source_translation
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
            self.transformation, distances, iterations, self.transformation_inverse, self.dest2source_translation = \
                icp(self.source, self.destination, dest2source_translation=self.dest2source_translation,
                    transformation_type=self.transformation_type)
        # elif method == 'manual'         : mapping_manual(source, destination)
        # elif method == 'automatic'      : mapping_automatic(source, destination)
        else: print('Method not found')

    def show_mapping_transformation(self, figure=None):
        if not figure: figure = plt.figure()
        destination_from_source = self.transform_coordinates(self.source)

        axis = figure.gca()

        axis.scatter(self.source[:, 0], self.source[:, 1], c='g')
        axis.scatter(self.destination[:, 0], self.destination[:, 1], c='r')
        axis.scatter(destination_from_source[:, 0], destination_from_source[:, 1], c='y')

    def transform_coordinates(self, coordinates, inverse = False):
        print(self.transformation)
<<<<<<< HEAD
        if self.transformation_type == 'linear':
            if inverse is False: return transform(coordinates, self.transformation)
            elif inverse is True: return transform(coordinates, self.transformation_inverse)

        elif self.transformation_type == 'nonlinear':
            if inverse is False: return polywarp_apply(self.transformation[0],self.transformation[1],coordinates)
            elif inverse is True: return polywarp_apply(self.transformation_inverse[0],self.transformation_inverse[1],coordinates)  
        else: print('Transformation type not found')
        #still to make nonlinear?? or use polywarp_apply    
        
        #    if inverse is False: return 
#=======
#        elif self.transformation_type=='nonlinear' : #still to be tested
#            from trace_analysis.image_adapt.polywarp import polywarp, polywarp_apply
#            if inverse is False: return polywarp_apply(self.transformation[0], self.transformation[1], coordinates)
#            elif inverse is True: 
#                 destinate=polywarp_apply(self.transformation[0], self.transformation[1], coordinates)
#                 kx, ky = polywarp(coordinates[0],coordinates[1],\
#                              destinate[0],destinate[1]) 
#                 return polywarp_apply(kx,ky, coordinates)
#                 
#        #still to make nonlinear?? or use polywarp_apply    
#        else: print('Transformation not found')
##            
#>>>>>>> da2376d467064b10880c3b6b836b42c93f56e740


# from trace_analysis.icp_nonrigid import icp_nonrigid
# from trace_analysis.image_adapt.polywarp import polywarp_apply
# kx, ky, distances, i = icp_nonrigid(donor,acceptor, tolerance=0.00000001, max_iterations=50)
# acceptor_calculated = polywarp_apply(kx, ky, donor)
=======
        if self.method == 'icp':
            return icp_apply_transform(coordinates, inverse, self.transformation,self.transformation_inverse, self.transformation_type, self.dest2source_translation)
                                     
        else: print('transform_coordinates only works for icp')
>>>>>>> 9cbf291dd367d0509b050101527471f809805c78
