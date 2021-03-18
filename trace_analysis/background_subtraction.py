# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 16:10:42 2021

input = file, method, *args
methods+args:
    Channel_mean
    Channel_median
    ROI_min, args= filter_neighbourhood_size ** DEFAULT
    ROI_median, args= filter_neighbourhood_size
    Surrounding_mean, filter=radius
    Surrounding_median, filter=radius
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage.filters as filters

def extract_background(file, method='ROI_min',  *args):
    coordinates=file.coordinates
    coordinates[:,0]+=150
    average_image=file.average_image
    #TODO: later on extend it to a loop over all individual images, for now stick to average_image
    donor_image = file.movie.get_channel(image=average_image, channel='d')
    acceptor_image = file.movie.get_channel(image=average_image, channel='a')
    background=np.zeros(len(coordinates))
    
    default_filtersize=10
    default_radius=5
    if method=='Channel_mean':
        background[0::2]= np.mean(donor_image)
        background[1::2]= np.mean(acceptor_image)
        # TODO: extend to multiple channels
    elif method=='Channel_median': 
        background[0::2]= np.median(donor_image)
        background[1::2]= np.median(acceptor_image)
    elif method=='ROI_min':
        try: 
            filter_neighbourhood_size=args[0]
        except: 
            filter_neighbourhood_size=default_filtersize
        image_min = filters.minimum_filter(average_image, filter_neighbourhood_size) 
        #TODO: improve edge between donor&acceptor channel
        plt.imshow(image_min)
        for ii in range(file.number_of_channels*file.number_of_molecules):
            xpos=int(coordinates[ii][0])
            ypos=int(coordinates[ii][1])
            background[ii]=image_min[ypos,xpos]
            
    elif method=='ROI_median':
        try: filter_neighbourhood_size=args[0]
        except: filter_neighbourhood_size=default_filtersize
        image_median = filters.median_filter(average_image, filter_neighbourhood_size) 
        #TODO: improve edge between donor&acceptor channel
        for ii in range(file.number_of_channels*file.number_of_molecules):
            xpos=int(coordinates[ii][0])
            ypos=int(coordinates[ii][1])
            background[ii]=image_median[ypos,xpos]
            
    elif method=='Surrounding_mean': #not optimal when multiple spots are close, better to use median
        from trace_analysis.coordinate_optimization import crop, circle
        try: radius=args[0]
        except: radius=default_radius
        circle_matrix = circle(radius)
        for ii in range(file.number_of_channels*file.number_of_molecules):
            cropped_peak = crop(average_image, coordinates[ii], radius*2+1)
            AA=cropped_peak * circle_matrix
            background[ii]=np.mean(AA[np.nonzero(AA)])
            
    elif method=='Surrounding_median': #similar to coordinate_optimization, coordinates_without_intensity_at_radius
        from trace_analysis.coordinate_optimization import crop, circle
        radius=args[0]
        circle_matrix = circle(radius)
        for ii in range(file.number_of_channels*file.number_of_molecules):
            cropped_peak = crop(average_image, coordinates[ii], radius*2+1)
            AA=cropped_peak * circle_matrix
            background[ii]=np.median(AA[np.nonzero(AA)])
            
    return background