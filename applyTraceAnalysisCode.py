# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 10:50:59 2018

@author: ivoseverins
"""


import os

import numpy as np
from matplotlib import pyplot as plt

# Change path to traceAnalysis directory (assuming this file is in that directory)
#os.chdir(os.path.dirname(os.path.abspath(__file__)))
from trace_analysis import Experiment
from trace_analysis.mapping.mapping import Mapping2

# Define path to data, replace by your own directory
# (Now it is set to the twoColourExampleData folder in the traceAnalysis directory)
mainPath = r'./twoColourExampleData'
mainPath = r'O:\Ivo\20190710 - Single-molecule setup (TIR-I)'
mainPath = r'D:\ivoseverins\Desktop\sifx testdata'
mainPath = r'D:\ivoseverins\Desktop\pma testdata\Newly analyzed'
mainPath = r'D:\ivoseverins\SURFdrive\Promotie\Code\Python\traceAnalysis\twoColourExampleData\sifx'
mainPath = r'D:\ivoseverins\Desktop\20190820 drift test microfluidic for multiple imaging rounds of same FOV\1. First test run for drift'

# Initialize an experiment
exp = Experiment(mainPath, colours=('r'))
mov = exp.files[0].movie

image = mov.get_channel(channel='a')
c = mov.find_peaks(image=image, method = 'local-maximum', threshold=8)

mov.write_coordinates_to_pks_file(c)
traces = mov.get_all_traces(c, channel='a')
mov.write_traces_to_traces_file(traces)
#exp.files[0].select()



#Mapping(mappingFilePath)
#mov.make_average_tif(write = True)


#exp.files[-1].use_for_mapping()

#mov.mapping.show_mapping_transformation()


#
# channel = 'donor'
#
# image = mov.get_channel(mov.average_image, channel)
# plt.imshow(image, vmin = 0, vmax = 200)
# coordinates = mov.find_peaks(image, method='local-maximum')
# coordinates = coordinates[mov.is_within_margin(coordinates)]
#
#
# plt.scatter(coordinates[:,0],coordinates[:,1])



# mov.show_coordinates(image, coordinates)
#
#
#
#
#
#
#
#





# from trace_analysis.coordinate_transformations import transform
# from trace_analysis.icp import icp
#
# T, distances, i = icp(donor,acceptor, tolerance=0.0001)
#
# acceptor_calculated = transform(donor, T)
# #



#plt.ioff()
#for i, file in enumerate(exp.files):
#    mov = file.movie
#    #image_mean = mov.make_average_tif()
#    #image_mean_corrected = mov.subtract_background(image_mean)
#
#    im = mov.get_channel(mov.average_image, 'a')
#    coordinates = mov.find_peaks(im , method='local-maximum', threshold = 125)
#
#    coordinates = coordinates[mov.is_within_margin(coordinates, edge = np.array([[0,512],[0,1024]]), margin = 20)]
#
#    mov.show_coordinates(im,coordinates,vmin=450,vmax=850)
#
#    mov.write_coordinates_to_pks_file(coordinates)
#    print(f'{i+1} out of {len(exp.files)}')



#%%
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#X = Y = np.arange(1024)
#X, Y = np.meshgrid(X, Y)
#surf = ax.plot_surface(X,Y,test,cmap=cm.coolwarm,
#                       linewidth=0, antialiased=False)

#%%
#plt.figure()
#plt.imshow()
#plt.imshow(mov.get_channel(image_mean,'acceptor'), vmax=200)
#plt.scatter(pts[:,0],pts[:,1], marker = 'o', facecolors='none', edgecolors='r')
#plt.show()
