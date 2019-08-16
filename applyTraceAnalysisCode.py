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

# Define path to data, replace by your own directory
# (Now it is set to the twoColourExampleData folder in the traceAnalysis directory)
mainPath = r'./twoColourExampleData'
mainPath = r'O:\Ivo\20190710 - Single-molecule setup (TIR-I)'
mainPath = r'D:\ivoseverins\Desktop\sifx testdata'

# Initialize an experiment
exp = Experiment(mainPath)

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