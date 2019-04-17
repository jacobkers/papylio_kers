# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 13:50:42 2019

@author: https://stackoverflow.com/questions/35854197/how-to-use-opencvs-connected-components-with-stats-in-python
"""
import cv2
import numpy as np

def analyze(src):
    ret, thresh = cv2.threshold(src.astype(np.uint8),0,1,cv2.THRESH_BINARY)
    # You need to choose 4 or 8 for connectivity type
    connectivity = 4  
    # Perform the operation
    output = cv2.connectedComponentsWithStats(thresh, connectivity, cv2.CV_32S)
    # Get the results
    # The first cell is the number of labels
    num_labels = output[0]
    # The second cell is the label matrix
    labels = output[1]
    # The third cell is the stat matrix
    #stats = output[2]
    # The fourth cell is the centroid matrix
    centroids = output[3]
    
#    for ii in num_labels:
#        LL=labels.copy()
#        LL[LL!=ii]=0
#        xx,yy=np.where(LL>0)
#        centroids[ii,0]=np.mean(xx)
#        centroids[ii,1]=np.mean(yy)
    #apparentely centroids is already np.mean
        
    return num_labels,labels,centroids