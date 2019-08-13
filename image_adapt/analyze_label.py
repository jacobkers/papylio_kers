# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 13:50:42 2019

@author: https://stackoverflow.com/questions/35854197/how-to-use-opencvs-connected-components-with-stats-in-python

returns the number of spots, and per spot its centroid and the number of pixels (can be used to discard too large spots)
"""
import cv2
import numpy as np
from ..autopick.pick_spots_akaze_ALL import enhance_blobies_single
from .find_threshold import get_threshold #,remove_background
 
def analyze(src, threshold = None):
#    if 0:
#        # old version
#        ret, thresh = cv2.threshold(src.astype(np.uint8),0,1,cv2.THRESH_BINARY)
#    
#        # You need to choose 4 or 8 for connectivity type
#        connectivity = 4  
#        # Perform the operation
#        output = cv2.connectedComponentsWithStats(thresh, connectivity, cv2.CV_32S)
#        # Get the results
#        # The first cell is the number of labels
#        num_labels = output[0]
#        # The second cell is the label matrix
#        labels = output[1]
#        
#        size_label=np.zeros(num_labels)
#        for jj in range(num_labels):
#            size_label[jj]=np.sum(labels==jj)
#        
#        # The third cell is the stat matrix
#        #stats = output[2]
#        # The fourth cell is the centroid matrix
#        ctrd = output[3]
        
    #    for ii in num_labels:
    #        LL=labels.copy()
    #        LL[LL!=ii]=0
    #        xx,yy=np.where(LL>0)
    #        centroids[ii,0]=np.mean(xx)
    #        centroids[ii,1]=np.mean(yy)
        #apparentely centroids is already np.mean
        
#    elif 0:
#        #also old version, fL does not depend on image at all
#        ctrd=[]
#        if len(src)>512: LL=50*16
#        else: LL=200 
#        fL=50000
#        tol=1
#        while np.shape(ctrd)[0]<LL: 
#            # while loop to lower f and increase the number of spots found
#            gray1= enhance_blobies_single(src[:, :src.shape[1]//2],fL, tol)
#            # initialize the AKAZE descriptor, then detect keypoints and extract
#            # local invariant descriptors from the image
#            detector = cv2.AKAZE_create()
#            (kps1, descs1) = detector.detectAndCompute(gray1, None)
#            ctrd=cv2.KeyPoint_convert(kps1);
#            fL=fL*0.9
#    else:
#   
    if not threshold:
        fL=get_threshold(src)
    else:
        fL = threshold
    gray1=enhance_blobies_single(src,fL,1) #remove_background(src, fL)
    detector = cv2.AKAZE_create()
    (kps1, descs1) = detector.detectAndCompute(gray1, None)
    ctrd=np.array(cv2.KeyPoint_convert(kps1))
                
    # remove all pixels at the edge (within 10 pix) 
    num_labels=   len(ctrd) 
    dim1,dim0=np.shape(src)
    for ii in range(num_labels-1,-1,-1):
        discard=ctrd[ii,0]<10 or ctrd[ii,1]<10 or ctrd[ii,0]>dim0-10 or ctrd[ii,1]>dim1-10 or src[int(ctrd[ii,1]),int(ctrd[ii,0])]==0# or size_label[ii]>100
        #for some reason also spots are found on immean20 with no intensity --> discard
        if discard:
            ctrd=np.delete(ctrd,ii, axis=0)
   #         size_label=np.delete(size_label,ii, axis=0)
    num_labels=   len(ctrd) 
    
      
    return num_labels,0, ctrd

 