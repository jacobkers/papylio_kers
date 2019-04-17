# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 16:40:50 2019

@author: mwdocter
"""

# applyMapping: file in other folder which needs to improt mapping from pick_spots_akaze

import sys
from autopick.do_before import clear_all

from autopick.pick_spots_akaze_final import mapping
from load_file import read_one_page_pma
from rolling_ball import rollingball
from find_xy_position.Gaussian import makeGaussian
import numpy as np
import matplotlib.pyplot as plt
from find_threshold import remove_background
import time
import analyze_label
import cv2

clear_all()

#
sys.path.insert(0,"E://CMJ trace analysis") 
#sys.path.append("E://CMJ trace analysis/autopick") 
sys.path.insert(0,"E://CMJ trace analysis/autopick") 

#
#!@# first mapping
#
tetra_image_location="E://CMJ trace analysis/autopick/tetraspeck.tif"

transform_matrix=mapping(tetra_image_location,0)[0]
# alternatives: 
if 0:
    transform_matrix,_=mapping("E://CMJ trace analysis/autopick/tetraspeck.tif",0)
    transform_matrix,pts=mapping("E://CMJ trace analysis/autopick/tetraspeck.tif",0)
    transform_matrix,pts=mapping("E://CMJ trace analysis/autopick/tetraspeck.tif",1)
    
#later on to be used on acceptor files: 
# im4=cv2.warpPerspective(gray2, transformation_matrix, array_size[::-1] )

#
#!@# from data: find average over first 20 frames
#
root='H:\\projects\\research practicum\\single molecule fluorescence\\Matlab\\HJA-data from Ivo'
name='hel6.pma'

#follow IDL path: sum 20 image, find spots& background. Then loop over all image, do background subtract+ extract traces 
#note: a second path would be localization sub pixel, which includes local background estimation/subtraction
_,hdim,vdim,nImages=(read_one_page_pma(root,name, pageNb=0))
im_sum=(read_one_page_pma(root,name, pageNb=0 )[0]).astype(float)
for ii in range (1,20):
    im=(read_one_page_pma(root,name, pageNb=ii)[0]).astype(float)
    im[im<0]=0
    im_sum=im_sum+im
im_mean20=(im_sum/20).astype(int)

#
#!@# from data: find background, to later subtract
#
im_bg=rollingball(im_mean20)[1]

im_correct=im_mean20-im_bg
im_correct[im_correct<0]=0

im_correct2,threshold=remove_background(im_correct)
#
#!@# find positions of corresponding points
#
###_,ptsG,ptsR=mapping(im_correct2,0)
###pts_number= len(ptsG)
pts_number,labels,ptsG=analyze_label.analyze(im_correct2[:,0:int(vdim/2)])

#
#!@# extract the data of all points in a loop
#
#
circle=np.zeros(shape=(11,11))
circle[0] = [ 0,0,0,0,0,0,0,0,0,0,0]# a typical way of defining a vector
circle[1] = [ 0,0,0,0,1,1,1,0,0,0,0]
circle[2] = [ 0,0,0,1,0,0,0,1,0,0,0]
circle[3] = [ 0,0,1,0,0,0,0,0,1,0,0]
circle[4] = [ 0,1,0,0,0,0,0,0,0,1,0]
circle[5] = [ 0,1,0,0,0,0,0,0,0,1,0]
circle[6] = [ 0,1,0,0,0,0,0,0,0,1,0]
circle[7] = [ 0,0,1,0,0,0,0,0,1,0,0]
circle[8] = [ 0,0,0,1,0,0,0,1,0,0,0]
circle[9] = [ 0,0,0,0,1,1,1,0,0,0,0]
circle[10]= [ 0,0,0,0,0,0,0,0,0,0,0]

# double check these spots are good (well separated)

good_spots=np.zeros(pts_number)

for jj in range(pts_number):
    labels_bin=np.zeros(np.shape(labels))
    labels_bin[labels==jj]=1
    IM_in=im_correct2[:,0:int(vdim/2)]*labels_bin.astype(float)
    xpix=ptsG[jj][1]
    xpix_int=int(xpix)
    ypix=ptsG[jj][0]
    ypix_int=int(ypix)
    if xpix>10 and xpix<hdim-10 and ypix>10 and ypix<vdim/2-10:
        # extract 7x7 pixels around xy location
        
            impix=IM_in[ (xpix_int-5) : (xpix_int+6) , (ypix_int-5) : (ypix_int+6)]
        # check whether intensity is really above background
       
#        if np.amax(circle*impix)> threshold+0.65*im[xpix_int,ypix_int]:
#        if np.amax(impix)==im[xpix_int,ypix_int]:
            good_spots[jj]=1
            # make gaussian and multiply
            GG=makeGaussian(11, fwhm=3, center=(ypix-ypix_int+5,xpix-xpix_int+5))
            multip=impix*GG
            plt.subplot(2,2,1)
            plt.imshow(impix)
            plt.subplot(2,2,2)
            plt.imshow(GG)
            plt.subplot(2,2,3)
            plt.imshow(multip)
            plt.title(jj)
            plt.waitforbuttonpress()
#            plt.pause(0.3)
#            break  
        # sum all intensities and store in acceptro donor array
        
# save all traces in pks and trace file