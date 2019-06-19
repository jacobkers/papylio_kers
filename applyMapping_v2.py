# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 16:40:50 2019

@author: mwdocter
order = transform - find spot locations - extract traces - save all

"""

# applyMapping: file in other folder which needs to improt mapping from pick_spots_akaze

import sys
import time
import analyze_label
import cv2

from autopick.do_before import clear_all
from autopick.pick_spots_akaze_final import mapping
from load_file import read_one_page_pma
from rolling_ball import rollingball
from find_threshold import remove_background
from find_xy_position.Gaussian import makeGaussian

import numpy as np
import matplotlib.pyplot as plt


clear_all()
show=0
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
# same as in K:\bn\cmj\Shared\Ivo\Voorbeelddata Holiday junction
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

donor=np.zeros((pts_number,nImages ))
acceptor=np.zeros(( pts_number,nImages))

size_label=np.zeros(shape=(pts_number))
for jj in range(pts_number):
    size_label[jj]=np.sum(labels==jj)

t = time.time()  
for ii in range(0,nImages):
    im=(read_one_page_pma(root,name, pageNb=ii))[0] # this one opens, closes all the time
    im_correct=im-im_bg
    im_correct[im_correct<0]=0
    im_correct2=remove_background(im_correct, threshold)[0]
    print(ii)
    IM_donor=im_correct2[:,0:int(vdim/2)]
    IM_acceptor=im_correct2[:,int(vdim/2):]
    array_size=np.shape(IM_acceptor)
    imA=cv2.warpPerspective(IM_acceptor.astype(float), transform_matrix,array_size[::-1])
    
#    elapsed = time.time() - t ;   print('time before for loop {0:f}'.format(elapsed))
    for jj in range(0,pts_number):
        
        xpix=ptsG[jj][1]
        ypix=ptsG[jj][0]
        if xpix>10 and xpix<hdim-10 and ypix>10 and ypix<vdim/2-10 and size_label[jj]<100:
#            t = time.time()
            xpix_int=int(xpix)
            ypix_int=int(ypix)
            #first crop around spot, then do multiplication
            impix=IM_donor[ (xpix_int-5) : (xpix_int+6) , (ypix_int-5) : (ypix_int+6)]
#            elapsed = time.time() - t; print('{0:f}'.format(elapsed))  
            lab=labels[ (xpix_int-5) : (xpix_int+6) , (ypix_int-5) : (ypix_int+6)]
#            elapsed = time.time() - t; print('{0:f}'.format(elapsed))  
            impix[lab!=jj]=0
#            elapsed = time.time() - t; print('{0:f}'.format(elapsed))   
           
            GG=makeGaussian(11, fwhm=3, center=(ypix-ypix_int+5,xpix-xpix_int+5))
#            elapsed = time.time() - t; print('{0:f}'.format(elapsed))  
            multip=impix*GG
#            elapsed = time.time() - t; print('{0:f}'.format(elapsed))   
           
            if show:
                    plt.subplot(1,3,1)
                    plt.imshow(impix)
                    plt.subplot(1,3,2)
                    plt.imshow(GG)
                    plt.subplot(1,3,3)
                    plt.imshow(multip)
                    plt.title(jj)
            donor[jj,ii]=np.sum(multip)

             ###  
#            elapsed = time.time() - t; print('{0:f}'.format(elapsed))   
            impix=imA[ (xpix_int-5) : (xpix_int+6) , (ypix_int-5) : (ypix_int+6)]
#            elapsed = time.time() - t; print('{0:f}'.format(elapsed))   
            impix[lab!=jj]=0
#            elapsed = time.time() - t; print('{0:f}'.format(elapsed))   
            multip=impix*GG
#            elapsed = time.time() - t; print('{0:f}'.format(elapsed))   
            acceptor[jj,ii]=np.sum(multip)
#            elapsed = time.time() - t; print('{0:f}'.format(elapsed))   
#            plt.waitforbuttonpress()
#            plt.pause(0.3)

            # sum all intensities and store in acceptro donor array
    # do stuff
elapsed = time.time() - t   
print('after extracting intensities of all positions all image {0:f}'.format(elapsed))
  
    # save all traces in pks and trace file
film_l=nImages
num_good=pts_number
    

with open(root+'\\'+name[:-4]+'a.traces', 'w') as outfile:
    off = np.array([film_l], dtype=np.int32)
    off.tofile(outfile)
    off = np.array([2*num_good], dtype=np.int16)
    off.tofile(outfile)
    for jj in range(pts_number):
        off = np.array(donor[jj,:], dtype=np.int16)
        off.tofile(outfile)
        off = np.array(acceptor[jj,:], dtype=np.int16)
        off.tofile(outfile)
elapsed = time.time() - t   
print('after saving traces to file {0:f}'.format(elapsed))
