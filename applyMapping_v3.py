# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 16:40:50 2019

@author: mwdocter
"""

# applyMapping: file in other folder which needs to improt mapping from pick_spots_akaze

import sys
import time
import analyze_label
import cv2

from autopick.do_before import clear_all
from autopick.pick_spots_akaze_v3 import mapping #akaze_final
from load_file import read_one_page_pma
from image_adapt.rolling_ball import rollingball
from image_adapt.find_threshold import remove_background
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
file_tetra="E://CMJ trace analysis/autopick/tetraspeck.tif" # most likely not related to the data
transform_matrix=mapping(file_tetra,show=1,bg=0)[0]
#
file_tetra='H://projects//research practicum//single molecule fluorescence//Matlab//HJA-data from Ivo//hel6_ave.tif'
transform_matrix=mapping(file_tetra,show=1,bg=None, tol=0)[0]
# alternatives: 
if 0:
    transform_matrix,_=mapping(file_tetra,0)
    transform_matrix,pts=mapping(file_tetra,0)
    transform_matrix,pts=mapping(file_tetra,1)

 
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
labels_ALL=analyze_label.analyze(im_correct2)[1]
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
acceptor=np.zeros((pts_number,nImages ))
acceptorB=np.zeros(( pts_number,nImages))

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
    # don't transform the full image, transform coordinates
    imA=cv2.warpPerspective(IM_acceptor.astype(float), transform_matrix,array_size[::-1])
    
    dstG = cv2.perspectiveTransform(ptsG.reshape(-1,1,2),transform_matrix) #approach2: transformcoordinates

    for jj in range(0,pts_number):
        xpix=ptsG[jj][1]
        ypix=ptsG[jj][0]
        
        if xpix>10 and xpix<hdim-10 and ypix>10 and ypix<vdim/2-10 and size_label[jj]<100:
            xpix_int=int(xpix)
            ypix_int=int(ypix)
            xf=dstG[jj,0,1]#approach2: transformcoordinates
            yf=dstG[jj][0][0]+256#approach2: transformcoordinates
            xf_int=int(xf)#approach2: transformcoordinates
            yf_int=int(yf)#approach2: transformcoordinates
            #first crop around spot, then do multiplication
            impixD=IM_donor[ (xpix_int-5) : (xpix_int+6) , (ypix_int-5) : (ypix_int+6)]
            labD=labels[ (xpix_int-5) : (xpix_int+6) , (ypix_int-5) : (ypix_int+6)]
            #impixD[labD!=jj]=0
          
            GG=makeGaussian(11, fwhm=3, center=(ypix-ypix_int+5,xpix-xpix_int+5))
            multipD=impixD*GG
            donor[jj,ii]=np.sum(multipD)
            if show:
                    plt.subplot(1,3,1)
                    plt.imshow(impixD)
                    plt.subplot(1,3,2)
                    plt.imshow(GG)
                    plt.subplot(1,3,3)
                    plt.imshow(multipD)
                    plt.title(jj)
           
             ###  
            impixA=imA[ (xpix_int-5) : (xpix_int+6) , (ypix_int-5) : (ypix_int+6)]
            labA=labels[ (xpix_int-5) : (xpix_int+6) , (ypix_int-5) : (ypix_int+6)]
            #impix[lab!=jj]=0
            multip=impixA*GG
            acceptor[jj,ii]=np.sum(multip)
            
            impixB=im_correct2[ (xf_int-5) : (xf_int+6) , (yf_int-5) : (yf_int+6)]#approach2: transformcoordinates
            labB=labels_ALL[ (xf_int-5) : (xf_int+6) , (yf_int-5) : (yf_int+6)]
            GGB=makeGaussian(11, fwhm=3, center=(yf-yf_int+5,xf-xf_int+5))#approach2: transformcoordinates
            multipB=impixB*GGB#approach2: transformcoordinates
            acceptorB[jj,ii]=np.sum(multipB)

elapsed = time.time() - t   
print('after extracting intensities of all positions all image {0:f}'.format(elapsed))
  
# save all traces file
# in comparison to IDL: film_l=nImages, num_good=pts_number
# in comparison to Trace.m: 
Nframes=nImages
Ntraces=pts_number
Ncolours=2
with open(root+'\\'+name[:-4]+'-P.traces', 'w') as outfile:
    off = np.array([Nframes], dtype=np.int32)
    off.tofile(outfile)
    off = np.array([Ntraces], dtype=np.int16)
    off.tofile(outfile)
    time_tr=np.zeros((Nframes,Ntraces))
    for jj in range(Ntraces//Ncolours):
        time_tr[:,jj*2] = donor[jj,:]
        time_tr[:,jj*2+1]=  acceptor[jj,:]
    off = np.array((time_tr), dtype=np.int16)
    off.tofile(outfile)
elapsed = time.time() - t   
print('after saving traces to file {0:f}'.format(elapsed))

#saving to pks file
with open(root+'\\'+name[:-4]+'a.pks', 'w') as outfile:
     for jj in range(0,pts_number):
         xpix=ptsG[jj][1]
         ypix=ptsG[jj][0]
         outfile.write('{0:4.0f}\t {1:4.2f}\t {2:4.2f}\t {3:4.2f}\t {4:4.2f}\n'.format(jj, xpix, ypix, 0, 0, width4=4, width6=6))
