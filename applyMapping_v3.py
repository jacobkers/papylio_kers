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
import autopick.pick_spots_akaze_manual 
import autopick.pick_spots_akaze_final 
from load_file import read_one_page_pma
from image_adapt.rolling_ball import rollingball
from image_adapt.find_threshold import remove_background
from find_xy_position.Gaussian import makeGaussian

import numpy as np
import matplotlib.pyplot as plt


clear_all()
plt.close("all")
show=0
#
sys.path.insert(0,"E://CMJ trace analysis") 
#sys.path.append("E://CMJ trace analysis/autopick") 
sys.path.insert(0,"E://CMJ trace analysis/autopick") 

#
#!@# first mapping
# Automatic mapping: use this if you have an accompanying tetra speck image RECOMMENDED
#file_tetra="E://CMJ trace analysis/autopick/tetraspeck.tif" # most likely not related to the data
file_tetra='N:\\tnw\\BN\\CMJ\\Shared\\Margreet\\20181203_HJ_training\\mapping\\rough_ave.tif'
transform_matrix=autopick.pick_spots_akaze_final.mapping(file_tetra,show=1,bg=None, tol=0,f=10000)[0]
plt.show()
plt.pause(0.05) # in hope the figures display with mapping images

print('Are you satisfied with the mapping (yes/no)?')
x = input()
if x[0]!='y': # do manual mapping
    transform_matrix=autopick.pick_spots_akaze_manual.mapping(file_tetra,show=1,bg=None, tol=0,f=10000)[0]
plt.show()
plt.pause(0.05) 

#example data (without mapping) in K:\bn\cmj\Shared\Ivo\Voorbeelddata Holiday junction

root='N:\\tnw\\BN\\CMJ\\Shared\\Margreet\\20181203_HJ_training\\#3.10_0.1mg-ml_streptavidin_50pM_HJC_G_movies'
name='hel4.pma'

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
if 0: #check whether rolling ball went well
    plt.figure()
    plt.subplot(1,3,1), plt.imshow(im_mean20)
    plt.subplot(1,3,2), plt.imshow(im_bg)
    plt.subplot(1,3,3), plt.imshow(im_correct2)
#!@# find positions of corresponding points
#

pts_number,label_size,ptsG=analyze_label.analyze(im_correct2[:,0:int(vdim/2)])
    
#labels_ALL=analyze_label.analyze(im_correct2)[1]
#
#!@# extract the data of all points in a loop
#
##
#circle=np.zeros(shape=(11,11))
#circle[0] = [ 0,0,0,0,0,0,0,0,0,0,0]# a typical way of defining a vector
#circle[1] = [ 0,0,0,0,1,1,1,0,0,0,0]
#circle[2] = [ 0,0,0,1,0,0,0,1,0,0,0]
#circle[3] = [ 0,0,1,0,0,0,0,0,1,0,0]
#circle[4] = [ 0,1,0,0,0,0,0,0,0,1,0]
#circle[5] = [ 0,1,0,0,0,0,0,0,0,1,0]
#circle[6] = [ 0,1,0,0,0,0,0,0,0,1,0]
#circle[7] = [ 0,0,1,0,0,0,0,0,1,0,0]
#circle[8] = [ 0,0,0,1,0,0,0,1,0,0,0]
#circle[9] = [ 0,0,0,0,1,1,1,0,0,0,0]
#circle[10]= [ 0,0,0,0,0,0,0,0,0,0,0]

# double check these spots are good (well separated)

donor=np.zeros((pts_number,nImages ))
acceptor=np.zeros((pts_number,nImages ))
acceptorB=np.zeros(( pts_number,nImages))
acceptorC=np.zeros(( pts_number,nImages))

#transform location acceptor xy before loop, plus discard if they go out of boundary
dstG= cv2.perspectiveTransform(ptsG.reshape(-1,1,2),np.linalg.inv(transform_matrix)) #reshape needed for persp.transform
dstG= dstG.reshape(-1,2)

for ii in range(pts_number-1,-1,-1):
        discard=dstG[ii,0]<10 or dstG[ii,1]<10 or dstG[ii,0]>hdim-10 or dstG[ii,1]>vdim-10
        if discard:
            ptsG=np.delete(ptsG,ii, axis=0)
            dstG=np.delete(dstG,ii, axis=0)
            label_size=np.delete(label_size,ii, axis=0)
pts_number= len(dstG) 
    

dstG2=np.array([[ii[0]+256,ii[1]] for ii in dstG])

# test to see you extract pixels correctly
if 0:
    IM_donor=im_correct2[:,0:int(vdim/2)]
    IM_acceptor=im_correct2[:,int(vdim/2):]
    array_size=np.shape(IM_acceptor)
    imA=cv2.warpPerspective(IM_acceptor.astype(float), transform_matrix,array_size[::-1])
    DD=10
    for jj in range(0,pts_number):
        xpix=ptsG[jj][1]
        ypix=ptsG[jj][0]
        xpix_int=int(xpix)
        ypix_int=int(ypix)
        impixD=IM_donor[ (xpix_int-DD) : (xpix_int+DD+1) , (ypix_int-DD) : (ypix_int+DD+1)]
        impixE=IM_acceptor[ (xpix_int-DD) : (xpix_int+DD+1) , (ypix_int-DD) : (ypix_int+DD+1)]
        xf2=dstG2[jj][1]#approach3,  NOTE: xy are swapped here
        yf2=dstG2[jj][0]#approach3
        xf2_int=int(xf2)#approach3
        yf2_int=int(yf2)#approach3
                        
        impixC=im_correct2[ (xf2_int-DD) : (xf2_int+DD+1) , (yf2_int-DD) : (yf2_int+DD+1)]
        
        plt.figure(20), 
        plt.subplot(1,3,1), plt.imshow(impixD), plt.title(jj)
        plt.subplot(1,3,2), plt.imshow(impixE), plt.title([xpix_int,ypix_int])
        plt.subplot(1,3,3), plt.imshow(impixC), plt.title([xf2_int,yf2_int])
        plt.pause(1)
        
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
    # approach 1: take location ptsG from IM_donor and imA (the transformed IM_acceptor)
    # approach2: take location ptsG from IM_donor and the transformed coordinates ptsG from IM_acceptor
    # approach 3, most similar to IDL: take ptsG and dstG from im_correct2
    # for later concern: compare the three outcomes of the three approaches and compare to find the best
    imA=cv2.warpPerspective(IM_acceptor.astype(float), transform_matrix,array_size[::-1])

       
    for jj in range(0,pts_number):
        xpix=ptsG[jj][1]
        ypix=ptsG[jj][0]
        
        if 1: #already incorporated in reshaping ptsG xpix>10 and xpix<hdim/2-10 and ypix>10 and ypix<vdim-10 and label_size[jj]<100 and dstG[jj,1]>10 and dstG[jj,1]<hdim-10 and dstG[jj,0]>10 and dstG[jj,0]<vdim/2-10: #also include dstG here!!!!
            xpix_int=int(xpix)
            ypix_int=int(ypix)
            
            #first crop around spot, then do multiplication
            impixD=IM_donor[ (xpix_int-5) : (xpix_int+6) , (ypix_int-5) : (ypix_int+6)]
            GG=makeGaussian(11, fwhm=3, center=(ypix-ypix_int+5,xpix-xpix_int+5))
            multipD=impixD*GG
            donor[jj,ii]=np.sum(multipD)
#            if jj==82: #testing purposes
#                print(ii,np.sum(multipD),donor[jj,ii])
            if show:
                    plt.subplot(1,3,1)
                    plt.imshow(impixD)
                    plt.subplot(1,3,2)
                    plt.imshow(GG)
                    plt.subplot(1,3,3)
                    plt.imshow(multipD)
                    plt.title(jj)
            ###  
            # approach 1, find same coordinates in transformed image
            impixA=imA[ (xpix_int-5) : (xpix_int+6) , (ypix_int-5) : (ypix_int+6)]
            multip=impixA*GG
            acceptor[jj,ii]=np.sum(multip)
            
            if 1: #testing approaches
                # approach 2, find transformed coordinates in individual images
                xf=dstG[jj][1]#approach2: transformcoordinates, NOTE: xy are swapped here
                yf=dstG[jj][0]#approach2: transformcoordinates
                xf_int=int(xf)#approach2: transformcoordinates
                yf_int=int(yf)#approach2: transformcoordinates
                
                impixB=IM_acceptor[ (xf_int-5) : (xf_int+6) , (yf_int-5) : (yf_int+6)]#approach2: transformcoordinates
                GGB=makeGaussian(11, fwhm=3, center=(yf-yf_int+5,xf-xf_int+5))#approach2: transformcoordinates
                multipB=impixB*GGB#approach2: transformcoordinates
                acceptorB[jj,ii]=np.sum(multipB)
                
                # approach 3, find transformed coordinates in im_correct (full image)
                xf2=dstG2[jj][1]#approach3
                yf2=dstG2[jj][0]#approach3
                xf2_int=int(xf2)#approach3
                yf2_int=int(yf2)#approach3
                
                impixC=im_correct2[ (xf2_int-5) : (xf2_int+6) , (yf2_int-5) : (yf2_int+6)]
                GGC=makeGaussian(11, fwhm=3, center=(yf2-yf2_int+5,xf2-xf2_int+5))#approach3
                multipC=impixC*GGC#approach3
                acceptorC[jj,ii]=np.sum(multipC)
       # else: print([jj,ptsG[jj][1],ptsG[jj][0],dstG[jj,0],dstG[jj,1],label_size[jj]])

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
    off = np.array([2*Ntraces], dtype=np.int16)
    off.tofile(outfile)
    time_tr=np.zeros((Nframes,2*Ntraces))
    for jj in range(2*Ntraces//Ncolours):
            time_tr[:,jj*2] = donor[jj,:]
            time_tr[:,jj*2+1]=  acceptor[jj,:]
    off = np.array((time_tr), dtype=np.int16)
    off.tofile(outfile)
elapsed = time.time() - t   
print('after saving traces to file {0:f}'.format(elapsed))

if 1: #testing approaches
     
    with open(root+'\\'+name[:-4]+'-PB.traces', 'w') as outfile:
        off = np.array([Nframes], dtype=np.int32)
        off.tofile(outfile)
        off = np.array([2*Ntraces], dtype=np.int16)
        off.tofile(outfile)
        time_tr=np.zeros((Nframes,2*Ntraces))
        for jj in range(2*Ntraces//Ncolours):
            time_tr[:,jj*2] = donor[jj,:]
            time_tr[:,jj*2+1]=  acceptorB[jj,:]
        off = np.array((time_tr), dtype=np.int16)
        off.tofile(outfile)
    elapsed = time.time() - t   
    print('after saving traces to file {0:f}'.format(elapsed))
    
    with open(root+'\\'+name[:-4]+'-PC.traces', 'w') as outfile:
        off = np.array([Nframes], dtype=np.int32)
        off.tofile(outfile)
        off = np.array([2*Ntraces], dtype=np.int16)
        off.tofile(outfile)
        time_tr=np.zeros((Nframes,2*Ntraces))
        for jj in range(2*Ntraces//Ncolours):
            time_tr[:,jj*2] = donor[jj,:]
            time_tr[:,jj*2+1]=  acceptorC[jj,:]
        off = np.array((time_tr), dtype=np.int16)
        off.tofile(outfile)
    elapsed = time.time() - t   
    print('after saving traces to file {0:f}'.format(elapsed))

if 1: #testing approaches
    #saving to pks file
    with open(root+'\\'+name[:-4]+'-PC.pks', 'w') as outfile:
         for jj in range(0,pts_number):
             pix0=ptsG[jj][0]
             pix1=ptsG[jj][1]
             outfile.write(' {0:4.0f} {1:4.4f} {2:4.4f} {3:4.4f} {4:4.4f} \n'.format((jj*2)+1, pix0, pix1, 0, 0, width4=4, width6=6))
             pix0=dstG2[jj][0]
             pix1=dstG2[jj][1]
             outfile.write(' {0:4.0f} {1:4.4f} {2:4.4f} {3:4.4f} {4:4.4f} \n'.format((jj*2)+2, pix0, pix1, 0, 0, width4=4, width6=6))

#saving to .map file
P=np.zeros((4,4))   
tm=np.linalg.inv(transform_matrix)      
P[0,0]=tm[0,2]
P[0,1]=tm[0,1]
P[1,0]=tm[0,0]
Q=np.zeros((4,4))         
Q[0,0]=tm[1,2]
Q[0,1]=tm[1,1]
Q[1,0]=tm[1,0]
with open(root+'\\'+name[:-4]+'-P.map', 'w') as outfile:
   for ii in range (P.size):
       outfile.write('{0:4.10e}\n'.format(np.hstack(P)[ii]))
   for ii in range (Q.size):
       outfile.write('{0:4.10e}\n'.format(np.hstack(Q)[ii]))

#saving to .coeff file:
with open(root+'\\'+name[:-4]+'-P.coeff', 'w') as outfile:
    outfile.write('{0:4.10e}\n'.format(transform_matrix[0,2]+256))
    outfile.write('{0:4.10e}\n'.format(transform_matrix[0,0]))
    outfile.write('{0:4.10e}\n'.format(transform_matrix[0,1]))
    outfile.write('{0:4.10e}\n'.format(transform_matrix[1,2]))
    outfile.write('{0:4.10e}\n'.format(transform_matrix[1,0]))
    outfile.write('{0:4.10e}\n'.format(transform_matrix[1,1]))
                    
       
    