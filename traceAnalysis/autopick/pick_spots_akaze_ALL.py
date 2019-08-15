# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 13:08:00 2019

@author: mwdocter

find spots and make transformation, see Shirani's report
version 4

warning: black box numbers: manual input: 4 points, minimum # points to be found: 50 per 512*512 pixels
"""
# # NOTES # #
# Install the packages opencv and opencv-contrib.
# Both versions of OpenCV have to be lower than 3.4.3.
# This is because SIFT and SURF algorithms are both patented and removed from newer versions
# Range of uint16: [0, 65535]
# code via Shirani Bisnajak (BEP 2018-2019)

#https://opencv-python-tutroals.readthedocs.io/en/latest/py_tutorials/py_feature2d/py_matcher/py_matcher.html
from autopick.do_before import clear_all
clear_all()
import cv2 #computer vision?
#from PIL import Image # Python Imaging Library (PIL)
import tifffile as TIFF
import matplotlib.pyplot as plt
#from skimage import filters #  Image processing in Python — scikit-image
import numpy as np
import bisect #This module provides support for maintaining a list in sorted order without having to sort the list after each insertion.
from image_adapt.find_threshold import remove_background, get_threshold
import os
from image_adapt.polywarp import polywarp, polywarp_apply
from image_adapt.distance_calculation import dist_calc
#import image_adapt

def imghist(img): 
#img is a np.array, imghist makes a histogram
    binrange = [np.min(img), np.max(img)]
    binlength = binrange[1] - binrange[0]
    hist,bins = np.histogram(img.flatten(),binlength, binrange) #img.flatten changes RGB into one channel
    cdf = hist.cumsum()
    cdf_normalized = cdf * hist.max()/ cdf.max()
    plt.plot(cdf_normalized, color = 'b')
    plt.hist(img.flatten(), binlength, binrange, color = 'r')
    plt.xlim(binrange)
    plt.legend(('cdf','histogram'), loc = 'upper left')
    plt.show()
    
    
def imadjust(src, tol=1, vout=(0,255)): 
#imadjust scales the gray scale in the image to 0-255
# src : input one-layer image (numpy array)
# tol : tolerance, from 0 to 100.
# vin  : src image bounds
# vout : dst image bounds
# return : output img

    assert len(src.shape) == 2 ,'Input image should be 2-dims'

    tol = max(0, min(100, tol))

    vin = [np.min(src), np.max(src)]
    vout = [0, 65535] # 65535=16 bits
    if tol > 0:
# Compute in and out limits
# Histogram
        hist = np.histogram(src,bins=list(range(vin[1] - vin[0])),range=tuple(vin))[0]

# Cumulative histogram
        cum = hist.copy()
        for i in range(0, vin[1]-vin[0]-1): cum[i] = cum[i - 1] + hist[i] # why not hist.cumsum() here?

# Compute bounds
        total = src.shape[0] * src.shape[1]
        low_bound = total * tol / 100
        upp_bound = total * (100 - tol) / 100
        vin[0] = bisect.bisect_left(cum, low_bound)
        vin[1] = bisect.bisect_left(cum, upp_bound)
        if vin[0]==0 and vin[1]==0: #do not allow vin=0,0
             vin = [np.min(src), np.max(src)]
# Stretching
    scale = (vout[1] - vout[0]) / (vin[1] - vin[0])
    vs = src-vin[0]
    vs[src<vin[0]]=0 #everything below zero becomes 0
    vd = vs*scale+0.5 + vout[0] # why +0.5?
    vd[vd>vout[1]] = vout[1]
    dst = vd

    return dst.astype(np.uint16), scale


def im_binarize(img, f): 
# makes all values below zero 0
    temp = img.copy()
    temp[temp<f] = 0
    return temp.astype(np.uint8)


def enhance_blobies(image, f, tol): 
# sequential scaling (imadjust) and removing <0 (imbinarize)
    l, r = image[:, :image.shape[1]//2], image[:, image.shape[1]//2:]
    l_adj, scale_l,r_adj, scale_r = imadjust(l.copy(), tol), imadjust(r.copy(), tol)
    l_bin, r_bin = im_binarize(l_adj, f*scale_l).astype(np.uint8), im_binarize(r_adj,f*scale_r).astype(np.uint8)
    return l, r, l_bin, r_bin

def enhance_blobies_single(image, f, tol): 
# sequential scaling (imadjust) and removing <0 (imbinarize)
    l_adj,scale = imadjust(image.copy(), tol)
    l_bin = im_binarize(l_adj, f*scale).astype(np.uint8)
    return l_bin

def mapping_manual(file_tetra, show=0, f=None,bg=None, tol=1):
# let you pick four points in left image, adjust the position with zoomed in versions, 
# and calculates Transformation matrix (to be later adjusted in automatic mapping)
#'E:\CMJ trace analysis\\autopick\\tetraspeck.tif'
    import matplotlib.pyplot as plt
    plt.close('all') 
    print('manual aligning')
    root, name = os.path.split(file_tetra)
    save_fn=os.path.join(root,name[:-4]+'-P.coeff') 
    
    global points_right # globals in nested function, need to be defined here as well
    global points_left
    global ii
# Open image
    if type(file_tetra)==str:
        import re 
        from image_adapt.load_file import read_one_page, read_header
        if re.search('tif',file_tetra)!=None: 
            image_tetra_raw = TIFF.imread(file_tetra)
        elif re.search('sifx',file_tetra)!=None: 
            hdim, vdim, n_images,A = read_header(file_tetra)
            im_array = np.dstack([(read_one_page(file_tetra, pageNb=jj,A=A,ii=np.array(range(1024*1024)))).astype(float) for jj in range(20)])
            image_tetra_raw = np.mean(im_array, axis=2).astype(np.uint16)
            TIFF.imwrite(file_tetra[:-5]+'.tif',image_tetra_raw.astype(np.uint16))
    else: #assume you are passing an image
        image_tetra_raw=file_tetra
        
# default threshold for 16 bits 50000, for 8 bits 200 (=256*50000/64000)
    if f==None:
            f=50000
        
# calculate background if not supplied. different value for dondor&acceptor channel
    if type(bg)==type(None): # take two different backgrounds, one for donor, one for acceptor channel
            sh=np.shape(image_tetra_raw)
            thr_donor=get_threshold(image_tetra_raw[:,1:sh[0]//2])
            thr_acceptor=get_threshold(image_tetra_raw[:,sh[0]//2:])
            bg=np.zeros(sh)
            bg[:,1:sh[0]//2]=thr_donor
            bg[:,sh[0]//2:]=thr_acceptor
    image_tetra=remove_background(image_tetra_raw.astype(float),bg)
    image_tetra=image_tetra.astype(np.uint16)    
    position1=[]
    position2=[]
# for later keypoint detection, you need a threshold, which depends on the image
#start f depends on max (image_tetra), only interesting if you are going to not load both manual&automatic alignments
    if not( os.path.isfile(save_fn) and os.path.isfile(os.path.join(root,name[:-4]+'-P.map')) ):
        position1 = image_adapt.analyze_label.analyze(image_tetra[:, :image_tetra.shape[1]//2])[2]  
        position2=  image_adapt.analyze_label.analyze(image_tetra[:, image_tetra.shape[1]//2:])[2]  
        fL=get_threshold(image_tetra[:, :image_tetra.shape[1]//2])
        fR=get_threshold(image_tetra[:, image_tetra.shape[1]//2:])
#
        gray1= enhance_blobies_single(image_tetra[:, :image_tetra.shape[1]//2],fL, tol)
##                # initialize the AKAZE descriptor, then detect keypoints and extract
##                # local invariant descriptors from the image
#        detector = cv2.AKAZE_create()
#        (kps1, descs1) = detector.detectAndCompute(gray1, None)
#        position1=cv2.KeyPoint_convert(kps1);
#
        gray2= enhance_blobies_single(image_tetra[:, image_tetra.shape[1]//2:],fR, tol)
##                # initialize the AKAZE descriptor, then detect keypoints and extract
##                # local invariant descriptors from the image
#        detector = cv2.AKAZE_create()
#        (kps2, descs2) = detector.detectAndCompute(gray2, None)
#        position2=cv2.KeyPoint_convert(kps2);
        
    transformation_matrixC=np.zeros((3,3))
    transformation_matrixC[2,2]=1
# if the transformation matrix has been produced previously, and saved to .coeff file, load it    
    if os.path.isfile(save_fn):
        #import from .coeff file:
        transformation_matrixC=load_coeff(save_fn) 
        # also load the detected points            
        points_left,points_right=load_coeffpoints(os.path.join(root,name[:-4]+'-P.coeffpoints'),image_tetra.shape[1] )
        
        try: fL
        except NameError: fL,fR=f,f
        return  transformation_matrixC,points_right,points_left,fL,fR
#else calculate the l(linear) transformation matrix
    else:
    # click on 4 points in donor image
        from pynput.keyboard import Key, Listener
        import matplotlib.pyplot as plt
        
        fig = plt.figure(10, figsize=(18,9))
        ax1 = fig.add_subplot(1,2,1)  
        ax1.imshow(image_tetra_raw[:, :image_tetra.shape[1]//2],vmin=thr_donor-30,vmax=thr_donor+5)#ax1.imshow(gray1, vmax=50)
        for ii in range((np.amax(np.shape(position1)))): 
            ax1.plot(position1[ii][0],position1[ii][1],'wo',markerfacecolor='none', markersize=4)
  
      
        ax2 = fig.add_subplot(1,2,2) 
        ax2.imshow(image_tetra_raw[:, image_tetra.shape[1]//2:],vmin=thr_acceptor-30,vmax=thr_acceptor+5)#ax2.imshow(gray2,vmax=50)
        for ii in range((np.amax(np.shape(position2)))): 
            ax2.plot(position2[ii][0],position2[ii][1],'wo',markerfacecolor='none', markersize=4)
#$$$$$$ BLACK BOX NUMBER #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$        
        if len(gray2)<=512: 
            points_left=plt.ginput(4) # Python needs 4 points. NOTE: for better alignment use more =10  points
            ax1.set_title('click on bright spots in the left image \n select four corners')
        else: 
            points_left=plt.ginput(10)
            ax1.set_title('click on bright spots in the left image \n select TEN points ')
#$$$$$$ BLACK BOX NUMBER #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        points_right=points_left.copy()
        
    # mark/overlay all four points in the image
        xpL=[xx for xx,yy in points_left]
        ypL=[yy for xx,yy in points_left]
        lineL,=ax1.plot(xpL,ypL,markersize=10, c='w', marker='o', fillstyle='none', linestyle='none') 
                
        xpR=[xx for xx,yy in points_right]
        ypR=[yy for xx,yy in points_right]
        lineR,=ax2.plot(xpR,ypR,markersize=10, c='w', marker='o', fillstyle='none', linestyle='none') 
        fig.canvas.draw()
        fig.canvas.flush_events()
        plt.pause(1)
        
        ax1.set_title('')
        ax2.set_title('move point with arrows, press esc when it matches the location in other channel')
    
    # this definition is to let the position of the selected spot move with the key presses 
        def on_release(key):
            global points_right
            global points_left
            global ii
            if key == Key.up:
                points_right[ii]=[points_right[ii][0], points_right[ii][1]+1]
            elif key == Key.down:
                points_right[ii]=[points_right[ii][0], points_right[ii][1]-1]
            elif key == Key.right:
                points_right[ii]=[points_right[ii][0]+1, points_right[ii][1]]
            elif key == Key.left:
                points_right[ii]=[points_right[ii][0]-1, points_right[ii][1]]
            elif key == Key.f1:
                points_left[ii]=[points_left[ii][0], points_left[ii][1]+1]
            elif key == Key.f2:
                points_left[ii]=[points_left[ii][0], points_left[ii][1]-1]
            elif key == Key.f4:
                points_left[ii]=[points_left[ii][0]+1, points_left[ii][1]]
            elif key == Key.f3:
                points_left[ii]=[points_left[ii][0]-1, points_left[ii][1]]                   
            
            if key == Key.esc or (key==Key.up) or (key == Key.down) or (key == Key.right) or (key == Key.left) or (key == Key.left) or (key == Key.f1) or (key == Key.f2) or (key == Key.f3) or (key == Key.f4):
                # Stop listener
                return False
        
    # loop over 4 spots, move then to corresponding features in both channels, press esc to go to the next point
        matches=len(points_left)   
        for ii in range(matches):
            xpR_new=0*xpR
            ypR_new=0*ypR
            xpL_new=0*xpL
            ypL_new=0*ypL
            
            ax1.set_xlim(points_left[ii][0]-50, points_left[ii][0]+50)
            ax1.set_ylim(points_left[ii][1]-50, points_left[ii][1]+50)
            ax2.set_xlim(points_left[ii][0]-50, points_left[ii][0]+50)
            ax2.set_ylim(points_left[ii][1]-50, points_left[ii][1]+50)
            line1=ax1.plot(points_left[ii][0],points_left[ii][1],markersize=10, c='y', marker='+', fillstyle='none')      
            line2=ax2.plot(points_right[ii][0],points_right[ii][1],markersize=10, c='y', marker='+', fillstyle='none')      
            fig.canvas.draw()
            fig.canvas.flush_events()
            plt.pause(0.1)  
            while xpR!=xpR_new or ypR!=ypR_new or xpL!=xpL_new or ypL!=ypL_new:
                xpR=[xx for xx,yy in points_right]
                ypR=[yy for xx,yy in points_right]
                xpL=[xx for xx,yy in points_left]
                ypL=[yy for xx,yy in points_left]
                
                # Collect events until released
                with Listener(on_release=on_release) as listener:
                    listener.join()
            
                xpR_new=[xx for xx,yy in points_right]
                ypR_new=[yy for xx,yy in points_right]
                xpL_new=[xx for xx,yy in points_left]
                ypL_new=[yy for xx,yy in points_left]
                
                lineR.set_xdata(xpR_new)
                lineR.set_ydata(ypR_new)
                lineL.set_xdata(xpL_new)
                lineL.set_ydata(ypL_new)
                fig.canvas.draw()
                fig.canvas.flush_events()
                plt.pause(0.1)  
           
            line1[-1].remove()
            line2[-1].remove()
    
            fig.canvas.draw()
            fig.canvas.flush_events()
            plt.pause(0.1)  
               
    # find the transformation matrix with cv2. findhomography
        points_right=np.array(points_right)
        points_left=np.array(points_left)
        if show: print('the number of right points are  ' );    print(points_right)
        if show: print('the number of left points are  ' );     print(points_left)
        transformation_matrixC, mask = cv2.findHomography(points_right, points_left, cv2.RANSAC,20)
        #transformation_matrixC, mask = cv2.getPerspectiveTransform(points_right, points_left)         
    # produce an image in which the overlay between two channels is shown
        array_size=np.shape(gray2)
        imC=cv2.warpPerspective(gray2, transformation_matrixC, array_size[::-1] )
        
        if show:
            plt.figure(11, figsize=(18,9))
            plt.subplot(1,1,1)
            plt.subplot(1,6,1),
            plt.imshow(gray1, extent=[0,array_size[1],0,array_size[0]], aspect=1)
            plt.title('green channel')
                
            plt.subplot(1,6,2),
            plt.imshow(gray2, extent=[0,array_size[1],0,array_size[0]], aspect=1)
            plt.title('red channel')    
            
            plt.subplot(1,6,3),
            plt.imshow(imC, extent=[0,array_size[1],0,array_size[0]], aspect=1)
            plt.title('red transformed')
            plt.show()
        
            plt.subplot(1,6,4),
            A=(gray1>0)+2*(gray2>0)
            plt.imshow(A, extent=[0,array_size[1],0,array_size[0]], aspect=1)
           # plt.colorbar()
            plt.title( 'unaligned #(yellow) \nspots overlap {:d}'.format(np.sum(A==3))   ) 
                
            plt.subplot(1,6,5),
            AA=(gray1>0)+2*(imC>0)
            plt.imshow((gray1>0)+2*(imC>0), extent=[0,array_size[1],0,array_size[0]], aspect=1)
            #plt.colorbar()
            plt.title(  'manual align \n#spots overlap y{:d}'.format(np.sum(AA==3)) )    
            plt.show()
            plt.pause(0.05)
    
        plt.figure(10), plt.close()            
               
    #saving to .coeff file:
        save_coeff(save_fn,transformation_matrixC)
    # also save the found points
        
        save_coeffpoints(os.path.join(root,name[:-4]+'-P.coeffpoints'), points_left,points_right)
          
    return  transformation_matrixC,points_right,points_left, fL,fR


def mapping_automatic(file_tetra, tf1_matrix,show=0, fL=None,fR=None,bg=None, tol=1 ): #'E:\CMJ trace analysis\\autopick\\tetraspeck.tif'
# this function finds all spots, and alignes them with a polywarp function (higher order, not only linear transformation)    
    root, name = os.path.split(file_tetra)
    save_fn=os.path.join(root,name[:-4]+'-P.map') 
    print('automatic aligning')
    
# if auto mapping has been done, load the transformation matrix and found points, no need to recalculate
    if os.path.isfile(save_fn ):
        P,Q,transformation_matrixC=load_map(save_fn)
        
           
    #load found points, matched between channels: 
        #pts1 is location found in donor, pts2 is location found in acceptor channel
        #dstG is transformed pts1 location, to hopefully match with pts2 location
        
        pts1,pts2,dst2=load_mappoints(os.path.join(root,name[:-4]+'-P mappoints'))
        deg=len(P)-1
        P21,Q21=polywarp(pts1[:,0],pts1[:,1],pts2[:,0],pts2[:,1],degree=deg) 
    # load found positions (unmatched between channels)
        position1,position2=load_mapposition12(os.path.join(root,name[:-4]+'-P mapposition1'))

    else: 
# calculate the automatic mapping
    # Open image
        if type(file_tetra)==str:
            import re 
            from image_adapt.load_file import read_one_page, read_header
            if re.search('tif',file_tetra)!=None: 
                image_tetra_raw = TIFF.imread(file_tetra)
            elif re.search('sifx',file_tetra)!=None: 
                try: image_tetra_raw = TIFF.imread(file_tetra[:-5]+'.tif')
                except:
                    hdim, vdim, n_images,A = read_header(file_tetra)
                    im_array = np.dstack([(read_one_page(file_tetra, pageNb=jj,A=A,ii=np.array(range(1024*1024)))).astype(float) for jj in range(20)])
                    image_tetra_raw = np.mean(im_array, axis=2).astype(np.uint16)
                    TIFF.imwrite(file_tetra[:-5]+'.tif',image_tetra_raw.astype(np.uint16))
            
        else: #assume you are passing an image
            image_tetra_raw=file_tetra
        
    # take two different backgrounds, one for donor, one for acceptor channel
        if bg==None:
            sh=np.shape(image_tetra_raw)
            thr_donor=get_threshold(image_tetra_raw[:,1:sh[0]//2])
            thr_acceptor=get_threshold(image_tetra_raw[:,sh[0]//2:])
            bg=np.zeros(sh)
            bg[:,1:sh[0]//2]=thr_donor
            bg[:,sh[0]//2:]=thr_acceptor
        PL=plt.figure(1,figsize=(40,40)); plt. subplot(1,1,1)
        plt.imshow(image_tetra_raw, vmin=np.amin(image_tetra_raw), vmax=np.amin(image_tetra_raw)+200)
        PL.savefig(file_tetra[:-4]+'-P image_raw.tif')    
        
        image_tetra=remove_background(image_tetra_raw.astype(float),bg).astype(float)
        PL=plt.figure(1,figsize=(40,40)); plt. subplot(1,1,1)
        type(image_tetra)
        plt.imshow(image_tetra)
        plt.imshow(image_tetra, vmin=np.amin(image_tetra), vmax=np.amin(image_tetra)+50)
        PL.savefig(file_tetra[:-4]+'-P image_tetra.tif')      
        image_tetra=image_tetra.astype(np.uint16)    

        fL=get_threshold(image_tetra[:, :image_tetra.shape[1]//2])
        fR=get_threshold(image_tetra[:, image_tetra.shape[1]//2:])

        gray1= enhance_blobies_single(image_tetra[:, :image_tetra.shape[1]//2],fL, tol)
        gray2= enhance_blobies_single(image_tetra[:, image_tetra.shape[1]//2:],fR, tol) 
        
        position1 = image_adapt.analyze_label.analyze(image_tetra[:, :image_tetra.shape[1]//2])[2]  
        position2=  image_adapt.analyze_label.analyze(image_tetra[:, image_tetra.shape[1]//2:])[2]  

        for ii in range(np.amax(np.shape(position2))):  #adapt position2 to right channel asap
            position2[ii][0]=position2[ii][0]+ np.shape(image_tetra)[0]//2
            
        PL=plt.figure(1,figsize=(40,40)); plt. subplot(1,1,1)
        plt.imshow(image_tetra, vmin=np.amin(image_tetra), vmax=np.amin(image_tetra)+20)
        for ii in range((np.amax(np.shape(position1)))): 
            plt.plot(position1[ii,0],position1[ii,1], 'wo',markerfacecolor='none', markersize=5)
        for ii in range((np.amax(np.shape(position2)))): 
            plt.plot(position2[ii,0],position2[ii,1], 'ws',markerfacecolor='none', markersize=5)
        PL.savefig(file_tetra[:-4]+'-P positions.tif')
        if show: print('length position 1 is '); print(len(position1))
        if show: print('length position 2 is '); print(len(position2))
         
# find matching points (no other features) based on the manual mapping
        dst= cv2.perspectiveTransform(position1.reshape(-1,1,2),np.linalg.inv(tf1_matrix)) #reshape needed for persp.transform
        dst= dst.reshape(-1,2)
        for ii in range(np.amax(np.shape(dst))): #adapt dst to right channel asap
            dst[ii][0]=dst[ii][0]+ np.shape(image_tetra)[0]//2
        
        if show: print('shape position1 is '); print(np.shape(position1))   
        if show: print('shape position2 is '); print(np.shape(position2)) 
        if len(gray1)<=512: 
            LEN=4
        else:
            LEN=20 # for CMOS
        BB, dist=dist_calc(position2,dst,LEN)
        pts1, pts2 = [], []    
        for jj in range(len(BB)):
            pts1.append(position1[BB[jj,1]])
            pts2.append(position2[BB[jj,0]])
            
        if show: print([np.shape(ii),np.shape(jj)])
        pts1 = np.array(pts1).astype(np.float32)
        pts2 = np.array(pts2).astype(np.float32)             
        if 1: print('shape pts1 is '); print(np.shape(pts1))            
        if show: print('shape pts2 is '); print(np.shape(pts2))
# find the transformation matrix linear          
        transformation_matrixC,mask=cv2.findHomography(pts2,pts1, cv2.RANSAC,20)
        if show:
            # produce an image in which the overlay between two channels is shown
            array_size=np.shape(gray2)
            plt.imshow(gray2),
            if show: print('transformation matrix is '); print(transformation_matrixC)
            if show: print('array size is '); print(array_size)
            imC=cv2.warpPerspective(gray2, transformation_matrixC, array_size[::-1] )
            
            plt.figure(11,figsize=(18,9))
            plt.subplot(1,6,6),
            AA=(gray1>0)+2*(imC>0)
            plt.imshow((gray1>0)+2*(imC>0), extent=[0,array_size[1],0,array_size[0]], aspect=1)
            #plt.colorbar()
            plt.title(  'automatic align \n#spots overlap {:d}'.format(np.sum(AA==3))   )    
            plt.show()
            plt.pause(0.05)
            
        
        #the nonlinear polywarp
        deg=3
        if 1: #works well
            P,Q=polywarp(pts2[:,0],pts2[:,1],pts1[:,0],pts1[:,1],degree=deg) 
            P21,Q21=polywarp(pts1[:,0],pts1[:,1],pts2[:,0],pts2[:,1],degree=deg) 
        else:#backup if linear transformation was not correctly found
            P,Q=polywarp(pts2[:,0]-np.shape(image_tetra)[0]//2 ,pts2[:,1],pts1[:,0],pts1[:,1],degree=deg) 
            P21,Q21=polywarp(pts1[:,0],pts1[:,1],pts2[:,0]-np.shape(image_tetra)[0]//2 ,pts2[:,1],degree=deg)     
            
            P[0,0]=P[0,0]+np.shape(image_tetra)[0]//2
        # maybe later on you might find a need for orders higher than 3
            
        save_map(save_fn,P,Q)
        
        dst2=polywarp_apply(P,Q,pts1)# [np.sum([P[ii,jj]*pts1[kk,0]**ii * pts1[kk,1]**jj for ii in range(deg+1) for jj in range(deg+1)]) for kk in range(len(pts1))]
                                         #[np.sum([Q[ii,jj]*pts1[kk,0]**ii * pts1[kk,1]**jj for ii in range(deg+1) for jj in range(deg+1)]) for kk in range(len(pts1))]
            
# double check selected points are not within 10 pixels of the boundary               
        for ii in range(len(dst2)-1,-1,-1): # range(5,-1,-1)=5,4,3,2,1,0
            discard_pts1=pts1[ii,0]<10 or pts1[ii,1]<10 or pts1[ii,0]>np.shape(image_tetra)[0]//2-10 or pts1[ii,1]>np.shape(image_tetra)[1]-10
            discard_pts2=pts2[ii,0]<10+np.shape(image_tetra)[0]//2 or pts2[ii,1]<10 or pts2[ii,0]>np.shape(image_tetra)[0]-10 or pts2[ii,1]>np.shape(image_tetra)[1]-10
            discard_dst2=dst2[ii,0]<10+np.shape(image_tetra)[0]//2 or dst2[ii,1]<10 or dst2[ii,0]>np.shape(image_tetra)[0]-10 or dst2[ii,1]>np.shape(image_tetra)[1]-10
            discard=discard_dst2+discard_pts1+discard_pts2
            if discard:
                pts1=np.delete(pts1,ii, axis=0)
                pts2=np.delete(pts2,ii, axis=0)
                dst2=np.delete(dst2,ii, axis=0)

##actual locations of pts2,ptsG,position2 are in right channel, so add half of the horizontal dimension
#        for ii in range(np.amax(np.shape(pts2))): 
#            pts2[ii][0]=pts2[ii][0]+ np.shape(image_tetra)[0]//2
#            dst2[ii][0]=dst2[ii][0]+ np.shape(image_tetra)[0]//2

                         
# save selected points to file            
        save_mappoints(os.path.join(root,name[:-4]+'-P mappoints'),pts1,pts2,dst2)
#save all found positions            
        save_mapposition12(os.path.join(root,name[:-4]+'-P mapposition1'),position1,position2 )  
            
# plot markers around the found spots     and save to dile
        if not(os.path.isfile(file_tetra[:-4]+'-P selected.tif')):
              # Open image
                if not('image_tetra' in locals()): # another way of try:image_tetra, except NameError:
                    if type(file_tetra)==str:
                        image_tetra = TIFF.imread(file_tetra)
                    else: #assume you are passing an image
                        image_tetra=file_tetra
                
                PL=plt.figure(1,figsize=(40,40)); plt. subplot(1,1,1)
                #plt.imshow(image_tetra, vmin=np.amin(image_tetra), vmax=np.amin(image_tetra)+20)
                if show: print(pts1)
                for ii in range((np.amax(np.shape(pts1)))): 
                    plt.plot(pts1[ii][0],pts1[ii][1], 'wo',markerfacecolor='none', markersize=8)
                                
                if show:
                    print('number of found matching points is '); print(np.shape(pts1))
                    print('pts1 is '); print(pts1)
                    print('number of found matching points is '); print(np.shape(pts2))
                    print('pts2 is '); print(pts2)
                    print('number of found matching points is '); print(np.shape(dst2))
                    print('dstG is '); print(dst2)
          
#                pts2 = np.array([[ii[0] , ii[1]] for ii in pts2])
#                dst2 = np.array([[ii[0] , ii[1]] for ii in dst2])
                
                for ii in range(np.amax(np.shape(pts1))): 
                    plt.plot(pts2[ii][0],pts2[ii][1], 'ws',markerfacecolor='none', markersize=8)
                    plt.plot(dst2[ii][0],dst2[ii][1], 'wd',markerfacecolor='none', markersize=8)
                    
                PL.savefig(file_tetra[:-4]+'-P selected.tif')
        else:  # or load the image if it already exists
                im = TIFF.imread(file_tetra[:-4]+'-P selected.tif')
                PL=plt.figure(1), plt. subplot(1,1,1)
                plt.imshow(im)
            
    #print('done with automatic align')
    
    return  P,Q, transformation_matrixC, position1, position2, pts1, pts2, dst2, P21, Q21

def load_coeff(filename):
    tf1_matrix =np.zeros((3,3))
    tf1_matrix[2,2]=1
    with open(filename, 'r') as infile:
        tf1_matrix[0,2]    =float(infile.readline())
        tf1_matrix[0,0]    =float(infile.readline())
        tf1_matrix[0,1]    =float(infile.readline())
        tf1_matrix[1,2]    =float(infile.readline())
        tf1_matrix[1,0]    =float(infile.readline())
        tf1_matrix[1,1]    =float(infile.readline()) 
    return tf1_matrix    

def save_coeff(filename,transformation_matrixC):  
    #saving to .coeff file:
    with open(filename, 'w') as outfile:
        outfile.write('{0:4.10e}\n'.format(transformation_matrixC[0,2]))
        outfile.write('{0:4.10e}\n'.format(transformation_matrixC[0,0]))
        outfile.write('{0:4.10e}\n'.format(transformation_matrixC[0,1]))
        outfile.write('{0:4.10e}\n'.format(transformation_matrixC[1,2]))
        outfile.write('{0:4.10e}\n'.format(transformation_matrixC[1,0]))
        outfile.write('{0:4.10e}\n'.format(transformation_matrixC[1,1]))
            
def load_coeffpoints(filename,LEN):
    with open(filename) as infile:
        if LEN<=512: rm=4
        else: rm=10
        points_right=np.zeros((rm,2))   
        points_left=np.zeros((rm,2))
        for ii in range(0,rm):
            A=infile.readline().split()
            points_left[ii,0]=float(A[0])
            points_left[ii,1]=float(A[1])
            A=infile.readline().split()
            points_right[ii,0]=float(A[0])
            points_right[ii,1]=float(A[1])
    return points_left, points_right

def save_coeffpoints(filename, points_left,points_right):
     with open (filename,'w' )    as outfile:
            for ii in range (len(points_left)):
               outfile.write('{0:4.10e} {1:4.10e}\n'.format(points_left[ii,0],points_left[ii,1]))
               outfile.write('{0:4.10e} {1:4.10e}\n'.format(points_right[ii,0],points_right[ii,1]))
        
def load_map(filename):
    P=np.zeros((4,4))   
    Q=np.zeros((4,4))         
    tf2_matrix=np.zeros((3,3))
    with open(filename, 'r') as infile:
        for ii in range(0,16):
            P[ii//4,ii%4]=float(infile.readline())
        for ii in range(0,16):
            Q[ii//4,ii%4]=float(infile.readline())
                
    tf2_matrix[0,2]=P[0,0]
    tf2_matrix[0,1]=P[0,1]
    tf2_matrix[0,0]=P[1,0]
    tf2_matrix[1,2]=Q[0,0]
    tf2_matrix[1,1]=Q[0,1]
    tf2_matrix[1,0]=Q[1,0]
    tf2_matrix[2,2]=1
    tf2_matrix=np.linalg.inv(tf2_matrix)   ## is this necessary?  
    return P,Q,tf2_matrix
    
def save_map(filename,P,Q):
    with open(filename, 'w') as outfile:
           for ii in range (P.size):
               outfile.write('{0:4.10e}\n'.format(np.hstack(P)[ii]))
           for ii in range (Q.size):
               outfile.write('{0:4.10e}\n'.format(np.hstack(Q)[ii]))
        
def load_mappoints(filename):
    with open(filename , 'r') as infile:
            pts1=[]  
            pts2=[]
            dst2=[]
            for ii in range(0,20000): ## is 20000 enough?
                ii
                try:
                    A=infile.readline().split()
                    pts1.append([float(A[0]),float(A[1])])
                    A=infile.readline().split()
                    pts2.append([float(A[0]),float(A[1])])
                    A=infile.readline().split()
                    dst2.append([float(A[0]),float(A[1])])
                except:
                    break
    pts1=  np.array(pts1)      
    pts2=  np.array(pts2)      
    dst2=  np.array(dst2) 
    return pts1,pts2,dst2          
        
def save_mappoints(filename,pts1,pts2,dst2):
    with open (filename,'w' )    as outfile:
            for ii in range (len(pts1)):
               outfile.write('{0:4.10e} {1:4.10e}\n'.format(pts1[ii,0],pts1[ii,1]))
               outfile.write('{0:4.10e} {1:4.10e}\n'.format(pts2[ii,0],pts2[ii,1]))
               outfile.write('{0:4.10e} {1:4.10e}\n'.format(dst2[ii,0],dst2[ii,1]))
               
def load_mapposition12(filename):
    with open (filename,'r' )    as infile:
            position1=[]  
            for ii in range (20000):
                try:
                    A=infile.readline().split()
                    position1.append([float(A[0]),float(A[1])])
                except:
                  break
    with open (filename[:-1]+'2','r' )    as infile:
            position2=[]  
            for ii in range (20000):
                try:
                    A=infile.readline().split()
                    position2.append([float(A[0]),float(A[1])])
                except:
                  break   
    return position1, position2   

def save_mapposition12(filename, position1, position2):
    with open (filename,'w' )    as outfile:
            for ii in range (len(position1)):
               outfile.write('{0:4.10e} {1:4.10e}\n'.format(position1[ii,0],position1[ii,1]))
    with open (filename[:-1]+'2','w' )    as outfile:
            for ii in range (len(position2)):
               outfile.write('{0:4.10e} {1:4.10e}\n'.format(position2[ii,0],position2[ii,1]))              