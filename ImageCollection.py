# -*- coding: utf-8 -*-
"""
Created on Jun 2019
@author: carlos

some features are shared within an image collection, like the file name, mapping, background, threshold, ...)
so imagecollection is a class ;)
"""
from load_file import read_one_page#_pma, read_one_page_tif
from image_adapt.rolling_ball import rollingball

from image_adapt.find_threshold import remove_background
from image_adapt.find_threshold import get_threshold
import numpy as np
from cached_property import cached_property
from Mapping import Mapping
from Image import Image
import analyze_label
import cv2
import os
from find_xy_position.Gaussian import makeGaussian

class ImageCollection(object):
    def __init__(self, tetra_fn, image_fn, **kwargs):
        self.image_fn = image_fn
        self.mapping = Mapping(tetra_fn)
        #self.set_background_and_transformation() # Carlos: isn't this double, you call set_background twice, Margreet: this is the bg of the image, not the tetra_image. Same formulas though
        self.pks_fn=kwargs.get('pks_fn',image_fn) # default pks_fn=image_fn, but you could give in a different name; ImageCollection(name, name, pks_fn='   ')
        self.choice_channel=kwargs.get('choice_channel','d') #default 'd', 'd' for donor, 'a' for acceptor, 'da' for the sum

        (self.background,
         self.threshold,
         self.pts_number,
         self.dstG,
         self.ptsG,
         self.im_mean20_correct,
         self.n_images,
         self.hdim,
         self.Gauss) = self.set_background_and_transformation()
 
#    @cached_property
#    def read_one_page(self):
#        return read_one_page ##normally this funciton needs two inputs, why not here?
#        if '.pma' in self.image_fn:
#            return read_one_page_pma
#        elif '.tif' in self.image_fn:
#            return read_one_page_tif

    def set_background_and_transformation(self):
        """
        sum 20 image, find spots& background. Then loop over all image, do background subtract+ extract traces
        :return:
        """
        _, hdim, vdim, n_images = read_one_page(self.image_fn, pageNb=0)
        im_array = np.dstack([(read_one_page(self.image_fn, pageNb=ii)[0]).astype(float) for ii in range(20)])
        im_mean20 = np.mean(im_array, axis=2).astype(int)
        bg = rollingball(im_mean20)[1]
        im_mean20_correct = im_mean20 - bg
        im_mean20_correct[im_mean20_correct < 0] = 0       
        threshold = get_threshold(im_mean20_correct)
        im_mean20_correct=remove_background(im_mean20_correct,threshold)        

        
        #note: optionally a fixed threshold can be set, like with IDL
        # note 2: do we need a different threshold for donor and acceptor?
        
        root, name = os.path.split(self.pks_fn)
        pks_fn=os.path.join(root,name[:-4]+'-P.pks') 
        if os.path.isfile(pks_fn): # if you can load the pks data, load it
             ptsG=[]
             dstG=[]
             with open(pks_fn, 'r') as infile:
                 for jj in range(0,10000):
                     A=infile.readline()
                     if A=='':
                         break
                     ptsG.append([float(A.split()[1]),float(A.split()[2])])
                     A=infile.readline()
                     dstG.append([float(A.split()[1]),float(A.split()[2])])
             ptsG=np.array(ptsG)
             dstG=np.array(dstG)
             pts_number =len(ptsG)
            
        else: # if you cannot load the pks data, calculate it
            if self.pks_fn== self.image_fn: # if they are the same, reuse im_mean20_correct; im_mean20_correctA is not stored/exported
                im_mean20_correctA=im_mean20_correct
            else: #otherwise make your own im_mean correct for pks detection
                _, hdim, vdim, n_images = read_one_page(self.image_fn, pageNb=0)
                im_array = np.dstack([(read_one_page(self.image_fn, pageNb=ii)[0]).astype(float) for ii in range(20)])
                im_mean20 = np.mean(im_array, axis=2).astype(int)
                bg = rollingball(im_mean20)[1]
                im_mean20_correctA = im_mean20 - bg
                im_mean20_correctA[im_mean20_correctA < 0] = 0       
                threshold = get_threshold(im_mean20_correctA)
                im_mean20_correctA=remove_background(im_mean20_correctA,threshold)      
                
            if self.choice_channel=='d': # with donor channel, 
                pts_number, label_size, ptsG = analyze_label.analyze(im_mean20_correctA[:, 0:vdim//2])       
                dstG = cv2.perspectiveTransform(ptsG.reshape(-1, 1, 2), np.linalg.inv(self.mapping._tf2_matrix))#transform_matrix))
                dstG = dstG.reshape(-1, 2)
                dstG = np.array([[ii[0] + 256, ii[1]] for ii in dstG])
            
            elif self.choice_channel=='a':
                pts_number, label_size, dstG = analyze_label.analyze(im_mean20_correctA[:, vdim//2:])     
                ptsG = cv2.perspectiveTransform(dstG.reshape(-1, 1, 2),(self.mapping._tf2_matrix))#transform_matrix))
                ptsG = ptsG.reshape(-1, 2)
                ptsG = np.array([[ii[0] + 256, ii[1]] for ii in ptsG])
            
            elif self.choice_channel=='da':
                print('I have no clue yet how to do this')
                #pts_number, label_size, ptsG = analyze_label.analyze(im_mean20_correctA[:, 0:vdim//2+im_mean20_correctA[:, vdim//2:]])                                
            else:
                print('make up your mind, choose wisely d/a/da')
          
            #saving to pks file
            with open(pks_fn, 'w') as outfile:
                for jj in range(0,pts_number):
                         pix0=ptsG[jj][0]
                         pix1=ptsG[jj][1]
                         outfile.write(' {0:4.0f} {1:4.4f} {2:4.4f} {3:4.4f} {4:4.4f} \n'.format((jj*2)+1, pix0, pix1, 0, 0, width4=4, width6=6))
                         pix0=dstG[jj][0]
                         pix1=dstG[jj][1]
                         outfile.write(' {0:4.0f} {1:4.4f} {2:4.4f} {3:4.4f} {4:4.4f} \n'.format((jj*2)+2, pix0, pix1, 0, 0, width4=4, width6=6))
        
        ALL_GAUSS=makeGaussian(11, fwhm=3, center=(5, 5))          
          
        return bg, threshold, pts_number, dstG, ptsG, im_mean20_correct, n_images, hdim, ALL_GAUSS

    def subtract_background(self, im):
        im_correct = im - self.background
        im_correct[im_correct < 0] = 0
        return remove_background(im_correct, self.threshold)

    def get_image(self, idx):
        img, hdim, vdim, n_images = read_one_page(self.image_fn, pageNb=idx)
        img = self.subtract_background(img)
        return Image(img, vdim, self.mapping._tf2_matrix, self.ptsG, self.dstG, self.pts_number, self.Gauss)
    
    def get_all_traces(self):
        root, name = os.path.split(self.image_fn)
        traces_fn=os.path.join(root,name[:-4]+'-P.traces') 
        Ncolours=2
        if os.path.isfile(traces_fn):
             with open(traces_fn, 'r') as infile:
                 Nframes = np.fromfile(infile, dtype = np.int32, count = 1).item()
                 Ntraces = np.fromfile(infile, dtype = np.int16, count = 1).item()
                 rawData = np.fromfile(infile, dtype = np.int16, count = Ncolours*Nframes * Ntraces)
             orderedData = np.reshape(rawData.ravel(), (Ncolours, Ntraces//Ncolours, Nframes), order = 'F') 
             donor=orderedData[0,:,:]   
             acceptor=orderedData[1,:,:]
             donor=np.transpose(donor)
             acceptor=np.transpose(acceptor)
        else:
            donor=np.zeros(( self.n_images,self.pts_number))
            acceptor=np.zeros((self.n_images,self.pts_number))
            import time
            t0 = time.time()  
            for ii in range(0,self.n_images):
                img=self.get_image(ii)
                donor[ii,:]=img.donor
                acceptor[ii,:]=img.acceptor
            t1=time.time()
            elapsed_time=t1-t0; print(elapsed_time)    
            
            #root, name = os.path.split(self.image_fn)
            
            #if os.path.isfile(trace_fn):
               
            with open(traces_fn, 'w') as outfile:
                 off = np.array([self.n_images], dtype=np.int32)
                 off.tofile(outfile)
                 off = np.array([2*self.pts_number], dtype=np.int16)
                 off.tofile(outfile)
                 time_tr=np.zeros((self.n_images,2*self.pts_number))
                 Ncolours=2
                 for jj in range(2*self.pts_number//Ncolours):
                     time_tr[:,jj*2] = donor[:,jj]
                     time_tr[:,jj*2+1]=  acceptor[:,jj]
                 off = np.array((time_tr), dtype=np.int16)
                 off.tofile(outfile)
        
        return donor, acceptor