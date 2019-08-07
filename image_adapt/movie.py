# -*- coding: utf-8 -*-
"""
Created on Jun 2019
@author: carlos

some features are shared within an image collection, like the file name, mapping, background, threshold, ...)
so imagecollection is a class ;)

warning blackboax number: (size+fwhm) gauss for extracting donor&acceptor
"""
from image_adapt.load_file import read_one_page#_pma, read_one_page_tif
from image_adapt.load_file import read_header
from image_adapt.rolling_ball import rollingball

from image_adapt.find_threshold import remove_background
from image_adapt.find_threshold import get_threshold
import matplotlib.pyplot as plt
import numpy as np
#from cached_property import cached_property
from image_adapt.Mapping import Mapping
from image_adapt.Image import Image
import image_adapt.analyze_label # note analyze label is differently from the approach in pick spots
#import cv2
import os
from find_xy_position.Gaussian import makeGaussian
import time
from image_adapt.polywarp import polywarp, polywarp_apply
 
class Movie(object):
    def __init__(self, filepath):#, **kwargs):
        self.filepath = filepath
       
        #self.set_background_and_transformation() # Carlos: isn't this double, you call set_background twice, Margreet: this is the bg of the image, not the tetra_image. Same formulas though
        
        #Ivo: This is a good idea I think, but I think we should put this option in the method that produces the pks file.
        #self.pks_fn=kwargs.get('pks_fn',image_fn) # default pks_fn=image_fn, but you could give in a different name; ImageCollection(name, name, pks_fn='   ')
        
        #Ivo: same here
        #self.choice_channel=kwargs.get('choice_channel','d') #default 'd', 'd' for donor, 'a' for acceptor, 'da' for the sum
        
        #Ivo: what are these?
        #self.generic_map=kwargs.get('generic',0)
        #self.ii=np.array(range(1024*1024))
    
        #Ivo: Commented this because I would like to be able to instantiate the object without doing this initially
        #self.mapping = Mapping(tetra_fn,generic=self.generic_map)
#        (self.background,
#         self.pts_number,
#         self.dstG,
#         self.ptsG,
#         self.im_mean20_correct,
#         self.n_images,
#         self.hdim,
#         self.vdim,
#         self.Gauss,
#         self.ptsG2) = self.set_background_and_transformation()
#        self.show_selected_spots()
 
#    @cached_property
#    def read_one_page(self):
#        return read_one_page ##normally this funciton needs two inputs, why not here?
#        if '.pma' in self.image_fn:
#            return read_one_page_pma
#        elif '.tif' in self.image_fn:
#            return read_one_page_tif

    def __repr__(self):
        return(f'{self.__class__.__name__}({str(self.filename)})')

    def set_background_and_transformation(self):
        """
        sum 20 image, find spots& background. Then loop over all image, do background subtract+ extract traces
        :return:
        """
        hdim, vdim, n_images,self.A = read_header(self.image_fn)
        im_array = np.dstack([(read_one_page(self.image_fn, pageNb=jj,A=self.A,ii=self.ii)).astype(float) for jj in range(20)])
        im_mean20 = np.mean(im_array, axis=2).astype(int)
        if 0:
            bg = rollingball(im_mean20,hdim/10)[1] # this one is not used in pick_spots_akaze
            im_mean20_correct = im_mean20 - bg
            im_mean20_correct[im_mean20_correct < 0] = 0       
            threshold = get_threshold(im_mean20_correct)
            im_mean20_correct=remove_background(im_mean20_correct,threshold)  
        else:
            sh=np.shape(im_mean20)
            thr_donor=get_threshold(im_mean20[:,1:sh[0]//2])
            thr_acceptor=get_threshold(im_mean20[:,sh[0]//2:])
            bg=np.zeros(sh)
            bg[:,1:sh[0]//2]=thr_donor
            bg[:,sh[0]//2:]=thr_acceptor
            im_mean20_correct=remove_background(im_mean20,bg)   
        
        #note: optionally a fixed threshold can be set, like with IDL
        # note 2: do we need a different threshold for donor and acceptor?
        
        root, name = os.path.split(self.pks_fn)
        pks_fn=os.path.join(root,name[:-4]+'-P.pks') 
        pks2_fn=os.path.join(root,name[:-4]+'-P2.pks') 
        if os.path.isfile(pks2_fn): 
        # if you can load the pks data, load it
             ptsG=[]
             dstG=[]
             with open(pks_fn, 'r') as infile:
                 for jj in range(0,10000): # there will be a time when more than 10000 frames are generated
                     A=infile.readline()
                     if A=='':
                         break
                     ptsG.append([float(A.split()[1]),float(A.split()[2])])
                     A=infile.readline()
                     dstG.append([float(A.split()[1]),float(A.split()[2])])
                     
             ptsG=np.array(ptsG)
             dstG=np.array(dstG)
             pts_number =len(ptsG)
             im_mean20_correctA=im_mean20_correct
             
             ptsG2=[]
             with open(pks2_fn, 'r') as infile:
                for jj in range(0,pts_number):
                    for jj in range(0,10000): # there will be a time when more than 10000 frames are generated
                     A=infile.readline()
                     if A=='':
                         break
                     ptsG2.append([float(A.split()[1]),float(A.split()[2])])
        else: 
        # if you cannot load the pks data, calculate it
            if self.pks_fn== self.image_fn: 
            # if they are the same, reuse im_mean20_correct; im_mean20_correctA is not stored/exported
                im_mean20_correctA=im_mean20_correct
            else: 
            #otherwise make your own im_mean correct for pks detection
                hdim, vdim, n_images,A = read_header(self.pks_fn)
                im_array = np.dstack([read_one_page(self.pks_fn, pageNb=jj,A=self.A,ii=self.ii).astype(float) for jj in range(20)])
                im_mean20 = np.mean(im_array, axis=2).astype(int)
                if 0: #older version, not matching pick spots
                    bg = rollingball(im_mean20)[1]
                    im_mean20_correctA = im_mean20 - bg
                    im_mean20_correctA[im_mean20_correctA < 0] = 0       
                    threshold = get_threshold(im_mean20_correctA)
                    im_mean20_correctA=remove_background(im_mean20_correctA,threshold)      
                else:
                    sh=np.shape(im_mean20)
                    thr_donor=get_threshold(im_mean20[:,1:sh[0]//2])
                    thr_acceptor=get_threshold(im_mean20[:,sh[0]//2:])
                    bg=np.zeros(sh)
                    bg[:,1:sh[0]//2]=thr_donor
                    bg[:,sh[0]//2:]=thr_acceptor
                    im_mean20_correctA=remove_background(im_mean20,bg)   
            if self.choice_channel=='d': 
            # with donor channel ptsG, calculate position in acceptor dstG
                pts_number, label_size, ptsG = image_adapt.analyze_label.analyze(im_mean20_correctA[:, 0:vdim//2])       
                ptsG2 = image_adapt.analyze_label.analyze(im_mean20_correctA[:, vdim//2:])[2] 
                ptsG2 = np.array([[ii[0] + hdim/2, ii[1]] for ii in ptsG2])

                dstG=polywarp_apply(self.mapping.P,self.mapping.Q,ptsG)
           #discard point close to edge image
                for ii in range(pts_number-1,-1,-1): # range(5,-1,-1)=5,4,3,2,1,0
                    discard_dstG=dstG[ii,0]<hdim//2-10 or dstG[ii,1]<10 or dstG[ii,0]>hdim-10 or dstG[ii,1]>vdim-10
                    discard_ptsG=ptsG[ii,0]<10 or ptsG[ii,1]<10 or ptsG[ii,0]>hdim/2-10 or ptsG[ii,1]>vdim-10
                    discard=discard_dstG+discard_ptsG
                    if discard:
                        ptsG=np.delete(ptsG,ii, axis=0)
                        dstG=np.delete(dstG,ii, axis=0)
                pts_number=len(ptsG) 
                
                print(pts_number)
            
            elif self.choice_channel=='a':
            # with acceptor dstG, calculate position in donor channel ptsG
                pts_number, label_size, dstG = image_adapt.analyze_label.analyze(im_mean20_correctA[:, vdim//2:])   
                ptsG2 = image_adapt.analyze_label.analyze(im_mean20_correctA[:, :vdim//2])[2]  
#                ptsG = cv2.perspectiveTransform(dstG.reshape(-1, 1, 2),(self.mapping._tf2_matrix))#transform_matrix))
#                ptsG = ptsG.reshape(-1, 2)
                ptsG=polywarp_apply(self.mapping.P21,self.mapping.Q21,dstG)
               
            #discard point close to edge image
                for ii in range(pts_number-1,-1,-1): # range(5,-1,-1)=5,4,3,2,1,0
                    discard_dstG=dstG[ii,0]<hdim//2-10 or dstG[ii,1]<10 or dstG[ii,0]>hdim-10 or dstG[ii,1]>vdim-10
                    discard_ptsG=ptsG[ii,0]<10 or ptsG[ii,1]<10 or ptsG[ii,0]>hdim/2-10 or ptsG[ii,1]>vdim-10
                    discard=discard_dstG+discard_ptsG
                    if discard:
                        ptsG=np.delete(ptsG,ii, axis=0)
                        dstG=np.delete(dstG,ii, axis=0)
                pts_number=   len(ptsG) 
                print(pts_number)
           
            elif self.choice_channel=='da':
                print('I have no clue yet how to do this')
                # most likely they do not overlap before finding transformation, so what is the point of doing D+A?
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
            
            root, name = os.path.split(self.pks_fn)
            pks_fn=os.path.join(root,name[:-4]+'-P2.pks') 
            with open(pks_fn, 'w') as outfile:
                for jj in range(len(ptsG2)):
                    pix0=ptsG2[jj][0]
                    pix1=ptsG2[jj][1]
                    outfile.write(' {0:4.0f} {1:4.4f} {2:4.4f} {3:4.4f} {4:4.4f} \n'.format((jj*2)+1, pix0, pix1, 0, 0, width4=4, width6=6))
        
#$$$$$$ BLACK BOX NUMBER #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$        
        sizeGauss=11
        ALL_GAUSS=makeGaussian(sizeGauss, fwhm=3, center=(sizeGauss//2, sizeGauss//2))          
#$$$$$$ BLACK BOX NUMBER #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
          
        return bg, pts_number, dstG, ptsG, im_mean20_correct, n_images, hdim, vdim, ALL_GAUSS, ptsG2

    def show_selected_spots(self):
    #make a plot with selected spots
        #make image with found spots
        PL=plt.figure(14,figsize=(40,40))    
        plt.imshow(self.im_mean20_correct,vmax=np.amin(self.im_mean20_correct)+5)
        for ii in range((np.amax(np.shape(self.dstG)))): 
            plt.plot(self.ptsG[ii][0],self.ptsG[ii][1], 'wo',markerfacecolor='none', markersize=8)
            plt.plot(self.dstG[ii][0],self.dstG[ii][1], 'wv',markerfacecolor='none', markersize=8)
        for ii in range((np.amax(np.shape(self.ptsG2)))):
            plt.plot(self.ptsG2[ii][0],self.ptsG2[ii][1],'y^',markerfacecolor='none', markersize=8)
        PL.savefig(self.image_fn[:-4]+'-P data found spots.tif')
        
        PL=plt.figure(15,figsize=(40,40))
        if self.choice_channel=='d': 
            for ii in range((np.amax(np.shape(self.ptsG2)))):
                plt.plot(self.ptsG2[ii][0]-len(self.im_mean20_correct)//2,self.ptsG2[ii][1], 'r^')#,markerfacecolor='none', markersize=8)
        else:
            for ii in range((np.amax(np.shape(self.ptsG2)))):
                plt.plot(self.ptsG2[ii][0],self.ptsG2[ii][1], 'r^')#,markerfacecolor='none', markersize=8)
        for ii in range((np.amax(np.shape(self.dstG)))): 
            plt.plot(self.ptsG[ii][0],self.ptsG[ii][1], 'ko',markerfacecolor='none', markersize=8)
            plt.plot(self.dstG[ii][0]-len(self.im_mean20_correct)//2,self.dstG[ii][1], 'kv',markerfacecolor='none', markersize=8)
       
        
        PL.savefig(self.image_fn[:-4]+'-P data location spots.tif')  
        
    def subtract_background(self, im):
        im_correct = im - self.background
        im_correct[im_correct < 0] = 0
        return remove_background(im_correct, self.threshold)

    def get_image(self, idx):
        img= read_one_page(self.image_fn, idx,self.A, self.ii)
        #img = self.subtract_background(img)
        return Image(img, self.vdim, self.mapping._tf2_matrix, self.ptsG, self.dstG, self.pts_number, self.Gauss)
    
    def get_image_show(self, idx, hs,ws,siz): # example hs=650,ws=950,siz=20
        img= read_one_page(self.image_fn, idx,self.A,self.ii)
        plt.figure(idx)
        ax1=plt.subplot(1,2,1)
        ax1.imshow(img)
        ax1.set_xlim(hs,hs+siz)
        ax1.set_ylim(ws, ws+siz)
        
        img = self.subtract_background(img)
        ax2=plt.subplot(1,2,2)
        ax2.imshow(img)
        ax2.set_xlim(hs,hs+siz)
        ax2.set_ylim(ws, ws+siz)
        return Image(img, self.vdim, self.mapping._tf2_matrix, self.ptsG, self.dstG, self.pts_number, self.Gauss)
    
    def get_all_traces(self):
    # reutnr donor and acceptor for the full data set
        root, name = os.path.split(self.image_fn)
        traces_fn=os.path.join(root,name[:-4]+'-P.traces') 
        Ncolours=2
        if os.path.isfile(traces_fn):
        # load if traces file already exist
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
        # go through all images, extract donor and acceptor signal
            donor=np.zeros(( self.n_images,self.pts_number))
            acceptor=np.zeros((self.n_images,self.pts_number))
           
            t0 = time.time()  
            for ii in range(self.n_images): #self.n_images also works for pm, len(self.A.filelist) not
                print(ii)
                img=self.get_image(ii)
                donor[ii,:]=img.donor # will multiply with gaussian, spot location is not drift compensated
                acceptor[ii,:]=img.acceptor # will multiply with gaussian, spot location is not drift compensated
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
    