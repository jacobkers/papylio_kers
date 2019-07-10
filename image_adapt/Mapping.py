# -*- coding: utf-8 -*-
"""
Created on Jun 2019
@author: carlos

a seperate file for doing mapping. 
this should be adapted to start with manual 3 points, then automated numerous points.
like the idl code, and save .coeff and .map files
"""
import autopick.pick_spots_akaze_ALL
import numpy as np
import os
import tifffile as tiff
import matplotlib.pyplot as plt
import cv2
from image_adapt.find_threshold import remove_background, get_threshold

class Mapping:#(object):
    def __init__(self, tetra_fname,f=10000, generic=0):
        self._tetra_fn = tetra_fname
        if generic==1:
            print('generic')
           
            self._tf1_matrix =np.zeros((3,3))
            self._tf1_matrix[2,2]=1
            save_fn=os.path.join(os.path.split(tetra_fname)[0],'generic.coeff') 
            with open(save_fn, 'r') as infile:
                self._tf1_matrix[0,2]    =float(infile.readline())-256
                self._tf1_matrix[0,0]    =float(infile.readline())
                self._tf1_matrix[0,1]    =float(infile.readline())
                self._tf1_matrix[1,2]    =float(infile.readline())
                self._tf1_matrix[1,0]    =float(infile.readline())
                self._tf1_matrix[1,1]    =float(infile.readline())
                
            P=np.zeros((4,4))   
            Q=np.zeros((4,4))         
            self._tf2_matrix=np.zeros((3,3))
            save_fn=os.path.join(os.path.split(tetra_fname)[0],'generic.map')  
            with open(save_fn, 'r') as infile:
                for ii in range(0,16):
                    P[ii//4,ii%4]=float(infile.readline())
                for ii in range(0,16):
                    Q[ii//4,ii%4]=float(infile.readline())
                
            self._tf2_matrix[0,2]=P[0,0]
            self._tf2_matrix[0,1]=P[0,1]
            self._tf2_matrix[0,0]=P[1,0]
            self._tf2_matrix[1,2]=Q[0,0]
            self._tf2_matrix[1,1]=Q[0,1]
            self._tf2_matrix[1,0]=Q[1,0]
            self._tf2_matrix[2,2]=1
            self._tf2_matrix=np.linalg.inv(self._tf2_matrix)
            
        else:
            self.manual_align()
            self.automatic_align()
        self.visualize()
                

    @property
    def tetra_fn(self):
        return self._tetra_fn

    def manual_align(self):
        self._tf1_matrix, self.points_right,self.points_left, self.fL,self.fR= autopick.pick_spots_akaze_ALL.mapping_manual(self._tetra_fn,
                                                                             show=1,
                                                                             bg=None,
                                                                             tol=0, f=None)
        
        print(self.fL,self.fR)
        return self._tf1_matrix
    
    #@tetra_fn.setter
    def automatic_align(self):
        self._tf2_matrix,self.position1, self.position2, self.pts1, self.pts2,self.dstG = autopick.pick_spots_akaze_ALL.mapping_automatic(self._tetra_fn,self._tf1_matrix,
                                                                        show=1,
                                                                        bg=None,
                                                                        tol=0, fL=self.fL,fR=self.fR)
        return self._tf2_matrix
    
    
    def visualize(self):
#        plt.close('all')
#        image_tetra_raw = tiff.imread(self._tetra_fn)
#        PL=plt.figure(1,figsize=(40,40)); plt. subplot(1,1,1)
#        plt.imshow(image_tetra_raw, vmin=np.amin(image_tetra_raw), vmax=np.amin(image_tetra_raw)+200)
#        PL.savefig(self._tetra_fn[:-4]+' VIS image_raw.tif')    
        
#        sh=np.shape(image_tetra_raw)
#        thr_donor=get_threshold(image_tetra_raw[:,1:sh[0]//2])
#        thr_acceptor=get_threshold(image_tetra_raw[:,sh[0]//2:])
#        bg=np.zeros(sh)
#        bg[:,1:sh[0]//2]=thr_donor
#        bg[:,sh[0]//2:]=thr_acceptor
#        image_tetra=remove_background(image_tetra_raw.astype(float),bg)
#        
#        PL=plt.figure(2,figsize=(40,40)); plt. subplot(1,1,1)
#        plt.imshow(image_tetra, vmin=np.amin(image_tetra), vmax=np.amin(image_tetra)+200)
#        dstGm = cv2.perspectiveTransform(self.points_right.reshape(-1, 1, 2), np.linalg.inv(self._tf1_matrix))#transform_matrix))
#        dstGm = dstGm.reshape(-1, 2)
#        for ii in range((np.amax(np.shape(self.points_left)))): 
#            plt.plot(self.points_left[ii][0],self.points_left[ii][1], 'wo',markerfacecolor='none', markersize=50)
#            plt.plot(self.points_right[ii][0],self.points_right[ii][1], 'ws',markerfacecolor='none', markersize=70)
#            plt.plot(self.points_right[ii][0]+len(image_tetra_raw)//2,self.points_right[ii][1], 'ws',markerfacecolor='none', markersize=70)
#            plt.plot(dstGm[ii][0]+len(image_tetra_raw)//2,dstGm[ii][1], 'wd',markerfacecolor='none', markersize=50)
#        PL.savefig(self._tetra_fn[:-4]+' VIS manual mapping.tif')   
            
       
#    
#        PL=plt.figure(3,figsize=(40,40)); plt. subplot(1,1,1)
#        plt.imshow(image_tetra, vmin=np.amin(image_tetra), vmax=np.amin(image_tetra)+20)
#        for ii in range((np.amax(np.shape(self.position1)))): 
#            plt.plot(self.position1[ii][0],self.position1[ii][1], 'wo',markerfacecolor='none', markersize=5)
#        for ii in range((np.amax(np.shape(self.position2)))): 
#            plt.plot(self.position2[ii][0]+len(image_tetra)//2,self.position2[ii][1], 'ws',markerfacecolor='none', markersize=5)
#        PL.savefig(self._tetra_fn[:-4]+' VIS positions.tif')
#
#       # dstG = cv2.perspectiveTransform(self.pts1.reshape(-1, 1, 2), np.linalg.inv(self._tf2_matrix))#transform_matrix))
#      #  dstG = dstG.reshape(-1, 2)
#        PL=plt.figure(4,figsize=(40,40)); plt. subplot(1,1,1)
#        plt.imshow(image_tetra, vmin=np.amin(image_tetra), vmax=np.amin(image_tetra)+20)
#        for ii in range((np.amax(np.shape(self.position1)))): 
#            plt.plot(self.position1[ii][0],self.position1[ii][1], 'wo',markerfacecolor='none', markersize=5)
#        for ii in range((np.amax(np.shape(self.position2)))): 
#            plt.plot(self.position2[ii][0]+len(image_tetra)//2,self.position2[ii][1], 'ws',markerfacecolor='none', markersize=5)
#        PL.savefig(self._tetra_fn[:-4]+' VIS positions.tif')
#        for ii in range((np.amax(np.shape(self.pts1)))): 
#            plt.plot(self.pts1[ii][0],self.pts1[ii][1], 'wo',markerfacecolor='none', markersize=8)
#        for ii in range(np.amax(np.shape(self.pts1))): 
#            plt.plot(self.pts2[ii][0]+len(image_tetra)//2,self.pts2[ii][1], 'ws',markerfacecolor='none', markersize=8)
#            plt.plot(self.dstG[ii][0]+len(image_tetra)//2,self.dstG[ii][1], 'wd',markerfacecolor='none', markersize=8)
#                    
#        PL.savefig(self._tetra_fn[:-4]+' VIS selected.tif')
#                
        PL=plt.figure(5,figsize=(20,40)); plt. subplot(1,1,1) 
        for ii in range((np.amax(np.shape(self.pts1)))): 
            plt.plot(self.pts1[ii][0],self.pts1[ii][1], 'ko',markerfacecolor='none', markersize=8)
            plt.plot(self.pts2[ii][0],self.pts2[ii][1], 'ks',markerfacecolor='none', markersize=8)
            plt.plot(self.dstG[ii][0],self.dstG[ii][1], 'kd',markerfacecolor='none', markersize=8)
        PL.savefig(self._tetra_fn[:-4]+' VIS select.tif')  
        PL.set_size_inches(5,10, forward=True)
        
        PL=plt.figure(6,figsize=(20,40)); plt. subplot(1,1,1) 
        for ii in range((np.amax(np.shape(self.position1)))): 
            plt.plot(self.position1[ii][0],self.position1[ii][1], 'ko',markerfacecolor='none', markersize=8)
        for ii in range((np.amax(np.shape(self.position2)))): 
            plt.plot(self.position2[ii][0],self.position2[ii][1], 'ks',markerfacecolor='none', markersize=6)
        PL.savefig(self._tetra_fn[:-4]+' VIS pos_overlap.tif')  
        PL.set_size_inches(5,10, forward=True)