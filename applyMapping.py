# -*- coding: utf-8 -*-
"""
Created on Jun 2019
@author: carlos
meant for mapping, but now with class implementation
"""
import os
from image_adapt.ImageCollection import ImageCollection
from image_adapt. Mapping import Mapping
import numpy as np
import matplotlib.pyplot as plt
import time

if 0: # use for trial on loadig data
    tetra_fn=os.path.normpath('N:/tnw/BN/CMJ/Shared/Margreet/20181203_HJ_training/mapping/rough_ave.tif')
    root=os.path.normpath("N:/tnw/BN/CMJ/Shared/Margreet/20181203_HJ_training/#3.10_0.1mg-ml_streptavidin_50pM_HJC_G_movies")
    name='hel4.pma'
elif 0: # remove -P.pks and -P.traces in the root folder, to test analysis on pms file
    tetra_fn=os.path.normpath('N:/tnw/BN/CMJ/Shared/Margreet/20181203_HJ_training/mapping/rough_ave.tif')
    tmp=Mapping(tetra_fn)
    root=os.path.normpath(r"N:\tnw\BN\CMJ\Shared\Margreet\20181203_HJ_training\#1.10_0.1mg-ml_streptavidin_50pM_HJA_G_movies")
    name='hel2.pma'                          
elif 1:
    tetra_fn=os.path.normpath(r'M:\tnw\bn\cmj\Shared\margreet\Cy3 G50\cmos tetra G2,4 v2.tif')
    tmp=Mapping(tetra_fn,5000) # if you want to test the mapping seperately, it is included in ImageCollection
    root=os.path.normpath(r'M:\tnw\bn\cmj\Shared\margreet\Cy3 G50')
    name='Spooled files.sifx'
else:
    tetra_fn=os.path.normpath(r'N:\tnw\BN\CMJ\Shared\Margreet\181218 - First single-molecule sample (GattaQuant)\Signal  - Image Series-16 bit.tif')
    tmp=Mapping(tetra_fn,5000) # if you want to test the mapping seperately, it is included in ImageCollection
    root=os.path.normpath("N:/tnw/BN/CMJ/Shared/Margreet/181218 - First single-molecule sample (GattaQuant)/RawData") 
    name='Spooled files.sifx'
    # 
image_fn=os.path.join(root,name)
imc = ImageCollection(tetra_fn, image_fn)  # **kwargs, 'pks_fn'= ; choice_channel='d'
#

###
if 1:
    from image_adapt.load_file import read_header
    hdim,vdim,nImages,A=read_header(image_fn)
    donor,acceptor=imc.get_all_traces(A) # 40 seconds for 166 traces from 200 images of 512x512, reduced to 24 sec
    ##donor,acceptor=imc.get_all_traces_sifx_fast()
#
##example codes to return data from mapping, image collection
#    img = imc.get_image(0,A)
#    #tmp._tf2_matrix
#    #imc.mapping._tf2_matrix
#    #imc.pts_number # is too 19 for the pma
#    #import matplotlib.pyplot as plt; plt.imshow(imc.im_mean20_correct)
#
#imc.get_image(1,A)
#imc.get_image(2,A)
#imc.get_image(3,A)
#imc.get_image(4,A)
#plt.figure(10), plt.plot(donor[:,1478])
##
if 0:
    for ii in range(imc.pts_number):
        if np.amax(donor[:,ii])>200:
            plt.plot(donor[:,ii])
            plt.ylim([0, 300])
            plt.title([ii,imc.ptsG[ii,:]])
            plt.pause(0.001)
            plt.waitforbuttonpress()
            plt.cla()