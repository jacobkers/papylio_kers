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
    tmp=Mapping(tetra_fn)
    root=os.path.normpath("N:/tnw/BN/CMJ/Shared/Margreet/20181203_HJ_training/#3.10_0.1mg-ml_streptavidin_50pM_HJC_G_movies")
    name='hel4.pma'
    image_fn=os.path.join(root,name)
    imc = ImageCollection(tetra_fn, image_fn)  # **kwargs, 'pks_fn'= ; choice_channel='d'

elif 0: # remove -P.pks and -P.traces in the root folder, to test analysis on pms file
    tetra_fn=os.path.normpath('N:/tnw/BN/CMJ/Shared/Margreet/20181203_HJ_training/mapping/rough_ave.tif')
    root=os.path.normpath(r"N:\tnw\BN\CMJ\Shared\Margreet\20181203_HJ_training\#1.10_0.1mg-ml_streptavidin_50pM_HJA_G_movies")
    name='hel2.pma'                          
    image_fn=os.path.join(root,name)
    imc = ImageCollection(tetra_fn, image_fn)  # **kwargs, 'pks_fn'= ; choice_channel='d'

elif 0:
    tetra_fn=os.path.normpath(r'M:\tnw\bn\cmj\Shared\margreet\Cy3 G50\cmos tetra G2,4 v2.tif')
    root=os.path.normpath(r'M:\tnw\bn\cmj\Shared\margreet\Cy3 G50')
    name='Spooled files.sifx'
    image_fn=os.path.join(root,name)
    imc = ImageCollection(tetra_fn, image_fn)  # **kwargs, 'pks_fn'= ; choice_channel='d'

elif 0:
    tetra_fn=os.path.normpath(r'M:\tnw\bn\cmj\Shared\margreet\Cy3 G50\cmos tetra G2,4 v2.tif')
    root=os.path.normpath(r"N:\tnw\BN\CMJ\Shared\Margreet\190628 HJ long TIR-Ivo\HJ G50")
    name='Spooled files.sifx'
    image_fn=os.path.join(root,name)
    imc = ImageCollection(tetra_fn, image_fn)  # **kwargs, 'pks_fn'= ; choice_channel='d'

    
elif 0:
    tetra_fn=os.path.normpath(r'N:\tnw\BN\CMJ\Shared\Margreet\181218 - First single-molecule sample (GattaQuant)\Signal  - Image Series-16 bit.tif')
    root=os.path.normpath("N:/tnw/BN/CMJ/Shared/Margreet/181218 - First single-molecule sample (GattaQuant)/RawData") 
    name='Spooled files.sifx'
    image_fn=os.path.join(root,name)
    imc = ImageCollection(tetra_fn, image_fn)  # **kwargs, 'pks_fn'= ; choice_channel='d'
    
    # 
elif 1: 
    tetra_fn=os.path.normpath(r'E:\CMJ trace analysis\development_not_to_be_included_in_git\cmos tetra G2,4 v2.tif')
    tmp=Mapping(tetra_fn)
else: 
    tetra_fn=os.path.normpath(r'N:\tnw\BN\CMJ\Shared\Margreet\20181203_HJ_training\mapping\rough_ave.tif')
    tmp=Mapping(tetra_fn)
    
#

###
if 0:
    donor,acceptor=imc.get_all_traces() # 40 seconds for 166 traces from 200 images of 512x512, reduced to 24 sec
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
    for ii in range(1731,1732):#imc.pts_number):
        if np.amax(donor[:,ii])>200:
            plt.plot(donor[:,ii])
#            plt.ylim([0, 300])
            plt.title([ii,imc.ptsG[ii,:]])
            plt.pause(0.1)
            #plt.cla()
            
if 0:
    import numpy as np
    import tifffile as TIFF
    for ii in range(5):
        im=imc.get_image(ii,A).image
        print(ii)
        naam=r'M:\tnw\bn\cmj\Shared\margreet\Cy3 G50\ModifiedData\Python'+'{:03d}'.format(ii)+'.tif'
        TIFF.imwrite(naam, np.uint16(im))
        
if 0:
 
    BB=[]
    for ii in range(imc.pts_number):
        if np.amax(acceptor[:,ii])-np.amin(acceptor[:,ii])>100:
            BB.append(ii)
 
            bg_add_subtr=1
            ii=677
            plt.subplot(2,1,2), plt.cla()
            plt.subplot(2,1,1), plt.cla()
            plt.subplot(2,1,1),
            if bg_add_subtr:
                plt.plot(donor[:,ii]-np.amin(donor[:,ii]),'g')
                plt.plot(acceptor[:,ii]-np.amin(acceptor[:,ii]),'r')
            else:
                plt.plot(donor[:,ii],'g')
                plt.plot(acceptor[:,ii],'r')
 
#            plt.ylim([0, 300])
            plt.title([ii,imc.ptsG[ii,:]])
            plt.subplot(2,1,2)
            if bg_add_subtr:
                fret=(acceptor[:,ii]-np.amin(acceptor[:,ii]))/(donor[:,ii]+acceptor[:,ii]-np.amin(donor[:,ii])-np.amin(acceptor[:,ii]))
            else:
                fret= acceptor[:,ii] /(donor[:,ii]+acceptor[:,ii])
            plt.plot(fret,'b')
            plt.pause(0.1)
            plt.waitforbuttonpress(2)
            
         