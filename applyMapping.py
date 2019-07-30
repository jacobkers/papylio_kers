# -*- coding: utf-8 -*-
"""
Created on Jun 2019
@author: carlos
meant for mapping, but now with class implementation
"""
import os
from image_adapt.ImageCollection import ImageCollection
#from image_adapt. Mapping import Mapping
import numpy as np
import matplotlib.pyplot as plt
#import time

plt.close('all')
if 0: # use for trial on HJ example data for RP practicum
    tetra_fn=os.path.normpath('N:/tnw/BN/CMJ/Shared/Margreet/20181203_HJ_training/mapping/rough_ave.tif')
   # image_fn=os.path.normpath(r"N:\tnw\BN\CMJ\Shared\Margreet\20181203_HJ_training\#1.10_0.1mg-ml_streptavidin_50pM_HJA_G_movies\hel2.pma")
    # do NOT use #1 pma, there are hardly any anticorrelated traces to be found.
    image_fn=os.path.normpath("N:/tnw/BN/CMJ/Shared/Margreet/20181203_HJ_training/#3.10_0.1mg-ml_streptavidin_50pM_HJC_G_movies\hel4.pma")
   
    imc = ImageCollection(tetra_fn, image_fn)  # **kwargs, 'pks_fn'= ; choice_channel='d'

elif 0: # use for trial on CMOS data, only Cy3 here, so no anticorrelation
    ## tetra_fn=os.path.normpath(r'N:\tnw\BN\CMJ\Shared\Margreet\190628 HJ long TIR-Ivo\Cy3 G50\cmos tetra G2,4 v2.tif') 
    ## unfortunatelymapiing goes well for the above image, but this map does not seem to fit on the data below
    root=os.path.normpath(r'N:\tnw\BN\CMJ\Shared\Margreet\190628 HJ long TIR-Ivo\Cy3 G50')
    name='Spooled files.sifx'
    image_fn=os.path.join(root,name)
    tetra_fn=image_fn ##here having the same tetra as image file does not work
    imc = ImageCollection(tetra_fn, image_fn)  # **kwargs, 'pks_fn'= ; choice_channel='d'

elif 1: # this is a great data set: Holliday Junction (so test anticorrelation) on Ivo's setup
    tetra_fn=os.path.normpath(r'N:\tnw\BN\CMJ\Shared\Margreet\190628 HJ long TIR-Ivo\Cy3 G50\cmos tetra G2,4 v2.tif') 
    ## unfortunatelymapiing goes well for the aboe image, but this map does not seem to fit on the data below
    root=os.path.normpath(r'N:\tnw\BN\CMJ\Shared\Margreet\190628 HJ long TIR-Ivo\HJ 120mW_3_long')
    name='Spooled files.sifx'
    image_fn=os.path.join(root,name)
    tetra_fn=image_fn # used the data itself to do mappping, since for the Cy3 G50 this was also better.
    imc = ImageCollection(tetra_fn, image_fn)  # **kwargs, 'pks_fn'= ; choice_channel='d'

else:
    tetra_fn=os.path.normpath(r'N:\tnw\BN\CMJ\Shared\Margreet\20190620 - Demo Chirlmin\Mapping (old files)\rough_ave.tif')
    root=os.path.normpath(r'N:\tnw\BN\CMJ\Shared\Margreet\20190620 - Demo Chirlmin\HJ 50mM Mg 2x intensity (20mW)')
    name='hel20.pma'# or 21,22,23
    image_fn=os.path.join(root,name)
    imc = ImageCollection(tetra_fn, image_fn) 
    
###
if 1:
    donor,acceptor=imc.get_all_traces() # 18 seconds for 222 traces from 2000 images of 512x512
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
        im=imc.get_image(ii).image
        print(ii)
        naam=r'M:\tnw\bn\cmj\Shared\margreet\Cy3 G50\ModifiedData\Python'+'{:03d}'.format(ii)+'.tif'
        TIFF.imwrite(naam, np.uint16(im))
        
if 0:    #show traces for which a matching acceptor point is found
    BB=[]
    if imc.choice_channel=='d': 
        dist=np.zeros((len(imc.ptsG2),len(imc.dstG)))
        for ii in range(0, len(imc.ptsG2)):
                for jj in range(0, len(imc.dstG)):
                    dist[ii,jj]=np.sqrt((imc.ptsG2[ii][0]-imc.dstG[jj][0])**2+(imc.ptsG2[ii][1]-imc.dstG[jj][1])**2)
    elif imc.choice_channel=='a': 
        dist=np.zeros((len(imc.ptsG2),len(imc.ptsG)))
        for ii in range(0, len(imc.ptsG2)):
                for jj in range(0, len(imc.ptsG)):
                    dist[ii,jj]=np.sqrt((imc.ptsG2[ii][0]-imc.ptsG[jj][0])**2+(imc.ptsG2[ii][1]-imc.ptsG[jj][1])**2)
    for ii in range(0,len(imc.ptsG2)):
        jj=np.where(dist[ii,:]==min(dist[ii,:]))[0][0] # use only 0th element
        if dist[ii,jj]<4 and ii==np.where(dist[:,jj]==min(dist[:,jj]))[0][0]: 
            BB.append(jj)
    BB.sort()                
    plt.figure(20)

    
    for jj in range(len(BB)):
            ii=BB[jj]
            bg_add_subtr=1
            ###ii=220
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
            plt.pause(5)
#            plt.waitforbuttonpress(5)
if 1:
    plt.close('all')
    BB=range(0,imc.pts_number)
    selected=np.zeros(np.shape(BB))
#    BB=[]
#    for ii in range(imc.pts_number):
#        if np.amax(acceptor[:,ii])-np.amin(acceptor[:,ii])>150:
#            BB.append(ii)
    for jj in range(0,len(BB)):
            ii=BB[jj]
            
            ###ii=220
            for kk in range(1,5):
                plt.subplot(4,1,kk), plt.cla()
            
            plt.subplot(4,1,1)
            plt.plot(donor[:,ii]-np.amin(donor[:,ii]),'g')
            plt.plot(acceptor[:,ii]-np.amin(acceptor[:,ii]),'r')
            plt.xlim((0,650))
            plt.title([ii,imc.ptsG[ii,:]])
            
            plt.subplot(4,1,3)
            plt.plot(donor[:,ii],'g')
            plt.plot(acceptor[:,ii],'r')
            plt.xlim((0,650))
           
            plt.subplot(4,1,2)
            fret=(acceptor[:,ii]-np.amin(acceptor[:,ii]))/(donor[:,ii]+acceptor[:,ii]-np.amin(donor[:,ii])-np.amin(acceptor[:,ii]))
            plt.plot(fret,'b')
            plt.xlim((0,650))
            
            plt.subplot(4,1,4)
            fret= acceptor[:,ii] /(donor[:,ii]+acceptor[:,ii])
            plt.plot(fret)
            plt.xlim((0,650))
            plt.pause(0.1)
            A=input()
            try: selected[jj]=int(A)
            except: 
                if A=='x': break
            
            if np.mod(jj,100)==0:
                name='Spooled files.sifx'
                np.savetxt(os.path.join(root,'traces_select'), selected  )
                np.savetxt(os.path.join(root,'traces_select.txt'), selected  ,fmt='%1d')
     np.savetxt(os.path.join(root,'traces_select'), selected  )
     np.savetxt(os.path.join(root,'traces_select.txt'), selected  ,fmt='%1d')
             