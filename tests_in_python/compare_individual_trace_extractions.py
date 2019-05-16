# -*- coding: utf-8 -*-
"""
Created on Tue May 14 11:12:14 2019

@author: mwdocter
"""
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

plt.close('all')
root='N:\\tnw\\BN\\CMJ\\Shared\\Margreet\\20181203_HJ_training\\#3.10_0.1mg-ml_streptavidin_50pM_HJC_G_movies'
name='hel4.pma'

fn_IDL='N:\\tnw\\BN\\CMJ\\Shared\\Margreet\\20181203_HJ_training\\mapping\\rough.map';
fn_Python='N:\\tnw\\BN\\CMJ\\Shared\\Margreet\\20181203_HJ_training\\#3.10_0.1mg-ml_streptavidin_50pM_HJC_G_movies\\hel4-P.map'
fn_Python_COEFF='N:\\tnw\\BN\\CMJ\\Shared\\Margreet\\20181203_HJ_training\\#3.10_0.1mg-ml_streptavidin_50pM_HJC_G_movies\\hel4-P.coeff'

## copied from applyMapping_v3
_,hdim,vdim,nImages=(read_one_page_pma(root,name, pageNb=0))
im_sum=(read_one_page_pma(root,name, pageNb=0 )[0]).astype(float)
for ii in range (1,20):
    im=(read_one_page_pma(root,name, pageNb=ii)[0]).astype(float)
    im[im<0]=0
    im_sum=im_sum+im
im_mean20=(im_sum/20).astype(int)
im_bg=rollingball(im_mean20)[1]
##

transform_matrix=np.zeros((3,3))
transform_matrix[2][2]=1
## read in coeff into transform_matrix
with open(root+'\\'+name[:-4]+'-P.coeff', 'r') as infile:
    transform_matrix[0,2]=float(infile.readline())-256
    transform_matrix[0,0]=float(infile.readline())   
    transform_matrix[0,1]=float(infile.readline())
    transform_matrix[1,2]=float(infile.readline())  
    transform_matrix[1,0]=float(infile.readline())
    transform_matrix[1,1]=float(infile.readline())  
    
ii=1
jj=82


im=(read_one_page_pma(root,name, pageNb=ii))[0]
plt.close('all')
plt.figure(1), plt.imshow(im),plt.colorbar()

im_correct=im-im_bg
im_correct[im_correct<0]=0
im_correct2,threshold=remove_background(im_correct, show=1)
plt.figure(1), plt.subplot(1,1,1), plt.imshow(im_correct2),plt.colorbar()

IM_donor=im_correct2[:,0:int(vdim/2)]
IM_acceptor=im_correct2[:,int(vdim/2):]

array_size=np.shape(IM_acceptor)
imA=cv2.warpPerspective(IM_acceptor.astype(float), transform_matrix,array_size[::-1])

plt.figure(10), 
ax=plt.subplot(1,3,1); ax.imshow(IM_donor), ax.set_title('donor')
ax=plt.subplot(1,3,2); ax.imshow(IM_acceptor),ax.set_title('acceptor, not shifted')
ax=plt.subplot(1,3,3); ax.imshow(imA),ax.set_title('acceptor, shifted')

A=(IM_donor>0)+2*(IM_acceptor>0)
B=(IM_donor>0)+2*(imA>0)
plt.figure(11), plt.subplot(1,2,1), plt.imshow(A), plt.title([(A==1).sum(),(A==2).sum(),(A==3).sum()])
plt.figure(11), plt.subplot(1,2,2), plt.imshow(B), plt.title([(B==1).sum(),(B==2).sum(),(B==3).sum()])
# display the size of non overlapping channels, and overlap

#initiating plots, plus put on right position in screen
fig=plt.figure(2, figsize=(6.5,2.5)) #figsize in inches
fig.canvas.manager.window.move(0,0)

fig=plt.figure(3, figsize=(6.5,2.5)) #figsize in inches
fig.canvas.manager.window.move(0,300)

fig=plt.figure(4, figsize=(6.5,2.5)) #figsize in inches
fig.canvas.manager.window.move(0,600)

fig=plt.figure(5, figsize=(6.5,2.5)) #figsize in inches
fig.canvas.manager.window.move(1000,0)

fig=plt.figure(6, figsize=(6.5,2.5)) #figsize in inches
fig.canvas.manager.window.move(1000,300)

fig=plt.figure(7, figsize=(6.5,2.5)) #figsize in inches
fig.canvas.manager.window.move(1000,600)

## donor intensity ##

##LOAD VALUES FROM PKS pts_number,labels,ptsG=analyze_label.analyze(im_correct2[:,0:int(vdim/2)])
##also dstG
ptsG=[]
dstG2=[]
with open(root+'\\'+name[:-4]+'-PC.pks', 'r') as infile:
    cc=0
    for line in infile:
        cc=cc+1
        if cc%2 != 0: # cc%2!= 0 means odd, so ptsG, cc%2!= 1 means even, so dstG
            A=[float(ii) for ii in line.split()]
            ptsG.append(A[1:3])
        else:
            B=[float(ii) for ii in line.split()]
            dstG2.append(B[1:3])
pts_number=len(ptsG)
ptsG=np.array(ptsG)
dstG2=np.array(dstG2)

tmp = np.genfromtxt(fn_IDL) 
Pidl=tmp.reshape(8,4)[0:4,:]
Qidl=tmp.reshape(8,4)[4:,:]
dstGidl=np.zeros(np.shape(ptsG))
    
tmp = np.genfromtxt(fn_Python) 
Ppyt=tmp.reshape(8,4)[0:4,:]
Qpyt=tmp.reshape(8,4)[4:,:]
dstGpyt=np.zeros(np.shape(ptsG))

## make a loop over jj to watch each selected spot individually             
for jj in range(pts_number):
    print(jj)
    xpix=ptsG[jj][1] 
    ypix=ptsG[jj][0]
    
    xpix_int=int(xpix)
    ypix_int=int(ypix)
    DD=10
    impixD=IM_donor[ (xpix_int-DD) : (xpix_int+DD+1) , (ypix_int-DD) : (ypix_int+DD+1)]
    GG=makeGaussian(2*DD+1, fwhm=3, center=(ypix-ypix_int+DD,xpix-xpix_int+DD))
    multipD=impixD*GG
    donor=np.sum(multipD)
    
    plt.figure(2), 
    plt.subplot(1,3,1)
    plt.imshow(impixD), plt.title('donor')
    plt.subplot(1,3,2)
    plt.imshow(GG), plt.title([jj,xpix_int,ypix_int] )
    plt.subplot(1,3,3)
    plt.imshow(multipD)
    plt.title(sum(sum(multipD)))
    plt.pause(0.05)
    
    
    ## acceptor intensity, without warping ##
    
    impixE=IM_acceptor[ (xpix_int-DD) : (xpix_int+DD+1) , (ypix_int-DD) : (ypix_int+DD+1)]
    GGE=makeGaussian(2*DD+1, fwhm=3, center=(ypix-ypix_int+DD,xpix-xpix_int+DD))
    multipE=impixE*GGE
    acceptor0=np.sum(multipE)
    
    plt.figure(3)
    plt.subplot(1,3,1)
    plt.imshow(impixE), plt.title('acceptor, original')
    plt.subplot(1,3,2)
    plt.imshow(GG), plt.title([jj,xpix_int,ypix_int] )
    plt.subplot(1,3,3)
    plt.imshow(multipE)
    plt.title(sum(sum(multipE)))
    plt.pause(0.05)
    
    ## acceptor intensity, warping the image ##
    
    impixA=imA[ (xpix_int-DD) : (xpix_int+DD+1) , (ypix_int-DD) : (ypix_int+DD+1)]
    
    multip=impixA*GG
    acceptor=np.sum(multip)
    
    plt.figure(4)
    plt.subplot(1,3,1)
    plt.imshow(impixA), plt.title('acceptor, image warped')
    plt.subplot(1,3,2)
    plt.imshow(GG), plt.title([jj,xpix_int,ypix_int] )
    plt.subplot(1,3,3)
    plt.imshow(multip)
    plt.title(sum(sum(multip)))
    plt.pause(0.05)
   
    # approach 3, find transformed coordinates in im_correct (full image)
    xf2=dstG2[jj][1]#approach3,  NOTE: xy are swapped here
    yf2=dstG2[jj][0]#approach3
    xf2_int=int(xf2)#approach3
    yf2_int=int(yf2)#approach3
                    
    impixC=im_correct2[ (xf2_int-DD) : (xf2_int+DD+1) , (yf2_int-DD) : (yf2_int+DD+1)]
    GGC=makeGaussian(2*DD+1, fwhm=3, center=(yf2-yf2_int+DD,xf2-xf2_int+DD))#approach3
    multipC=impixC*GGC#approach3
    acceptorC=np.sum(multipC)
    
    plt.figure(5)
    plt.subplot(1,3,1)
    plt.imshow(impixC), plt.title('warped positions')
    plt.subplot(1,3,2)
    plt.imshow(GGC), plt.title([jj,xf2_int,yf2_int] )
    plt.subplot(1,3,3)
    plt.imshow(multipC)
    plt.title(sum(sum(multipC)))
    plt.pause(0.05)
                
    ## acceptor intensity, warping the position, with the IDL transform matrix  ##
    
    ### compare to extraction IDL

    #for LL in range( len(ptsG)):
    new0=0
    new1=0
    old0=ptsG[jj][0]
    old1=ptsG[jj][1]
    for iii in range(4):
            for jjj in range(4):
                new0=new0+Pidl[iii,jjj]* (old0**iii)*(old1**jjj)
                new1=new1+Qidl[iii,jjj]* (old0**iii)*(old1**jjj)
    dstGidl[jj][0]=new0
    dstGidl[jj][1]=new1
    
    xnew=dstGidl[jj][1]
    ynew=dstGidl[jj][0]+256
    xnew_int=int(xnew)
    ynew_int=int(ynew)
    
    fig=plt.figure(6,figsize=(6.5,2.5))
    fig.canvas.manager.window.move(1000,300)
    if (xnew_int>10 and ynew>10 and xnew<hdim/2-10 and ynew<vdim-10):
        impixG=im_correct2[ (xnew_int-DD) : (xnew_int+DD+1) , (ynew_int-DD) : (ynew_int+DD+1)]
        GGC=makeGaussian(2*DD+1, fwhm=3, center=(ynew-ynew_int+DD,xnew-xnew_int+DD))#approach3
        multipG=impixG*GGC#approach3
        acceptorG=np.sum(multipG)
        
        plt.figure(6)
        plt.subplot(1,3,1)
        plt.imshow(impixG), plt.title('warped positions with IDL P&Q')
        plt.subplot(1,3,2)
        plt.imshow(GG),plt.title([jj,xnew_int,ynew_int] )
        plt.subplot(1,3,3)
        plt.imshow(multipG)
        plt.title(sum(sum(multipG)))
        plt.pause(0.05)
    
    ## acceptor intensity, warping the position, with the Python transform matrix  ##
    ### compare to extraction IDL
    

    #for LL in range( len(ptsG)):
    new0=0
    new1=0
    old0=ptsG[jj][0]
    old1=ptsG[jj][1]
    for iii in range(4):
            for jjj in range(4):
                new0=new0+Ppyt[iii,jjj]* (old0**iii)*(old1**jjj)
                new1=new1+Qpyt[iii,jjj]* (old0**iii)*(old1**jjj)
    dstGpyt[jj][0]=new0
    dstGpyt[jj][1]=new1
    
    xnew=dstGpyt[jj][1]
    ynew=dstGpyt[jj][0]+256
    xnew_int=int(xnew)
    ynew_int=int(ynew)
    
    plt.figure(7)
    plt.cla()
    if (xnew_int>10 and ynew>10 and xnew<hdim/2-10 and ynew<vdim-10):
        impixF=im_correct2[ (xnew_int-DD) : (xnew_int+DD+1) , (ynew_int-DD) : (ynew_int+DD+1)]
        GGF=makeGaussian(2*DD+1, fwhm=3, center=(ynew-ynew_int+DD,xnew-xnew_int+DD))#approach3
        multipF=impixF*GGF#approach3
        acceptorF=np.sum(multipD)
        
        plt.figure(7)
        plt.subplot(1,3,1)
        plt.imshow(impixF), plt.title('warped positions with Python P&Q')
        plt.subplot(1,3,2)
        plt.imshow(GG),plt.title([jj,xnew_int,ynew_int] )
        plt.subplot(1,3,3)
        plt.imshow(multipF)
        plt.title(sum(sum(multipF)))
        plt.pause(0.05)
    
    plt.pause(0.1)
    print(jj)

### from the above conclude that the mapping + alignment is working pretty well. 
### from below: 
########################################################################
#
from traceAnalysisCode_MD import Experiment

mainPath=r'N:\\tnw\\BN\\CMJ\\Shared\\Margreet\\20181203_HJ_training\\#3.10_0.1mg-ml_streptavidin_50pM_HJC_G_movies'

exp1 = Experiment(mainPath)
#

for ii in range(len(exp1.files)):
    print([ii,exp1.files[ii].name])
jjM=82
jj=82
plt.figure(1), plt.subplot(1,1,1)
plt.subplot(4,3,1), plt.plot(donor[jj,:])
plt.subplot(4,3,2), plt.plot(acceptor[jj,:])
#plt.subplot(2,3,3), plt.plot(acceptor[jj,:]/(acceptor[jj,:]+donor[jj,:]))

donor_out3=exp1.files[0].molecules[jjM].intensity[0,:]
acceptor_out3=exp1.files[0].molecules[jjM].intensity[1,:]
plt.subplot(4,3,4), plt.plot(donor_out3)
plt.subplot(4,3,5), plt.plot(acceptor_out3)   
#plt.subplot(2,3,6), plt.plot(acceptor_out/(acceptor_out+donor_out))

donor_out4=exp1.files[4].molecules[jj].intensity[0,:]
acceptor_out4=exp1.files[4].molecules[jj].intensity[1,:]
plt.subplot(4,3,7), plt.plot(donor_out4)
plt.subplot(4,3,8), plt.plot(acceptor_out4)   

donor_out5=exp1.files[5].molecules[jj].intensity[0,:]
acceptor_out5=exp1.files[5].molecules[jj].intensity[1,:]
plt.subplot(4,3,10), plt.plot(donor_out5)
plt.subplot(4,3,11), plt.plot(acceptor_out5)
   
if 0:
        
    plt.figure(15) # from this I can tell acceptor 
    cc1=0
    FF=0
    print(exp1.files[FF].name)
    for jj in range(len(exp1.files[FF].molecules)):
            plt.cla()
            exp1.files[FF].molecules[jj].plot()
            plt.title(jj)
            plt.pause(0.1)
            if np.sum(exp1.files[FF].molecules[jj].intensity)>0:
                cc1=cc1+1
                
                
      
            
    plt.figure(16)
    FF=5
    cc2=0
    for jj in range(len(exp1.files[FF].molecules)):
            plt.cla()
            exp1.files[FF].molecules[jj].plot()
            plt.title(jj)
            plt.pause(0.1)
            if np.sum(exp1.files[FF].molecules[jj].intensity)>0:
             cc2=cc2+1  
        
