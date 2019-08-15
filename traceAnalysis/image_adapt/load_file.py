# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 10:45:41 2019

@author: mwdocter

main function = read_one_page, which then refers to reading a tif/pma/sifx file
"""
import os
import re
import tifffile
import numpy as np 
import matplotlib.pyplot as plt
import time
#from sifreader.sifreader import SIFFile
from image_adapt.sifreaderA import SIFFile

from autopick.do_before import clear_all
clear_all()

def read_one_page(image_fn,pageNb,A=None,ii=None): # distributes to specific readers (sif/tif/pma)
    root, name = os.path.split(image_fn)
    for string2compare in ['.pma$','.tif$','.sifx$','.sif$']:
        AA=re.search(string2compare,name)
        im0=0
        if AA!=None: 
            #print(name, string2compare)
            if string2compare=='.pma$':
                im0=read_one_page_pma(root,name,pageNb,A)
            elif string2compare=='.tif$':
                im0=read_one_page_tif(root,name,pageNb)
            elif string2compare=='.sifx$':
                im0=read_one_page_sifx(root,name,pageNb,A,ii)
#            elif A!=None and string2compare=='.sif$':
#                im0,hdim,vdim,nImages=read_one_page_sif(root,name,pageNb)
            break
    return im0

def read_header(image_fn): # distributes to specific readers (sif/tif/pma)
    root, name = os.path.split(image_fn)
    for string2compare in ['.pma$','.tif$','.sifx$','.sif$']:
        AA=re.search(string2compare,name)
        if AA!=None: 
            #print(name, string2compare)
            if string2compare=='.pma$':
                hdim,vdim,nImages,A=read_header_pma(root,name)
            elif string2compare=='.tif$':
                hdim,vdim,nImages=read_header_tif(root,name)
                A=None
            elif string2compare=='.sifx$':
                hdim,vdim,nImages,A=read_header_sifx(root,name)
#            elif A!=None and string2compare=='.sif$':
#                im0,hdim,vdim,nImages=read_one_page_sif(root,name,pageNb)
            break
    return hdim,vdim,nImages,A
    
def read_header_pma(root,name):
    statinfo = os.stat(root+'\\'+name)       
    #determine 8 bits or 16 bits
    A=re.search('_16.pma$',name)
    with open(root+'\\'+name, 'rb') as fid:
        hdim = np.fromfile(fid, np.int16,count=1)
        vdim=  np.fromfile(fid, np.int16,count=1)
        hdim=int(hdim[0])
        vdim=int(vdim[0])
    
    nImages=int((statinfo.st_size-4)/(hdim*vdim))
        
    return hdim,vdim,nImages,A

def read_one_page_pma(root,name, pageNb=0,A=None):
     
    with open(root+'\\'+name, 'rb') as fid:
        hdim = np.fromfile(fid, np.int16,count=1)
        vdim=  np.fromfile(fid, np.int16,count=1)
        hdim=int(hdim[0])
        vdim=int(vdim[0])
    
        if A==None: #8 bits
#        with open(root+'\\'+name, 'rb') as fid: #did the offset reset?    # is already open
            # for image pageNb, 4 for skipping header, plus certain amount of images to read image pageNb
            fid.seek(4+ (pageNb*(hdim*vdim)), os.SEEK_SET)
            im=np.reshape(np.fromfile(fid,np.int8,count=hdim*vdim),(hdim,vdim))
        else:
#        with open(root+'\\'+name, 'rb') as fid: #did the offset reset?  #is already open
            fid.seek(4+ 2*pageNb*(hdim*vdim), os.SEEK_SET)
            msb=np.reshape(np.fromfile(fid,np.int8,count=(hdim*vdim)),(hdim,vdim))
            lsb=np.reshape(np.fromfile(fid,np.int8,count=(hdim*vdim)),(hdim,vdim))
#            msb = np.core.records.fromfile(fid, 'int8', offset=4+ 2*pageNb*(hdim*vdim), shape=(hdim,vdim)) # for first image
#            lsb = np.core.records.fromfile(fid, 'int8', offset=4+ (1+2*pageNb)*(hdim*vdim), shape=(hdim,vdim)) # for first image
            im=256*msb+lsb;
    return im # still need to convert im

def read_header_tif(root,name):
    with tifffile.TiffFile(root+'\\'+name) as tif:
        tifpage=tif.pages
    nImages=(len(tifpage))   
    im = tifpage[0].asarray()
    hdim,vdim=np.shape(im)
    return hdim,vdim,nImages
    
def read_one_page_tif(root,name, pageNb=0):
    t = time.time()

    with tifffile.TiffFile(root+'\\'+name) as tif:
#         tif_tags = {}
#         for tag in tif.pages[0].tags.values():
#            name, value = tag.name, tag.value
#            tif_tags[name] = value
#         hdim=tif_tags['ImageWidth']
#         vdim=tif_tags['ImageLength']   
        ### hdim,vdim=tif.pages[0].shape
        
         tifpage=tif.pages
         nImages=(len(tifpage))
  
         if nImages==1:
             #return -1,0,0,0
             im=tifpage[0].asarray()
         elif (nImages-1)>=pageNb:
              im = tifpage[pageNb].asarray()
         else:
              im = tifpage[nImages-1].asarray()
              print('pageNb out of range, printed image {0} instead'.format(nImages))
         hdim,vdim=np.shape(im)
    elapsed = time.time() - t
    print(elapsed),
    return im

    
def read_one_page_sif(root,name,pageNb=0):
    file = SIFFile(root+'\\'+name)
    im = file.read_block(0) 
    hdim,vdim=np.shape(im)
    
    print('sif not working yet')
    return im, 

def read_header_sifx(root,name):
 #   filelist = [x for x in os.listdir(root) if x.endswith("spool.dat")]
    A=SIFFile(root+'\\Spooled files.sifx')
    hh=A.height
    ww=A.width
    nImages=A.stacksize
    LL=len(A.filelist)
    new_filelist=[]
    
    #making new_list with correct numerical image name
    for ii in range(LL):
        filename=A.filelist[ii]
        new_filelist.append(filename[9::-1])
        
    #merge the old and new filenames into 2D list    
    def merge(list1, list2): 
        merged_list = [(list1[i], list2[i]) for i in range(0, len(list1))] 
        return merged_list 
    AA=merge(A.filelist,new_filelist)
    
    #sort the 2d list to the second element=new filename
    BB=AA.copy()
    def takeSecond(elem):
        return elem[1]
    BB.sort(key=takeSecond)
    
    for ii in range(LL):
        #do something with ....
        #print([ii,BB[ii][0],BB[ii][1]])
        A.filelist[ii]=BB[ii][0]
  
    return hh,ww, nImages,A
   
def read_one_page_sifx(root,name, pageNb,A,ii=0):
    if (A.xbin == 2) and (A.ybin == 2):
         count=A.height*A.width*4
         #name should follow from A.filelist
         with open(root+'\\'+A.filelist[pageNb], 'rb') as fid:
             raw=np.uint16(np.fromfile(fid,np.uint8,count))
         ALL = raw[0::4]+raw[1::4]*256
        
    else:
         count=A.height*A.width*3//2
         #name should follow from A.filelist
         with open(root+'\\'+A.filelist[pageNb], 'rb') as fid:
              raw=np.uint16(np.fromfile(fid,np.uint8,count))
         
         
         ii = np.array(range(int(A.width/2)*int(A.height/2)))
            
         #print([A.height,A.width,A.stacksize,np.shape(raw)])        
         AA=raw[ii*6+0]*16 + (raw[ii*6+1]%16)
         BB=raw[ii*6+2]*16 + (raw[ii*6+1]//16)
         CC=raw[ii*6+3]*16 + (raw[ii*6+4]%16)
         DD=raw[ii*6+5]*16 + (raw[ii*6+4]//16) 
              
         ALL=np.uint16(np.zeros(A.height*A.width))
         ALL[0::4] = AA
         ALL[1::4] = BB
         ALL[2::4] = CC
         ALL[3::4] = DD
              
    im=np.reshape(ALL,(A.height, A.width))
    im=np.rot90(im)     
    if 0: # for testing match real data
        plt.imshow(im)
        tifffile.imwrite(root+"Python image.tif" , im ,  photometric='minisblack')
    
    return im   
 
    
# MATLAB code    
#    filename='0000000000spool.dat'
#    % filedir='D:\data\personal data\20190110 testing matlab image conversion\spool\';
#    
#    file = fopen([filedir,filename],'rb')
#    [data] = fread(file,'uint8=>double',0);
#    fclose(file);
#    %
#    dd=0; %images for dd 0:3:39 look reasonable, but not identical
#    ii=1:(1024*1024); %floor(numel(data)/6)=1048582, 1024*1024=1048576
#    AA=data(dd+(ii-1)*6+1)*(2^4) +rem(dd+data((ii-1)*6+2),2^4);
#    BB=data(dd+(ii-1)*6+3)*(2^4) +floor(dd+data((ii-1)*6+2)/(2^4));
#    CC=data(dd+(ii-1)*6+4)*(2^4) +rem(dd+data((ii-1)*6+5),2^4);
#    DD=data(dd+(ii-1)*6+6)*(2^4) +floor(dd+data((ii-1)*6+5)/(2^4));
#    
#    ALL=zeros(2048);
#    ALL(1:4:end)=AA(1:end);
#    ALL(2:4:end)=BB(1:end);
#    ALL(3:4:end)=CC(1:end);
#    ALL(4:4:end)=DD(1:end);
#        
#    ALLm=ALL(end:-1:1,:);    
    
# BELOW A TEST IS DESCRIBED
#if 0: #test pma
#    mainPath='H:\\projects\\research practicum\\single molecule fluorescence\\Matlab\\HJA-data from Ivo'
#elif 0: #test multipage tif
#    mainPath='N:\\tnw\\BN\\CMJ\\Shared\\Margreet\\180920 super resolution folder\\190326 data superres C1 D1-210'
#elif 1:
#    #mainPath='N:\\tnw\\BN\\CMJ\\Shared\\Margreet\\older\\170111 cmos camera\\data\\170228 cy37'
#    mainPath='E:\\CMJ trace analysis\\test data'
#elif 0: 
#    mainPath='N:\\tnw\\BN\\CMJ\\Shared\\Margreet\\181218 - First single-molecule sample (GattaQuant)\\RawData'
##read in pma data
#for root, dirs, fileNames in os.walk(mainPath, topdown=False):
#    for name in fileNames:
#        # detect whether you would like to read a pma file, sif(x), or multipage tif
#        im0,hdim,vdim,nImages=read_one_page(root,name,0)
#        im1=read_one_page(root,name,1)[0]
#        im100=read_one_page(root,name,100)[0]
#        if type(im0)!=int:
#            plt.figure(1)
#            plt.subplot(1,3,1)
#            plt.imshow(im0)
#            plt.subplot(1,3,2)
#            plt.imshow(im1)
#            plt.subplot(1,3,3)
#            plt.imshow(im100)
#            plt.show()
#        
#            break #for testing do only the first one in the file directory
#        