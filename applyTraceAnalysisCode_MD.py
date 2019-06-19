# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 10:50:59 2018

@author: ivoseverins, adapted by margreet to play with my own code

for example: 
    make histogram of fret values
    read in pks value and make scatter plot
    plot individual traces

"""
# this code is used to analyse the difference between data analyzed through IDL
# and data analyzed through Python.
# currently the FRET signal is too low of the Python part. Why? 
# - too much background subtraction

from traceAnalysisCode_MD import Experiment
import os
import matplotlib.pyplot as plt #Provides a MATLAB-like plotting framework
import numpy as np
#from pick_spots_akaze import mapping

#mainPath = r'D:\ivoseverins\SURFdrive\Promotie\Code\Python\traceAnalysis\twoColourExampleData\HJ A'
mainPath=r'E:\CMJ trace analysis\test data'
os.chdir(os.path.dirname(os.path.abspath(__file__)))

#mainPath = './twoColourExampleData/HJ A'

#T,fig=mapping()

plt.close("all")

exp1 = Experiment(r'E:\CMJ trace analysis\test data\voorbeeld data door Python')
exp2 = Experiment(r'E:\CMJ trace analysis\test data\voorbeeld data door IDL')

plt.figure(1)
exp1.histogram()

plt.figure(2)
exp2.histogram()


# compare .pks files
filename="E:\\CMJ trace analysis\\test data\\voorbeeld data door IDL\\hel6.pks";
pks = np.genfromtxt(filename) 

filename2='E:\\CMJ trace analysis\\test data\\voorbeeld data door Python\\hel6-P.pks';
pks2=np.genfromtxt(filename2) 

filename3='E:\\CMJ trace analysis\\test data\\voorbeeld data door Python\\hel6-PC.pks';
pks3=np.genfromtxt(filename3) 

#im1=np.zeros((512,512));
#for ii in range (len(pks)):
#    xp=int(np.min([np.max([1,np.round(pks[ii,2])]),512]))
#    yp=int(np.min([np.max([1,np.round(pks[ii,3])]),512]))
#    im1[xp,yp]=1;
#
#im2=np.zeros((512,512));
#for ii in range (len(pks2)):
#    xp=int(np.min([np.max([1,np.round(pks2[ii,2])]),512]))
#    yp=int(np.min([np.max([1,np.round(pks2[ii,3])]),512]))
#    im2[xp,yp]=1;
#
#im3=np.zeros((512,512));
#for ii in range (len(pks3)):
#    xp=int(np.min([np.max([1,np.round(pks3[ii,2])]),512]))
#    yp=int(np.min([np.max([1,np.round(pks3[ii,3])]),512]))
#    im3[xp,yp]=1;
#
#plt.figure(11), plt.imshow(im1)
#plt.figure(12), plt.imshow(im2)
#plt.figure(13), plt.imshow(im3)

plt.figure(11), plt.scatter(pks[:,1],pks[:,2])
#plt.figure(12), plt.scatter(pks2[:,1],pks2[:,2])
plt.figure(13), plt.scatter(pks3[:,1],pks3[:,2])
plt.figure(14), plt.scatter(pks[:,1],pks[:,2],c='r', marker='o'), 
plt.scatter(pks3[:,1],pks3[:,2],c='b', marker='x')

# manually matched points
plt.figure(), exp2.files[1].molecules[125].plot() # pks 251&252
plt.figure(), exp1.files[2].molecules[794].plot() # pks-PC 1589&1590

plt.figure(5) # from this I can tell acceptor 
cc1=0
FF=0
for ii in range(len(exp1.files[FF].molecules)):
    if np.sum(exp1.files[FF].molecules[ii].intensity)>0:
        plt.cla()
        exp1.files[FF].molecules[ii].plot()
        plt.title(ii)
        cc1=cc1+1
        plt.pause(0.1)
        
plt.figure(6)
cc2=0
for ii in range(len(exp2.files[1].molecules)):
    if np.sum(exp2.files[1].molecules[ii].intensity)>0:
        plt.cla()
        exp2.files[1].molecules[ii].plot()
        plt.title(ii)
        cc2=cc2+1  
        plt.pause(0.1)
    