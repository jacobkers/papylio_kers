# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 12:09:49 2019

@author: mwdocter
"""

# distance calculation
import numpy as np

def dist_calc(posA, posB,LEN=4):
    if type(posA)==list:
        posA=np.array(posA)
        posB=np.array(posB)
    BB=[]
    dist=np.zeros((len(posA),len(posB)))
    for ii in range(0, len(posA)):
        for jj in range(0, len(posB)):
            dist[ii,jj]=np.sqrt((posA[ii,0]-posB[jj,0])**2+(posA[ii,1]-posB[jj,1])**2)
    for ii in range(0,len(posA)):
        jj=np.where(dist[ii,:]==min(dist[ii,:]))[0][0] # use only 0th element
        if dist[ii,jj]<LEN and ii==np.where(dist[:,jj]==min(dist[:,jj]))[0][0]: 
            BB.append([ii,jj])
    BB=np.array(BB)
    
    return BB,dist