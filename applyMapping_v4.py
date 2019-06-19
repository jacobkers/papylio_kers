# -*- coding: utf-8 -*-
"""
Created on Jun 2019
@author: carlos
meant for mapping, but now with class implementation
"""
import os
from ImageCollection import ImageCollection
from Mapping import Mapping
import numpy as np
import matplotlib.pyplot as plt
import time

tetra_fn=os.path.normpath('N:/tnw/BN/CMJ/Shared/Margreet/20181203_HJ_training/mapping/rough_ave.tif')

#tetra_fn='/home/carlos/PycharmProjects/margreet_code_review/rough_ave.tif'
#image_fn ='/home/carlos/PycharmProjects/margreet_code_review/hel4.pma'

tmp=Mapping(tetra_fn) # if you want to test the mapping seperately, it is included in ImageCollection

if 1:
    root=os.path.normpath("N:/tnw/BN/CMJ/Shared/Margreet/20181203_HJ_training/#3.10_0.1mg-ml_streptavidin_50pM_HJC_G_movies")
    name='hel4.pma'
    root=os.path.normpath(r"N:\tnw\BN\CMJ\Shared\Margreet\20181203_HJ_training\#1.10_0.1mg-ml_streptavidin_50pM_HJA_G_movies")
    name='hel2.pma'                          
else:
    root=os.path.normpath("N:/tnw/BN/CMJ/Shared/Margreet/181218 - First single-molecule sample (GattaQuant)/RawData") 
    name='Spooled files.sifx'
 
image_fn=os.path.join(root,name)
imc = ImageCollection(tetra_fn, image_fn)  # **kwargs, 'pks_fn'= ; choice_channel='d'
# better: give image directory, and search for specific name or extension
img = imc.get_image(1)

t1=time.time()
donor,acceptor=imc.get_all_traces()
t2=time.time()-t1;
print(t2)
#example codes to return data from mapping, image collection
    #tmp._tf2_matrix
    #imc.mapping._tf2_matrix
    #imc.pts_number # is too 19 for the pma
    #import matplotlib.pyplot as plt; plt.imshow(imc.im_mean20_correct)


         