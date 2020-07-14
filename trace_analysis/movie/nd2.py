# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 15:50:57 2020

@author: mwdoc
https://github.com/soft-matter/pims_nd2
read_header and read_frame adapted, def __init__ unchanged
"""


from pathlib import Path
import os, sys
   

import time
import numpy as np
import matplotlib.pyplot as plt
from nd2reader import ND2Reader

from trace_analysis.movie.movie import Movie


class ND2Movie(Movie):
    def __init__(self, arg, *args, **kwargs):
        super().__init__(arg, *args, **kwargs)
        
        self.writepath = self.filepath.parent
        self.name = self.filepath.with_suffix('').name
        

        self.threshold = {  'view':             (0,200),
                            'point-selection':  (45,25)
                            }


        self.read_header()

#
    def read_header(self):
        with ND2Reader(self.filepath) as images:
            self.width = images.metadata['width']
            self.height = images.metadata['height']
            self.number_of_frames = images.metadata['num_frames']
            self.bitdepth = 16
            self.exp_time=images.metadata['experiment']['loops'][0]['sampling_interval']
#            self.exp_time_start=images.metadata['experiment']['loops'][0]['start']
#            self.exp_time_duration=images.metadata['experiment']['loops'][0]['duration']
#            self.pixelmicron=images.metadata['experiment']['pixel_microns']
            ##Note: nothing is said 
            


    def read_frame(self, frame_number):
        # t = time.time()
        with ND2Reader(self.filepath) as images:
            images.iter_axes = 'z' #maybe later on you want to iterate over other axis of 6D measurement
            if self.number_of_frames == 1:
                # return -1,0,0,0
                im =images[0]
            elif (self.number_of_frames - 1) >= frame_number:
                im = images[frame_number]
            else:
                im = images[self.number_of_frames - 1]
                print('pageNb out of range, printed image {0} instead'.format(self.number_of_frames))
        # note: im is a Frame, which is pims.frame.Frame, a np. array with additiona frame umber and metadata
        
        
        im1=np.rot90(im) 
        # note this is required for Sung Hyun's setup, rotation can also be done in NIS Elements.
        # to be double checked whether left if Cy3 or Cy5; otherwise rotate + fliplr or flipud
     #   plt.figure(109),  plt.subplot(1,2,1),  plt.imshow(im1)
        print('!@#$%^&*(*&^%$#@#$%^&*(*&^%$##$%^&*()*&^%$#$%^&*(*&^%$#')
        return im1


if __name__ == "__main__":
    print('test')
    
    
## from the examples:
#    
#from pims import ND2_Reader
#frames = ND2_Reader('some_movie.nd2')
#frames[82]  # display frame 82
#frames.close()
#
#from pims import ND2_Reader
#with ND2_Reader('cluster.nd2') as frames:
#	frames.iter_axes = 't'  # 't' is the default already
#	frames.bundle_axes = 'zyx'  # when 'z' is available, this will be default
#	frames.default_coords['c'] = 1  # 0 is the default setting
#	for frame in frames[:3]:
#		# do something with 3D frames in channel 1
#        
#frames.metadata['mpp']  # calibration in microns per pixel
#frames[0].metadata['t_ms']  # time of frame in milliseconds


    