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
    extensions = ['.nd2']

    def __init__(self, arg, *args, **kwargs):
        super().__init__(arg, *args, **kwargs)

        self.writepath = self.filepath.parent
        self.name = self.filepath.with_suffix('').name

        self.threshold = {'view': (0, 200),
                          'point-selection': (45, 25)
                          }

        # We should probably put this in the configuration file
        # SHK: self.rot90 should be set before reading the header.
        self.rot90 = 1

        self.f_obj = ND2Reader(self.filepath)
        self.read_header()
        self.create_frame_info()  # Possibly move to Movie later on


    def _read_header(self):
        with ND2Reader(self.filepath) as images:
            self.width = images.metadata['width']
            self.height = images.metadata['height']
            # self.number_of_frames = images.metadata['num_frames']
            self.number_of_frames = len(images)
            # self.bitdepth = 16
            self.exp_time = images.metadata['experiment']['loops'][0]['sampling_interval']
           # self.exp_time_start=images.metadata['experiment']['loops'][0]['start']
           # self.exp_time_duration=images.metadata['experiment']['loops'][0]['duration']
           # self.pixelmicron=images.metadata['experiment']['pixel_microns']

    def _read_frame(self, frame_number):
        #   images.iter_axes = 'z' #maybe later on you want to iterate over other axis of 6D measurement
        if self.number_of_frames == 1:
            im = self.f_obj[0]
        elif (self.number_of_frames - 1) >= frame_number:
            im = self.f_obj[frame_number]
        else:
            im = self.f_obj[self.number_of_frames - 1]
            print('pageNb out of range, printed image {0} instead'.format(self.number_of_frames))
        # note: im is a Frame, which is pims.frame.Frame, a np. array with additiona frame umber and metadata

        return im


if __name__ == "__main__":
    print('test')
