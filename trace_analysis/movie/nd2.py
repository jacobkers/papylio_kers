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
import xarray as xr
from nd2reader import ND2Reader

from trace_analysis.movie.movie import Movie, Illumination


class ND2Movie(Movie):
    extensions = ['.nd2']

    def __init__(self, arg, *args, **kwargs):
        super().__init__(arg, *args, **kwargs)
        # super().__init__(arg)

        self.writepath = self.filepath.parent
        self.name = self.filepath.with_suffix('').name

        if 'fov' in self.name:
            token_position = self.name.find('_fov')
            self.name = self.name[:token_position]

        # setting for multi fov measurement
        # TODO: Move fov info / fov selection here if possible
        self.fov_info = None
        if 'fov_info' in kwargs:
            self.fov_info = kwargs['fov_info']  # fov=Field of View

        self.threshold = {'view': (0, 200),
                          'point-selection': (45, 25)
                          }

        # We should probably put this in the configuration file
        # SHK: self.rot90 should be set before reading the header.
        self.rot90 = 1

        self.f_obj = ND2Reader(str(self.filepath))
        if 'c' in self.f_obj.axes:
            self.f_obj.iter_axes = 'tc'
        else:
            self.f_obj.iter_axes = 't'
        # self.read_header()

        # self.time = self.time[self.fov_info['first_frame_of_each_fov'][self.fov_info['fov_chosen']]:(self.fov_info['last_frame_of_each_fov'][self.fov_info['fov_chosen']]+1)]
        # self.illumination = self.illumination[self.fov_info['first_frame_of_each_fov'][self.fov_info['fov_chosen']]:(self.fov_info['last_frame_of_each_fov'][self.fov_info['fov_chosen']]+1)]
        self.create_frame_info()  # Possibly move to Movie later on

        # self._initialized = True

    def _read_header(self):
        with ND2Reader(str(self.filepath)) as images:
            self.width = images.metadata['width']
            self.height = images.metadata['height']
            if 'c' in images.axes:
                images.iter_axes = 'tc'
            else:
                images.iter_axes = 't'
            # self.number_of_fields_of_view = len(images.metadata["experiment"]["loops"])  # number of fov is now available from self.fov_info
            self.number_of_frames = len(images)

            self.illuminations = [Illumination(self, name) for name in images.metadata["channels"]]

            if self.fov_info:
                self.frame_offset = self.fov_info['first_frame_of_each_fov'][self.fov_info['fov_chosen']] * self.number_of_illuminations
                frame_end = (self.fov_info['last_frame_of_each_fov'][self.fov_info['fov_chosen']]+1) * self.number_of_illuminations
                self.number_of_frames = frame_end - self.frame_offset
            else:
                self.frame_offset = 0
                frame_end = self.number_of_frames

            self.illumination_arrangement = np.arange(len(self.illuminations))

            if self.fov_info:
                self.time = xr.DataArray(np.repeat(
                    images.timesteps[self.fov_info['first_frame_of_each_fov'][self.fov_info['fov_chosen']]:(self.fov_info['last_frame_of_each_fov'][self.fov_info['fov_chosen']]+1)],
                    self.number_of_illuminations)/1000, dims='frame', coords={}, attrs={'units': 's'})
            else:
                self.time = xr.DataArray(np.repeat(images.timesteps, self.number_of_illuminations)/1000, dims='frame',
                                         coords={}, attrs={'units': 's'})

            # self.exp_time = images.metadata['experiment']['loops'][0]['sampling_interval']
            # self.exp_time_start=images.metadata['experiment']['loops'][0]['start']
            # self.exp_time_duration=images.metadata['experiment']['loops'][0]['duration']
            # self.pixelmicron=images.metadata['experiment']['pixel_microns']

    def _read_frame(self, frame_number):
        if frame_number < self.number_of_frames:
            im = self.f_obj[frame_number + self.frame_offset]
        else:
            im = self.f_obj[self.number_of_frames + self.frame_offset - 1]
            print(f'pageNb out of range. The last frame (fr#{self.number_of_frames + self.frame_offset - 1}) is loaded instead')
        # note: im is a Frame, which is pims.frame.Frame, a np. array with additional frame number and metadata
        return im


if __name__ == "__main__":
    print('test')
