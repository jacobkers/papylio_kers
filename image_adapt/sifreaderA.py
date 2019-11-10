# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 11:18:26 2019

@author: https://github.com/lightingghost/sifreader/blob/master/sifreader/sifreader.py, mwdocter

Andor does not provide a sif reader for Python. Many details can be found by trial and error in sifx
in the metadata file you can actually find the pixel encoding
"""

import os
import time
import numpy as np

#example filename=r'N:\tnw\BN\CMJ\Shared\Margreet\181218 - First single-molecule sample (GattaQuant)\RawData\Spooled files.sifx';
# A=SIFFile(filename)
# A.__repr__
# A.width,A.height,A.stacksize returns  (2048, 2048, 5000)
class SIFFile(object):
    """
    A class that reads the contents and metadata of an Andor .sif file. Compatible with images as well as spectra.
    Exports data as numpy array or xarray.DataArray.

    Example: SIFFile('my_spectrum.sif').read_all()

    In addition to the raw data, SIFFile objects provide a number of meta data variables:
    :ivar x_axis: the horizontal axis (can be pixel numbers or wavelength in nm)
    :ivar original_filename: the original file name of the .sif file
    :ivar date: the date the file was recorded
    :ivar model: camera model
    :ivar temperature: sensor temperature in degrees Celsius
    :ivar exposuretime: exposure time in seconds
    :ivar cycletime: cycle time in seconds
    :ivar accumulations: number of accumulations
    :ivar readout: pixel readout rate in MHz
    :ivar xres: horizontal resolution
    :ivar yres: vertical resolution
    :ivar width: image width
    :ivar height: image height
    :ivar xbin: horizontal binning
    :ivar ybin: vertical binning
    :ivar gain: EM gain level
    :ivar vertical_shift_speed: vertical shift speed
    :ivar pre_amp_gain: pre-amplifier gain
    :ivar stacksize: number of frames
    :ivar filesize: size of the file in bytes
    :ivar m_offset: offset in the .sif file to the actual data
    """

    def __init__(self, filepath):
        self.filepath = filepath
        self._read_header(filepath)
        self.find_filelist(filepath) 

    def find_filelist(self,filepath):
        root=self.filepath.replace('Spooled files.sifx','')
        files = os.listdir(root)
        files_spool = [i for i in files if i.endswith('spool.dat')]
        self.filelist=files_spool
       
    def __repr__(self):
        info = (('Original Filename', self.original_filename),
                ('Date', self.date),
                ('Camera Model', self.model),
                ('Temperature (deg.C)', '{:f}'.format(self.temperature)),
                ('Exposure Time', '{:f}'.format(self.exposuretime)),
                ('Cycle Time', '{:f}'.format(self.cycletime)),
                ('Number of accumulations', '{:d}'.format(self.accumulations)),
                ('Pixel Readout Rate (MHz)', '{:f}'.format(self.readout)),
                ("Horizontal Camera Resolution", '{:d}'.format(self.xres)),
                ("Vertical Camera Resolution", '{:d}'.format(self.yres)),
                ("Image width", '{:d}'.format(self.width)),
                ("Image Height", '{:d}'.format(self.height)),
                ("Horizontal Binning", '{:d}'.format(self.xbin)),
                ("Vertical Binning", '{:d}'.format(self.ybin)),
                ("EM Gain level", '{:f}'.format(self.gain)),
                ("Vertical Shift Speed", '{:f}'.format(self.vertical_shift_speed)),
                ("Pre-Amplifier Gain", '{:f}'.format(self.pre_amp_gain)),
                ("Stacksize", '{:d}'.format(self.stacksize)),
                ("Filesize", '{:d}'.format(self.filesize)),
                ("Offset to Image Data", '{:f}'.format(self.m_offset)))
        desc_len = max([len(d) for d in list(zip(*info))[0]]) + 3
        res = ''
        for description, value in info:
            res += ('{:' + str(desc_len) + '}{}\n').format(description + ': ', value)

        res = super().__repr__() + '\n' + res
        return res

    def _read_header(self, filepath):
        f = open(filepath, 'rb')
        headerlen = 32
    #    spool = 0
        for ii in range(50):#headerlen + spool):
            line = f.readline().strip()
         #   print(ii,line)
            if ii == 0:
                if line != b'Andor Technology Multi-Channel File':
                    f.close()
                    raise Exception('{} is not an Andor SIF file'.format(filepath))
            # elif ii==1: # line=b'65538 1' , no clue what this means. 2048**2/64=65536
            elif ii == 2:
                tokens = line.split()
                self.temperature = float(tokens[5])
                self.date = time.strftime('%c', time.localtime(float(tokens[4])))
                self.exposuretime = float(tokens[12])
                self.cycletime = float(tokens[13])
                self.accumulations = int(tokens[15])
                self.readout = 1 / float(tokens[18]) / 1e6
                self.gain = float(tokens[21])
                self.vertical_shift_speed = float(tokens[41])
                self.pre_amp_gain = float(tokens[43])
            elif ii == 3:
                self.model = line.decode('utf-8')
            elif ii==4: #nImages is wrong, for test file it should be 5000, Python returns 40
                self.width,self.height,_=[int(ii) for ii in line.decode('utf-8').split()]
            elif ii == 5:
                self.original_filename = line.decode('utf-8') # not so useful if you copy to a different computer
            #elif ii==6: # b'65538 2048'
#            elif ii == 7: # lots of characters, tokens[0] is not Spooled, maybe wrong line?
#                tokens = line.split()
#                if len(tokens) >= 1 and tokens[0] == 'Spooled':
#                    spool = 1
            #elif ii==8: #same junk as ii==7
            #elif i == 9: #similar junk characters as ii==7
            #    wavelength_info = line.split() # example [b'\x00\x00\x00\x00\x00\x00\x00\xfb\x97\x1e2\x05\x00\x00\x00\x00\x00\x00\x00']
#                self.center_wavelength = float(wavelength_info[3])
#                self.grating = float(wavelength_info[6])
#                self.grating_blaze = float(wavelength_info[7])
            #elif 13: # (b'65538 \x01 \x02 \x03 \x00 0 0',)
            #elif ii==14: # (b'65540 0 0 500 0 0 1200 1200',)
            #elif ii==17: # (b'0 SR303i',) 
#            elif i == 19: #(b'0 10',)
#                self.wavelength_coefficients = [float(num) for num in line.split()][::-1]
#                self.wavelength_coefficients = 0
            #elif 21: # (b'65537 1 500 200',)
            #elif 7 < ii < headerlen - 12:
###            elif len(line) == 37 and line[0:6] == b'65539 ':#len(line) == 17
###                   # and line[7] == b'x01' and line[8] == b'x20' \
###                   # and line[9] == b'x00':
###                    headerlen = headerlen + 12
#                   
###            elif ii == 43: #headerlen - 2:
            elif line[:12] == b'Pixel number' and len(line)>14:
                line = line[12:]
                tokens = line.split()
                if len(tokens) < 6:
                    raise Exception('Not able to read stacksize.')
                self.yres = int(tokens[2])
                self.xres = int(tokens[3])
                self.stacksize = int(tokens[5])
###            elif ii == 44: #headerlen - 1:  ( b'65538 1 2048 2048 1 1 1 0')
                #continue with next line
                line = f.readline().strip()
    #            print(ii),print(line)
                tokens = line.decode('utf-8').split()
#                print(tokens)
                if len(tokens) < 7:
                   raise Exception("Not able to read Image dimensions.")
                self.left = int(tokens[1])
                self.top = int(tokens[2])
                self.right = int(tokens[3])
                self.bottom = int(tokens[4])
                self.xbin = int(tokens[5])
                self.ybin = int(tokens[6])
#                 self.left=0
#                 self.right=self.left+self.width
#                 self.bottom=0
#                 self.top=self.bottom+self.height
#                 self.xbin=1
#                 self.ybin=1
                break
     
        f.close()

#        width = self.right - self.left + 1
        width=self.width
        mod = width % self.xbin
        self.width = int((width - mod) / self.ybin)
#        height = self.top - self.bottom + 1
        height=self.height
        mod = height % self.ybin
        self.height = int((height - mod) / self.xbin)

        self.filesize = os.path.getsize(filepath)
        self.datasize = self.width * self.height * 4 * self.stacksize
        self.m_offset = self.filesize - self.datasize - 8

#        self.x_axis = np.polyval(self.wavelength_coefficients, np.arange(self.left, self.right + 1))
# from here, the rest is most likely for .sif and not .sifx
#    def read_block(self, num=0):
#        """
#        Returns a specific block (i.e. frame) in the .sif file as a numpy array.
#        :param num: block number
#        :return: a numpy array with shape (y, x)
#        """
#        f = open(self.filepath, 'rb')
#        f.seek(self.m_offset + num * self.width * self.height * 4)
#        block = f.read(self.width * self.height * 4)
#        data = np.fromstring(block, dtype=np.float32)
#        f.close()
#        return data.reshape(self.height, self.width)
#
#    def read_all(self):
#        """
#        Returns all blocks (i.e. frames) in the .sif file as a numpy array.
#        :return: a numpy array with shape (blocks, y, x)
#        """
#        f = open(self.filepath, 'rb')
#        f.seek(self.m_offset)
#        block = f.read(self.width * self.height * self.stacksize * 4)
#        data = np.fromstring(block, dtype=np.float32)
#        f.close()
#        return data.reshape(self.stacksize, self.height, self.width)
#
#    def as_xarray(self, x_axis_quantity='wavelength'):
#        """
#        Returns an xarray.DataArray object containing all frames, all metadata and all coordinate axis values. If the
#        file contains a spectrum, the wavelength axis can be converted to wavenumbers or photon energy.  This method
#        requires the xarray package to be installed.
#        :param x_axis_quantity: Only relevant for spectra. Can be either 'wavelength' (default), 'wavenumber' or \
#        'photon energy'.
#        :return: An xarray.DataArray containing the entire contents of the .sif file.
#        """
#        try:
#            import xarray as xr
#            data = self.read_all()
#            y_axis = np.arange(data.shape[1])
#            # determine if it's an image or a spectrum: check if the x_axis spacing is always one
#            if (np.abs((np.diff(self.x_axis) - 1)) < 1e-5).all():
#                # it's an image
#                x_axis = self.x_axis.astype(np.uint16)
#                x_axis_quantity = 'x'
#                x_unit = 'px'
#                x_name = 'x'
#            else:
#                # it's a spectrum
#                if x_axis_quantity == 'wavelength':
#                    x_axis = self.x_axis
#                    x_unit = 'nm'
#                    x_name = 'Wavelength'
#                elif x_axis_quantity == 'wavenumber':
#                    x_axis = (1e7 / self.x_axis)[::-1]
#                    data = np.flip(data, 2)
#                    x_unit = 'cm^-1'
#                    x_name = 'Wavenumber'
#                elif x_axis_quantity == 'photon_energy':
#                    x_axis = (1239.84 / self.x_axis)[::-1]
#                    data = np.flip(data, 2)
#                    x_unit = 'eV'
#                    x_name = 'Photon energy'
#                else:
#                    raise RuntimeError('X-axis quantity "{}" not recognized!'.format(x_axis_quantity))
#
#            if data.shape[0] == 1:
#                # Only one frame
#                data = np.transpose(data[0])
#                data_array = xr.DataArray(data, coords=[(x_axis_quantity, x_axis), ('y', y_axis)],
#                                    name='intensity')
#            else:
#                # multiple frames
#                frame_axis = np.arange(data.shape[0])
#                data = np.transpose(data, [2, 1, 0])
#                data_array = xr.DataArray(data, coords=[(x_axis_quantity, x_axis), ('y', y_axis),
#                                                  ('frames', frame_axis)], name='intensity')
#                data_array.frames.attrs['long_name'] = 'Frame number'
#
#            data_array.attrs['long_name'] = 'Intensity'
#            data_array.attrs['units'] = 'arb. u.'
#            data_array.y.attrs['units'] = 'px'
#            data_array[x_axis_quantity].attrs['long_name'] = x_name
#            data_array[x_axis_quantity].attrs['units'] = x_unit
#            data_array.attrs['sif_metadata'] = str(self)
#
#            return data_array
#
#        except ImportError:
#            raise RuntimeError("xarray package required for this method!")