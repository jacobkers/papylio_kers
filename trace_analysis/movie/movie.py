# -*- coding: utf-8 -*-
"""
Created on Jun 2019
@author: carlos

some features are shared within an image collection, like the file name, mapping, background, threshold, ...)
so imagecollection is a class ;)

warning blackboax number: (size+fwhm) gauss for extracting donor&acceptor
"""

import os
import time
from pathlib import Path
import tifffile as TIFF
import numpy as np
import matplotlib.pyplot as plt

# from trace_analysis.image_adapt.load_file import read_one_page#_pma, read_one_page_tif
# from trace_analysis.image_adapt.load_file import read_header
from trace_analysis.image_adapt.rolling_ball import rollingball
from trace_analysis.image_adapt.find_threshold import remove_background, get_threshold
# from trace_analysis.image_adapt.Mapping import Mapping
from trace_analysis.image_adapt.Image import Image

# from cached_property import cached_property
from trace_analysis.mapping.mapping import Mapping2

class Movie:
    def __init__(self, filepath):  # , **kwargs):
        self.filepath = Path(filepath)
        self._average_image = None
        self._maximum_projection_image = None
        self.is_mapping_movie = False

        self.channels = [['green', 'g', 'donor', 'd'],
                         ['red', 'r', 'acceptor', 'a']]
        self._channel_grid = np.array([2, 1])  # (x,y)
        self._number_of_channels = 2
        self.rot90 = 0

        if not self.filepath.suffix == '.sifx':
            self.writepath = self.filepath.parent
            self.name = self.filepath.with_suffix('').name

        self.read_header()

    def __repr__(self):
        return (f'{self.__class__.__name__}({str(self.filepath)})')

    @property
    def average_image(self):
        if self._average_image is None: self.make_average_image(write=True)
        return self._average_image

    @property
    def maximum_projection_image(self):
        if self._maximum_projection_image is None: self.make_maximum_projection(write=True)
        return self._maximum_projection_image

    @property
    def channel_grid(self):
        """ numpy.array : number of channels in the horizontal and vertical dimension

        Setting the channel_grid variable will assume equally spaced channels
        """
        return self._channel_grid

    @channel_grid.setter
    def channel_grid(self, channel_grid):
        channel_grid = np.array(channel_grid)
        # Possibly support multiple cameras by adding a third dimension
        if len(channel_grid) == 2 and np.all(np.array(channel_grid) > 0):
            self._channel_grid = channel_grid
            self._number_of_channels = np.product(channel_grid)

    @property
    def number_of_channels(self):
        """ int : number of channels in the movie

        Setting the number of channels will divide the image horizontally in equally spaced channels.
        """
        return self._number_of_channels

    @number_of_channels.setter
    def number_of_channels(self, number_of_channels):
        if number_of_channels > 0:
            self._number_of_channels = number_of_channels
            self._channel_grid = (number_of_channels, 1)
        else:
            raise ValueError('Number of channels should be at least 1')

    def read_header(self):
        # self.width_pixels, self.height_pixels, self.number_of_frames, self.movie_file_object = read_header(self.filepath)
        self._read_header()
        if not (self.rot90 % 2 == 0):
            width = self.width
            height = self.height
            self.width = height
            self.height = width

    def read_frame(self, frame_number):
        frame = self._read_frame(frame_number)
        return np.rot90(frame, self.rot90)

    def get_channel(self, image=None, channel='d'):
        if image is None: image = self.average_image
        channel_boundaries = self.channel_boundaries(channel)
        #
        #
        #     return image
        # else
        return image[channel_boundaries[0, 1]:channel_boundaries[1, 1],
               channel_boundaries[0, 0]:channel_boundaries[1, 0]]

    def get_channel_number(self, channel):
        """Get the channel number belonging to a specific channel (name)
        If

        Parameters
        ----------
        channel : str or int
            The name or number of a channel

        Returns
        -------
        i: int
            The index of the channel to which the channel name belongs

        """
        if isinstance(channel, int):
            # We should probably integrate this into the for loop
            if channel < self._number_of_channels:
                return channel
        for i, channel_names in enumerate(self.channels):
            if channel in channel_names:
                return i

    def channel_boundaries(self, channel):
        """Get the x and y boundaries of the channel within the movie

        Parameters
        ----------
        channel : str
            Name of a channel or 'all'

        Returns
        -------
        channel_boundaries : np.array
            Formatted as two coordinates, with the lowest and highest x and y values respectively
        """
        if channel == 'all':
            horizontal_boundaries = [0, self.width]
            vertical_boundaries = [0, self.height]
        else:
            channel_number = self.get_channel_number(channel)

            channel_width = self.width // self.channel_grid[0]
            horizontal_boundaries = np.array([0, channel_width]) + \
                                    channel_width * (channel_number % self.channel_grid[0])

            channel_height = self.height // self.channel_grid[1]
            vertical_boundaries = np.array([0, channel_height]) + \
                                  channel_height * (channel_number // self.channel_grid[0])

        return np.vstack([horizontal_boundaries, vertical_boundaries]).T

        # if channel is 'd':
        #     return np.array([[0, self.width // 2],[0,self.height]])
        # elif channel is 'a':
        #     return np.array([[self.width // 2, self.width], [0, self.height]])

    def channel_vertices(self, channel):
        """Get the vertices of the channel within the movie

        Parameters
        ----------
        channel : str
            Name of a channel or 'all'

        Returns
        -------
        channel_vertices : np.array
            Four coordinates giving the four corners of the channel
            Coordinates form a closed shape
        """
        if channel == 'all':
            channel_width = self.width
            channel_height = self.height
            channel_origin = [0, 0]
        else:
            channel_number = self.get_channel_number(channel)
            channel_width = self.width // self.channel_grid[0]
            channel_height = self.height // self.channel_grid[1]
            channel_origin = [channel_width * (channel_number % self.channel_grid[0]),
                              channel_height * (channel_number // self.channel_grid[0])]

        channel_vertices = np.array([channel_origin, ] * 4)
        channel_vertices[[1, 2], 0] += channel_width
        channel_vertices[[2, 3], 1] += channel_height

        return channel_vertices

    def saveas_tif(self):
        tif_filepath = self.writepath.joinpath(self.name + '.tif')
        try:
            tif_filepath.unlink()
            # When upgraded to python 3.8 the try-except struture can be replaced by
            # tif_filepath.unlink(missing_ok=True)
        except OSError:
            pass

        for i in range(self.number_of_frames):
            # frame = self.get_image(ii).image
            # frame = read_one_page(self.filepath, pageNb=i, A = self.movie_file_object)
            frame = self.read_frame(frame_number=i)
            print(i)
            # naam=r'M:\tnw\bn\cmj\Shared\margreet\Cy3 G50\ModifiedData\Python'+'{:03d}'.format(ii)+'.tif'
            TIFF.imwrite(tif_filepath, np.uint16(frame), append=True)

    def make_projection_image(self, type='average', start_frame=0, number_of_frames=20, write=False):
        """ Construct a projection image
        Determine a projection image for a number_of_frames starting at start_frame.
        i.e. [start_frame, start_frame + number_of_frames)

        Parameters
        ----------
        type : str
            'average' for average image
            'maximum' for maximum projection image
        start_frame : int
            Frame to start with
        number_of_frames : int
            Number of frames to average over
        write : bool
            If true, a tif file will be saved in writepath

        Returns
        -------
        np.ndarray
            2d image array with the projected image
        """
        # Check and specify number of frames
        if number_of_frames == 'all':
            number_of_frames = self.number_of_frames
        elif (self.number_of_frames - start_frame) < number_of_frames:
            print('Number of frames entered exceeds size movie')
            number_of_frames = (self.number_of_frames - start_frame)

        # Calculate sum of frames and find mean
        image = np.zeros((self.height, self.width))

        if type == 'average':
            for i in range(start_frame, start_frame + number_of_frames):
                frame = self.read_frame(frame_number=i).astype(float)
                image = image + frame
            image = (image / number_of_frames).astype(int)
            self._average_image = image
            if write:
                self.write_image(image, '_ave.tif')
        elif type == 'maximum':
            for i in range(start_frame, start_frame + number_of_frames):
                frame = self.read_frame(frame_number=i)
                image = np.maximum(image, frame)
            self._maximum_projection_image = image
            if write:
                self.write_image(image, '_max.tif')

        return image

    def make_average_image(self, start_frame=0, number_of_frames=20, write=False):
        """ Construct an average image
        Determine average image for a number_of_frames starting at start_frame.
        i.e. [start_frame, start_frame + number_of_frames)

        Parameters
        ----------
        start_frame : int
            Frame to start with
        number_of_frames : int
            Number of frames to average over
        write : bool
            If true, the a tif file will be saved in the writepath

        Returns
        -------
        np.ndarray
            2d image array with the average image

        """
        return self.make_projection_image('average', start_frame, number_of_frames, write)

    def make_maximum_projection(self, start_frame=0, number_of_frames=20, write=False):
        """ Construct a maximum projection image
        Determine maximum projection image for a number_of_frames starting at start_frame.
        i.e. [start_frame, start_frame + number_of_frames)

        Parameters
        ----------
        start_frame : int
            Frame to start with
        number_of_frames : int
            Number of frames to average over
        write : bool
            If true, the a tif file will be saved in the writepath

        Returns
        -------
        np.ndarray
            2d image array with the maximum projection image
        """

        return self.make_projection_image('maximum', start_frame, number_of_frames, write)

    def write_image(self, image, extension):
        """Write an image to the tif file format

        Parameters
        ----------
        image : np.ndarray
            2d image array
        extension : str
            String to add after the filename
        """
        if '.tif' not in extension:
            'Only tif export is supported (at the moment)'

        tif_filepath = self.writepath.joinpath(self.name + extension).with_suffix('.tif')
        if self.bitdepth == 16:
            TIFF.imwrite(tif_filepath, np.uint16(image))
        elif self.bitdepth == 8:
            TIFF.imwrite(tif_filepath, np.uint8(image))

    def show(self):
        return MoviePlotter(self)

    # Moved to file, can probably be removed
    def show_average_image(self, mode='2d', figure=None):
        if not figure: figure = plt.figure()  # Or possibly e.g. plt.figure('Movie')
        if mode == '2d':
            axis = figure.gca()
            axis.imshow(self.average_image)
        if mode == '3d':
            from matplotlib import cm
            axis = figure.gca(projection='3d')
            X = np.arange(self.average_image.shape[1])
            Y = np.arange(self.average_image.shape[0])
            X, Y = np.meshgrid(X, Y)
            axis.plot_surface(X, Y, self.average_image, cmap=cm.coolwarm,
                              linewidth=0, antialiased=False)
        # plt.show()

    def subtract_background(self, image, method='per_channel'):
        if method == 'rollingball':
            background = rollingball(image, self.width_pixels / 10)[1]  # this one is not used in pick_spots_akaze
            image_correct = image - background
            image_correct[image_correct < 0] = 0
            threshold = get_threshold(image_correct)
            return remove_background(image_correct, threshold)
        elif method == 'per_channel':  # maybe there is a better name
            sh = np.shape(image)
            threshold_donor = get_threshold(self.get_channel(image, 'donor'))
            threshold_acceptor = get_threshold(self.get_channel(image, 'acceptor'))
            background = np.zeros(np.shape(image))
            background[:, 0:sh[0] // 2] = threshold_donor
            background[:, sh[0] // 2:] = threshold_acceptor
            return remove_background(image, background)

        # note: optionally a fixed threshold can be set, like with IDL
        # note 2: do we need a different threshold for donor and acceptor?

    def show_coordinates(self, image, coordinates, figure=None, **kwargs):
        if not figure: figure = plt.gcf()  # Or possibly e.g. plt.figure('Movie')
        #        sorted_intensities = np.sort(image)
        #        vmin = np.percentile(sorted_intensities, 5)
        #        vmax = np.percentile(sorted_intensities, 99)
        plt.imshow(image, **kwargs)
        plt.scatter(coordinates[:, 0], coordinates[:, 1], marker='o', facecolors='none', edgecolors='r')
        # plt.show()

        plt.savefig(self.writepath.joinpath(self.name + '_ave_circles.png'), dpi=600)

    # Moved to coordinate_optimalization, so can probably be removed [IS 01-11-2019]
    def is_within_margin(self, coordinates,
                         edge=None,
                         margin=10):

        if edge is None: edge = np.array([[0, self.width // 2], [0, self.height]])
        if coordinates.size == 0: return np.array([])

        criteria = np.array([(coordinates[:, 0] > edge[0, 0] + margin),
                             (coordinates[:, 0] < edge[0, 1] - margin),
                             (coordinates[:, 1] > edge[1, 0] + margin),
                             (coordinates[:, 1] < edge[1, 1] - margin)
                             ])

        return criteria.all(axis=0)

    def get_channel_coordinates(self, within_margin=True, show=False, channel='donor', threshold=100):
        image = self.get_channel(self.average_image, channel)
        coordinates = self.find_peaks(image, method='local-maximum', threshold=threshold)
        coordinates = coordinates[self.is_within_margin(coordinates)]

        if show == True: self.show_coordinates(image, coordinates)

        return coordinates

    def calculate_acceptor_coordinates(self, donor_coordinates):
        # acceptor_coordinates = polywarp_apply(self.mapping.P,self.mapping.Q,donor_coordinates)
        acceptor_coordinates = self.mapping.transform_coordinates(donor_coordinates)
        return acceptor_coordinates

    def calculate_donor_coordinates(self, acceptor_coordinates):
        # donor_coordinates = polywarp_apply(self.mapping.P21,self.mapping.Q21,acceptor_coordinates)
        donor_coordinates = self.mapping.transform_coordinates(acceptor_coordinates)
        return donor_coordinates

    # Will be removed, as it is moved to file
    def write_coordinates_to_pks_file(self, coordinates):
        pks_filepath = self.writepath.joinpath(self.name + '.pks')
        with pks_filepath.open('w') as pks_file:
            for i, coordinate in enumerate(coordinates):
                # outfile.write(' {0:4.0f} {1:4.4f} {2:4.4f} {3:4.4f} {4:4.4f} \n'.format(i, coordinate[0], coordinate[1], 0, 0, width4=4, width6=6))
                pks_file.write('{0:4.0f} {1:4.4f} {2:4.4f} \n'.format(i + 1, coordinate[0], coordinate[1]))

    def generate_pks_file(self, channel):

        image_mean = self.make_average_image(number_of_frames=20, write=True)

        # image_mean_corrected = self.subtract_background(image_mean, method = 'per_channel')

        if channel == 'd':
            # donor_coordinates = self.find_peaks(self.get_channel(image_mean_corrected, 'donor'))
            donor_coordinates = self.find_peaks(self.get_channel(image_mean, 'donor'),
                                                method='local-maximum', threshold=self.threshold['point-selection'][0])
            acceptor_coordinates = self.calculate_acceptor_coordinates(donor_coordinates)
            acceptor_coordinates[:, 0] = acceptor_coordinates[:, 0] + self.width // 2

        elif channel == 'a':
            # acceptor_coordinates = self.find_peaks(self.get_channel(image_mean_corrected, 'acceptor'))
            acceptor_coordinates = self.find_peaks(self.get_channel(image_mean, 'acceptor'),
                                                   method='local-maximum',
                                                   threshold=self.threshold['point-selection'][0])
            donor_coordinates = self.calculate_donor_coordinates(acceptor_coordinates)
            acceptor_coordinates[:, 0] = acceptor_coordinates[:, 0] + self.width // 2

        elif channel == 'da':
            print('I have no clue yet how to do this')
            # most likely they do not overlap before finding transformation, so what is the point of doing D+A?
            # pts_number, label_size, ptsG = analyze_label.analyze(im_mean20_correctA[:, 0:self.height_pixels//2+im_mean20_correctA[:, self.height_pixels//2:]])

            # Ivo: I think they always transform the entire image of the acceptor channel using the mapping file,
            # so that they do overlap. From what I can see in IDL at least.

        else:
            print('make up your mind, choose wisely d/a/da')

            # Discard point close to edge image
        donor_edge = np.array([[0, self.width // 2], [0, self.height]])
        acceptor_edge = np.array([[self.width // 2, self.width], [0, self.height]])
        margin = 10

        both_within_margin = (self.is_within_margin(donor_coordinates, donor_edge, margin) &
                              self.is_within_margin(acceptor_coordinates, acceptor_edge, margin))

        donor_coordinates = donor_coordinates[both_within_margin]
        acceptor_coordinates = acceptor_coordinates[both_within_margin]

        all_coordinates = np.array([donor_coordinates, acceptor_coordinates])
        s = all_coordinates.shape
        all_coordinates = np.reshape(all_coordinates.T, (s[0], s[1] * s[2])).T

        self.write_coordinates_to_pks_file(all_coordinates)

    # Can likely be removed
    def use_for_mapping(self):
        self.is_mapping_movie = True
        donor = self.get_channel_coordinates(channel='donor',
                                             threshold=self.threshold['point-selection'][0])
        acceptor = self.get_channel_coordinates(channel='acceptor',
                                                threshold=self.threshold['point-selection'][1])
        self.mapping = Mapping2(donor, acceptor)

        return self.mapping

    def show_selected_spots(self):  # Not yet checked IS 09-09-2019
        # make a plot with selected spots
        # make image with found spots
        PL = plt.figure(14, figsize=(40, 40))
        plt.imshow(self.im_mean20_correct, vmax=np.amin(self.im_mean20_correct) + 5)
        for ii in range((np.amax(np.shape(self.dstG)))):
            plt.plot(self.ptsG[ii][0], self.ptsG[ii][1], 'wo', markerfacecolor='none', markersize=8)
            plt.plot(self.dstG[ii][0], self.dstG[ii][1], 'wv', markerfacecolor='none', markersize=8)
        for ii in range((np.amax(np.shape(self.ptsG2)))):
            plt.plot(self.ptsG2[ii][0], self.ptsG2[ii][1], 'y^', markerfacecolor='none', markersize=8)
        PL.savefig(self.filepath[:-4] + '-P data found spots.tif')

        PL = plt.figure(15, figsize=(40, 40))
        if self.choice_channel == 'd':
            for ii in range((np.amax(np.shape(self.ptsG2)))):
                plt.plot(self.ptsG2[ii][0] - len(self.im_mean20_correct) // 2, self.ptsG2[ii][1],
                         'r^')  # ,markerfacecolor='none', markersize=8)
        else:
            for ii in range((np.amax(np.shape(self.ptsG2)))):
                plt.plot(self.ptsG2[ii][0], self.ptsG2[ii][1], 'r^')  # ,markerfacecolor='none', markersize=8)
        for ii in range((np.amax(np.shape(self.dstG)))):
            plt.plot(self.ptsG[ii][0], self.ptsG[ii][1], 'ko', markerfacecolor='none', markersize=8)
            plt.plot(self.dstG[ii][0] - len(self.im_mean20_correct) // 2, self.dstG[ii][1], 'kv',
                     markerfacecolor='none', markersize=8)

        PL.savefig(self.filepath[:-4] + '-P data location spots.tif')

    def get_image(self, idx):  # Not yet checked IS 09-09-2019
        img = self.read_frame(idx)
        # img = self.subtract_background(img)
        return Image(img, self.height_pixels, self.mapping._tf2_matrix, self.ptsG, self.dstG, self.pts_number,
                     self.Gauss)

    def get_image_show(self, idx, hs, ws, siz):  # example hs=650,ws=950,siz=20  # Not yet checked IS 09-09-2019
        img = self.read_frame(idx)
        plt.figure(idx)
        ax1 = plt.subplot(1, 2, 1)
        ax1.imshow(img)
        ax1.set_xlim(hs, hs + siz)
        ax1.set_ylim(ws, ws + siz)

        img = self.subtract_background(img)
        ax2 = plt.subplot(1, 2, 2)
        ax2.imshow(img)
        ax2.set_xlim(hs, hs + siz)
        ax2.set_ylim(ws, ws + siz)
        return Image(img, self.height_pixels, self.mapping._tf2_matrix, self.ptsG, self.dstG, self.pts_number,
                     self.Gauss)

    # Will be removed, as it is moved to file
    def write_traces_to_traces_file(self, traces):
        traces_filepath = self.writepath.joinpath(self.name + '.traces')
        with traces_filepath.open('w') as traces_file:
            np.array([traces.shape[2]], dtype=np.int32).tofile(traces_file)
            np.array([traces.shape[0] * traces.shape[1]], dtype=np.int16).tofile(traces_file)
            # time_tr = np.zeros((self.number_of_frames, 2 * self.pts_number))
            # Ncolours=2
            # for jj in range(2*self.pts_number//Ncolours):
            #     time_tr[:,jj*2] = donor[:,jj]
            #     time_tr[:,jj*2+1]=  acceptor[:,jj]
            np.array(traces.T, dtype=np.int16).tofile(traces_file)


class MoviePlotter:
    # Adapted from Matplotlib Image Slices Viewer
    def __init__(self, movie):
        fig, ax = plt.subplots(1, 1)
        fig.canvas.mpl_connect('scroll_event', self.on_scroll)
        plt.show()

        self.ax = ax
        ax.set_title('use scroll wheel to navigate images')

        self.movie = movie
        self.slices, rows, cols = (movie.number_of_frames, movie.height, movie.width)
        self.ind = self.slices // 2

        self.im = ax.imshow(self.movie.read_frame(self.ind))
        self.update()

    def on_scroll(self, event):
        print("%s %s" % (event.button, event.step))
        if event.button == 'up':
            self.ind = (self.ind + 1) % self.slices
        else:
            self.ind = (self.ind - 1) % self.slices
        self.update()

    def update(self):
        self.im.set_data(self.movie.read_frame(self.ind))
        self.ax.set_ylabel('slice %s' % self.ind)
        self.im.axes.figure.canvas.draw()
