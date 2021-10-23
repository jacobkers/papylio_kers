import sys
from pathlib import Path
import pandas as pd
import tifffile as TIFF
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import xarray as xr

from trace_analysis.image_adapt.rolling_ball import rollingball
from trace_analysis.image_adapt.find_threshold import remove_background, get_threshold
from trace_analysis.timer import Timer

class Movie:
    def __init__(self, filepath):  # , **kwargs):
        self.filepath = Path(filepath)
        # self.filepaths = [Path(filepath) for filepath in filepaths] # For implementing multiple files, e.g. two channels over two files
        self._average_image = None
        self._maximum_projection_image = None
        self.is_mapping_movie = False

        self.illuminations = [Illumination(self, 'green', 'g')]
        self.illumination_arrangement = np.array([0])

        self.channels = [Channel(self, 'green', 'g', other_names=['donor', 'd']),
                         Channel(self, 'red', 'r', other_names=['acceptor', 'a'])]


        self.channel_arrangement = np.array([[[0, 1]]]) #[[[0,1]]] # First level: frames, second level: y within frame, third level: x within frame

        self.rot90 = 0

        self._data_type = np.dtype(np.uint16)
        self.intensity_range = (np.iinfo(self.data_type).min, np.iinfo(self.data_type).max)
        self.illumination_correction = None

        if not self.filepath.suffix == '.sifx':
            self.writepath = self.filepath.parent
            self.name = self.filepath.with_suffix('').name

        self.header_is_read = False

        # self.create_frame_info()

    def __repr__(self):
        return (f'{self.__class__.__name__}({str(self.filepath)})')

    @property
    def pixels_per_frame(self):
        return self.width*self.height

    @property
    def bitdepth(self):
        return self.data_type.itemsize*8 # 8 bits in a byte

    @property
    def bytes_per_frame(self):
        return self.data_type.itemsize * self.pixels_per_frame

    @property
    def average_image(self):
        if self._average_image is None:
            if self.filepath.with_name(self.name+'_ave.tif')
            self.make_average_image(write=True)
        return self._average_image

    @property
    def maximum_projection_image(self):
        if self._maximum_projection_image is None:
            self.make_maximum_projection(write=True)
        return self._maximum_projection_image

    # @property
    # def channel_grid(self):
    #     """ numpy.array : number of channels in the horizontal and vertical dimension
    #
    #     Setting the channel_grid variable will assume equally spaced channels
    #     """
    #     return self._channel_grid
    #
    # @channel_grid.setter
    # def channel_grid(self, channel_grid):
    #     channel_grid = np.array(channel_grid)
    #     # Possibly support multiple cameras by adding a third dimension
    #     if len(channel_grid) == 2 and np.all(np.array(channel_grid) > 0):
    #         self._channel_grid = channel_grid
    #         self._number_of_channels = np.product(channel_grid)
    @property
    def number_of_illuminations(self):
        """ int : number of channels in the movie

        Setting the number of channels will divide the image horizontally in equally spaced channels.
        """
        return len(self.illuminations)

    @property
    def number_of_channels(self):
        """ int : number of channels in the movie

        Setting the number of channels will divide the image horizontally in equally spaced channels.
        """
        return len(self.channels)

    @property
    def data_type(self):
        return self._data_type

    @data_type.setter
    def data_type(self, data_type):
        self._data_type = data_type
        self.intensity_range = (np.iinfo(self.data_type).min, np.iinfo(self.data_type).max)

    def create_frame_info(self):
        # files = [0] # For implementing multiple files
        frames = range(self.number_of_frames)

        index = pd.Index(data=frames, name='frame')
        self.frame_info = pd.DataFrame(index=index, columns=['time', 'illumination', 'channel'])
        # self.frame_info['file'] = len(self.frame_info) * [list(range(2))] # For implementing multiple files
        # self.frame_info = self.frame_info.explode('file') # For implementing multiple files
        self.frame_info['time'] = self.frame_info.index.to_frame()['frame'].values
        self.frame_info['illumination'] = self.illumination_arrangement.tolist() * (self.number_of_frames // self.illumination_arrangement.shape[0])
        self.frame_info['channel'] = self.channel_arrangement.tolist() * (self.number_of_frames // self.channel_arrangement.shape[0])

        self.frame_info = self.frame_info.explode('channel').explode('channel')

        categorical_columns = ['illumination', 'channel']
        self.frame_info[categorical_columns] = self.frame_info[categorical_columns].astype('category')

    def read_header(self):
        self._read_header()
        if not (self.rot90 % 2 == 0):
            width = self.width
            height = self.height
            self.width = height
            self.height = width

        self.header_is_read = True

    def read_frame_raw(self, frame_number):
        if not self.header_is_read:
            self.read_header()
        frame = self._read_frame(frame_number)
        return np.rot90(frame, self.rot90)

    def read_frame(self, frame_number, channel=None):
        frame = self.read_frame_raw(frame_number)
        frame = xr.DataArray(np.dstack([channel.crop_image(frame) for channel in self.channels]),
                             dims=('x', 'y', 'channel'),
                             coords={'channel': [channel.index for channel in self.channels]})

        channels = self.get_channels_from_names(channel)
        channel_indices = self.get_channel_indices_from_names(channel)
        frame = frame.sel(channel=channel_indices)

        if self.illumination_correction is not None:
            for channel in channels:
                frame = frame.astype(float)
                frame.loc[{'channel': channel.index}] *= self.illumination_correction[frame_number, channel.index]

        frame_out = np.zeros((self.height, self.width))
        for channel in channels:
            frame_out[channel.boundaries[0, 1]:channel.boundaries[1, 1],
                      channel.boundaries[0, 0]:channel.boundaries[1, 0]] = frame.sel(channel=channel.index).values

        return frame_out

    def get_channel(self, image=None, channel='d'):
        if image is None:
            image = self.average_image
        if channel in [None, 'all']:
            return image

        if not isinstance(channel, Channel):
            channel = self.get_channel_from_name(channel)

        return channel.crop_image(image)

    def get_channel_from_name(self, channel_name):
        """Get the channel index belonging to a specific channel (name)
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
        for channel in self.channels:
            if channel_name in channel.names:
                return channel
        else:
            raise ValueError('Channel name not found')

    def get_channels_from_names(self, channel_names):
        """Get the channel index belonging to a specific channel (name)
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
        if channel_names in [None, 'all']:
            return self.channels

        return [self.get_channel_from_name(channel_name) for channel_name in channel_names]

    def get_channel_indices_from_names(self, channel_names):
        channels = self.get_channels_from_names(channel_names)
        return [channel.index for channel in channels]

    def saveas_tif(self):
        tif_filepath = self.writepath.joinpath(self.name + '.tif')
        try:
            tif_filepath.unlink()
            # When upgraded to python 3.8 the try-except structure can be replaced by
            # tif_filepath.unlink(missing_ok=True)
        except OSError:
            pass

        for i in range(self.number_of_frames):
            frame = self.read_frame(frame_number=i)
            TIFF.imwrite(tif_filepath, np.uint16(frame), append=True)

            #     tifffile.imwrite(self.writepath.joinPath(f'{self.name}_fr{frame_number}.tif'), image,  photometric='minisblack')

    def make_projection_image(self, projection_type='average', start_frame=0, number_of_frames=20, illumination=None,
                              channel=None, write=False):
        """ Construct a projection image
        Determine a projection image for a number_of_frames starting at start_frame.
        i.e. [start_frame, start_frame + number_of_frames)

        Parameters
        ----------
        projection_type : str
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
        if not self.header_is_read:
            self.read_header()

        frames = self.frame_info
        frames = frames.loc[start_frame:]
        filename_addition = ''

        if illumination is not None:
            frames = frames.query(f'illumination=={illumination}')
            if self.number_of_illuminations > 1:
                filename_addition += f'_i{illumination}'
        if channel is not None:
            if not isinstance(channel, Channel):
                channel = self.get_channel_from_name(channel)
            frames = frames.query(f'channel=={channel.index}')
            filename_addition += f'_c{channel.index}'

        # Determine frame indices to be used
        frame_indices = frames.index.unique().values
        if number_of_frames == 'all':
            pass
        elif type(number_of_frames) is int:
            if number_of_frames > len(frame_indices):
                print(f'Number of frames entered exceeds available frames, used {len(frame_indices)} instead of {number_of_frames} frames')
            frame_indices = frame_indices[:number_of_frames]
        else:
            raise ValueError('Incorrect value for number_of_frames')

        # Calculate sum of frames and find mean
        image = self.get_channel(np.zeros((self.height, self.width)), channel=channel)
        number_of_frames = len(frame_indices)

        if projection_type == 'average':
            if len(frame_indices) > 100:
                print(f'\n Making average image of {self.name}')
            for i, frame_index in enumerate(frame_indices):
                if len(frame_indices) > 100 and i % 13 == 0:
                    sys.stdout.write(f'\r   Processing frame {frame_index} in {frame_indices[0]}-{frame_indices[-1]}')
                frame = self.read_frame(frame_number=frame_index, channel=channel).astype(float)
                image = image + frame
            image = (image / number_of_frames).astype(self.data_type)
            self._average_image = image # Possibly we should save only the overview image not the last image [IS: 20-04-2021]
        elif projection_type == 'maximum':
            print(f'\n Making maximum projection image of {self.name}')
            for i, frame_index in enumerate(frame_indices):
                if i % 13 == 0:
                    sys.stdout.write(f'\r   Processing frame {frame_index} in {frame_indices[0]}-{frame_indices[-1]}')
                frame = self.read_frame(frame_number=frame_index)
                image = np.maximum(image, frame)
            sys.stdout.write(f'\r   Processed frames {frame_indices[0]}-{frame_indices[-1]}\n')
            self._maximum_projection_image = image # Possibly we should save only the overview image not the last image [IS: 20-04-2021]

        if write:
            filename = self.name + '_'+projection_type[:3]+f'_{number_of_frames}fr'+filename_addition
            filepath = self.writepath.joinpath(filename)
            TIFF.imwrite(filepath.with_suffix('.tif'), image)
            # plt.imsave(filepath.with_suffix('.tif'), image, format='tif', cmap=colour_map, vmin=self.intensity_range[0], vmax=self.intensity_range[1])
            # plt.imsave(filepath.with_suffix('.png'), image, cmap=colour_map, vmin=self.intensity_range[0], vmax=self.intensity_range[1])

        # return image

    def make_projection_images(self, projection_type='average', start_frame=0, number_of_frames=20):
        illumination_indices, channel_indices = \
            self.frame_info[['illumination','channel']].drop_duplicates().sort_values(by=['channel','illumination']).values.T

        # for illumination_index in np.unique(illumination_indices):
        #     self.make_projection_image(projection_type, start_frame, number_of_frames, illumination_index, write=True)

        images = []
        for illumination_index, channel_index in zip(illumination_indices, channel_indices):
            self.make_projection_image(projection_type, start_frame, number_of_frames,
                                       illumination_index, channel_index, write=True)

            image = self.make_projection_image(projection_type, start_frame, number_of_frames,
                                               illumination_index, channel_index)
            image = (image - self.intensity_range[0]) / (self.intensity_range[1]-self.intensity_range[0])
            images.append(self.channels[channel_index].colour_map(image, bytes=True))

        images_combined = np.hstack(images)
        filepath = self.writepath.joinpath(self.name + '_' + projection_type[:3] + f'_{number_of_frames}fr')
        plt.imsave(filepath.with_suffix('.png'), images_combined)

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
        return self.make_projection_image('average', start_frame=start_frame, number_of_frames=number_of_frames,
                                          write=write)

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

        return self.make_projection_image('maximum', start_frame=start_frame, number_of_frames=number_of_frames,
                                          write=write)

    def show(self):
        return MoviePlotter(self)

    # Moved to file, can probably be removed
    def show_projection_image(self, projection_type='average', figure=None):
        # This is only the last saved projection image, which may not be what we want here
        # We may want an standard overview image [IS: 20-04-2021]
        if not figure:
            figure = plt.figure()
        axis = figure.gca()

        if projection_type == 'average':
            image = self.average_image
            axis.set_title('Average image')
        elif projection_type == 'maximum':
            image = self.maximum_projection_image
            axis.set_title('Maximum projection')

        axis.imshow(image)

    def show_average_image(self, figure=None):
        self.show_projection_image(projection_type='average', figure=figure)

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

    def determine_illumination_correction(self, filter_neighbourhood_size=10):
        import scipy.ndimage.filters as filters
        illumination_intensity = np.zeros((self.number_of_frames, self.number_of_channels))

        for i in range(self.number_of_frames):
            frame = self.read_frame_raw(i)

            filtered_frame = filters.minimum_filter(frame, filter_neighbourhood_size)
            illumination_intensity[i, 0] = np.sum(self.get_channel(filtered_frame, 'g'))
            illumination_intensity[i, 1] = np.sum(self.get_channel(filtered_frame, 'r'))

        # figure = plt.figure()
        # axis = figure.gca()
        # axis.plot(illumination_intensity[:, 0], 'g', illumination_intensity[:, 1], 'r')
        # axis.set_title('Minimum filtered intensity sum per frame')
        # axis.set_xlabel('Frame (0.1s)')
        # axis.set_ylabel('Minimum filtered intensity sum')
        # axis.set_ylim((0, None))

        self.illumination_correction = illumination_intensity.max(axis=0) / illumination_intensity


class Channel:
    def __init__(self, movie, name, short_name, other_names=None, colour_map=None):
        self.movie = movie
        self.name = name
        self.short_name = short_name
        self.other_names = other_names
        if colour_map is None:
            channel_colour = list({'green', 'red', 'blue'}.intersection(self.names))[0]
            self.colour_map = make_colour_map(channel_colour)

    def __repr__(self):
        return (f'{self.__class__.__name__}({self.name})')

    @property
    def names(self):
        return [self.index, self.name, self.short_name] + self.other_names

    @property
    def index(self):
        try:
            return self.movie.channels.index(self)
        except:
            pass

    @property
    def location(self):
        return [int(i) for i in np.where(self.movie.channel_arrangement == self.index)]

    @property
    def width(self):
        return self.movie.width // self.movie.channel_arrangement.shape[2]

    @property
    def height(self):
        return self.movie.height // self.movie.channel_arrangement.shape[1]
        # for frame_index, frame in enumerate(self.channel_arrangement):
        #     for y_index, y in enumerate(frame):
        #         try:
        #             x_index = y.index(channel_index)
        #             return frame_index, y_index, x_index
        #         except ValueError:
        #             pass

    @property
    def origin(self):
        return [self.width * self.location[2],
                self.height * self.location[1]]

    @property
    def boundaries(self):
        # channel_boundaries: np.array
        # #         Formatted as two coordinates, with the lowest and highest x and y values respectively
        horizontal_boundaries = np.array([0, self.width]) + self.width * self.location[2]
        vertical_boundaries = np.array([0, self.height]) + self.height * self.location[1]
        return np.vstack([horizontal_boundaries, vertical_boundaries]).T

    @property
    def vertices(self):
        #     channel_vertices : np.array
        #         Four coordinates giving the four corners of the channel
        #         Coordinates form a closed shape
        channel_vertices = np.array([self.origin, ] * 4)
        channel_vertices[[1, 2], 0] += self.width
        channel_vertices[[2, 3], 1] += self.height
        return channel_vertices

    def crop_image(self, image):
        return image[self.boundaries[0, 1]:self.boundaries[1, 1],
                     self.boundaries[0, 0]:self.boundaries[1, 0]]

class Illumination:
    def __init__(self, movie, name, short_name, other_names=None):
        self.movie = movie
        self.name = name
        self.short_name = short_name
        self.other_names = other_names

    def __repr__(self):
        return (f'{self.__class__.__name__}({self.name})')

    @property
    def names(self):
        return [self.index, self.name, self.short_name] + self.other_names

    @property
    def index(self):
        try:
            return self.movie.illuminations.index(self)
        except:
            pass

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


def make_colour_map(colour, N=256):
    values = np.zeros((N, 3))
    if colour == 'grey':
        values[:, 0] = values[:, 1] = values[:, 2] = np.linspace(0, 1, N)
    elif colour == 'red':
        values[:, 0] = np.linspace(0, 1, N)
    elif colour == 'green':
        values[:, 1] = np.linspace(0, 1, N)
    elif colour == 'blue':
        values[:, 2] = np.linspace(0, 1, N)
    else:
        values[:, 0] = values[:, 1] = values[:, 2] = np.linspace(0, 1, N)

    return ListedColormap(values)

