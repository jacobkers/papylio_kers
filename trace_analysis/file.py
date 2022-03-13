if __name__ == '__main__':
    import sys
    from pathlib import Path
    p = Path(__file__).parents[1]
    sys.path.insert(0, str(p))

from pathlib import Path # For efficient path manipulation
import numpy as np #scientific computing with Python
import pandas as pd
import matplotlib.pyplot as plt #Provides a MATLAB-like plotting framework
import skimage.io as io
import xarray as xr
import skimage as ski
import warnings
import sys
# from trace_analysis.molecule import Molecule
from trace_analysis.movie.sifx import SifxMovie
from trace_analysis.movie.pma import PmaMovie
from trace_analysis.movie.nsk import NskMovie
from trace_analysis.movie.tif import TifMovie
from trace_analysis.movie.nd2 import ND2Movie
from trace_analysis.movie.binary import BinaryMovie
from trace_analysis.plotting import histogram
from trace_analysis.mapping.mapping import Mapping2
from trace_analysis.peak_finding import find_peaks
from trace_analysis.coordinate_optimization import  coordinates_within_margin, \
                                                    coordinates_after_gaussian_fit, \
                                                    coordinates_without_intensity_at_radius, \
                                                    merge_nearby_coordinates, \
                                                    set_of_tuples_from_array, array_from_set_of_tuples, \
                                                    coordinates_within_margin_selection
from trace_analysis.trace_extraction import extract_traces
from trace_analysis.coordinate_transformations import translate, transform # MD: we don't want to use this anymore I think, it is only linear
                                                                           # IS: We do! But we just need to make them usable with the nonlinear mapping
from trace_analysis.background_subtraction import extract_background
# from trace_analysis.plugin_manager import PluginManager
# from trace_analysis.plugin_manager import PluginMetaClass
from trace_analysis.plugin_manager import plugins


@plugins
class File:
    # plugins = []
    # _plugin_mixin_class = None
    #
    # @classmethod
    # def add_plugin(cls, plugin_class):
    #     cls.plugins.append(plugin_class)
    #     cls._plugin_mixin_class = type(cls.__name__, (cls,) + tuple(cls.plugins), {})
    #
    # def __new__(cls, *args, **kwargs):
    #     if not cls._plugin_mixin_class:
    #         return super().__new__(cls)
    #     else:
    #         return super().__new__(cls._plugin_mixin_class)

    def __init__(self, relativeFilePath, experiment, *args, **kwargs):
        relativeFilePath = Path(relativeFilePath)
        self.experiment = experiment

        self.relativePath = relativeFilePath.parent
        self.name = relativeFilePath.name
        self.extensions = set()

        # self.molecules = Molecules()

        self.exposure_time = None  # Found from log file or should be inputted

        self.log_details = None  # a string with the contents of the log file
        self.number_of_frames = None

        self.isSelected = False
        self.is_mapping_file = False

        self.movie = None
        self.mapping = None
        self._average_image = None
        self._maximum_projection_image = None

        if 'fov_info' in kwargs:
            self.nd2_fov_info = kwargs['fov_info']  # fov = Field of View

        self.dataset_variables = ['molecule', 'coordinates', 'background', 'intensity', 'FRET', 'selected', 'molecule_in_file']

        # I think it will be easier if we have import functions for specific data instead of specific files.
        # For example. the sifx, pma and tif files can better be handled in the Movie class. Here we then just have a method import_movie.
        # [IS 10-08-2020]
        # TODO: Make an import_movie method and move the specific file type handling to the movie class (probably this should also include the log file)
        # TODO: Make an import_mapping method and move the specific mapping type handling (.map, .coeff) to the mapping class.

        self.importFunctions = {'.sifx': self.import_sifx_file,
                                '.pma': self.import_pma_file,
                                '.nsk': self.import_nsk_file,
                                '.nd2': self.import_nd2_file,
                                '.tif': self.import_tif_file,
                                '.tiff': self.import_tif_file,
                                '.bin': self.import_bin_file,
                                '.TIF': self.import_tif_file,
                                '.TIFF': self.import_tif_file,
                                '_ave.tif': self.import_average_tif_file,
                                '_max.tif': self.import_maximum_projection_tif_file,
                                '.coeff': self.import_coeff_file,
                                '.map': self.import_map_file,
                                '.mapping': self.import_mapping_file,
                                '.pks': self.import_pks_file,
                                '.traces': self.import_traces_file,
                                '.log': self.import_log_file,
                                '_steps_data.xlsx': self.import_excel_file,
                                '_selected_molecules.txt': self.import_selected,
                                '.nc': self.noneFunction
                                }

        print(self)

        super().__init__()

        # if self.experiment.import_all is True:
        #     self.findAndAddExtensions()


    def __repr__(self):
        return (f'{self.__class__.__name__}({self.relativePath.joinpath(self.name)})')

    @property
    def relativeFilePath(self):
        return self.relativePath.joinpath(self.name)

    @property
    def absoluteFilePath(self):
        return self.experiment.main_path.joinpath(self.relativeFilePath)

    @property
    def number_of_molecules(self):
        return len(self.molecule)

    # @property
    # def molecule(self):
    #     with xr.open_dataset(self.relativeFilePath.with_suffix('.nc'), engine='h5netcdf') as dataset:
    #         return dataset['molecule'].load()
    # @number_of_molecules.setter
    # def number_of_molecules(self, number_of_molecules):
    #     if not self.molecules:
    #         for molecule in range(0, number_of_molecules):
    #             self.addMolecule()
    #     elif number_of_molecules != self.number_of_molecules:
    #         raise ValueError(f'Requested number of molecules ({number_of_molecules}) differs from existing number of '
    #                          f'molecules ({self.number_of_molecules}) in {self}. \n'
    #                          f'If you are sure you want to proceed, empty the molecules list file.molecules = [], or '
    #                          f'possibly delete old pks or traces files')

    @property
    def number_of_channels(self):
        return self.experiment.number_of_channels

    @property
    def selected_molecules(self):
        return self.molecule[self.selected]

    @property
    def average_image(self):
        if self._average_image is None:
            # Refresh configuration
            self.experiment.import_config_file()
            # number_of_frames = self.experiment.configuration['compute_image']['number_of_frames']
            configuration_show_movie = self.experiment.configuration['show_movie']
            if configuration_show_movie['illumination'] == 'None' or configuration_show_movie['illumination'] == 'none':
                illumination = None
            else:
                illumination = int(configuration_show_movie['illumination'])
            first_frame = int(configuration_show_movie['frames_for_show_movie']['first_frame'])
            if configuration_show_movie['frames_for_show_movie']['last_frame'] == 'last':
                last_frame = self.number_of_frames-1
            else:
                last_frame = int(configuration_show_movie['frames_for_show_movie']['last_frame'])
            number_of_frames = last_frame-first_frame + 1

            self._average_image = self.movie.make_average_image(start_frame=first_frame,
                                                                number_of_frames=number_of_frames,
                                                                illumination=illumination, write=True)
        return self._average_image

    @property
    def maximum_projection_image(self):
        if self._maximum_projection_image is None:
            # Refresh configuration
            self.experiment.import_config_file()
            # number_of_frames = self.experiment.configuration['compute_image']['number_of_frames']
            configuration_show_movie = self.experiment.configuration['show_movie']
            if configuration_show_movie['illumination'] == 'None' or configuration_show_movie['illumination'] == 'none':
                illumination = None
            else:
                illumination = int(configuration_show_movie['illumination'])
            first_frame = int(configuration_show_movie['frames_for_show_movie']['first_frame'])
            if configuration_show_movie['frames_for_show_movie']['last_frame'] == 'last':
                last_frame = self.number_of_frames-1
            else:
                last_frame = int(configuration_show_movie['frames_for_show_movie']['last_frame'])
            number_of_frames = last_frame-first_frame + 1

            self._maximum_projection_image = self.movie.make_maximum_projection(start_frame=first_frame,
                                                                                number_of_frames=number_of_frames,
                                                                                illumination=illumination, write=True)
        return self._maximum_projection_image

    # @property
    # def coordinates(self):
    #     with xr.open_dataset(self.relativeFilePath.with_suffix('.nc'), engine='h5netcdf') as dataset:
    #         #.set_index({'molecule': ('molecule_in_file','file')})
    #         return dataset['coordinates'].load()
    #
    # @coordinates.setter

    def set_coordinates_of_channel(self, coordinates, channel):
        # TODO: make this usable for more than two channels
        channel_index = self.movie.get_channel_number(channel)  # Or possibly make it self.channels
        if channel_index == 0:
            coordinates_in_main_channel = coordinates
        elif channel_index == 1:
            coordinates_in_main_channel = self.mapping.transform_coordinates(coordinates, inverse=True)
        else:
            raise NotImplementedError('File.set_coordinates_of_channel not implemented for more than two channels')

        # for channel in self.channels or for mapping in self.mappings
        coordinates_list = [coordinates_in_main_channel]
        for i in range(1, self.number_of_channels):
            if self.number_of_channels > 2:
                raise NotImplementedError('File.set_coordinates_of_channel not implemented for more than two channels')
            coordinates_in_other_channel = self.mapping.transform_coordinates(coordinates_in_main_channel, direction='Donor2Acceptor')
            coordinates_list.append(coordinates_in_other_channel)
        coordinates = np.hstack(coordinates_list).reshape((-1, 2))

        self.coordinates = coordinates

    def coordinates_from_channel(self, channel):
        # if not self._pks_file:
        #     _pks_file = PksFile(self.absoluteFilePath.with_suffix('.pks'))

        #return np.concatenate([[molecule.coordinates[0, :] for molecule in self.molecules]])
        if type(channel) is str:
            channel = {'d': 0, 'a': 1, 'g':0, 'r':1}[channel]

        return self.coordinates.sel(channel=channel)

    def __getattr__(self, item):
        if item in self.dataset_variables:
            with xr.open_dataset(self.relativeFilePath.with_suffix('.nc'), engine='h5netcdf') as dataset:
                return dataset[item].load()
        # else:
        #     super().__getattribute__(item)
        raise AttributeError

    def get_data(self, key):
        with xr.open_dataset(self.relativeFilePath.with_suffix('.nc'), engine='h5netcdf') as dataset:
            return dataset[key].load()

    @property
    def dataset(self):
        with xr.open_dataset(self.relativeFilePath.with_suffix('.nc'), engine='h5netcdf') as dataset:
            return dataset.load()

    # def get_coordinates(self, selected=False):
    #     if selected:
    #         molecules = self.selectedMolecules
    #     else:
    #         molecules = self.molecules
    #
    #     return np.vstack([molecule.coordinates for molecule in molecules])

    #in analogy with coordinates, also background:
    # @background.setter
    # def background(self, background, number_of_channels = None):
    #     if number_of_channels is None:
    #         number_of_channels = self.number_of_channels
    #     self.number_of_molecules = np.shape(background)[0]//number_of_channels
    #
    #     for i, molecule in enumerate(self.molecules):
    #         molecule.background = background[(i * number_of_channels):((i + 1) * number_of_channels)]
    #
    # def background_from_channel(self, channel):
    #     if type(channel) is str:
    #         channel = {'d': 0, 'a': 1, 'g':0, 'r':1}[channel]
    #
    #     return np.vstack([molecule.background[channel] for molecule in self.molecules])
    #

    @property
    def time(self):  # the time axis of the experiment, if not found in log it will be asked as input
        if self.exposure_time is None:
            # SHK modification for debuging. Should be set back to original code later
            print('exposure_time is set to 0.1s (see file.time())')
            self.exposure_time = 0.1
            # original code:
            # self.exposure_time = float(input(f'Exposure time for {self.name}: '))
        return np.arange(0, self.number_of_frames)*self.exposure_time

    def _init_dataset(self, number_of_molecules):
        dataset = xr.Dataset(
            {
                'selected': ('molecule', xr.DataArray(False, coords=[range(number_of_molecules)]))#,
                # 'x': (('molecule', 'channel'), xr.DataArray(np.nan, coords=[molecule_multiindex, []])),
                # 'y': (('molecule', 'channel'), xr.DataArray(np.nan, coords=[molecule_multiindex, []])),
                # 'background': (('molecule', 'channel'), xr.DataArray(np.nan, coords=[molecule_multiindex, []])),
                # 'intensity': (
                # ('molecule', 'channel', 'frame'), xr.DataArray(np.nan, coords=[molecule_multiindex, [], []]))
            },
            coords=
            {
                'molecule': ('molecule', range(number_of_molecules)),
                # # pd.MultiIndex.from_tuples([], names=['molecule_in_file', 'file'])),
                # 'frame': ('frame', np.array([], dtype=int)),
                # 'channel': ('channel', np.array([], dtype=int))
            }
        )
        dataset = dataset.reset_index('molecule').rename(molecule_='molecule_in_file')
        dataset = dataset.assign_coords({'file': ('molecule', [str(self.relativeFilePath)]*number_of_molecules)})

        dataset.to_netcdf(self.relativeFilePath.with_suffix('.nc'), engine='h5netcdf', mode='w')
        self.extensions.add('.nc')

    def findAndAddExtensions(self):
        file_keyword = self.name
        # special treatment for multi-fov nd2 file
        if hasattr(self, 'nd2_fov_info'):
            token_position = self.name.find('_fov')
            file_keyword = self.name[:token_position]

        foundFiles = [file.name for file in self.experiment.main_path.joinpath(self.relativePath).glob(file_keyword + '*')]
        # foundExtensions = [file[len(self.name):] for file in foundFiles]
        foundExtensions = [file[len(file_keyword):] for file in foundFiles]

        # For the special case of a sifx file, which is located inside a folder
        if '' in foundExtensions: foundExtensions[foundExtensions.index('')] = '.sifx'

        newExtensions = [extension for extension in foundExtensions if extension not in self.extensions]
        # self.extensions = self.extensions + newExtensions
        for extension in newExtensions: self.importExtension(extension)

    def importExtension(self, extension):

        # print(f.relative_to(self.experiment.main_path))

        # if extension not in self.extensions: # better to use sets here
        #     self.extensions.append(extension)

        # print(extension)

        self.importFunctions.get(extension, self.noneFunction)()
        if extension in self.importFunctions.keys():
            self.extensions.add(extension)

    def noneFunction(self):
        return

    def import_log_file(self):
        self.exposure_time = np.genfromtxt(f'{self.relativeFilePath}.log', max_rows=1)[2]
        print(f'Exposure time set to {self.exposure_time} sec for {self.name}')
        self.log_details = open(f'{self.relativeFilePath}.log').readlines()
        self.log_details = ''.join(self.log_details)

    def import_sifx_file(self):
        imageFilePath = self.absoluteFilePath.joinpath('Spooled files.sifx')
        self.movie = SifxMovie(imageFilePath)
        # self.movie.number_of_channels = self.experiment.number_of_channels
        self.number_of_frames = self.movie.number_of_frames

    def import_pma_file(self):
        imageFilePath = self.absoluteFilePath.with_suffix('.pma')
        self.movie = PmaMovie(imageFilePath)
        # self.movie.number_of_channels = self.experiment.number_of_channels
        self.number_of_frames = self.movie.number_of_frames

    def import_nsk_file(self):
        imageFilePath = self.absoluteFilePath.with_suffix('.nsk')
        self.movie = NskMovie(imageFilePath)
        # self.movie.number_of_channels = self.experiment.number_of_channels
        self.number_of_frames = self.movie.number_of_frames

    def import_tif_file(self):
        #TODO: Pass all image files to the Movie class and let the Movie class decide what to do
        imageFilePath = self.absoluteFilePath.with_suffix('.tif')
        self.movie = TifMovie(imageFilePath)
        # self.movie.number_of_channels = self.experiment.number_of_channels
        self.number_of_frames = self.movie.number_of_frames

    def import_nd2_file(self):
        if hasattr(self, 'nd2_fov_info'):
            # special treatment for multi-fov nd2 file
            token_position = self.name.find('_fov')
            file_keyword = self.name[:token_position]
            imageFilePath = self.absoluteFilePath.parent.joinpath(file_keyword).with_suffix('.nd2')
            self.movie = ND2Movie(imageFilePath, fov_info=self.nd2_fov_info)
        else:
            imageFilePath = self.absoluteFilePath.with_suffix('.nd2')
            self.movie = ND2Movie(imageFilePath)
        self.number_of_frames = self.movie.number_of_frames

    def import_bin_file(self):
        imageFilePath = self.absoluteFilePath.with_suffix('.bin')
        self.movie = BinaryMovie(imageFilePath)
        self.number_of_frames = self.movie.number_of_frames

    def import_average_tif_file(self):
        averageTifFilePath = self.absoluteFilePath.with_name(self.name+'_ave.tif')
        self._average_image = io.imread(averageTifFilePath, as_gray=True)

    def import_maximum_projection_tif_file(self):
        maxTifFilePath = self.absoluteFilePath.with_name(self.name+'_max.tif')
        self._maximum_projection_image = io.imread(maxTifFilePath, as_gray=True)

    def import_coeff_file(self):
        from skimage.transform import AffineTransform
        if self.mapping is None: # the following only works for 'linear'transformation_type
            file_content=np.genfromtxt(str(self.relativeFilePath) + '.coeff')
            if len(file_content)==12:
                [coefficients, coefficients_inverse] = np.split(file_content,2)
            elif len(file_content)==6:
                coefficients = file_content
            else:
                raise TypeError('Error in importing coeff file, wrong number of lines')

            self.mapping = Mapping2(transformation_type='linear')

            transformation = np.zeros((3,3))
            transformation[2,2] = 1
            transformation[[0,0,0,1,1,1],[2,0,1,2,0,1]] = coefficients
            self.mapping.transformation = AffineTransform(matrix=transformation)

            if len(file_content)==6:
                self.mapping.transformation_inverse=\
                    AffineTransform(matrix=np.linalg.inv(self.mapping.transformation.params))
            else:
                transformation_inverse = np.zeros((3,3))
                transformation_inverse[2,2] = 1
                transformation_inverse[[0,0,0,1,1,1],[2,0,1,2,0,1]] = coefficients_inverse
                self.mapping.transformation_inverse = AffineTransform(matrix=transformation_inverse)

            self.mapping.file = self
            self.mapping.source_name = 'Donor'
            self.mapping.destination_name = 'Acceptor'

    def export_coeff_file(self):
        warnings.warn('The export_coeff_file method will be depricated, use export_mapping instead')
        self.export_mapping(filetype='classic')

    def export_mapping(self, filetype='yml'):
        self.mapping.save(self.absoluteFilePath, filetype)

    def import_map_file(self):
        from trace_analysis.mapping.polywarp import PolywarpTransform
        #coefficients = np.genfromtxt(self.relativeFilePath.with_suffix('.map'))
        file_content=np.genfromtxt(self.relativeFilePath.with_suffix('.map'))
        if len(file_content) == 64:
            [coefficients, coefficients_inverse] = np.split(file_content, 2)
        elif len(file_content) == 32:
            coefficients = file_content
        else:
            raise TypeError('Error in import map file, incorrect number of lines')

        degree = int(np.sqrt(len(coefficients) // 2) - 1)
        P = coefficients[:len(coefficients) // 2].reshape((degree + 1, degree + 1))
        Q = coefficients[len(coefficients) // 2 : len(coefficients)].reshape((degree + 1, degree + 1))

        self.mapping = Mapping2(transformation_type='nonlinear')
        self.mapping.transformation = PolywarpTransform(params=(P,Q)) #{'P': P, 'Q': Q}
        #self.mapping.file = self

        if len(file_content)==64:
            degree = int(np.sqrt(len(coefficients_inverse) // 2) - 1)
            Pi = coefficients_inverse[:len(coefficients_inverse) // 2].reshape((degree + 1, degree + 1))
            Qi = coefficients_inverse[len(coefficients_inverse) // 2 : len(coefficients_inverse)].reshape((degree + 1, degree + 1))
        else:
            grid_range = 500 # in principle the actual image size doesn't matter
            # image_height = self._average_image.shape[0]

            # Can't we make this independent of the image?
            grid_coordinates = np.array([(a,b) for a in np.arange(0, grid_range//2, 5) for b in np.arange(0, grid_range, 5)])
            from trace_analysis.mapping.polywarp import polywarp, polywarp_apply
            transformed_grid_coordinates = polywarp_apply(P, Q, grid_coordinates)
            # plt.scatter(grid_coordinates[:, 0], grid_coordinates[:, 1], marker='.')
            # plt.scatter(transformed_grid_coordinates[:,0], transformed_grid_coordinates[:,1], marker='.')
            Pi, Qi = polywarp(grid_coordinates, transformed_grid_coordinates)
            # transformed_grid_coordinates2 = polywarp_apply(Pi, Qi, transformed_grid_coordinates)
            # plt.scatter(transformed_grid_coordinates2[:, 0], transformed_grid_coordinates2[:, 1], marker='.')
            # plt.scatter(grid_coordinates[:, 0], grid_coordinates[:, 1], marker='.', facecolors='none', edgecolors='r')
       # self.mapping = Mapping2(transformation_type='nonlinear')
        self.mapping.transformation_inverse = PolywarpTransform(params=(Pi,Qi)) # {'P': Pi, 'Q': Qi}
        self.mapping.file = self
        self.mapping.source_name = 'Donor'
        self.mapping.destination_name = 'Acceptor'

    def export_map_file(self):
        warnings.warn('The export_map_file method will be depricated, use export_mapping instead')
        self.export_mapping(filetype='classic')

    def import_mapping_file(self):
        self.mapping = Mapping2(load=self.absoluteFilePath.with_suffix('.mapping'))

    def find_coordinates(self, configuration=None):
        '''
        This function finds and sets the locations of all molecules within the movie's images.

        Main part of the work is done within 'find_molecules_utils.py'(imported as findMols)

        findMols.find_unique_molecules(), will loop through frames of the movie:
        1. creates an average image or maximum projection (for every N frames)
        2. finds the peaks in the new image (wraps around built-in method)
        3. keeps only the unique positions/pixels of peak locations ()
             3B allows to select only 'FRET pairs': intensity peak is seen in both donor/acceptor
        4. Finally, correct for an uncertainty in exact location due to photon shot noise etc. to remove 'functional duplicates'
        5. generate the coordinates array that is compatible with the Molecule() properties:
            [ [array1: donor coordinates],[array2: acceptor coordinates]  ]
            First array has all the coordinates (either found directly or inferred) of the donors:
            donor_coordinates = [[x1,y1],[x2,y2], etc.]
            Similar for acceptor coordinates



        configurations to be set by user:
        (within the find_coordinates section of the configuration file)
        --------------------------------
        channel : choose 'd'/'donor' or 'a'/'acceptor' for using only one of the two channels. alternatively
                  choose 'total_intensity' to use the sum of donor and acceptor signals
                  choose 'FRET pair' to only keep peak intensities found in both channels
                  choose 'both channels' to keep all peaks found in the entire image

        method:  choose 'average_image' or 'maximum_projection_image' to
                         set type of image used to find peak intensities.

        uncertainty_pixels: set number of pixels within which two peak intensities
                            should be considered to correspond to the same molecule

        img_per_N_pixels: use the image type of 'method' for every N frames of the movie (sliding window).
                         If no sliding window is used, the first N frames are used (a 'single window')

        use_sliding_window: type 'True' or 'False' to activate sliding window


        Additional configurations to be set
        (within 'peak_finding' section of configuration file)
        ------------------------------------
        ADD DESCRIPTIONS HERE!!!
        '''

        # --- Refresh configuration ----
        if not configuration:
            self.experiment.import_config_file()
            configuration = self.experiment.configuration['find_coordinates']

        # --- Get settings from configuration file ----
        channels = configuration['channels']
        method = configuration['method']
        peak_finding_configuration = configuration['peak_finding']
        projection_image_type = configuration['projection_image_type']

        sliding_window = configuration['sliding_window']
        use_sliding_window = bool(sliding_window['use_sliding_window'])
        window_size = sliding_window['window_size']
        minimal_point_separation = sliding_window['minimal_point_separation']

        frames_for_peak_finding = configuration['frames_for_peak_finding']
        if frames_for_peak_finding['last_frame'] == 'last':
            frames_for_peak_finding['last_frame'] = self.movie.number_of_frames-1

        # --- set illumination configuration
        #  An integer number for choosing one of the laser lines (the order of it first appeared)
        #  ex. Two laser lines (532 and 640) in Alex mode starting with 532 excitation: 0 for green and 1 for red
        #  None for simple average of the frames regardless of the order of illumination profile.
        illumination = None
        if 'illumination' in configuration:
            if configuration['illumination'] == 'None':
                illumination = None
            else:
                illumination = configuration['illumination']
                print(f'  frames with "{self.movie.illuminations[illumination]}" were chosen for peak finding')


        # --- make the windows
        # (if no sliding windows, just a single window is made to make it compatible with next bit of code) ----
        if use_sliding_window:
            window_start_frames = [i * window_size for i in range(self.number_of_frames // window_size)]
        else:
            window_start_frames = [int(frames_for_peak_finding['first_frame'])]
            window_size = int(frames_for_peak_finding['last_frame']) - frames_for_peak_finding['first_frame'] + 1

        # coordinates = set()
        if method == 'by_channel':
            coordinate_sets = [set() for channel in channels]
        elif method in ('average_channels', 'sum_channels'):
            # SHK: following two lines are not needed for average/sum channels
            # if len(channels) < 2:
            #     raise ValueError('No channels to overlay')
            # END SHK
            coordinate_sets = [set()]

        # coordinates_sets = dict([(channel, set()) for channel in channels])
        # coordinate_sets = [set() for channel in channels]

        # --- Loop over all frames and find unique set of molecules ----
        for window_start_frame in window_start_frames:

            # --- allowed to apply sliding window to either the max projection OR the averages ----
            image = self.movie.make_projection_image(projection_type=projection_image_type, start_frame=window_start_frame,
                                                     number_of_frames=window_size, illumination=illumination)

            # Do we need a separate image?
            # # --- we output the "sum of these images" ----
            # find_coords_img += image

            if method == 'by_channel':
                # coordinates_per_channel = dict([(channel, set()) for channel in channels])
                channel_images = [self.movie.get_channel(image=image, channel=channel) for channel in channels]

            elif method in ('average_channels', 'sum_channels'):
                # Possibly we can move making the overlayed image to the Movie class.
                # TODO: make this usable for any number of channels
                donor_image = self.movie.get_channel(image=image, channel='d')
                # acceptor_image = self.movie.get_channel(image=image, channel='a')
                image_transformed = self.mapping.transform_image(image, direction='Acceptor2Donor')
                acceptor_image_transformed = self.movie.get_channel(image=image_transformed, channel='d')

                if method == 'average_channels':
                    channel_images = [(donor_image + acceptor_image_transformed) / 2]
                elif method == 'sum_channels':
                    channel_images = [(donor_image + acceptor_image_transformed)]
                channels = ['d'] # When number of channels can be > 2 this should probably be the channel with the lowest number

                # TODO: Make this a separate plotting function, possibly in Movie
                # plt.imshow(np.stack([donor_image.astype('uint8'),
                #                      acceptor_image_transformed.astype('uint8'),
                #                      np.zeros((self.movie.height,
                #                                self.movie.width // 2)).astype('uint8')],
                #                     axis=-1))
            else:
                raise ValueError(f'"{method}" is not a valid method.')

            print(f' Finding molecules in {self}')
            for i, channel_image in enumerate(channel_images):
                channel_coordinates = find_peaks(image=channel_image, **peak_finding_configuration)  # .astype(int)))

                # ---- optimize / fine-tune the coordinate positions ----
                coordinate_optimization_functions = \
                    {'coordinates_within_margin': coordinates_within_margin,
                     'coordinates_after_gaussian_fit': coordinates_after_gaussian_fit,
                     'coordinates_without_intensity_at_radius': coordinates_without_intensity_at_radius}
                for f, kwargs in configuration['coordinate_optimization'].items():
                    channel_coordinates = coordinate_optimization_functions[f](channel_coordinates, channel_image, **kwargs)

                channel_coordinates = set_of_tuples_from_array(channel_coordinates)

                coordinate_sets[i].update(channel_coordinates)

        # Check whether points are found
        for coordinate_set in coordinate_sets:
            if len(coordinate_set) == 0:
                self._init_dataset(0)  # SHK: Creating a dummy dataset tp avoid errors in the downstream analysis
                print('no peaks found')
                return

        # --- correct for photon shot noise / stage drift ---
        # Not sure whether to put this in front of combine_coordinate_sets/detect_FRET_pairs or behind [IS: 12-08-2020]
        # I think before, as you would do it either for each window, or for the combined windows.
        # Transforming the coordinate sets for each window will be time consuming and changes the distance_threshold.
        # And you would like to combine the channel sets on the merged coordinates.
        for i in range(len(coordinate_sets)):
            # --- turn into array ---
            coordinate_sets[i] = array_from_set_of_tuples(coordinate_sets[i])

            if use_sliding_window: # Do we actually need to put this if statement here [IS: 31-08-2020]
                                   # If not, in default configuration take minimal_point_separation outside sliding_window
                coordinate_sets[i] = merge_nearby_coordinates(coordinate_sets[i], distance_threshold=minimal_point_separation)

            # Map coordinates to main channel in movie
            # TODO: make this usable for any number of channels

            coordinate_sets[i] = coordinate_sets[i]+self.movie.get_channel_from_name(channels[i]).boundaries[0]
            if channels[i] in ['a', 'acceptor']:
            # if i > 0: #i.e. if channel is not main channel # this didn't work when selecting only the acceptor channel
                # Maybe we can do this earlier, right after point detection, then we need only a single coordinate_set
                coordinate_sets[i] = self.mapping.transform_coordinates(coordinate_sets[i],
                                                                        direction='Acceptor2Donor')

        # TODO: make this usable for any number of channels
        if len(coordinate_sets) == 1:
            coordinates = coordinate_sets[0]
        elif len(coordinate_sets) > 1:
            raise NotImplementedError('Assessing found coordinates in multiple channels does not work properly yet')
            # TODO: Make this function.
            #  This can easily be done by creating a cKDtree for each coordinate set and
            #  by finding the points close to each other
            coordinates = combine_coordinate_sets(coordinate_sets, method='and')  # the old detect_FRET_pairs

        # TODO: make this usable for more than two channels
        # TODO: Use set_coordinates_of_channel
        coordinates_in_main_channel = coordinates
        coordinates_list = [coordinates]
        for i in range(self.number_of_channels)[1:]: # This for loop will only be useful once we make this usable for more than two channels
            if self.number_of_channels > 2:
                raise NotImplementedError()
            coordinates_in_other_channel = self.mapping.transform_coordinates(coordinates_in_main_channel, direction='Donor2Acceptor')
            coordinates_list.append(coordinates_in_other_channel)

        coordinates = np.hstack(coordinates_list)

        coordinates_selections = [coordinates_within_margin_selection(coordinates, bounds=self.movie.channels[i].boundaries)
                                  for i, coordinates in enumerate(coordinates_list)]
        selection = np.vstack(coordinates_selections).all(axis=0)
        coordinates = coordinates[selection]

        coordinates = coordinates.reshape((-1, 2))

        sys.stdout.write('\r')
        print(f'   {coordinates.shape[0]} molecules found')

        # should also have incorporated check coordinatesDA_within_margin from MD_check_boundaries
        # --- finally, we set the coordinates of the molecules ---
        # self.coordinates = coordinates

        peaks = xr.DataArray(coordinates, dims=("peak", 'dimension'),
                     coords={'peak': range(len(coordinates)), 'dimension': ['x', 'y']}, name='coordinates')

        coordinates = split_dimension(peaks, 'peak', ('molecule', 'channel'), (-1, 2)).reset_index('molecule', drop=True)
        # file = str(self.relativeFilePath)
        # #coordinates = split_dimension(coordinates, 'molecule', ('molecule_in_file', 'file'), (-1, 1), (-1, [file]), to='multiindex')
        # coordinates = coordinates.reset_index('molecule').rename(molecule_='molecule_in_file')
        # self.experiment.dataset.drop_sel(file=str(self.relativeFilePath), errors='ignore')

        # if len(coordinates) !=0:

        self.coordinates = coordinates

        # self.molecules.export_pks_file(self.relativeFilePath.with_suffix('.pks'))

    @property
    def coordinates(self):
        return self.dataset.coordinates

    @coordinates.setter
    def coordinates(self, coordinates):
        # Reset current .nc file
        self._init_dataset(len(coordinates.molecule))

        coordinates.to_netcdf(self.relativeFilePath.with_suffix('.nc'), engine='h5netcdf', mode='a')
        self.extract_background()

    def extract_background(self, configuration=None):
        sys.stdout.write(f' Calculating background in {self}')
        # --- Refresh configuration ----
        if not configuration:
            self.experiment.import_config_file()
            configuration = self.experiment.configuration['background']

        # --- get settings from configuration file ----
        if 'frames_for_background' in configuration.keys():
            frames_for_background = configuration['frames_for_background']
            if frames_for_background['last_frame'] == 'last':
                frames_for_background['last_frame'] = self.movie.number_of_frames-1
            frames_for_background = [int(frames_for_background['first_frame']),
                                     int(frames_for_background['last_frame']) - frames_for_background['first_frame'] + 1]
        else:
            print('configurations for the background subtraction is not found. All the frames will be used.')
            frames_for_background = [0, self.movie.number_of_frames-1]

        # --- get the averaged images for background extraction per illumination profile
        # todo: In the following code, let's use 'illumination' instead of 'illuminations', similar to channel.
        if self.movie.illumination_arrangement is not None:
            image_for_background = [None] * len(self.movie.illumination_arrangement)
            illuminations_to_use = self.movie.illumination_arrangement
        else:
            image_for_background = [None]
            illuminations_to_use = [None]
        for illumination_id, illumination in enumerate(illuminations_to_use):
            image_for_background[illumination_id] = \
                        self.movie.make_projection_image(projection_type='average',
                                                         start_frame=frames_for_background[0],
                                                         number_of_frames=frames_for_background[1],
                                                         illumination=illumination)

        # START of original code
        # background_list = []
        # for i, channel in enumerate(self.movie.channels):
        #     channel_image = self.movie.get_channel(image_for_background[illumination_id], i)
        #     channel_coordinates = self.coordinates_from_channel(i).values-self.movie.channels[i].vertices[0]
        #     background_list.append(extract_background(channel_image, channel_coordinates, method=configuration['method']))
        #
        # background = xr.DataArray(np.vstack(background_list).T, dims=['molecule','channel'], name='background')

        # END of original, START modified
        background_list = []
        for illumination_id, illumination in enumerate(illuminations_to_use):
            tmp_background_list = []
            for i, channel in enumerate(self.movie.channels):
                channel_image = self.movie.get_channel(image_for_background[illumination_id], i)
                channel_coordinates = self.coordinates_from_channel(i).values-self.movie.channels[i].vertices[0]
                tmp_background_list.append(extract_background(channel_image, channel_coordinates, method=configuration['method']))
            background_list.append(np.vstack(tmp_background_list).T)
        background = xr.DataArray(background_list, dims=['illuminations', 'molecule', 'channel'], name='background')
        # END modified

        background.to_netcdf(self.relativeFilePath.with_suffix('.nc'), engine='h5netcdf', mode='a')
        sys.stdout.write(f'\r   background calculated {self}\n')

    def import_excel_file(self, filename=None):
        if filename is None:
            filename = f'{self.relativeFilePath}_steps_data.xlsx'
        try:
            steps_data = pd.read_excel(filename, index_col=[0,1],
                                       dtype={'kon':np.str})       # reads from the 1st excel sheet of the file

            print(f'imported steps data from excel file for {self.name}')
        except FileNotFoundError:
            print(f'No saved analysis for {self.name} as {filename}')
            return
        molecules = steps_data.index.unique(0)
        indices = [int(m.split()[-1]) for m in molecules]
        for mol in self.molecules:
            if mol.index + 1 not in indices:
                continue
            mol.steps = steps_data.loc[f'mol {mol.index + 1}']
            # if saved steps are found for molecule it is assumed selected
            # mol.isSelected = True
            if 'kon' in mol.steps.columns:
                k = [int(i) for i in mol.steps.kon[0]]
                mol.kon_boolean = np.array(k).astype(bool).reshape((4,3))
        return steps_data

    def import_selected(self):
        '''
        Imports the selected molecules stored in {filename}_selected_molecules.txt
        '''
        pass
        # try:
        #     filename = f'{self.relativeFilePath}_selected_molecules.txt'
        #     selected = np.atleast_1d(np.loadtxt(filename, dtype=int))
        # except FileNotFoundError:
        #     return
        # # print(selected, type(selected))
        # for i in list(selected):
        #     self.molecules[i-1].isSelected = True

    def extract_traces(self, configuration=None):
        if self.number_of_molecules == 0:
            print('   no traces available!!')
            return

        # Refresh configuration
        self.experiment.import_config_file()

        if self.movie is None: raise FileNotFoundError('No movie file was found')

        print(f'  Extracting traces in {self}')

        if configuration is None: configuration = self.experiment.configuration['trace_extraction']
        channel = configuration['channel']  # Default was 'all'
        mask_size = configuration['mask_size']  # Default was 11
        neighbourhood_size = configuration['neighbourhood_size']  # Default was 11
        subtract_background = configuration['subtract_background']

        if mask_size == 'TIR-T' or mask_size == 'TIR-V':
            mask_size = 1.291
        elif mask_size == 'TIR-S 1.5x 2x2':
            mask_size = 0.8
        elif mask_size == 'TIR-S 1x 2x2':
            mask_size = 0.55

        if subtract_background:
            background = self.background.stack(peak=('molecule', 'channel'))
        else:
            background = None

        coordinates = self.coordinates.stack(peak=('molecule', 'channel')).T

        if self.movie.illumination_arrangement is not None:
            number_illumination = len(self.movie.illumination_arrangement)
        else:
            number_illumination = 1

        # todo: we do not have to pass "number_illumination" as an argument, because "self.movie" is also passed to the function.
        traces = extract_traces(self.movie, coordinates.values, background=background.values, channel=channel,
                                mask_size=mask_size, neighbourhood_size=neighbourhood_size,
                                number_illumination=number_illumination)

        intensity = xr.DataArray(traces, dims=['peak', 'frame'], name='intensity')\
            .assign_coords({'peak': coordinates.peak.to_index()})\
            .unstack('peak').reset_index('molecule', drop=True)\
            .assign_coords(self.coordinates.sel(dimension='x', drop=True).coords)\
            .transpose('molecule', 'channel', 'frame')
        # number_of_molecules = len(traces) // self.number_of_channels
        # traces = traces.reshape((number_of_molecules, self.number_of_channels, self.movie.number_of_frames)).swapaxes(0, 1)

        if hasattr(self.movie, 'time'):
            intensity = intensity.assign_coords(time=self.movie.time)

        if self.movie.illumination is not None:
            intensity = intensity.assign_coords(illumination=self.movie.illumination)

        intensity.to_netcdf(self.relativeFilePath.with_suffix('.nc'), engine='h5netcdf', mode='a')

        self.calculate_FRET()

        #self.export_traces_file()

    def calculate_FRET(self):
        # TODO: Make suitable for mutliple colours
        # TODO: Implement corrections
        donor = self.intensity.sel(channel=0, drop=True)
        acceptor = self.intensity.sel(channel=1, drop=True)
        FRET = acceptor/(donor+acceptor)
        FRET.name = 'FRET'
        FRET.to_netcdf(self.relativeFilePath.with_suffix('.nc'), engine='h5netcdf', mode='a')


    def import_pks_file(self):
        peaks = import_pks_file(self.relativeFilePath.with_suffix('.pks'))
        peaks = split_dimension(peaks, 'peak', ('molecule', 'channel'), (-1, 2)).reset_index('molecule', drop=True)
        # peaks = split_dimension(peaks, 'molecule', ('molecule_in_file', 'file'), (-1, 1), (-1, [file]), to='multiindex')

        if not self.relativeFilePath.with_suffix('.nc').is_file():
            self._init_dataset(len(peaks.molecule))

        coordinates = peaks.sel(parameter=['x', 'y']).rename(parameter='dimension')
        background = peaks.sel(parameter='background', drop=True)

        xr.Dataset({'coordinates': coordinates, 'background': background})\
            .to_netcdf(self.relativeFilePath.with_suffix('.nc'), engine='h5netcdf', mode='a')

    def export_pks_file(self):
        peaks = xr.merge([self.coordinates.to_dataset('dimension'), self.background.to_dataset()])\
            .stack(peaks=('molecule', 'channel')).to_array(dim='parameter').T
        export_pks_file(peaks, self.relativeFilePath.with_suffix('.pks'))
        self.extensions.add('.pks')

    def import_traces_file(self):
        traces = import_traces_file(self.relativeFilePath.with_suffix('.traces'))
        intensity = split_dimension(traces, 'trace', ('molecule', 'channel'), (-1, 2))\
            .reset_index(['molecule','frame'], drop=True)

        if not self.relativeFilePath.with_suffix('.nc').is_file():
            self._init_dataset(len(intensity.molecule))

        xr.Dataset({'intensity': intensity}).to_netcdf(self.relativeFilePath.with_suffix('.nc'), engine='h5netcdf', mode='a')

    def export_traces_file(self):
        traces = self.intensity.stack(trace=('molecule', 'channel')).T
        export_traces_file(traces, self.relativeFilePath.with_suffix('.traces'))
        self.extensions.add('.traces')


    # def addMolecule(self):
    #     index = len(self.molecules) # this is the molecule number
    #     self.molecules.append(Molecule(self))
    #     self.molecules[-1].index = index

    def histogram(self, axis=None, bins=100, parameter='E', molecule_averaging=False,
                  makeFit=False, export=False, **kwargs):
        histogram(self.molecules, axis=axis, bins=bins, parameter=parameter, molecule_averaging=molecule_averaging, makeFit=makeFit, collection_name=self, **kwargs)
        if export:
            plt.savefig(self.absoluteFilePath.with_name(f'{self.name}_{parameter}_histogram').with_suffix('.png'))



    def savetoExcel(self, filename=None, save=True):
        if filename is None:
            filename = f'{self.relativeFilePath}_steps_data.xlsx'

        # Find the molecules for which steps were selected
        molecules_with_data = [mol for mol in self.molecules if mol.steps is not None]


        # Concatenate all steps dataframes that are not None
        mol_data = [mol.steps for mol in molecules_with_data]
        if not mol_data:
            print(f'no data to save for {self.name}')
            return
        keys = [f'mol {mol.index + 1}' for mol in molecules_with_data]

        steps_data = pd.concat(mol_data, keys=keys, sort=False)
        # drop duplicate columns
        steps_data = steps_data.loc[:,~steps_data.columns.duplicated()]
        if save:
            print("data saved in: " + filename)
            writer = pd.ExcelWriter(filename)
            steps_data.to_excel(writer, self.name)
            writer.save()
        return steps_data

    def autoThreshold(self, trace_name, threshold=100, max_steps=20,
                      only_selected=False, kon_str='000000000'):
        nam = trace_name
        for mol in self.molecules:

            trace = mol.I(0)*int((nam == 'green')) + \
                    mol.I(1)*int((nam == 'red')) +\
                     mol.E()*int((nam == 'E'))  # Here no offset corrections are applied yet

            d = mol.find_steps(trace)
            frames = d['frames']
            times = frames*self.exposure_time
            times = np.sort(times)
            mol.steps = pd.DataFrame({'time': times, 'trace': nam,
                                  'state': 1, 'method': 'thres',
                                'thres': threshold, 'kon': kon_str})
        filename = self.name+'_steps_data.xlsx'
        data = self.savetoExcel(filename)
        return data

    def select(self, figure=None):
        # iasonas: I think this function is not needed anymore
        plt.ion()
        for index, molecule in enumerate(self.molecules):
            molecule.plot(figure=figure)
            plt.title('Molecule ' + str(index), y=-0.01)
            plt.show()
            plt.pause(0.001)
            print('Molecule ' + str(index))
            input("Press enter to continue")

    def perform_mapping(self, configuration = None):
        # Refresh configuration
        if not configuration:
            self.experiment.import_config_file()

        image = self.average_image
        if configuration is None:
            configuration = self.experiment.configuration['mapping']

        transformation_type = configuration['transformation_type']
        print(transformation_type)
        method = configuration['method']

        donor_image = self.movie.get_channel(image=image, channel='d')
        acceptor_image = self.movie.get_channel(image=image, channel='a')
        donor_coordinates = find_peaks(image=donor_image,
                                       **configuration['peak_finding']['donor'])
        if donor_coordinates.size == 0: #should throw a error message to warm no acceptor molecules found
            print('No donor molecules found')
        acceptor_coordinates = find_peaks(image=acceptor_image,
                                          **configuration['peak_finding']['acceptor'])
        if acceptor_coordinates.size == 0: #should throw a error message to warm no acceptor molecules found
            print('No acceptor molecules found')
        acceptor_coordinates = transform(acceptor_coordinates, translation=[image.shape[1]//2, 0])
        # print(acceptor_coordinates.shape, donor_coordinates.shape)
        print(f'Donor: {donor_coordinates.shape[0]}, Acceptor: {acceptor_coordinates.shape[0]}')
        coordinates = np.append(donor_coordinates, acceptor_coordinates, axis=0)

        # coordinate_optimization_functions = \
        #     {'coordinates_within_margin': coordinates_within_margin,
        #      'coordinates_after_gaussian_fit': coordinates_after_gaussian_fit,
        #      'coordinates_without_intensity_at_radius': coordinates_without_intensity_at_radius}
        #
        # for f, kwargs in configuration['coordinate_optimization'].items():
        #     coordinates = coordinate_optimization_functions[f](coordinates, image, **kwargs)

        if 'coordinates_after_gaussian_fit' in configuration['coordinate_optimization']:
            gaussian_width = configuration['coordinate_optimization']['coordinates_after_gaussian_fit']['gaussian_width']
            coordinates = coordinates_after_gaussian_fit(coordinates, image, gaussian_width)

        if 'coordinates_without_intensity_at_radius' in configuration['coordinate_optimization']:
            coordinates = coordinates_without_intensity_at_radius(coordinates, image,
                                                                  **configuration['coordinate_optimization']['coordinates_without_intensity_at_radius'])
                                                                  # radius=4,
                                                                  # cutoff=np.median(image),
                                                                  # fraction_of_peak_max=0.35) # was 0.25 in IDL code

        if 'coordinates_within_margin' in configuration['coordinate_optimization']:
            margin = configuration['coordinate_optimization']['coordinates_within_margin']['margin']
        else:
            margin = 0

        donor_coordinates = coordinates_within_margin(coordinates,
                                                      bounds=self.movie.get_channel_from_name('d').boundaries, margin=margin)
        acceptor_coordinates = coordinates_within_margin(coordinates,
                                                         bounds=self.movie.get_channel_from_name('a').boundaries, margin=margin)

        # TODO: put overlapping coordinates in file.coordinates for mapping file
        # Possibly do this with mapping.nearest_neighbour match
        # self.coordinates = np.hstack([donor_coordinates, acceptor_coordinates]).reshape((-1, 2))

        if ('initial_translation' in configuration) and (configuration['initial_translation'] == 'width/2'):
            # initial_transformation = {'translation': [image.shape[0] // 2, 0]}
            initial_transformation = {'translation': [image.shape[1] // 2, 0]}
        else:
            if configuration['initial_translation'][0] == '[':  # remove brackets
                arr = [float(x) for x in configuration['initial_translation'][1:-1].split(' ')]
                initial_transformation = {'translation': arr}
            else:
                arr = [float(x) for x in configuration['initial_translation'].split(' ')]
                initial_transformation = {'translation': configuration['initial_translation']}

        # Obtain specific mapping parameters from configuration file
        additional_mapping_parameters = {key: configuration[key]
                                         for key in (configuration.keys() and {'distance_threshold'})}

        self.mapping = Mapping2(source_name='Donor',
                                source=donor_coordinates,
                                destination_name='Acceptor',
                                destination=acceptor_coordinates,
                                method=method,
                                transformation_type=transformation_type,
                                initial_transformation=initial_transformation)
        self.mapping.perform_mapping(**additional_mapping_parameters)
        self.mapping.file = self

        # self.export_mapping(filetype='classic')
        self.export_mapping()

        self.use_mapping_for_all_files()

    def copy_coordinates_to_selected_files(self):
        for file in self.experiment.selectedFiles:
            if file is not self:
                file._init_dataset(len(self.molecule))
                self.coordinates.to_netcdf(file.relativeFilePath.with_suffix('.nc'), engine='h5netcdf')

    def use_mapping_for_all_files(self):
        print(f"\n{File} used as mapping")
        self.is_mapping_file = True
        #mapping = self.movie.use_for_mapping()
        for file in self.experiment.files:
            if file is not self:
                file.mapping = self.mapping
                file.is_mapping_file = False

    def show_image(self, image_type='default', figure=None, **kwargs):
        # Refresh configuration
        if image_type == 'default':
            self.experiment.import_config_file()
            image_type = self.experiment.configuration['show_movie']['image']

        if figure is None:
            figure = plt.figure()
        axis = figure.gca()

        # Choose method to plot
        if image_type == 'average_image':
            image = self.average_image
            axis.set_title('Average image')
        elif image_type == 'maximum_image':
            image = self.maximum_projection_image
            axis.set_title('Maximum projection')

        image_handle = axis.imshow(image)

        # process keyword arguments
        colorscale = list(image_handle.get_clim())
        if 'vmin' in kwargs:
            colorscale[0] = kwargs['vmin']
        if 'vmax' in kwargs:
            colorscale[1] = kwargs['vmax']
        image_handle.set_clim(colorscale)

    def show_average_image(self, figure=None):
        self.show_image(image_type='average_image', figure=figure)

    def show_coordinates(self, figure=None, annotate=None, **kwargs):
        # Refresh configuration
        self.experiment.import_config_file()

        if not figure:
            figure = plt.figure()

        if annotate is None:
            annotate = self.experiment.configuration['show_movie']['annotate']

        if self.coordinates is not None:
            axis = figure.gca()
            coordinates = self.coordinates.stack({'peak': ('molecule', 'channel')}).T.values
            sc_coordinates = axis.scatter(coordinates[:, 0], coordinates[:, 1], facecolors='none', edgecolors='red', **kwargs)
            # marker='o'

            selected_coordinates = self.coordinates[self.selected].stack({'peak': ('molecule', 'channel')}).T.values
            axis.scatter(selected_coordinates[:, 0], selected_coordinates[:, 1], facecolors='none', edgecolors='green', **kwargs)

            if annotate:
                annotation = axis.annotate("", xy=(0, 1.03), xycoords=axis.transAxes) # x in data units, y in axes fraction
                annotation.set_visible(False)

                #molecule_indices = np.repeat(np.arange(0, self.number_of_molecules), self.number_of_channels)
                molecule_indices = np.repeat(self.molecule_in_file.values, self.number_of_channels)
                # sequence_indices = np.repeat(self.sequence_indices, self.number_of_channels)

                def update_annotation(ind):
                    # print(ind)

                    # text = "Molecule number: {} \nSequence: {}".format(" ".join([str(indices[ind["ind"][0]])]),
                    #                        " ".join([str(sequences[ind["ind"][0]].decode('UTF-8'))]))
                    plot_index = ind["ind"]
                    molecule_index = molecule_indices[plot_index]
                    text = f'Molecule number: {", ".join(map(str, molecule_index))}'

                    if hasattr(self, 'sequences'):
                        sequence_names = [str(self.sequence_name[index]) for index in molecule_index]
                        sequences = [str(self.sequence[index]) for index in molecule_index]
                        text += f'\nSequence name: {", ".join(sequence_names)}'
                        text += f'\nSequence: {", ".join(sequences)}'

                    annotation.set_text(text)

                def hover(event):
                    vis = annotation.get_visible()
                    if event.inaxes == axis:
                        cont, ind = sc_coordinates.contains(event)
                        if cont:
                            update_annotation(ind)
                            annotation.set_visible(True)
                            figure.canvas.draw_idle()
                        else:
                            if vis:
                                annotation.set_visible(False)
                                figure.canvas.draw_idle()

                figure.canvas.mpl_connect("motion_notify_event", hover)

            plt.show()

    def show_coordinates_in_image(self, figure=None, **kwargs):
        if not figure:
            figure = plt.figure()

        self.show_image(figure=figure, **kwargs)
        self.show_coordinates(figure=figure)
        # plt.savefig(self.writepath.joinpath(self.name + '_ave_circles.png'), dpi=600)


def import_pks_file(pks_filepath):
    pks_filepath = Path(pks_filepath)
    data = np.atleast_2d(np.genfromtxt(pks_filepath)[:,1:])
    if data.shape[1] == 2:
        data = np.hstack([data, np.zeros((len(data),1))])
    return xr.DataArray(data, dims=("peak",'parameter'),
                        coords={'peak': range(len(data)), 'parameter': ['x', 'y', 'background']})


def export_pks_file(peaks, pks_filepath):
    pks_filepath = Path(pks_filepath)
    with pks_filepath.open('w') as pks_file:
        for i, (x, y, background) in enumerate(peaks.values):
            # outfile.write(' {0:4.0f} {1:4.4f} {2:4.4f} {3:4.4f} {4:4.4f} \n'.format(i, coordinate[0], coordinate[1], 0, 0, width4=4, width6=6))
            # pks_file.write('{0:4.0f} {1:4.4f} {2:4.4f} \n'.format(i + 1, coordinate[0], coordinate[1]))
            pks_file.write('{0:4.0f} {1:4.4f} {2:4.4f} {3:4.4f}\n'.format(i + 1, x, y, background))


def import_traces_file(traces_filepath):
    traces_filepath = Path(traces_filepath)
    with traces_filepath.open('r') as traces_file:
        number_of_frames = np.fromfile(traces_file, dtype=np.int32, count=1).item()
        number_of_traces = np.fromfile(traces_file, dtype=np.int16, count=1).item()
        # number_of_molecules = number_of_traces // number_of_channels
        raw_data = np.fromfile(traces_file, dtype=np.int16, count=number_of_frames * number_of_traces)
    # traces = np.reshape(rawData.ravel(),
    #                         (number_of_channels, number_of_molecules, number_of_frames),
    #                         order='F')  # 3d array of trace # 2d array of traces
    traces = np.reshape(raw_data.ravel(), (number_of_traces, number_of_frames), order='F') # 2d array of traces
    traces = xr.DataArray(traces, dims=("trace", "frame"), coords=(range(number_of_traces), range(number_of_frames)))
    return traces


def export_traces_file(traces, traces_filepath):
    traces_filepath = Path(traces_filepath)
    with traces_filepath.open('w') as traces_file:
        # Number of frames
        np.array([len(traces.frame)], dtype=np.int32).tofile(traces_file)
        # Number of traces
        np.array([len(traces.trace)], dtype=np.int16).tofile(traces_file)
        traces.values.T.astype(np.int16).tofile(traces_file)

def split_dimension(data_array, old_dim, new_dims, new_dims_shape=None, new_dims_coords=None, to='dimensions'):
    all_dims = list(data_array.dims)
    old_dim_index = all_dims.index(old_dim)
    all_dims[old_dim_index:old_dim_index + 1] = new_dims
    new_dims_shape = np.array(new_dims_shape)
    if sum(new_dims_shape == -1) == 1:
        fixed_dim_prod = np.prod(new_dims_shape[new_dims_shape!=-1])
        old_len = data_array.shape[data_array.dims==old_dim]
        if old_len % fixed_dim_prod != 0:
            raise ValueError('Incorrect dimension shape')
        new_dims_shape[new_dims_shape == -1] = old_len // fixed_dim_prod
    elif sum(new_dims_shape == -1) > 1:
        raise ValueError

    if new_dims_coords is None:
        new_dims_coords = [-1]*len(new_dims_shape)
    new_dims_coords = (range(new_dims_shape[i]) if new_dims_coord == -1 else new_dims_coord
                       for i, new_dims_coord in enumerate(new_dims_coords))

    new_index = pd.MultiIndex.from_product(new_dims_coords, names=new_dims)
    data_array = data_array.assign_coords(**{old_dim: new_index})

    if to == 'dimensions':
        return data_array.unstack(old_dim).transpose(*all_dims)
    elif to == 'multiindex':
        return data_array
    else:
        raise ValueError
