# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 15:24:46 2018

@author: ivoseverins
"""
# Use the following lines on Mac
# from sys import platform
# if platform == "darwin":
#     from matplotlib import use
#     use('WXAgg')

import os  # Miscellaneous operating system interfaces - to be able to switch from Mac to Windows
from pathlib import Path  # For efficient path manipulation
import yaml
import numpy as np
import pandas as pd
# import matplotlib
# matplotlib.use('WXAgg')

import matplotlib.pyplot as plt  # Provides a MATLAB-like plotting framework
import xarray as xr

from trace_analysis.file import File
# from trace_analysis.molecule import Molecules
from trace_analysis.collection import Collection
from trace_analysis.plotting import histogram
# from trace_analysis.plugin_manager import PluginManager
# from trace_analysis.plugin_manager import PluginMetaClass
from trace_analysis.plugin_manager import plugins

import re  # Regular expressions
import warnings
from nd2reader import ND2Reader


# import matplotlib.pyplot as plt #Provides a MATLAB-like plotting framework
# import itertools #Functions creating iterators for efficient looping
# np.seterr(divide='ignore', invalid='ignore')
# import pandas as pd
# from threshold_analysis_v2 import stepfinder
# import pickle


@plugins
class Experiment:
    """ Main experiment class

    Class containing all the files in an experiment.
    In fact it can contain any collection of files.

    .. warning:: Only works with one or two channels.

    Attributes
    ----------
    name : str
        Experiment name based on the name of the main folder
    main_path : str
        Absolute path to the main experiment folder
    files : list of :obj:`File`
        Files
    import_all : bool
        If true, then all files in the main folder are automatically imported. \n
        If false, then files are detected, but not imported.
    """
    # TODO: Add presets for specific microscopes
    def __init__(self, main_path, channels=['g', 'r'], import_all=True):
        """Init method for the Experiment class

        Loads config file if it locates one in the main directory, otherwise it exports the default config file to the main directory.
        Scans all directory in the main directory recursively and imports all found files (if import_all is set to `True`).

        Parameters
        ----------
        main_path : str
            Absolute path to the main experiment folder
        channels : list of str
            Channels used in the experiment
        import_all : bool
            If true, then all files in the main folder are automatically imported. \n
            If false, then files are detected, but not imported.
        """

        self.name = os.path.basename(main_path)
        self.main_path = Path(main_path).absolute()
        self.files = Collection()
        self.import_all = import_all

        self._channels = np.atleast_1d(np.array(channels))
        self._number_of_channels = len(channels)
        self._pairs = [[c1, c2] for i1, c1 in enumerate(channels) for i2, c2 in enumerate(channels) if i2 > i1]

        # Load custom config file or otherwise load the default config file
        if self.main_path.joinpath('config.yml').is_file():
            self.import_config_file()
        else:
            with Path(__file__).with_name('default_configuration.yml').open('r') as yml_file:
                self.configuration = yaml.load(yml_file, Loader=yaml.SafeLoader)
            try:
                self.export_config_file()
            except:
                FileNotFoundError
                pass

        os.chdir(main_path)

        # file_paths = self.find_file_paths()
        # self.add_files(file_paths, test_duplicates=False)

        self.add_files(self.main_path, test_duplicates=False)

        # Find mapping file
        for file in self.files:
            if file.mapping is not None:
                file.use_mapping_for_all_files()
                break

        print('\nInitialize experiment: \n' + str(self.main_path))

    def __getstate__(self):
        d = self.__dict__.copy()
        # d.pop('files')
        d['files'] = []
        d['_do_not_update'] = None # This is for parallelization in Collection
        return d

    def __setstate__(self, dict):
        self.__dict__.update(dict)

    def __repr__(self):
        return (f'{self.__class__.__name__}({self.name})')

    @property
    def channels(self):
        """list of str : Channels used in the experiment.

        Setting the channels will automatically update pairs.
        """
        return self._channels

    @channels.setter
    def channels(self, channels):
        self._channels = np.atleast_1d(np.array(channels))
        self._number_of_channels = len(channels)
        self._pairs = [[c1, c2] for i1, c1 in enumerate(channels) for i2, c2 in enumerate(channels) if i2 > i1]

    @property
    def number_of_channels(self):
        """int : Number of channels used in the experiment (read-only)"""
        return self._number_of_channels

    @property
    def pairs(self):
        """list of list of str : List of channel pairs"""
        return self._pairs

    # @property
    # def molecules(self):
    #     """list of Molecule : List of all molecules in the experiment"""
    #     return Molecules.sum([file.molecules for file in self.files])

    @property
    def selectedFiles(self):
        """list of File : List of selected files"""
        return [file for file in self.files if file.isSelected]

    # @property
    # def selectedMoleculesInSelectedFiles(self):
    #     """list of Molecule : List of selected molecules in selected files"""
    #     return [molecule for file in self.selectedFiles for molecule in file.selectedMolecules]

    # @property
    # def selectedMoleculesInAllFiles(self):
    #     """list of Molecule : List of selected molecules in all files"""
    #     return [molecule for file in self.files for molecule in file.selectedMolecules]

    @property
    def mapping_file(self):
        for file in self.files:
            if file.is_mapping_file:
                return file

    @property
    def file_paths(self):
        return [file.relativeFilePath for file in self.files]

    @property
    def nc_file_paths(self):
        return [file.relativeFilePath.with_suffix('.nc') for file in self.files if '.nc' in file.extensions]

    def import_config_file(self):
        """Import configuration file from main folder into the configuration property."""
        with self.main_path.joinpath('config.yml').open('r') as yml_file:
            self.configuration = yaml.load(yml_file, Loader=yaml.SafeLoader)

    def export_config_file(self):
        """Export from the configuration property into the configuration file in main folder"""
        with self.main_path.joinpath('config.yml').open('w') as yml_file:
            yaml.dump(self.configuration, yml_file, sort_keys=False)

    def find_file_paths(self):
        """Find unique files in all subfolders and add them to the experiment

        Get all files in all subfolders of the main_path and remove their suffix (extensions), and add them to the experiment.

        Note
        ----
        Non-relevant files are excluded e.g. files with underscores or 'Analysis' in their name, or files with dat, db,
        ini, py and yml extensions.

        Note
        ----
        Since sifx files made using spooling are all called 'Spooled files' the parent folder is used as file instead of the sifx file

        """

        file_paths = [p.relative_to(self.main_path).with_suffix('') for p in self.main_path.glob('**/*')
                      if (
                          # Use only files
                              p.is_file() &
                              # Exclude stings in filename
                              all(name not in p.with_suffix('').name for name in
                                  self.configuration['files']['excluded_names']) &
                              # Exclude strings in path
                              all(path not in str(p.relative_to(self.main_path).parent) for path in
                                  self.configuration['files']['excluded_paths']) &
                              # Exclude hidden folders
                              ('.' not in [s[0] for s in p.parts]) &
                              # Exclude file extensions
                              (p.suffix[1:] not in self.configuration['files']['excluded_extensions'])
                      )
                      ]

        for i, file_path in enumerate(file_paths):
            if (file_path.name == 'Spooled files'):
                file_path[i] = file_path[i].parent

        unique_file_paths = np.unique(file_paths)

        return unique_file_paths

    def add_files(self, paths, test_duplicates=True):
        """Find unique files in all subfolders and add them to the experiment

        Get all files in all subfolders of the main_path and remove their suffix (extensions), and add them to the experiment.

        Note
        ----
        Non-relevant files are excluded e.g. files with underscores or 'Analysis' in their name, or files with dat, db,
        ini, py and yml extensions.

        Note
        ----
        Since sifx files made using spooling are all called 'Spooled files' the parent folder is used as file instead of the sifx file

        """

        if isinstance(paths, str) or isinstance(paths, Path):
            paths = paths.glob('**/*')

        # TODO: Merge with method find_file_paths
        file_paths_and_extensions = \
            [[p.relative_to(self.main_path).with_suffix(''), p.suffix]
             for p in paths
             if (
                 # Use only files
                     p.is_file() &
                     # Exclude stings in filename
                     all(name not in p.with_suffix('').name for name in
                         self.configuration['files']['excluded_names']) &
                     # Exclude strings in path
                     all(path not in str(p.relative_to(self.main_path).parent) for path in
                         self.configuration['files']['excluded_paths']) &
                     # Exclude hidden folders
                     ('.' not in [s[0] for s in p.parts]) &
                     # Exclude file extensions
                     (p.suffix[1:] not in self.configuration['files']['excluded_extensions'])
             )
             ]

        # TODO: Test spooled file and nd2 file import
        new_file_paths_and_extensions = []
        for i, (file_path, extensions) in enumerate(file_paths_and_extensions):
            if (file_path.name == 'Spooled files'):
                new_file_paths_and_extensions.append([file_path.parent, extensions])
                # file_paths_and_extensions[i, 0] = file_paths_and_extensions[i, 0].parent
            elif '.nd' in extensions and not 'fov' in file_path:
                fov_info = self.get_fov_from_nd2(file_path)
                if fov_info['number_of_fov'] > 1:  # if the file is nd2 with multiple field of views
                    for fov_id in range(fov_info['number_of_fov']):
                        new_path = Path(str(file_path) + f'_fov{fov_id:03d}')
                        new_file_paths_and_extensions.append([new_path, extensions])
                        # fov_info['fov_chosen'] = fov_id
                        # new_file = File(new_path, self, fov_info=fov_info.copy())
                        # if new_file.extensions:
                        #     self.files.append(new_file)
                else:
                    new_file_paths_and_extensions.append([file_path, extensions])
            else:
                new_file_paths_and_extensions.append([file_path, extensions])
        file_paths_and_extensions = new_file_paths_and_extensions

        file_paths_and_extensions = np.array(file_paths_and_extensions)

        file_paths_and_extensions = file_paths_and_extensions[file_paths_and_extensions[:, 0].argsort()]
        unique_file_paths, indices = np.unique(file_paths_and_extensions[:, 0], return_index=True)
        extensions_per_filepath = np.split(file_paths_and_extensions[:, 1], indices[1:])

        for file_path, extensions in zip(unique_file_paths, extensions_per_filepath):
            if not test_duplicates or (file_path.absolute().relative_to(self.main_path) not in self.file_paths):
                self.files.append(File(file_path, extensions, self))
            else:
                i = self.file_paths.find(file_path.absolute().relative_to(self.main_path))
                self.files[i].add_extensions(extensions)

        # nd2_file = list(self.main_path.glob(str(relativeFilePath) + '.nd2'))
        # fov_info = {'number_of_fov': 1}  # fov=Field of View
        # if nd2_file:  # check if the nd2 file has multiple fov data,
        #     fov_info = self.get_fov_from_nd2(nd2_file[0])
        #
        # if fov_info['number_of_fov'] > 1:  # if the file is nd2 with multiple field of views
        #     for fov_id in range(fov_info['number_of_fov']):
        #         new_path = Path(str(relativeFilePath) + f'_fov{fov_id:03d}')
        #         fov_info['fov_chosen'] = fov_id
        #         new_file = File(new_path, self, fov_info=fov_info.copy())
        #         if new_file.extensions:
        #             self.files.append(new_file)



    # def add_files(self, file_paths, test_duplicates=True):
    #     for file_path in file_paths:
    #         self.add_file(file_path, test_duplicates)
    #
    #     for file in self.files:
    #         if file.mapping is not None:
    #             file.use_mapping_for_all_files()
    #             break

    # def add_file(self, file_path, test_duplicates=True):
    #     """Add a file to the experiment
    #
    #     Add the file to the experiment only if the file object has found and imported relevant extensions .
    #     If the file is already present in experiment, then try to find and import new extensions.
    #
    #     Parameters
    #     ----------
    #     relativeFilePath : pathlib.Path or str
    #         Path with respect to the main experiment path
    #
    #     """
    #     # Perhaps move this conversion to relative file path to File
    #     relative_file_path = Path(file_path).absolute().relative_to(self.main_path)
    #
    #     # if there is no extension, add all files with the same name with all extensions
    #     # if there is an extension just add that file if the filename is the same
    #
    #     # Test whether file is already in experiment
    #     # for file in self.files:
    #     #     if file.relativeFilePath == relativeFilePath:
    #     #         file.findAndAddExtensions()
    #     #         break
    #     # else:
    #     #     new_file = File(relativeFilePath, self)
    #     #     if new_file.extensions:
    #     #         self.files.append(new_file)
    #
    #
    #     if not test_duplicates or (relative_file_path not in self.file_paths):
    #         self.files.append(File(relative_file_path, self))
    #     else:
    #         i = self.file_paths.find(file_path.absolute().relative_to(self.main_path))
    #         self.files[i].findAndAddExtensions()

    def get_fov_from_nd2(self, nd2_fullpath):
        images = ND2Reader(str(nd2_fullpath))
        y_positions = images._parser._raw_metadata.y_data   # nikon sample stage position
        x_positions = images._parser._raw_metadata.x_data   # nikon sample stage position

        # set the image data order in the nd2 file
        if 'c' in images.axes:
            images.iter_axes = 'tc'  # for alex measurements
        else:
            images.iter_axes = 't'

        n_illumination = len(images.metadata["channels"])
        n_frames = len(x_positions)
        position_tolerance = 10  # xy tol = tolerance in um
        first_frame_of_each_fov = [0]
        last_frame_of_each_fov = []
        for fri in range(n_frames - 1):
            if abs(x_positions[fri] - x_positions[fri + 1]) > position_tolerance or abs(y_positions[fri] - y_positions[fri + 1]) > position_tolerance:
                first_frame_of_each_fov.append(fri + 1)
                last_frame_of_each_fov.append(fri)
        last_frame_of_each_fov.append(n_frames-1)
        fov_info = {'number_of_fov': len(first_frame_of_each_fov),
                    'first_frame_of_each_fov': first_frame_of_each_fov,
                    'last_frame_of_each_fov': last_frame_of_each_fov}
        return fov_info


    def histogram(self, axis=None, bins=100, parameter='E', molecule_averaging=False,
                  fileSelection=False, moleculeSelection=False, makeFit=False, export=False, **kwargs):
        """FRET histogram of all molecules in the experiment or a specified selection

        Parameters
        ----------
        axis : matplotlib.axis
            Axis to use for histogram plotting
        bins : int
            Number of bins
        parameter : str
            Parameter to be used for histogram I or E
        molecule_averaging : bool
            If True an time average of the trace is used
        fileSelection : bool
            If True the histogram is made only using selected files.
        moleculeSelection : bool
            If True the histogram is made only using selected molecules.
        makeFit : bool
            If True perform Gaussian fitting.
        export : bool
            If True the graph is exported.
        **kwargs
            Arbitrary keyword arguments.

        """
        # files = [file for file in exp.files if file.isSelected]
        # files = self.files

        if (fileSelection & moleculeSelection):
            molecules = [molecule for file in self.selectedFiles for molecule in file.selectedMolecules]
        elif (fileSelection & (not moleculeSelection)):
            molecules = [molecule for file in self.selectedFiles for molecule in file.molecules]
        elif ((not fileSelection) & moleculeSelection):
            molecules = [molecule for file in self.files for molecule in file.selectedMolecules]
        else:
            molecules = [molecule for file in self.files for molecule in file.molecules]

        histogram(molecules, axis=axis, bins=bins, parameter=parameter, molecule_averaging=molecule_averaging,
                  makeFit=makeFit, collection_name=self, **kwargs)
        if export: plt.savefig(self.main_path.joinpath(f'{self.name}_{parameter}_histogram').with_suffix('.png'))

    def boxplot_number_of_molecules(self):
        """Boxplot of the number of molecules in each file"""
        fig, ax = plt.subplots(figsize=(8, 1.5))
        pointCount = [len(file.molecules) for file in self.files]
        plt.boxplot(pointCount, vert=False, labels=[''], widths=(0.8))
        plt.xlabel('Count')
        plt.title('Molecules per file')
        plt.tight_layout()

        fig.savefig(self.main_path.joinpath('number_of_molecules.pdf'), bbox_inches='tight')
        fig.savefig(self.main_path.joinpath('number_of_molecules.png'), bbox_inches='tight')

    def select(self):
        """Simple method to look through all molecules in the experiment

        Plots a molecule. If enter is pressed the next molecule is shown.

        """
        for molecule in self.molecules:
            molecule.plot()
            input("Press enter to continue")

    def print_files(self):
        for i, file in enumerate(self.files):
            print(f"{i:3d}.  {file.relativeFilePath}")

    def plot_trace(self, files=None, query={}):
        from trace_analysis.trace_plot import TraceAnalysisFrame
        import wx

        if files is None:
            files = self.files

        file_paths = [file.relativeFilePath.with_suffix('.nc') for file in files if '.nc' in file.extensions]

        with xr.open_mfdataset(file_paths, concat_dim='molecule', combine='nested') as ds:
            ds_sel = ds.query(query)  # HJ1_WT, HJ7_G116T
            app = wx.App(False)
            # app = wit.InspectableApp()
            frame = TraceAnalysisFrame(None, ds_sel, "Sample editor")
            # frame.molecules = exp.files[1].molecules
            # print('test')
            # import wx.lib.inspection
            # wx.lib.inspection.InspectionTool().Show()
            app.MainLoop()

    def export_number_of_molecules_per_file(self):
        df = pd.DataFrame(columns=['Number of molecules'])
        for i, file in enumerate(self.files):
            n = str(file.relativeFilePath)
            try:
                nms = file.number_of_molecules
            except FileNotFoundError:
                nms = -1
            df.loc[n] = nms
        df.to_excel(self.main_path.joinpath('number_of_molecules'))
