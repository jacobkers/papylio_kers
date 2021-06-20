import numpy as np #scientific computing with Python
import matplotlib.pyplot as plt #Provides a MATLAB-like plotting framework
import pandas as pd
import xarray as xr
from pathlib import Path
from trace_analysis.analysis.autoThreshold import stepfinder
# from trace_analysis.plugin_manager import PluginManager
# from trace_analysis.plugin_manager import PluginMetaClass
from trace_analysis.plugin_manager import plugins
from trace_analysis.trace_extraction import make_gaussian
import copy

class Molecules:
    def load(filepath):
        return Molecules().load(filepath)

    def __init__(self, name=''):
        # for value in molecules:
        #     if not isinstance(value, Molecule):
        #         raise TypeError('MoleculeList can only contain Molecule objects')
        # self.molecules = molecules
        self.name = name

        self._dataset = xr.Dataset(
            coords=
            {
                'molecule':     ('molecule', pd.MultiIndex.from_tuples([],names=['molecule_in_file','file'])),
                'frame':        ('frame', np.array([], dtype=int)),
                'channel':      ('channel', np.array([], dtype=int))
            }
        )


    @property
    def dataset(self):
        return self._dataset

    @dataset.setter
    def dataset(self, dataset):
        if not self._dataset:
            self._dataset = dataset
            if 'selected' not in self._dataset.keys():
                self.dataset['selected'] = xr.DataArray(False, coords=[dataset.molecule])
        else:
            self._dataset = dataset

        # self.empty_parameter_index = pd.MultiIndex.from_arrays([['is_selected'],[]], names=['Parameter','Channel'])
        # self.empty_traces_index = pd.MultiIndex.from_arrays([[],[]], names=['Frame','Channel'])
        # self._traces = None
        # self._parameters = None
        # self.channels = [0,1]

    # @property
    # def traces(self):
    #     return self._traces
    #
    # @traces.setter
    # def traces(self, traces):
    #     if self._traces is not None and (len(traces.columns) != len(self.molecule_index)):
    #         raise ValueError('Number of molecules does not match current number of molecules, to proceed reset the object first ')
    #     self._traces = traces
    #     self.molecule_index = traces.columns
    #
    # @property
    # def parameters(self):
    #     return self._parameters
    #
    # @parameters.setter
    # def parameters(self, parameters):
    #     if self._parameters is not None and (len(parameters.columns) != len(self.molecule_index)):
    #         raise ValueError('Number of molecules does not match current number of molecules, to proceed reset the object first ')
    #     self._parameters = parameters
    #     self.molecule_index = parameters.columns
    #
    # @property
    # def molecule_index(self):
    #     if self._traces is not None:
    #         return self.traces.columns
    #
    # @molecule_index.setter
    # def molecule_index(self, molecule_index):
    #     if self._traces is None:
    #         self._traces = pd.DataFrame(index=self.empty_traces_index, columns=molecule_index)
    #     else:
    #         self._traces.columns = molecule_index
    #     if self._parameters is None:
    #         self._parameters = pd.DataFrame(index=self.empty_parameter_index, columns=molecule_index)
    #     else:
    #         self._parameters.columns = self.molecule_index
    #
    # def reset(self):
    #     self._traces = None
    #     self._parameters = None

    def __len__(self):
        return len(self.dataset.molecule)

    def __getitem__(self, item):
        new = Molecules()
        new.dataset = self.dataset.isel(molecule=item)
        return new

    def __add__(self, other):
        return xr.concat([self.dataset, other.data_set], dim='molecule')

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)

    def __getattr__(self, name):
        try:
            return self.dataset[name]
        except KeyError:
            super().__getattribute__(name)

    # def __setattr__(self, key, value):
    #     if key in self.__dict__.keys():
    #         super().__setattr__(key, value)
    #     else:
    #         self._parameters.loc[key,:] = value

    # def add_parameters(self, added_parameters):
    #     parameters = self.parameters.append(added_parameters)
    #     self.parameters = parameters[~parameters.index.duplicated(keep='last')]
    #
    # def add_parameter_from_list(self, name, parameter_list, channel=''):
    #     self.parameters.loc[(name, channel), :] = parameter_list

    def import_file(self, filepath):
        filepath = Path(filepath)
        if filepath.suffix == '.traces':
            self.import_traces_file(filepath)
        elif filepath.suffix == '.pks':
            self.import_traces_file(filepath)
        else:
            raise FileNotFoundError(filepath)

    def export_file(self, filepath):
        filepath = Path(filepath)
        if filepath.suffix == '.traces':
            self.export_traces_file(filepath)
        elif filepath.suffix == '.pks':
            self.export_traces_file(filepath)
        else:
            raise FileNotFoundError(filepath)

    def import_pks_file(self, pks_filepath):
        peaks = import_pks_file(pks_filepath)
        peaks = split_dimension(peaks, 'peak', ('molecule', 'channel'), (-1, 2))
        peaks = split_dimension(peaks, 'molecule', ('molecule_in_file', 'file'), (-1, 1), to='multiindex')
        self.dataset = self.dataset.merge(peaks.to_dataset('parameter'))

    def export_pks_file(self, pks_filepath):
        peaks = self.dataset[['x','y','background']].reset_index('molecule')\
            .stack(peaks=('molecule', 'channel')).to_array(dim='parameter').T
        export_pks_file(peaks, pks_filepath)

    def import_traces_file(self, traces_filepath):
        traces = import_traces_file(traces_filepath)
        traces = split_dimension(traces, 'trace', ('molecule', 'channel'), (-1, 2))
        traces = split_dimension(traces, 'molecule', ('molecule_in_file', 'file'), (-1, 1), to='multiindex')
        dataset = xr.Dataset({'traces': traces})
        self.dataset = self.dataset.merge(dataset)

    def export_traces_file(self, traces_filepath):
        traces = self.traces.reset_index('molecule').stack(trace=('molecule', 'channel')).T
        export_traces_file(traces, traces_filepath)

    def save(self, filepath):
        filepath = Path(filepath)
        self.dataset.to_netcdf(filepath.with_suffix('.nc'))

    def load(self, filepath):
        filepath = Path(filepath)
        self.dataset = xr.open_dataset(filepath.with_suffix('.nc'))


def import_pks_file(pks_filepath):
    pks_filepath = Path(pks_filepath)
    data = np.atleast_2d(np.genfromtxt(pks_filepath)[:,1:])
    if data.shape[1] == 2:
        data = np.hstack([data, np.zeros((len(data),1))])
    return xr.DataArray(data, dims=("peak",'parameter'),
                        coords={'peak': range(len(data)), 'parameter': ['x', 'y', 'background']})


    #
    # number_of_molecules = len(data) // number_of_channels
    # index = pd.MultiIndex.from_product([[pks_filepath.with_suffix('').name],
    #                                     np.arange(number_of_molecules),
    #                                     np.arange(number_of_channels)],
    #                                     names=['File','Molecule', 'Channel'])
    # columns = pd.Index(['x', 'y', 'background'], name='Parameter')
    # return pd.DataFrame(data=data, index=index, columns=columns).T.stack('Channel')

# def export_pks_file(coordinates_and_background, pks_filepath):
#     coordinates_and_background = coordinates_and_background.unstack('Channel').T
#     pks_filepath = Path(pks_filepath)
#     with pks_filepath.open('w') as pks_file:
#         for i, (background, x, y) in enumerate(coordinates_and_background.values):
#             # outfile.write(' {0:4.0f} {1:4.4f} {2:4.4f} {3:4.4f} {4:4.4f} \n'.format(i, coordinate[0], coordinate[1], 0, 0, width4=4, width6=6))
#             # pks_file.write('{0:4.0f} {1:4.4f} {2:4.4f} \n'.format(i + 1, coordinate[0], coordinate[1]))
#             pks_file.write('{0:4.0f} {1:4.4f} {2:4.4f} {3:4.4f}\n'.format(i + 1, x, y, background))

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
    #                         order='F')  # 3d array of traces
    traces = np.reshape(raw_data, (number_of_frames, number_of_traces)).T  # 2d array of traces
    traces = xr.DataArray(traces, dims=("trace", "frame"), coords=(range(number_of_traces), range(number_of_frames)))
    return traces


    # column_index = pd.MultiIndex.from_product([[traces_filepath.with_suffix('').name],
    #                                     np.arange(number_of_molecules),
    #                                     np.arange(number_of_channels)],
    #                                    names=['File', 'Molecule', 'Channel'])
    # index = pd.Index(data=np.arange(number_of_frames), name='Frame')
    # return pd.DataFrame(traces.T, index=index, columns=column_index).stack('Channel')

def split_dimension(data_array, old_dim, new_dims, new_dims_shape, to='dimensions'):
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

    new_index = pd.MultiIndex.from_product((range(dim_len) for dim_len in new_dims_shape), names=new_dims)
    data_array = data_array.assign_coords(**{old_dim: new_index})

    if to == 'dimensions':
        return data_array.unstack(old_dim).transpose(*all_dims)
    elif to == 'multiindex':
        return data_array
    else:
        raise ValueError

# test.reset_index('trace').set_index(trace=['file','molecule','channel'])
# test3.reset_index('trace').reset_coords('file', drop=True).assign_coords({'file': ('trace',[1]*278)}).set_index(trace=['file','molecule','channel'])
# number_of_channels = 2
# number_of_molecules = a.shape[a.dims=='trace']
# traces = traces.unstack('Channel')
def export_traces_file(traces, traces_filepath):
    traces_filepath = Path(traces_filepath)
    with traces_filepath.open('w') as traces_file:
        # Number of frames
        np.array([len(traces.frame)], dtype=np.int32).tofile(traces_file)
        # Number of traces
        np.array([len(traces.trace)], dtype=np.int16).tofile(traces_file)
        traces.values.T.astype(np.int16).tofile(traces_file)



@plugins
class Molecule:
    pass
#     slots = ('file', 'index', '_coordinates', 'intensity', '_background', 'is_selected', 'steps', 'kon_boolean')
#
#     def __init__(self, file):
#         self.file = file
#         self.index = None
#         self._coordinates = None
#         self.intensity = None
#         self._background = None
#
#         self.is_selected = False
#
#         self.steps = None  #Defined in other classes as: pd.DataFrame(columns=['frame', 'trace', 'state', 'method','thres'])
#         self.kon_boolean = None  # 3x3 matrix that is indicates whether the kon will be calculated from the beginning, in-between molecules or for the end only
#         #self.bg_scale=np.sum(make_gaussian(self.file.experiment.configuration['find_coordinates']['coordinate_optimization']['coordinates_after_gaussian_fit']['gaussian_width']))
#
#     @property
#     def coordinates(self):
#         return self._coordinates
#
#     @property
#     def background(self):
#         return self._background
#
#     @coordinates.setter
#     def coordinates(self, coordinates):
#         self._coordinates = np.atleast_2d(coordinates)
#
#     def background(self, background):
#         self.background=background # should be dependent on emission channel as well
#
#     @property  # this is just for the stepfinder to be called through Molecule. Maybe not needed
#     def find_steps(self):
#         return stepfinder
#
#     def I(self, emission, Ioff=0):
#         return self.intensity[emission, :] - Ioff # - self.background[emission] * self.bg_scale #this number comes from sum(make_gaussian) in trace_extraction
#
#     def E(self, Imin=0, Iroff=0, Igoff=0, alpha=0):
#         red = np.copy(self.I(1, Ioff=Iroff))
#         green = self.I(0, Ioff=Igoff)
#         np.putmask(green, green < 0, 0) # green < 0 is taken as 0
#         np.putmask(red, red < Imin, 0)  # the mask makes all elements of acceptor that are below the Imin zero, for E caclulation
#         E =  (red - alpha*green) / (green + red - alpha*green)
#         E = np.nan_to_num(E)  # correct for divide with zero = None values
#         return E
#
#     def plot(self, ylim=(0, 500), xlim=(), Ioff=[],  save=False, **fretkwargs):
#         plt.style.use('seaborn-dark')
#         plt.style.use('seaborn-colorblind')
#         figure = plt.figure(f'{self.file.name}_mol_{self.index}', figsize=(7,4))
#         if len(self.file.experiment.pairs) > 0:
#             axis_I = figure.add_subplot(211)
#         else:
#             axis_I = figure.gca()
#
#         axis_I.set_ylabel('Intensity (a.u.)')
#         axis_I.set_ylim(ylim[0], ylim[1])
#         if xlim == ():
#             axis_I.set_xlim(0, self.file.time.max()+1)
#         else:
#             axis_I.set_xlim(xlim[0], xlim[1])
#
#         axis_I.set_title(f'Molecule {self.index} /{len(self.file.molecules)}')
#         if Ioff == []:
#             Ioff = [0]*self.file.number_of_channels
#         for i, channel in enumerate(self.file.experiment.channels):
#             axis_I.plot(self.file.time, self.I(i, Ioff=Ioff[i]), channel)
#
#         if len(self.file.experiment.pairs) > 0:
#             axis_E = figure.add_subplot(212, sharex=axis_I)
#             axis_E.set_xlabel('Time (s)')
#             axis_E.set_ylabel('FRET')
#             axis_E.set_ylim(0,1.1)
#             for i, pair in enumerate(self.file.experiment.pairs):
#                 axis_E.plot(self.file.time, self.E(**fretkwargs), 'b')
#
#         plt.tight_layout()
#         if save:
#             plt.savefig(f'{self.file.relativeFilePath}_mol_{self.index}.eps', transparent=True)
#             plt.savefig(f'{self.file.relativeFilePath}_mol_{self.index}.png', facecolor='white', dpi=300, transparent=True)
#





if __name__ == '__main__':
    filepath = r'D:\SURFdrive\Promotie\Code\Python\traceAnalysis\twoColourExampleData\20141017 - Holliday junction - Copy\test'
    # traces_filepath = r'P:\SURFdrive\Promotie\Code\Python\traceAnalysis\twoColourExampleData\20141017 - Holliday junction - Copy\test'
    filepath = Path(filepath)

    test = Molecules()
    test.import_traces_file(filepath.with_suffix('.traces'))
    test.import_pks_file(filepath.with_suffix('.pks'))

    # pks_filepath = Path(traces_filepath).with_suffix('.pks')
    # test.import_pks_file(pks_filepath, 2)
    # test.export_pks_file(pks_filepath.with_name('test2.pks'))
    #
    #
    # test[5:10]
    #
    # test[5]+test[10]


    test.export_traces_file(filepath.with_name('test2.traces'))
    test.export_pks_file(filepath.with_name('test2.pks'))