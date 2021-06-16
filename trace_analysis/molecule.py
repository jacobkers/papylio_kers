import numpy as np #scientific computing with Python
import matplotlib.pyplot as plt #Provides a MATLAB-like plotting framework
import pandas as pd
from pathlib import Path
from trace_analysis.analysis.autoThreshold import stepfinder
# from trace_analysis.plugin_manager import PluginManager
# from trace_analysis.plugin_manager import PluginMetaClass
from trace_analysis.plugin_manager import plugins
from trace_analysis.trace_extraction import make_gaussian
import copy

class Molecules:
    def __init__(self, molecules=[], name=''):
        # for value in molecules:
        #     if not isinstance(value, Molecule):
        #         raise TypeError('MoleculeList can only contain Molecule objects')
        # self.molecules = molecules
        self.name = name
    #
    # def __len__(self):
    #     return self.sequence.shape[0]

    def __getitem__(self, item):
        # new = copy.copy(self)
        new = Molecules(name=self.name)
        for attribute, value in vars(self).items():
            if isinstance(value, pd.DataFrame):
                #new_value = value.xs(0, axis=1, level='Molecule', drop_level=False)
                if isinstance(item, int):
                    item = [item]
                new.__setattr__(attribute, value.iloc[:, item])
            elif isinstance(value, list):
                new.__setattr__(attribute, value[item])
            else:
                new.__setattr__(attribute, value)
        return new

    def __add__(self, other):
        new = Molecules(name=self.name)
        for attribute, value in vars(self).items():
            if isinstance(value, pd.DataFrame):
                #new_value = value.xs(0, axis=1, level='Molecule', drop_level=False)
                new_value = pd.concat([value, other.__getattribute__(attribute)], axis=1)
                new.__setattr__(attribute, new_value)
            elif isinstance(value, list):
                new.__setattr__(attribute, value+other.__getattribute__(attribute))
            else:
                new.__setattr__(attribute, value)
        return new

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)
    # def __getattribute__(self, item):
    #     try:
    #         datatype = type(self.molecules[0].item)
    #         if isinstance(datatype, pd.DataFrame):
    #             [getattr(molecule, item) for molecule in self.molecules]
    #     except AttributeError
    #         super().__getattribute__(item)
    #
    #
    # def __setattr__(self, key, value):
    #     if isinstance(value, pd.DataFrame):
    #         number_of_molecules = value.index.levshape[[name == 'Molecule' for name in value.index.names].index(True)]
    #         if number_of_molecules != self.__len__():
    #             raise ValueError('Wrong number of molecules')
    #         for molecule in self.molecules:
    #             molecule.key = value.loc[]

    def import_traces_file(self, traces_filepath, number_of_channels):
        self.traces = import_traces_file(traces_filepath, number_of_channels)

    def export_traces_file(self, traces_filepath):
        export_traces_file(self.traces, traces_filepath)



def import_traces_file(traces_filepath, number_of_channels):
    traces_filepath = Path(traces_filepath)
    with traces_filepath.open('r') as traces_file:
        number_of_frames = np.fromfile(traces_file, dtype=np.int32, count=1).item()
        number_of_traces = np.fromfile(traces_file, dtype=np.int16, count=1).item()
        number_of_molecules = number_of_traces // number_of_channels
        rawData = np.fromfile(traces_file, dtype=np.int16, count=number_of_frames * number_of_traces)
    # traces = np.reshape(rawData.ravel(),
    #                         (number_of_channels, number_of_molecules, number_of_frames),
    #                         order='F')  # 3d array of traces
    traces = np.reshape(rawData.ravel(), (number_of_channels * number_of_molecules, number_of_frames),
                        order='F')  # 2d array of traces

    column_index = pd.MultiIndex.from_product([[traces_filepath.with_suffix('')],
                                        np.arange(number_of_molecules),
                                        np.arange(number_of_channels)],
                                       names=['File', 'Molecule', 'Channel'])
    index = pd.Index(data=np.arange(number_of_frames), name='Frame')
    return pd.DataFrame(traces.T, index=index, columns=column_index).stack('Channel')

def export_traces_file(traces, traces_filepath):
    traces_filepath = Path(traces_filepath)
    traces = traces.unstack('Channel')
    with traces_filepath.open('w') as traces_file:
        # Number of frames
        np.array([len(traces.index)], dtype=np.int32).tofile(traces_file)
        # Number of traces
        np.array([len(traces.columns)], dtype=np.int16).tofile(traces_file)
        traces.values.astype(np.int16).tofile(traces_file)

@plugins
class Molecule:
    slots = ('file', 'index', '_coordinates', 'intensity', '_background', 'is_selected', 'steps', 'kon_boolean')

    def __init__(self, file):
        self.file = file
        self.index = None
        self._coordinates = None
        self.intensity = None
        self._background = None

        self.is_selected = False

        self.steps = None  #Defined in other classes as: pd.DataFrame(columns=['frame', 'trace', 'state', 'method','thres'])
        self.kon_boolean = None  # 3x3 matrix that is indicates whether the kon will be calculated from the beginning, in-between molecules or for the end only
        #self.bg_scale=np.sum(make_gaussian(self.file.experiment.configuration['find_coordinates']['coordinate_optimization']['coordinates_after_gaussian_fit']['gaussian_width']))

    @property
    def coordinates(self):
        return self._coordinates
    
    @property
    def background(self):
        return self._background
    
    @coordinates.setter
    def coordinates(self, coordinates):
        self._coordinates = np.atleast_2d(coordinates)

    def background(self, background):
        self.background=background # should be dependent on emission channel as well
        
    @property  # this is just for the stepfinder to be called through Molecule. Maybe not needed
    def find_steps(self):
        return stepfinder

    def I(self, emission, Ioff=0):
        return self.intensity[emission, :] - Ioff # - self.background[emission] * self.bg_scale #this number comes from sum(make_gaussian) in trace_extraction

    def E(self, Imin=0, Iroff=0, Igoff=0, alpha=0):
        red = np.copy(self.I(1, Ioff=Iroff))
        green = self.I(0, Ioff=Igoff)
        np.putmask(green, green < 0, 0) # green < 0 is taken as 0
        np.putmask(red, red < Imin, 0)  # the mask makes all elements of acceptor that are below the Imin zero, for E caclulation
        E =  (red - alpha*green) / (green + red - alpha*green)
        E = np.nan_to_num(E)  # correct for divide with zero = None values
        return E

    def plot(self, ylim=(0, 500), xlim=(), Ioff=[],  save=False, **fretkwargs):
        plt.style.use('seaborn-dark')
        plt.style.use('seaborn-colorblind')
        figure = plt.figure(f'{self.file.name}_mol_{self.index}', figsize=(7,4))
        if len(self.file.experiment.pairs) > 0:
            axis_I = figure.add_subplot(211)
        else:
            axis_I = figure.gca()

        axis_I.set_ylabel('Intensity (a.u.)')
        axis_I.set_ylim(ylim[0], ylim[1])
        if xlim == ():
            axis_I.set_xlim(0, self.file.time.max()+1)
        else:
            axis_I.set_xlim(xlim[0], xlim[1])

        axis_I.set_title(f'Molecule {self.index} /{len(self.file.molecules)}')
        if Ioff == []:
            Ioff = [0]*self.file.number_of_channels
        for i, channel in enumerate(self.file.experiment.channels):
            axis_I.plot(self.file.time, self.I(i, Ioff=Ioff[i]), channel)

        if len(self.file.experiment.pairs) > 0:
            axis_E = figure.add_subplot(212, sharex=axis_I)
            axis_E.set_xlabel('Time (s)')
            axis_E.set_ylabel('FRET')
            axis_E.set_ylim(0,1.1)
            for i, pair in enumerate(self.file.experiment.pairs):
                axis_E.plot(self.file.time, self.E(**fretkwargs), 'b')

        plt.tight_layout()
        if save:
            plt.savefig(f'{self.file.relativeFilePath}_mol_{self.index}.eps', transparent=True)
            plt.savefig(f'{self.file.relativeFilePath}_mol_{self.index}.png', facecolor='white', dpi=300, transparent=True)






if __name__ == '__main__':
    traces_filepath = r'D:\SURFdrive\Promotie\Code\Python\traceAnalysis\twoColourExampleData\20141017 - Holliday junction - Copy\test.traces'
    test = Molecules()
    test.import_traces_file(traces_filepath, 2)
    test[5:10]

    test[5]+test[10]