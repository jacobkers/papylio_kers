from pathlib import Path # For efficient path manipulation
import numpy as np #scientific computing with Python
import matplotlib.pyplot as plt #Provides a MATLAB-like plotting framework
from trace_analysis.molecule import Molecule
from trace_analysis.image_adapt.sifx_file import SifxFile
from trace_analysis.image_adapt.pma_file import PmaFile
from trace_analysis.plotting import histogram
from trace_analysis.mapping.mapping import Mapping2

class File:
    def __init__(self, relativeFilePath, experiment):
        relativeFilePath = Path(relativeFilePath)
        self.experiment = experiment

        self.relativePath = relativeFilePath.parent
        self.name = relativeFilePath.name
        self.extensions = list()

        self.molecules = list()
        self.exposure_time = None  # Here the exposure time is given but it should be found from the log file if possible
        self.number_of_frames = None

        self.background = np.array([0, 0])

        self.isSelected = False
        self.is_mapping_file = False

        self.movie = None
        self.mapping = None

        if self.experiment.import_all is True:
            self.findAndAddExtensions()

    def __repr__(self):
        return (f'{self.__class__.__name__}({self.relativePath.joinpath(self.name)})')

    @property
    def relativeFilePath(self):
        return self.relativePath.joinpath(self.name)

    @property
    def absoluteFilePath(self):
        return self.experiment.mainPath.joinpath(self.relativeFilePath)

    @property
    def number_of_molecules(self):
        return len(self.molecules)

    @number_of_molecules.setter
    def number_of_molecules(self, number_of_molecules):
        if not self.molecules:
            for molecule in range(0, number_of_molecules):
                self.addMolecule()
        elif number_of_molecules != self.number_of_molecules:
            raise ValueError('Requested number of molecules differs from existing number of molecules')

    @property
    def number_of_colours(self):
        return self.experiment.Ncolours

    @property
    def selectedMolecules(self):
        return [molecule for molecule in self.molecules if molecule.isSelected]

    @property
    def coordinates(self):
        # if not self._pks_file:
        #     _pks_file = PksFile(self.absoluteFilePath.with_suffix('.pks'))

        return np.concatenate([[molecule.coordinates[0, :] for molecule in self.molecules]])

    @coordinates.setter
    def coordinates(self, coordinates, number_of_colours = None):
        if number_of_colours is None: number_of_colours = self.number_of_colours
        self.number_of_molecules = np.shape(coordinates)[0]//number_of_colours

        for i, molecule in enumerate(self.molecules):
            molecule.coordinates = coordinates[(i * number_of_colours):((i + 1) * number_of_colours), :]

    @property
    def traces(self):
        np.dstack([molecule.intensity for molecule in self.molecules]).swapaxes(1, 2) # 3d array of traces
        # np.concatenate([molecule.intensity for molecule in self.molecules]) # 2d array of traces

    @traces.setter
    def traces(self, traces):
        for i, molecule in enumerate(self.molecules):
            molecule.intensity = traces[:, i, :] # 3d array of traces
            # molecule.intensity = traces[(i * self.number_of_colours):((i + 1) * self.number_of_colours), :] # 2d array of traces

    def findAndAddExtensions(self):
        foundFiles = [file.name for file in self.experiment.mainPath.joinpath(self.relativePath).glob(self.name + '*')]
        foundExtensions = [file[len(self.name):] for file in foundFiles]

        # For the special case of a sifx file, which is located inside a folder
        if '' in foundExtensions: foundExtensions[foundExtensions.index('')] = '.sifx'

        newExtensions = [extension for extension in foundExtensions if extension not in self.extensions]
        self.extensions = self.extensions + newExtensions
        for extension in newExtensions: self.importExtension(extension)

    def importExtension(self, extension):

        # print(f.relative_to(self.experiment.mainPath))

        if extension not in self.extensions:
            self.extensions.append(extension)

        # print(extension)
        importFunctions = { '.sifx': self.import_sifx_file,
                            '.pma': self.import_pma_file,
                            '.coeff': self.import_coeff_file,
                            '.map': self.import_map_file,
                            '.pks': self.import_pks_file,
                            '.traces': self.import_traces_file
                            #                           '.sim'        : self.importSimFile
                            }

        importFunctions.get(extension, self.noneFunction)()

    #        if extension == '.pks':
    # self.importPksFile()

    def noneFunction(self):
        return

    def import_sifx_file(self):
        imageFilePath = self.absoluteFilePath.joinpath('Spooled files.sifx')
        self.movie = SifxFile(imageFilePath)

    def import_pma_file(self):
        imageFilePath = self.absoluteFilePath.with_suffix('.pma')
        self.movie = PmaFile(imageFilePath)

    def import_coeff_file(self):
        if self.mapping is None:
            self.mapping = Mapping2(transformation_type='linear')
            self.mapping.transformation = np.zeros((3,3))
            self.mapping.transformation[2,2] = 1
            self.mapping.transformation[[0,0,0,1,1,1],[2,0,1,2,0,1]] = \
                np.genfromtxt(str(self.relativeFilePath) + '.coeff')

    def import_map_file(self):
        coefficients = np.genfromtxt(p.joinpath('rough').with_suffix('.map'))
        degree = int(np.sqrt(len(coefficients) // 2) - 1)
        P = coefficients[:len(coefficients) // 2].reshape((degree + 1, degree + 1))
        Q = coefficients[len(coefficients) // 2 : len(coefficients)].reshape((degree + 1, degree + 1))

        self.mapping = Mapping2(transformation_type='polynomial')
        self.mapping.transformation = {'P': P, 'Q': Q}

    def import_pks_file(self):
        # Background value stored in pks file is not imported yet
        coordinates = np.genfromtxt(str(self.relativeFilePath) + '.pks')
        coordinates = np.atleast_2d(coordinates)[:,1:3]

        self.coordinates = coordinates

    def export_pks_file(self):
        pks_filepath = self.absoluteFilePath.with_suffix('.pks')
        with pks_filepath.open('w') as pks_file:
            for i, coordinate in enumerate(self.coordinates):
                # outfile.write(' {0:4.0f} {1:4.4f} {2:4.4f} {3:4.4f} {4:4.4f} \n'.format(i, coordinate[0], coordinate[1], 0, 0, width4=4, width6=6))
                pks_file.write('{0:4.0f} {1:4.4f} {2:4.4f} \n'.format(i + 1, coordinate[0], coordinate[1]))

    def import_traces_file(self):
        traces_filepath = self.absoluteFilePath.with_suffix('.traces')
        with traces_filepath.open('r') as traces_file:
            self.number_of_frames = np.fromfile(traces_file, dtype=np.int32, count=1).item()
            number_of_traces = np.fromfile(traces_file, dtype=np.int16, count=1).item()
            self.number_of_molecules = number_of_traces // self.number_of_colours
            rawData = np.fromfile(traces_file, dtype=np.int16, count=self.number_of_frames * number_of_traces)
        self.traces = np.reshape(rawData.ravel(), (self.number_of_colours, self.number_of_molecules, self.number_of_frames), order='F')  # 3d array of traces
        #self.traces = np.reshape(rawData.ravel(), (self.number_of_colours * self.number_of_molecules, self.number_of_frames), order='F') # 2d array of traces

    def export_traces_file(self):
        traces_filepath = self.writepath.joinpath(self.name + '.traces')
        with traces_filepath.open('w') as traces_file:
            np.array([traces.shape[2]], dtype=np.int32).tofile(traces_file)
            np.array([traces.shape[0]*traces.shape[1]], dtype=np.int16).tofile(traces_file)
            # time_tr = np.zeros((self.number_of_frames, 2 * self.pts_number))
            # Ncolours=2
            # for jj in range(2*self.pts_number//Ncolours):
            #     time_tr[:,jj*2] = donor[:,jj]
            #     time_tr[:,jj*2+1]=  acceptor[:,jj]
            np.array(traces.T, dtype=np.int16).tofile(traces_file)

    #    def importSimFile(self):
    #        file = open(str(self.relativeFilePath) + '.sim', 'rb')
    #        self.data = pickle.load(file)
    #        red, green  = self.data['red'], self.data['green']
    #        Ntraces = red.shape[0]
    #        self.Nframes = red.shape[1]
    #
    #        if not self.molecules:
    #            for molecule in range(0, Ntraces):
    #                self.addMolecule()
    #
    #        for i, molecule in enumerate(self.molecules):
    #            molecule.intensity = np.vstack((green[i], red[i]))
    #        file.close()

    def addMolecule(self):
        self.molecules.append(Molecule(self))
        self.molecules[-1].index = len(self.molecules)  # this is the molecule number

    def histogram(self):
        histogram(self.molecules)

    #    def importExcel(self, filename=None):
    #        if filename is None:
    #            filename = self.name+'_steps_data.xlsx'
    #        try:
    #            steps_data = pd.read_excel(filename, index_col=[0,1],
    #                                            dtype={'kon':np.str})       # reads from the 1st excel sheet of the file
    #        except FileNotFoundError:
    #            return
    #        molecules = steps_data.index.unique(0)
    #        indices = [int(m.split()[-1]) for m in molecules]
    #        for mol in self.molecules:
    #            if mol.index not in indices:
    #                continue
    #            mol.steps = steps_data.loc[f'mol {mol.index}']
    #            k = [int(i) for i in mol.steps.kon[0]]
    #            mol.kon_boolean = np.array(k).astype(bool).reshape((3,3))

    def select(self, figure=None):
        plt.ion()
        for molecule in self.molecules:
            molecule.plot(figure=figure)
            plt.show()
            plt.pause(0.001)
            input("Press enter to continue")

    def use_for_mapping(self):
        self.is_mapping_file = True
        mapping = self.movie.use_for_mapping()
        for file in self.experiment.files:
            file.movie.mapping = mapping