# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 15:24:46 2018

@author: ivoseverins
"""

import os
import re # Regular expressions
import warnings
import numpy as np

# Use the following instead of: import matplotlib as mpl
#from matplotlib import use
#use('WXAgg')
from matplotlib import pyplot as plt

import itertools
from .threshold_analysis import stepfinder, plot_steps
import multiprocessing as mp
from functools import wraps, partial

class Experiment:
    def __init__(self, mainPath, exposure_time=0.1):
        self.name = os.path.basename(mainPath)
        self.mainPath = mainPath
        self.files = list()
        self.Ncolours = 2
        self.exposure_time = exposure_time

        os.chdir(mainPath)

        self.addAllFilesInMainPath()

    @property
    def molecules(self):
        return [molecule for file in self.files for molecule in file.molecules]
    @property
    def selectedFiles(self):
        return [file for file in self.files if file.isSelected]
    @property
    def selectedMoleculesInSelectedFiles(self):
        return [molecule for file in self.selectedFiles for molecule in file.selectedMolecules]
    @property
    def selectedMoleculesInAllFiles(self):
        return [molecule for file in self.files for molecule in file.selectedMolecules]

    def addAllFilesInMainPath(self):
        # Only works/tested for files in the main folder for now

        for root, dirs, fileNames in os.walk('.', topdown=False):

            for fileName in fileNames:
                    print(root)
                    print(fileName)

                    self.addFile(root, fileName)

            for name in dirs:
                warnings.warn('Import of files in subfolders does not work properly yet')
                print(os.path.join(root, name))

    def addFile(self, relativePath, fileName):
        fileNameSearch = re.search('(hel[0-9]*)(.*\..*)', fileName)

        if fileNameSearch:
            name, extension = fileNameSearch.groups()

            fullPath = os.path.join(relativePath, name)

            # Test whether file is already in experiment
            for file in self.files:
                if file.fullPath == fullPath:
                    foundFile = file
                    break
            else:
                foundFile = None

            # If not found: add file and extension, or if the file is already there then add the extention to it.
            # If the file and extension are already imported, display a warning message.
            if not foundFile:
                self.files.append(File(relativePath, name, self, self.exposure_time))
                self.files[-1].addExtension(extension)
            else:
                if extension not in foundFile.extensions:
                    foundFile.addExtension(extension)
                else:
                    warnings.warn(fullPath + extension + ' already imported')
        else:
            warnings.warn('FileName ' + fileName + ' not added, filename should contain hel...')


    def histogram(self, axis = None, fileSelection = False, moleculeSelection = False, makeFit = False):
        #files = [file for file in exp.files if file.isSelected]
        #files = self.files

        if (fileSelection & moleculeSelection):
            histogram([molecule for file in self.selectedFiles for molecule in file.selectedMolecules], axis, makeFit)
        elif (fileSelection & (not moleculeSelection)):
            histogram([molecule for file in self.selectedFiles for molecule in file.molecules], axis, makeFit)
        elif ((not fileSelection) & moleculeSelection):
            histogram([molecule for file in self.files for molecule in file.selectedMolecules], axis, makeFit)
        else:
            histogram([molecule for file in self.files for molecule in file.molecules], axis, makeFit)


    def select(self):
        for molecule in self.molecules:
            molecule.plot()
            input("Press enter to continue")

class File:
    def __init__(self, relativePath, name, experiment, exposure_time):
        self.relativePath = relativePath
        self.name = name
        self.extensions = list()
        self.experiment = experiment
        self.molecules = list()
        self.exposure_time = exposure_time

        self.isSelected = False

    @property
    def fullPath(self):
        return os.path.join(self.relativePath, self.name)

    @property
    def coordinates(self):
        return np.concatenate([[molecule.coordinates[0, :]
                                for molecule in self.molecules]])

    @property
    def selectedMolecules(self):
        return [molecule for molecule in self.molecules if molecule.isSelected]

    def addExtension(self, extension):
        self.extensions.append(extension)
        print(extension)
        importFunctions = {'.pks'        : self.importPksFile,
                           '.traces'     : self.importTracesFile}

        importFunctions.get(extension, print)()
#        if extension == '.pks':
            #self.importPksFile()

    def importPksFile(self):
        # Background value stored in pks file is not imported yet
        Ncolours = self.experiment.Ncolours
        pks = np.genfromtxt(self.name + '.pks', delimiter='      ')
        Ntraces = np.shape(pks)[0]

        if not self.molecules:
            for molecule in range(0, Ntraces, Ncolours):
                self.addMolecule()

#        for trace in range(0, Ntraces, Ncolours):
#            coordinates = pks[trace:(trace+Ncolours),1:3]
#            self.addMolecule(coordinates)

        for i, molecule in enumerate(self.molecules):
            molecule.coordinates = pks[(i*Ncolours):((i+1)*Ncolours), 1:3]

    def importTracesFile(self):
        Ncolours = self.experiment.Ncolours

        file = open(self.name + '.traces', 'r')

        self.Nframes = np.fromfile(file, dtype=np.int32, count=1).item()
        Ntraces = np.fromfile(file, dtype=np.int16, count=1).item()

        if not self.molecules:
            for molecule in range(0, Ntraces, Ncolours):
                self.addMolecule()

        rawData = np.fromfile(file, dtype=np.int16, count=self.Nframes * Ntraces)
        orderedData = np.reshape(rawData.ravel(), (Ncolours, Ntraces//Ncolours, self.Nframes), order = 'F')

        for i, molecule in enumerate(self.molecules):
            molecule.intensity = orderedData[:,i,:]

    def addMolecule(self):
        self.molecules.append(Molecule(self, self.exposure_time))

    def histogram(self):
        histogram(self.molecules)

    def find_dwell_times(self, trace, **params):
        dwell_times = [mol.find_steps(trace, **params)["dwell_times"]
                       for mol in self.molecules]
        return dwell_times

#    def find_molecule_dwell_times(self, trace, **params):
#        @wraps (self.find_molecule_dwell_times)
#        def inner(mol):
#            if mol.dwell_times(trace, **params) != {}:  # this is where the dwell-time calculation is performed
#                    dwell_times = (mol.dwell_times(trace, **params)["dwell_times"])
#            else:
#                    dwell_times = []
#            return dwell_times
##        inner.__module__ = "__main__"
#        return inner

#    def find_molecule_dwell_times(self,mol, trace, **params):
#        if mol.dwell_times(trace, **params) != {}:  # this is where the dwell-time calculation is performed
#                dwell_times = (mol.dwell_times(trace, **params)["dwell_times"])
#        else:
#                dwell_times = []
#        return dwell_times
#
#
#    def find_dwell_times_mp(self, trace, **params):
#        pool = mp.Pool(mp.cpu_count() - 1)
##        finders = [self.find_molecule_dwell_times(mol) for mol in self.molecules]
#        func = self.find_molecule_dwell_times(mol, trace, **params)
##        func.module__ = "__main__"
#        trace_list = [trace for i in len(self.molecules)]
#        params_list = [params for i in len(self.molecules)]
#        dwell_times = pool.starmap(func, (self.molecules, trace_list, params_list))
#        pool.stop()
#        pool.join()
#        return dwell_times



class Molecule:

    def __init__(self, file, exposure_time):
        self.file = file
        self.coordinates = None
        self.intensity = None
        self.exposure_time = exposure_time

        self.isSelected = False

        self.dwell_times_auto = {"red":{}, "green":{},
                                 "E":{}}

    def I(self, emission):
        return self.intensity[emission, :]

    def E(self):
        return self.I(1) / ( self.I(0) + self.I(1) )

    def plot(self):
        plt.plot(self.intensity[0,:], 'g')
        plt.plot(self.intensity[1,:], 'r')
        plt.show()


    def __apply_to_trace(self, function):
        def inner (trace, **params):
            if trace in ["green"]:
                return function(self.intensity[0, :], exposure_time=self.exposure_time, **params)
            if trace in ["red"]:
                return function(self.intensity[1, :], exposure_time=self.exposure_time, **params)
            if trace in ["E"]:
                return function(self.E(), exposure_time=self.exposure_time, **params)
            else:
                print("You didn't input correct key for trace")

        return inner

    def find_steps(self, trace, **params):
        dwells = self.__apply_to_trace(stepfinder)(trace, **params)
        self.dwell_times_auto[trace] = dwells
        return dwells

    def plot_trace(self, trace, **params):
        return self.__apply_to_trace(plot_steps)(trace, **params)


def histogram(input, axis, makeFit = False):
    if not input: return None
    if not axis: axis = plt.gca()
    #    if not isinstance(input,list): input = [input]
#
#    molecules = list()
#
#    for i in input:
#        if isinstance(i, Molecule):
#            molecules.append(i)
#        else:
#            molecules.append(i.molecules)
    molecules = input
    data = np.concatenate([molecule.E() for molecule in molecules])
    axis.hist(data,100, range = (0,1))

    if makeFit:
        fit(data, axis)


def fit(data, axis):

    hist, bin_edges = np.histogram(data,100, range = (0,1))
    bin_centers = (bin_edges[0:-1]+bin_edges[1:])/2

    #plt.plot(bin_centers,hist)

    from scipy.signal import butter
    from scipy.signal import filtfilt
    b, a = butter(2, 0.2, 'low')
    output_signal = filtfilt(b, a, hist)
    plt.plot(bin_centers, output_signal)

    from scipy.signal import find_peaks
    peaks, properties = find_peaks(output_signal, prominence=5, width=7) #prominence=1
    plt.plot(bin_centers[peaks], hist[peaks], "x")


    def func(x, a, b, c, d, e, f):
        return a * np.exp(-(x-b)**2/(2*c**2)) + d * np.exp(-(x-e)**2/(2*f**2))

    from scipy.optimize import curve_fit
    popt, pcov = curve_fit(func, bin_centers, hist, method = 'trf',
                           p0 = [hist[peaks[0]],bin_centers[peaks[0]],0.1,hist[peaks[1]],bin_centers[peaks[1]],0.1], bounds = (0,[np.inf,1,1,np.inf,1,1]))

    axis.plot(bin_centers,func(bin_centers, *popt))
    #plt.plot(bin_centers,func(bin_centers, 10000,0.18,0.1,5000,0.5,0.2))

#uniqueFileNames = list(set([re.search('hel[0-9]*',fileName).group() for fileName in fileNames]))

"""

for root, dirs, files in os.walk('..', topdown=False):
	for name in files:
		print(os.path.join(root, name))
		print(root)
		print(name)

	for name in dirs:
		print(os.path.join(root, name))



for root, dirs, files in os.walk('..', topdown=False):
	print(files)

"""
