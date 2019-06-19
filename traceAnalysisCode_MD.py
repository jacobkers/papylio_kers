# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 15:24:46 2018

@author: ivoseverins
"""

import os # Miscellaneous operating system interfaces - to be able to switch from Mac to Windows
import re # Regular expressions
import warnings
import numpy as np #scientific computing with Python
import matplotlib.pyplot as plt #Provides a MATLAB-like plotting framework
import itertools #Functions creating iterators for efficient looping

class Experiment:
    def __init__(self, mainPath):
        self.name = os.path.basename(mainPath)
        self.mainPath = mainPath
        self.files = list()
        self.Ncolours = 2
    		
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
                  #  print(root) #suppress printing skipped files and directories
                  #  print(fileName) #suppress printing skipped files and directories
                    
                    self.addFile(root, fileName)
	
            for name in dirs:
                warnings.warn('Import of files in subfolders does not work properly yet')
                print(os.path.join(root, name))
			
    def addFile(self, relativePath, fileName):
        
        #fileNameSearch = re.search('(hel[0-9]*)(.*\..*)',fileName) # here you search for hel + numerics, everything beyond is supposed to be extension
        fileNameSearch = re.search('(hel[0-9].*)((?=.).*\..*)',fileName) # here you search for anything that starts with hel, and you find the extension after the dot
        
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
                self.files.append(File(relativePath, name, self))
                self.files[-1].addExtension(extension)
            else:
                if extension not in foundFile.extensions:
                    foundFile.addExtension(extension)
                else:
                    warnings.warn(fullPath + extension + ' already imported')
        #else: #suppress printing skipped files and directories
          #  warnings.warn('FileName ' + fileName + ' not added, filename should contain hel...')
	

    def histogram(self, axis = None, fileSelection = False, moleculeSelection = False):
        #files = [file for file in exp.files if file.isSelected]
        #files = self.files
        
        if (fileSelection & moleculeSelection): 
            histogram([molecule for file in self.selectedFiles for molecule in file.selectedMolecules], axis)
        elif (fileSelection & (not moleculeSelection)):
            histogram([molecule for file in self.selectedFiles for molecule in file.molecules], axis)
        elif ((not fileSelection) & moleculeSelection):
            histogram([molecule for file in self.files for molecule in file.selectedMolecules], axis)
        else:
            histogram([molecule for file in self.files for molecule in file.molecules], axis)
	


class File:
    def __init__(self, relativePath, name, experiment):
        self.relativePath = relativePath
        self.name = name
        self.extensions = list()
        self.experiment = experiment
        self.molecules = list()
        
        self.isSelected = False
        
    @property
    def fullPath(self):
        return os.path.join(self.relativePath, self.name)
    
    @property
    def coordinates(self):
        return np.concatenate([[molecule.coordinates[0,:] for molecule in self.molecules]])
    
    @property
    def selectedMolecules(self):
        return [molecule for molecule in self.molecules if molecule.isSelected]
    
    def addExtension(self, extension):
        self.extensions.append(extension)
        #print(extension) #suppress printing skipped files and directories
        
        importFunctions = {'.pks'        : self.importPksFile,
                           '.traces'     : self.importTracesFile}
                
        importFunctions.get(extension, print)()
        
        
#        if extension == '.pks':
            #self.importPksFile()

   
    
    def importPksFile(self):
        # Background value stored in pks file is not imported yet
        Ncolours = self.experiment.Ncolours
        
   #     pks = np.genfromtxt(self.name + '.pks', delimiter='      ')  #MD190104 you get an error when using 6 spaces for tab
        pks = np.genfromtxt(self.name + '.pks')  #MD190104 By default, any consecutive whitespaces act as delimiter.
        print(self.name + '.pks') 
        Ntraces = np.shape(pks)[0]
        
        if not self.molecules:
            for molecule in range(0, Ntraces, Ncolours):
                self.addMolecule()
        
#        for trace in range(0, Ntraces, Ncolours):
#            coordinates = pks[trace:(trace+Ncolours),1:3]
#            self.addMolecule(coordinates)
        
        for i, molecule in enumerate(self.molecules):
            molecule.coordinates = pks[(i*Ncolours):((i+1)*Ncolours),1:3]
    
    def importTracesFile(self):
        Ncolours = self.experiment.Ncolours
        
        file = open(self.name + '.traces', 'r')
        
        self.Nframes = np.fromfile(file, dtype = np.int32, count = 1).item()
        
        Ntraces = np.fromfile(file, dtype = np.int16, count = 1).item()
        
        print(self.name + '.traces  ' + str(self.Nframes) +' ' + str(Ntraces ))
        
        if not self.molecules:
            for molecule in range(0, Ntraces, Ncolours):
                self.addMolecule()
        
        rawData = np.fromfile(file, dtype = np.int16, count = self.Nframes * Ntraces)
        orderedData = np.reshape(rawData.ravel(), (Ncolours, Ntraces//Ncolours, self.Nframes), order = 'F')
        
        for i, molecule in enumerate(self.molecules):
            molecule.intensity = orderedData[:,i,:]
            
            
    def addMolecule(self):
        self.molecules.append(Molecule(self))
        
    
    def histogram(self, axis=None):
        histogram(self.molecules, axis)
    
    
		
class Molecule:
    def __init__(self, file):
        self.file = file
        self.coordinates = None
        self.intensity = None
        
        self.isSelected = False
    

    def I(self, emission):
        return self.intensity[emission,:]
    
    def E(self):
        # A= self.I(1) / ( self.I(0) + self.I(1) )
        with np.errstate(divide='ignore'): #divide errors are ignored, taken care of with np.where
            #A= np.where( (self.I(0) + self.I(1)) ==0 ,  0,  self.I(1) / ( self.I(0) + self.I(1) ))
            A= np.where( (np.array(self.I(0)) + np.array(self.I(1))) ==0 ,  0,  np.array(self.I(1)) / ( np.array(self.I(0)) + np.array(self.I(1)) ))
            return A
        
    
    def plot(self):
        plt.plot(self.intensity[0,:], 'g')
        plt.plot(self.intensity[1,:], 'r')
        plt.show()
#MD190104: why not add a subplot with FRET here as well, to match with display Matlab?

    #def dwelltime


def histogram(input, axis):
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
    
    #data = np.concatenate([molecule.intensity[0,:] for molecule in molecules])
    data = np.concatenate([molecule.E() for molecule in molecules])
    
    #plt.hist(data, 100)
    axis.hist(data,100, range = (0,1))


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