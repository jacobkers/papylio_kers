# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 15:24:46 2018

@author: ivoseverins
"""

import os
import re # Regular expressions
import warnings
import numpy as np
import matplotlib.pyplot as plt

class Experiment:
    def __init__(self, mainPath):
        self.name = os.path.basename(mainPath)
        self.mainPath = mainPath
        self.files = list()
        self.Ncolours = 2
    		
        os.chdir(mainPath)
        
        self.addAllFilesInMainPath()
		
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
        
        fileNameSearch = re.search('(hel[0-9]*)(.*\..*)',fileName)
                
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
        else:
            warnings.warn('FileName ' + fileName + ' not added, filename should contain hel...')
	

    def histogram(self):
        histogram([molecule for file in self.files for molecule in file.molecules])	
	


class File:
    def __init__(self, relativePath, name, experiment):
        self.relativePath = relativePath
        self.name = name
        self.extensions = list()
        self.experiment = experiment
        self.molecules = list()
        
    @property
    def fullPath(self):
        return os.path.join(self.relativePath, self.name)
    
    @property
    def coordinates(self):
        return np.concatenate([[molecule.coordinates[0,:] for molecule in self.molecules]])
    
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
        
        for trace in range(0, Ntraces, Ncolours):
            coordinates = pks[trace:(trace+Ncolours),1:3]
            self.addMolecule(coordinates)
    
    def importTracesFile(self):
        Ncolours = self.experiment.Ncolours
        
        file = open(self.name + '.traces', 'r')
        
        self.Nframes = np.fromfile(file, dtype = np.int32, count = 1).item()
        Ntraces = np.fromfile(file, dtype = np.int16, count = 1).item()
        
        rawData = np.fromfile(file, dtype = np.int16, count = self.Nframes * Ntraces)
        orderedData = np.reshape(rawData.ravel(), (Ncolours, Ntraces//Ncolours, self.Nframes), order = 'F')
        
        for i, molecule in enumerate(self.molecules):
            molecule.intensity = orderedData[:,i,:]
            
            
    
    def addMolecule(self, coordinates):
        self.molecules.append(Molecule(coordinates, self))
        
    
    def histogram(self):
        histogram(self.molecules)
    
    
		
class Molecule:
    def __init__(self, coordinates, file):
        self.coordinates = coordinates
        self.file = file
        self.intensity = None
    
    def plot(self):
        plt.plot(self.intensity[0,:], 'g')
        plt.plot(self.intensity[1,:], 'r')
        plt.show()




def histogram(molecules):
    data = np.concatenate([molecule.intensity[0,:] for molecule in molecules])
    
    plt.hist(data, 100)



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