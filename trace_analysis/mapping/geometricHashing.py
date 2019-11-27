# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 08:59:08 2019

@author: ivoseverins
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle
from pathlib2 import Path
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.pyplot import cm
import scipy.stats
from trace_analysis.mapping.mapping import Mapping2
from trace_analysis.coordinate_transformations import translate, rotate, magnify, reflect, transform
from trace_analysis.plotting import plot_match


def mapToPoint(pointSet,startPoints,endPoints,returnTransformationMatrix = False, tr = None, di = None, ro = None):
    startPoints = np.atleast_2d(startPoints)
    endPoints = np.atleast_2d(endPoints)
    if len(startPoints) == 1 & len(endPoints) == 1: 
        tr = True; ro = False; di = False
       
    elif len(startPoints) == 2 & len(endPoints) == 2:
        if tr is None: tr = True
        if di is None: di = True
        if ro is None: ro = True
    
    transformationMatrix = np.identity(3)
    
    if tr:
        translationMatrix = translate(endPoints[0] - startPoints[0])
        transformationMatrix = translationMatrix @ transformationMatrix
    
    if di or ro:
        diffs = np.array([startPoints[0]-startPoints[1],endPoints[0]-endPoints[1]])
        diffLengths = np.linalg.norm(diffs,axis=1,keepdims=True)
        unitDiffs = diffs/diffLengths
        
        if di:
            dilationMatrix = magnify(diffLengths[1]/diffLengths[0], endPoints[0])
            transformationMatrix = dilationMatrix @ transformationMatrix
    
        if ro:
            angle = -np.arctan2(np.linalg.det(unitDiffs),np.dot(unitDiffs[0],unitDiffs[1]))
            #angle = np.arccos(np.dot(diffs[0]/endLength,diffs[1]/startLength))
            rotationMatrix = rotate(angle,endPoints[0])
            transformationMatrix = rotationMatrix @ transformationMatrix
    
    
    pointSet = np.append(pointSet,np.ones((pointSet.shape[0],1)),axis=1)
    transformedPointSet = (transformationMatrix @ pointSet.T)[0:2,:].T
    
    if returnTransformationMatrix:
        return transformedPointSet, transformationMatrix
    else:
        return transformedPointSet

def plotPointSet(pointSet1, pointSet2 = None):
    fig, ax = plt.subplots()
    #ax.cla()
   
    #colour=iter(cm.rainbow(np.linspace(0,1,len(pointSet2))))
    
    ax.scatter(pointSet1[:,0],pointSet1[:,1], marker = 'o', facecolors = 'none', edgecolors='b')
    if pointSet2 is not None:
        if not isinstance(pointSet2, list):
            ax.scatter(pointSet2[:,0],pointSet2[:,1],c='r',marker = 'x')
        else:
            for pointSet in pointSet2:
                #c = next(colour)
                ax.scatter(pointSet[:,0],pointSet[:,1],marker = 'x')


def pointHash(pointSet, bases='all', mode='similarity', hashTableRange = np.array([-1,1]), nBins = 100,
              rotationRange = None, magnificationRange = None):
#    basis0, basis1 = np.mgrid[0:pointSet.shape[0],0:pointSet.shape[0]]
#    basis0 = basis0.flatten()
#    basis1 = basis1.flatten()
#
    #pointSet = similarityTransformation(pointSet)
    
    #nBins = 200
    #nBins = 100

    binEdges = np.linspace(hashTableRange[0],hashTableRange[1],nBins+1)
    
    
    hashTable = {}
    for binX in (np.arange(nBins)+1):
        for binY in (np.arange(nBins)+1):
            hashTable[(binX,binY)] = []
  
    
#    np.histogram2d(normalizedPointSet[:,0],normalizedPointSet[:,1],[binEdges,binEdges])
#    
#    statistic, xedges, yedges, binnumber = scipy.stats.binned_statistic_2d(
#             normalizedPointSet[:,0],
#             normalizedPointSet[:,1],
#             values=np.ones(normalizedPointSet.shape[0]),
#             bins=[binEdges,binEdges])
#    
    if mode == 'similarity':
        endPoints = np.array([[-0.5,0],[0.5,0]])
        NbasisPoints = 2
    elif mode == 'translation':
        endPoints = np.array([[0,0]])
        NbasisPoints = 1
    
    
    if bases is 'all':
        #bases = np.arange(testPointSet.shape[0])
        
        baseNumbersPerAxis = np.mgrid[[slice(pointSet.shape[0]) for i in range(NbasisPoints)]]
        bases = np.column_stack([numbersOnAxis.flatten() for numbersOnAxis in baseNumbersPerAxis])
        bases = bases[[len(np.unique(basis)) == NbasisPoints for basis in bases]]
    
    print('0 of ' + str(len(bases)) + ' bases')
    
    for i, base in enumerate(bases):
        
        ## This part is temporary
        if (mode == 'similarity') and ((rotationRange is not None) or (magnificationRange is not None)):
            startPoints = pointSet[base]
            diffs = np.array([startPoints[0]-startPoints[1],endPoints[0]-endPoints[1]])
            diffLengths = np.linalg.norm(diffs,axis=1,keepdims=True)
            unitDiffs = diffs/diffLengths

            magnification = diffLengths[1] / diffLengths[0]
            rotation = -np.arctan2(np.linalg.det(unitDiffs), np.dot(unitDiffs[0], unitDiffs[1]))

            if not ((rotationRange[0] < rotation < rotationRange[1]) and \
                    (magnificationRange[0] < 1/magnification < magnificationRange[1])):
                #print('skipped')
                continue

        # dC = pointSet[base[1]] - pointSet[base[0]]
        # length = np.linalg.norm(dC,axis=0,keepdims=True)
        # #print(length)
        #
        # if length < 1500: continue
        # if length > 2500: continue
        
        ## Part above is temporary
        
        normalizedPointSet = mapToPoint(pointSet,pointSet[base],endPoints)
        
        #normalizedPointSet = similarityTransformation(normalizedPointSet)
        
        normalizedPointSet = np.round(normalizedPointSet,10)
        #print(normalizedPointSet)
        
        normalizedPointSetWithoutBasisPoints = np.delete(normalizedPointSet,base,0)
        
        #normalizedPointSetWithoutBasisPoints = similarityTransformation(normalizedPointSetWithoutBasisPoints)
        
        bins = np.digitize(normalizedPointSetWithoutBasisPoints,binEdges)
        #print(bins)
        bins = bins[~np.any(np.logical_or(bins == 0, bins == nBins+1),axis=1)]
        if bins.size > 0:
            bins = np.unique(bins, axis=0)
            #print(bins.shape)
            
            for aBin in bins:
                hashTable[tuple(aBin)].append(base)
        #else:
            #print('empty')
            
        print(str(i) + ' of ' + str(len(bases)) + ' bases')    
        
    plt.figure()
    hashTableImage = np.zeros((nBins,nBins))
    for key, value in hashTable.items():
        hashTableImage[key[0]-1,key[1]-1] = len(value)
    plt.imshow(np.rot90(hashTableImage), vmin=0)
    
    fig = plt.figure()
    X = np.arange(nBins)
    Y = np.arange(nBins)
    X, Y = np.meshgrid(X, Y)
    Z = hashTableImage
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X,Y,Z)
    ax.set_zlim(0, 500)
    
    return hashTable



def findMatch(testPointSet, hashTable, bases = 'all', returnMatchedBases = False, 
              mode='similarity', hashTableRange = np.array([-1,1]), nBins = 100,
              rotationRange = None, magnificationRange = None):
    
    #nBins = 200
    #nBins = 100
    
    binEdges = np.linspace(hashTableRange[0],hashTableRange[1],nBins+1)

    matchedBases = []
    allBestBaseMatches = {'hashTableBases': [], 'testBases': [], 'counts': []}
    bestBasis = {'hashTableBasis': [], 'testBasis': [], 'counts': 0}
    
    if mode == 'similarity':
        endPoints = np.array([[-0.5,0],[0.5,0]])
        NbasisPoints = 2
    elif mode == 'translation':
        endPoints = np.array([[0,0]])
        NbasisPoints = 1
    
    if bases is 'all':
        #bases = np.arange(testPointSet.shape[0])
        
        baseNumbersPerAxis = np.mgrid[[slice(testPointSet.shape[0]) for i in range(NbasisPoints)]]
        bases = np.column_stack([numbersOnAxis.flatten() for numbersOnAxis in baseNumbersPerAxis])
        bases = bases[[len(np.unique(basis)) == NbasisPoints for basis in bases]]
    
    elif isinstance(bases, int):
        Nbases = bases
        baseNumbersPerAxis = np.mgrid[[slice(testPointSet.shape[0]) for i in range(NbasisPoints)]]
        bases = np.column_stack([numbersOnAxis.flatten() for numbersOnAxis in baseNumbersPerAxis])
        bases = bases[[len(np.unique(basis)) == NbasisPoints for basis in bases]]
        np.random.shuffle(bases)
        bases = bases[:Nbases]
    
    for base in bases:

        if (mode == 'similarity') and ((rotationRange is not None) or (magnificationRange is not None)):
            startPoints = testPointSet[base]
            diffs = np.array([startPoints[0]-startPoints[1],endPoints[0]-endPoints[1]])
            diffLengths = np.linalg.norm(diffs,axis=1,keepdims=True)
            unitDiffs = diffs/diffLengths

            magnification = diffLengths[1] / diffLengths[0]
            rotation = -np.arctan2(np.linalg.det(unitDiffs), np.dot(unitDiffs[0], unitDiffs[1]))

            if not ((rotationRange[0] < rotation < rotationRange[1]) and \
                    (magnificationRange[0] < 1/magnification < magnificationRange[1])):
                continue

        normalizedPointSet = mapToPoint(testPointSet,testPointSet[base],endPoints)
        
        #normalizedPointSet = similarityTransformation(normalizedPointSet)
        
        normalizedPointSet = np.round(normalizedPointSet,10)
        #print(normalizedPointSet)
        normalizedPointSetWithoutBasisPoints = np.delete(normalizedPointSet,base,0)
        
        bins = np.digitize(normalizedPointSetWithoutBasisPoints,binEdges)
        #print(bins)
        bins = bins[~np.any(np.logical_or(bins == 0, bins == nBins+1),axis=1)]
        
        #print(bins)
        
        foundBases = [element for aBin in bins for element in hashTable[tuple(aBin)]]
        
        if foundBases:
            uniqueBases, counts = np.unique(foundBases, return_counts=True, axis = 0)
            
            matchedBases.append([base,uniqueBases,counts])
            
            
            allBestBaseMatches['hashTableBases'].append(uniqueBases[np.argmax(counts)])
            allBestBaseMatches['testBases'].append(np.array(base))
            allBestBaseMatches['counts'].append(np.max(counts))
            
            currentBestBasis = {'hashTableBasis': uniqueBases[np.argmax(counts)],
                                                              'testBasis': np.array(base),
                                                              'counts': np.max(counts)}
            if currentBestBasis['counts'] > bestBasis['counts']:
                bestBasis = currentBestBasis
                print('New best basis')
                print(bestBasis['counts'])
        else:
            allBestBaseMatches['hashTableBases'].append(None)
            allBestBaseMatches['testBases'].append(np.array(base))
            allBestBaseMatches['counts'].append(0)
            #print('No match found')
        
    if returnMatchedBases:
        return bestBasis, matchedBases, allBestBaseMatches
    else:
        return bestBasis


def similarityTransformation(pointSet):
    u = pointSet[:,0]
    v = pointSet[:,1]
    
    
    rho = (0.5-3/(4*(u**2+v**2)+3))*2
    #rho = (0.5-3/(2*(u**2+v**2)+3))*2
    phi = np.arctan2(v,u)/np.pi
    
    return np.array([rho,phi]).T
    
#    x = rho * np.cos(phi)
#    y = rho * np.sin(phi)
#    
#    return np.array([x,y]).T



def pointHashTranslation(pointSet, bases='all'):
#    basis0, basis1 = np.mgrid[0:pointSet.shape[0],0:pointSet.shape[0]]
#    basis0 = basis0.flatten()
#    basis1 = basis1.flatten()
#
    #pointSet = similarityTransformation(pointSet)
    pointSet = np.append(pointSet,np.ones((pointSet.shape[0],1)),axis=1)
    
    nBins = 100
    
    binEdges = np.linspace(-10000,10000,nBins+1)
    
    
    hashTable = {}
    for binX in (np.arange(nBins)+1):
        for binY in (np.arange(nBins)+1):
            hashTable[(binX,binY)] = []

    #endPoints = np.array([[-0.5,0],[0.5,0]])
    
    if bases is 'all':
        bases = np.arange(pointSet.shape[0])
        
#        basesX, basesY = np.mgrid[0:pointSet.shape[0],0:pointSet.shape[0]]
#        bases = np.column_stack((basesX.flatten(),basesY.flatten()))
#        bases = bases[~(bases[:,0]==bases[:,1])]
    
    for base in bases:
      
        #normalizedPointSet = mapToPoint(pointSet,pointSet[base],endPoints)
        
        normalizedPointSet = (translate(-pointSet[base]) @ pointSet.T)[0:2,:].T
        
        
        normalizedPointSet = np.round(normalizedPointSet,10)
        #print(normalizedPointSet)
        
        normalizedPointSetWithoutBasisPoints = np.delete(normalizedPointSet,base,0)
        
        #normalizedPointSetWithoutBasisPoints = similarityTransformation(normalizedPointSetWithoutBasisPoints)
        
        bins = np.digitize(normalizedPointSetWithoutBasisPoints,binEdges)
        #print(bins)
        bins = bins[~np.any(np.logical_or(bins == 0, bins == nBins+1),axis=1)]
        if bins.size > 0:
            bins = np.unique(bins, axis=0)
            #print(bins.shape)
            
            for aBin in bins:
                hashTable[tuple(aBin)].append(base)
        #else:
            #print('empty')

    hashTableImage = np.zeros((nBins,nBins))
    for key, value in hashTable.items():
        hashTableImage[key[0]-1,key[1]-1] = len(value)
    plt.imshow(np.rot90(hashTableImage), vmin=0)
    
    fig = plt.figure()
    X = np.arange(nBins)
    Y = np.arange(nBins)
    X, Y = np.meshgrid(X, Y)
    Z = hashTableImage
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X,Y,Z)
    ax.set_zlim(0, 10000)
    
    return hashTable



def findMatchTranslation(testPointSet, hashTable, bases = 'all', returnMatchedBases = False):
    nBins = 200
    
    binEdges = np.linspace(-10000,10000,nBins+1)
    
    testPointSet = np.append(testPointSet,np.ones((testPointSet.shape[0],1)),axis=1)

    matchedBases = []
    allBestBaseMatches = {'hashTableBases': [], 'testBases': [], 'counts': []}
    bestBasis = {'hashTableBasis': [], 'testBasis': [], 'counts': 0}
    
    endPoints = np.array([[-0.5,0],[0.5,0]])
    
    if bases is 'all':
        bases = np.arange(testPointSet.shape[0])

    
    elif isinstance(bases, int):
        Nbases = bases
        bases = np.arange(testPointSet.shape[0])
        np.random.shuffle(bases)
        bases = bases[:Nbases]
    
    for base in bases:
        #normalizedPointSet = mapToPoint(testPointSet,testPointSet[base],endPoints)
        normalizedPointSet = (translate(-testPointSet[base]) @ testPointSet.T)[0:2,:].T
        
        normalizedPointSet = np.round(normalizedPointSet,10)
        #print(normalizedPointSet)
        normalizedPointSetWithoutBasisPoints = np.delete(normalizedPointSet,base,0)
        
        bins = np.digitize(normalizedPointSetWithoutBasisPoints,binEdges)
        #print(bins)
        bins = bins[~np.any(np.logical_or(bins == 0, bins == nBins+1),axis=1)]
        
        #print(bins)
        
        foundBases = [element for aBin in bins for element in hashTable[tuple(aBin)]]
        
        if foundBases:
            uniqueBases, counts = np.unique(foundBases, return_counts=True, axis = 0)
            
            matchedBases.append([base,uniqueBases,counts])
            
            
            allBestBaseMatches['hashTableBases'].append(uniqueBases[np.argmax(counts)])
            allBestBaseMatches['testBases'].append(np.array(base))
            allBestBaseMatches['counts'].append(np.max(counts))
            
            currentBestBasis = {'hashTableBasis': uniqueBases[np.argmax(counts)],
                                                              'testBasis': np.array(base),
                                                              'counts': np.max(counts)}
            if currentBestBasis['counts'] > bestBasis['counts']:
                bestBasis = currentBestBasis
                print('New best basis')
                print(bestBasis['counts'])
        else:
            allBestBaseMatches['hashTableBases'].append(None)
            allBestBaseMatches['testBases'].append(np.array(base))
            allBestBaseMatches['counts'].append(0)
            #print('No match found')
        
    if returnMatchedBases:
        return bestBasis, matchedBases, allBestBaseMatches
    else:
        return bestBasis


class SequencingDataMapping:
    def __init__(self, tile, files, mode, dataPath, nBins=100, rotationRange=None, magnificationRange=None):
        self.tile = tile
        self.files = files # List of coordinate sets
        self.dataPath = Path(dataPath)

        self.initial_image_transformation = {'reflection': 0, 'rotation': 0, 'magnification': 1} # This reflects with respect to axis 0.

        self.mode = mode
        if mode == 'translation': self.hashTableRange = [-10000, 10000]
        else: self.hashTableRange = [-1,1]
        self.nBins = nBins

        self.bases_hashTable = 'all'
        self.bases_findMatch = 20

        self._hashTable = None
        self._matches = None

        self.rotationRange = rotationRange
        self.magnificationRange = magnificationRange


    @property
    def hashTable(self):
        if self._hashTable is None:
            if self.dataPath.joinpath(self.tile.name+'.ht').is_file():
                with self.dataPath.joinpath(self.tile.name+'.ht').open('rb') as f:
                    self._hashTable = pickle.load(f)
            else:
                self._hashTable = pointHash(self.tile.coordinates, bases=self.bases_hashTable, mode=self.mode,
                                            hashTableRange=self.hashTableRange, nBins=self.nBins,
                                            rotationRange=self.rotationRange, magnificationRange=self.magnificationRange)
                with self.dataPath.joinpath(self.tile.name+'.ht').open('wb') as f:
                    pickle.dump(self._hashTable, f)
        return self._hashTable

    @property
    def matches(self):
        if self._matches is None:
            if self.dataPath.joinpath(self.tile.name+'.matches').is_file():
                with self.dataPath.joinpath(self.tile.name+'.matches').open('rb') as f:
                    self._matches = pickle.load(f)
            else:
                self._matches = self.findMatches()
                with self.dataPath.joinpath(self.tile.name+'.matches').open('wb') as f:
                    pickle.dump(self._matches, f)
        return self._matches

    def save(self):
        if self.dataPath.joinpath(self.tile.name+'.sm').is_file():
            print('File already exists')
        else:
            with self.dataPath.joinpath(self.tile.name+'.sm').open('wb') as f:
                pickle.dump(self, f)

    @staticmethod
    def load(path):
        with path.open('rb') as f:
            return pickle.load(f)

    def findMatches(self):
        matches = []
        for file in self.files:
            print(file.name)
            coordinates, initial_transformation = transform(file.coordinates, returnTransformationMatrix = True,
                                                            **self.initial_image_transformation)

            [bestBasis, matchedBases, allBestBaseMatches] = findMatch(coordinates, self.hashTable,
                                                                      bases=self.bases_findMatch,
                                                                      returnMatchedBases=True,
                                                                      mode=self.mode,
                                                                      nBins=self.nBins,
                                                                      hashTableRange=self.hashTableRange,
                                                                      rotationRange=self.rotationRange,
                                                                      magnificationRange=self.magnificationRange)

            coordinates_transformed, transformation_matrix = \
                mapToPoint(coordinates, coordinates[bestBasis['testBasis']],
                            self.tile.coordinates[bestBasis['hashTableBasis']],
                            returnTransformationMatrix=True)

            transformation_matrix = transformation_matrix @ initial_transformation

            match = Mapping2(file.coordinates.copy(), self.tile.coordinates.copy(),
                                  method='geometric-hashing',
                                  transformation_type='linear')

            # match.name =
            match.transformation = transformation_matrix
            match.best_image_basis = bestBasis
            match.count = bestBasis['counts']
            #match.image_coordinates_transformed = coordinates_transformed
            match.meanError = np.mean(
                [np.min(np.linalg.norm(self.tile.coordinates - row, axis=1)) for row in coordinates_transformed])
            match.setWidth = np.linalg.norm(np.max(coordinates_transformed, axis=0) - np.min(coordinates_transformed, axis=0))
            match.percentMatch = bestBasis['counts'] / coordinates_transformed.shape[0]

            matches.append(match)

        return matches

    def histogram_matches(self, export = False):
        counts = [match.count for match in self.matches]
        fig, ax = plt.subplots(figsize = (6,3))
        ax.hist(counts, np.arange(0.5, np.max(counts) + 1.5, 1), facecolor = 'k', rwidth = 0.5)

        ax.set_xlim([0,np.max(counts)+1])

        ax.set_xlabel('Number of matches for best matching basis')
        ax.set_ylabel('Count')
        plt.tight_layout()

        if export:
            fig.savefig(self.dataPath.joinpath('histogramNumberOfMatches.pdf'), bbox_inches='tight')
            fig.savefig(self.dataPath.joinpath('histogramNumberOfMatches.png'), bbox_inches='tight')

    def show_match(self, match, figure = None, view='destination'):
        if not figure: figure = plt.gcf()
        figure.clf()
        ax = figure.gca()

        #ax.scatter(ps2[:,0],ps2[:,1],c='g',marker = '+')

        ax.scatter(match.destination[:,0],match.destination[:,1], marker = '.', facecolors = 'k', edgecolors='k')
        ax.scatter(match.transform_source_to_destination[:,0],match.transform_source_to_destination[:,1],c='r',marker = 'x')

        destination_basis_index = match.best_image_basis['hashTableBasis']
        source_basis_index = match.best_image_basis['testBasis']
        ax.scatter(match.destination[destination_basis_index, 0], match.destination[destination_basis_index, 1], marker='.', facecolors='g', edgecolors='g')
        ax.scatter(match.transform_source_to_destination[source_basis_index, 0], match.transform_source_to_destination[source_basis_index, 1], c='g',
                   marker='x')

        ax.set_aspect('equal')
        ax.set_title('Tile:' + self.tile.name +', File: ' + str(self.files[self.matches.index(match)].relativeFilePath))

        if view == 'source':
            maxs = np.max(match.transform_source_to_destination, axis=0)
            mins = np.min(match.transform_source_to_destination, axis=0)
            ax.set_xlim([mins[0], maxs[0]])
            ax.set_ylim([mins[1], maxs[1]])
        elif view == 'destination':
            maxs = np.max(match.destination, axis=0)
            mins = np.min(match.destination, axis=0)
            ax.set_xlim([mins[0], maxs[0]])
            ax.set_ylim([mins[1], maxs[1]])
            # ax.set_xlim([0, 31000])
            # ax.set_ylim([0, 31000])

        name = str(self.files[self.matches.index(match)].relativeFilePath)
        print(name)
        n = name.replace('\\', '_')

        figure.savefig(self.dataPath.joinpath(n + '_raw.pdf'), bbox_inches='tight')
        figure.savefig(self.dataPath.joinpath(n + '_raw.png'), bbox_inches='tight', dpi=1000)

    def plot_match(self, match):
        name = str(self.files[self.matches.index(match)].relativeFilePath)
        print(name)
        plot_match(match, self.dataPath, name, unit='um')




#fig, ax = plt.subplots()
#
#nBins = 50
#lims = np.array([-1,1])
#binEdges = np.linspace(lims[0],lims[1],nBins+1)
#
#from matplotlib.ticker import MultipleLocator
#ax.cla()
#endPoints = np.array([[-0.5,0],[0.5,0]])
#nps1 = mapToPoint(ps1,ps1[[5,31]],endPoints)
#nps2 = mapToPoint(ps2,ps2[[0,1]],endPoints)
#
#nps1 = similarityTransformation(nps1)
#nps2 = similarityTransformation(nps2)
#
#ax.scatter(nps1[:,0],nps1[:,1],c='g',marker = '+')
#ax.scatter(nps2[:,0],nps2[:,1],c='k',marker = 'x')
#
#ax.set_xlim([-2,2])
#ax.set_ylim([-2,2])
#ax.set_xticks(binEdges, minor = True)
#ax.set_yticks(binEdges, minor = True)
#
#plt.grid(True, which='minor')
#
#
#fig, ax = plt.subplots()
#for t in nps1:
#    plt.polar(t[1],t[0],'ro')
#    
#plt.show()

    
#    
#ps1 = np.array([[1,1],[1,10],[10,1],[10,10]])
#ps2 = transform(ps1, t=[4,2],r=np.pi/4, m=1.2)

#if __name__== "__main__":
#    ps1 = S
#    ht = pointHash(ps1)
#    
#    
#    ps2o = M
#    ps2 = M
#    #ps2 = M.astype('uint8').astype('float')
#    
#    indices = np.arange(ps2.shape[0])
#    np.random.shuffle(indices)
#    ps2 = ps2[indices[10:]]
#    
#    
#    ps2 = ps2 + np.random.rand(ps2.shape[0],ps2.shape[1])
#    
#    
#    [bestBasis, matchedBases, allBestBaseMatches] = findMatch(ps2, ht, bases = 'all', returnMatchedBases = True)
#    
#    ps3 = mapToPoint(ps2,ps2[bestBasis['testBasis']],ps1[bestBasis['hashTableBasis']])
#    
#    fig, ax = plt.subplots()
#    ax.cla()
#       
#    ax.scatter(ps1[:,0],ps1[:,1], marker = 'o', facecolors = 'none', edgecolors='b')
#    ax.scatter(ps2o[:,0],ps2o[:,1],c='g',marker = '+')
#    ax.scatter(ps2[:,0],ps2[:,1],c='g',marker = 'x')
#    ax.scatter(ps3[:,0],ps3[:,1],c='k',marker = 'x')
#    
#    
#    ax.set_aspect('equal')
#    
#    
#    meanError = []
#    plt.ion()
#    fig, ax = plt.subplots()
#    for i in np.arange(len(allBestBaseMatches['counts'])):
#        if not np.any(allBestBaseMatches['hashTableBases'][i] == None):
#            
#            ps3 = mapToPoint(ps2,ps2[allBestBaseMatches['testBases'][i]],ps1[allBestBaseMatches['hashTableBases'][i]])
#            
#            
#            meanError.append(np.mean([np.min(np.linalg.norm(Sselection - row,axis=1)) for row in ps3]))
#            
#            
#            ax.cla()
#            
#            
#            ax.scatter(ps1[:,0],ps1[:,1], marker = 'o', facecolors = 'none', edgecolors='b')
#            ax.scatter(ps2[:,0],ps2[:,1],c='g',marker = '+')
#            ax.scatter(ps3[:,0],ps3[:,1],c='k',marker = 'x')
#            
#            
#            ax.set_aspect('equal')
#            
#            #fig.canvas.mpl_connect('key_press_event',print('test'))
#            #plt.draw()
#            plt.pause(0.1)
#            plt.waitforbuttonpress()
#            #plt.show()
#            
#            #input()
#            
#    
#    np.sum(np.array(meanError)<1)/len(meanError)
#    
    
    
    
    
    






#endPoints = np.array([[-0.5,0],[0.5,0]])
#fig, ax = plt.subplots()
#ax.cla()
#
#tps1 = similarityTransformation(mapToPoint(ps1,ps1[bestBasis['hashTableBasis']],endPoints))
#tps2 = similarityTransformation(mapToPoint(ps2,ps2[bestBasis['testBasis']],endPoints))
#   
#ax.scatter(tps1[:,0],tps1[:,1], marker = 'o', facecolors = 'none', edgecolors='b')
#ax.scatter(tps2[:,0],tps2[:,1],c='g',marker = '+')
##ax.scatter(ps3[:,0],ps3[:,1],c='k',marker = 'x')
#
#nBins = 50
#lims = np.array([-1,1])
#binEdges = np.linspace(lims[0],lims[1],nBins+1)
#
#
#ax.set_xlim([-2,2])
#ax.set_ylim([-2,2])
#ax.set_xticks(binEdges, minor = True)
#ax.set_yticks(binEdges, minor = True)
#
#plt.grid(True, which='minor')
#
#ax.set_aspect('equal')



#fig, ax = plt.subplots()


