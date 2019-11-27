# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 13:55:53 2019

@author: ivoseverins
"""

import numpy as np
import matplotlib.pyplot as plt

from trace_analysis.mapping.icp import icp
from trace_analysis.mapping.icp_nonrigid import icp_nonrigid
from trace_analysis.coordinate_transformations import transform
from trace_analysis.plotting import plot_match

class Mapping2:
    def __init__(self, source = None, destination = None, method = None,
                 transformation_type = 'linear', initial_translation = None):
        self.source = source
        self.destination = destination
        self.method = method
        self.transformation_type = transformation_type
        self.transformation = None
        self.transformation_inverse = None
        self.initial_translation = initial_translation

        # if (source is not None) and (destination is not None):
        #     if self.method is None: self.method = 'icp'
        #     self.perform_mapping()

    @property
    def translation(self):
        return self.transformation[0:2,2]

    @property
    def magnification(self):
        return np.linalg.norm(self.transformation[:,0:2],axis=0)

    @property
    def rotation(self):
        rotation_matrix = self.transformation[0:2, 0:2]/self.magnification[0]
        return np.arctan2(rotation_matrix[0,1]-rotation_matrix[1,0],rotation_matrix[0,0]+rotation_matrix[1,1])/(2*np.pi)*360

    @property
    def transform_source_to_destination(self):
        return self.transform_coordinates(self.source)

    def perform_mapping(self):
        if self.method == 'icp':
            self.transformation, distances, iterations = icp(self.source, self.destination, initial_translation=self.initial_translation)
        elif self.method == 'icp-non-rigid':
            self.transformation, distances, iterations = icp_nonrigid(self.source, self.destination)
        # elif method == 'manual'         : mapping_manual(source, destination)
        # elif method == 'automatic'      : mapping_automatic(source, destination)
        else: print('Method not found')

        self.transformation_inverse = np.linalg.inv(self.transformation)

    def show_mapping_transformation(self, figure=None):
        if not figure: figure = plt.figure()

        destination_from_source = self.transform_coordinates(self.source)

        plt.scatter(self.source[:, 0], self.source[:, 1], c='b')
        plt.scatter(self.destination[:, 0], self.destination[:, 1], c='r')
        plt.scatter(destination_from_source[:, 0], destination_from_source[:, 1], c='g')

    def transform_coordinates(self, coordinates, inverse = False):
        if self.transformation_type == 'linear':
            if inverse is False: return transform(coordinates, self.transformation)
            elif inverse is True: return transform(coordinates, self.transformation_inverse)
        else: print('Transformation not found')


# from trace_analysis.icp_nonrigid import icp_nonrigid
# from trace_analysis.image_adapt.polywarp import polywarp_apply
# kx, ky, distances, i = icp_nonrigid(donor,acceptor, tolerance=0.00000001, max_iterations=50)
# acceptor_calculated = polywarp_apply(kx, ky, donor)

from trace_analysis.mapping.geometricHashing import pointHash, findMatch, mapToPoint
import pickle
import json
from pathlib2 import Path

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
            match.image_coordinates_transformed = coordinates_transformed
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

