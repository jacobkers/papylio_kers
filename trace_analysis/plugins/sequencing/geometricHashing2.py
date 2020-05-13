import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import KDTree
import random
from trace_analysis.mapping.geometricHashing import mapToPoint
from trace_analysis.plotting import scatter_coordinates
from skimage.transform import AffineTransform
import time

# Make random source and destination dataset
np.random.seed(42)
destination = np.random.rand(1000,2)*1000
source_bounds_in_destination = np.array([[300, 300], [450, 600]])



def crop_coordinates(coordinates, bounds):
    bounds.sort(axis=0)
    selection = (coordinates[:, 0] > bounds[0, 0]) & (coordinates[:, 0] < bounds[1, 0]) & \
                (coordinates[:, 1] > bounds[0, 1]) & (coordinates[:, 1] < bounds[1, 1])
    return coordinates[selection]

#
# selection = (destination[:,0] > source_bounds[0,0]) & (destination[:,0] < source_bounds[1,0]) & \
#             (destination[:,1] > source_bounds[0,1]) & (destination[:,1] < source_bounds[1,1])
#source = destination[selection]
transformation = AffineTransform(scale=(0.1, 0.1), rotation=np.pi, shear=None, translation=(-100,350))
source = transformation(crop_coordinates(destination, source_bounds_in_destination))
source_bounds = transformation(source_bounds_in_destination)
plt.figure()
plt.scatter(destination[:,0],destination[:,1])
plt.scatter(source[:,0],source[:,1])

#
# t0 = time.time()
# from trace_analysis.mapping.geometricHashing import pointHash, findMatch
# ht = pointHash(destination, bases='all', magnificationRange=[0,10], rotationRange=[-np.pi,np.pi])
# matched_bases = findMatch(source, ht, bases='all', magnificationRange=[0,1], rotationRange=[-np.pi,np.pi])
# source_coordinate_tuple = source[matched_bases['testBasis']]
# destination_coordinate_tuple = destination[matched_bases['hashTableBasis']]
# t1 = time.time()
# plt.close('all')


# centers = np.array([(destination[pair[0]]+destination[pair[1]])/2 for pair in pairs])
#
# centers_KDTree = KDTree(centers)
#
# # Points within the circle containing the two original points in the pair
# indices_in_between_pairs = centers_KDTree.query_ball_tree(destination_KDTree, r=distance/2*0.8)


def geometric_hash(point_set, maximum_distance=100, point_set_size=4):
    # TODO: Add minimum_distance and implement
    # TODO: Make invariant to mirroring

    point_set_KDTree = KDTree(point_set)

    pairs = list(point_set_KDTree.query_pairs(maximum_distance))

    point_tuples = []
    hash_table = []
    for pair in pairs:
        pair_coordinates = point_set[list(pair)]
        center = (pair_coordinates[0] + pair_coordinates[1]) / 2
        distance = np.linalg.norm(pair_coordinates[0]-pair_coordinates[1])
        internal_points = point_set_KDTree.query_ball_point(center, distance/2*0.99)
        if len(internal_points) >= (point_set_size-2):
            random.shuffle(internal_points)
            internal_points = internal_points[0:(point_set_size-2)]

            internal_coordinates = point_set[internal_points]
            # Sort internal points based on x coordinate
            # internal_point_order = np.argsort(internal_coordinates[:, 0])
            # internal_points = [internal_points[i] for i in internal_point_order]
            # internal_coordinates = internal_coordinates[internal_point_order]

            end_points = np.array([[0, 0], [1, 1]])
            hash_coordinates = mapToPoint(internal_coordinates, pair_coordinates, end_points)
            # Break similarity of pair
            if np.sum(hash_coordinates[:,0]) > ((point_set_size-2)/2):
                pair = pair[::-1]
                pair_coordinates = pair_coordinates[::-1]
                hash_coordinates = mapToPoint(internal_coordinates, pair_coordinates, end_points)

            # Break similarity of internal points
            internal_point_order = np.argsort(hash_coordinates[:, 0])
            internal_points = [internal_points[i] for i in internal_point_order]
            internal_coordinates = internal_coordinates[internal_point_order]
            hash_coordinates = hash_coordinates[internal_point_order]

            hash_code = hash_coordinates.flatten()
            hash_table.append(hash_code)

            point_tuples.append(pair + tuple(internal_points))

            # plot_tuple(np.vstack([pair_coordinates, internal_coordinates])
            # plot_tuple(np.vstack([end_points, hash_coordinates]))

    hash_table_KDTree = KDTree(np.array(hash_table))

    return point_set_KDTree, point_tuples, hash_table_KDTree

def find_match_after_hashing(source_KDTree, source_tuples, source_hash_table_KDTree, source_bounds,
                             destination_KDTree, destination_tuples, destination_hash_table_KDTree):

    for source_tuple_index in np.arange(len(source_tuples)):
        distance, destination_tuple_index = destination_hash_table_KDTree.query(
            source_hash_table_KDTree.data[source_tuple_index])
        # We can also put a threshold on the distance here possibly
        # if distance < 0.01:
        print(distance, source_tuple_index)
        # source_coordinate_tuple = source[list(source_tuples[source_tuple_index])]
        # destination_coordinate_tuple = destination[list(destination_tuples[destination_tuple_index])]

        source_tuple = source_tuples[source_tuple_index]
        destination_tuple = destination_tuples[destination_tuple_index]
        match = tuple_match(source_KDTree, destination_KDTree, source_bounds, source_tuple, destination_tuple)
        if match:
            return match

def tuple_match(source_KDTree, destination_KDTree, source_bounds, source_tuple, destination_tuple):
    source_coordinate_tuple = source_KDTree.data[list(source_tuple)]
    destination_coordinate_tuple = destination_KDTree.data[list(destination_tuple)]

    source_transformed, transformation_matrix = mapToPoint(source_KDTree.data, source_coordinate_tuple[:2], destination_coordinate_tuple[:2], returnTransformationMatrix=True)
    # scatter_coordinates([source_transformed])

    found_transformation = AffineTransform(transformation_matrix)
    source_bounds_transformed = found_transformation(source_bounds)
    destination_cropped = crop_coordinates(destination, source_bounds_transformed)

    source_transformed_area = np.linalg.norm(source_bounds_transformed[0] - source_bounds_transformed[1])
    pDB = 1 / source_transformed_area

    alpha = 0
    test_radius = 5
    K=1
    K_threshold = 10e9
    for coordinate in source_transformed:
        points_within_radius = destination_KDTree.query_ball_point(coordinate, test_radius)
        pDF = alpha/source_transformed_area + (1-alpha)*len(points_within_radius)/len(destination_cropped)
        K = K * pDF/pDB
        if K > K_threshold:
            return found_transformation

    t1 = time.time()



t0 = time.time()

def find_match(source, destination, source_bounds):
    # destination_KDTree, destination_tuples, destination_hash_table_KDTree = geometric_hash(destination, 40, 4)
    # source_KDTree, source_tuples, source_hash_table_KDTree = geometric_hash(source, 4, 4)
    # 200 points 200,20
    # 10000 points 10,1
    return find_match_after_hashing(*geometric_hash(source, 4, 4), source_bounds, *geometric_hash(destination, 40, 4))

match = find_match(source, destination, source_bounds)
scatter_coordinates([source,destination,match(source)])


#
# destination_KDTree, destination_tuples, destination_hash_table_KDTree = geometric_hash(destination, 40, 4)
# source_KDTree, source_tuples, source_hash_table_KDTree = geometric_hash(source, 4, 4)





import pickle
from pathlib2 import Path
from trace_analysis.mapping.mapping import Mapping2
from trace_analysis.coordinate_transformations import translate, rotate, magnify, reflect, transform
from trace_analysis.plotting import plot_match

class SequencingDataMapping:
    def __init__(self, tile, files, dataPath):
        self.tile = tile
        self.files = files # List of coordinate sets
        self.dataPath = Path(dataPath)

        self.initial_image_transformation = {'reflection': 0, 'rotation': 0, 'magnification': 1} # This reflects with respect to axis 0.

        # self.mode = mode
        # if mode == 'translation': self.hashTableRange = [-10000, 10000]
        # else: self.hashTableRange = [-1,1]
        # self.nBins = nBins

        self.bases_hashTable = 'all'
        self.bases_findMatch = 20

        self.tuple_size = 4
        self.maximum_distance_tile = 40
        self.maximum_distance_file = 4

        self._hashTable = None
        self._matches = None

        # self.rotationRange = rotationRange
        # self.magnificationRange = magnificationRange


    @property
    def hashTable(self):
        if self._hashTable is None:
            if self.dataPath.joinpath(self.tile.name+'.ht').is_file():
                with self.dataPath.joinpath(self.tile.name+'.ht').open('rb') as f:
                    self._hashTable = pickle.load(f)
            else:
                # self._hashTable = pointHash(self.tile.coordinates, bases=self.bases_hashTable, mode=self.mode,
                #                             hashTableRange=self.hashTableRange, nBins=self.nBins,
                #                             rotationRange=self.rotationRange, magnificationRange=self.magnificationRange)
                self._hashTable = geometric_hash(self.tile.coordinates, self.maximum_distance_tile, self.tuple_size)
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

            file_hash = geometric_hash(coordinates, self.maximum_distance_file, self.tuple_size)

            #TODO: Make file.bounds use consistent
            match = find_match_after_hashing(*file_hash, file.bounds, *self.hashTable)
            # match = findMatch(coordinates, self.hashTable,
            #                                                           bases=self.bases_findMatch,
            #                                                           returnMatchedBases=True,
            #                                                           mode=self.mode,
            #                                                           nBins=self.nBins,
            #                                                           hashTableRange=self.hashTableRange,
            #                                                           rotationRange=self.rotationRange,
            #                                                           magnificationRange=self.magnificationRange)

            transformation_matrix = match.params

            transformation_matrix = transformation_matrix @ initial_transformation

            match = Mapping2(file.coordinates.copy(), self.tile.coordinates.copy(),
                                  method='geometric-hashing',
                                  transformation_type='linear')

            # match.name =
            match.transformation = transformation_matrix
            # match.best_image_basis = bestBasis
            # match.count = bestBasis['counts']
            # #match.image_coordinates_transformed = coordinates_transformed
            # match.meanError = np.mean(
            #     [np.min(np.linalg.norm(self.tile.coordinates - row, axis=1)) for row in coordinates_transformed])
            # match.setWidth = np.linalg.norm(np.max(coordinates_transformed, axis=0) - np.min(coordinates_transformed, axis=0))
            # match.percentMatch = bestBasis['counts'] / coordinates_transformed.shape[0]

            matches.append(match)

        return matches
    #
    # def histogram_matches(self, export = False):
    #     counts = [match.count for match in self.matches]
    #     fig, ax = plt.subplots(figsize = (6,3))
    #     ax.hist(counts, np.arange(0.5, np.max(counts) + 1.5, 1), facecolor = 'k', rwidth = 0.5)
    #
    #     ax.set_xlim([0,np.max(counts)+1])
    #
    #     ax.set_xlabel('Number of matches for best matching basis')
    #     ax.set_ylabel('Count')
    #     plt.tight_layout()
    #
    #     if export:
    #         fig.savefig(self.dataPath.joinpath('histogramNumberOfMatches.pdf'), bbox_inches='tight')
    #         fig.savefig(self.dataPath.joinpath('histogramNumberOfMatches.png'), bbox_inches='tight')

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

    def loop_through_matches(self, figure=plt.figure()):
        plt.ion()
        for match in self.matches:
            self.show_match(match, figure=figure)
            plt.show()
            plt.pause(0.001)
            input("Press enter to continue")

    def give_matches_to_files(self, match_threshold = 0):
        for file, match in zip(self.files, self.matches):
            if match.count >= match_threshold:
                file.sequence_match = match
                file.sequence_match.tile = self.tile





# scatter_coordinates([destination,source_transformed])
# print(t1-t0)



# destination_coordinate_tuples = [destination[list(t)] for t in destination_tuples]
# source_coordinate_tuples = [source[list(t)] for t in source_tuples]
# source_coordinate_tuple = source_coordinate_tuples[source_tuple_index]
# destination_coordinate_tuple = destination_coordinate_tuples[destination_tuple_index]


#
# scatter_coordinates([source_coordinate_tuple, destination_coordinate_tuple])
#
# def connect_pairs(pairs):
#     for pair in pairs:
#         plt.plot(pair[:,0], pair[:,1], color='r')
#
# def plot_tuple(tuple_coordinates):
#     pair_coordinates = tuple_coordinates[:2]
#     internal_coordinates = tuple_coordinates[2:]
#     plt.figure()
#     center = (pair_coordinates[0] + pair_coordinates[1]) / 2
#     distance = np.linalg.norm(pair_coordinates[0] - pair_coordinates[1])
#     connect_pairs([pair_coordinates])
#     scatter_coordinates([pair_coordinates, np.atleast_2d(center), internal_coordinates])
#     plt.gca().set_aspect('equal')
#     circle = plt.Circle(center, distance / 2, fill=False)
#     plt.gcf().gca().add_artist(circle)
#     axis_limits = np.array([center-distance/2*1.2, center+distance/2*1.2])
#     plt.xlim(axis_limits[:, 0])
#     plt.ylim(axis_limits[:, 1])
#     plt.xlabel('x')
#     plt.ylabel('y')







