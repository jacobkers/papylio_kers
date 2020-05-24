import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as pth
from scipy.spatial import KDTree
import random
from trace_analysis.mapping.geometricHashing import mapToPoint
from trace_analysis.plotting import scatter_coordinates
from skimage.transform import AffineTransform
import time


def crop_coordinates(coordinates, vertices):
    return coordinates[pth.Path(vertices).contains_points(coordinates)]

    # bounds.sort(axis=0)
    # selection = (coordinates[:, 0] > bounds[0, 0]) & (coordinates[:, 0] < bounds[1, 0]) & \
    #             (coordinates[:, 1] > bounds[0, 1]) & (coordinates[:, 1] < bounds[1, 1])
    # return coordinates[selection]

#
# selection = (destination[:,0] > source_bounds[0,0]) & (destination[:,0] < source_bounds[1,0]) & \
#             (destination[:,1] > source_bounds[0,1]) & (destination[:,1] < source_bounds[1,1])
#source = destination[selection]


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

def find_match_after_hashing(source_KDTree, source_tuples, source_hash_table_KDTree, source_vertices,
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
        match = tuple_match(source_KDTree, destination_KDTree, source_vertices, source_tuple, destination_tuple)
        if match:
            return match

def tuple_match(source_KDTree, destination_KDTree, source_vertices, source_tuple, destination_tuple):
    source_coordinate_tuple = source_KDTree.data[list(source_tuple)]
    destination_coordinate_tuple = destination_KDTree.data[list(destination_tuple)]

    source_transformed, transformation_matrix = mapToPoint(source_KDTree.data, source_coordinate_tuple[:2], destination_coordinate_tuple[:2], returnTransformationMatrix=True)
    # scatter_coordinates([source_transformed])

    found_transformation = AffineTransform(transformation_matrix)
    source_vertices_transformed = found_transformation(source_vertices)
    destination_cropped = crop_coordinates(destination_KDTree.data, source_vertices_transformed)

    #source_transformed_area = np.linalg.norm(source_vertices_transformed[0] - source_vertices_transformed[1])
    source_transformed_area = np.abs(np.cross(source_vertices_transformed[1] - source_vertices_transformed[0],
                                              source_vertices_transformed[3] - source_vertices_transformed[0]))
    pDB = 1 / source_transformed_area

    alpha = 0
    # test_radius = 5
    test_radius = 50
    K=1
    K_threshold = 10e9
    for coordinate in source_transformed:
        points_within_radius = destination_KDTree.query_ball_point(coordinate, test_radius)
        pDF = alpha/source_transformed_area + (1-alpha)*len(points_within_radius)/len(destination_cropped)
        K = K * pDF/pDB
        if K > K_threshold:
            return found_transformation



def find_match(source, destination, source_vertices):
    # destination_KDTree, destination_tuples, destination_hash_table_KDTree = geometric_hash(destination, 40, 4)
    # source_KDTree, source_tuples, source_hash_table_KDTree = geometric_hash(source, 4, 4)
    # 200 points 200,20
    # 10000 points 10,1
    return find_match_after_hashing(*geometric_hash(source, 4, 4), source_vertices, *geometric_hash(destination, 40, 4))



if __name__ == '__main__':
    # # Make random source and destination dataset
    # np.random.seed(42)
    # destination = np.random.rand(1000,2)*1000
    # source_vertices_in_destination = np.array([[300, 300], [450, 300], [450, 600], [300,600]])
    #
    #
    # transformation = AffineTransform(scale=(0.1, 0.1), rotation=np.pi, shear=None, translation=(-100,350))
    # source = transformation(crop_coordinates(destination, source_vertices_in_destination))
    # source_vertices = transformation(source_vertices_in_destination)
    # plt.figure()
    # plt.scatter(destination[:,0],destination[:,1])
    # plt.scatter(source[:,0],source[:,1])
    # scatter_coordinates([source,destination,crop_coordinates(destination, source_vertices_in_destination)])
    #
    # destination_KDTree, destination_tuples, destination_hash_table_KDTree = geometric_hash(destination, 40, 4)
    # source_KDTree, source_tuples, source_hash_table_KDTree = geometric_hash(source, 4, 4)
    #
    # match = find_match(source, destination, source_vertices)
    # scatter_coordinates([source, destination, match(source), source_vertices,source_vertices_in_destination])



    file_coordinates = np.loadtxt(
        r'C:\Users\Ivo Severins\Desktop\seqdemo\20190924 - Single-molecule setup (TIR-I)\16L\spool_6.pks')[:, 1:3]
    tile_coordinates = np.loadtxt(r'C:\Users\Ivo Severins\Desktop\seqdemo\20190926 - Sequencer (MiSeq)\2102.loc')

    source = file_coordinates
    destination = tile_coordinates
    source_vertices = np.array([[1024,0],[2048,0],[2048,2048],[1024,2048]])

    # source = file_coordinates[[1,5,9,7]]
    # destination = tile_coordinates[[73,82,85,84]]

    source[:,0] = -source[:,0]
    source_vertices[:,0] = -source_vertices[:,0]

    # destination_KDTree, destination_tuples, destination_hash_table_KDTree = geometric_hash(destination, 3000, 4)
    # source_KDTree, source_tuples, source_hash_table_KDTree = geometric_hash(source, 100, 4)

    source_hash_data = geometric_hash(source, 1000, 4)
    destination_hash_data = geometric_hash(destination, 3000, 4)

    source_KDTree, source_tuples, source_hash_table_KDTree = source_hash_data
    destination_KDTree, destination_tuples, destination_hash_table_KDTree = destination_hash_data

    match = find_match_after_hashing(*source_hash_data, source_vertices, *destination_hash_data)

    #match = find_match(source, destination, source_vertices)
    #scatter_coordinates([source, destination, match(source), source_vertices, match(source_vertices)])



    def scatter_with_number(coordinates):
        # plot the chart
        plt.scatter(coordinates[:,0], coordinates[:,1])

        # zip joins x and y coordinates in pairs
        for i, coordinate in enumerate(coordinates):
            label = str(i)

            # this method is called for each point
            plt.annotate(label,  # this is the text
                         (coordinate[0], coordinate[1]),  # this is the point to label
                         textcoords="offset points",  # how to position the text
                         xytext=(0, 10),  # distance from text to points (x,y)
                         ha='center')

    scatter_with_number(source)
    scatter_with_number(destination)

    for i, t in enumerate(source_tuples):
        if ((1 in t) & (5 in t) & (9 in t) & (7 in t)):
            print(i)


    #
    #
    # def show_match(self, match, figure = None, view='destination'):
    #     if not figure: figure = plt.gcf()
    #     figure.clf()
    #     ax = figure.gca()
    #
    #     #ax.scatter(ps2[:,0],ps2[:,1],c='g',marker = '+')
    #
    #     ax.scatter(match.destination[:,0],match.destination[:,1], marker = '.', facecolors = 'k', edgecolors='k')
    #     ax.scatter(match.transform_source_to_destination[:,0],match.transform_source_to_destination[:,1],c='r',marker = 'x')
    #
    #     destination_basis_index = match.best_image_basis['hashTableBasis']
    #     source_basis_index = match.best_image_basis['testBasis']
    #     ax.scatter(match.destination[destination_basis_index, 0], match.destination[destination_basis_index, 1], marker='.', facecolors='g', edgecolors='g')
    #     ax.scatter(match.transform_source_to_destination[source_basis_index, 0], match.transform_source_to_destination[source_basis_index, 1], c='g',
    #                marker='x')
    #
    #     ax.set_aspect('equal')
    #     ax.set_title('Tile:' + self.tile.name +', File: ' + str(self.files[self.matches.index(match)].relativeFilePath))
    #
    #     if view == 'source':
    #         maxs = np.max(match.transform_source_to_destination, axis=0)
    #         mins = np.min(match.transform_source_to_destination, axis=0)
    #         ax.set_xlim([mins[0], maxs[0]])
    #         ax.set_ylim([mins[1], maxs[1]])
    #     elif view == 'destination':
    #         maxs = np.max(match.destination, axis=0)
    #         mins = np.min(match.destination, axis=0)
    #         ax.set_xlim([mins[0], maxs[0]])
    #         ax.set_ylim([mins[1], maxs[1]])
    #         # ax.set_xlim([0, 31000])
    #         # ax.set_ylim([0, 31000])
    #
    #     name = str(self.files[self.matches.index(match)].relativeFilePath)
    #     print(name)
    #     n = name.replace('\\', '_')
    #
    #     figure.savefig(self.dataPath.joinpath(n + '_raw.pdf'), bbox_inches='tight')
    #     figure.savefig(self.dataPath.joinpath(n + '_raw.png'), bbox_inches='tight', dpi=1000)
    #
    # def loop_through_matches(self, figure=plt.figure()):
    #     plt.ion()
    #     for match in self.matches:
    #         self.show_match(match, figure=figure)
    #         plt.show()
    #         plt.pause(0.001)
    #         input("Press enter to continue")





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







