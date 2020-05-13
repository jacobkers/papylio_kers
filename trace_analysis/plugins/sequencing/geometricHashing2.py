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
source_bounds = np.array([[300, 300], [450, 600]])



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
source = transformation(crop_coordinates(destination, source_bounds))
source_bounds_transformed = transformation(source_bounds)
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

t0 = time.time()
destination_KDTree, destination_tuples, destination_hash_table_KDTree = geometric_hash(destination, 40, 4)

source_KDTree, source_tuples, source_hash_table_KDTree = geometric_hash(source, 4, 4)

# 200 points 200,20
# 10000 points 10,1

plt.figure()
scatter_coordinates([destination])

source_tuple_index = 0
for source_tuple_index in np.arange(len(source_tuples)):
    distance, destination_tuple_index = destination_hash_table_KDTree.query(source_hash_table_KDTree.data[source_tuple_index])
    # We can also put a threshold on the distance here possibly
    # if distance < 0.01:
    print(distance, source_tuple_index)
    source_coordinate_tuple = source[list(source_tuples[source_tuple_index])]
    destination_coordinate_tuple = destination[list(destination_tuples[destination_tuple_index])]

    source_transformed, transformation_matrix = mapToPoint(source, source_coordinate_tuple[:2], destination_coordinate_tuple[:2], returnTransformationMatrix=True)
    scatter_coordinates([source_transformed])

    found_transformation = AffineTransform(transformation_matrix)
    source_bounds_in_destination = found_transformation(source_bounds_transformed)
    destination_cropped = crop_coordinates(destination, source_bounds_in_destination)

    source_in_destination_area = np.linalg.norm(source_bounds_in_destination[0] - source_bounds_in_destination[1])
    pDB = 1 / source_in_destination_area

    alpha = 0
    test_radius = 5
    K=1
    K_threshold = 10e9
    is_match = False
    for coordinate in source_transformed:
        points_within_radius = destination_KDTree.query_ball_point(coordinate, test_radius)
        pDF = alpha/source_in_destination_area + (1-alpha)*len(points_within_radius)/len(destination_cropped)
        K = K * pDF/pDB
        if K > K_threshold:
            is_match = True
            break
    if is_match: break

    t1 = time.time()








scatter_coordinates([destination,source_transformed])
print(t1-t0)



#destination_coordinate_tuples = [destination[list(t)] for t in destination_tuples]
#source_coordinate_tuples = [source[list(t)] for t in source_tuples]
#source_coordinate_tuple = source_coordinate_tuples[source_tuple_index]
#destination_coordinate_tuple = destination_coordinate_tuples[destination_tuple_index]



scatter_coordinates([source_coordinate_tuple, destination_coordinate_tuple])

def connect_pairs(pairs):
    for pair in pairs:
        plt.plot(pair[:,0], pair[:,1], color='r')

def plot_tuple(tuple_coordinates):
    pair_coordinates = tuple_coordinates[:2]
    internal_coordinates = tuple_coordinates[2:]
    plt.figure()
    center = (pair_coordinates[0] + pair_coordinates[1]) / 2
    distance = np.linalg.norm(pair_coordinates[0] - pair_coordinates[1])
    connect_pairs([pair_coordinates])
    scatter_coordinates([pair_coordinates, np.atleast_2d(center), internal_coordinates])
    plt.gca().set_aspect('equal')
    circle = plt.Circle(center, distance / 2, fill=False)
    plt.gcf().gca().add_artist(circle)
    axis_limits = np.array([center-distance/2*1.2, center+distance/2*1.2])
    plt.xlim(axis_limits[:, 0])
    plt.ylim(axis_limits[:, 1])
    plt.xlabel('x')
    plt.ylabel('y')