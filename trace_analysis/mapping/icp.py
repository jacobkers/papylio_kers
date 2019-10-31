import numpy as np
from sklearn.neighbors import NearestNeighbors
from trace_analysis.coordinate_transformations import transform

def best_fit_transform(A, B):
    '''
    Calculates the least-squares best-fit transform that maps corresponding points A to B in m spatial dimensions
    Input:
      A: Nxm numpy array of corresponding points
      B: Nxm numpy array of corresponding points
    Returns:
      T: (m+1)x(m+1) homogeneous transformation matrix that maps A on to B
      R: mxm rotation matrix
      t: mx1 translation vector
    '''

    assert A.shape == B.shape

    # get number of dimensions
    m = A.shape[1]

    # translate points to their centroids
    centroid_A = np.mean(A, axis=0)
    centroid_B = np.mean(B, axis=0)
    AA = A - centroid_A
    BB = B - centroid_B

    # rotation matrix
    H = np.dot(AA.T, BB)
    U, S, Vt = np.linalg.svd(H)
    R = np.dot(Vt.T, U.T)

    # special reflection case
    if np.linalg.det(R) < 0:
       Vt[m-1,:] *= -1
       R = np.dot(Vt.T, U.T)

    # translation
    t = centroid_B.T - np.dot(R,centroid_A.T)

    # homogeneous transformation
    T = np.identity(m+1)
    T[:m, :m] = R
    T[:m, m] = t

    return T, R, t



def least_squares_fit(src, dst):
    T, res, rank, s = np.linalg.lstsq(src, dst, rcond=None)
    plt.figure()
    plt.scatter(src[:,0],src[:,1], color = 'b')
    plt.scatter(dst[:,0],dst[:,1], color = 'r')

   # src_transformed = ((T) @ (src.T)).T
    plt.scatter(src_transformed[:, 0], src_transformed[:, 1], color = 'g')



def nearest_neighbor(src, dst):
    '''
    Find the nearest (Euclidean) neighbor in dst for each point in src
    Input:
        src: Nxm array of points
        dst: Nxm array of points
    Output:
        distances: Euclidean distances of the nearest neighbor
        indices: dst indices of the nearest neighbor
    '''

    #assert src.shape == dst.shape

    neigh = NearestNeighbors(n_neighbors=1)
    neigh.fit(dst)
    distances, indices = neigh.kneighbors(src, return_distance=True)
    return distances.ravel(), indices.ravel()

def nearest_neighbor_pair(pointset1, pointset2):
    distances2, indices2 = nearest_neighbor(pointset1, pointset2)
    distances1, indices1 = nearest_neighbor(pointset2, pointset1)

    i1 = indices1[indices2[indices1] == np.arange(len(indices1))]
    i2 = np.where(indices2[indices1] == np.arange(len(indices1)))[0]
    # s2 = indices1[indices2] == np.arange(len(indices2))
    #
    # i2 = indices2[indices1[indices2] == np.arange(len(indices2))]

    return distances1[i2], i1, i2

import matplotlib.pyplot as plt
def scatter_coordinates(pointsets):
    for pointset in pointsets:
        plt.scatter(pointset[:,0], pointset[:,1])

def show_point_connections(pointset1,pointset2):
    for coordinate1, coordinate2 in zip(pointset1, pointset2):
        plt.plot([coordinate1[0],coordinate2[0]],[coordinate1[1],coordinate2[1]], color='r')

def icp(source, destination, max_iterations=20, tolerance=0.001):
    '''
    The Iterative Closest Point method: finds best-fit transform that maps points A on to points B
    Input:
        A: Nxm numpy array of source mD points
        B: Nxm numpy array of destination mD point
        init_pose: (m+1)x(m+1) homogeneous transformation
        max_iterations: exit algorithm after max_iterations
        tolerance: convergence criteria
    Output:
        T: final homogeneous transformation that maps A on to B
        distances: Euclidean distances (errors) of the nearest neighbor
        i: number of iterations to converge
    '''


    source = np.hstack([source, np.ones((len(source), 1))])
    destination = np.hstack([destination, np.ones((len(destination), 1))])

    print(source)

    source_dummy = source.copy()
    #transformation_final = np.identity(3)

    # Initial translation to overlap both point-sets
    initial_translation = np.identity(3)
    initial_translation[0:2,2] = (np.mean(destination, axis=0) - np.mean(source, axis=0))[0:2]
    # Possibly add initial rotation and reflection as well, see best_fit_transform?

    source_dummy = (initial_translation @ source_dummy.T).T
    #transformation_final = transformation_final @ initial_translation

    prev_error = 0

    for i in range(max_iterations):
        # Find the nearest neighbors between the current source and destination points
        #distances, indices = nearest_neighbor(source_dummy[:,:2], destination[:,:2])
        distances, source_indices, destination_indices = \
            nearest_neighbor_pair(source_dummy[:, :2], destination[:, :2])

        plt.figure()
        scatter_coordinates([source_dummy,destination])
        show_point_connections(source_dummy[source_indices],destination[destination_indices])

        # compute the transformation between the current source and nearest destination points
        # T,_,_ = best_fit_transform(src[:m,:].T, dst[:m,indices].T)

        T, res, rank, s = np.linalg.lstsq(source_dummy[source_indices], destination[destination_indices], rcond=None)

        transformation = T.T
        # Update the source dummy
        source_dummy = (transformation @ source_dummy.T).T
        #transformation_final = transformation @ transformation_final

        # Check error
        mean_error = np.mean(distances)
        print(mean_error)
        if np.abs(prev_error - mean_error) < tolerance:
            break
        prev_error = mean_error

    # Calculate final transformation
    T, res, rank, s = np.linalg.lstsq(source[source_indices], destination[destination_indices], rcond=None)
    transformation = T.T

    return transformation, distances, i


if __name__ == '__main__':
    Npoints = 40

    plt.close('all')
    np.random.seed(32)
    source = np.random.rand(Npoints, 2) * 1000
    destination = source.copy()
    #destination = destination * 1.20 - 100 +
    destination = transform(source, r=np.pi/180*10, m=1.2, t=[0, 0]) #+ np.random.uniform(-20, 20, (Npoints, 2))
    np.random.shuffle(destination)
    destination = destination[:30]

    transformation, distances, i = icp(source,destination)

