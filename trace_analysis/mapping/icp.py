import numpy as np
from sklearn.neighbors import NearestNeighbors
from trace_analysis.image_adapt.polywarp import polywarp, polywarp_apply #required for nonlinear
import cv2 #required for nonlinear 
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

    return distances1[i2], i1, i2

import matplotlib.pyplot as plt
def scatter_coordinates(pointsets, marker=['+','x'], c=['g','r']):
    if type(pointsets)==np.ndarray:
        plt.scatter(pointsets[:,0], pointsets[:,1], marker=marker[0], c=c[0])
    if type(pointsets)==list:
        if len(pointsets)==2:
            for ii in range(len(pointsets)): plt.scatter(pointsets[ii][:,0], pointsets[ii][:,1], marker=marker[ii], c=c[ii])
        else: 
            for ii in range(len(pointsets)): plt.scatter(pointsets[ii][:,0], pointsets[ii][:,1])

def show_point_connections(pointset1,pointset2):
    for coordinate1, coordinate2 in zip(pointset1, pointset2):
        plt.plot([coordinate1[0],coordinate2[0]],[coordinate1[1],coordinate2[1]], color='gray')

def icp(source, destination, max_iterations=20, tolerance=0.001, dest2source_translation = None, transformation_type = 'linear'):
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

    source = np.hstack([source, np.ones((len(source), 1))]) # source=donor=left side of the image
    destination = np.hstack([destination, np.ones((len(destination), 1))]) #destination=acceptor=right side of the image

    destination_dummy = destination.copy() # destination_dummy=destination (=acceptor) transformed to donor location, left side image
    #transformation_final = np.identity(3)
    if dest2source_translation is None:
        # Initial translation to overlap both point-sets
        dest2source_translation = np.identity(3)
        dest2source_translation[0:2,2] = (np.mean(source, axis=0) - np.mean(destination, axis=0))[0:2] # need to be set back to mapping2
        # Possibly add initial rotation and reflection as well, see best_fit_transform?  
    '''destination_dummy is the destination, moved to source location'''
    destination_dummy = (dest2source_translation @ destination_dummy.T).T
    #transformation_final = transformation_final @ initial_translation

    prev_error = 0

    for i in range(max_iterations):
        print(i)
        # Find the nearest neighbors between the current source and destination_dummy points (which are in the same left part of the image now)
        # distances, indices = nearest_neighbor(source_dummy[:,:2], destination[:,:2])
        distances, source_indices, destination_indices = \
            nearest_neighbor_pair(source[:, :2], destination_dummy[:, :2])

#        plt.figure() # don't plot for every iteration --> move to after the lop
#        scatter_coordinates([source_dummy,destination])
#        show_point_connections(source_dummy[source_indices],destination[destination_indices])
        if transformation_type=='nonlinear':
            # Might be useful as well for linear. with cutoff you remove the entries with outlier distances
            cutoff = np.median(distances)+np.std(distances) #changed
            source_indices = source_indices[distances<cutoff]
            destination_indices = destination_indices[distances<cutoff]

   #         kx_s2d, ky_s2d = polywarp(destination_dummy[destination_indices,:], source[source_indices,:]) 
            kx_d2s, ky_d2s = polywarp(source[source_indices,:],destination_dummy[destination_indices,:]) 
            # these are the values needed to transform source into destination_dummy
            
            destination_dummy = polywarp_apply(kx_d2s, ky_d2s, destination_dummy)


        elif transformation_type=='linear':
			# compute the transformation between the current source and nearest destination points
			#T,_,_ = best_fit_transform(src[:m,:].T, dst[:m,indices].T)
            T, res, rank, s = np.linalg.lstsq(source[source_indices], destination_dummy[destination_indices], rcond=None)
            transformation_s2d = T.T
			# Update the source dummy
            destination_dummy = (transformation_s2d @ destination_dummy.T).T

        plt.figure(i) # don't plot for every iteration --> move to after the lop
        scatter_coordinates([source,destination_dummy])
        show_point_connections(source[source_indices],destination_dummy[destination_indices])
        # Check error
        mean_error = np.mean(distances)
        print(i, mean_error, prev_error, len(distances))
        if np.abs(prev_error - mean_error) < tolerance:
            break
        prev_error = mean_error
    # continue the loop until the improvement in match is limited. The outcome are a set of matching locations, afterwards do one more transform to find final transform matrix
        
    print(i)
    plt.figure() # don't plot for every iteration --> move to after the lop
    scatter_coordinates([source,destination_dummy])
    show_point_connections(source[source_indices],destination_dummy[destination_indices])
 
    # Calculate final transformation, need to be redone, since above you retrieve kx,ky per iteration, now you want the overall one
    if transformation_type == 'nonlinear': ## zit hier de initiele translatie nog in??
        kx_d2s, ky_d2s = polywarp(source[source_indices,:],destination[destination_indices,:]) ## or destination_dummy, really depends on how you do apply transform
        kx_s2d, ky_s2d = polywarp(destination[destination_indices,:],source[source_indices,:])
        transformation = (kx_s2d, ky_s2d) # should be renamed to transformL2R?
        transformation_inverse = (kx_d2s,ky_d2s) #  should be renamed to transformR2L?
        
    elif transformation_type=='linear': # replace with transform 
        T, res, rank, s = np.linalg.lstsq(source[source_indices], destination[destination_indices], rcond=None)
        transformation = T.T
        transformation_inverse= np.linalg.inv(transformation)
    return transformation, distances, i, transformation_inverse, dest2source_translation

def icp_apply_transform  (coordinates, inverse, self_transformation,self_transformation_inverse, self_transformation_type,dest2source_translation):
 
    # first move destination (=acceptor) to source
    coords_dummy= coordinates.copy()
    
    if self_transformation_type == 'linear': #$$$$$$$$$$$$$$$$$
        # Maybe we should rename transform to linear_transform [IS 05-03-2020]
        if inverse is False: return transform([source,coords_dummy[:,:2]], self_transformation)
        elif inverse is True: return transform([source,coords_dummy[:,:2]], self_transformation_inverse)

    elif self_transformation_type == 'nonlinear':
            if inverse is False: return polywarp_apply(self_transformation[0],self_transformation[1],coords_dummy)
            elif inverse is True: return polywarp_apply(self_transformation_inverse[0],self_transformation_inverse[1],coords_dummy)
    
    return coords_dummy[:,:1]
            
if __name__ == '__main__':
    Npoints = 40

    np.random.seed(32)
    source = np.random.rand(Npoints, 2) * 1000
    destination = source.copy()
    #destination = destination * 1.20 - 100 +
    destination = transform(source, r=np.pi/180*10, m=1.2, t=[0, 0]) #+ np.random.uniform(-20, 20, (Npoints, 2))
    np.random.shuffle(destination)
    destination = destination[:30]

    transformation, distances, i = icp(source,destination)

