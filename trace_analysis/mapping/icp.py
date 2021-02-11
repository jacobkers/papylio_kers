import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors

from trace_analysis.plotting import scatter_coordinates, show_point_connections
from trace_analysis.image_adapt.polywarp import polywarp, polywarp_apply #required for nonlinear
import cv2 #required for nonlinear
from skimage.transform import AffineTransform, PolynomialTransform
from trace_analysis.mapping.polywarp import PolywarpTransform
# from trace_analysis.coordinate_transformations import transform, translate
import matplotlib.pyplot as plt

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

# Do we use this somewhere? I guess it can be removed, as it doesn't return anything. [IS: 29-11-2020]
def least_squares_fit(src, dst):
    T, res, rank, s = np.linalg.lstsq(src, dst, rcond=None)
    plt.figure()
    plt.scatter(src[:,0],src[:,1], color = 'b')
    plt.scatter(dst[:,0],dst[:,1], color = 'r')

   # src_transformed = ((T) @ (src.T)).T
   # plt.scatter(src_transformed[:, 0], src_transformed[:, 1], color = 'g')

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


def icp(source, destination, max_iterations=20, tolerance=0.001, cutoff=None, initial_transformation=None,
        transform=AffineTransform, transform_final=None, **kwargs):
    '''Iterative closest point algorithm for point-set registration
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

    if transform_final is None:
        transform_final=transform

    if cutoff == 'auto':
        auto_cutoff = True
    else:
        auto_cutoff = False

    plot = icp_plot()
    plot.append_data(source, destination, title='Start')

    source_moving_to_destination = source.copy() # destination_moved2source=destination transformed to source location, left side image

    #transformation_final = np.identity(3)
    if initial_transformation is None:
        # Initial translation to overlap both point-sets
        initial_transformation = AffineTransform(translation=(np.mean(destination, axis=0) - np.mean(source, axis=0))) # need to be set back to mapping2
        # Possibly add initial rotation and reflection as well, see best_fit_transform?  
    '''destination_moved2source is the destination, moved to source location'''
    source_moving_to_destination = initial_transformation(source_moving_to_destination)
    #transformation_final = transformation_final @ initial_translation

    plot.append_data(source_moving_to_destination, destination, title='Initial transformation')

    previous_error = 0


    for i in range(max_iterations):
        # print(i)
        # Find the nearest neighbors between the current source and destination_moved2source points (which are in the same left part of the image now)
        # use = nearest_neighbor()
        distances, source_indices, destination_indices = \
            nearest_neighbor_pair(source_moving_to_destination, destination)

        if auto_cutoff:
            cutoff = np.median(distances) + np.std(distances)

        if type(cutoff) in (float, int):
            source_indices = source_indices[distances < cutoff]
            destination_indices = destination_indices[distances < cutoff]

        transformation_step = transform()
        transformation_step.estimate(source_moving_to_destination[source_indices], destination[destination_indices], **kwargs)
        source_moving_to_destination = transformation_step(source_moving_to_destination)

        # if transformation_type=='nonlinear':
        #     kx, ky = polywarp(destination[destination_indices,0:2], source_moving_to_destination[source_indices,0:2])
        #     # these are the values needed to transform source into destination_moved2source
        #     source_moving_to_destination = polywarp_apply(kx, ky, source_moving_to_destination)
        #
        # elif transformation_type=='linear':
        #     # compute the transformation between the current source and nearest destination points
        #     #T,_,_ = best_fit_transform(src[:m,:].T, dst[:m,indices].T)
        #     T, res, rank, s = np.linalg.lstsq(source_moving_to_destination[source_indices], destination[destination_indices], rcond=None)
        #     transformation = T.T
        #     source_moving_to_destination = (transformation @ source_moving_to_destination.T).T

        plot.append_data(source_moving_to_destination, destination, source_indices, destination_indices,
                         title=f'Step {i}')
        # Check error
        mean_error = np.mean(distances) # These are not the distances after cutoff
        mean_squared_error = np.sqrt(np.mean(distances**2))
        print(f'Iteration: {i} \t Mean_squared_error: {mean_squared_error} \t Number of pairs: {len(distances)}')
        if np.abs(previous_error - mean_squared_error) < tolerance:
            break
        # previous_error = mean_error
        previous_error = mean_squared_error
    # continue the loop until the improvement in match is limited. The outcome are a set of matching locations, afterwards do one more transform to find final transform matrix

    transformation_final = transform_final()
    transformation_final.estimate(source[source_indices], destination[destination_indices], **kwargs)

    transformation_final_inverse = transform_final()
    transformation_final_inverse.estimate(destination[destination_indices], source[source_indices], **kwargs)

    distances, _, _ = \
        nearest_neighbor_pair(transformation_final(source[source_indices]), destination[destination_indices])
    mean_squared_error = np.sqrt(np.mean(distances ** 2))
    print(f'Final \t\t Mean squared error: {mean_squared_error} \t Number of pairs: {len(distances)}')

    plot.append_data(transformation_final(source), destination, source_indices, destination_indices,
                     title=f'Final')
    plot.plot()

    # # Calculate final transformation, need to be redone, since above you retrieve kx,ky per iteration, now you want the overall one
    # if transformation_type == 'nonlinear': ## zit hier de initiele translatie nog in??
    #     kx_inv, ky_inv = polywarp(source[source_indices,:],destination[destination_indices,:])
    #     kx, ky = polywarp(destination[destination_indices,:],source[source_indices,:])
    #     transformation = (kx, ky) # should be renamed to transformL2R?
    #     transformation_inverse = (kx_inv,ky_inv) #  should be renamed to transformR2L?
    #
    # elif transformation_type=='linear': # replace with transform
    #     T, res, rank, s = np.linalg.lstsq(source[source_indices], destination[destination_indices], rcond=None)
    #     transformation = T.T
    #     transformation_inverse = np.linalg.inv(transformation)

    return transformation_final, transformation_final_inverse, distances, i


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

class icp_plot:
    def __init__(self):
        self.data = []
        self.xlim = [0,512]
        self.ylim = [0,512]

    def append_data(self, source, destination, source_indices=None, destination_indices=None,
                    title=None):
        self.data.append((source, destination, source_indices, destination_indices, title))

    def set_step(self, step_index):
        step_index = int(step_index)
        self.axis.cla()
        self.axis.set_xlim(self.xlim)
        self.axis.set_ylim(self.ylim)
        plot_icp_step(*self.data[step_index], self.axis)
        self.figure.canvas.draw_idle()
        # plt.title('#donor=' + str(len(source)) + '#acceptor=' + str(len(destination)) + '#overlap=' + str(
        #     len(source_indices)))

    def plot(self):
        self.figure, self.axis = plt.subplots()
        plt.subplots_adjust(bottom=0.15)
        self.axis.set_aspect(1)
        axcolor = 'lightgoldenrodyellow'
        slider_axis = plt.axes([0.1, 0.05, 0.8, 0.025], facecolor=axcolor)

        self.slider = Slider(slider_axis, 'Step', 0, len(self.data)-1,
                             valinit=len(self.data), valstep=1)

        self.slider.on_changed(self.set_step)

        self.set_step(len(self.data)-1)

def plot_icp_step(source, destination, source_indices=None, destination_indices=None,
                  title=None, axis=None):
    if not axis:
        axis = plt.gca()
    if (source_indices is not None) and (destination_indices is not None):
        show_point_connections(source[source_indices], destination[destination_indices], axis)
    scatter_coordinates([destination, source], axis)
    axis.set_title(title)


if __name__ == '__main__':
    # Npoints = 40
    #
    # np.random.seed(32)
    # source = np.random.rand(Npoints, 2) * 1000
    # destination = source.copy()
    # #destination= destination * 1.20 - 100 +
    # destination = transform(source, r=np.pi/180*3, m=1.1, t=[0, 0]) #+ np.random.uniform(-20, 20, (Npoints, 2))
    # np.random.shuffle(destination)
    # destination = destination[:30]

    # plt.figure(i)  # don't plot for every iteration --> move to after the lop
    # scatter_coordinates([source_moving_to_destination, destination])
    #
    # transformation, distances, i, transformation_inverse = icp(source,destination,transformation_type='linear')#, initial_translation=translate([0,0]))
    # transformation, distances, i, transformation_inverse = icp(source, destination,
    #                                                            transformation_type='nonlinear')  # , initial_translation=translate([0,0]))

    from trace_analysis.plugins.sequencing.point_set_simulation import simulate_mapping_test_point_set


    number_of_source_points = 40
    transformation = AffineTransform(translation=[256,0], rotation=5/360*2*np.pi, scale=[0.98, 0.98])
    source_bounds = ([0, 0], [256, 512])
    source_crop_bounds = None
    fraction_missing_source = 0.1
    fraction_missing_destination = 0.1
    maximum_error_source = 4
    maximum_error_destination = 4
    shuffle = True

    source, destination = simulate_mapping_test_point_set(number_of_source_points, transformation,
                                                          source_bounds, source_crop_bounds,
                                                          fraction_missing_source, fraction_missing_destination,
                                                          maximum_error_source, maximum_error_destination, shuffle)

    max_iterations = 20
    tolerance = 0.0000001
    cutoff = None
    initial_transformation = None
    transform = AffineTransform
    transform_final = AffineTransform

    transformation, transformation_inverse, distances, i = icp(source, destination, max_iterations, tolerance,
                                                               cutoff, initial_transformation,
                                                               transform, transform_final)

