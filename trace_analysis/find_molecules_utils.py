'''
Functions called in File().find_molecules()




'''
import numpy as np
from trace_analysis import peak_finding
from trace_analysis import coordinate_transformations as ctrans
import skimage as ski

import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

'''
Main function called directly within 'file.py'--> 'find_molecules()'
'''
def find_unique_molecules(file,
                          sliding_window=True,
                          projection_image_type='average',
                          method='by_channel',
                          channels=['donor'],
                          window_size=20,
                          uncertainty_pixels=2,
                          configs_peak_finding=None):
    '''
    this function will loop through frames of the movie:
    1. creates an average image or maximum projection (for every N frames)
    2. finds the peaks in the new image (wraps around built-in method)
    3. keeps only the unique positions/pixels of peak locations ()
         3B. allows to select only 'FRET pairs': intensity peak is seen in both donor/acceptor
    4. Finally, correct for an uncertainty in exact location due to photon shot noise etc. to remove 'functional duplicates'

    '''

    # --- fatch values from config file and objects ---
    nmbr_pixels = file.movie.width
    number_of_frames = file.number_of_frames


    # --- initialize output image ---
    find_coords_img = np.zeros((nmbr_pixels, nmbr_pixels))


    # --- make the windows
    # (if no sliding windows, just a single window is made to make it compatible with next bit of code) ----
    if sliding_window:
        window_start_frames = [i * window_size for i in range(number_of_frames // window_size)]
    else:
        window_start_frames = [0]

    # coordinates = set()
    if method == 'by_channel':
        coordinate_sets = [set() for channel in channels]
    elif method == 'overlay_channels':
        if len(channels) < 2:
            raise ValueError('No channels to overlay')
        coordinate_sets = [set()]

    #coordinates_sets = dict([(channel, set()) for channel in channels])
    # coordinate_sets = [set() for channel in channels]

    # --- Loop over all frames and find unique set of molecules ----
    for window_start_frame in window_start_frames:

        # --- allowed to apply sliding window to either the max projection OR the averages ----
        image = file.movie.make_projection_image(type=projection_image_type, start_frame=window_start_frame, number_of_frames=window_size)

        # --- we output the "sum of these images" ----
        find_coords_img += image

        if method == 'by_channel':
            coordinates_per_channel = dict([(channel, set()) for channel in channels])
            channel_images = [file.movie.get_channel(image=image, channel='donor') for channel in channels]

        if method == 'overlay_channels':
            # TODO: make this usable for any number of channels
            donor_image = file.movie.get_channel(image=image, channel='d')
            acceptor_image = file.movie.get_channel(image=image, channel='a')

            image_transformation = ctrans.translate([-file.movie.width / 2, 0]) @ file.mapping.transformation
            acceptor_image_transformed = ski.transform.warp(acceptor_image, image_transformation, preserve_range=True)
            # MD: problem: this is a linear transform, while yo u might have found a nonlinear transform; is nonlinear transform of image available?
            channel_images = [(donor_image + acceptor_image_transformed) / 2]
            channels = ['d']

        for i, channel_image in enumerate(channel_images):
            channel_coordinates = peak_finding.find_peaks(image=channel_image, **configs_peak_finding)#.astype(int)))
            channel_coordinates = set_of_tuples_from_array(channel_coordinates)

            coordinate_optimization_functions = \
                {'coordinates_within_margin': coordinates_within_margin,
                 'coordinates_after_gaussian_fit': coordinates_after_gaussian_fit,
                 'coordinates_without_intensity_at_radius': coordinates_without_intensity_at_radius}
            for f, kwargs in configuration['coordinate_optimization'].items():
                coordinates = coordinate_optimization_functions[f](coordinates, image, **kwargs)

            coordinate_sets[i].update(channel_coordinates)

    # --- correct for photon shot noise / stage drift ---
    # Not sure whether to put this in front of combine_coordinate_sets/detect_FRET_pairs or behind [IS: 12-08-2020]
    # I think before, as you would do it either for each window, or for the combined windows.
    # Transforming the coordinate sets for each window will be time consuming and changes the distance_threshold.
    # And you would like to combine the channel sets on the merged coordinates.
    for i in range(len(coordinate_sets)):
        if sliding_window:
            coordinate_sets[i] = merge_nearby_coordinates(coordinate_sets[i], distance_threshold=uncertainty_pixels)

        # Map coordinates to main channel in movie
        # TODO: make this usable for any number of channels
        coordinate_sets[i] = transform(coordinate_sets[i], translation=file.movie.channel_boundaries('a')[:,0])
        if channel=='a':
            coordinate_sets[i] = self.mapping.transform_coordinates(coordinate_sets[i], direction='destination2source')

    #TODO: make this work properly using a the mapping transformation
    #TODO: make this usable for any number of channels
    if len(coordinate_sets) == 1:
        coordinates = coordinate_sets[0]
    if len(coordinate_sets) > 1:
        raise NotImplementedError('Assessing found coordinates in multiple channels does not work properly yet')
        coordinates = combine_coordinate_sets(coordinate_sets, method='and') # the old detect_FRET_pairs

    donor_coordinates = coordinates
    acceptor_coordinates = self.mapping.transform_coordinates(coordinates, direction='source2destination')
    coordinates = np.hstack([donor_coordinates, acceptor_coordinates]).reshape((-1, 2))

    # --- turn into array ---
    coordinates = array_from_set_of_tuples(coordinates)

    return coordinates, find_coords_img


def merge_nearby_coordinates(coordinates, distance_threshold=2, plot=False):
    """Merge nearby coordinates to a single coordinate

    Coordinates are stored in a KD-tree.
    Each pair of points with a distance smaller than the distance threshold is obtained
    Pairs are chained to obtain groups of points
    For each group find the center coordinate and use that as a new coordinate
    (do this only if each member of the group is within the distance threshold from the center coordinate).
    Add individual points to the new coordinate list, i.e. points that do not have other points within the distance threshold.

    Parameters
    ----------
    coordinates : numpy.ndarray of ints or floats OR set of tuples
        Array with each row a set of coordinates
    distance_threshold : int or float
        Points closer than this distance are considered belonging to the same molecule.
    plot : bool
        If True shows a scatter plot of the coordinates and the new coordinates on top. (Only for 2D coordinates)

    Returns
    -------
    new_coordinates : numpy.ndarray of floats
        Coordinate array after merging nearby coordinates

    """

    # Convert to numpy array in case the coordinates are given as a set of tuples
    coordinates = array_from_set_of_tuples(coordinates)

    # Put coordinates in KD-tree for fast nearest-neighbour finding
    coordinates_KDTree = cKDTree(coordinates)

    # Determine pairs of points closer than the distance_threshold
    close_pairs = coordinates_KDTree.query_pairs(r=distance_threshold)
    close_pairs = [set(pair) for pair in close_pairs] # Convert to list of sets

    # Chain the pairs to obtain groups (or clusters) of points
    groups_of_points = combine_overlapping_sets(close_pairs)

    # Calculate the new coordinates by taking the center of all the neighbouring points.
    # A threshold for the total group is applied, i.e. all points must lie within the distance_threshold
    # from the center coordinate.
    new_coordinates = []
    for group in groups_of_points:
        group_coordinates = coordinates[list(group)]
        center_coordinate = np.mean(group_coordinates, axis=0)
        distances_to_center = np.sqrt(np.sum((group_coordinates-center_coordinate)**2, axis=1))
        if not (np.max(distances_to_center) > distance_threshold): # This could be another threshold
            new_coordinates.append(center_coordinate)

    # Obtain individual points, i.e. those that do not have another point within the distance_threshold.
    # This is done by taking the difference from all points and the ones that are present in any of the groups.
    all_points_in_groups = set(point for group in groups_of_points for point in group)
    all_points = set(range(len(coordinates)))
    individual_points = all_points.difference(all_points_in_groups)

    # Add individual points to new_coordinates list
    for point in individual_points:
        new_coordinates.append(coordinates[point])

    # Convert to numpy array
    new_coordinates = np.array(new_coordinates)

    if plot:
        axis = plt.figure().gca()
        axis.scatter(coordinates[:,0],coordinates[:,1])
        axis.scatter(new_coordinates[:,0],new_coordinates[:,1])

    return new_coordinates


def combine_overlapping_sets(old_list_of_sets):
    """ Combine sets that have overlap

    Go through each set, if it has overlap with one of the sets in the new list of sets, then combine it with this set
    If there is no overlap, append the set to the new list of sets.
    Perform this function recursively until the new_list_of_sets does not change anymore.

    Parameters
    ----------
    old_list_of_sets : list of sets
        List of sets of which overlapping ones should be combined

    Returns
    -------
    new_list_of_sets : list of sets
        Combined list of sets

    """

    # test_set1 = [set((1,2)),set((3,4)),set((5,2)),set((5,6)),set((4,10))]
    # test_set2 = [set((1,2)),set((3,4)),set((2,3)),set((5,6)),set((4,10))]

    new_list_of_sets = []
    for old_set in old_list_of_sets:
        append = True
        for new_set in new_list_of_sets:
            if not (old_set.isdisjoint(new_set)):
                new_set.update(old_set)
                append = False
        if append:
            new_list_of_sets.append(old_set.copy())

    if not (new_list_of_sets == old_list_of_sets):
        new_list_of_sets = combine_overlapping_sets(new_list_of_sets)

    return new_list_of_sets


def set_of_tuples_from_array(array):
    return set([tuple(a) for a in array])


def array_from_set_of_tuples(set_of_tuples):
    return np.array([t for t in set_of_tuples])

'''
utility functions (called within find_unique_molecules())
'''
# def check_new_molecules(image_to_check, unique_set, configs_peak_finding):
#     '''
#     from an input image this function finds the 'peaks'/molecule locations (built by others)
#     and keeps the unique set of molecules by comparing with the already found molecules
#
#
#     THE trick:
#     use the 'set operations' to get the union: "A U B"
#     '''
#     # --- use peak finding alogithm to find the locations of peaks/molecules in current image ---
#     new_molecule_locations = peak_finding.find_peaks(image=image_to_check, **configs_peak_finding)
#
#
#     # --- store the newly found molecules into a set ----
#     NewMolecules = set()
#     for new_loc in new_molecule_locations:
#         NewMolecules.add(tuple(new_loc))#.astype(int)))
#
#     # --- use set operations (like in statistics books) to get unique set of molecules w. Python ---
#     # ------ this next line just allows me to use the more intuitive name as the input parameter w/o
#     #  having to carry it  all the way through the function  -----------
#     UniqueMolecules = unique_set
#
#     # --- THE update: simply keep the union "A U B" ----
#     UniqueMolecules = NewMolecules.union(UniqueMolecules)
#     return UniqueMolecules

def coordinates_from_set(Set, channel, nmbr_pixels=512):
    '''
    my function used the 'set' type of Python to find the unique postions.
    this small function turns this back into a useable array

    note: in case of acceptor channel only coordinates --> shift to get coordinates w.r.t complete image

    '''
    coordinates = []
    for member in Set:
        coordinates.append(list(member))

    coordinates = np.array(coordinates)
    if channel in ['a', 'acceptor']:
        coordinates[:, 0] += (nmbr_pixels // 2)

    return coordinates




def find_FRET_pairs(Acceptor, Donor, uncertainty, nmbr_pixels):
    '''
    This algorithm is probably very inefficient, but it will get the job done ;)

    Should ask JooCs to help brainstorm of getting a different way of finding all the FRET pairs
    after including a tolerance of some pixels in both x and y directions

    first version: Do the brute force comparisson of one element in set A to an element in set D
    '''
    # --- the following just allows for the nicer name of the input variable
    # ------- without me needing to carry this over into the whole code  ---------
    u = uncertainty

    # --- transform/map donor coordinates onto acceptor channel ---
    A = coordinates_from_set(Acceptor, channel='donor', nmbr_pixels=nmbr_pixels)
    acceptor_transformed_coordinates = ctrans.transform(A)

    # --- turn into Python set-type to easily perform intersection operation ---
    # ---  set of points/peaks to find back in donor channel ----
    Acceptor_mapped = set()
    for coordinate in acceptor_transformed_coordinates:
        Acceptor_mapped.add(tuple(coordinate.astype(int)))

    # ---- using a tolerance of N pixels (2 by default), search for pairs in Donor channel ---
    # --- PROBABLY THE VERY INEFFICIENT PART -------
    FRETpairs = set()

    for a in Acceptor_mapped:
        for d in Donor:
            # --- now use the following to keep all useable pairs ---
            condition_x_coordinate = ((a[0] - u) <= d[0]) and (d[0] <= (a[0] + u))
            condition_y_coordinate = ((a[1] - u) <= d[1]) and (d[1] <= (a[1] + u))
            if condition_x_coordinate and condition_y_coordinate:
                FRETpairs.add(a)
                break  # No need to continue checking the Donor channel, already found a match

    return FRETpairs


def get_coordinate_pairs(file, peak_coordinates, channel):
    '''
    after having found the peak intensity locations in the image(s), we set Molecule() coordinates using the
    following structure:
    [ [array1: donor coordinates],[array2: acceptor coordinates]  ]
    First array has all the coordinates (either found directly or inferred) of the donors:
    donor_coordinates = [[x1,y1],[x2,y2], etc.]
    Similar for acceptor coordinates
    '''

    # --- Turn set of Peaks into an array
    # (translate the horizontal coordinates if peaks have only been found in acceptor channel) ----
    coordinates_identified = coordinates_from_set(peak_coordinates, channel=channel, nmbr_pixels=file.movie.width_pixels)


    # --- initialize ----
    Donors = []
    Acceptors = []


    for (x,y) in coordinates_identified:
        if x > (file.movie.width_pixels//2):
            '''
            this is a location in donor channel

            add this coordinate to Donors. 
            transform onto acceptor Channel and add to Acceptors
            '''
            donor = [x,y]
            acceptor  =  file.mapping.transform_coordinates(donor,direction='source2destination')


        else:
            '''
            this is a location in the acceptor channel
            
            add this location to Acceptors. 
            apply inverse transform to get location in donor channel and add to Donors
            '''
            acceptor = [x,y]
            donor = file.mapping.transform_coordinates(acceptor, direction='destination2source')

        # --- add to the lists ----
        Donors.append(donor)
        Acceptors.append(acceptor)

    # --- turn into arrays -------
    donor_coordinates = np.array(Donors)
    acceptor_coordinates = np.array(Acceptors)

    # --- output the coordinate array that is compatible with Molecule() class ---
    coordinates = np.hstack([donor_coordinates, acceptor_coordinates]).reshape((-1, 2))
    return coordinates