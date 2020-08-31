'''
Functions called in File().find_molecules()




'''
import numpy as np
from trace_analysis import peak_finding
from trace_analysis import coordinate_transformations as ctrans
import skimage as ski


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

        for i, channel_image in enumerate(channel_images):
            channel_coordinates = peak_finding.find_peaks(image=channel_image, **configs_peak_finding)#.astype(int)))
            channel_coordinates = set_of_tuples_from_array(channel_coordinates)
            coordinate_sets[i].update(channel_coordinates)


    #TODO: make this work properly using a the mapping transformation
    #TODO: make this usable for any number of channels
    if len(channels) > 1:
        raise NotImplementedError('Assessing found coordinates in multiple channels does not work properly yet')
        coordinates = combine_coordinate_sets(coordinate_sets, method='and') # the old detect_FRET_pairs

    # --- correct for photon shot noise / stage drift ---
    # Not sure whether to put this in front of combine_coordinate_sets/detect_FRET_pairs or behind [IS: 12-08-2020]
    coordinates = apply_uncertainty(coordinates, uncertainty=uncertainty_pixels, nmbr_pixels=nmbr_pixels)

    #TODO: map coordinates to main channel in movie

    # --- turn into array ---
    coordinates = array_from_set_of_tuples(coordinates)

    return coordinates, find_coords_img

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

def apply_uncertainty(Peaks, uncertainty=2, nmbr_pixels=512):
    '''
    due to photon shot noise + stage drift, the exact molecule position on the image can shift a bit

    first version:
    * pick a point from the full set of points found
    * determine difference in x and y coordinates to all other points
    * if ANY of the dx,dy pairs lies within the allowed uncertainty...
    * keep the test (x,y) and remove the ones that lie too close
    * In this fashion, the total loop should not need to run over all points.
    Once a peak is discarded because it is concidered "belonging to the same molecule", we don't need to use it as the
    reference (x,y) no longer.
    '''
    # --- the following just allows for the nicer name of the input variable
    # ------- without me needing to carry this over into the whole code  ---------
    u = uncertainty

    # --- make array of peak locations ---
    NewPeaks = list(Peaks.copy())
    PeakArray = coordinates_from_set(Peaks, channel='donor', nmbr_pixels=nmbr_pixels)
    xcoords = PeakArray[:, 0]
    ycoords = PeakArray[:, 1]

    removedIDs = []
    # --- choose a reference peak/molecule ---
    # --- needed to loop over the original set of peaks, otherwise loop ends before checking all pairs
    # --- because you reduce the size of the new set ----
    for ID, (x, y) in enumerate(Peaks):

        # --- If this pixel location allready has been discarded, no need to continue ---
        if ID in removedIDs:
            continue

        # --- determine if ANY of the remainig peaks can be concidered identical to the reference location ---
        xdiffs = np.abs(xcoords - x)
        ydiffs = np.abs(ycoords - y)
        condition_x_coordinate = xdiffs <= u
        condition_y_coordinate = ydiffs <= u

        # -- return the 'molecule'/'peak' indices that represent the same molecule/fluoresence spot ---
        no_need_these_peaks = np.where(condition_x_coordinate & condition_y_coordinate)[0]

        # --- remove the duplicates (keep the one used as reference) ---
        NoNeeds = no_need_these_peaks[no_need_these_peaks != ID]

        # --- to prevent the algorithm of trying to remove things twice ---
        NoNeeds = NoNeeds[~np.isin(NoNeeds, removedIDs)]

        for p in PeakArray[NoNeeds]:
            xNo = p[0]
            yNo = p[1]
            NewPeaks.remove((xNo, yNo))


        # --- to prevent the loop of checking something we allready know is going to be removed ---
        removedIDs.extend(NoNeeds)

    return set(NewPeaks)


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