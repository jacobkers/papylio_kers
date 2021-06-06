import numpy as np

# def make_gaussian(size, fwhm = 3, center=None):
#     # From https://stackoverflow.com/questions/7687679/how-to-generate-2d-gaussian-with-python
#     """ Make a square gaussian kernel.
#
#     size is the length of a side of the square
#     fwhm is full-width-half-maximum, which
#     can be thought of as an effective radius.
#     """
#
#     x = np.arange(0, size, 1, float)
#     y = x[:,np.newaxis]
#
#     if center is None:
#         x0 = y0 = size // 2
#     else:
#         x0 = center[0]
#         y0 = center[1]
#
#     return np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)

def make_gaussian(size, center=None, offset=[0,0]):
    x = np.arange(0, size, 1, float)
    y = x[:,np.newaxis]

    if center is None: center = [size // 2, size //2]

    return 2.0 * np.exp(- 0.3 * ((x-center[0]-offset[0])**2 + (y-center[0]-offset[1])**2))

def extract_trace_values_from_image(image, coordinates, background, twoD_gaussians):  # extract traces
    coordinates = np.atleast_2d(coordinates)

    # Probably indeed better to get this outside of the function, so that it is not redefined every time.
    half_size_Gaussian = len(twoD_gaussians[0]) // 2

    # This should likely be put on a central place in selection of locations
    # coordinates = coordinates[self.is_within_margin(coordinates, edge = None, margin = half_size_Gaussian + 1)]

    coordinates = np.floor(coordinates).astype(int)

    trace_values = np.zeros(len(coordinates))

    for i, coordinate in enumerate(coordinates):
        # First crop around spot, then do multiplication
        intensities = image[(coordinate[1] - half_size_Gaussian):(coordinate[1] + half_size_Gaussian + 1),
                      (coordinate[0] - half_size_Gaussian):(coordinate[0] + half_size_Gaussian + 1)
                      ]

        intensities = intensities - background[i]

        weighted_intensities = intensities * twoD_gaussians[i]
        trace_values[i] = np.sum(weighted_intensities)
        #trace_values[i]=np.sum(intensities) # MD testing
    return trace_values


def extract_traces(movie, coordinates, background=None, channel='all', gauss_width=4):
    # return donor and acceptor for the full data set
    #     root, name = os.path.split(self.filepath)
    #     traces_fn=os.path.join(root,name[:-4]+'-P.traces')

    # if os.path.isfile(traces_fn):
    # # load if traces file already exist
    #      with open(traces_fn, 'r') as infile:
    #          Nframes = np.fromfile(infile, dtype = np.int32, count = 1).item()
    #          Ntraces = np.fromfile(infile, dtype = np.int16, count = 1).item()
    #          rawData = np.fromfile(infile, dtype = np.int16, count = self.number_of_channels*Nframes * Ntraces)
    #      orderedData = np.reshape(rawData.ravel(), (self.number_of_channels, Ntraces//self.number_of_channels, Nframes), order = 'F')
    #      donor=orderedData[0,:,:]
    #      acceptor=orderedData[1,:,:]
    #      donor=np.transpose(donor)
    #      acceptor=np.transpose(acceptor)
    # else:

    # go through all images, extract donor and acceptor signal

    # This should likely be put on a central place in selection of locations
    # coordinates = coordinates[self.is_within_margin(coordinates, edge = None, margin = self.gauss_width // 2 + 1)]

    # donor=np.zeros(( self.number_of_frames,self.pts_number))
    # acceptor=np.zeros((self.number_of_frames,self.pts_number))

    traces = np.zeros((len(coordinates), movie.number_of_frames))

    if background is None:
        background = np.zeros(len(coordinates))

    # t0 = time.time()

    #twoD_gaussian = make_gaussian(gauss_width, fwhm=3, center=(gauss_width // 2, gauss_width // 2))

    offsets = coordinates % 1

    twoD_gaussians = [make_gaussian(gauss_width, offset=offsets[i]) for i in range(len(coordinates))]

    for frame_number in range(movie.number_of_frames):  # self.number_of_frames also works for pm, len(self.movie_file_object.filelist) not
        print(frame_number)
        image = movie.read_frame(frame_number)
        image = movie.get_channel(image, channel)
        trace_values_in_frame = extract_trace_values_from_image(image, coordinates, background, twoD_gaussians)

        traces[:,frame_number] = trace_values_in_frame  # will multiply with gaussian, spot location is not drift compensated
    # t1=time.time()
    # elapsed_time=t1-t0; print(elapsed_time)

    # root, name = os.path.split(self.filepath)

    # if os.path.isfile(trace_fn):

    return traces