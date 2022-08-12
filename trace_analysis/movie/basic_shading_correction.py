#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from pathlib import Path
import cv2
import matplotlib.pyplot as plt
import numpy as np
import tqdm

from trace_analysis.correction.shading_correction import BaSiC


def BaSiCfun(directory, output_directory, extension='.tif', estimate_darkfield=False, apply_correction=False,
             use_flatfield=None, use_darkfield=None, epsilon=1e-6, l_s=None, l_d=None,
             output_flatfield_filename=None, output_darkfield_filename=None, verbose=False):
    # Create the BaSiC shading correction object
    optimizer = BaSiC(directory, estimate_darkfield=estimate_darkfield, extension=extension, verbose=verbose)

    # Set some optimizer parameters
    optimizer.l_s = l_s
    optimizer.l_d = l_d
    # Prepare the optimization
    optimizer.prepare()

    # Extract the flat and dark fields
    perform_estimation = True
    if use_flatfield is not None:
        img = cv2.imread(use_flatfield, cv2.IMREAD_ANYDEPTH)
        optimizer.set_flatfield(img)
        perform_estimation = False

    if use_darkfield is not None:
        img = cv2.imread(use_darkfield, cv2.IMREAD_ANYDEPTH)
        optimizer.set_flatfield(img)
        perform_estimation = False

    # Perform the estimation
    if perform_estimation:
        optimizer.run()

        # Save the estimated fields (only if the profiles were estimated)
        directory = Path(output_directory)
        directory.mkdir(parents=True, exist_ok=True)
        if output_flatfield_filename is not None:
            flatfield_name = Path(output_flatfield_filename).resolve()
            flatfield_name.parent.mkdir(parents=True, exist_ok=True)
        else:
            flatfield_name = directory / "flatfield.tif"
        if output_darkfield_filename is not None:
            darkfield_name = Path(output_darkfield_filename).resolve()
            darkfield_name.parent.mkdir(parents=True, exist_ok=True)
        else:
            darkfield_name = directory / "darkfield.tif"

        cv2.imwrite(str(flatfield_name), optimizer.flatfield_fullsize.astype(np.float32))
        cv2.imwrite(str(darkfield_name), optimizer.darkfield_fullsize.astype(np.float32))

    # Apply shading correction.
    if apply_correction:
        optimizer.write_images(output_directory, epsilon=epsilon)




def squeeze_channel_from_frames(frames):
    return xr.combine_by_coords(
        [frames.sel(channel=channel).set_index(x='x_pixel', y='y_pixel') for channel in frames.channel])

    # xr.combine_by_coords([frames.sel(channel=channel.index).set_index(x='x_pixel', y='y_pixel') for channel in self.channels])
    #
    # test = frames.assign_coords(channel_index_x=('channel',[0,1]),channel_index_y=('channel',[0,0]))
    # test2 = test.drop('channel').set_index(channel=('channel_index_x','channel_index_y'))
    # test3 = test2.unstack('channel').stack(x_pixel=('channel_index_x','x')).stack(y_pixel=('channel_index_y','y'))
    #
    # test2.unstack('channel').stack(channel=('channel_index_x','channel_index_y'))
    # test2.unstack('channel').stack(x_pixel=('channel_index_x','x'), y_pixel=('channel_index_y','y'))

def spatial_shading_correction(movies):
    frames = xr.concat([movie.read_frames_raw([0]) for movie in tqdm.tqdm(movies)], dim='frame')

    flatfield = xr.ones_like(frames.sel(frame=0), dtype=float)
    darkfield = xr.zeros_like(frames.sel(frame=0), dtype=float)

    for channel in frames.channel:
        optimizer = BaSiC(frames.sel(channel=channel).values, estimate_darkfield=True, extension=None, verbose=True)
        optimizer.prepare()
        optimizer.run()
        flatfield[dict(channel=channel)] = optimizer.flatfield_fullsize
        darkfield[dict(channel=channel)] = optimizer.darkfield_fullsize

    flatfield = squeeze_channel_from_frames(flatfield)
    darkfield = squeeze_channel_from_frames(darkfield)

    return darkfield, flatfield

#
# movies = files_green_laser[0::25].movie
# darkfield, flatfield = spatial_shading_correction(movies)
#
# import tifffile
# save_path = Path(r'N:\tnw\BN\CMJ\Shared\Ivo\PhD_data\20220602 - Objective-type TIRF (BN)')
# tifffile.imwrite(save_path.joinpath('flatfield.tif'), flatfield)
# tifffile.imwrite(save_path.joinpath('darkfield.tif'), darkfield)

# plt.figure()
# plt.imshow(flatfield)
#
# plt.figure()
# plt.imshow(darkfield)

#
# frames = files_green_laser[0].movie.read_frames_raw()
# frames = squeeze_channel_from_frames(frames)

if __name__ == '__main__':
    import tifffile
    from pathlib2 import Path
    import numpy as np

    pth = Path(r'N:\tnw\BN\CMJ\Shared\Ivo\PhD_data\20220602 - Objective-type TIRF (BN)\Analysis\Test')

    img_stack = tifffile.imread(pth.joinpath('TIRF 561 0001.tif'))
    img_stack = np.rot90(img_stack, 1, axes=(1,2))
    flatfield = tifffile.imread(pth.joinpath('flatfield.tif'))
    darkfield = tifffile.imread(pth.joinpath('darkfield.tif'))

    from trace_analysis.correction.shading_correction import get_photobleach

    import time
    start = time.time()
    c1 = get_photobleach(img_stack[:,:,0:256], flatfield[:,0:256], darkfield[:,0:256], size=(32,64))
    print(time.time()-start)
    start = time.time()
    c2 = get_photobleach(img_stack[:,:,256:512], flatfield[:,256:512], darkfield[:,256:512], size=(32,64))
    print(time.time()-start)

    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(c1.T/c1.max())
    plt.plot(c1.T*c2.max()/c2.T)
    plt.plot(c2.T-c1.T*c2.max()/c2.T)
    plt.plot(c2.max()/c2.T)
    plt.plot(c2.T/c2.max())
    plt.plot(test2/test2.max())




from trace_analysis.correction.shading_correction import get_photobleach
import time
start = time.time()
test = get_photobleach(frames.values, flatfield, darkfield=darkfield)
print(time.time()-start)

import scipy.ndimage
start = time.time()
# test2 = np.array([scipy.ndimage.minimum_filter(((frame)), size=15, mode='wrap').sum() for frame in img_stack])
test2 = np.array([scipy.ndimage.gaussian_filter(((frame-darkfield)/flatfield)[:,0:256], sigma=50, mode='wrap').mean() for frame in img_stack])
test2 = np.array([scipy.ndimage.minimum_filter(((frame-darkfield)/flatfield)[:,0:256], size=15, mode='wrap').mean() for frame in img_stack])
# test2 = np.array([frame.sum() for frame in frames])
# test2 = np.array([scipy.ndimage.minimum_filter(frame, size=15, mode='wrap').sum() for frame in img_stack])
print(time.time()-start)

test = test.squeeze()
test = test/test.sum()
test2 = test2/test2.sum()

plt.figure()
plt.plot(test)
plt.plot(test2)





