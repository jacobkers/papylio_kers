#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from pathlib import Path
import cv2
import matplotlib.pyplot as plt
import numpy as np
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


import tqdm

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

    flatfield = xr.ones_like(frames.sel(frame=0))
    darkfield = xr.zeros_like(frames.sel(frame=0))

    for channel in frames.channel:
        optimizer = BaSiC(frames.sel(channel=channel).values, estimate_darkfield=True, extension=None, verbose=True)
        optimizer.prepare()
        optimizer.run()
        flatfield[dict(channel=channel)] = optimizer.flatfield_fullsize
        darkfield[dict(channel=channel)] = optimizer.darkfield_fullsize

    flatfield = squeeze_channel_from_frames(flatfield)
    darkfield = squeeze_channel_from_frames(darkfield)

    return darkfield, flatfield


movies = files_green_laser.movie
darkfield, flatfield = spatial_shading_correction(movies)



plt.figure()
plt.imshow(flatfield)

plt.figure()
plt.imshow(darkfield)




OPTIONS = {'size': 128}

img_stack = np.zeros((OPTIONS["size"], OPTIONS["size"], frames.shape[0]))

for i, frame in enumerate(frames):
    img_stack[:, :, i] = cv2.resize(frame, (OPTIONS["size"], OPTIONS["size"]), interpolation=cv2.INTER_LINEAR,).astype(np.float64)

flatfield_small = cv2.resize(flatfield, (OPTIONS["size"], OPTIONS["size"]), interpolation=cv2.INTER_LINEAR,).astype(np.float64)
darkfield_small = cv2.resize(darkfield, (OPTIONS["size"], OPTIONS["size"]), interpolation=cv2.INTER_LINEAR,).astype(np.float64)

from numba import njit


@njit
def mean_axis_0(arr):
    return np.array([arr[..., i].mean() for i in range(arr.shape[1])])

def get_photobleach(imgflt_stack, flatfield, darkfield=None):
    size=128
    imgflt_stack = np.reshape(
        imgflt_stack, (size * size, -1)
    ).astype(np.float64)

    imgflt_stack_svd = np.linalg.svd(imgflt_stack, full_matrices=False, compute_uv=False)

    return _get_photobleach(imgflt_stack, imgflt_stack_svd, flatfield, darkfield=darkfield)

# From https://github.com/PolusAI/polus-plugins/tree/dev/regression/polus-basic-flatfield-correction-plugin
@njit
def _get_photobleach(imgflt_stack, imgflt_stack_svd, flatfield, darkfield=None):
    """Calculate the global effect of photobleaching for each image
    Using the original data, flatfield, and darkfield images, estimate the total
    contribution of photobleaching to an image in a series of images.
    Inputs:
        imgflt_stack - Numpy stack of images
        flatfield - numpy floating precision matrix containing flatfield values
        darkfield - numpy floating precision matrix containing darkfield values
    Outputs:
        A_coeff - A 1xn matrix of photobleaching offsets, where n is the number
            of input images
    """
    size=128

    # Initialize matrices
    # Initialize matrices
    # imgflt_stack = np.reshape(
    #     imgflt_stack, (size * size, -1)
    # ).astype(np.float64)
    if darkfield is None:
        darkfield = np.zeros(flatfield.shape, dtype=np.float64)

    # Initialize weights and tolerances
    weights = np.ones(imgflt_stack.shape, dtype=np.float64)
    epsilon = np.float64(0.1)
    tol = np.float64(10 ** -6)

    # Run optimization exactly 5 times
    for r in range(5):
        print(r)
        # Calculate weights, offsets and coefficients
        W_idct_hat = np.reshape(flatfield, (-1, 1))
        A_offset = np.reshape(darkfield, (-1, 1))
        A_coeff = np.reshape(mean_axis_0(imgflt_stack), (1, -1))

        # Initialization values and learning rates
        # temp = np.linalg.svd(imgflt_stack, full_matrices=False, compute_uv=False)
        temp = imgflt_stack_svd
        norm_two = np.float64(temp[0])
        mu = np.float64(12.5) / norm_two
        mu_bar = mu * 10 ** 7
        rho = np.float64(1.5)

        ent1 = 1
        wem = weights / (ent1 * mu)

        # Normalization factors
        d_norm = np.linalg.norm(imgflt_stack)#, "fro")

        # Initialize augmented representation and error
        A = np.zeros(imgflt_stack.shape, dtype=np.float64)
        E1 = np.zeros(imgflt_stack.shape, dtype=np.float64)
        # Y1 = np.float64(0)
        Y1 = np.zeros(imgflt_stack.shape)

        # Run optimization
        iternum = 0
        converged = False
        while not converged:
            iternum += 1
            print(iternum)
            # Calculate augmented representation
            # A = np.matmul(W_idct_hat, A_coeff) + A_offset
            A = (W_idct_hat * A_coeff) + A_offset

            # Calculate errors
            # E1 = E1 + np.divide(imgflt_stack - A - E1 + np.multiply(1 / mu, Y1), ent1)
            E1 = imgflt_stack - A + np.multiply(1 / mu, Y1)

            # E1 = np.max(
            #     np.reshape(
            #         E1 - weights / (ent1 * mu),
            #         (imgflt_stack.shape[0], imgflt_stack.shape[1], 1),
            #     ),
            #     -1,
            #     initial=10 ** -6,
            # ) + np.min(
            #     np.reshape(
            #         E1 + weights / (ent1 * mu),
            #         (imgflt_stack.shape[0], imgflt_stack.shape[1], 1),
            #     ),
            #     -1,
            #     initial=0,
            # )

            E1 = np.maximum(E1 - wem, 10 ** -6) + np.minimum(E1 + wem, 0)

            # Calculate coefficients
            R1 = imgflt_stack - E1
            # A_coeff = np.reshape(np.mean(R1, axis=0), (1, -1)) - np.mean(A_offset)
            A_coeff = np.reshape(mean_axis_0(R1), (1, -1)) - np.mean(A_offset)
            # A_coeff[A_coeff < 0] = 0  # pixel values are never negative
            A_coeff = np.maximum(A_coeff, 0)

            # Loss
            Z1 = imgflt_stack - A - E1

            # Error updates
            Y1 = Y1 + mu * Z1

            # Update learning rate
            # mu = np.min(mu * rho, initial=mu_bar)
            mu = np.minimum(mu * rho, mu_bar)

            # Stop if below threshold
            #stopCriterion = np.linalg.norm(Z1, "fro") / d_norm
            stopCriterion = np.linalg.norm(Z1) / d_norm

            if stopCriterion < tol:
                converged = True

        # Update weights
        # XE_norm = np.reshape(np.mean(A, axis=0), (1, -1)) / E1
        XE_norm = np.reshape( mean_axis_0(A), (1, -1)) / E1
        weights = 1 / np.abs(XE_norm + epsilon)
        weights = weights * weights.size / np.sum(weights)

    return A_coeff

import time
start = time.time()
test = get_photobleach(img_stack, flatfield_small, darkfield=darkfield_small)
print(time.time()-start)

import scipy.ndimage
start = time.time()
test2 = np.array([scipy.ndimage.minimum_filter(np.array(frame), size=15, mode='wrap').sum() for frame in frames])
# test2 = np.array([scipy.ndimage.gaussian_filter(np.array(frame), sigma=15, mode='wrap').sum() for frame in frames])
# test2 = np.array([frame.sum() for frame in frames])
print(time.time()-start)

test = test.squeeze()
test = test/test.sum()
test2 = test2/test2.sum()

plt.figure()
plt.plot(test)
plt.plot(test2)





