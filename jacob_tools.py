""" Here I collect some basic tools that might be of use in FRET analysis. JacobKers'25 """
# -*- coding: utf-8 -*-
"""
spot image
Created on Sep 5 2023
@author: jkerssemakers
"""

import matplotlib.pyplot as plt
import numpy as np
import math as mt


def sorted_pixels_treshold(im, demo=0):
    rr,cc=im.shape
    # sort and scale on number of pixels (to equalize axes)
    impixels=im.flatten()
    impixels_sorted=np.sort(impixels)
    Npix=len(impixels)
    pix_ax=np.arange(0, Npix, 1)
    Ipix=np.max(impixels)
    impixels_sorted=impixels_sorted/Ipix*Npix

    #fit on lower half of N:
    lowerhalf_N=impixels_sorted[0:int(Npix/2)]
    lowerhalf_pix_ax=pix_ax[0:int(Npix/2)]
    lowerfit_p=np.polyfit(lowerhalf_pix_ax,lowerhalf_N,1)
    lowerfit=np.polyval(lowerfit_p, pix_ax)

    #fit on higher half of I:
    upperhalf_I=impixels_sorted[impixels_sorted>Npix/2]
    upperhalf_pix_ax=pix_ax[impixels_sorted>Npix/2]
    upperfit_p=np.polyfit(upperhalf_pix_ax,upperhalf_I,1)
    upperfit=np.polyval(upperfit_p, pix_ax)
    #get cross-point
    xc=(lowerfit_p[1]-upperfit_p[1])/(upperfit_p[0]-lowerfit_p[0])
    yc=np.polyval(lowerfit_p,xc)

    #get 'knee'
    #rr=np.hypot((1:length(sim))-xc).'.^2, (sim-yc).^2);
    rr=np.hypot(pix_ax-xc,impixels_sorted-yc)
    x_kn=pix_ax[(rr== min(rr))]
    y_kn=impixels_sorted[(rr== min(rr))]
    #scale value back
    treshold=y_kn/Npix*Ipix
    msk=im>treshold
    spot_tres=msk*im
    
    if demo==1:
        fig, ax = plt.subplots(1, 3)
        ax[0].imshow(im)
        ax[1].plot(pix_ax,impixels_sorted)
        ax[1].plot(lowerhalf_pix_ax,lowerhalf_N,'o')
        ax[1].plot(upperhalf_pix_ax,upperhalf_I, 'o')
        ax[1].plot(pix_ax,lowerfit, 'r-')
        ax[1].plot(pix_ax,upperfit, 'r-')
        ax[1].plot(xc,yc, 'ko')
        ax[1].plot(x_kn,y_kn, 'mo')
        ax[1].set_ylim(np.min(impixels_sorted),np.max(impixels_sorted) )
        ax[2].imshow(spot_tres)
        fig.tight_layout()
        return fig,ax
    else:
        return spot_tres, msk, treshold 


def outlier_flag(data=0, tolerance=2.5, sig_change=0.7, how=1, demo=0):
    """
    An iterative tool to separate a distribution from its outliers.
    An initial estimate of average 'mu' and standard deviation 'sigma' is used to identify outliers.
    These are removed and [mu,sigma] is re-determnined] after wchich the sequence is repeated
    until sigma does not change much anymore
    Input:
        data: single array of values
        tolerance: outliers are points more than [tolerance] standard deviations away from the average.
        sigchange: iteration stops if sigma is changed less than a fraction  'sigchange'.
        how
    sho
    demo
    Output:
        1) array of outliers
        2) array of inliers
        3) flags: binary trace indicating which points where outliers in the original data.
    Jacob Kers '2022
    """
    
    sig_ratio = 0
    sigma_nw = 1e20
    flags = 0 * data + 1
    while sig_ratio < sig_change:
        sigma_old = sigma_nw
        ix_in = np.ndarray.nonzero(flags == 1)
        ix_out = np.ndarray.nonzero(flags == 0)
        inliers = data[ix_in]
        outliers = data[ix_out]
        av = np.median(inliers)
        sigma_nw = np.std(inliers)
        sig_ratio = sigma_nw / sigma_old
        if how == 1:
            flags = (data - av) < tolerance * sigma_nw
        elif how == 0:
            flags = abs(data - av) < tolerance * sigma_nw
        elif how == -1:
            flags = (data - av) > -tolerance * sigma_nw
        if demo:
            lo = np.min(inliers)
            hi = np.max(inliers)
            bins = np.linspace(lo, hi, 40)
            fig2, ax2 = plt.subplots(1, 1)
            ax2.hist(inliers, bins, histtype="bar")
            fig2.show()
        else:
            return inliers, outliers, flags



def generate_spot(
    psf=2,  # point spread
    size=50,  # image size (square shape)
    amplitude=1,  # amplitude
    SN_ratio=0.1,  # signal/noise ratio
    demo=0,  # demo plot or not
    node_x=0,  # relative center position of the spot (0=center)
    node_y=0,  # relative center position of the spot (1=edge)
):
    """
    This generates a spot for quick testing.
    coordinates image run from -10 to 10
    """
    x, y = np.linspace(-10 - node_x * 10, 10 - node_x * 10, size), np.linspace(
        -10 - node_y * 10, 10 - node_y * 10, size
    )
    X, Y = np.meshgrid(x, y)
    rr = np.hypot(X, Y)
    frame = amplitude * np.exp((-(rr/ (2*psf))** 2))
    noise = SN_ratio * amplitude * np.random.randn(size, size)
    spot_im = frame + noise
    # demo section start ------------------------
    if demo > 0:
        fig, ax = plt.subplots(1, 1)
        ax.imshow(spot_im)
        ax.set_title("pattern")
        ax.set_xlabel("x (pixels)")
        ax.set_ylabel("y (pixels)")
        return fig, ax
    # demo section stop  ------------------------
    else:
        return spot_im

def random_data_twopeaks(demo=0):
    # build trace with 2 distributions
    N_pts = 2000
    s1 = 10
    u1 = 0
    N_otl = 200
    s2 = 10
    u2 = 100
    temp_ax = np.arange(0, N_pts, 1)
    data = s1 * np.random.randn(N_pts, 1) + u1
    for ii in range(N_otl):
        randii = int((N_pts - 1) * np.random.rand(1, 1))
        data[randii] = s2 * np.random.randn(1, 1) + u2
    if demo > 0:
        fig, ax = plt.subplots(1, 1)
        ax.plot(data)
        ax.set_title("pattern")
        ax.set_xlabel("x (a.u)")
        ax.set_ylabel("y (a.u)")
        return fig, ax
    # demo section stop  ------------------------
    else:
        return data    
