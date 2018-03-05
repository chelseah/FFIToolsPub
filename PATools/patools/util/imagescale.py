# 
# Copyright (C) 2017 - Massachusetts Institute of Technology (MIT) 
# 
# This program is free software: you can redistribute it and/or modify 
# it under the terms of the GNU General Public License as published by 
# the Free Software Foundation, either version 3 of the License, or 
# (at your option) any later version. 
# 
# This program is distributed in the hope that it will be useful, 
# but ITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details. 
# 
# You should have received a copy of the GNU General Public License 
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 



import numpy as np


def zscale(image, nsamples=1000, contrast=0.25, max_reject=0.5, min_npixs=5, max_iters=5, krej=2.5):
    """
    implementation of IRAF zscale algorithm by module numdisplay from stsdas.stsci.edu
    """
    samples = zsc_sample(image, nsamples)
    npix = len(samples)
    center_pixel = (npix - 1) // 2 
    zmin = np.nanmin(samples)
    zmax = np.nanmax(samples)
    median = np.nanmedian(samples)

    samples.sort()

    minpix = max(min_npixs, int(npix*max_reject))
    ngrow = max (1, int (npix * 0.01))
    ngoodpix, zstart, zslope = zsc_fit_line (samples, npix, krej, ngrow, max_iters, minpix)
    if ngoodpix < minpix: 
        z1 = zmin 
        z2 = zmax 
    else: 
        if contrast > 0: zslope = zslope / contrast 
        z1 = max (zmin, median - (center_pixel - 1) * zslope) 
        z2 = min (zmax, median + (npix - center_pixel) * zslope) 
    return z1, z2

def zsc_sample (image, maxpix): 
    """ Figure out which pixels to use for the zscale algorithm 
    Returns the 1-d array samples.  
    Sample in a square grid, and return the first maxpix in the sample """
    nc = image.shape[0] 
    nl = image.shape[1] 
    stride = max (1.0, np.sqrt((nc - 1) * (nl - 1) / float(maxpix))) 
    stride = int (stride) 
    samples = image[::stride,::stride].flatten() 
    return samples[:maxpix]

def zsc_fit_line (samples, npix, krej, ngrow, maxiter, minpix): 
    # 
    # First re-map indices from -1.0 to 1.0 
    xscale = 2.0 / (npix - 1) 
    xnorm = np.arange(npix) 
    xnorm = xnorm * xscale - 1.0 
  
    ngoodpix = npix 
    last_ngoodpix = npix + 1 
  
    # This is the mask used in k-sigma clipping.  0 is good, 1 is bad 
    goodpixels = samples>0 
  
    # 
    #  Iterate 
  
    for niter in range(maxiter): 
  
        if (ngoodpix >= last_ngoodpix) or (ngoodpix < minpix): 
            break 
         
        # Accumulate sums to calculate straight line fit 
        # goodpixels = np.where(badpix == 0) 
        sumx = xnorm[goodpixels].sum() 
        sumxx = (xnorm[goodpixels]*xnorm[goodpixels]).sum() 
        sumxy = (xnorm[goodpixels]*samples[goodpixels]).sum() 
        sumy = samples[goodpixels].sum() 
        length = len(xnorm[goodpixels]) 
  
        delta = length * sumxx - sumx * sumx 
        # Slope and intercept 
        intercept = (sumxx * sumy - sumx * sumxy) / delta 
        slope = (length * sumxy - sumx * sumy) / delta 
         
        # Subtract fitted line from the data array 
        fitted = xnorm*slope + intercept 
        flat = samples - fitted 
  
        # Compute the k-sigma rejection threshold
        mean = np.nanmean(flat[goodpixels])
        sigma = np.nanstd(flat[goodpixels])
  
        threshold = sigma * krej 
  
        # Detect and reject pixels further than k*sigma from the fitted line 
        goodpixels = np.abs(flat)<threshold 
         
        niter += 1 
  
    # Transform the line coefficients back to the X range [0:npix-1] 
    zstart = intercept - slope 
    zslope = slope * xscale 
  
    return ngoodpix, zstart, zslope 
