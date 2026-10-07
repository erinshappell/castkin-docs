################################################################################
# castkin.py
# Written by Erin Shappell for Lu Lab
# 
# This module contains all functions used to generate continuous estimates
# of the angular velocity due to head casting in freely moving C. elegans
#
################################################################################

# Imports
import os
import cv2
import glob

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path
from IPython.display import clear_output
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import splev, splrep, interp1d
from scipy.stats import levene, median_abs_deviation

################################################################################
def remove_nans(arr):
    """
    Remove NaN values from an array.

    Inputs:
        arr (array-like): Input array that may contain NaN values.

    Output:
        ndarray: An array containing only the non-NaN values from the input array.
        
    """
    return arr[~np.isnan(arr)]

################################################################################
def smooth(data, sigma):
    """
    Apply one-dimensional Gaussian smoothing to a signal.

    Inputs:
        data (array-like): Input signal to smooth.
        sigma (float):     Standard deviation of the Gaussian smoothing kernel.

    Output:
        ndarray: The smoothed version of the input signal.
        
    """
    data_smoothed = gaussian_filter1d(data, sigma, mode='nearest')
    return data_smoothed

################################################################################
def filter_data(data, threshold=100, interpolate_gaps=False, recovery_window=5,
                confirm_steps=2, anchor_window=10, anchor_min_count=3):
    """
    Remove large discontinuities and tracking errors from a signal.

    Signal values that differ from the most recent valid value by more than
    the specified threshold are marked as invalid. Optionally, invalid
    regions can be filled using interpolation.

    Inputs:
        data (array-like):                 Input signal.
        threshold (float, optional):       Maximum allowed change between consecutive
                                           valid samples.
        interpolate_gaps (bool, optional): If True, interpolate across invalid regions.
        recovery_window (int, optional):   Number of samples ahead to search for signal
                                           recovery.
        confirm_steps (int, optional):     Number of consecutive stable samples required
                                           to confirm recovery.
        anchor_window (int, optional):     Number of samples used when searching for a
                                           valid starting anchor.
        anchor_min_count (int, optional):  Minimum number of nearby points required to
                                           establish an anchor.

    Output:
        ndarray: Filtered signal with invalid regions removed or interpolated.
        
    """
    data = np.asarray(data, dtype=float)
    data_filtered = data.copy()

    ### Find a trusted starting anchor via consensus
    # Look for a cluster of `anchor_min_count` points within `anchor_window`
    # that are all mutually within threshold of a candidate anchor value.
    start        = None
    anchor_value = None

    for i in range(min(anchor_window, len(data))):
        candidate = data[i]
        # Count how many points in the window are close to this candidate
        window = data[i: i + anchor_window]
        close_count = np.sum(np.abs(window - candidate) <= threshold)
        if close_count >= anchor_min_count:
            start = i
            anchor_value = candidate
            break

    if start is None: return np.full_like(data, np.nan)

    data_filtered[:start] = np.nan
    last_valid = anchor_value
    i = start + 1

    while i < len(data_filtered):
        if abs(data[i] - last_valid) > threshold:
            # Mark as invalid and search ahead for recovery
            data_filtered[i] = np.nan
            recovered = False

            for offset in range(1, recovery_window + 1):
                j = i + offset
                if j >= len(data): break

                # Candidate must be close to last known valid value
                if abs(data[j] - last_valid) > threshold: continue

                # Confirm recovery: next `confirm_steps` must also be stable
                stable = True
                for k in range(1, confirm_steps + 1):
                    if j + k >= len(data): break
                    if abs(data[j + k] - data[j + k - 1]) > threshold:
                        stable = False
                        break

                if stable:
                    data_filtered[i:j] = np.nan
                    last_valid = data[j]
                    i = j
                    recovered = True
                    break

            if not recovered: pass
        else: last_valid = data[i]
        i += 1

    # Optional interpolation
    if interpolate_gaps:
        indices = np.arange(len(data_filtered))
        valid = ~np.isnan(data_filtered)
        if valid.sum() > 1:
            data_filtered = np.interp(indices, indices[valid], data_filtered[valid])

    return data_filtered

################################################################################
def spatial_curv(x, y, anginc):
    """
    Calculate curvature and bending angle along a body contour across multiple frames.

    Curvature is estimated by fitting circles through neighboring points along
    the contour and computing the inverse radius of curvature.
    
    NOTE: this function was converted to Python from MATLAB and originally written
    by Chris Pierce and Lucinda Peng for Lu Lab. 

    Inputs:
        x (ndarray):  X-coordinates of body points with shape (num_frames, num_points).
        y (ndarray):  Y-coordinates of body points with shape (num_frames, num_points).
        anginc (int): Point spacing used when calculating local curvature and bending angle.

    Outputs:
        ndarray: Curvature values for each body point and frame.
        ndarray: Bending angles for each body point and frame.
        
    """
    # Transpose the data
    x = x.T
    y = y.T
    
    # Get dimnesions from data
    [numpts,numframes] = x.shape 

    Curvature = np.zeros((numpts - anginc, numframes))
    Angles    = np.zeros((numpts - anginc, numframes))
    
    # Loop for calculating curvature
    for m in range(numframes):
        tempx = x[:,m]
        tempy = y[:,m]

        # angle, radius and curvature all along the track
        # Calculate using formula for a fit circle
        SplineAngles    = np.zeros((numpts-anginc))
        SplineRadius    = np.zeros((numpts-anginc))
        SplineCurvature = np.zeros((numpts-anginc))
        
        for i in range(anginc,numpts-anginc):
            pasttopresentvector = np.array([tempx[i]-tempx[i-anginc],tempy[i]-tempy[i-anginc]])
            presenttonextvector = np.array([tempx[i+anginc]-tempx[i],tempy[i+anginc]-tempy[i]])

            dotprod          = np.dot(pasttopresentvector, presenttonextvector)
            dotprod          = dotprod / (np.linalg.norm(pasttopresentvector) * np.linalg.norm(presenttonextvector))
            SplineAngles[i]  = np.arccos(dotprod)
            crossprod        = np.cross(pasttopresentvector,presenttonextvector)

            pastpresslope = (tempy[i] - tempy[i-anginc]) / (tempx[i] - tempx[i-anginc])
            presfutslope  = (tempy[i+anginc]-tempy[i])/(tempx[i+anginc]-tempx[i]);
            
            circx = (pastpresslope * presfutslope * (tempy[i-anginc] - tempy[i+anginc])+\
                     presfutslope  * (tempx[i-anginc] + tempx[i])-\
                     pastpresslope * (tempx[i] + tempx[i+anginc])) / (2*(presfutslope-pastpresslope))
            circy = (-1*(circx-(tempx[i-anginc] + tempx[i])/2) / pastpresslope) + (tempy[i-anginc]+tempy[i])/2
            
            SplineRadius[i]    = np.power(np.power(tempx[i]-circx,2)+np.power(tempy[i]-circy,2),0.5)
            SplineCurvature[i] = 1 / SplineRadius[i]
            
            if crossprod < 0:
                SplineRadius[i]    = -1*SplineRadius[i]
                SplineCurvature[i] = -1*SplineCurvature[i]
                SplineAngles[i]    = -1*SplineAngles[i]
                
        Curvature[:,m] = SplineCurvature
        Angles[:,m]    = SplineAngles
    
    return [Curvature, Angles]

################################################################################
def calc_angle(v1, v2, origin=None):
    """
    Calculate the angle between two vectors in degrees.

    Vectors may share an arbitrary origin point. If an origin is supplied,
    the vectors are first translated so that the origin becomes (0, 0, ...).

    Inputs:
        v1 (array-like):               Endpoint of the first vector.
        v2 (array-like):               Endpoint of the second vector.
        origin (array-like, optional): Origin point shared by both vectors.
                                       If None, assumed to be (0, 0, ...).

    Output:
        float: Angle between v1 and v2 in degrees.
        
    """
    v1 = np.array(v1, dtype=float)
    v2 = np.array(v2, dtype=float)

    if origin is not None:
        origin = np.array(origin, dtype=float)
        v1 = v1 - origin
        v2 = v2 - origin

    dot_product = np.dot(v1, v2)
    norm_product = np.linalg.norm(v1) * np.linalg.norm(v2)

    cos_theta = np.clip(dot_product / norm_product, -1.0, 1.0)
    angle_rad = np.arccos(cos_theta)

    # Use 2D cross product (z-component) to determine direction
    cross_z = v1[0] * v2[1] - v1[1] * v2[0]

    if cross_z < 0: angle_rad = 2 * np.pi - angle_rad

    angle_deg = np.degrees(angle_rad)

    return angle_deg