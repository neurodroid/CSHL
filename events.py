import numpy as np
import sys

import filter

def find_peaks(data, dt, threshold, mode='epsc', f_lo=0.001, f_hi=0.02):
    
    databp = filter.bandpass(data, dt, f_lo, f_hi)
    if mode is 'epsc':
        crossings_start = np.where(
            np.diff((databp < threshold).astype(np.float)) == 1)[0]
        crossings_end = np.where(
            np.diff((databp < threshold).astype(np.float)) == -1)[0]
    else:
        crossings_start = np.where(
            np.diff((databp > threshold).astype(np.float)) == 1)[0]
        crossings_end = np.where(
            np.diff((databp > threshold).astype(np.float)) == -1)[0]


    if not len(crossings_start) or not len(crossings_end):
        sys.stderr.write("No threshold crossings detected. Try to lower threshold.")
        return None

    if crossings_start[0] > crossings_end[0]:
        crossings_start = crossings_start[1:]

    if len(crossings_start) > len(crossings_end):
        crossings_start = crossings_start[:-1]

    if len(crossings_end) > len(crossings_start):
        crossings_end = crossings_end[:-1]

    ioffset = 20 # account for delay of butterworth filter
    if mode is 'epsc':
        return np.array([
            np.argmin(data[start-ioffset:end]) + start-ioffset
            for start, end in zip(crossings_start, crossings_end)])
    else:
        return np.array([
            np.argmax(data[start-ioffset:end]) + start-ioffset
            for start, end in zip(crossings_start, crossings_end)])
