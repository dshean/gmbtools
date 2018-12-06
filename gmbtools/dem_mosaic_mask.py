#! /usr/bin/env python

"""
Script to mask output from make_mos.sh and dem_mosaic_validtiles.py to isolate good elevation values
"""

import sys
import os
import glob
from pygeotools.lib import iolib
import numpy as np

def domask(tile_fn): 
    prefix = '-'.join(tile_fn.split('-')[:-1])
    print("\nLoading: %s" % tile_fn)
    #dem_fn = prefix+'-median.tif'
    dem_fn = tile_fn
    dem_ds = iolib.fn_getds(dem_fn)
    dem = iolib.ds_getma(dem_ds)

    #Get original mask, True where masked
    mask = np.ma.getmaskarray(dem)
    valid_px_count = (~mask).sum()
    if valid_px_count < 1:
        sys.exit("No valid pixels remain")
    print("Valid pixel count: %i" % valid_px_count)

    min_count = 2
    count_fn = prefix+'-count.tif'
    print("Loading: %s" % count_fn)
    count = iolib.fn_getma(count_fn)
    print("min: %i max: %i" % (count.min(), count.max()))
    print("Masking: (count < %i)" % min_count)
    mask = np.logical_or(mask, (count < min_count))
    valid_px_count = (~mask).sum()
    if valid_px_count < 1:
        sys.exit("No valid pixels remain")
    print("Valid pixel count: %i" % valid_px_count)

    max_std = 3.0
    #std_fn = prefix+'-std.tif'
    std_fn = prefix+'-nmad.tif'
    print("Loading: %s" % std_fn)
    std = iolib.fn_getma(std_fn)
    print("min: %i max: %i" % (std.min(), std.max()))
    print("Masking: (std/nmad >= %i)" % max_std)
    mask = np.logical_or(mask, (std >= max_std))
    valid_px_count = (~mask).sum()
    if valid_px_count < 1:
        sys.exit("No valid pixels remain")
    print("Valid pixel count: %i" % valid_px_count)

    print("Applying mask")
    dem_masked = np.ma.array(dem, mask=mask)
    out_fn = os.path.splitext(dem_fn)[0]+'_masked.tif'
    print("Writing: %s" % out_fn)
    iolib.writeGTiff(dem_masked, out_fn, dem_ds)

def main():
    if len(sys.argv) != 2:
        sys.exit("Usage: '%s tile_fn'" % sys.argv[0])
    else:
        tile_fn = sys.argv[1]
        if os.path.exists(tile_fn):
            domask(tile_fn)
        else:
            sys.exit("Unable to find input: %s" % tile_fn)

if __name__== "__main__":
    main()
