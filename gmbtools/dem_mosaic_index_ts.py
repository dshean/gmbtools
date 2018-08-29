#! /usr/bin/env python

"""
Create timestamp arrays from dem_mosaic filename index arrays
"""

import os
import sys
from pygeotools.lib import iolib, timelib
import numpy as np

def make_dem_mosaic_index_ts(index_tif_fn):
    index_txt_fn = index_tif_fn+'-index-map.txt'
    out_fn = os.path.splitext(index_tif_fn)[0]+'_ts.tif'

    if not os.path.exists(index_tif_fn):
        print("Unable to find input file: %s" % index_tif_fn)
        return False 
    if not os.path.exists(index_txt_fn):
        print("Unable to find input file: %s" % index_txt_fn)
        return False
    if os.path.exists(out_fn):
        print("Existing ts output found: %s" % out_fn)
        return True

    #Load dem_mosaic index tif
    index_tif_ds = iolib.fn_getds(index_tif_fn)
    index_tif = iolib.ds_getma(index_tif_ds)

    #Create output array
    #index_ts_tif = np.zeros(index_tif.shape, dtype=np.float32)
    #index_ts_tif = np.ma.masked_all_like(index_tif)
    index_ts_tif = np.ma.masked_all(index_tif.shape, dtype=np.float32)

    #Read in dem_mosaic index txt file - should be "fn, index" for each record
    index_txt = np.genfromtxt(index_txt_fn, dtype=str)

    #Extract Python datetime object from filename (should pull out YYYYMMDD_HHMM)
    index_ts = [timelib.fn_getdatetime(fn) for fn in index_txt[:,0]]

    #Convert to desired output timestamp format
    #Python ordinal
    #index_ts = timelib.dt2o(index_ts)
    #YYYYMMDD integer
    #index_ts = [ts.strftime('%Y%m%d') for ts in index_ts]
    #Decimal year
    index_ts = [timelib.dt2decyear(ts) for ts in index_ts]

    for n, dt in enumerate(index_ts):
        index_ts_tif[index_tif == n] = dt

    iolib.writeGTiff(index_ts_tif, out_fn, index_tif_ds, ndv=0) 
    return True

def main():
    index_tif_fn = sys.argv[1]
    if not os.path.exists(index_tif_fn):
        sys.exit("Unable to find input file: %s" % index_tif_fn)
    result = make_dem_mosaic_index_ts(index_tif_fn)

if __name__ == "__main__":
    main()
