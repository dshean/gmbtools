#! /usr/bin/env python

"""
Convert GSFC hdf GRACE mascon data to geopandas dataframe
Adapted from: https://github.com/NASA-Planetary-Science/HiMAT/blob/master/himatpy/GRACE_MASCON/pygrace.py
One-time script to prepare much smaller gpkg of mascons boundaries
"""

import sys
import os
import h5py
import pandas as pd
import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon
from pygeotools.lib import geolib 

def polygeom(mascon_s):
    x = np.array([-1,1,1,-1]) * mascon_s['lon_span'] / 2 + mascon_s['lon_center']
    y = np.array([-1,-1,1,1]) * mascon_s['lat_span'] / 2 + mascon_s['lat_center']
    return Polygon(list(zip(x, y)))

fn = sys.argv[1]

f = h5py.File(fn)
mascon_ds = f['mascon']

mascon_dct = {}
poly_geom = []

dataset_list = list(mascon_ds.keys())
for d in dataset_list:
    mascon_dct.update({ d: mascon_ds[d][0, :]})

mascon_df = pd.DataFrame.from_dict(mascon_dct)
for k, m in mascon_df.iterrows():
    poly_geom.append(polygeom(m))

#Need to clean up - the global df generates errors during spatial join, likely near 0,360 degrees
#Clipping to HMA resolves the issue

CRS = geolib.wgs_srs.ExportToProj4()
mascon_gdf = gpd.GeoDataFrame(mascon_df, crs=CRS, geometry=poly_geom)
mascon_gdf.index = mascon_gdf.index + 1
out_fn = os.path.splitext(fn)[0]+'.gpkg'
mascon_gdf.to_file(out_fn, driver='GPKG')
