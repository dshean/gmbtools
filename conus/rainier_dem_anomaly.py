#! /usr/bin/env python

"""
Create anomaly map time series
"""

import os
import sys
import glob

import matplotlib.pyplot as plt

from pygeotools.lib import iolib, timelib
from imview.lib import pltlib

def makefig(dem, hs, anomaly, ds, title=None):
    f,axa = plt.subplots(1,2,figsize=(10,5))
    dem_clim = (2300, 4200)
    hs_clim = (1, 255)
    anomaly_clim = (-15, 15)
    hs_im = axa[0].imshow(hs, vmin=hs_clim[0], vmax=hs_clim[1], cmap='gray')
    dem_im = axa[0].imshow(dem, vmin=dem_clim[0], vmax=dem_clim[1], cmap='cpt_rainbow', alpha=0.5)
    res = 8
    pltlib.add_scalebar(axa[0], res=res)
    pltlib.add_cbar(axa[0], dem_im, label='Elevation (m WGS84)')
    anomaly_im = axa[1].imshow(anomaly, vmin=anomaly_clim[0], vmax=anomaly_clim[1], cmap='RdBu')
    pltlib.add_cbar(axa[1], anomaly_im, label='Elevation Anomaly (m)')
    pltlib.shp_overlay(axa[1], ds, shp_fn, color='darkgreen')
    plt.tight_layout()
    for ax in axa:
        pltlib.hide_ticks(ax)
        ax.set_facecolor('k')
        if title is not None:
            ax.set_title(title)
    return f

dem_ref_fn = 'rainier_allgood_mos-tile-0_warp.tif'
dem_ref = iolib.fn_getma(dem_ref_fn)
dem_fn_list = glob.glob('*8m_trans_warp.tif')
shp_fn = '/Volumes/SHEAN_1TB_SSD/usgs_dems/rainier/final_clip/rainier_24k_1970-2015_mb_lines.shp'

outdir = 'movie_2panel'
if not os.path.exists(outdir):
    os.makedirs(outdir)

for dem_fn in [dem_ref_fn]+dem_fn_list:
    print(dem_fn)
    dem_ds = iolib.fn_getds(dem_fn)
    dem = iolib.ds_getma(dem_ds)
    dem_hs_fn = os.path.splitext(dem_fn)[0]+'_hs_az315.tif'
    dem_hs = iolib.fn_getma(dem_hs_fn)
    anomaly = dem - dem_ref
    dt = timelib.fn_getdatetime(dem_fn)
    if dt is not None:
        title = dt.strftime('%Y-%m-%d')
    else: 
        title = 'Reference (2014-2016 DEM Mean)'
    f = makefig(dem, dem_hs, anomaly, ds=dem_ds, title=title)
    out_fn = os.path.join(outdir, os.path.splitext(dem_fn)[0]+'_fig.png')
    f.savefig(out_fn, bbox_inches='tight', dpi=150)

#~/src/demtools/make_movie.py rainier_allgood_mos-tile-0_warp_fig.png rainier_allgood_mos-tile-0_warp_fig.png rainier_allgood_mos-tile-0_warp_fig.png rainier_allgood_mos-tile-0_warp_fig.png rainier_allgood_mos-tile-0_warp_fig.png 2*png
#mv movie.mp4 rainier_summit_2014-2016_DEM_anomaly_movie_2panel.mp4

