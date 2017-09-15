#! /usr/bin/env python

"""
Create anomaly map time series
"""

import os
import sys
import glob

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid

from pygeotools.lib import iolib, timelib, geolib, malib
from imview.lib import pltlib

add_cbar = False

#dem_fn_list = glob.glob('*32m.tif')
#dem_ref_fn = 'rainier_allgood_mos-tile-0_warp.tif'
dem_fn_list = sys.argv[1:]

dems = np.ma.array([iolib.fn_getma(fn) for fn in dem_fn_list])
dem_clim = malib.calcperc(dems, (2,98))

w=10.0
h=7.5
f = plt.figure(figsize=(w,h))

n = len(dem_fn_list)
ncols = int(np.ceil(np.sqrt((float(n)*w)/h)))
nrows = int(np.ceil(float(n)/ncols))
#nrows = int(np.sqrt(n))+1
#ncols = nrows

if add_cbar:
    grid = ImageGrid(f, 111, nrows_ncols=(nrows, ncols), axes_pad=(0.05,0.2), share_all=True, cbar_mode='single', cbar_pad=0.1, cbar_set_cax=False)
else:
    grid = ImageGrid(f, 111, nrows_ncols=(nrows, ncols), axes_pad=(0.05,0.2), share_all=True, cbar_mode='None')

#f, axa = plt.subplots(nrows, ncols, sharex='all', sharey='all', figsize=(10,10))
#plt.subplots_adjust(wspace=0, hspace=0.2)
#for ax in axa.ravel():
#    pltlib.hide_ticks(ax)
#    ax.set_aspect('equal')

#Should extract this automatically
#Rainier
#dem_clim = (1300, 3700)
#SCG
#dem_clim = (760, 2270)
#Baker
#dem_clim = (550, 2650)
#Ngozumpa
#dem_clim = (4500, 7400)
#GM
#dem_clim = (1766, 3247)
#SBB
#dem_clim = (2934, 3983)
hs_clim = (1, 255)

for i,dem_fn in enumerate(dem_fn_list):
    ax = grid[i]
    print(dem_fn)
    dem_ds = iolib.fn_getds(dem_fn)
    dem = iolib.ds_getma(dem_ds)
    #dem_hs_fn = os.path.splitext(dem_fn)[0]+'_hs_az315.tif'
    #dem_hs = iolib.fn_getma(dem_hs_fn)
    dem_hs = geolib.gdaldem_mem_ds(dem_ds, 'hillshade', returnma=True)
    dt = timelib.fn_getdatetime(dem_fn)
    if dt is not None:
        title = dt.strftime('%Y-%m-%d')
        t = ax.set_title(title, fontdict={'fontsize':8})
        #t.set_position([0.5, 0.95])
    hs_im = ax.imshow(dem_hs, vmin=hs_clim[0], vmax=hs_clim[1], cmap='gray')
    dem_im = ax.imshow(dem, vmin=dem_clim[0], vmax=dem_clim[1], cmap='cpt_rainbow', alpha=0.5)
    ax.set_facecolor('k') 
    pltlib.hide_ticks(ax)

for ax in grid[i+1:]:
    ax.axis('off')

#for i in range(nrows*ncols):
#    ax = grid[i]

if add_cbar:
    cbar_lbl = 'Elevation (m WGS84)'
    cbar_kwargs = {'extend':'both', 'alpha':1.0}
    cbar = grid.cbar_axes[0].colorbar(dem_im, **cbar_kwargs) 
    cbar.update_bruteforce(dem_im)
    cbar.set_label_text(cbar_lbl)

#res = geolib.get_res(dem_ds)[0]
#pltlib.add_scalebar(grid[-1], res=res)

#out_fn = os.path.join(outdir, os.path.splitext(dem_fn)[0]+'_fig.png')
out_fn = os.path.join('dem_gallery.png')
f.savefig(out_fn, bbox_inches='tight', dpi=300)
