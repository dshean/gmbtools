#! /usr/bin/env python

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

from imview.lib import pltlib
from pygeotools.lib import warplib, geolib, iolib, malib, timelib

def plothist(ax, x,y,xlim,ylim):
    bins = (100, 100)
    H, xedges, yedges = np.histogram2d(x,y,range=[xlim,ylim],bins=bins)
    H = np.rot90(H)
    H = np.flipud(H)
    #Hmasked = np.ma.masked_where(H==0,H)
    Hmasked = H
    H_clim = malib.calcperc(Hmasked, (0,99))
    ax.pcolormesh(xedges,yedges,Hmasked,cmap='hot',vmin=H_clim[0], vmax=H_clim[1])

prism_fn = '/Users/dshean/data/PRISM_ppt_30yr_normal_800mM2_10-05_winter_cum.tif'
dem_fn = sys.argv[1]
#hs_fn = os.path.splitext(dem_fn)[0]+'_hs_az315.tif'
dz_fn = sys.argv[2]
dem_ts = timelib.fn_getdatetime(dem_fn)

hs = geolib.gdaldem_wrapper(dem_fn, product='hs')
slope = geolib.gdaldem_wrapper(dem_fn, product='slope')
aspect = geolib.gdaldem_wrapper(dem_fn, product='aspect')

dem_ds, dz_ds, prism_ds = warplib.memwarp_multi_fn([dem_fn, dz_fn, prism_fn], extent='first', res='first', t_srs='first')

res = geolib.get_res(dem_ds)[0]
dem = iolib.ds_getma(dem_ds)
dem_clim = malib.calcperc(dem, (1,99))
#dem_clim = (1200, 2950)

#hs = iolib.ds_getma(hs_ds)
hs_clim = (1,255)

rho_s = 0.5
dz = iolib.ds_getma(dz_ds)
swe = dz * rho_s
#swe_clim = malib.calcperc(swe, (1,99))
swe_clim = (0,6)

prism = iolib.ds_getma(prism_ds)/1000.
prism = np.ma.array(prism, mask=np.ma.getmaskarray(dz))

#Map plots
f, axa = plt.subplots(1, 3, figsize=(10,5), sharex=True, sharey=True, subplot_kw={'aspect':'equal', 'adjustable':'box-forced'})
hs_im = axa[0].imshow(hs, vmin=hs_clim[0], vmax=hs_clim[1], cmap='gray')
dem_im = axa[0].imshow(dem, vmin=dem_clim[0], vmax=dem_clim[1], cmap='cpt_rainbow', alpha=0.5)
axa[0].set_facecolor('k')
swe_im = axa[1].imshow(swe, vmin=swe_clim[0], vmax=swe_clim[1], cmap='inferno')
axa[1].set_facecolor('0.5')
prism_im = axa[2].imshow(prism, vmin=swe_clim[0], vmax=swe_clim[1], cmap='inferno')
axa[2].set_facecolor('0.5')
for ax in axa:
    pltlib.hide_ticks(ax)

axa[0].set_title('Summer %i' % dem_ts.year)
pltlib.add_cbar(axa[0], dem_im, label='Elevation (m WGS84)')
pltlib.add_scalebar(axa[0], res)
axa[1].set_title('WY%i' % dem_ts.year)
pltlib.add_cbar(axa[1], swe_im, label='SWE (m w.e.)')
axa[2].set_title('PRISM Oct-May Cumulative Precip')
pltlib.add_cbar(axa[2], swe_im, label='Precip (m w.e.)')
plt.tight_layout()

#Histogram plots
mask = malib.common_mask([dem, swe, slope, aspect])
dem = np.ma.array(dem, mask=mask).compressed()
#dem_clim = malib.calcperc(dem, (1,99))
swe = np.ma.array(swe, mask=mask).compressed()
#swe_clim = malib.calcperc(swe, (1,99))
prism = np.ma.array(prism, mask=mask).compressed()
slope = np.ma.array(slope, mask=mask).compressed()
slope_clim = malib.calcperc(slope, (0,99))
aspect = np.ma.array(aspect, mask=mask).compressed()
aspect_clim = (0., 360.)

f, axa = plt.subplots(1, 4, figsize=(10,5))

plothist(axa[0], swe, dem, swe_clim, dem_clim)
axa[0].set_xlabel('SWE (m w.e.)')
axa[0].set_ylabel('Elevation (m WGS84)')
plothist(axa[1], swe, slope, swe_clim, slope_clim)
axa[1].set_xlabel('SWE (m w.e.)')
axa[1].set_ylabel('Slope (deg)')
plothist(axa[2], swe, aspect, swe_clim, aspect_clim)
axa[2].set_xlabel('SWE (m w.e.)')
axa[2].set_ylabel('Aspect (deg)')
plothist(axa[3], swe, prism, swe_clim, swe_clim)
axa[3].set_xlabel('SWE (m w.e.)')
axa[3].set_ylabel('PRISM Precip (m w.e.)')
plt.tight_layout()

plt.show()
