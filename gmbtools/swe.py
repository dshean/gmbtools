#! /usr/bin/env python

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats

from imview.lib import pltlib
from pygeotools.lib import warplib, geolib, iolib, malib, timelib, filtlib

def plothist(ax, x,y,xlim,ylim):
    bins = (100, 100)
    H, xedges, yedges = np.histogram2d(x,y,range=[xlim,ylim],bins=bins)
    H = np.rot90(H)
    H = np.flipud(H)
    Hmasked = np.ma.masked_where(H==0,H)
    #Hmasked = H
    H_clim = malib.calcperc(Hmasked, (0,99))
    ax.pcolormesh(xedges,yedges,Hmasked,cmap='inferno',vmin=H_clim[0], vmax=H_clim[1])

def get_snotel_sites(dem_ds):
    snotel_fn = '/Users/dshean/src/conus/conus/snotel_latlon.csv'
    snotel_srs = geolib.wgs_srs
    sites = np.loadtxt(snotel_fn, delimiter=',', dtype=None)
    dem_extent = geolib.ds_extent(dem_ds, snotel_srs)
    valid_idx = ((sites[:,2] > dem_extent[0]) & (sites[:,2] < dem_extent[2]) & (sites[:,1] > dem_extent[1]) & (sites[:,1] < dem_extent[3]))
    valid_sites = sites[valid_idx]
    return valid_sites
    
def plot_snotel(ax, dem_ds):
    sites = get_snotel_sites(dem_ds)
    if sites.size > 0:
        mX, mY, dummy = geolib.cT_helper(sites[:,2], sites[:,1], 0, geolib.wgs_srs, geolib.get_ds_srs(dem_ds))
        pX, pY = geolib.mapToPixel(mX, mY, dem_ds.GetGeoTransform())
        ax.scatter(pX, pY, s=16, facecolors='w', edgecolors='k')
        for i, lbl in enumerate(sites[:,0]):
            bbox=dict(boxstyle='round,pad=0.1', fc='k', alpha=0.7)
            ax.annotate(str(int(lbl)), xy=(pX[i], pY[i]), xytext=(0, 4), textcoords='offset points', fontsize=8, color='w', bbox=bbox)

site='baker'

prism_fn = '/Users/dshean/data/PRISM_ppt_30yr_normal_800mM2_10-05_winter_cum.tif'
dem_fn = sys.argv[1]
#hs_fn = os.path.splitext(dem_fn)[0]+'_hs_az315.tif'
dz_fn = sys.argv[2]
dem_ts = timelib.fn_getdatetime(dem_fn)
wy = dem_ts.year + 1

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
swe_f = filtlib.rolling_fltr(swe, size=5)
swe_f = filtlib.gauss_fltr_astropy(swe, size=9)
swe = swe_f
swe_clim = np.array((0,6))

prism = iolib.ds_getma(prism_ds)/1000.
prism = np.ma.array(prism, mask=np.ma.getmaskarray(swe))

if True:
    #Map plots
    f, axa = plt.subplots(1, 3, figsize=(10,4), sharex=True, sharey=True, subplot_kw={'aspect':'equal', 'adjustable':'box-forced'})
    hs_im = axa[0].imshow(hs, vmin=hs_clim[0], vmax=hs_clim[1], cmap='gray')
    dem_im = axa[0].imshow(dem, vmin=dem_clim[0], vmax=dem_clim[1], cmap='cpt_rainbow', alpha=0.5)
    axa[0].set_facecolor('k')
    swe_im = axa[1].imshow(swe, vmin=swe_clim[0], vmax=swe_clim[1], cmap='inferno')
    axa[1].set_facecolor('0.3')
    prism_im = axa[2].imshow(prism, vmin=swe_clim[0], vmax=swe_clim[1], cmap='inferno')
    axa[2].set_facecolor('0.3')
    for ax in axa:
        pltlib.hide_ticks(ax)

    axa[0].set_title('Late Summer %i' % dem_ts.year, fontdict={'fontsize':8})
    pltlib.add_cbar(axa[0], dem_im, label='Elevation (m WGS84)')
    pltlib.add_scalebar(axa[0], res)
    axa[1].set_title('WY%i (Summer %i to Spring %i Elev. Diff.)' % (wy, dem_ts.year, (dem_ts.year+1)), fontdict={'fontsize':8})
    pltlib.add_cbar(axa[1], swe_im, label=r'SWE Estimate (m w.e., $\rho_s$=0.5)')
    axa[2].set_title('~30-year PRISM Normal: Oct-May Precip', fontdict={'fontsize':8})
    pltlib.add_cbar(axa[2], swe_im, label='Cumulative Precip (m w.e.)')

    #plot_snotel(axa[2], dem_ds)
    plot_snotel(axa[0], dem_ds)

    plt.tight_layout()
    fig_fn = '%s_WY%i_SWE_maps.png' % (site, wy)
    plt.savefig(fig_fn, dpi=300, bbox_inches='tight')

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

swe_dem_slope, swe_dem_intercept, r_value, p_value, std_err = scipy.stats.linregress(swe, dem)
swe_dem_f = swe_dem_slope * swe_clim + swe_dem_intercept

if True:
    f, axa = plt.subplots(1, 4, figsize=(10,2.5))
    dem_clim = malib.calcperc(dem, (5,99.9))
    plothist(axa[0], dem, swe, dem_clim, swe_clim)
    axa[0].set_ylabel('SWE (m w.e.)')
    axa[0].set_xlabel('Elevation (m WGS84)')
    axa[0].plot(swe_dem_f, swe_clim, color='limegreen', ls='--', lw=0.5)
    plothist(axa[1], slope, swe, slope_clim, swe_clim)
    axa[1].set_ylabel('SWE (m w.e.)')
    axa[1].set_xlabel('Slope (deg)')
    plothist(axa[2], aspect, swe, aspect_clim, swe_clim)
    axa[2].set_ylabel('SWE (m w.e.)')
    axa[2].set_xlabel('Aspect (deg)')
    plothist(axa[3], prism, swe, swe_clim, swe_clim)
    axa[3].set_ylabel('SWE (m w.e.)')
    axa[3].set_xlabel('PRISM Precip (m w.e.)')
    axa[3].plot(swe_clim, swe_clim, color='limegreen', ls='-', lw=0.5)
    for ax in axa:
        ax.set_facecolor('0.3')
    plt.tight_layout()
    fig_fn = '%s_WY%i_SWE_WV_plots.png' % (site, wy)
    plt.savefig(fig_fn, dpi=300, bbox_inches='tight')

if True:
    swe = prism
    f, axa = plt.subplots(1, 4, figsize=(10,2.5), facecolor='w')
    dem_clim = malib.calcperc(dem, (5,99.9))
    plothist(axa[0], dem, swe, dem_clim, swe_clim)
    axa[0].set_ylabel('PRISM Precip (m w.e.)')
    axa[0].set_xlabel('Elevation (m WGS84)')
    axa[0].plot(swe_dem_f, swe_clim, color='limegreen', ls='--', lw=0.5)
    plothist(axa[1], slope, swe, slope_clim, swe_clim)
    axa[1].set_ylabel('PRISM Precip (m w.e.)')
    axa[1].set_xlabel('Slope (deg)')
    plothist(axa[2], aspect, swe, aspect_clim, swe_clim)
    axa[2].set_ylabel('PRISM Precip (m w.e.)')
    axa[2].set_xlabel('Aspect (deg)')
    plothist(axa[3], prism, swe, swe_clim, swe_clim)
    axa[3].set_ylabel('PRISM Precip (m w.e.)')
    axa[3].set_xlabel('PRISM Precip (m w.e.)')
    axa[3].plot(swe_clim, swe_clim, color='limegreen', ls='-', lw=0.5)
    for ax in axa:
        ax.set_facecolor('0.3')
    plt.tight_layout()
    fig_fn = '%s_WY%i_SWE_PRISM_plots.png' % (site, wy)
    plt.savefig(fig_fn, dpi=300, bbox_inches='tight')

plt.show()
