#! /usr/bin/env python

"""
Generate plots for CONUS glacier mass balance
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap as Basemap
from osgeo import osr

from pygeotools.lib import malib
from pygeotools.lib import geolib 
from imview.lib import pltlib

def get_equal_vmin_vmax(x):
    out_stats = malib.print_stats(x)
    if out_stats[5] < 0:
        vmin = out_stats[5]-out_stats[6]*2
        vmax = -vmin 
    else:
        vmax = out_stats[5]-out_stats[6]*2
        vmin = -vmax
    return vmin, vmax

def add_legend(ax, scale_factor=16, loc='upper right'):
    ax.autoscale(False)
    leg_s = np.array([0.1, 0.5, 1.0, 5.0, 10.0])
    leg_x = np.full(leg_s.size, -999999999)
    leg_y = np.full(leg_s.size, -999999999)
    #leg_sc = ax.scatter(leg_x, leg_y, c='0.8', s=leg_s)
    #ax.legend(leg_sc, ['%0.1f km^2' % s for s in leg_s], scatterpoints=1, loc='upper right')
    for i, s in enumerate(leg_s):
        lbl = r'$%0.1f\/km^2$' % s
        ax.scatter(leg_x[i], leg_y[i], s=s*scale_factor, c='gray', label=lbl)
    return ax.legend(title='Glacier Area', scatterpoints=1, loc=loc, prop={'size':8})

def mapplot(out, field):
    f, ax = plt.subplots(figsize=(8,8),dpi=300)
    #Get this from glacier shp or DEM mosaic
    extent = [-674693.945810, -729687.166879, 795021.159447, 688725.556370]
    extent = geolib.pad_extent(geolib.extent_round(extent, precision=1000), width=60000)
    w = extent[2] - extent[0]
    h = extent[3] - extent[1]

    conus_aea_proj = geolib.conus_aea_proj
    lon_0, lat_0 = (-115, 43)
    m = Basemap(width=w, height=h, resolution='h', area_thresh=10000, \
            projection='aea',lat_1=36,lat_2=49,lon_0=-115,lat_0=43,ellps='WGS84', ax=ax)
    xoff, yoff = m(lon_0, lat_0)
    x = out['x'] + xoff
    y = out['y'] + yoff
    m.fillcontinents(color='0.9',zorder=0)
    m.drawcoastlines(linewidth=0.25)
    m.drawcountries(linewidth=0.25)
    m.drawstates(linewidth=0.25)
    m.drawrivers(color='lightblue',linewidth=0.25,zorder=1)
    #m.shadedrelief()
    m.drawmapscale(-115,37,lon_0,lat_0,400,barstyle='fancy',zorder=10)
    parallels = np.arange(0.,91.,5.)
    m.drawparallels(parallels,linewidth=0.5,labels=[True,False,False,False])
    meridians = np.arange(-180.,181.,10.)
    m.drawmeridians(meridians,linewidth=0.5,labels=[False,False,False,True])
    if field == 'mb_mwea':
        cmap = 'RdBu'
        label = 'Mass balance (mwe/yr)'
        vmin,vmax = get_equal_vmin_vmax(out[field])
        lw = 0.2
    elif field == 't1':
        cmap = 'inferno'
        label = None 
        vmin = out[field].min()
        vmax = out[field].max()
        lw = 0.0
    sc = m.scatter(x, y, c=out[field], cmap=cmap, s=out['area_km2']*16, \
            edgecolor='k', lw=lw, vmin=vmin, vmax=vmax)
    #ax.set(adjustable='box-forced', aspect='equal')
    #ax.set_facecolor('0.5')
    cbar = pltlib.add_cbar(ax, sc, label=label)
    #plt.tight_layout()
    #plt.savefig('conus_mb.png', dpi=300)
    leg = add_legend(ax)
    return f

#out_fn = 'conus_mb_20170204.csv'
out_fn = 'conus_mb_summer2014-2016_20170205.csv'
#Load into structured array
out = np.genfromtxt(out_fn, delimiter=',', dtype=None, names=True)
#Sort by area
out = np.sort(out, order='area_km2')[::-1]

#Compute lat, lon
utm_srs = osr.SpatialReference()
utm_srs.ImportFromEPSG(32610)

lon, lat, dummy = geolib.cT_helper(out['x'],out['y'],0,geolib.conus_aea_srs,geolib.wgs_srs)
x_utm, y_utm, dummy = geolib.cT_helper(out['x'],out['y'],0,geolib.conus_aea_srs,utm_srs)
vmin, vmax = get_equal_vmin_vmax(out['mb_mwea'])

f = mapplot(out, field='mb_mwea')
fig_fn = 'conus_mb_map.png'
plt.title("Long-term (~30-60 year) Geodetic Mass Balance")
plt.savefig(fig_fn, dpi=300, bbox_inches='tight')

f = mapplot(out, field='t1')
fig_fn = 'conus_nedyear_map.png'
plt.title("NED Source Date")
plt.savefig(fig_fn, dpi=300, bbox_inches='tight')

f, ax = plt.subplots()
ax.set_facecolor('0.8')
sc = ax.scatter(lat, out['z_med'], c=out['mb_mwea'], s=out['area_km2']*16, \
        edgecolor='k', lw='0.2', cmap='RdBu', vmin=vmin, vmax=vmax)
leg = add_legend(ax, loc='lower left')
ax.set_xlabel('Latitude')
ax.set_ylabel('Elevation (m WGS84)')
ax.minorticks_on()
#ax.tick_params(left=True, right=True, bottom=True, top=True)
cbar = pltlib.add_cbar(ax, sc, label='Mass balance (mwe/yr)')
fig_fn = 'conus_mb_elev_lat.png'
plt.title("Long-term (~30-60 year) Geodetic Mass Balance")
plt.savefig(fig_fn, dpi=300, bbox_inches='tight')

plt.show()
