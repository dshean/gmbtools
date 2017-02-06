#! /usr/bin/env python

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap as Basemap

from pygeotools.lib import malib
from pygeotools.lib import geolib 
from imview.lib import pltlib

out_fn = 'conus_mb_20170204.csv'
out = np.loadtxt(out_fn, delimiter=',')

f, ax = plt.subplots(figsize=(10,10),dpi=300)
#Get this from glacier shp or DEM mosaic
extent = [-674693.945810, -729687.166879, 795021.159447, 688725.556370]
extent = geolib.pad_extent(geolib.extent_round(extent, precision=1000), width=60000)
w = extent[2] - extent[0]
h = extent[3] - extent[1]

conus_aea_proj = geolib.conus_aea_proj
lon_0, lat_0 = (-115, 43)
m = Basemap(width=w, height=h, resolution='h', area_thresh=10000, projection='aea',lat_1=36,lat_2=49,lon_0=-115,lat_0=43,ellps='WGS84')
xoff, yoff = m(lon_0, lat_0)
out[:,0] += xoff
out[:,1] += yoff
m.fillcontinents(color='lightyellow',zorder=0)
m.drawcoastlines(linewidth=0.25)
m.drawcountries(linewidth=0.25)
m.drawstates(linewidth=0.25)
m.drawrivers(color='lightblue',linewidth=0.25)
#m.shadedrelief()
m.drawmapscale(-115,37,lon_0,lat_0,400,barstyle='fancy')
parallels = np.arange(0.,91.,5.)
m.drawparallels(parallels,linewidth=0.5,labels=[True,False,False,False])
meridians = np.arange(-180.,181.,10.)
m.drawmeridians(meridians,linewidth=0.5,labels=[False,False,False,True])

out_mb = out[:,2]
out_stats = malib.print_stats(out_mb)
if out_stats[5] < 0:
    vmin = out_stats[5]-out_stats[6]*2
    vmax = -vmin 
else:
    vmax = out_stats[5]-out_stats[6]*2
    vmin = -vmax
#f, ax = plt.subplots(figsize=(10,10),dpi=300)
ax = plt.gca()
sc = m.scatter(out[:,0], out[:,1], c=out_mb, cmap='RdBu', s=out[:,3]*16, edgecolor='k', lw='0.2', vmin=vmin, vmax=vmax)
#ax.set(adjustable='box-forced', aspect='equal')
#ax.set_facecolor('0.5')
cbar = pltlib.add_cbar(ax, sc, label='Long-term mass balance (mwe/yr)')
plt.tight_layout()
plt.savefig('conus_mb.png', dpi=300)

out_dt = out[:,4]
vmin = out_dt.min()
vmax = out_dt.max()
f, ax = plt.subplots(figsize=(10,10),dpi=300)
sc = ax.scatter(out[:,0], out[:,1], c=out_dt, cmap='inferno', s=out[:,3]*16, vmin=vmin, vmax=vmax)
ax.set(adjustable='box-forced', aspect='equal')
ax.set_facecolor('0.5')
#cbar = pltlib.add_cbar(ax, sc, label='Time interval (yr)')
cbar = pltlib.add_cbar(ax, sc, label='NED Source Date (yr)')
plt.savefig('conus_dt.png', dpi=300)
