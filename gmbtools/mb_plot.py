#! /usr/bin/env python

"""
Generate plots for regional glacier mass balance
"""

import os
import sys
from datetime import datetime

import numpy as np
from numpy.lib import recfunctions
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap as Basemap
from mpl_toolkits.mplot3d import Axes3D
from osgeo import osr

from pygeotools.lib import malib
from pygeotools.lib import geolib 
from imview.lib import pltlib

glacier_dict = {}
#glacier_dict[10480] = 'SouthCascadeGlacier'
#glacier_dict[(-447928.98,614122.05)] = 'SCG'
#glacier_dict[10480] = 'SCG'

def get_equal_vmin_vmax(x):
    a_stats = malib.print_stats(x)
    if a_stats[5] < 0:
        vmin = a_stats[5]-a_stats[6]*2
        vmax = -vmin 
    else:
        vmax = a_stats[5]+a_stats[6]*2
        vmin = -vmax
    return vmin, vmax

def add_legend(ax, sf=16, loc='upper right'):
    """
    Create legend for scaled scatterplot markers
    """
    ax.autoscale(False)
    #leg_s = np.array([0.1, 0.5, 1.0, 5.0, 10.0])
    leg_s = np.array([0.1, 1.0, 10.0, 100.0])
    leg_x = np.full(leg_s.size, -999999999)
    leg_y = np.full(leg_s.size, -999999999)
    #leg_sc = ax.scatter(leg_x, leg_y, c='0.8', s=leg_s)
    #ax.legend(leg_sc, ['%0.1f km^2' % s for s in leg_s], scatterpoints=1, loc='upper right')
    for i, s in enumerate(leg_s):
        lbl = r'$%0.1f\/km^2$' % s
        ax.scatter(leg_x[i], leg_y[i], s=s*sf, c='gray', label=lbl)
    return ax.legend(title='Glacier Area', scatterpoints=1, loc=loc, prop={'size':8})

def add_lbl(ax, a):
    for pt,lbl in glacier_dict.iteritems():
        ax.annotate(lbl,xy=pt)

def scplot_fields(a, fx, fy, fc, fs, sf=16, ax=None, clim=None):
    ax = scplot(a[fx], a[fy], a[fc], a[fs], sf, ax, clim)
    for glacnum,lbl in glacier_dict.iteritems():
        b = a[a['glacnum'] == float(glacnum)]
        ax.annotate(lbl, xy=(b[fx], b[fy]), fontsize=6)
    return ax

def scplot(x, y, c, s, sf=16, ax=None, clim=None):
    """
    Create scatter plot with scaled circles
    """
    if clim is None:
        vmin, vmax = get_equal_vmin_vmax(c)
    else:
        vmin, vmax = clim
    if ax is None:
        ax = plt.gca()
    ax.set_facecolor('0.8')
    sc = ax.scatter(x, y, c=c, s=s*sf, edgecolor='k', lw='0.2', cmap='RdBu', vmin=vmin, vmax=vmax)
    cbar = pltlib.add_cbar(ax, sc, label='Mass balance (m we/yr)')
    #leg = add_legend(ax, sf=sf, loc='lower left')
    #leg = add_legend(ax, sf=sf, loc='upper right')
    leg = add_legend(ax, sf=sf, loc='lower right')
    ax.minorticks_on()
    #ax.tick_params(left=True, right=True, bottom=True, top=True)
    return ax

def scplot3D(x, y, z, c, s, sf=4, ax=None, clim=None):
    """
    Create scatter plot with scaled circles
    """
    if clim is None:
        vmin, vmax = get_equal_vmin_vmax(c)
    else:
        vmin, vmax = clim
    ax.set_facecolor('0.8')
    kwargs = {'edgecolor':'k', 'lw':'0.2', 'cmap':'RdBu', 'vmin':vmin, 'vmax':vmax}
    #kwargs = {'cmap':'RdBu', 'vmin':vmin, 'vmax':vmax}
    sc = ax.scatter(x, y, z, c=c, s=s*sf, depthshade=True, **kwargs)
    cbar = pltlib.add_cbar(ax, sc, label='Mass balance (m we/yr)')
    leg = add_legend(ax, sf=sf, loc='lower left')
    ax.minorticks_on()
    #ax.tick_params(left=True, right=True, bottom=True, top=True)
    return ax

def mapplot(a, field, srs, sf=16, ax=None):
    if ax is None:
        f, ax = plt.subplots()
    #Get this from glacier shp or DEM mosaic
    #CONUS
    #extent = [-674693.945810, -729687.166879, 795021.159447, 688725.556370]
    #extent = geolib.pad_extent(geolib.extent_round(extent, precision=1000), width=60000)
    extent = [a['x'].min(), a['y'].min(), a['x'].max(), a['y'].max()]
    print(extent)
    extent = geolib.pad_extent(geolib.extent_round(extent, precision=1000), width=100000)
    print(extent)
    #w = extent[2] - extent[0]
    #h = extent[3] - extent[1]
    w = 2*np.abs([extent[0], extent[2]]).max()
    h = 2*np.abs([extent[1], extent[3]]).max()
    lon_0, lat_0 = (srs.GetProjParm("longitude_of_center"), srs.GetProjParm("latitude_of_center"))
    lat_1, lat_2 = (srs.GetProjParm("standard_parallel_1"), srs.GetProjParm("standard_parallel_2"))
    proj_kwargs = {'projection':'aea','lat_1':lat_1,'lat_2':lat_2,'lon_0':lon_0,'lat_0':lat_0,'ellps':'WGS84'}
    print(proj_kwargs)
    m = Basemap(width=w, height=h, resolution='h', area_thresh=10000, ax=ax, **proj_kwargs)
    #clon, clat, dummy = geolib.cT_helper(extent[0]+w/2.0,extent[1]+h/2.0,0,srs,geolib.wgs_srs)
    xoff, yoff = m(lon_0, lat_0)
    print(xoff, yoff)
    x = a['x'] + xoff
    y = a['y'] + yoff
    m.fillcontinents(color='0.9',zorder=0)
    m.drawcoastlines(linewidth=0.25,zorder=1)
    m.drawcountries(linewidth=0.25,zorder=1)
    m.drawstates(linewidth=0.25,zorder=1)
    m.drawrivers(color='lightblue',linewidth=0.25,zorder=2)
    #m.shadedrelief()
    #m.bluemarble()
    #m.drawmapscale(-115,37,lon_0,lat_0,400,barstyle='fancy',zorder=10)
    parallels = np.arange(0.,91.,5.)
    m.drawparallels(parallels,linewidth=0.5,labels=[True,False,False,False],zorder=2)
    meridians = np.arange(-180.,181.,10.)
    m.drawmeridians(meridians,linewidth=0.5,labels=[False,False,False,True],zorder=2)
    if field == 'mb_mwea':
        cmap = 'RdBu'
        label = 'Mass balance (m we/yr)'
        vmin,vmax = get_equal_vmin_vmax(a[field])
        lw = 0.1
    elif field == 't1':
        cmap = 'inferno'
        label = None 
        vmin = a[field].min()
        vmax = a[field].max()
        lw = 0.0
    print sf
    sc = m.scatter(x, y, c=a[field], cmap=cmap, s=a['area_km2']*sf, \
            edgecolor='k', lw=lw, vmin=vmin, vmax=vmax)
    #ax.set(adjustable='box-forced', aspect='equal')
    #ax.set_facecolor('0.5')
    cbar = pltlib.add_cbar(ax, sc, label=label)
    #plt.tight_laya()
    #plt.savefig('conus_mb.png', dpi=300)
    leg = add_legend(ax, sf=sf, loc='upper right')
    return ax

csv_fn = sys.argv[1]

if 'conus' in csv_fn or 'other' in csv_fn:
    site = 'conus'
elif 'hma' in csv_fn:
    site = 'hma'

#Load into structured array
a = np.genfromtxt(csv_fn, delimiter=',', dtype=None, names=True)
#Sort by area
a = np.sort(a, order='area_km2')[::-1]

ts = datetime.now().strftime('%Y%m%d_%H%M')

if site == 'conus':
    aea_srs = geolib.conus_aea_srs
    title = "CONUS Long-term Glacier Mass Balance (~1950s-1980s to 2015)"
    sf = 16
elif site == 'hma':
    aea_srs = geolib.hma_aea_srs
    title = "HMA ~15-year Glacier Mass Balance (2000 to 2015)"
    sf = 2 

#Compute lat, lon
lon, lat, dummy = geolib.cT_helper(a['x'],a['y'],0,aea_srs,geolib.wgs_srs)
a = recfunctions.append_fields(a, 'lon', lon, dtypes=None, usemask=False)
a = recfunctions.append_fields(a, 'lat', lat, dtypes=None, usemask=False)

#utm_srs = osr.SpatialReference()
#utm_srs.ImportFromEPSG(32610)
#x_utm, y_utm, dummy = geolib.cT_helper(a['x'],a['y'],0,aea_srs,utm_srs)

vmin, vmax = get_equal_vmin_vmax(a['mb_mwea'])

if True:
    #f, ax = plt.subplots(figsize=(10,8))
    f, ax = plt.subplots()
    ax = mapplot(a, field='mb_mwea', srs=aea_srs, sf=sf, ax=ax)
    fig_fn = '%s_mb_map_%s.png' % (site, ts)
    ax.set_title(title)
    plt.savefig(fig_fn, dpi=300, bbox_inches='tight')

if False:
    #f, ax = plt.subplots(figsize=(10,8))
    f, ax = plt.subplots()
    ax = mapplot(a, field='t1', srs=aea_srs, sf=sf, ax=ax)
    fig_fn = '%s_nedyear_map_%s.png' % (site, ts)
    plt.title("NED Source Date")
    plt.savefig(fig_fn, dpi=300, bbox_inches='tight')

if True:
    f, ax = plt.subplots()
    #ax = scplot(lat, a['z_med'], a['mb_mwea'], a['area_km2'], sf=sf, ax=ax, clim=(vmin, vmax))
    ax = scplot_fields(a, 'lat', 'z_med', 'mb_mwea', 'area_km2', sf=sf, ax=ax, clim=(vmin, vmax))
    ax.set_title(title)
    ax.set_xlabel('Latitude')
    ax.set_ylabel('Elevation (m WGS84)')
    fig_fn = '%s_mb_elev_lat_%s.png' % (site, ts)
    plt.savefig(fig_fn, dpi=300, bbox_inches='tight')

if True:
    f, ax = plt.subplots()
    #ax = scplot(lon, a['z_med'], a['mb_mwea'], a['area_km2'], sf=sf, ax=ax, clim=(vmin, vmax))
    ax = scplot_fields(a, 'lon', 'z_med', 'mb_mwea', 'area_km2', sf=sf, ax=ax, clim=(vmin, vmax))
    ax.set_title(title)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Elevation (m WGS84)')
    fig_fn = '%s_mb_elev_lon_%s.png' % (site, ts)
    plt.savefig(fig_fn, dpi=300, bbox_inches='tight')

#3D plot - pretty slow
if False:
    f = plt.figure()
    ax = f.add_subplot(111, projection='3d')
    #ax = scplot3D(lon, lat, a['z_med'], a['mb_mwea'], a['area_km2']*4, ax=ax, clim=(vmin, vmax))
    ax = scplot3D(a['x'], a['y'], a['z_med'], a['mb_mwea'], a['area_km2'], sf=sf, ax=ax, clim=(vmin, vmax))
    ax.set_title(title)
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_zlabel('Elevation (m WGS84)')
    #fig_fn = '%s_mb_elev_lat_%s.png' % (site, ts)
    #plt.savefig(fig_fn, dpi=300, bbox_inches='tight')

if 'ppt_a' in a.dtype.names:
    if True:
        f, ax = plt.subplots()
        #ax = scplot(a['ppt_a'], a['tmean_a'], a['mb_mwea'], a['area_km2'], sf=sf, ax=ax, clim=(vmin, vmax))
        ax = scplot_fields(a, 'ppt_a', 'tmean_a', 'mb_mwea', 'area_km2', sf=sf, ax=ax, clim=(vmin, vmax))
        ax.set_title(title)
        ax.set_xlabel('Mean Annual Precip (m we)')
        ax.set_ylabel('Mean Annual Temp (C)')
        ax.set_ylim(-6,6)
        ax.axhline(0, ls=':', color='k', lw=0.5)
        fig_fn = '%s_mb_ppt_a_tmean_a_%s.png' % (site, ts)
        plt.savefig(fig_fn, dpi=300, bbox_inches='tight')

    if True:
        f, ax = plt.subplots()
        #ax = scplot(a['ppt_w'], a['tmean_w'], a['mb_mwea'], a['area_km2'], sf=sf, ax=ax, clim=(vmin, vmax))
        ax = scplot_fields(a, 'ppt_w', 'tmean_a', 'mb_mwea', 'area_km2', sf=sf, ax=ax, clim=(vmin, vmax))
        ax.set_title(title)
        ax.set_xlabel('Mean Winter Precip (m we)')
        ax.set_ylabel('Mean Annual Temp (C)')
        ax.set_ylim(-6, 6)
        ax.axhline(0, ls=':', color='k', lw=0.5)
        fig_fn = '%s_mb_ppt_w_tmean_a_%s.png' % (site, ts)
        plt.savefig(fig_fn, dpi=300, bbox_inches='tight')

    if True:
        f, ax = plt.subplots()
        #ax = scplot(a['ppt_w'], a['tmean_w'], a['mb_mwea'], a['area_km2'], sf=sf, ax=ax, clim=(vmin, vmax))
        ax = scplot_fields(a, 'ppt_w', 'tmean_w', 'mb_mwea', 'area_km2', sf=sf, ax=ax, clim=(vmin, vmax))
        ax.set_title(title)
        ax.set_xlabel('Mean Winter Precip (m we)')
        ax.set_ylabel('Mean Winter Temp (C)')
        ax.set_ylim(-9,3)
        ax.axhline(0, ls=':', color='k', lw=0.5)
        fig_fn = '%s_mb_ppt_w_tmean_w_%s.png' % (site, ts)
        plt.savefig(fig_fn, dpi=300, bbox_inches='tight')

    if True:
        f, ax = plt.subplots()
        #ax = scplot(a['ppt_w'], a['tmean_s'], a['mb_mwea'], a['area_km2'], sf=sf, ax=ax, clim=(vmin, vmax))
        ax = scplot_fields(a, 'ppt_w', 'tmean_s', 'mb_mwea', 'area_km2', sf=sf, ax=ax, clim=(vmin, vmax))
        ax.set_title(title)
        ax.set_xlabel('Mean Winter Precip (m we)')
        ax.set_ylabel('Mean Summer Temp (C)')
        ax.set_ylim(1,13)
        ax.axhline(0, ls=':', color='k', lw=0.5)
        fig_fn = '%s_mb_ppt_w_tmean_s_%s.png' % (site, ts)
        plt.savefig(fig_fn, dpi=300, bbox_inches='tight')

plt.show()
