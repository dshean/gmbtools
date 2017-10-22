#! /usr/bin/env python

"""
Combine glacier polygons with geodetic mass balance output
Uses geopandas to join, compute stats for differen regions and plot
"""

import os, sys
import pandas as pd
import matplotlib.pyplot as plt

import cartopy.crs as ccrs 
import cartopy.feature as cfeature
#cartopy.crs.AlbersEqualArea(central_longitude=0.0, central_latitude=0.0, false_easting=0.0, false_northing=0.0, standard_parallels=(20.0, 50.0), globe=None)
crs = ccrs.AlbersEqualArea(central_longitude=-115, central_latitude=43, standard_parallels=(36, 49))

def cartopy_extent(extent):
    return [extent[0], extent[2], extent[1], extent[3]]

f = plt.figure()
ax = plt.axes(projection=crs)
#Doesn't work with current version of cartopy
#gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')

from osgeo import gdal
from pygeotools.lib import iolib, geolib
print("Loading image")
#hs_fn = 'conus_20171018_mos_32m_trans_100m_hs_az315.tif'
hs_fn = '/scr/mb/gpd/conus_20171018_mos_32m_trans_100m_hs_az315_1km.tif'
hs_ds = gdal.Open(hs_fn)
hs = iolib.ds_getma(hs_ds)
hs_extent = geolib.ds_extent(hs_ds)
hs_extent_cartopy = cartopy_extent(hs_extent)
print("Plotting image")
ax.imshow(hs, cmap='gray', origin='upper', extent=hs_extent_cartopy, transform=crs)

import geopandas as gpd

print("Loading glacier polygons")
glac_shp_fn = '/Users/dshean/data/rgi60/regions/rgi60_merge_CONUS_aea.shp'
glac_df = gpd.read_file(glac_shp_fn)

print("Loading mb")
mb_csv_fn = sys.argv[1]
mb_df = pd.read_csv(mb_csv_fn)
#mb_df['RGIId'] = 'RGI60-0'+mb_df['RGIId'].astype(str)
mb_df['RGIId'] = 'RGI60-'+mb_df['RGIId'].map('{:08.5f}'.format)
mergefield = 'RGIId'

"""
glac_shp_fn = 'conus_glacierpoly_24k_aea.shp'
mb_csv_fn = 'conus_mb_20171018_1434.csv'
glac_df = gpd.read_file(glac_shp_fn)
glac_df['GLACNUM'] = glac_df['GLACNUM'].astype(long)
mb_df = pd.read_csv(mb_csv_fn)
mb_df.rename(index=str, columns={"glacnum":"GLACNUM"}, inplace=True)
mergefield = 'GLACNUM'
"""

#mb_df['mb_mwea_area'] = mb_df['mb_mwea'] * mb_df['area_km2']*1E6

print("Merging")
glac_df_mb = glac_df.merge(mb_df, on=mergefield)
out_fn = os.path.splitext(mb_csv_fn)[0]+'_'+os.path.splitext(os.path.split(glac_shp_fn)[-1])[0]+'.shp'
glac_df_mb.to_file(out_fn)

#glac_df_mb_o = glac_df_mb[[mergefield,'mb_mwea']]

area_filter=True
if area_filter:
    #glac_df = glac_df_mb[glac_df_mb['area_km2'] > 0.1]
    glac_df = glac_df_mb[glac_df_mb['area_m2'] > 1E5]

print("Loading regions")
region_shp_fn = '/Users/dshean/Documents/UW/CONUS/regions/conus_mb_regions_aea.shp'
region_df = gpd.read_file(region_shp_fn)
region_df.crs = glac_df_mb.crs
region_df['region_centroid'] = region_df.centroid

print("Spatial join by region")
glac_df_mb_region = gpd.sjoin(glac_df_mb, region_df, how="inner", op="intersects")

#Watershed

print("Dissolve by region")
glac_df_mb_region_mean = glac_df_mb_region.dissolve(by='region', aggfunc='mean')
glac_df_mb_region_sum = glac_df_mb_region.dissolve(by='region', aggfunc='sum')
#glac_df_mb_region_sum = glac_df_mb_region.dissolve(by='region', aggfunc='sum').sort_values('area_km2')
#glac_df_mb_region_mean.reindex(glac_df_mb_region_sum.index)

#fig, ax = plt.subplots()
#ax.set_aspect('equal')

print("Loading states")
state_shp_fn = '/Users/dshean/Documents/UW/GIS_data/cb_2016_us_state_500k/cb_2016_us_state_500k_CONUS.shp'
state_df = gpd.read_file(state_shp_fn)

print("Plotting states")
state_df = state_df.to_crs(glac_df.crs)
style_kwd = {'facecolor':'none','edgecolor':'k', 'linewidth':0.5}
state_df.plot(ax=ax, **style_kwd)

#This is currently broken
#states = cfeature.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='50m', facecolor='none')
#ax.add_feature(states, edgecolor='k')

print("Plotting glacier polygons")
clim = (-1.0, 1.0)
style_kwd = {'edgecolor':'k', 'linewidth':0.5}
glac_df_mb.plot(ax=ax, column='mb_mwea', cmap='RdBu', vmin=clim[0], vmax=clim[1], **style_kwd)

#Plot scaled circles for each region at centroid x and y
#ax.scatter()

glac_df_extent = glac_df_mb.total_bounds
glac_df_extent_cartopy = cartopy_extent(glac_df_extent)
#Also currently broken
#ax.set_extent(glac_df_extent_cartopy, crs=crs)

f2, ax2 = plt.subplots(figsize=(8,8))
ax2.set_aspect('equal')
style_kwd = {'facecolor':'0.9','edgecolor':'k', 'linewidth':0.5}
state_df.plot(ax=ax2, **style_kwd)
#ax2.imshow(hs, cmap='gray', origin='upper', extent=hs_extent_cartopy)

scaling_f = 5
clim=(-0.4, 0.4)
x = glac_df_mb_region_mean['x']
y = glac_df_mb_region_mean['y']
s = scaling_f*glac_df_mb_region_sum['area_m2']/1E6
c = glac_df_mb_region_mean['mb_mwea']
#c = glac_df_mb_region_sum['mb_m3wea']/glac_df_mb_region_sum['area_m2']
sc = ax2.scatter(x, y, s, c, vmin=clim[0], vmax=clim[1], cmap='RdBu', edgecolor='k', linewidth=0.5)

from imview.lib import pltlib
pltlib.add_cbar(ax2, sc, label='Mass Balance (m we/yr)')
style_kwd = {'facecolor':'none', 'edgecolor':'k', 'linewidth':0.5}
glac_df_mb.plot(ax=ax2, **style_kwd)
pltlib.add_scalebar(ax2, res=1)
ax2.set_xlim(-856800, 910700)
ax2.set_ylim(-789000, 839400)
pltlib.hide_ticks(ax2)

for k, v in glac_df_mb_region_mean.iterrows():
    lbl = '%0.2f +/- %0.2f' % (v['mb_mwea'], v['mb_mwea_sigma'])
    ax2.annotate(lbl, xy=(v['x'],v['y']), xytext=(1,0), textcoords='offset points', family='sans-serif', fontsize=8, color='k')

plt.show()
