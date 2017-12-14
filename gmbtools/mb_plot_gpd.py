#! /usr/bin/env python

"""
Combine glacier polygons with geodetic mass balance output
Uses geopandas to join, compute stats for differen regions and plot
"""

import os, sys
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from imview.lib import pltlib

from osgeo import gdal
from pygeotools.lib import iolib, geolib

import geopandas as gpd

def cartopy_extent(extent):
    return [extent[0], extent[2], extent[1], extent[3]]

def append_centroid_xy(df):
    xlist = []
    ylist = []
    for index, row in df.centroid.iteritems():
        xlist.append(row.x)
        ylist.append(row.y)
    df['centroid_x'] = xlist
    df['centroid_y'] = ylist

site = 'hma'
#site = 'conus'

area_filter=True
#min_area_m2 = 1E5
min_area_m2 = 1E6

mb_clim = (-1.0, 1.0)

extent = None
if site == 'hma':
    hs_fn = '/Users/dshean/Documents/UW/HMA/mos/hma_20170716_mos/mos_32m/hma_20170716_mos_32m_100m_hs_az315.tif'
    #Subset of Himalayas for testing
    #glac_shp_fn = '/Users/dshean/data/rgi60/regions/rgi60_merge_HMA_aea_test.shp'
    glac_shp_fn = '/Users/dshean/data/rgi60/regions/rgi60_merge_HMA_aea.shp'
    border_shp_fn = '/Users/dshean/data/NaturalEarth/10m_cultural/10m_cultural/ne_10m_admin_0_countries_lakes.shp'
    basin_shp_fn = '/Users/dshean/data/HydroBASINS/hybas_lake_as_lev01-12_v1c/hybas_lake_as_lev04_v1c.shp'
    basin_col = 'HYBAS_ID'
    #region_shp_fn = None
    #region_shp_fn = basin_shp_fn
    #region_col = 'HYBAS_ID'
    #region_shp_fn = '/Users/dshean/data/geocells/geocell_1deg_HMA_rgi60_intersect.shp'
    #region_col = 'name'
    region_shp_fn = '/Users/dshean/Documents/UW/HMA/Kaab_regions/regions_from_kaab2015_merged.shp'
    #Already have name in glac_shp_df
    region_col = 'Name'
    glac_crs = {u'datum':u'WGS84',u'lat_0':36,u'lat_1':25,u'lat_2':47,u'lon_0':85,u'no_defs':True,u'proj':u'aea',u'units':u'm',u'x_0':0,u'y_0':0}
    extent = [-1610900, -1142400, 1767400, 1145700]
elif site == 'conus':
    #hs_fn = '/scr/mb/gpd/conus_20171018_mos_32m_trans_100m_hs_az315_1km.tif'
    glac_shp_fn = '/Users/dshean/data/rgi60/regions/rgi60_merge_CONUS_aea.shp'
    region_shp_fn = '/Users/dshean/Documents/UW/CONUS/regions/conus_mb_regions_aea.shp'
    region_col = 'region'
    #border_shp_fn = '/Users/dshean/Documents/UW/GIS_data/cb_2016_us_state_500k/cb_2016_us_state_500k_CONUS.shp'
    border_shp_fn = '/Users/dshean/data/NaturalEarth/10m_cultural/10m_cultural/ne_10m_admin_1_states_provinces_lakes.shp'
    basin_shp_fn = '/Users/dshean/data/HydroBASINS/hybas_lake_na_lev01-12_v1c/hybas_lake_na_lev04_v1c.shp'
    #minx, miny, maxx, maxy
    extent = [-856800, -789000, 910700, 839400]
else:
    sys.exit("Site not currently supported")

#Input csv from mb_parallel.py
mb_csv_fn = sys.argv[1]

merge_fn = os.path.splitext(mb_csv_fn)[0]+'_'+os.path.splitext(os.path.split(glac_shp_fn)[-1])[0]+'.geojson'
#driver = 'ESRI Shapefile'
driver = 'GeoJSON'

print("Loading glacier polygons")
glac_shp_join_fn = os.path.splitext(merge_fn)[0]+'_join.geojson'
if os.path.exists(glac_shp_join_fn):
    glac_df = gpd.read_file(glac_shp_join_fn)
    #This is a hack, as geojson doesn't properly preserve custom aae proj
    glac_df.crs = glac_crs
else:
    glac_df = gpd.read_file(glac_shp_fn)

if region_shp_fn is not None:
    print("Loading regions")
    region_df = gpd.read_file(region_shp_fn)
    #Convert to glac crs
    region_df = region_df.to_crs(glac_df.crs)

if basin_shp_fn is not None:
    print("Loading basins")
    basin_df = gpd.read_file(basin_shp_fn)
    #Convert to glac crs
    basin_df = basin_df.to_crs(glac_df.crs)

#Add region and basin fields to RGI polygons
if not os.path.exists(glac_shp_join_fn):
    if region_shp_fn is not None:
        print("One-time spatial join by region")
        glac_df = gpd.sjoin(glac_df, region_df, how="inner", op="intersects")
        glac_df.rename(index=str, columns={region_col+"_right": "region", 'index_right':'region_id'}, \
                inplace=True)

    if basin_shp_fn is not None:
        print("One-time spatial join by basin")
        glac_df = gpd.sjoin(glac_df, basin_df, how="inner", op="intersects")
        glac_df.rename(index=str, columns={"HYBAS_ID": "basin"}, inplace=True)

    print("Writing out: %s" % glac_shp_join_fn)
    glac_df.to_file(glac_shp_join_fn, driver=driver)

region_col = 'region'
basin_col = 'basin'

print("Loading mb")
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

if os.path.exists(merge_fn):
    print("Loading merged polygons and mb")
    glac_df_mb = gpd.read_file(merge_fn)
else:
    print("Merging glacier polygons and mb results")
    glac_df_mb = glac_df.merge(mb_df, on=mergefield)
    glac_df_mb.to_file(merge_fn, driver=driver)

if area_filter:
    print("Filtering by glacier polygon area (min %0.2f km^2)" % (min_area_m2/1E6))
    orig_count = glac_df_mb.shape[0]
    glac_df_mb = glac_df_mb[glac_df_mb['area_m2'] > min_area_m2]
    print("%i of %i records preserved" % (glac_df_mb.shape[0], orig_count))

#glac_df_mb_o = glac_df_mb[[mergefield,'mb_mwea']]

if region_shp_fn is not None:
    glac_df_mb_region_fn = os.path.splitext(merge_fn)[0]+'_region_dissolve.geojson'
    if os.path.exists(glac_df_mb_region_fn):
        print("Loading region dissolve")
        glac_df_mb_region = gpd.read_file(glac_df_mb_region_fn)
    else:
        print("Dissolve by region")
        aggfunc = {'area_m2':[np.mean, np.sum], 'mb_mwea':[np.mean, np.std, np.sum], 'mb_mwea_sigma':np.mean}
        #glac_df_mb_region = glac_df_mb.dissolve(by=region_col+'_right', aggfunc=aggfunc)
        glac_df_mb_region = glac_df_mb.dissolve(by=region_col, aggfunc=aggfunc)
        append_centroid_xy(glac_df_mb_region)
        #glac_df_mb_region_sum = glac_df_mb_region.dissolve(by=region_col, aggfunc='sum').sort_values('area_km2')
        #glac_df_mb_region_mean.reindex(glac_df_mb_region_sum.index)
        #No need for tuple column names
        glac_df_mb_region.columns = ['_'.join(x) if isinstance(x, tuple) else x \
                for x in glac_df_mb_region.columns]
        glac_df_mb_region.to_file(glac_df_mb_region_fn, driver=driver)

if basin_shp_fn is not None:
    glac_df_mb_basin_fn = os.path.splitext(merge_fn)[0]+'_basin_dissolve.geojson'
    if os.path.exists(glac_df_mb_basin_fn):
        print("Loading basin dissolve")
        glac_df_mb_basin = gpd.read_file(glac_df_mb_basin_fn)
    else:
        print("Dissolve by basin")
        aggfunc = {'area_m2':[np.mean, np.sum], 'mb_mwea':[np.mean, np.std, np.sum], 'mb_mwea_sigma':np.mean}
        glac_df_mb_basin = glac_df_mb.dissolve(by=basin_col, aggfunc=aggfunc)
        append_centroid_xy(glac_df_mb_basin)
        glac_df_mb_basin.columns = ['_'.join(x) if isinstance(x, tuple) else x \
                for x in glac_df_mb_basin.columns]
        glac_df_mb_basin.to_file(glac_df_mb_basin_fn, driver=driver)

if True:
    import cartopy.crs as ccrs 
    #cartopy.crs.AlbersEqualArea(central_longitude=0.0, central_latitude=0.0, false_easting=0.0, false_northing=0.0, standard_parallels=(20.0, 50.0), globe=None)
    crs = ccrs.AlbersEqualArea(central_longitude=-115, central_latitude=43, standard_parallels=(36, 49))
    fig = plt.figure()
    ax = plt.axes(projection=crs)
    #Currently unsupported for AEA
    #gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')

#fig, ax = plt.subplots()
#ax.set_aspect('equal')

if False:
    print("Loading shaded relief map")
    hs_ds = gdal.Open(hs_fn)
    hs = iolib.ds_getma(hs_ds)
    hs_extent = geolib.ds_extent(hs_ds)
    hs_extent_cartopy = cartopy_extent(hs_extent)
    print("Plotting image")
    ax.imshow(hs, cmap='gray', origin='upper', extent=hs_extent_cartopy, transform=crs, alpha=0.6)

"""
#This is currently broken
import cartopy.feature as cfeature
borders = cfeature.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='50m', facecolor='none')
ax.add_feature(borders, edgecolor='k')
"""

if border_shp_fn is not None:
    print("Loading borders")
    #Load local copy of border polygons
    border_df = gpd.read_file(border_shp_fn)
    border_df = border_df.to_crs(glac_df.crs)
    print("Plotting borders")
    border_style = {'facecolor':'0.9','edgecolor':'k', 'linewidth':0.5}
    border_df.plot(ax=ax, **border_style)

if basin_shp_fn is not None:
    basin_style = {'facecolor':'none','edgecolor':'b', 'linewidth':0.5, 'alpha':0.5}
    basin_df.plot(ax=ax, **basin_style)

print("Plotting glacier polygons")
glac_style = {'edgecolor':'k', 'linewidth':0.5}
#glac_df_mb.plot(ax=ax, column='mb_mwea', cmap='RdBu', vmin=mb_clim[0], vmax=mb_clim[1], **glac_style)

#This is minx, miny, maxx, maxy
if extent is None:
    glac_df_extent = glac_df_mb.total_bounds
    extent = glac_df_extent

glac_df_extent_cartopy = cartopy_extent(extent)
mb_sc = glac_df_mb.plot(ax=ax, column='mb_mwea', cmap='RdBu', vmin=mb_clim[0], vmax=mb_clim[1], **glac_style)

ax.set_extent(glac_df_extent_cartopy, crs=crs)
pltlib.add_scalebar(ax, res=1)

"""
f3, (t1_ax, mb_ax) = plt.subplots(1,2,figsize=(12,8), sharex=True, sharey=True)
t1_clim = (glac_df_mb['t1'].min(), glac_df_mb['t1'].max()) 
for ax in (t1_ax, mb_ax):
    ax.set_aspect('equal')
    border_df.plot(ax=ax, **border_style)
    pltlib.hide_ticks(ax)
t1_sc = glac_df_mb.plot(ax=t1_ax, column='t1', cmap='inferno', vmin=t1_clim[0], vmax=t1_clim[1], **glac_style)
mb_sc = glac_df_mb.plot(ax=mb_ax, column='mb_mwea', cmap='RdBu', vmin=mb_clim[0], vmax=mb_clim[1], **glac_style)
#This is a hack
#https://stackoverflow.com/questions/36008648/colorbar-on-geopandas
t1_sm = plt.cm.ScalarMappable(cmap='inferno', norm=matplotlib.colors.Normalize(vmin=t1_clim[0], vmax=t1_clim[1]))
t1_sm._A = []
pltlib.add_cbar(t1_ax, t1_sm, label='Source Date (yr)')
mb_sm = plt.cm.ScalarMappable(cmap='RdBu', norm=matplotlib.colors.Normalize(vmin=mb_clim[0], vmax=mb_clim[1]))
mb_sm._A = []
pltlib.add_cbar(mb_ax, mb_sm, label='Mass Balance (m we/yr)')
#t1_ax.set_xlim(extent[0], extent[2])
#t1_ax.set_ylim(extent[1], extent[3])
"""

if region_shp_fn is not None:
    f2, ax2 = plt.subplots(figsize=(8,8))
    ax2.set_aspect('equal')
    border_df.plot(ax=ax2, **border_style)
    #ax2.imshow(hs, cmap='gray', origin='upper', extent=hs_extent_cartopy)

    #This is original region_col
    region_style = {'column':'Name', 'cmap':'cpt_rainbow', 'edgecolor':'k', 'linewidth':0.5, 'alpha':0.2}
    #region_style = {'column':'Name', 'cmap':'gray', 'edgecolor':'k', 'linewidth':0.5, 'alpha':0.3}
    region_df.plot(ax=ax2, **region_style)

    scaling_f = 1
    #Update these with centroids
    x = glac_df_mb_region['centroid_x']
    y = glac_df_mb_region['centroid_y']
    s = scaling_f*glac_df_mb_region[('area_m2_sum')]/1E6
    c = glac_df_mb_region[('mb_mwea_mean')]
    sc = ax2.scatter(x, y, s, c, vmin=mb_clim[0], vmax=mb_clim[1], cmap='RdBu', edgecolor='k', linewidth=0.5)

    pltlib.add_cbar(ax2, sc, label='Mass Balance (m we/yr)')
    style_kwd = {'facecolor':'none', 'edgecolor':'k', 'linewidth':0.5}
    #glac_df_mb.plot(ax=ax2, **style_kwd)
    pltlib.add_scalebar(ax2, res=1)
    ax2.set_xlim(extent[0], extent[2])
    ax2.set_ylim(extent[1], extent[3])
    pltlib.hide_ticks(ax2)

    for k, v in glac_df_mb_region.iterrows():
        lbl = '%0.2f +/- %0.2f' % (v[('mb_mwea_mean')], v[('mb_mwea_sigma_mean')])
        ax2.annotate(lbl, xy=(v['centroid_x'],v['centroid_y']), xytext=(1,0), textcoords='offset points', family='sans-serif', fontsize=8, color='k')

plt.show()
