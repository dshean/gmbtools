#! /usr/bin/env python

"""
Combine glacier polygons with geodetic mass balance output
Uses geopandas to join, compute stats for differen regions and plot
"""

#Todo
#Fix issue with int32 conversion after join
#Rivers
#Use 'Name' and HYBAS_ID for index
#Output shp, clean up field names <10 char
#Fix issue with geojson crs definition in output files (currently not written)
#Sum area for all glaciers in each region, not just those with mb numbers

import os, sys
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from imview.lib import pltlib

from osgeo import gdal
from pygeotools.lib import iolib, geolib

import geopandas as gpd
import cartopy.crs as ccrs 

site = 'hma'
#site = 'conus'

area_filter=False
min_area_m2 = 1E6

#Default mb clim
mb_clim = (-1.0, 1.0)
#mb_clim = (-0.7, 0.7)

#HMA
scaling_f = 0.2 
#CONUS
#scaling_f = 3 

region_col = 'region'
basin_col = 'basin'
#basin_col = 'HYBAS_ID'
qdgc_col = 'qdgc'

extent = None
crs = None
if site == 'hma':
    hs_fn = '/Users/dshean/Documents/UW/HMA/mos/hma_mos_32m_20180723/hma_mos_32m_100m_hs_az315.tif'
    glac_shp_fn = '/Users/dshean/data/rgi60/regions/rgi60_merge_HMA_aea.shp'
    #glac_shp_fn = '/Users/dshean/data/rgi60/regions/rgi60_merge_HMA_aea_0.1km.shp'
    #glac_shp_fn = '/Users/dshean/data/rgi60/regions/rgi60_merge_HMA_aea_1km.shp'
    border_shp_fn = '/Users/dshean/data/NaturalEarth/10m_cultural/10m_cultural/ne_10m_admin_0_countries_lakes.shp'
    basin_shp_fn = '/Users/dshean/data/HydroBASINS/hybas_lake_as_lev01-12_v1c/hybas_lake_as_lev04_v1c.shp'
    region_shp_fn = '/Users/dshean/Documents/UW/HMA/Kaab_regions/regions_from_kaab2015_merged.shp'
    qdgc_shp_fn = '/Users/dshean/data/qdgc/qdgc_asia/qdgc_01_asia.shp'
    #This is geopandas crs format
    glac_crs = {u'datum':u'WGS84',u'lat_0':36,u'lat_1':25,u'lat_2':47,u'lon_0':85,u'no_defs':True,u'proj':u'aea',u'units':u'm',u'x_0':0,u'y_0':0}
    #minx, miny, maxx, maxy
    extent = [-1610900, -1142400, 1767400, 1145700]
elif site == 'conus':
    #hs_fn = '/scr/mb/gpd/conus_20171018_mos_32m_trans_100m_hs_az315_1km.tif'
    glac_shp_fn = '/Users/dshean/data/rgi60/regions/rgi60_merge_CONUS_aea.shp'
    region_shp_fn = '/Users/dshean/Documents/UW/CONUS/regions/conus_mb_regions_aea.shp'
    border_shp_fn = '/Users/dshean/data/NaturalEarth/10m_cultural/10m_cultural/ne_10m_admin_1_states_provinces_lakes.shp'
    basin_shp_fn = '/Users/dshean/data/HydroBASINS/hybas_lake_na_lev01-12_v1c/hybas_lake_na_lev07_v1c.shp'
    glac_crs = {u'datum':u'WGS84',u'lat_0':43,u'lat_1':36,u'lat_2':49,u'lon_0':-115,u'no_defs':True,u'proj':u'aea',u'units':u'm',u'x_0':0,u'y_0':0}
    #This is cartopy
    crs = ccrs.AlbersEqualArea(central_longitude=-115, central_latitude=43, standard_parallels=(36, 49))
    extent = [-856800, -789000, 910700, 839400]
else:
    sys.exit("Site not currently supported")

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

def make_map(mb_dissolve_df=None, glac_df_mb=None, region_df=None, col='mb_mwea', border_df=None, \
        basin_df=None, crs=crs, extent=None, hs=None, hs_extent=None, clim=None):

    fig, ax = plt.subplots(figsize=(10,8))
    ax.set_aspect('equal')

    cmap = 'RdBu'
    label = 'Mass Balance (m we/yr)'
    if 't1' in col:
        cmap = 'inferno'
        label = 'Source Date (year)'

    if clim is None:
        clim = (glac_df_mb[col].min(), glac_df_mb[col].max())

    #This is cartopy-enabled axes
    #ax = plt.axes(projection=crs)

    #Currently unsupported for AEA
    #gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')

    if hs is not None:
        print("Plotting image")
        hs_style = {'cmap':'gray', 'origin':'upper', 'extent':cartopy_extent(hs_extent), 'transform':crs}
        ax.imshow(hs, **hs_style)

    if border_df is not None:
        print("Plotting borders")
        border_style = {'facecolor':'0.95','edgecolor':'k', 'linewidth':0.7}
        border_df.plot(ax=ax, **border_style)

    if region_df is not None:
        #This is original region_col
        #region_style = {'column':col_name, 'cmap':'cpt_rainbow', 'edgecolor':'k', 'linewidth':0.5, 'alpha':0.2}
        #region_style = {'column':'Name', 'cmap':'gray', 'edgecolor':'k', 'linewidth':0.5, 'alpha':0.3}
        #region_style = {'cmap':'cpt_rainbow', 'edgecolor':'k', 'linewidth':0.5, 'alpha':0.05}
        region_style = {'facecolor':'none','edgecolor':'k', 'linewidth':0.3, 'alpha':0.4}
        region_df.plot(ax=ax, **region_style)

    if basin_df is not None:
        basin_style = {'facecolor':'none','edgecolor':'k', 'linewidth':0.3, 'alpha':0.4}
        basin_df.plot(ax=ax, **basin_style)

    #https://stackoverflow.com/questions/36008648/colorbar-on-geopandas
    # fake up the array of the scalar mappable. Urgh...
    sc = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=clim[0], vmax=clim[1]))
    sc._A = []

    if mb_dissolve_df is not None:
        #Plot single values for region or basin
        x = mb_dissolve_df['centroid_x']
        y = mb_dissolve_df['centroid_y']
        #s = scaling_f*mb_dissolve_df[('area_m2_sum')]/1E6
        #s = scaling_f*mb_dissolve_df[('Area_sum')]
        s = scaling_f*mb_dissolve_df[('Area_all', '')]
        #c = mb_dissolve_df[('mb_mwea_mean')]
        c = mb_dissolve_df['mb_mwea', 'mean']
        sc_style = {'cmap':cmap, 'edgecolor':'k', 'linewidth':0.5}
        sc = ax.scatter(x, y, s, c, vmin=clim[0], vmax=clim[1], **sc_style) 
        #Add labels
        for k, v in mb_dissolve_df.iterrows():
            #lbl = '%0.2f +/- %0.2f' % (v[('mb_mwea_mean')], v[('mb_mwea_sigma_mean')])
            lbl = '%0.2f +/- %0.2f' % (v[col, 'mean'], v[col+'_sigma','mean'])
            ax.annotate(lbl, xy=(v['centroid_x'],v['centroid_y']), xytext=(1,0), textcoords='offset points', family='sans-serif', fontsize=8, color='k')

    if glac_df_mb is not None:
        print("Plotting glacier polygons")
        glac_style = {'edgecolor':'k', 'linewidth':0.5}
        glac_ax = glac_df_mb.plot(ax=ax, column=col, cmap=cmap, vmin=clim[0], vmax=clim[1], **glac_style)

    #This is minx, miny, maxx, maxy
    if extent is None:
        if glac_df_mb is not None:
            extent = glac_df_mb.total_bounds
        else:
            extent = mb_dissolve_df.total_bounds
    #For cartopy axes
    #ax.set_extent(cartopy_extent(extent), crs=crs)
    ax.set_xlim(extent[0], extent[2])
    ax.set_ylim(extent[1], extent[3])

    #Adding colorbar doesn't work with the cartopy axes
    pltlib.add_cbar(ax, sc, label=label)
    pltlib.add_scalebar(ax, res=1)
    pltlib.hide_ticks(ax)

    plt.tight_layout()

    return fig

#Input csv from mb_parallel.py
mb_csv_fn = sys.argv[1]

driver = 'GeoJSON'
ext = 'geojson'
#driver = 'ESRI Shapefile'
#ext = 'shp'

merge_fn = os.path.splitext(mb_csv_fn)[0]+'_'+os.path.splitext(os.path.split(glac_shp_fn)[-1])[0]+'.'+ext

print("Loading glacier polygons")
glac_shp_join_fn = os.path.splitext(merge_fn)[0]+'_join.'+ext
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

if qdgc_shp_fn is not None:
    print("Loading qdgc")
    qdgc_df = gpd.read_file(qdgc_shp_fn)
    #Convert to glac crs
    qdgc_df = qdgc_df.to_crs(glac_df.crs)

#Add region and basin fields to RGI polygons
#There's a bug here
#https://github.com/Toblerity/Shapely/issues/553
#pip3 uninstall shapely fiona ; pip3 install --no-binary shapely shapely ; pip3 install fiona
#if not os.path.exists(glac_shp_join_fn):
if qdgc_shp_fn is not None:
    print("One-time spatial join by qdgc")
    glac_df = gpd.sjoin(glac_df, qdgc_df, how="inner", op="intersects")

if region_shp_fn is not None:
    print("One-time spatial join by region")
    glac_df = gpd.sjoin(glac_df, region_df, how="inner", op="intersects")
    glac_df.rename(index=str, columns={region_col+"_right": "region", 'index_right':'region_id'}, \
            inplace=True)

if basin_shp_fn is not None:
    print("One-time spatial join by basin")
    glac_df = gpd.sjoin(glac_df, basin_df, how="inner", op="intersects")
    glac_df.rename(index=str, columns={"HYBAS_ID": "basin"}, inplace=True)

#These 'id' are all 0
#glac_df.drop('id', 1)

print("Writing out: %s" % glac_shp_join_fn)
glac_df.to_file(glac_shp_join_fn, driver=driver)

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

print("%i records loaded" % (glac_df_mb.shape[0]))

if area_filter:
    print("Filtering by glacier polygon area (min %0.2f km^2)" % (min_area_m2/1E6))
    orig_count = glac_df_mb.shape[0]
    glac_df_mb = glac_df_mb[glac_df_mb['area_m2'] > min_area_m2]
    print("%i of %i records preserved" % (glac_df_mb.shape[0], orig_count))

#Define aggregation function for the dissolve
aggfunc = {'area_m2':[np.mean, np.sum], 'mb_mwea':[np.mean, np.std, np.sum], 'mb_mwea_sigma':np.mean, 'Area':np.sum}

if region_shp_fn is not None:
    glac_df_region_sum = glac_df.groupby(region_col).sum()
    glac_df_mb_region = glac_df_mb.groupby(region_col).agg(aggfunc)
    glac_df_mb_region['Area_all'] = glac_df_region_sum['Area']
    #glac_df_mb_region = gpd.DataFrame(glac_df_mb_region)
    append_centroid_xy(region_df)
    #region_df.rename(index=str, columns={"Name": "region"}, inplace=True)
    glac_df_mb_region = glac_df_mb_region.merge(region_df[['region', 'centroid_x', 'centroid_y']], left_index=True, right_on='region')

if basin_shp_fn is not None:
    glac_df_basin_sum = glac_df.groupby(basin_col).sum()
    glac_df_mb_basin = glac_df_mb.groupby(basin_col).agg(aggfunc)
    glac_df_mb_basin['Area_all'] = glac_df_basin_sum['Area']
    append_centroid_xy(basin_df)
    basin_df.rename(index=str, columns={"HYBAS_ID": "basin"}, inplace=True)
    glac_df_mb_basin = glac_df_mb_basin.merge(basin_df[['basin', 'centroid_x', 'centroid_y']], left_index=True, right_on='basin')

if False:
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    print("Loading shaded relief map")
    hs_ds = gdal.Open(hs_fn)
    hs = iolib.ds_getma(hs_ds)
    hs_extent = geolib.ds_extent(hs_ds)
    hs_extent_cartopy = cartopy_extent(hs_extent)
    print("Plotting image")
    ax.imshow(hs, cmap='gray', origin='upper', extent=hs_extent_cartopy, transform=crs, alpha=0.6)
else:
    hs = None

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

#Create source date map for CONUS
date_fig = make_map(col='t1', glac_df_mb=glac_df_mb, border_df=border_df, clim=None, crs=crs, extent=extent)

if region_shp_fn is not None:
    region_fig = make_map(col='mb_mwea', mb_dissolve_df=glac_df_mb_region, glac_df_mb=glac_df_mb, region_df=region_df, \
            border_df=border_df, clim=mb_clim, crs=crs, extent=extent)
    region_fig_fn = os.path.splitext(glac_shp_join_fn)[0]+'_region_fig.png'
    region_fig.savefig(region_fig_fn, dpi=300, bbox_inches='tight', pad_inches=0) 
if basin_shp_fn is not None:
    basin_fig = make_map(col='mb_mwea', mb_dissolve_df=glac_df_mb_basin, glac_df_mb=glac_df_mb, basin_df=basin_df, \
            border_df=border_df, clim=mb_clim, crs=crs, extent=extent)
    basin_fig_fn = os.path.splitext(glac_shp_join_fn)[0]+'_basin_fig.png'
    basin_fig.savefig(basin_fig_fn, dpi=300, bbox_inches='tight', pad_inches=0) 

plt.show()
