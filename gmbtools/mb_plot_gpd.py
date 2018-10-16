#! /usr/bin/env python

"""
Combine glacier polygons with geodetic mass balance output
Uses geopandas to join, compute stats for differen regions and plot
"""

#Todo
#Plot rivers

import os, sys
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from imview.lib import pltlib

from osgeo import gdal
from pygeotools.lib import iolib, geolib, malib

import geopandas as gpd
import cartopy.crs as ccrs 

pd.options.display.float_format = '{:,.2f}'.format

site = 'hma'
#site = 'conus'

area_filter = False
min_area_m2 = 1E6

outlier_removal = True
outlier_perc = (0.001, 0.999)

#Ocean area
#3.625×108 km2 (Cogley et al., 2011)
Gt2mm = 362.5

plot = False
map_plots = False

#Default mb clim
mb_clim = (-1.0, 1.0)
#mb_clim = (-1.2, 1.2)
#mb_clim = (-0.7, 0.7)

#suptitle = "Glacier Mass Balance (ASTER 2000–2009)"
suptitle = "Glacier Mass Balance (ASTER 2000–2018)"
#suptitle = "Glacier Mass Balance (ASTER 2009–2018)"
#suptitle = "Glacier Mass Balance (SRTM 2000 to WV/GE median)"

#HMA
scaling_f = 0.2 
#CONUS
#scaling_f = 3 

rgi_col = 'RGIId'
region_col = 'region'
basin_col = 'HYBAS_ID'
qdgc_col = 'qdgc'
mascon_col = 'mascon'

extent = None
crs = None
if site == 'hma':
    basin_col = 'basin_name'
    hs_fn = '/Users/dshean/Documents/UW/HMA/mos/hma_mos_32m_20180723/hma_mos_32m_100m_hs_az315.tif'
    glac_shp_fn = '/Users/dshean/data/rgi60/regions/rgi60_merge_HMA_aea.shp'
    #glac_shp_fn = '/Users/dshean/data/rgi60/regions/rgi60_merge_HMA_aea_0.1km.shp'
    #glac_shp_fn = '/Users/dshean/data/rgi60/regions/rgi60_merge_HMA_aea_1km.shp'
    border_shp_fn = '/Users/dshean/data/NaturalEarth/10m_cultural/10m_cultural/ne_10m_admin_0_countries_lakes.shp'
    #basin_shp_fn = '/Users/dshean/data/HydroBASINS/hybas_lake_as_lev01-12_v1c/hybas_lake_as_lev04_v1c.shp'
    basin_shp_fn = '/Users/dshean/data/HydroBASINS/HiMAT_full_210_IDs_subset_merged_clip_names.gpkg'
    region_shp_fn = '/Users/dshean/Documents/UW/HMA/Kaab_regions/regions_from_kaab2015_merged_clean.shp'
    #http://www.mindland.com/wp/download-qdgc-continents/
    qdgc_shp_fn = '/Users/dshean/data/qdgc/qdgc_asia/qdgc_01_asia.shp'
    mascon_shp_fn = '/Users/dshean/data/grace_mascons/GSFC.glb.200301_201607_v02.4_clip.gpkg'
    #This is geopandas crs format
    glac_crs = {u'datum':u'WGS84',u'lat_0':36,u'lat_1':25,u'lat_2':47,u'lon_0':85,u'no_defs':True,u'proj':u'aea',u'units':u'm',u'x_0':0,u'y_0':0}
    #minx, miny, maxx, maxy
    extent = [-1610900, -1142400, 1767400, 1145700]
elif site == 'conus':
    #hs_fn = '/scr/mb/gpd/conus_20171018_mos_32m_trans_100m_hs_az315_1km.tif'
    glac_shp_fn = '/Users/dshean/data/rgi60/regions/rgi60_merge_CONUS_aea.shp'
    region_shp_fn = '/Users/dshean/Documents/UW/CONUS/regions/conus_mb_regions.shp'
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

def aggregate(glac_df, glac_df_mb, agg_df, col):
    #Define aggregation function for the dissolve
    aggfunc = {'area_m2':[np.mean, np.sum], 'mb_mwea':[np.mean, np.std, np.sum, np.size], 'mb_mwea_sigma':np.mean, 'mb_m3wea':[np.mean,np.sum], 'mb_Gta':[np.sum], 'Area':[np.mean, np.sum], 't1':[np.mean, np.min, np.max], 't2':[np.mean, np.min, np.max], 'dt':[np.mean, np.min, np.max]}
    #This is for all glaciers - mostly just to get total area
    glac_df_agg_sum = glac_df.groupby(col).sum()
    glac_df_agg_mean = glac_df.groupby(col).mean()
    #glac_df_agg_mean = glac_df.groupby(col).median()
    #Perform the aggregation
    glac_df_mb_agg = glac_df_mb.groupby(col).agg(aggfunc)
    #This is count of number of glaciers
    glac_df_mb_agg[('mb_mwea', 'count')] = glac_df_mb_agg[('mb_mwea', 'size')].astype('int64') 
    #These are values for all glaciers in each region (not just those with mb numbers)
    glac_df_mb_agg[('Area_all', 'mean')] = glac_df_agg_mean['Area']
    glac_df_mb_agg[('Area_all', 'sum')] = glac_df_agg_sum['Area']
    #Percent coverage of mb numbers
    glac_df_mb_agg[('Area','perc')] = 100*glac_df_mb_agg[('Area','sum')]/glac_df_mb_agg[('Area_all', 'sum')]
    glac_df_mb_agg[('mb_mwea', 'total_m3a')] = glac_df_mb_agg[('mb_mwea', 'mean')] * glac_df_mb_agg[('Area_all', 'sum')] * 1E6
    glac_df_mb_agg[('mb_mwea', 'total_Gta')] = glac_df_mb_agg[('mb_mwea', 'total_m3a')]/1E9 
    glac_df_mb_agg[('mb_mwe_cum', 'mean')] = glac_df_mb_agg[('mb_mwea','mean')] * glac_df_mb_agg[('dt', 'mean')]
    glac_df_mb_agg[('mb_mwe_cum', 'total_m3')] = glac_df_mb_agg[('mb_mwea', 'total_m3a')] * glac_df_mb_agg[('dt', 'mean')] 
    glac_df_mb_agg[('mb_mwe_cum', 'total_Gt')] = glac_df_mb_agg[('mb_mwea', 'total_Gta')] * glac_df_mb_agg[('dt', 'mean')] 
    #Compute numbers for meltwater (polygons with mb < 0)
    #glac_df_mb_agg_meltwater = glac_df_mb[glac_df_mb['mb_mwea' < 0]].sum()
    glac_df_mb_agg_meltwater = glac_df_mb[glac_df_mb['mb_mwea'] < 0].groupby(col).agg(aggfunc)
    glac_df_mb_agg[('meltwater', 'count')] = glac_df_mb_agg_meltwater[('mb_mwea', 'size')].astype('int64') 
    glac_df_mb_agg[('meltwater', 'total_m3a')] = glac_df_mb_agg_meltwater[('mb_m3wea', 'sum')]
    glac_df_mb_agg[('meltwater', 'total_Gta')] = glac_df_mb_agg[('meltwater', 'total_m3a')]/1E9 
    glac_df_mb_agg[('meltwater', 'total_mmSLEa')] = glac_df_mb_agg[('meltwater', 'total_Gta')]/Gt2mm
    glac_df_mb_agg[('meltwater_cum', 'total_m3')] = glac_df_mb_agg[('meltwater', 'total_m3a')] * glac_df_mb_agg[('dt', 'mean')] 
    glac_df_mb_agg[('meltwater_cum', 'total_Gt')] = glac_df_mb_agg[('meltwater', 'total_Gta')] * glac_df_mb_agg[('dt', 'mean')] 
    glac_df_mb_agg[('meltwater_cum', 'total_mmSLE')] = glac_df_mb_agg[('meltwater', 'total_Gta')] * glac_df_mb_agg[('dt', 'mean')] 
    append_centroid_xy(agg_df)
    if 'basin' in col:
        #Preserve basin attributes (endorheic flag, discharge)
        glac_df_mb_agg = glac_df_mb_agg.merge(agg_df[['ENDO', 'centroid_x', 'centroid_y']], left_index=True, right_index=True)
    else:
        glac_df_mb_agg = glac_df_mb_agg.merge(agg_df[['centroid_x', 'centroid_y']], left_index=True, right_index=True)
    glac_df_mb_agg.sort_values(by=('Area_all', 'sum'), ascending=False, inplace=True)
    glac_df_mb_agg.df_name = col
    return glac_df_mb_agg

def add_legend(ax, sf=16, loc='upper right'):
    """
    Create legend for scaled scatterplot markers
    """
    ax.autoscale(False)
    #CONUS
    #leg_s = np.array([0.1, 0.5, 1.0, 5.0, 10.0])
    #HMA
    #leg_s = np.array([0.1, 1.0, 10.0, 100.0])
    leg_s = np.array([1.0, 10.0, 100.0, 1000.0, 10000.0])
    #Spoof dummy coordinates way off map
    leg_x = np.full(leg_s.size, -999999999)
    leg_y = np.full(leg_s.size, -999999999)
    for i, s in enumerate(leg_s):
        #lbl = r'$%0.1f\/km^2$' % s
        lbl = '%i' % s
        ax.scatter(leg_x[i], leg_y[i], s=s*sf, c='gray', label=lbl)
    legend = ax.legend(title=r'$Glacier\/Area\/(km^2)$', scatterpoints=1, loc=loc, prop={'size':7}, ncol=leg_s.size)
    legend.get_title().set_fontsize('8')
    return legend

def make_map(mb_dissolve_df=None, glac_df_mb=None, agg_df=None, col=('mb_mwea', 'mean'), border_df=None, crs=crs, extent=None, hs=None, hs_extent=None, clim=None, labels='val', title=None):

    fig, ax = plt.subplots(figsize=(10,8))
    ax.set_aspect('equal')
    legend = add_legend(ax, sf=scaling_f)
    if title is not None:
        ax.set_title(title)

    if clim is None:
        #clim = (glac_df_mb[col].min(), glac_df_mb[col].max())
        clim = malib.calcperc_sym(mb_dissolve_df[col], perc=(1,99))

    cmap = 'RdBu'
    if 'mb_mwea' in col:
        label = 'Mass Balance (m we/yr)'
    elif 'mb_Gta' in col: 
        label = 'Mass Balance (Gt/yr)'
    elif 'meltwater' in col: 
        label = 'Excess Meltwater Runoff (Gt/yr)'
        #Reverse, as these are negative values
        cmap = 'YlOrRd_r'
        #cmap = 'inferno'
        clim = malib.calcperc(mb_dissolve_df[col], perc=(0,99))
    elif 't1' in col:
        cmap = 'inferno'
        label = 'Source Date (year)'

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
        border_style = {'facecolor':'0.65','edgecolor':'k', 'linewidth':0.7}
        border_df.plot(ax=ax, **border_style)

    if agg_df is not None:
        print("Plotting agg boundaries")
        #This was to get colored regions
        #agg_style = {'cmap':'cpt_rainbow', 'edgecolor':'none', 'linewidth':0, 'alpha':0.05}
        agg_style = {'cmap':'summer', 'edgecolor':'none', 'linewidth':0, 'alpha':0.05}
        #agg_style = {'facecolor':'0.95','edgecolor':'k', 'linewidth':0.3, 'alpha':0.2}
        agg_df.plot(ax=ax, **agg_style)

    if glac_df_mb is not None:
        print("Plotting glacier polygons")
        glac_style = {'edgecolor':'k', 'linewidth':0.1, 'alpha':0.2}
        #This plots mb color ramp for each glacier polygon
        #glac_ax = glac_df_mb.plot(ax=ax, column=col[0], cmap=cmap, vmin=clim[0], vmax=clim[1], **glac_style)
        #This plots outlines
        glac_ax = glac_df_mb.plot(ax=ax, facecolor='none', **glac_style)

    if agg_df is not None:
        agg_style = {'facecolor':'none', 'edgecolor':'w', 'linewidth':0.5}
        agg_df.plot(ax=ax, **agg_style)

    #https://stackoverflow.com/questions/36008648/colorbar-on-geopandas
    # fake up the array of the scalar mappable so we can plot colorbar. Urgh...
    sc = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=clim[0], vmax=clim[1]))
    sc._A = []

    if mb_dissolve_df is not None:
        print("Plotting scatterplot of %s values" % (col, ))
        #Plot single values for region or basin
        x = mb_dissolve_df['centroid_x']
        y = mb_dissolve_df['centroid_y']
        #Scale by total glacier area in each polygon 
        s = scaling_f*mb_dissolve_df[('Area_all', 'sum')]
        c = mb_dissolve_df[col]
        sc_style = {'cmap':cmap, 'edgecolor':'k', 'linewidth':0.5, 'alpha':0.8}
        sc = ax.scatter(x, y, s, c, vmin=clim[0], vmax=clim[1], **sc_style) 
        #Add labels
        text_kw = {'family':'sans-serif', 'fontsize':8, 'color':'k'}
        if labels is not None:
            print("Adding annotations")
            for k, v in mb_dissolve_df.iterrows():
                #lbl = '%0.2f +/- %0.2f' % (v[col], v[(col[0]+'_sigma',col[1])])
                if labels == 'name+val':
                    lbl = '%s\n%+0.2f' % (k, v[col])
                else:
                    lbl = '%+0.2f' % v[col]
                #ax.annotate(lbl, xy=(v['centroid_x'],v['centroid_y']), xytext=(1,0), textcoords='offset points', family='sans-serif', fontsize=6, color='darkgreen')
                txt = ax.annotate(lbl, xy=(v['centroid_x'],v['centroid_y']), ha='center', va='center', **text_kw)
                txt.set_path_effects([path_effects.Stroke(linewidth=0.75, foreground='w'),path_effects.Normal()])

    #This is minx, miny, maxx, maxy
    if extent is None:
        #if glac_df_mb is not None:
        #    extent = glac_df_mb.total_bounds
        #else:
        extent = mb_dissolve_df.total_bounds

    #For cartopy axes
    #ax.set_extent(cartopy_extent(extent), crs=crs)
    #Pad extent so labels fit within map
    #extent = geolib.pad_extent(extent, perc=0.01, uniform=True)
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

#driver = 'GeoJSON'
#ext = 'geojson'
#driver = 'ESRI Shapefile'
#ext = 'shp'
driver = 'GPKG'
ext= 'gpkg'

merge_fn = os.path.splitext(mb_csv_fn)[0]+'_'+os.path.splitext(os.path.split(glac_shp_fn)[-1])[0]+'.'+ext

glac_shp_join_fn = os.path.splitext(merge_fn)[0]+'_join.'+ext

if os.path.exists(glac_shp_join_fn):
    print("Loading glaciers polygons joined with regions, basins, qdgc, mascons")
    glac_df = gpd.read_file(glac_shp_join_fn)
    #This is a hack, as geojson doesn't properly preserve custom aae proj
    glac_df.crs = glac_crs
else:
    print("Loading glacier polygons")
    glac_df = gpd.read_file(glac_shp_fn)
    glac_df.set_index(rgi_col, inplace=True)
    #Add centroid field - needed for proper partitioning in spatial join
    glac_df['centroid_geom'] = gpd.GeoSeries(glac_df.centroid)
    glac_df['polygon_geom'] = gpd.GeoSeries(glac_df.geometry)
    glac_df.set_geometry('centroid_geom', inplace=True, drop=True)

#Area in km2
#glac_df.geometry.area.sum()/1E6
#glac_df['Area'].sum()

if region_shp_fn is not None:
    print("Loading regions")
    region_df = gpd.read_file(region_shp_fn)
    #Convert to glac crs
    region_df = region_df.to_crs(glac_df.crs)
    region_df.set_index(region_col, inplace=True)

if basin_shp_fn is not None:
    print("Loading basins")
    basin_df = gpd.read_file(basin_shp_fn)
    #There are issues with fiona writing large integers in basin IDs, convert to string to be safe
    #The dtype is int64, but finoa writes int32, and they all get truncated to 2147483647
    if basin_col == 'HYBAS_ID':
        for f in ['HYBAS_ID', 'NEXT_DOWN', 'NEXT_SINK', 'MAIN_BAS']:
            basin_df[f] = basin_df[f].astype('str')
    #Convert to glac crs
    basin_df = basin_df.to_crs(glac_df.crs)
    basin_df.set_index(basin_col, inplace=True)

if mascon_shp_fn is not None:
    print("Loading mascon")
    mascon_df = gpd.read_file(mascon_shp_fn)
    #Convert to glac crs
    mascon_df = mascon_df.to_crs(glac_df.crs)
    #Add a unique identifier
    mascon_df[mascon_col] = mascon_df.lat_center.map('{:,.0f}N'.format) + mascon_df.lon_center.map('{:,.0f}E'.format)
    mascon_df.set_index(mascon_col, inplace=True)

if qdgc_shp_fn is not None:
    print("Loading qdgc")
    qdgc_df = gpd.read_file(qdgc_shp_fn)
    #Convert to glac crs
    qdgc_df = qdgc_df.to_crs(glac_df.crs)
    qdgc_df.set_index(qdgc_col, inplace=True)

print(glac_df.shape)
#Add region, basin, etc fields to RGI polygons
if not os.path.exists(glac_shp_join_fn):
    if mascon_shp_fn is not None:
        print("One-time spatial join by mascon")
        glac_df = gpd.sjoin(glac_df, mascon_df, how="inner", op="intersects")
        glac_df.rename(index=str, columns={'index_right':mascon_col}, inplace=True)
        print(glac_df.shape)

    if qdgc_shp_fn is not None:
        print("One-time spatial join by qdgc")
        glac_df = gpd.sjoin(glac_df, qdgc_df, how="inner", op="intersects")
        glac_df.rename(index=str, columns={'index_right':qdgc_col}, inplace=True)
        print(glac_df.shape)

    if region_shp_fn is not None:
        print("One-time spatial join by region")
        glac_df = gpd.sjoin(glac_df, region_df, how="inner", op="intersects")
        glac_df.rename(index=str, columns={'index_right':region_col}, inplace=True)
        print(glac_df.shape)

    if basin_shp_fn is not None:
        print("One-time spatial join by basin")
        glac_df = gpd.sjoin(glac_df, basin_df, how="inner", op="intersects")
        glac_df.rename(index=str, columns={'index_right':basin_col}, inplace=True)
        print(glac_df.shape)

    glac_df.set_geometry('polygon_geom', inplace=True, drop=True)
    print("Writing out: %s" % glac_shp_join_fn)
    #With index set to RGIId, it is not written out, hack to create new column
    glac_df.reset_index().rename(columns={'index':rgi_col}).to_file(glac_shp_join_fn, driver=driver)

print("Loading mb")
mb_df = pd.read_csv(mb_csv_fn)
mb_df[rgi_col] = 'RGI60-'+mb_df[rgi_col].map('{:08.5f}'.format)
mb_df.set_index(rgi_col, inplace=True)
mb_df['mb_Gta'] = mb_df['mb_m3wea']/1E9
mb_df_t1 = mb_df['t1'].mean()
mb_df_t2 = mb_df['t2'].mean()
print(mb_df_t1, mb_df_t2)
dt_str = '%.0f-%.0f' % (mb_df_t1, mb_df_t2)

#mb_df['mb_mwea_area'] = mb_df['mb_mwea'] * mb_df['area_km2']*1E6

if os.path.exists(merge_fn):
    print("Loading merged polygons and mb")
    glac_df_mb = gpd.read_file(merge_fn)
    glac_df_mb.set_index(rgi_col, inplace=True)
    print(glac_df_mb.shape)
else:
    print("Merging glacier polygons and mb results")
    glac_df_mb = glac_df.merge(mb_df, left_index=True, right_index=True)
    #With index set to RGIId, it is not written out, hack to create new column
    print(glac_df_mb.shape)
    print("Writing out: %s" % merge_fn)
    glac_df_mb.reset_index().rename(columns={'index':rgi_col}).to_file(merge_fn, driver=driver)

print("%i merged records loaded" % (glac_df_mb.shape[0]))

if area_filter:
    print("Filtering by glacier polygon area (min %0.2f km^2)" % (min_area_m2/1E6))
    orig_count = glac_df_mb.shape[0]
    glac_df_mb = glac_df_mb[glac_df_mb['area_m2'] > min_area_m2]
    print("%i of %i records preserved" % (glac_df_mb.shape[0], orig_count))

if outlier_removal:
    print("Removing outliers")
    outlier_clim = (glac_df_mb['mb_mwea'].quantile(outlier_perc[0]), glac_df_mb['mb_mwea'].quantile(outlier_perc[1]))

    if plot:
        plt.figure()
        ax = plt.gca()
        ax.set_xlabel('Mass balance (m we/yr)')
        ax.set_ylabel('Number of glaciers')
        hist_clim = (-1.0, 1.0)
        ax.set_xlim(*hist_clim)
        print("%i records before outlier removal" % (glac_df_mb.shape[0]))
        glac_df_mb['mb_mwea'].hist(range=hist_clim, bins=256, label='Before outlier filter')
        inlier_idx = np.abs(glac_df_mb['mb_mwea'] - glac_df_mb['mb_mwea'].mean()) <= (3*glac_df_mb['mb_mwea'].std())
        glac_df_mb = glac_df_mb[inlier_idx]
        print("%i records after outlier removal" % (glac_df_mb.shape[0]))
        glac_df_mb['mb_mwea'].hist(range=hist_clim, bins=256, label='After outlier filter')
        ax.axvline(0, linewidth=0.5, color='k')
        ax.legend()

if mascon_shp_fn is not None:
    glac_df_mb_mascon = aggregate(glac_df, glac_df_mb, mascon_df, mascon_col)

if qdgc_shp_fn is not None:
    glac_df_mb_qdgc = aggregate(glac_df, glac_df_mb, qdgc_df, qdgc_col)

if region_shp_fn is not None:
    glac_df_mb_region = aggregate(glac_df, glac_df_mb, region_df, region_col)

if basin_shp_fn is not None:
    glac_df_mb_basin = aggregate(glac_df, glac_df_mb, basin_df, basin_col)
    glac_df_mb_basin_exo = glac_df_mb_basin[glac_df_mb_basin['ENDO'] == 0]

#Compile stats for all glaciers
all_stats = glac_df_mb[['mb_mwea','mb_m3wea']].mean()
#all_stats = glac_df_mb[['mb_mwea','mb_m3wea']].median()
all_stats['mb_mwea_std'] = glac_df_mb[['mb_mwea','mb_mwea']].std()
all_stats['mb_m3wea_std'] = glac_df_mb[['mb_mwea','mb_m3wea']].std()
#glac_df_mb[['mb_mwea','mb_m3wea']].apply(malib.mad, axis=0)
#all_stats = glac_df_mb[['mb_mwea']].mean()
#Total sampled area
all_stats['sample_area'] = glac_df_mb['Area'].sum()
#All area in RGI db
all_stats['RGI_total_area'] = glac_df['Area'].sum()
all_stats['sample_area_perc'] = all_stats['sample_area']/all_stats['RGI_total_area']
all_stats['mb_m3wea_sum'] = glac_df_mb['mb_m3wea'].sum()
all_stats['mb_m3wea_all'] = all_stats['mb_mwea'] * all_stats['RGI_total_area'] * 1E6
all_stats['mb_Gta'] = all_stats['mb_m3wea']/1E9
all_stats['mb_Gta_sum'] = all_stats['mb_m3wea_sum']/1E9
all_stats['mb_Gta_all'] = all_stats['mb_m3wea_all']/1E9
dt = glac_df_mb['dt'].mean()
all_stats_cum = all_stats * dt

print("All glaciers, rate")
print(all_stats)
#print("All glaciers, cumulative")
#print(all_stats_cum)

print("\nRegions")
region_summary = glac_df_mb_region[[('mb_mwea', 'count'),('mb_mwea', 'mean'),('mb_mwea', 'std'),('Area', 'sum'),('Area_all', 'sum'),('Area', 'perc'),('mb_mwea', 'total_Gta')]]
print(region_summary.to_string())
out_fn = os.path.splitext(glac_shp_join_fn)[0]+'__%s_region_summary.pkl' % dt_str
region_summary.to_pickle(out_fn)

print("\nBasins")
basin_summary = glac_df_mb_basin[[('meltwater', 'count'),('meltwater', 'total_m3a'),('meltwater', 'total_Gta'),('meltwater', 'total_mmSLEa')]]
print(basin_summary.to_string())
print(glac_df_mb_basin[[('meltwater', 'total_Gta'),('meltwater', 'total_mmSLEa')]].sum())
out_fn = os.path.splitext(glac_shp_join_fn)[0]+'__%s_basin_summary.pkl' % dt_str
basin_summary.to_pickle(out_fn)

print("\nExorheic Basins (SLR contribution)")
print(glac_df_mb_basin_exo[[('meltwater', 'count'),('meltwater', 'total_m3a'),('meltwater', 'total_Gta'),('meltwater', 'total_mmSLEa')]].to_string())
print(glac_df_mb_basin_exo[[('meltwater', 'total_Gta'),('meltwater', 'total_mmSLEa')]].sum())

#Compile stats for each division
print("\nTotal Gt/a for each aggregation")
print('agg', 'total_Gta', 'mb_mwea')
for i in [glac_df_mb_basin, glac_df_mb_region, glac_df_mb_qdgc]:
    print(i.df_name, i[('mb_mwea', 'total_Gta')].sum(), i['mb_mwea', 'mean'].mean())

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
#date_fig = make_map(col='t1', glac_df_mb=glac_df_mb, border_df=border_df, clim=None, crs=crs, extent=extent)

mascon_df_out = glac_df_mb_mascon[[('mb_mwea','count'),('mb_mwea', 'mean'), ('mb_mwea', 'std'), ('Area','sum'),('mb_Gta', 'sum'), ('Area_all', 'sum'), ('mb_mwea', 'total_Gta')]]
header = ('n_glaciers','mb_mwea_mean','mb_mwea_std','obs_glacier_area_km2','mb_Gta_sum','all_glacier_area_km2','mb_Gta_all')
mascon_csv_fn = os.path.splitext(glac_shp_join_fn)[0]+'_mascon.csv'
mascon_df_out.to_csv(mascon_csv_fn, float_format='%0.2f',header=header)

if map_plots:
    if basin_shp_fn is not None:
        basin_melt_gt_clim = (-5, 0)
        basin_fig_fn = os.path.splitext(glac_shp_join_fn)[0]+'_basin_excess_Gt_fig.png'
        print("Generating figure: %s" % basin_fig_fn)
        if basin_col == 'HYBAS_ID':
            title = suptitle + ": HydroBASINS level 4" 
        else:
            title = suptitle + ": Excess meltwater runoff from Glacier Mass Loss"
        #minx, miny, maxx, maxy
        basin_extent = [-2396534,-2619071, 3273634,2008000]
        #basin_extent = basin_df.total_bounds
        basin_fig = make_map(col=('meltwater', 'total_Gta'), mb_dissolve_df=glac_df_mb_basin, glac_df_mb=glac_df_mb, agg_df=basin_df, border_df=border_df, clim=basin_melt_gt_clim, crs=crs, extent=basin_extent, labels='name+val', title=title)
        #basin_fig = make_map(col=('mb_Gta', 'sum'), mb_dissolve_df=glac_df_mb_basin, glac_df_mb=glac_df_mb, agg_df=basin_df, border_df=border_df, clim=None, crs=crs, extent=extent, labels='val', title=title)
        print("Saving figure: %s" % basin_fig_fn)
        #basin_fig.savefig(basin_fig_fn, dpi=300, bbox_inches='tight', pad_inches=0) 
        basin_fig.savefig(basin_fig_fn, dpi=300, pad_inches=0) 

    if basin_shp_fn is not None:
        basin_mb_clim = (-0.5, 0.5)
        basin_fig_fn = os.path.splitext(glac_shp_join_fn)[0]+'_basin_mwe_fig.png'
        print("Generating figure: %s" % basin_fig_fn)
        if basin_col == 'HYBAS_ID':
            title = suptitle + ": HydroBASINS level 4" 
        else:
            title = suptitle + ": Basins"
        #minx, miny, maxx, maxy
        basin_extent = [-2396534,-2619071, 3273634,2008000]
        basin_fig = make_map(col=('mb_mwea', 'mean'), mb_dissolve_df=glac_df_mb_basin, glac_df_mb=glac_df_mb, agg_df=basin_df, border_df=border_df, clim=basin_mb_clim, crs=crs, extent=basin_extent, labels='name+val', title=title)
        print("Saving figure: %s" % basin_fig_fn)
        #basin_fig.savefig(basin_fig_fn, dpi=300, bbox_inches='tight', pad_inches=0) 
        basin_fig.savefig(basin_fig_fn, dpi=300, pad_inches=0) 

    if qdgc_shp_fn is not None:
        qdgc_fig_fn = os.path.splitext(glac_shp_join_fn)[0]+'_qdgc_mwe_fig.png'
        print("Generating figure: %s" % qdgc_fig_fn)
        #To plot grid cells, pass agg_df=qdgc_df
        title = suptitle + ": Quarter-degree Grid Cells"
        qdgc_fig = make_map(col=('mb_mwea', 'mean'), mb_dissolve_df=glac_df_mb_qdgc, glac_df_mb=glac_df_mb, agg_df=None, border_df=border_df, clim=mb_clim, crs=crs, extent=extent, labels=None, title=title)
        print("Saving figure: %s" % qdgc_fig_fn)
        qdgc_fig.savefig(qdgc_fig_fn, dpi=300, bbox_inches='tight', pad_inches=0) 

    if region_shp_fn is not None:
        region_gt_clim = (-3.0, 3.0)
        region_fig_fn = os.path.splitext(glac_shp_join_fn)[0]+'_region_Gt_fig.png'
        print("Generating figure: %s" % region_fig_fn)
        title = suptitle + ": Kaab Regions"
        #region_fig = make_map(col='mb_mwea', mb_dissolve_df=glac_df_mb_region, glac_df_mb=glac_df_mb, agg_df=region_df, border_df=border_df, clim=mb_clim, crs=crs, extent=extent, labels='name+val')
        region_fig = make_map(col=('mb_Gta', 'sum'), mb_dissolve_df=glac_df_mb_region, glac_df_mb=glac_df_mb, agg_df=region_df, border_df=border_df, clim=region_gt_clim, crs=crs, extent=extent, labels='name+val', title=title)
        print("Saving figure: %s" % region_fig_fn)
        region_fig.savefig(region_fig_fn, dpi=300, bbox_inches='tight', pad_inches=0) 

    if region_shp_fn is not None:
        region_mb_clim=(-0.5, 0.5)
        region_fig_fn = os.path.splitext(glac_shp_join_fn)[0]+'_region_mwe_fig.png'
        print("Generating figure: %s" % region_fig_fn)
        title = suptitle + ": Kaab Regions"
        #region_fig = make_map(col='mb_mwea', mb_dissolve_df=glac_df_mb_region, glac_df_mb=glac_df_mb, agg_df=region_df, border_df=border_df, clim=mb_clim, crs=crs, extent=extent, labels='name+val')
        region_fig = make_map(col=('mb_mwea', 'mean'), mb_dissolve_df=glac_df_mb_region, glac_df_mb=glac_df_mb, agg_df=region_df, border_df=border_df, clim=region_mb_clim, crs=crs, extent=extent, labels='name+val', title=title)
        print("Saving figure: %s" % region_fig_fn)
        region_fig.savefig(region_fig_fn, dpi=300, bbox_inches='tight', pad_inches=0) 

    if mascon_shp_fn is not None:
        mascon_fig_fn = os.path.splitext(glac_shp_join_fn)[0]+'_mascon_mwe_fig.png'
        print("Generating figure: %s" % mascon_fig_fn)
        #To plot grid cells, pass agg_df=mascon_df
        title = suptitle + ": GSFC GRACE Mascons"
        mascon_fig = make_map(col=('mb_mwea', 'mean'), mb_dissolve_df=glac_df_mb_mascon, glac_df_mb=glac_df_mb, agg_df=None, border_df=border_df, clim=mb_clim, crs=crs, extent=extent, labels=None, title=title)
        print("Saving figure: %s" % mascon_fig_fn)
        mascon_fig.savefig(mascon_fig_fn, dpi=300, bbox_inches='tight', pad_inches=0) 

    if mascon_shp_fn is not None:
        mascon_gt_clim = (-0.5, 0.5)
        mascon_fig_fn = os.path.splitext(glac_shp_join_fn)[0]+'_mascon_Gt_fig.png'
        print("Generating figure: %s" % mascon_fig_fn)
        #To plot grid cells, pass agg_df=mascon_df
        title = suptitle + ": GSFC GRACE Mascons"
        mascon_fig = make_map(col=('mb_Gta', 'sum'), mb_dissolve_df=glac_df_mb_mascon, glac_df_mb=glac_df_mb, agg_df=None, border_df=border_df, clim=mascon_gt_clim, crs=crs, extent=extent, labels=None, title=title)
        print("Saving figure: %s" % mascon_fig_fn)
        mascon_fig.savefig(mascon_fig_fn, dpi=300, bbox_inches='tight', pad_inches=0) 

