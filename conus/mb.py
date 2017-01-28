#! /usr/bin/env python
"""
Loop through input polygon features, compute dh/dt and mass balance
"""

import sys
import os

import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal, ogr

from pygeotools.lib import malib
from pygeotools.lib import warplib
from pygeotools.lib import geolib
from pygeotools.lib import iolib

#Glacier shp
#glac_shp_fn = '/nobackupp8/deshean/conus/shp/24k_selection_aea.shp'
#This has already been filtered by area
#glac_shp_fn = '/nobackupp8/deshean/conus/shp/24k_selection_aea_min0.1km2_rainier.shp'
glac_shp_fn = '/nobackupp8/deshean/conus/shp/24k_selection_aea_min0.1km2.shp'

#ogr2ogr -t_srs '+proj=aea +lat_1=36 +lat_2=49 +lat_0=43 +lon_0=-115 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ' 24k_selection_aea.shp 24k_selection_32610.shp

#Raster difference map between NED and WV mosaic
z1_fn = '/nobackupp8/deshean/rpcdem/ned1_2003/ned1_2003_adj.vrt'
z2_fn = '/nobackupp8/deshean/conus/dem2/conus_32611_8m/conus_32611_8m_mos.vrt'
dz_fn = '/nobackup/deshean/conus/dem2/conus_32611_8m/ned1_2003_adj_conus_32611_8m_mos_dz_eul_aea.tif'
#gdalwarp -co TILED=YES -co COMPRESS=LZW -co BIGTIFF=YES -r cubic -t_srs '+proj=aea +lat_1=36 +lat_2=49 +lat_0=43 +lon_0=-115 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ' ned1_2003_adj_conus_32611_8m_mos_dz_eul.tif ned1_2003_adj_conus_32611_8m_mos_dz_eul_aea.tif

#NED 2003 dates
z1_date_shp_fn = '/nobackupp8/deshean/rpcdem/ned1_2003/meta0306_PAL_24k_10kmbuffer_clean_dissolve_aea.shp'
#ogr2ogr -t_srs '+proj=aea +lat_1=36 +lat_2=49 +lat_0=43 +lon_0=-115 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ' meta0306_PAL_24k_10kmbuffer_clean_dissolve_aea.shp meta0306_PAL_24k_10kmbuffer_clean_dissolve_32611.shp

#Field name with scrubbed dates: S_DATE_CLN
z2_date_shp_fn = ''

#Make sure we are in equal area projection
aea_srs = geolib.conus_aea_srs

#Filter glacier poly - let's stick with big glaciers for now
min_glac_area = 0.1 #km^2

f, ax = plt.subplots(figsize=(10,10),dpi=300)
out = []

glac_shp_ds = ogr.Open(glac_shp_fn, 1)
glac_shp_lyr = glac_shp_ds.GetLayer()
feat_count = glac_shp_lyr.GetFeatureCount()
glac_shp_srs = glac_shp_lyr.GetSpatialRef()
glac_shp_lyr.ResetReading()

z1_date_shp_ds = ogr.Open(z1_date_shp_fn)
z1_date_shp_lyr = z1_date_shp_ds.GetLayer()
z1_date_shp_srs = z1_date_shp_lyr.GetSpatialRef()
z1_date_shp_lyr.ResetReading()

print("Processing %i features" % feat_count)

outdir = '/nobackupp8/deshean/conus/dem2/conus_32611_8m/glac_poly_out'
if not os.path.exists(outdir):  
    os.makedirs(outdir)
    
#out_glac_shp_fn = os.path.join(outdir, os.path.splitext(glac_shp_fn)[0]+'_mb.shp')
field_defn = ogr.FieldDefn("mb_mwe", ogr.OFTReal)
glac_shp_lyr.CreateField(field_defn)

for feat in glac_shp_lyr:
    glacname = feat.GetField("GLACNAME")
    glacnum = int(feat.GetField("GLACNUM"))
    if glacname is None:
        glacname = ""
    else:
        glacname = glacname.replace(" ", "")
    feat_fn = "%i_%s" % (glacnum, glacname)
    print("\n%s\n" % feat_fn)
    glac_geom = feat.GetGeometryRef()
    glac_geom.AssignSpatialReference(glac_shp_srs)
    glac_geom_extent = geolib.geom_extent(glac_geom)
    glac_area = glac_geom.GetArea()
    cx, cy = glac_geom.Centroid().GetPoint_2D()

    #ds_list = warplib.memwarp_multi_fn([dz_fn,], res='source', extent=glac_geom_extent, t_srs=aea_srs)
    ds_list = warplib.memwarp_multi_fn([z1_fn, z2_fn], res='min', extent=glac_geom_extent, t_srs=aea_srs, verbose=False)

    glac_geom_mask = geolib.geom2mask(glac_geom, ds_list[0])
    z1 = np.ma.array(iolib.ds_getma(ds_list[0]), mask=glac_geom_mask)
    z2 = np.ma.array(iolib.ds_getma(ds_list[1]), mask=glac_geom_mask)
    dz = z2 - z1
    #dz = np.ma.array(iolib.ds_getma(ds_list[0]), mask=glac_geom_mask)
    if dz.count() == 0:
        print("No valid dz pixels")
        continue

    ds_res = geolib.get_res(ds_list[0])
    valid_area = dz.count()*ds_res[0]*ds_res[1]
    min_valid_area_perc = 0.80
    if valid_area/glac_area < min_valid_area_perc:
        print("Not enough valid pixels. %0.1f%% percent of glacier polygon area")
        continue

    #Rasterize NED source dates
    z1_date_r_ds = iolib.mem_drv.CreateCopy('', ds_list[0])
    gdal.RasterizeLayer(z1_date_r_ds, [1], z1_date_shp_lyr, options=["ATTRIBUTE=S_DATE_CLN"])
    z1_date = iolib.ds_getma(z1_date_r_ds)

    writeout=False
    if writeout:
        dz_r_fn = os.path.join(outdir, feat_fn+'_dz.tif')
        dz_r_ds = iolib.gtif_drv.CreateCopy(dz_r_fn, ds_list[0], options=iolib.gdal_opt)
        iolib.writeGTiff(dz, dz_r_fn, ds_list[0])

        z1_date_r_fn = os.path.join(outdir, feat_fn+'_ned_date.tif')
        iolib.writeGTiff(z1_date, z1_date_r_fn, z1_date_r_ds)

    #Filter dz - throw out abs differences >150 m

    #Compute dz, volume change, mass balance and stats
    #dz_stats = malib.print_stats(dz)
    #z1_stats = malib.print_stats(z1)
    #z2_stats = malib.print_stats(z2)

    #These can be timestamp arrays or datetime objects
    t2 = 2015.0
    t1 = z1_date
    #This is decimal years
    dt = t2 - t1
    #m/yr
    dhdt = dz/dt
    #dhdt_stats = malib.print_stats(dhdt)
    #dhdt_mean = dhdt_stats[3]

    rho_i = 0.91
    rho_s = 0.50
    rho_is = 0.85
    #Can estimate ELA values computed from hypsometry and typical AAR
    ela = None
    if ela is None:
        mb = dhdt * rho_is
    else:
        mb = np.where(dhdt < ela, dhdt*rho_i, dhdt*rho_s)

    mb_stats = malib.print_stats(mb)
    mb_mean = mb_stats[3]
    dmbdt_total_myr = mb_mean*glac_area

    print('%0.2f mwe/yr\n' % mb_mean)
    print('-------------------------------')
    
    out.append([cx, cy, mb_mean, (glac_area/1E6)])
    glac_shp_lyr.SetFeature(feat)
    feat.SetField("mb_mwe", '%0.2f' % mb_mean)
    glac_shp_lyr.SetFeature(feat)

    #Do AED for all
    #Compute mb using scaled AED vs. polygon
    #Extract slope and aspect numbers for polygon
    #Check for valid pixel count vs. feature area, fill if appropriate

    #Error analysis assuming date is wrong by +/- 1-2 years

glac_shp_ds = None

out = np.array(out)
out_stats = malib.print_stats(out[:,2])
if out_stats[5] < 0:
    vmin = out_stats[5]-out_stats[6]*2
    vmax = -vmin 
else:
    vmax = out_stats[5]-out_stats[6]*2
    vmin = -vmax

ax.scatter(out[:,0], out[:,1], c=out[:,2], cmap='RdYlBu', s=out[:,3], vmin=vmin, vmax=vmax)
plt.savefig('test.png')

#Write out new shp with features containing stats
#One shp preserves input features, regardless of source date
#Another shp splits glacier poly based on NED source date

#Plots 
#Scaled circles based on size
#Time periods
