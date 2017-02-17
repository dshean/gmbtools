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
#z2_fn = '/nobackupp8/deshean/conus/dem2/conus_32611_8m/conus_32611_8m_mos.vrt'
#z2_fn = '/nobackupp8/deshean/conus/dem2/conus_8m_tile_coreg_round3/conus_8m_tile_coreg_round3.vrt'
z2_fn = '/nobackupp8/deshean/conus/dem2/conus_8m_tile_coreg_round3_summer2014-2016/conus_8m_tile_coreg_round3_summer2014-2016.vrt'
#dz_fn = '/nobackup/deshean/conus/dem2/conus_32611_8m/ned1_2003_adj_conus_32611_8m_mos_dz_eul_aea.tif'
#gdalwarp -co TILED=YES -co COMPRESS=LZW -co BIGTIFF=YES -r cubic -t_srs '+proj=aea +lat_1=36 +lat_2=49 +lat_0=43 +lon_0=-115 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ' ned1_2003_adj_conus_32611_8m_mos_dz_eul.tif ned1_2003_adj_conus_32611_8m_mos_dz_eul_aea.tif

#NED 2003 dates
z1_date_shp_fn = '/nobackupp8/deshean/rpcdem/ned1_2003/meta0306_PAL_24k_10kmbuffer_clean_dissolve_aea.shp'
#ogr2ogr -t_srs '+proj=aea +lat_1=36 +lat_2=49 +lat_0=43 +lon_0=-115 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ' meta0306_PAL_24k_10kmbuffer_clean_dissolve_aea.shp meta0306_PAL_24k_10kmbuffer_clean_dissolve_32611.shp

#PRISM climate data, 800-m 
prism_ppt_fn = '/nobackupp8/deshean/conus/prism/normals/annual/ppt/PRISM_ppt_30yr_normal_800mM2_annual_bil.bil'
prism_tmean_fn = '/nobackupp8/deshean/conus/prism/normals/annual/tmean/PRISM_tmean_30yr_normal_800mM2_annual_bil.bil'

#Field name with scrubbed dates: S_DATE_CLN
z2_date_shp_fn = ''

#Make sure we are in equal area projection
aea_srs = geolib.conus_aea_srs

#Filter glacier poly - let's stick with big glaciers for now
min_glac_area = 0.1 #km^2

#List to hold output
out = []

glac_shp_ds = ogr.Open(glac_shp_fn, 0)
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
    #ds_list = warplib.memwarp_multi_fn([z1_fn, z2_fn], res='min', extent=glac_geom_extent, t_srs=aea_srs, verbose=False)
    ds_list = warplib.memwarp_multi_fn([z1_fn, z2_fn, prism_ppt_fn, prism_tmean_fn], res='min', extent=glac_geom_extent, t_srs=aea_srs, verbose=False)

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
    valid_area_perc = valid_area/glac_area
    min_valid_area_perc = 0.80
    if valid_area_perc < min_valid_area_perc:
        print("Not enough valid pixels. %0.1f%% percent of glacier polygon area" % valid_area_perc)
        continue

    #Rasterize NED source dates
    z1_date_r_ds = iolib.mem_drv.CreateCopy('', ds_list[0])
    gdal.RasterizeLayer(z1_date_r_ds, [1], z1_date_shp_lyr, options=["ATTRIBUTE=S_DATE_CLN"])
    z1_date = np.ma.array(iolib.ds_getma(z1_date_r_ds), mask=glac_geom_mask)

    #Filter dz - throw out abs differences >150 m

    #Compute dz, volume change, mass balance and stats
    #dz_stats = malib.print_stats(dz)
    z1_stats = malib.print_stats(z1)
    z2_stats = malib.print_stats(z2)
    z2_elev_med = z2_stats[5]
    z2_elev_p16 = z2_stats[11]
    z2_elev_p84 = z2_stats[12]

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
    rho_f = 0.60
    rho_is = 0.85
    #Can estimate ELA values computed from hypsometry and typical AAR
    #For now, assume ELA is mean
    z1_ela = None
    z1_ela = z1_stats[3]
    z2_ela = z2_stats[3]
    #Note: in theory, the ELA should get higher with mass loss
    #In practice, using mean and same polygon, ELA gets lower as glacier surface thins
    print("ELA(t1): %0.1f" % z1_ela)
    print("ELA(t2): %0.1f" % z2_ela)
    if z1_ela is None:
        mb = dhdt * rho_is
    else:
        #Initiate with average density
        mb = dhdt*(rho_is + rho_f)/2.
        #Everything that is above ELA at t2 is elevation change over firn, use firn density
        accum_mask = (z2 > z2_ela).filled(0).astype(bool)
        mb[accum_mask] = (dhdt*rho_f)[accum_mask]
        #Everything that is below ELA at t1 is elevation change over ice, use ice density
        abl_mask = (z1 <= z1_ela).filled(0).astype(bool)
        mb[abl_mask] = (dhdt*rho_is)[abl_mask]
        #Everything in between, use average of ice and firn density
        #mb[(z1 > z1_ela) || (z2 <= z2_ela)] = dhdt*(rhois + rho_f)/2.
        #Linear ramp
        #rho_f + z2*((rho_is - rho_f)/(z2_ela - z1_ela))
        #mb = np.where(dhdt < ela, dhdt*rho_i, dhdt*rho_s)

    mb_stats = malib.print_stats(mb)
    mb_mean = mb_stats[3]
    dmbdt_total_myr = mb_mean*glac_area
    mb_sum = np.sum(mb)*ds_res[0]*ds_res[1]

    prism_ppt = np.ma.array(iolib.ds_getma(ds_list[2]), mask=glac_geom_mask)/1000.
    prism_ppt_stats = malib.print_stats(prism_ppt)
    prism_ppt_mean = prism_ppt_stats[3]
    prism_tmean = np.ma.array(iolib.ds_getma(ds_list[3]), mask=glac_geom_mask)
    prism_tmean_stats = malib.print_stats(prism_tmean)
    prism_tmean_mean = prism_tmean_stats[3]

    print('Mean mb: %0.2f mwe/yr' % mb_mean)
    print('Sum/Area mb: %0.2f mwe/yr' % (mb_sum/glac_area))
    print('Mean mb * Area: %0.2f mwe/yr' % dmbdt_total_myr)
    print('Sum mb: %0.2f mwe/yr' % mb_sum)
    print('-------------------------------')
    
    #out.append([cx, cy, mb_mean, (glac_area/1E6), t1.mean(), dt.mean()])
    #out.append([glacnum, cx, cy, z2_elev_med, z2_elev_p16, z2_elev_p84, mb_mean, (glac_area/1E6), t1.mean(), dt.mean()])
    out.append([glacnum, cx, cy, z2_elev_med, z2_elev_p16, z2_elev_p84, mb_mean, (glac_area/1E6), t1.mean(), dt.mean(), \
            prism_ppt_mean, prism_tmean_mean])

    writeout=False
    if writeout:
        out_dz_fn = os.path.join(outdir, feat_fn+'_dz.tif')
        #dz_ds = iolib.gtif_drv.CreateCopy(dz_fn, ds_list[0], options=iolib.gdal_opt)
        iolib.writeGTiff(dz, out_dz_fn, ds_list[0])

        out_z1_fn = os.path.join(outdir, feat_fn+'_z1.tif')
        iolib.writeGTiff(z1, out_z1_fn, ds_list[0])

        out_z2_fn = os.path.join(outdir, feat_fn+'_z2.tif')
        iolib.writeGTiff(z2, out_z2_fn, ds_list[0])

        out_z1_date_fn = os.path.join(outdir, feat_fn+'_ned_date.tif')
        iolib.writeGTiff(z1_date, out_z1_date_fn, ds_list[0])

        out_prism_ppt_fn = os.path.join(outdir, feat_fn+'_precip.tif')
        iolib.writeGTiff(prism_ppt, out_prism_ppt_fn, ds_list[0])

        out_prism_tmean_fn = os.path.join(outdir, feat_fn+'_tmean.tif')
        iolib.writeGTiff(prism_tmean, out_prism_tmean_fn, ds_list[0])

    add_fields = False
    if add_fields:
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
#Sort by area
out = out[out[:,3].argsort()[::-1]]
#out_fn = 'conus_mb_summer2014-2016_20170205.csv'
#out_fn = 'conus_mb_summer2014-2016_20170215.csv'
out_fn = 'conus_mb_summer2014-2016_ela_20170216.csv'
#out_header = 'x,y,mb_mwea,area_km2,t1,dt'
#out_header = 'glacnum,x,y,z_med,z_p16,z_p84,mb_mwea,area_km2,t1,dt'
out_header = 'glacnum,x,y,z_med,z_p16,z_p84,mb_mwea,area_km2,t1,dt,precip_mwe,temp'
np.savetxt(out_fn, out, fmt='%0.2f', delimiter=',', header=out_header)

#Write out new shp with features containing stats
#One shp preserves input features, regardless of source date
#Another shp splits glacier poly based on NED source date
