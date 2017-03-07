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
from pygeotools.lib import timelib

from imview.lib import pltlib

from datetime import datetime, timedelta
import time

def get_bins(dem, bin_width=100.0):
    #Define min and max elevation
    minz, maxz= list(malib.calcperc(dem, perc=(0.01, 99.99)))
    minz = np.floor(minz/bin_width) * bin_width
    maxz = np.ceil(maxz/bin_width) * bin_width
    #Compute bin edges and centers
    bin_edges = np.arange(minz, maxz + bin_width, bin_width)
    bin_centers = bin_edges[:-1] + np.diff(bin_edges)/2.0
    return bin_edges, bin_centers

topdir='/nobackup/deshean'
site='conus'
#site='hma'

#topdir='/Volumes/SHEAN_1TB_SSD/site_poly_highcount_rect3_rerun/rainier'
#site='rainier'

writeout=True
ts = datetime.now().strftime('%Y%m%d_%H%M')
out_fn = '%s_mb_%s.csv' % (site, ts)

if site == 'conus':
    #'+proj=aea +lat_1=36 +lat_2=49 +lat_0=43 +lon_0=-115 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs '
    #outdir = os.path.join(topdir,'conus/dem2/conus_32611_8m/glac_poly_out')
    outdir = os.path.join(topdir,'conus/dem2/glac_poly_out')
    aea_srs = geolib.conus_aea_srs

    #Glacier shp
    #ogr2ogr -t_srs '+proj=aea +lat_1=36 +lat_2=49 +lat_0=43 +lon_0=-115 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ' 24k_selection_aea.shp 24k_selection_32610.shp
    #glac_shp_fn = '/nobackupp8/deshean/conus/shp/24k_selection_aea.shp'
    #This has already been filtered by area
    glac_shp_fn = os.path.join(topdir,'conus/shp/24k_selection_aea_min0.1km2.shp')

    #Raster difference map between NED and WV mosaic
    z1_fn = os.path.join(topdir,'rpcdem/ned1_2003/ned1_2003_adj.vrt')
    #NED 2003 dates
    z1_date_shp_fn = os.path.join(topdir,'rpcdem/ned1_2003/meta0306_PAL_24k_10kmbuffer_clean_dissolve_aea.shp')
    #ogr2ogr -t_srs '+proj=aea +lat_1=36 +lat_2=49 +lat_0=43 +lon_0=-115 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ' meta0306_PAL_24k_10kmbuffer_clean_dissolve_aea.shp meta0306_PAL_24k_10kmbuffer_clean_dissolve_32611.shp
    z1_date_shp_ds = ogr.Open(z1_date_shp_fn)
    z1_date_shp_lyr = z1_date_shp_ds.GetLayer()
    z1_date_shp_srs = z1_date_shp_lyr.GetSpatialRef()
    z1_date_shp_lyr.ResetReading()

    z2_fn = os.path.join(topdir,'conus/dem2/conus_8m_tile_coreg_round3_summer2014-2016/conus_8m_tile_coreg_round3_summer2014-2016.vrt')
    z2_date = datetime(2015, 1, 1)
    z2_date = 2015.0
    #PRISM climate data, 800-m 
    prism_ppt_annual_fn = os.path.join(topdir,'conus/prism/normals/annual/ppt/PRISM_ppt_30yr_normal_800mM2_annual_bil.bil')
    prism_tmean_annual_fn = os.path.join(topdir,'conus/prism/normals/annual/tmean/PRISM_tmean_30yr_normal_800mM2_annual_bil.bil')
    prism_ppt_summer_fn = os.path.join(topdir,'conus/prism/normals/monthly/PRISM_ppt_30yr_normal_800mM2_06-09_summer_mean.tif')
    prism_ppt_winter_fn = os.path.join(topdir,'conus/prism/normals/monthly/PRISM_ppt_30yr_normal_800mM2_10-05_winter_mean.tif')
    prism_tmean_summer_fn = os.path.join(topdir,'conus/prism/normals/monthly/PRISM_tmean_30yr_normal_800mM2_06-09_summer_mean.tif')
    prism_tmean_winter_fn = os.path.join(topdir,'conus/prism/normals/monthly/PRISM_tmean_30yr_normal_800mM2_10-05_winter_mean.tif')

    #Noisy Creek, Silver
    glacier_dict = {}
    glacier_dict[6012] = 'EmmonsGlacier'
    glacier_dict[6096] = 'Nisqually-WilsonGlacier'
    glacier_dict[10480] = 'SouthCascadeGlacier'
    #Note: Sandalee has 3 records, 2693, 2695, 2696
    glacier_dict[2696] = 'SandaleeGlacier'
    glacier_dict[3070] = 'NorthKlawattiGlacier'
    glacier_dict[1969] = 'NoisyCreekGlacier'
    glacier_dict[2294] = 'SilverGlacier'
    glacier_dict[2500] = 'EastonGlacier'
    glacier_dict[5510] = 'BlueGlacier'
    glacier_dict[5603] = 'EelGlacier'
    glacier_dict[1589] = 'SperryGlacier'
    glacier_dict[1277] = 'GrinnellGlacier'
    glacier_dict[10490] = 'ConnessGlacier'
    glacier_dict[10071] = 'WheelerGlacier'
    glacier_dict[9129] = 'LyellGlacier'
    glacier_dict[9130] = 'LyellGlacier'

elif site == 'hma':
    #'+proj=aea +lat_1=25 +lat_2=47 +lat_0=36 +lon_0=85 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs '
    outdir = os.path.join(topdir,'hma/hma1_2016dec22/glac_poly_out')
    aea_srs = geolib.conus_aea_srs
    glac_shp_fn = os.path.join(topdir,'data/rgi50/regions/rgi50_hma_aea.shp')
    #SRTM
    z1_fn = os.path.join(topdir,'rpcdem/hma/srtm1/hma_srtm_gl1.vrt')
    z1_date = timelib.dt2decyear(datetime(2000,2,11))
    #z2_fn = '/nobackup/deshean/hma/hma1_2016dec22/hma_8m_tile/hma_8m.vrt'
    #z2_fn = os.path.join(topdir,'hma/hma1_2016dec22/hma_8m_tile/hma_8m.vrt')
    z2_fn = os.path.join(topdir,'hma/hma1_2016dec22/hma_8m_tile_round2_20170220/hma_8m_round2.vrt')
    z2_date = datetime(2015, 1, 1)
    z2_date = 2015.0
elif site == 'rainier':
    outdir = os.path.join(topdir,'mb')
    aea_srs = geolib.conus_aea_srs
    glac_shp_fn = '/Users/dshean/data/conus_glacierpoly_24k/rainier_24k_1970-2015_mb_aea.shp'
    z1_fn = sys.argv[1]
    z1_date = timelib.mean_date(timelib.fn_getdatetime_list(z1_fn))
    z2_fn = sys.argv[2]
    z2_date = timelib.mean_date(timelib.fn_getdatetime_list(z2_fn))
else:
    sys.exit()

#Filter glacier poly - let's stick with big glaciers for now
min_glac_area = 0.1 #km^2

#List to hold output
out = []

if '24k' in glac_shp_fn: 
    glacname_fieldname = "GLACNAME"
    glacnum_fieldname = "GLACNUM"
elif 'rgi' in glac_shp_fn:
    #Use RGI
    glacname_fieldname = "Name"
    #RGIId (String) = RGI50-01.00004
    glacnum_fieldname = "RGIId"
else:
    sys.exit('Unrecognized glacier shp filename')

glac_shp_ds = ogr.Open(glac_shp_fn, 0)
glac_shp_lyr = glac_shp_ds.GetLayer()
glac_shp_srs = glac_shp_lyr.GetSpatialRef()
feat_count = glac_shp_lyr.GetFeatureCount()
print("Input glacier polygon count: %i" % feat_count)

#Spatial filter
z1_ds = gdal.Open(z1_fn)
z2_ds = gdal.Open(z2_fn)
dz_int_geom = geolib.ds_geom_intersection([z1_ds, z2_ds], t_srs=glac_shp_srs)
glac_shp_lyr.SetSpatialFilter(dz_int_geom)
feat_count = glac_shp_lyr.GetFeatureCount()
print("Filtered glacier polygon count: %i" % feat_count)
glac_shp_lyr.ResetReading()

print("Processing %i features" % feat_count)

if not os.path.exists(outdir):  
    os.makedirs(outdir)

#Update shp
#out_glac_shp_fn = os.path.join(outdir, os.path.splitext(glac_shp_fn)[0]+'_mb.shp')
#field_defn = ogr.FieldDefn("mb_mwe", ogr.OFTReal)
#glac_shp_lyr.CreateField(field_defn)

for n, feat in enumerate(glac_shp_lyr):

    glacname = feat.GetField(glacname_fieldname)

    if glacname is None:
        glacname = ""
    else:
        glacname = glacname.replace(" ", "")
        glacname = glacname.replace("_", "")

    glacnum = feat.GetField(glacnum_fieldname)
    if '24k' in glac_shp_fn:
        glacnum = int(glacnum)
    else:
        #RGIId (String) = RGI50-01.00004
        glacnum = float(glacnum.split('-')[-1])*1000000

    #if glacname != "EmmonsGlacier":
    #if glacname != "Nisqually-WilsonGlacier":
    #if glacname != "SouthCascadeGlacier":
    if glacnum not in glacier_dict.keys():
        continue
    else:
        glacname = glacier_dict[glacnum]

    feat_fn = "%s_%s" % (glacnum, glacname)
    print("\n%i of %i: %s\n" % (n+1, feat_count, feat_fn))
    glac_geom = feat.GetGeometryRef()
    glac_geom.AssignSpatialReference(glac_shp_srs)
    glac_geom_extent = geolib.geom_extent(glac_geom)
    glac_area = glac_geom.GetArea()
    if glac_area/1E6 < min_glac_area:
        print("Glacier area below %0.1f km2 threshold" % min_glac_area)
        continue

    cx, cy = glac_geom.Centroid().GetPoint_2D()

    #Warp everything to common res/extent/proj
    ds_list = warplib.memwarp_multi_fn([z1_fn, z2_fn], res='min', \
            extent=glac_geom_extent, t_srs=aea_srs, verbose=False)

    if site == 'conus':
        #Add prism datasets
        prism_fn_list = [prism_ppt_annual_fn, prism_tmean_annual_fn]
        prism_fn_list.extend([prism_ppt_summer_fn, prism_ppt_winter_fn, prism_tmean_summer_fn, prism_tmean_winter_fn])
        ds_list.extend(warplib.memwarp_multi_fn(prism_fn_list, res=ds_list[0], extent=glac_geom_extent, t_srs=aea_srs, verbose=False))

    #Check to see if z2 is empty, as z1 should be continuous
    z2 = iolib.ds_getma(ds_list[1])
    if z2.count() == 0:
        print("No z2 pixels")
        continue

    glac_geom_mask = geolib.geom2mask(glac_geom, ds_list[0])
    z1 = np.ma.array(iolib.ds_getma(ds_list[0]), mask=glac_geom_mask)
    z2 = np.ma.array(z2, mask=glac_geom_mask)
    dz = z2 - z1
    if dz.count() == 0:
        print("No valid dz pixels")
        continue

    filter_outliers = False
    #Remove clearly bogus pixels
    if filter_outliers:
        bad_perc = (0.1, 99.9)
        rangelim = malib.calcperc(dz, bad_perc)
        dz = np.ma.masked_outside(dz, *rangelim)

    ds_res = geolib.get_res(ds_list[0])
    valid_area = dz.count()*ds_res[0]*ds_res[1]
    valid_area_perc = valid_area/glac_area
    min_valid_area_perc = 0.80
    if valid_area_perc < min_valid_area_perc:
        print("Not enough valid pixels. %0.1f%% percent of glacier polygon area" % (100*valid_area_perc))
        continue

    #Rasterize NED source dates
    if site == 'conus':
        z1_date_r_ds = iolib.mem_drv.CreateCopy('', ds_list[0])
        gdal.RasterizeLayer(z1_date_r_ds, [1], z1_date_shp_lyr, options=["ATTRIBUTE=S_DATE_CLN"])
        z1_date = np.ma.array(iolib.ds_getma(z1_date_r_ds), mask=glac_geom_mask)

    #Filter dz - throw out abs differences >150 m

    #Compute dz, volume change, mass balance and stats
    z1_stats = malib.print_stats(z1)
    z2_stats = malib.print_stats(z2)
    z2_elev_med = z2_stats[5]
    z2_elev_p16 = z2_stats[11]
    z2_elev_p84 = z2_stats[12]

    #These can be timestamp arrays or datetime objects
    t1 = z1_date
    t2 = z2_date
    #This is decimal years
    dt = t2 - t1
    if isinstance(dt, timedelta):
        dt = dt.total_seconds()/timelib.spy
    #m/yr
    dhdt = dz/dt
    #dhdt_stats = malib.print_stats(dhdt)
    #dhdt_mean = dhdt_stats[3]

    #Output mean values for timestamp arrays
    if site == 'conus':
        t1 = t1.mean()
        dt = dt.mean()

    if isinstance(t1, datetime):
        t1 = timelib.dt2decyear(t1)

    if isinstance(t2, datetime):
        t2 = timelib.dt2decyear(t2)

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

    if z1_ela > z2_ela:
        min_ela = z2_ela
        max_ela = z1_ela
    else:
        min_ela = z1_ela
        max_ela = z2_ela

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

    outlist = [glacnum, cx, cy, z2_elev_med, z2_elev_p16, z2_elev_p84, mb_mean, (glac_area/1E6), t1, t2, dt]

    if site == 'conus':
        prism_ppt_annual = np.ma.array(iolib.ds_getma(ds_list[2]), mask=glac_geom_mask)/1000.
        prism_ppt_annual_stats = malib.print_stats(prism_ppt_annual)
        prism_ppt_annual_mean = prism_ppt_annual_stats[3]

        prism_tmean_annual = np.ma.array(iolib.ds_getma(ds_list[3]), mask=glac_geom_mask)
        prism_tmean_annual_stats = malib.print_stats(prism_tmean_annual)
        prism_tmean_annual_mean = prism_tmean_annual_stats[3]

        outlist.extend([prism_ppt_annual_mean, prism_tmean_annual_mean])

        #This is mean monthly summer precip, need to multiply by nmonths to get cumulative
        n_summer = 4
        prism_ppt_summer = n_summer * np.ma.array(iolib.ds_getma(ds_list[4]), mask=glac_geom_mask)/1000.
        prism_ppt_summer_stats = malib.print_stats(prism_ppt_summer)
        prism_ppt_summer_mean = prism_ppt_summer_stats[3]

        n_winter = 8
        prism_ppt_winter = n_winter * np.ma.array(iolib.ds_getma(ds_list[5]), mask=glac_geom_mask)/1000.
        prism_ppt_winter_stats = malib.print_stats(prism_ppt_winter)
        prism_ppt_winter_mean = prism_ppt_winter_stats[3]

        prism_tmean_summer = np.ma.array(iolib.ds_getma(ds_list[6]), mask=glac_geom_mask)
        prism_tmean_summer_stats = malib.print_stats(prism_tmean_summer)
        prism_tmean_summer_mean = prism_tmean_summer_stats[3]

        prism_tmean_winter = np.ma.array(iolib.ds_getma(ds_list[7]), mask=glac_geom_mask)
        prism_tmean_winter_stats = malib.print_stats(prism_tmean_winter)
        prism_tmean_winter_mean = prism_tmean_winter_stats[3]

        outlist.extend([prism_ppt_summer_mean, prism_ppt_winter_mean, prism_tmean_summer_mean, prism_tmean_winter_mean])

    print('Mean mb: %0.2f mwe/yr' % mb_mean)
    print('Sum/Area mb: %0.2f mwe/yr' % (mb_sum/glac_area))
    print('Mean mb * Area: %0.2f mwe/yr' % dmbdt_total_myr)
    print('Sum mb: %0.2f mwe/yr' % mb_sum)
    print('-------------------------------')

    #Write to master list
    out.append(outlist)

    if writeout:
        out_dz_fn = os.path.join(outdir, feat_fn+'_dz.tif')
        iolib.writeGTiff(dz, out_dz_fn, ds_list[0])

        out_z1_fn = os.path.join(outdir, feat_fn+'_z1.tif')
        iolib.writeGTiff(z1, out_z1_fn, ds_list[0])

        out_z2_fn = os.path.join(outdir, feat_fn+'_z2.tif')
        iolib.writeGTiff(z2, out_z2_fn, ds_list[0])

        if site == 'conus':
            out_z1_date_fn = os.path.join(outdir, feat_fn+'_ned_date.tif')
            iolib.writeGTiff(z1_date, out_z1_date_fn, ds_list[0])

            out_prism_ppt_annual_fn = os.path.join(outdir, feat_fn+'_precip_annual.tif')
            iolib.writeGTiff(prism_ppt_annual, out_prism_ppt_annual_fn, ds_list[0])
            out_prism_tmean_annual_fn = os.path.join(outdir, feat_fn+'_tmean_annual.tif')
            iolib.writeGTiff(prism_tmean_annual, out_prism_tmean_annual_fn, ds_list[0])

            out_prism_ppt_summer_fn = os.path.join(outdir, feat_fn+'_precip_summer.tif')
            iolib.writeGTiff(prism_ppt_summer, out_prism_ppt_summer_fn, ds_list[0])
            out_prism_ppt_winter_fn = os.path.join(outdir, feat_fn+'_precip_winter.tif')
            iolib.writeGTiff(prism_ppt_winter, out_prism_ppt_winter_fn, ds_list[0])

            out_prism_tmean_summer_fn = os.path.join(outdir, feat_fn+'_tmean_summer.tif')
            iolib.writeGTiff(prism_tmean_summer, out_prism_tmean_summer_fn, ds_list[0])
            out_prism_tmean_winter_fn = os.path.join(outdir, feat_fn+'_tmean_winter.tif')
            iolib.writeGTiff(prism_tmean_winter, out_prism_tmean_winter_fn, ds_list[0])

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

    mb_plot = True
    if mb_plot:
        print("Generating histograms")
        z_bin_edges, z_bin_centers = get_bins(z1, 10.)
        z1_bin_counts, z1_bin_edges = np.histogram(z1, bins=z_bin_edges)
        z1_bin_areas = z1_bin_counts * ds_res[0] * ds_res[1] / 1E6
        z2_bin_counts, z2_bin_edges = np.histogram(z2, bins=z_bin_edges)
        z2_bin_areas = z2_bin_counts * ds_res[0] * ds_res[1] / 1E6

        #dz_bin_edges, dz_bin_centers = get_bins(dz, 1.)
        #dz_bin_counts, dz_bin_edges = np.histogram(dz, bins=dz_bin_edges)
        #dz_bin_areas = dz_bin_counts * ds_res * ds_res / 1E6
        dz_bin_med = np.ma.masked_all_like(z1_bin_areas)
        dz_bin_mad = np.ma.masked_all_like(z1_bin_areas)
        idx = np.digitize(z1, z_bin_edges)
        for n in range(z_bin_centers.size):
            dz_bin_samp = mb[(idx == n+1)]
            #dz_bin_samp = dhdt[(idx == n+1)]
            dz_bin_med[n] = malib.fast_median(dz_bin_samp)
            dz_bin_mad[n] = malib.mad(dz_bin_samp)
            dz_bin_med[n] = dz_bin_samp.mean()
            dz_bin_mad[n] = dz_bin_samp.std()

        print("Generating plot")
        f,axa = plt.subplots(1,3, figsize=(10,7.5))
        f.suptitle(feat_fn)
        alpha = 1.0
        hs = True
        if hs:
            z1_hs = geolib.gdaldem_wrapper(out_z1_fn, product='hs', returnma=True)
            z2_hs = geolib.gdaldem_wrapper(out_z2_fn, product='hs', returnma=True)
            hs_clim = malib.calcperc(z2_hs, (2,98))
            z1_hs_im = axa[0].imshow(z1_hs, cmap='gray', vmin=hs_clim[0], vmax=hs_clim[1])
            z2_hs_im = axa[1].imshow(z2_hs, cmap='gray', vmin=hs_clim[0], vmax=hs_clim[1])
            alpha = 0.5
        z1_im = axa[0].imshow(z1, cmap='cpt_rainbow', vmin=z_bin_edges[0], vmax=z_bin_edges[-1], alpha=alpha)
        z2_im = axa[1].imshow(z2, cmap='cpt_rainbow', vmin=z_bin_edges[0], vmax=z_bin_edges[-1], alpha=alpha)
        axa[0].contour(z1, [z1_ela,], linewidths=0.5, linestyles=':', colors='w')
        axa[1].contour(z2, [z2_ela,], linewidths=0.5, linestyles=':', colors='w')
        t1_title = int(np.round(t1))
        t2_title = int(np.round(t2))
        axa[0].set_title(t1_title)
        axa[1].set_title(t2_title)
        axa[2].set_title('%i to %i' % (t1_title, t2_title))
        #dz_clim = (-10, 10)
        dz_clim = (-2.0, 2.0)
        dz_im = axa[2].imshow(dhdt, cmap='RdBu', vmin=dz_clim[0], vmax=dz_clim[1])
        for ax in axa:
            pltlib.hide_ticks(ax)
            ax.set_facecolor('k')
        sb_loc = pltlib.best_scalebar_location(z1)
        pltlib.add_scalebar(axa[0], ds_res[0], location=sb_loc)
        pltlib.add_cbar(axa[0], z1_im, label='Elevation (m WGS84)')
        pltlib.add_cbar(axa[1], z2_im, label='Elevation (m WGS84)')
        pltlib.add_cbar(axa[2], dz_im, label='dh/dt (m/yr)')
        plt.tight_layout()
        #Make room for suptitle
        plt.subplots_adjust(top=0.90)
        fig_fn = os.path.join(outdir, feat_fn+'_mb_map.png')
        plt.savefig(fig_fn, dpi=300)

        f,axa = plt.subplots(1,2, figsize=(6, 6))
        f.suptitle(feat_fn)
        axa[0].plot(z1_bin_areas, z_bin_centers, label='%0.2f' % t1)
        axa[0].plot(z2_bin_areas, z_bin_centers, label='%0.2f' % t2)
        axa[0].axhline(z1_ela, ls=':', c='C0')
        axa[0].axhline(z2_ela, ls=':', c='C1')
        axa[0].legend(prop={'size':8}, loc='upper right')
        axa[0].set_ylabel('Elevation (m WGS84)')
        axa[0].set_xlabel('Area $\mathregular{km^2}$')
        axa[0].minorticks_on()
        axa[1].yaxis.tick_right()
        axa[1].yaxis.set_label_position("right")
        axa[1].axvline(0, lw=1.0, c='k')
        axa[1].axvline(mb_mean, lw=0.5, ls=':', c='k', label='%0.2f m w.e./yr' % mb_mean)
        axa[1].legend(prop={'size':8}, loc='upper right')
        axa[1].plot(dz_bin_med, z_bin_centers, color='k')
        axa[1].fill_betweenx(z_bin_centers, 0, dz_bin_med, where=(dz_bin_med<0), color='r', alpha=0.2)
        axa[1].fill_betweenx(z_bin_centers, 0, dz_bin_med, where=(dz_bin_med>0), color='b', alpha=0.2)
        #axa[1].set_ylabel('Elevation (m WGS84)')
        #axa[1].set_xlabel('dh/dt (m/yr)')
        axa[1].set_xlabel('mb (m w.e./yr)')
        axa[1].minorticks_on()
        axa[1].set_xlim(-2.0, 2.0)
        plt.tight_layout()
        #Make room for suptitle
        plt.subplots_adjust(top=0.95)
        fig_fn = os.path.join(outdir, feat_fn+'_mb_aed.png')
        plt.savefig(fig_fn, dpi=300)
    
#plt.show()

glac_shp_ds = None

out = np.array(out)
#Sort by area
out = out[out[:,3].argsort()[::-1]]

ts = datetime.now().strftime('%Y%m%d_%H%M')
out_fn = '%s_mb_%s.csv' % (site, ts)
out_fn = os.path.join(outdir, out_fn)

out_header = 'glacnum,x,y,z_med,z_p16,z_p84,mb_mwea,area_km2,t1,t2,dt'
if site == 'conus':
    out_header += ',ppt_a,tmean_a'
    out_header += ',ppt_s,ppt_w,tmean_s,tmean_w'

np.savetxt(out_fn, out, fmt='%0.2f', delimiter=',', header=out_header)

#Write out new shp with features containing stats
#One shp preserves input features, regardless of source date
#Another shp splits glacier poly based on NED source date
