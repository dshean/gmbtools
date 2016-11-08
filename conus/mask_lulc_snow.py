#! /usr/bin/env python

#David Shean
#dshean@gmail.com

#Commands and notes to automate control surface identification using LULC and snowcover grids
#Note: should combine with pc_align_wrapper and apply_mask functionality

#Pregenerate toa
#parallel --delay 3 -j 12 '~/src/demtools/toa.sh {}' ::: {WV,GE}*
#for i in */dem*/*32m.tif; do ~/src/demtools/conus/mask_lulc_snow.py $i; done
#Don't use parallel, as modscag fn could be overwrittent
##parallel --delay 3 -j 12 '~/src/demtools/conus/mask_lulc_snow.py {}' ::: */dem*/*32m.tif

import sys
import os
import subprocess
import glob

import gdal
import osr
import numpy as np

from datetime import datetime, timedelta

from lib import malib
from lib import warplib
from lib import geolib
from lib import timelib

topdir = os.getcwd()

nlcd_dir='/nobackup/deshean/rpcdem/nlcd'
if not os.path.exists(nlcd_dir):
    nlcd_dir='/scr/nlcd_2011_landcover_2011_edition_2014_10_10'

#Original has nan as ndv
#nlcd_fn = 'nlcd_2011_landcover_2011_edition_2014_10_10_pnw_clip.tif'
#Repaired ndv
#nlcd_fn = 'nlcd_2011_landcover_2011_edition_2014_10_10_pnw_clip_ndv.tif'
#Download http://www.landfire.gov/bulk/downloadfile.php?TYPE=nlcd2011&FNAME=nlcd_2011_landcover_2011_edition_2014_10_10.zip
#unzip fails to extract ige file
#7za x nlcd_2011_landcover_2011_edition_2014_10_10.zip
nlcd_fn = os.path.join(nlcd_dir, 'nlcd_2011_landcover_2011_edition_2014_10_10_warp.tif')

#Downloaded from http://www.glims.org/RGI/rgi50_files/02_rgi50_WesternCanadaUS.zip
#ogr2ogr -t_srs EPSG:32610 02_rgi50_WesternCanadaUS_32610.shp 02_rgi50_WesternCanadaUS.shp
#Manually clipped to lower48 in QGIS
#glac_shp_fn = '/Volumes/500GB_1/RGI/02_rgi50_WesternCanadaUS/02_rgi50_WesternCanadaUS_32610_lower48.shp'
#glac_shp_fn = '/Volumes/500GB_1/RGI/24kgrteqaul0.01km/24k_selection_32610.shp'
glac_shp_fn = os.path.join(nlcd_dir, '24k_selection_32610.shp')
if not os.path.exists(glac_shp_fn):
    glac_shp_fn = '/Volumes/500GB_1/RGI/24kgrteqaul0.01km/24k_selection_32610.shp'

#Create rockmask from LULC and remove glaciers
#Alternatively, can create vector mask
#gdal_polygonize.py rockmask.tif -f 'ESRI Shapefile' rockmask.shp
def rockmask(nlcd_fn, outmask_fn=os.path.join(nlcd_dir,'rockmask_conus.tif'), glac_shp_fn=None):

    print "Creating rockmask"
    print "Loading LULC"
    ds = gdal.Open(nlcd_fn)
    b = ds.GetRasterBand(1)
    l = b.ReadAsArray()

    print "Isolating rock"
    #Isolate rock

    #2011 LULC grids, 30 m
    #http://www.mrlc.gov/nlcd11_leg.php

    #12 - ice
    #31 - rock
    #11 - open water
    #52 - shrub, <5 m tall, >20%
    #42 - evergreeen forest

    mask = (l==31)
    #This includes river valleys
    #mask = ((l==31) or (l==11))

    #Remove glacier and perennial snowfields 
    if glac_shp_fn:
        print "Masking glaciers"
        icemask = geolib.shp2array(glac_shp_fn, r_ds=ds)
        mask *= icemask

    #This writes out 1 for rock, 0 for everything else (ndv)
    print "Writing out"
    malib.writeGTiff(mask, outmask_fn, src_ds=ds)

#This works with ftp
def getfile(url, outdir=None):
    fn = os.path.split(url)[-1]
    if outdir is not None:
        fn = os.path.join(outdir, fn)
    if not os.path.exists(fn):
        import urllib
        print "Retrieving: %s" % url
        #Add progress bar
        urllib.urlretrieve(url, fn)
    return fn

def getfile2(url, auth=None, outdir=None):
    import requests

    print "Retrieving: %s" % url
    #url = "http://download.thinkbroadband.com/10MB.zip"
    fn = os.path.split(url)[-1]
    if outdir is not None:
        fn = os.path.join(outdir, fn)
    if auth is not None:
        r = requests.get(url, stream=True, auth=auth)
    else:
        r = requests.get(url, stream=True)

    #from tqdm import tqdm
    #with open(fn, "wb") as handle:
    #    for data in tqdm(r.iter_content()):
    #        handle.write(data)
    chunk_size = 1000000
    with open(fn, 'wb') as fd:
        for chunk in r.iter_content(chunk_size):
            fd.write(chunk)

#Should use global MODIS 500 m snowcover grids, 8 day
#h09v04 should cover WA
#http://nsidc.org/data/docs/daac/modis_v5/mod10a2_modis_terra_snow_8-day_global_500m_grid.gd.html
#HDF4, sinusoidal
#Should be able to load up with warplib without issue

#Alternatively, use SNODAS data
#http://nsidc.org/data/docs/noaa/g02158_snodas_snow_cover_model/index.html
#1036 is snow depth
#Need to delete 'Created by module comment' line from Hdr, can contain too many characters
#snodas_fn = sys.argv[2]
#snodas_fn = 'us_ssmv11036tS__T0001TTNATS2015042205HP001.Hdr'

#Note: these are only available frmo 2010-present
def get_snodas(dem_dt, outdir=None):
    import tarfile
    import gzip

    masked = False
    if dem_dt < datetime(2010,1,1):
        masked = True

    if masked:
        snodas_url_str = 'ftp://sidads.colorado.edu/DATASETS/NOAA/G02158/masked/%Y/%m_%b/SNODAS_%Y%m%d.tar'
        tar_subfn_str_fmt = 'us_ssmv11036tS__T0001TTNATS%%Y%%m%%d05HP001.%s.gz'
    else:
        #snodas_url = 'ftp://sidads.colorado.edu/DATASETS/NOAA/G02158/unmasked/2016/07_Jul/SNODAS_unmasked_20160701.tar'
        snodas_url_str = 'ftp://sidads.colorado.edu/DATASETS/NOAA/G02158/unmasked/%Y/%m_%b/SNODAS_unmasked_%Y%m%d.tar'
        #'./zz_ssmv11036tS__T0001TTNATS2015041705HP001.Hdr.gz'
        tar_subfn_str_fmt = './zz_ssmv11036tS__T0001TTNATS%%Y%%m%%d05HP001.%s.gz'

    snodas_url = dem_dt.strftime(snodas_url_str)
    snodas_tar_fn = getfile(snodas_url, outdir=outdir)
    print "Unpacking"
    #gunzip extract both dat and Hdr files, tar.gz
    tar = tarfile.open(snodas_tar_fn)
    for ext in ('dat', 'Hdr'):
        tar_subfn_str = tar_subfn_str_fmt % ext
        tar_subfn_gz = dem_dt.strftime(tar_subfn_str)
        tar_subfn = os.path.splitext(tar_subfn_gz)[0]
        print tar_subfn
        if outdir is not None:
            tar_subfn = os.path.join(outdir, tar_subfn)
        if not os.path.exists(tar_subfn):
            #Should be able to do this without writing intermediate gz to disk
            #f = tar.extractfile(tar_subfn)
            tar.extract(tar_subfn_gz)
            with gzip.open(tar_subfn_gz, 'rb') as f:
                outf = open(tar_subfn, 'wb')
                outf.write(f.read())
                outf.close()
            os.remove(tar_subfn_gz)

    #Need to delete 'Created by module comment' line from Hdr, can contain too many characters
    bad_str = 'Created by module comment'
    snodas_fn = tar_subfn
    f = open(snodas_fn)
    output = []
    for line in f:
        if not bad_str in line:
            output.append(line)
    f.close()
    f = open(snodas_fn, 'w')
    f.writelines(output)
    f.close()

    snodas_ds = gdal.Open(snodas_fn)
    return snodas_ds

def get_modscag(dem_dt, outdir=None, pad_days=timedelta(days=7)):
    #https://snow-data.jpl.nasa.gov/modscag-historic/2015/001/MOD09GA.A2015001.h07v03.005.2015006001833.snow_fraction.tif
    tile_list = ('h08v04', 'h09v04', 'h10v04', 'h08v05', 'h09v05')
    uname = 'dshean'
    pw = 'yGmFbH_gDeSAwZ3x'
    #wget -A'h8v4*snow_fraction.tif' --user=dshean --password=yGmFbH_gDeSAwZ3x
    import re
    import requests
    from requests.auth import HTTPDigestAuth
    from bs4 import BeautifulSoup
    auth = HTTPDigestAuth(uname, pw)
    modscag_url_str = 'https://snow-data.jpl.nasa.gov/modscag-historic/%Y/%j/' 
    dt_list = timelib.dt_range(dem_dt-pad_days, dem_dt+pad_days+timedelta(1), timedelta(1))
    out_vrt_fn_list = []
    for dt in dt_list:
        out_vrt_fn = os.path.join(outdir, dt.strftime('%Y%m%d_snow_fraction.vrt'))
        out_tif_fn = os.path.splitext(out_vrt_fn)[0]+'.tif'
        if os.path.exists(out_vrt_fn):
            out_vrt_fn_list.append(out_vrt_fn)
        else:
            modscag_url_base = dt.strftime(modscag_url_str)
            print modscag_url_base
            r = requests.get(modscag_url_base, auth=auth)
            if not r.ok:
                #Try to use real-time products
                modscag_url_str = 'https://snow-data.jpl.nasa.gov/modscag/%Y/%j/' 
                modscag_url_base = dt.strftime(modscag_url_str)
                print modscag_url_base
                r = requests.get(modscag_url_base, auth=auth)
            if r.ok:
                parsed_html = BeautifulSoup(r.content, "html.parser")
                modscag_fn_list = []
                #import ipdb; ipdb.set_trace()
                for tile in tile_list:
                    #modscag_url_str = 'https://%s:%s@snow-data.jpl.nasa.gov/modscag-historic/%Y/%j/MOD09GA.A%Y%j.%s.005.2015006001833.snow_fraction.tif' % (uname, pw, tile)
                    #modscag_url_str = 'https://snow-data.jpl.nasa.gov/modscag-historic/%%Y/%%j/MOD09GA.A%%Y%%j.%s.005.*.snow_fraction.tif' % tile
                    modscag_url_fn = None
                    modscag_url_fn = parsed_html.findAll(text=re.compile('%s.*snow_fraction.tif' % tile))
                    if modscag_url_fn:
                        modscag_url_fn = modscag_url_fn[0]
                        modscag_url = os.path.join(modscag_url_base, modscag_url_fn)
                        print modscag_url
                        modscag_fn = os.path.join(outdir, os.path.split(modscag_url_fn)[-1])
                        if not os.path.exists(modscag_fn):
                            getfile2(modscag_url, auth=auth, outdir=outdir)
                        modscag_fn_list.append(modscag_fn)
                #Mosaic tiles
                cmd = ['gdalbuildvrt', '-vrtnodata', '255', out_vrt_fn]
                cmd.extend(modscag_fn_list)
                print cmd
                subprocess.call(cmd)
                #if not os.path.exists(out_tif_fn):
                    #Clip to extent of glaciers, EPSG:32611
                    #-256034 3923861 1600244 5726292
                    #cmd = ['gdalwarp', '-overwrite', '-r', 'cubic', '-t_srs', 'EPSG:32610', out_vrt_fn, out_tif_fn]
                    #print cmd
                    #subprocess.call(cmd)
                #out_vrt_fn_list.append(out_tif_fn)
                out_vrt_fn_list.append(out_vrt_fn)
    #s = malib.DEMStack(out_vrt_fn_list)
    #return s
    #modscag_ds = modscag_proc(out_vrt_fn_list)
    #return modscag_ds 
    return out_vrt_fn_list

def modscag_proc(fn_list, extent=None, t_srs=None):
    ds_list = warplib.memwarp_multi_fn(fn_list, extent=extent, t_srs=t_srs, r='cubicspline')
    dtype = np.uint8
    ma_stack = np.ma.array([np.ma.masked_greater(malib.ds_getma(ds), 100) for ds in np.array(ds_list)], dtype=dtype)
    stack_count = np.ma.masked_equal(ma_stack.count(axis=0), 0).astype(dtype)
    stack_count.set_fill_value(0)
    stack_min = ma_stack.min(axis=0).astype(dtype)
    stack_min.set_fill_value(0)
    stack_max = ma_stack.max(axis=0).astype(dtype)
    stack_max.set_fill_value(0)
    stack_med = np.ma.median(ma_stack, axis=0).astype(dtype)
    stack_med.set_fill_value(0)
    stack_fn = os.path.splitext(fn_list[0])[0] + '_' + os.path.splitext(os.path.split(fn_list[-1])[1])[0] + '_stack_%i' % len(fn_list) 
    out_fn = stack_fn + '_count.tif'
    malib.writeGTiff(stack_count, out_fn, ds_list[0])
    out_fn = stack_fn + '_max.tif'
    malib.writeGTiff(stack_max, out_fn, ds_list[0])
    out_fn = stack_fn + '_min.tif'
    malib.writeGTiff(stack_min, out_fn, ds_list[0])
    out_fn = stack_fn + '_med.tif'
    malib.writeGTiff(stack_med, out_fn, ds_list[0])
    ds = gdal.Open(out_fn)
    return ds

snodas_outdir = os.path.join(topdir, 'snodas')
if not os.path.exists(snodas_outdir):
    os.makedirs(snodas_outdir)

modscag_outdir = os.path.join(topdir, 'modscag')
if not os.path.exists(modscag_outdir):
    os.makedirs(modscag_outdir)

mask_fn=os.path.join(nlcd_dir,'rockmask_conus.tif')

if not os.path.exists(mask_fn):
    rockmask(nlcd_fn, mask_fn, glac_shp_fn)

mask_ds = gdal.Open(mask_fn)
dem_fn = sys.argv[1]
print dem_fn

#Extract DEM timestamp
dem_dt = timelib.fn_getdatetime(dem_fn)
dem_dir = os.path.split(os.path.split(dem_fn)[0])[0]
dem_ds = gdal.Open(dem_fn)

#Use top of atmosphere scaled reflectance values (0-1)
toa_fn = glob.glob(os.path.join(dem_dir, '*toa.tif'))
#import ipdb; ipdb.set_trace()
if not toa_fn:
    from subprocess import call
    cmd = ['toa.sh', dem_dir]
    call(cmd)
    toa_fn = glob.glob(os.path.join(dem_dir, '*toa.tif'))
toa_fn = toa_fn[0]

toa_ds = gdal.Open(toa_fn)

#Get SNODAS products for DEM timestamp
snowdas_ds = get_snodas(dem_dt, snodas_outdir)

#MODSCAG products
#t_srs = osr.SpatialReference()
#t_srs.ImportFromEPSG(32610)
modscag_fn_list = get_modscag(dem_dt, modscag_outdir)
modscag_ds = modscag_proc(modscag_fn_list, extent=dem_ds, t_srs=dem_ds)

#snowcover_fn = os.path.splitext(snowcover_stack.stack_fn)[0]+'_min.tif'
#snowcover_ds = gdal.Open(snowcover_fn)

#Warp all masks to DEM extent/res
#Note: use cubicspline here to avoid artifacts with negative values
#ds_list = warplib.memwarp_multi([dem_ds, mask_ds, snowcover_ds, toa_ds], res=dem_ds, extent=dem_ds, t_srs=dem_ds, r='cubicspline')
ds_list = warplib.memwarp_multi([dem_ds, mask_ds, snowdas_ds, modscag_ds, toa_ds], res=dem_ds, extent=dem_ds, t_srs=dem_ds, r='cubicspline')

#Load into memory 
dem = malib.ds_getma(ds_list[0])

#Rockmask
rockmask = ds_list[1].GetRasterBand(1).ReadAsArray()
#rockmask is 1 for valid rock, 0 for everything else (ndv)
rockmask = rockmask.astype(bool)

#SNODAS depth filter
snodas_max_depth = 0.2
#SNODAS snow depth values are mm, convert to meters
snodas_depth = malib.ds_getma(ds_list[2])/1000.
snodas_mask = np.ma.masked_greater(snodas_depth, snodas_max_depth)
#This should be 1 for valid surfaces with no snow, 0 for snowcovered surfaces
snodas_mask = ~(np.ma.getmaskarray(snodas_mask))

#MODSCAG
modscag_thresh = 50
modscag_perc = malib.ds_getma(ds_list[3])
modscag_mask = (modscag_perc.filled(0) >= modscag_thresh) 
#This should be 1 for valid surfaces with no snow, 0 for snowcovered surfaces
modscag_mask = ~(modscag_mask)

#Reflectance filter
toa_thresh = 0.4
toa = malib.ds_getma(ds_list[4])
toa_mask = np.ma.masked_greater(toa, toa_thresh)
#This should be 1 for valid surfaces, 0 for snowcovered surfaces
toa_mask = ~(np.ma.getmaskarray(toa_mask))

#Filter based on expected snowline
#Simplest approach is just an altitude cutoff
#max_elev = 1500 
#newdem = np.ma.masked_greater(newdem, max_elev)

#Identify snow-free rock surfaces
newmask = rockmask
#newmask = np.logical_and(snowcover_mask, newmask)
newmask = np.logical_and(toa_mask, newmask)
newmask = np.logical_and(modscag_mask, newmask)
#Now invert to use to create final masked array
newmask = ~newmask

#Check that we have enough pixels, good distribution

#Apply mask to original DEM - use these surfaces for co-registration
newdem = np.ma.array(dem, mask=newmask)

#Write out
out_fn = os.path.splitext(dem_fn)[0]+'_snodas_depth.tif'
print "Writing out %s" % out_fn
malib.writeGTiff(snodas_depth, out_fn, src_ds=dem_ds)
out_fn = os.path.splitext(dem_fn)[0]+'_modscag_perc.tif'
print "Writing out %s" % out_fn
malib.writeGTiff(modscag_perc, out_fn, src_ds=dem_ds)
out_fn = os.path.splitext(dem_fn)[0]+'_rockmask.tif'
print "Writing out %s" % out_fn
malib.writeGTiff(rockmask, out_fn, src_ds=dem_ds)
out_fn = os.path.splitext(dem_fn)[0]+'_snowdas_mask.tif'
print "Writing out %s" % out_fn
malib.writeGTiff(snodas_mask, out_fn, src_ds=dem_ds)
out_fn = os.path.splitext(dem_fn)[0]+'_modscag_mask.tif'
print "Writing out %s" % out_fn
malib.writeGTiff(modscag_mask, out_fn, src_ds=dem_ds)
out_fn = os.path.splitext(dem_fn)[0]+'_toamask.tif'
print "Writing out %s" % out_fn
malib.writeGTiff(toa_mask, out_fn, src_ds=dem_ds)
out_fn = os.path.splitext(dem_fn)[0]+'_ref.tif'
print "Writing out %s" % out_fn
malib.writeGTiff(newdem, out_fn, src_ds=dem_ds)

print 
#pc_align_wrapper.sh 
