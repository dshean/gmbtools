#! /usr/bin/env python
"""
Identify CONUS sites, isolate DEMs, and make stack for each site
"""

#mkdir site_poly_highcount_rect3
#lfs setstripe --count 32 site_poly_highcount_rect3

#Requires predefined polygons for sites with name field
#Requires shp of DEM footprints from raster2shp.py
#cd /nobackup/deshean/conus/dem2
#raster2shp.py -merge_fn shp/conus_32m_trans_n284_20170202.shp */*/*32m_trans.tif
#~/src/conus/conus/conus_site_poly.py ../shp/conus_site_poly_highcount_rect_32611.shp  shp/conus_32m_n366_20160928_1436.shp

import sys
import os
from datetime import datetime, timedelta
import subprocess
import time

import numpy as np
from osgeo import gdal, ogr, osr

from pygeotools.lib import malib
from pygeotools.lib import geolib
from pygeotools.lib import timelib
from pygeotools.lib import warplib

if len(sys.argv) != 3:
    sys.exit('Usage is %s sitepoly.shp dempoly.shp' % sys.argv[0])

#This contains polygons defining study areas
#Could also sort glacier polygons by area, then just use largest directly
site_shp_fn = sys.argv[1]
if not os.path.exists(site_shp_fn):
    sys.exit('Unable to find input shp: %s' % site_shp_fn)

site_shp_ds = ogr.Open(site_shp_fn)
site_shp_lyr = site_shp_ds.GetLayer()
site_shp_srs = site_shp_lyr.GetSpatialRef()

#This contains merged shp of outlines for available DEMs
dem_shp_fn = sys.argv[2]
if not os.path.exists(dem_shp_fn):
    sys.exit('Unable to find input shp: %s' % dem_shp_fn)

dem_shp_ds = ogr.Open(dem_shp_fn)
dem_shp_lyr = dem_shp_ds.GetLayer()
dem_shp_srs = dem_shp_lyr.GetSpatialRef()

#outdir = os.path.splitext(shp_fn)[0]+'_stack'
#outdir = 'range_poly'
outdir = 'site_poly_highcount_rect3_rerun'
if not os.path.exists(outdir):
    os.makedirs(outdir)

#Output stack parameters
res = 8 
#buffer = 1000
buffer = None
min_area = 1000000

#Use Albers Equal Area
dst_srs = geolib.conus_aea_srs

#For rectangles
site_name_field = 1 
#For 24k glacier polygons
#site_name_field = 4 
dem_name_field = 0 

mos_cmd_list = []
std_cmd_list = []
dz_cmd_list = []
for n,site_feat in enumerate(site_shp_lyr):
    site_name = site_feat.GetFieldAsString(site_name_field) 
    site_name = site_name.replace(" ", "_").lower()
    #Could output json here
    stackdir=os.path.join(outdir,site_name)
    if not os.path.exists(stackdir):
        os.makedirs(stackdir)
    out_fn = os.path.join(stackdir, 'site_%s.csv' % site_name)
    f = open(out_fn, 'w')

    site_geom_orig = site_feat.GetGeometryRef()
    site_geom_orig.AssignSpatialReference(site_shp_srs)
    geolib.geom_transform(site_geom_orig, dst_srs)
    
    if buffer is not None:
        site_geom_buff = site_geom_orig.Buffer(buffer)
        site_geom = site_geom_buff
    else:
        site_geom = site_geom_orig

    #This is xmin, xmax, ymin, ymax
    site_extent = geolib.geom_extent(site_geom)
    
    print("========================\n")
    print("Site name: %s" % site_name)
    print("Extent:", site_extent)
    print("Width/Height (km)", np.array(geolib.geom_wh(site_geom))/1000.)
    print("Area (km2)", site_geom.Area()/1E6)

    #Get output filename from feature, make sure datestr is included
    f.write(site_name+'\n')
    f.write(str(site_extent)+'\n')
    
    #Can also set spatial filter for each site
    #https://pcjericks.github.io/py-gdalogr-cookbook/vector_layers.html#spatial-filter

    dem_fn_list = []
    dem_shp_lyr.ResetReading()
    for m,dem_feat in enumerate(dem_shp_lyr):
        dem_name = dem_feat.GetFieldAsString(dem_name_field)
        #Date
        dem_geom = dem_feat.GetGeometryRef()
        dem_geom.AssignSpatialReference(dem_shp_srs)
        geolib.geom_transform(dem_geom, dst_srs)

        #Add interseciton area filter
        igeom = site_geom.Intersection(dem_geom)
        if igeom:
            if not igeom.IsEmpty():
                if igeom.Area() > min_area:
                    dem_fn_list.append(dem_name)
                #Write out interseciton geom to new file
        #dem_feat.Destroy()

    if dem_fn_list:
        #Sort, assuming date is first key in file name
        dem_fn_list = sorted(dem_fn_list, key=lambda x: os.path.split(x)[-1])
        print(dem_fn_list)

        #Hack to use different DEM resolutions
        dem_ext = dem_fn_list[0].split('-')[-1]
        dem_fn_list = [i.replace(dem_ext,'DEM_%im_trans.tif' % res) for i in dem_fn_list]
        for fn in dem_fn_list:
            f.write(fn+'\n')
        f = None

        #print "Generating stack"
        #s = malib.DEMStack(fn_list=dem_fn_list, outdir=stackdir, res=res, extent=site_extent, srs=dst_srs, trend=True)

        dem_fn_list = np.array(dem_fn_list)
        dem_dt_list = np.array([timelib.fn_getdatetime(i) for i in dem_fn_list])
    
        #These are OrderedDict
        summer_dict = timelib.dt_filter_rel_annual_idx(dem_dt_list, min_rel_dt=(8,1), max_rel_dt=(10,31))
        spring_dict = timelib.dt_filter_rel_annual_idx(dem_dt_list, min_rel_dt=(4,1), max_rel_dt=(6,15))

        mos_fn_dict = {}
        for mos_prefix, dt_dict in (('spring', spring_dict), ('summer', summer_dict)):
            for year, dt_idx in dt_dict.iteritems():
                #print(year, dt_idx)
                mos_fn_list = dem_fn_list[dt_idx]
                mos_dt_list = dem_dt_list[dt_idx]
                mos_fn = '%s_%i_%s_n%i_%s-%s' % (site_name, year, mos_prefix, mos_fn_list.size, \
                        min(mos_dt_list).strftime('%Y%m%d'), max(mos_dt_list).strftime('%Y%m%d'))
                mos_out_fn = os.path.join(stackdir, mos_fn) 
                print(mos_out_fn)
                mos_fn_dict[(year, mos_prefix)] = mos_out_fn+'-tile-0.tif'
                cmd = geolib.get_dem_mosaic_cmd(mos_fn_list, mos_out_fn, tr=res, t_srs=dst_srs, t_projwin=site_extent, threads=1)
                if not os.path.exists(mos_out_fn+'-tile-0.tif'):
                    mos_cmd_list.append(cmd)
                if len(mos_fn_list) > 1 and not os.path.exists(mos_out_fn+'-tile-0-stddev.tif'):
                    cmd2 = list(cmd)
                    cmd2.append('--stddev')
                    std_cmd_list.append(cmd2)

        #Generate difference map commands
        if len(mos_fn_dict) > 1:
            validyears = np.unique([i[0] for i in mos_fn_dict.keys()])
            for y in validyears:
                if (y, 'spring') in mos_fn_dict:
                    if (y, 'summer') in mos_fn_dict:
                        dz_cmd_list.append(['compute_dz.py', mos_fn_dict[(y, 'spring')], mos_fn_dict[(y, 'summer')]])
                if (y, 'summer') in mos_fn_dict:
                    if (y+1, 'spring') in mos_fn_dict:
                        dz_cmd_list.append(['compute_dz.py', mos_fn_dict[(y, 'summer')], mos_fn_dict[(y+1, 'spring')]])
                if (y, 'summer') in mos_fn_dict:
                    if (y+1, 'summer') in mos_fn_dict:
                        dz_cmd_list.append(['compute_dz.py', mos_fn_dict[(y, 'summer')], mos_fn_dict[(y+1, 'summer')]])
                if (y, 'spring') in mos_fn_dict:
                    if (y+1, 'spring') in mos_fn_dict:
                        dz_cmd_list.append(['compute_dz.py', mos_fn_dict[(y, 'spring')], mos_fn_dict[(y+1, 'spring')]])

            #Make stack of all year/season products 
            stack_cmd = ['make_stack.py', '--trend', '-outdir', os.path.join(stackdir, 'stack_all')]
            stack_cmd.extend(dem_fn_list)
            dz_cmd_list.append(stack_cmd)

            #Make stack of all year/season products 
            stack_cmd = ['make_stack.py', '--trend', '-outdir', os.path.join(stackdir, 'stack_seasonal_all')]
            stack_fn_list = mos_fn_dict.values()
            stack_fn_list.sort()
            stack_cmd.extend(stack_fn_list)
            dz_cmd_list.append(stack_cmd)

            #Make stack of all year/summer products
            stack_cmd = ['make_stack.py', '--trend', '-outdir', os.path.join(stackdir, 'stack_seasonal_summer')]
            stack_fn_list = [mos_fn_dict[k] for k in mos_fn_dict.keys() if 'summer' in k]
            if stack_fn_list:
                stack_fn_list.sort()
                stack_cmd.extend(stack_fn_list)
                dz_cmd_list.append(stack_cmd)

        """
                #Compute seasonal balance each year, composites each spring and each summer, difference

                neddir = '/nobackup/deshean/rpcdem/ned13'
                ned = os.path.join(neddir, 'ned13_tiles_glac24k_115kmbuff.vrt')

                warp_ds = warplib.memwarp_multi_fn([summer_fn, ned], res='first', extent='first', t_srs='first', r='cubic')
                warplib.writeout(warp_ds[1], os.path.join(stackdir, '%s_ned13_warp.tif' % site_name))

                #glac_shp = '/nobackup/deshean/rpcdem/nlcd/24k_selection_32610.shp' 
                #glac_shp = '/nobackup/deshean/rpcdem/nlcd/24k_selection_32610_dissolve.shp' 
                #glac_mask = geolib.shp2array(glac_shp, warp_ds[0])

                #Rasterize NED source date
                #Long-term elevation change NED
                #Clip to glacier extent
                #Stats
        """

    #site_feat.Destroy()

from concurrent.futures import ThreadPoolExecutor
threads = 8 
delay = 3.0
outf = open(os.devnull, 'w')

if mos_cmd_list:
    with ThreadPoolExecutor(max_workers=threads) as executor:
        for cmd in mos_cmd_list:
            #print(cmd)
            #executor.submit(subprocess.call, cmd, stdout=outf, stderr=subprocess.STDOUT)
            executor.submit(subprocess.call, cmd)
            time.sleep(delay)

if std_cmd_list:
    with ThreadPoolExecutor(max_workers=threads) as executor:
        for cmd in std_cmd_list:
            #print(cmd)
            #executor.submit(subprocess.call, cmd, stdout=outf, stderr=subprocess.STDOUT)
            executor.submit(subprocess.call, cmd)
            time.sleep(delay)

if dz_cmd_list:
    with ThreadPoolExecutor(max_workers=threads) as executor:
        for cmd in dz_cmd_list:
            #print(cmd)
            #executor.submit(subprocess.call, cmd, stdout=outf, stderr=subprocess.STDOUT)
            executor.submit(subprocess.call, cmd)
            time.sleep(delay)

outf = None

#Want to mask output to rock+ice
#parallel -j 32 'dem_mask.py --no_icemask {}' ::: */*tile-0-stddev.tif */*eul.tif */stack*/*trend.tif */stack*/*std.tif

#cd site_poly_highcount_rect/
#parallel 'compute_dh.py {}/{}_ned13_warp.tif {}/{}_summer-tile-0.tif' ::: *
#parallel 'clip_raster_by_shp.sh {} /nobackup/deshean/rpcdem/nlcd/24k_selection_32610_dissolve.shp' ::: */*eul.tif
#hs.sh */*tile-0.tif
#parallel 'imview.py {}/*clip.tif -overlay {}/*az315.tif -of png -cmap inferno_r -clim -50 0' ::: *
#parallel 'imview.py {}/*clip.tif -overlay {}/*az315.tif -of png -cmap RdYlBu -alpha 0.7 -clim -80 80' ::: *
