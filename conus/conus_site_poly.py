#! /usr/bin/env python

#David Shean
#10/4/16
#dshean@gmail.com

#Identify DEMs for each site polygon
#Generate DEM stack

import sys
import os
import numpy as np

import ogr
import osr
import gdal

from lib import malib
from lib import geolib
from lib import timelib

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
outdir = 'sites_poly'
if not os.path.exists(outdir):
    os.makedirs(outdir)

#Output stack parameters
stack_res = 2 
#stack_res = 8 
buffer = 1000
min_area = 1000000
#UTM 11N
#stack_epsg = 32611
#UTM 10N
stack_epsg = 32610
#This is HMA UTM
#stack_epsg = 32644

dst_srs = osr.SpatialReference()
#dst_srs = geolib.wgs_srs
dst_srs.ImportFromEPSG(stack_epsg)

#site_name_field = 1 
#For 24k glacier polygons
site_name_field = 4 
dem_name_field = 0 

for n,site_feat in enumerate(site_shp_lyr):
    site_name = site_feat.GetFieldAsString(site_name_field) 
    site_name = site_name.replace(" ", "_").lower()
    #Could output json here
    stackdir=os.path.join(outdir,site_name)
    if not os.path.exists(stackdir):
        os.makedirs(stackdir)
    out_fn = os.path.join(stackdir, 'site_%s_%i.csv' % (site_name, stack_epsg))
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
    site_extent = site_geom.GetEnvelope()
    stack_extent = [site_extent[0], site_extent[2], site_extent[1], site_extent[3]]
    
    print
    print site_name
    print site_extent
    print np.array(geolib.geom_wh(site_geom))/1000.
    print site_geom.Area()/1E6
    print

    #Get output filename from feature, make sure datestr is included
    print>>f, site_name
    #print>>f, dst_srs.ExportAsProj4()
    print>>f, site_extent 
    
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
        if not igeom.IsEmpty():
            if igeom.Area() > min_area:
                dem_fn_list.append(dem_name)
            #Write out interseciton geom to new file
        #dem_feat.Destroy()

    #Hack to use different DEM resolutions
    #dem_fn_list = ['%s.tif' % i for i in dem_fn_list]
    #dem_fn_list = [i.replace('DEM_32m_trans_32611_100m','DEM_8m_trans.tif') for i in dem_fn_list]
    dem_fn_list = [i.replace('DEM_32m_trans_32611_100m','DEM_%im_trans.tif' % stack_res) for i in dem_fn_list]
    print>>f, dem_fn_list
    print "Generating stack"
    s = malib.DEMStack(fn_list=dem_fn_list, outdir=stackdir, res=stack_res, extent=stack_extent, srs=dst_srs, trend=True)

    #temp_lyr = temp_ds.CreateLayer('out', None, ogr.wkbPolygon)
    #temp_lyr.CreateFeature(feat)

    #Mask stack with glacier poly

    dem_dt_list = [timelib.fn_getdatetime(i) for i in dem_fn_list]

    #Isolate summer DEMs
    #dem_fn_list_summer = 
    #Compute seasonal balance each year, composites each spring and each summer, difference
    #See stack_fltr
    #summer = timelib.get_dt_bounds(dem_dt_list, min_rel_dt=(8,15), max_rel_dt=(10,15))
    #spring = timelib.get_dt_bounds(dem_dt_list, min_rel_dt=(4,15), max_rel_dt=(6,15))
    #Alternatively, pick an ideal date, then find everything around
    #summer = timelib.get_closest_dt_padded_idx(datetime(

    #Long-term elevation change NED
    #Clip to glacier extent
    #Stats

    #site_feat.Destroy()
    f = None
