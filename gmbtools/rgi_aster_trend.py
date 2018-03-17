#! /usr/bin/env python

"""
Script to create DEMStack objects from ASTER DEMs for RGI polygons

To do:
Should run in parallel with multiprocessing, similar to mb_parallel
Implement RANSAC linear regression
Implement smoothness constraint to penalize large spatial gradients in trend

"""

import os
import sys
from osgeo import ogr
from pygeotools.lib import geolib, warplib, malib

topdir='/nobackup/deshean/'
asterdir = os.path.join(topdir, 'hma/aster/align')
stackdir = os.path.join(asterdir, 'stack_n5_1825days')
if not os.path.exists(stackdir):
    os.makedirs(stackdir)

def rgi_name(feat, glacname_fieldname='Name', glacnum_fieldname='RGIId'):
    #glacname = feat.GetField(feat.GetFieldIndex(glacname_fieldname))
    glacname = feat.GetField(glacname_fieldname)
    if glacname is None:
        glacname = ""
    else:
        #RGI has some nonstandard characters
        glacname = glacname.decode('unicode_escape').encode('ascii','ignore')
        glacname = glacname.replace(" ", "")
        glacname = glacname.replace("_", "")
        glacname = glacname.replace("/", "")

    #glacnum = feat.GetField(feat.GetFieldIndex(glacnum_fieldname))
    glacnum = feat.GetField(glacnum_fieldname)
    fn = feat.GetDefnRef().GetName()
    #RGIId (String) = RGI50-01.00004
    #glacnum = '%0.5f' % float(glacnum.split('-')[-1])
    glacnum = float(glacnum.split('-')[-1])

    if glacname:
        feat_fn = "%08.5f_%s" % (glacnum, glacname)
    else:
        feat_fn = "%08.5f" % glacnum
    return feat_fn

glac_shp_fn = os.path.join(topdir,'data/rgi60/regions/rgi60_merge_HMA_aea.shp')
glac_shp_ds = ogr.Open(glac_shp_fn, 0)
glac_shp_lyr = glac_shp_ds.GetLayer()
#This should be contained in features
glac_shp_srs = glac_shp_lyr.GetSpatialRef()
feat_count = glac_shp_lyr.GetFeatureCount()
print("Input glacier polygon count: %i" % feat_count)

#Area filter
min_glac_area = 10 
glac_shp_lyr.SetAttributeFilter("Area > %s" % min_glac_area)
feat_count = glac_shp_lyr.GetFeatureCount()
print("Min. Area filter glacier polygon count: %i" % feat_count)
glac_shp_lyr.ResetReading()

#Generate ASTER index
#gdaltindex -t_srs EPSG:4326 aster_align_index.shp *align/*align.tif
#ogr2ogr -t_srs '+proj=aea +lat_1=25 +lat_2=47 +lat_0=36 +lon_0=85 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ' aster_align_index_aea.shp aster_align_index.shp

aster_index_fn = os.path.join(asterdir, 'aster_align_index_aea.shp')
aster_index_ds = ogr.Open(aster_index_fn, 0)
aster_index_lyr = aster_index_ds.GetLayer()
#This should be contained in features
aster_index_srs = aster_index_lyr.GetSpatialRef()
feat_count = aster_index_lyr.GetFeatureCount()
print("Input ASTER count: %i" % feat_count)

min_aster_count = 5
min_dt_ptp = 1825

for n, feat in enumerate(glac_shp_lyr):
    #glac_geom_orig = geolib.geom_dup(feat.GetGeometryRef())
    feat_fn = rgi_name(feat)
    print(n, feat_fn)
    glac_geom = geolib.geom_dup(feat.GetGeometryRef())
    glac_geom_extent = geolib.geom_extent(glac_geom)

    #Spatial filter
    aster_index_lyr.SetSpatialFilter(glac_geom)
    aster_count = aster_index_lyr.GetFeatureCount()
    print("ASTER count after spatial filter: %i" % aster_count)
    if aster_count > min_aster_count:
        fn_list = []
        for aster_feat in aster_index_lyr:
            #Only 1 field from gdaltindex, 'location'
            fn = os.path.join(asterdir, aster_feat.GetField(0))
            fn_list.append(fn)
        stack = malib.DEMStack(fn_list, outdir=os.path.join(stackdir, feat_fn), \
                res='max', extent=glac_geom_extent, srs=aster_index_srs, mask_geom=glac_geom, n_thresh=min_aster_count, min_dt_ptp=min_dt_ptp)
        #if stack.ma_stack is not None:
            #sys.exit()
            #glac_geom_mask = geolib.geom2mask(glac_geom, stack.get_ds())
            #ds_list = warplib.memwarp_multi_fn(fn_list, res='max', extent=glac_geom_extent)
    aster_index_lyr.ResetReading()

#Generate plots
#hs.sh */*mean.tif
#for i in */*trend.tif; do  imviewer.py -clim -10 10 -cmap RdYlBu -label 'Trend (m/yr)' -of png -overlay $(echo $i | sed 's/_trend/_mean_hs_az315/') $i -scalebar -outsize 8 8 -dpi 100; done
#for i in */*count.tif; do imviewer.py -clim 0 20 -cmap inferno -label 'Count' -of png -overlay $(echo $i | sed 's/_count/_mean_hs_az315/') $i -scalebar -outsize 8 8 -dpi 100; done
#for i in */*[0-9]_std.tif; do imviewer.py -clim 0 30 -label 'Std (m)' -of png -overlay $(echo $i | sed 's/_std/_mean_hs_az315/') $i -scalebar -outsize 8 8 -dpi 100; done
