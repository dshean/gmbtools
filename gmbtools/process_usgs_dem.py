#! /usr/bin/env python

import sys, os
from datetime import datetime
import re
import subprocess

#Replace spaces with underscore
#for i in *\ *; do mv "$i" $(echo $i | sed 's/ //g'); done

#Replace ( in filesnames
#mv wolverine_17sep08\(2\).dem wolverine_17sep08_2.dem

#Replace _fall_99, assume oct 1
#mv blue_fall_99.dem blue_fall_01oct99.dem
#mv north_klawatti_fall_99.dem north_klawatti_fall_01oct99.dem
#mv south_cascade_fall_99.dem south_cascade_fall_01oct99.dem

#Replace year only
#mv south_cascade_2008_NTM.dem south_cascade_NTM_guess_01oct08.dem
#mv south_cascade_2008_WV.dem south_cascade_WV_guess_01oct08.dem

#Remove trailing _
#for i in *_.dem; do mv "$i" $(echo $i | sed 's/_.dem/.dem/'); done
#mv Columbia_29mar12_5m_DTM_.tif Columbia_29mar12_5m_DTM.tif

#Random
#mv gukana_02sep10.dem gulkana_02sep10.dem
#mv littlejarvis_24jul07.dem little_jarvis_24jul07.dem
#mv south_cascade_20010907.dem south_cascade_07sep01.dem
#mv weasel_collar24sep05.dem weasel_collar_24sep05.dem
#mv sandalee12may11.dem sandalee_12may11.dem

#Fix spelling of monthname
#ls *l[0-9]* | grep -v jul
#mv bear_lake_24augl13.dem bear_lake_24aug13.dem
#mv bear_lake_30augl14.dem bear_lake_30aug14.dem
#mv bear_lake_30augl14.dth bear_lake_30aug14.dth
#mv bear_lake_30augl14.dtp bear_lake_30aug14.dtp
#mv polychrome_24sepl13.dem polychrome_24sep13.dem
#mv south_cascade_17augl12.dem south_cascade_17aug12.dem

#Find redo or _2

in_fn = sys.argv[1]
#fn = fn.replace('_comm','')
fn = in_fn.replace('__','_')

d = re.findall(r'[0-9][0-9][a-z][a-z][a-z][0-9][0-9]', fn)
if not d:
    print("No datestring: %s" % fn)
else:
    d = d[0]
    dt = datetime.strptime(d,'%d%b%y')

d_list = os.path.splitext(fn)[0].split('_')
d_list = [x for x in d_list if x != d]

fn_prefix = '_'.join(map(str, d_list))
#fn_prefix = fn.split(d)[0].rstrip('_')

outdir = 'clean'
if not os.path.exists(outdir):
    os.makedirs(outdir)
out_fn = '%s_%s.tif' % (dt.strftime('%Y%m%d'), fn_prefix)
out_fn = out_fn.lower()
out_fn = os.path.join(outdir, out_fn)

cmd = ['gdal_translate', '-co', 'TILED=YES', '-co', 'COMPRESS=LZW', '-co', 'BIGTIFF=IF_SAFER', in_fn, out_fn]
print(cmd)
subprocess.call(cmd)

#Create new tif with clean fn
#parallel '~/src/conus/conus/process_usgs_dem.py {}' ::: *dem *tif
#Now separate alaska from conus
#mkdir alaska; mkdir conus
#for i in *tif; do lat=$(gdalinfo $i | grep 'Lower Right' | awk -F',' '{print $NF}' | awk -F'd' '{print $1}'); if [ "$lat" -gt "50" ]; then mv $i alaska/ ; else mv $i conus/ ; fi ; done
#cd conus

export PROJ_DEBUG=4

#Replace nodata values
cd conus
parallel 'replace_ndv.py {} -9999' ::: *tif
parallel "gdalwarp -r cubic -t_srs '+proj=aea +lat_1=36 +lat_2=49 +lat_0=43 +lon_0=-115 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84+units=m +no_defs ' {} {.}_aea.tif" ::: *_ndv.tif
raster2shp.py -merge_fn usgs_conus_dem_aea.shp *aea.tif

cd ../alaska
parallel 'replace_ndv.py {} -9999' ::: *tif
parallel 'gdalwarp -r cubic -t_srs EPSG:3338 {} {.}_aa.tif' ::: *ndv.tif
raster2shp.py -merge_fn usgs_alaska_dem_aa.shp *aa.tif

#Reproject and convert to WGS84
#Note datum_convert does not reproject, it just does horiz/vert offset
#conus_proj='+proj=aea +lat_1=36 +lat_2=49 +lat_0=43 +lon_0=-115 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs '
#parallel 'datum_convert --t_srs "$conus_proj" {} {.}_aea.tif' ::: *_ndv.tif
#parallel "gdalwarp -r cubic -t_srs '+proj=aea +lat_1=36 +lat_2=49 +lat_0=43 +lon_0=-115 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ' {} {.}_aea.tif; replace_ndv.py {.}_aea.tif 0; mv {.}_aea_ndv.tif {.}_aea.tif" ::: *_ndv.tif

