#! /bin/bash

cd /Volumes/d/lidar
list=final_products.txt
#fn_list=$(ls -L -Sr $(cat $list | grep -v ^#))
#fn_list=$(cat $list | grep -v -e '^$' | grep -v ^# | tail -n 3)
#fn_list=$(ls -L $(cat $list | grep -v -e '^$' | grep -v ^#) | grep Martel)
fn_list=$(ls -L $(cat $list | grep -v -e '^$' | grep -v ^#))

echo $fn_list

#ned=/Volumes/d/ned/tiles_1arcsec/ned1_tiles_glac24k_115kmbuff.vrt
ned=/Volumes/d/ned/NED_2003_1arcsec/ned1_2003_adj.vrt

gdal_opt="-co COMPRESS=LZW -co TILED=YES -co BIGTIFF=IF_SAFER"

#gdal_translate -outsize 3.125% 3.125% $gdal_opt {.}_32610_1m.tif {.}_32610_32m.tif; \

parallel -j 7 "echo {}; if [ ! -e {.}_32610_1m.tif ] ; then gdalwarp -tr 1 1 $gdal_opt -t_srs EPSG:32610 -r cubic -dstnodata -9999 {} {.}_32610_1m.tif; fi; \
if [ ! -e {.}_32610_1m.tif.ovr ] ; then ~/src/demtools/gdaladdo_ro.sh {.}_32610_1m.tif; fi; \
if [ ! -e {.}_32610_32m.tif ] ; then gdalwarp -overwrite -tr 32 32 -r cubic $gdal_opt {.}_32610_1m.tif {.}_32610_32m.tif; fi; \
if [ ! -e {.}_32610_32m_hs_az315.tif ] ; then ~/src/demtools/hs.sh {.}_32610_32m.tif; fi; \
if [ ! -e ned_diff/$(basename ${ned%.*})_$(basename {.}_32610_32m_dz_eul.tif) ] ; then compute_dz.py -outdir ned_diff $ned {.}_32610_32m.tif; fi;" ::: $fn_list

exit

parallel 'imviewer.py -clim -20 20 -cmap RdYlBu {} -of png' ::: ned_diff/ned*eul.tif

#list=$(find . -name '*32610_1m.tif')
#gdalbuildvrt conus_lidar_1m.vrt $list
list=$(find . -name '*32610_32m.tif')
gdalbuildvrt conus_lidar_32m.vrt $list
~/src/demtools/hs.sh conus_lidar_32m.vrt
~/src/demtools/gdaladdo_ro.sh conus_lidar_32m.vrt
~/src/demtools/gdaladdo_ro.sh conus_lidar_32m_hs_az315.tif

#didoscpput $list /nobackupp8/deshean/rpcdem/conus_lidar/
