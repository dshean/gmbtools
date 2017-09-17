#! /bin/bash

#This script can be used to prepare a new basemap for a second iteration of high-quality DEM generation

#This assumes output from conus_site_poly.py
#dem=rainier_stack_all-tile-0.tif

dem=$1
#ref_src=/nobackup/deshean/rpcdem/ned13/ned13_tiles_glac24k_115kmbuff.vrt 
ref_src=/nobackup/deshean/hma/mos/latest/mos_8m/*mos_8m.vrt
#ref_src=/nobackup/deshean/hma/mos/latest/mos_32m/*mos_32m.vrt

#Prepare reference DEM
ref_bn=$(basename $ref_src)
if [ ! -e ${ref_bn%.*}_warp.tif ] ; then 
    warptool.py -outdir . -te $dem -tr $dem -t_srs $dem $ref_src
fi
ref=${ref_bn%.*}_warp.tif

#Remove outliers based elev diff
#Note, NED is a DSM, so will be ~20-50 m offsets for forested areas
if [ ! -e ${dem%.*}_dzfilt_0.00-100.00.tif ] ; then 
    filter.py $dem -filt dz -param $ref 0 100 
fi
dem=${dem%.*}_dzfilt_0.00-100.00.tif
#filter.py $dem -filt med -param 5 
#dem=${dem%.*}_medfilt_5px.tif

#Should check maximum dimensions of remaining holes, determine size of filter necessary
gauss_fn=''
for s in 5 11 21
do
    if [ ! -e ${dem%.*}_gaussfilt_${s}px.tif ] ; then 
        filter.py $dem -filt gauss -param $s 
    fi
    gauss_fn+=" ${dem%.*}_gaussfilt_${s}px.tif"
done

#Merge with NED
#dem_mosaic --priority-blending-length 5 -o ${dem%.*} $dem $ref

#Try to fill holes with dem_mosaic
#Fails for large values, won't fill near edges
#dem_mosaic --hole-fill-length 9999 -o ${dem%.*} $dem

#gdal_fillnodata.py
#dem_downsample_fill.py
#inpaint_dem.py

#Merge gauss pyramid
dem_mosaic --priority-blending-length 5 -o ${dem%.*}_gaussfill $dem $gauss_fn $ref

echo ${dem%.*}_gaussfill-tile-0.tif
