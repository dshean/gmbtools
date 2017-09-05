#! /bin/bash

#This script can be used to prepare a new basemap for a second iteration of high-quality DEM generation

#dem=rainier_stack_all-tile-0.tif
dem=$1
ref_src=/nobackup/deshean/rpcdem/ned13/ned13_tiles_glac24k_115kmbuff.vrt 

#Prepare reference DEM
warptool.py -outdir . -te $dem -tr $dem -t_srs $dem $ref_src
ref=$(basename $ref_src)
ref=${ref%.*}_warp.tif

#Remove outliers based elev diff
#Note, NED is a DSM, so will be ~20-50 m offsets for forested areas
filter.py $dem -filt dz -param $ref 0 100 
dem=${dem%.*}_dzfilt*.tif
#filter.py $dem -filt med -param 5 
#dem=${dem%.*}_medfilt_5px.tif

#Should check maximum dimensions of remaining holes, determine size of filter necessary
gauss_fn=''
for s in 5 11 21
do
    filter.py $dem -filt gauss -param $s 
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
