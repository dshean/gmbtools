#! /bin/bash

#Many of the CONUS lidar datasets have horizontal and vertical units of US feet
#See doc/lidar_notes_20161024.txt for specific examples

gdal_opt="-co COMPRESS=LZW -co TILED=YES -co BIGTIFF=IF_SAFER"

in=$1
#Convert vertical units to meters
image_calc -o ${in%.*}_unitmeters.tif -c "var_0*0.3048" -d float32 --threads 7 $in 

#Most WA lidar is ESRI:102749 
#http://spatialreference.org/ref/esri/nad-1983-stateplane-washington-south-fips-4602-feet/
#dem_geoid has trouble with NAD83(HARN)
#Use EPSG:2285 for WA state plane north
#Use EPSG:2286 for WA state plane south
gdal_edit.py -a_srs EPSG:2285 ${in%.*}_unitmeters.tif 

#Now remove geoid offset
dem_geoid --reverse-adjustment ${in%.*}_unitmeters.tif

#Output values are meters relative to WGS84 ellipsoid
#This tif should be fine for comparison with other datasets, could also reproject to something more desirable
echo ${in%.*}_unitmeters-adj.tif

i=${in%.*}_unitmeters-adj

if [ ! -e ${i}_32610_1m.tif ] ; then
    gdalwarp -tr 1 1 $gdal_opt -t_srs EPSG:32610 -r cubic -dstnodata -9999 $i ${i}_32610_1m.tif
fi

if [ ! -e ${i}_32610_1m.tif.ovr ] ; then 
    ~/src/demtools/gdaladdo_ro.sh ${i}_32610_1m.tif
fi

if [ ! -e ${i}_32610_32m.tif ] ; then 
    gdalwarp -overwrite -tr 32 32 -r cubic $gdal_opt ${i}_32610_1m.tif ${i}_32610_32m.tif
fi

if [ ! -e ${i}_32610_32m_hs_az315.tif ] ; then 
    hs.sh ${i}_32610_32m.tif
fi

#ned=/Volumes/d/ned/NED_2003_1arcsec/ned1_2003_adj.vrt
#if [ ! -e ned_diff/$(basename ${ned%.*})_$(basename ${i}_32610_32m_dz_eul.tif) ] ; then 
#    compute_dz.py -outdir ned_diff $ned ${i}_32610_32m.tif
#fi

#Now build supermosaic
#list=$(find . -name '*32610_32m.tif')
#gdalbuildvrt conus_lidar_32m.vrt $list
#hs.sh conus_lidar_32m.vrt
#gdaladdo_ro.sh conus_lidar_32m.vrt
#gdaladdo_ro.sh conus_lidar_32m_hs_az315.tif
