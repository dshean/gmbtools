#! /bin/bash

#Many of the CONUS lidar datasets have horizontal and vertical units of US feet
#See doc/lidar_notes_20161024.txt for specific examples

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


