#! /bin/bash

#Create source date products for NED

#Downloaded all tiles
#Processed
cd /Volumes/d/ned/NED_2003_1arcsec

#Create mosaic in original WGS84 coord
fn_list=$(ls -d dem*/ | sed 's#/##g')
gdalbuildvrt ned1_2003_wgs84.vrt $fn_list
gdaladdo_ro.sh ned1_2003_wgs84.vrt

shp_fn=meta0306_PAL_24k_10kmbuffer_clean.shp
#Split using QGSI tool Vector -> Data Mangement Tools -> Split Vector Layer
#Split based on META0306 (seems to be tile ID)
d_list=$(ogrinfo -al $shp_fn | grep META0306 | awk '{print $NF}' | sed '1,2d' | sort -u)
#for d in $d_list ; do ogr2ogr -where "META0306=$d" meta0306_split/${d}_s_date.shp $shp_fn ; done
cd meta0306_split
shp_list=$(ls *.shp)
for i in $shp_list
do
    date=$(ogrinfo -al $i | grep S_DATE_CLN | awk '{print $NF}' | sed '1d')
    shp_rename.sh $i ${i%.*}_S_DATE_CLN_${date}.shp
done
shp_list=$(ls *.shp)
parallel 'clip_raster_by_shp.py -out_fn {.}.tif -extent shp ../ned1_2003_wgs84.vrt {}' ::: *shp  

exit

#Previous processing

epsg=26911
fn_list=$(ls -d dem*/)
parallel "gdalwarp -co COMPRESS=LZW -co TILED=YES -co BIGTIFF=IF_SAFER -overwrite -r cubic -t_srs EPSG:$epsg -dstnodata -9999 -tr 30 30 {} {.}_${epsg}.tif; dem_geoid --threads 1 --reverse-adjustment {.}_${epsg}.tif" ::: $fn_list
gdalbuildvrt -resolution highest -vrtnodata -9999 ned1_2003_adj.vrt *adj.tif
gdaladdo_ro.sh ned1_2003_adj.vrt

#Split the meta
ogr2ogr meta0306_PAL.shp meta0306/ PAL

#At this point, had to manually clean up
#Joined with tiles_1arcsec_meta_merge.shp based on quadname, selecting only S_DATE
#Created new fields OLD_S_DATE, NEW_S_DATE and computed difference
#Created new field S_DATE_CLN with clean dates

gdal_rasterize -co TILED=YES -co COMPRESS=LZW -tr 30 30 -a S_DATE_CLN -a_nodata 0 -ot Int16 meta0306_PAL_24k_10kmbuffer_clean_dissolve_32611.shp meta0306_PAL_24k_10kmbuffer_clean_dissolve_32611.tif

shp=meta0306_PAL_24k_10kmbuffer_clean_dissolve_32611.shp
mkdir s_date_shp
d_list=$(ogrinfo -al meta0306_PAL_24k_10kmbuffer_clean_dissolve_32611.shp | grep S_DATE_CLN | awk '{print $NF}' | sed '1d' | sort -u)
for d in $d_list 
do
    ogr2ogr -where "S_DATE_CLN=$d" s_date_shp/${d}_s_date.shp $shp
    #Make sure clip_raster_by_shp extent is set to shp
    clip_raster_by_shp.sh ned1_2003_adj.vrt s_date_shp/${d}_s_date.shp 
    mv ned1_2003_adj_shpclip.tif s_date_shp/${d}_s_date.tif
done
