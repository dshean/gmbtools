#! /bin/bash

cd srtm1

#Define projection
#epsg=32611
epsg=32644

urllist="hma_srtm_gl1_url.txt"

#Modify lat/lon bounds and generate with the get_srtm_tilelist.py script
get_srtm_tilelist.py >> hma_srtm_gl1_url.txt

#Old approach, from an existing list
#list='/scr2/ned/NED_1arcsec_grid_intersect_24k_115kmbuff_tilelist.txt'
#echo -n > $urllist
#for i in $(cat $list | tr '[:lower:]' '[:upper:]')
#do
#    echo https://cloud.sdsc.edu/v1/AUTH_opentopography/Raster/SRTM_GL1/SRTM_GL1_srtm/$i.hgt >> $urllist
#done

wget -nc -i $urllist

fn_list=$(ls *hgt)

#NOTE: SRTM GL1 are relative to EGM96
#https://lta.cr.usgs.gov/SRTM1Arc
#February 11-22, 2000

#Convert and adjust datum
parallel -j 16 "gdalwarp -co COMPRESS=LZW -co TILED=YES -co BIGTIFF=IF_SAFER -overwrite -r cubic -t_srs EPSG:$epsg -dstnodata -32768 -tr 30 30 {} {.}_${epsg}.tif; dem_geoid --threads 1 --reverse-adjustment {.}_${epsg}.tif" ::: $fn_list

#Create vrt mosaic
gdalbuildvrt -resolution highest -vrtnodata -32768 hma_srtm_gl1.vrt *_${epsg}-adj.tif
