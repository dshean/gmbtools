#! /bin/bash

#Download and wrangle NED tiles

#Create a list of necessary tiles (in this case, those that are within 115 km buffer from glacier polygons)

#shp=/Volumes/500GB_1/RGI/24kgrteqaul0.01km/24k_selection_32610_115kmbuff.shp
#ogrinfo -al NED_1arcsec_grid_intersect_24k_115kmbuff.shp | grep -i FILE_ID | awk '{print $NF}' | sort -u > NED_1arcsec_grid_intersect_24k_115kmbuff_tilelist.txt
list='/scr2/ned/NED_1arcsec_grid_intersect_24k_115kmbuff_tilelist.txt'

urllist=${list%.*}_url.txt
echo -n > $urllist

#Resolution (arcsec)
#r=1
r=13

tiledir=/scr2/ned/tiles/NED_${r}
cd $tiledir

for i in $(cat $list)
do
    if [ ! -d $tiledir/USGS_NED_${r}_${i}_ArcGrid ] ; then 
        echo https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/${r}/ArcGrid/USGS_NED_${r}_${i}_ArcGrid.zip >> $urllist
    fi
    if [ ! -d $i ] ; then
        echo https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/${r}/ArcGrid/${i}.zip >> $urllist
    fi
done

wget -nc -i $urllist

parallel -j 4 'unzip -d {.} {}' ::: *zip

mkdir duplicates
for i in $(cat $list)
do
    if [ -d USGS_NED_${r}_${i}_ArcGrid ] ; then
        mv $i duplicates
    fi
done

fn_list=$(ls -d */grd*${r})

#WGS84 UTM 11N
#epsg=32611
#NAD83 UTM 11N
epsg=26911

#parallel "gdalwarp -co COMPRESS=LZW -co TILED=YES -co BIGTIFF=IF_SAFER -overwrite -r cubic -t_srs EPSG:$epsg -dstnodata -9999 -tr 30 30 {} {.}_${epsg}.tif; dem_geoid --threads 1 --reverse-adjustment {.}_${epsg}.tif" ::: $fn_list
parallel "gdalwarp -co COMPRESS=LZW -co TILED=YES -co BIGTIFF=IF_SAFER -overwrite -r cubic -t_srs EPSG:$epsg -dstnodata -9999 -tr 10 10 {} {.}_${epsg}.tif; dem_geoid --threads 1 --reverse-adjustment {.}_${epsg}.tif" ::: $fn_list
#parallel "dem_geoid --reverse-adjustment {.}_${epsg}.tif" ::: $fn_list

gdalbuildvrt -resolution highest -vrtnodata -9999 ned${r}_tiles_glac24k_115kmbuff.vrt *adj.tif

