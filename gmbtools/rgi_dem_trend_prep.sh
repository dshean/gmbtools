#! /bin/bash

#Prepare to generate stacks for each RGI polygon using rgi_dem_trend.py
#Run after dem_post_parallel.pbs to align, then dem_align_post.py

#TODO
#Switch from shp to gpkg

#prefix=dem_align_aster
prefix=dem_align_noqb

#Throw out outliers
mkdir bad_align; for i in $(cat ${prefix}_bad_fn.txt); do mv -v $(echo $i | awk -F'/' '{print $1 "/" $2}') bad_align/; done

#Organize into annual subdir, important for many files and shell list length limitations
#valid_years=$(seq 2000 2018)
valid_years=$(cat ${prefix}_good_fn.txt | cut -c 1-4 | sort -u)
yr1=$(echo $valid_years | awk '{print $1}')
yr2=$(echo $valid_years | awk '{print $NF}')
#yr_mid=2009
#for y in $valid_years ; do if [ ! -d $y ] ; then mkdir $y ; fi ; mv AST_${y}* $y/ done
for y in $valid_years ; do if [ ! -d years/$y ] ; then mkdir -pv years/$y ; fi ; cd years/$y ; for i in ../../${y}*align/*align.tif; do ln -s $i . ; done; cd ../../ ; done

#Generate DEM index
#parallel 'gdaltindex -t_srs EPSG:4326 {}/${prefix}_index_{}.shp {}/*align/*align.tif' ::: $valid_years
parallel "gdaltindex -t_srs EPSG:4326 years/{}/${prefix}_index_{}.shp years/{}/*align.tif" ::: $valid_years
ogr_merge.sh ${prefix}_index_$yr1-$yr2.shp years/2*/${prefix}_index_*.shp
#ogr_merge.sh ${prefix}_index_$yr1-$yr_mid.shp 2*/${prefix}_index_200[0-9].shp
#ogr_merge.sh ${prefix}_index_2009-2018.shp 2*/${prefix}_index_2009.shp 2*/${prefix}_index_201[0-9].shp

#Now convert to aea projection
#shp_list="${prefix}_index_2000-2018.shp ${prefix}_index_2000-2009.shp ${prefix}_index_2009-2018.shp"
shp_list="${prefix}_index_$yr1-$yr2.shp"
parallel "ogr2ogr -t_srs '+proj=aea +lat_1=25 +lat_2=47 +lat_0=36 +lon_0=85 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ' {.}_aea.shp {}" ::: $shp_list 

#Now run rgi_dem_trend.py
#Now run dem_post_parallel.pbs
