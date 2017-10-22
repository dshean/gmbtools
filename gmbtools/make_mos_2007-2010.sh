#! /bin/bash

res=8
proj='+proj=aea +lat_1=36 +lat_2=49 +lat_0=43 +lon_0=-115 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs '
out=mos/mos_2007-2010/mos_2007-2010_8m

wv_dir=/nobackupp8/deshean/conus_combined
wv_fn=$(ls $wv_dir/*_200[89]0[89]*/dem*/*DEM_${res}m_trans.tif $wv_dir/*_200[89]1[01]*/dem*/*DEM_${res}m_trans.tif)

#lidar_dir=/nobackupp8/deshean/rpcdem/lidar/lidar_32610_1m
#lidar_fn=$(ls $lidar_dir/*200[789]*tif)
#lidar_fn="20070315-20071107_hood_3m_unitmeters-adj_32610_1m.tif 20070901-20081001_rainier_1m-adj_32610_1m.tif 20091009_sisters_3ft_unitmeters-adj_32610_1m.tif 20100902_shastalidar1m-adj_32610_1m.tif"
#temp=""; for i in $lidar_fn ; do temp+=" $lidar_dir/$i" ; done; lidar_fn=$temp

lidar_dir=/nobackupp8/deshean/rpcdem/lidar/lidar_32610_8m
lidar_fn=$(ls $lidar_dir/*tif)

list="$lidar_fn $wv_fn"
list=$(echo $list | tr ' ' '\n' | sort -n -t'/' -k 7)
echo $list
echo

cd $wv_dir
#dem_mosaic_validtiles.py --t_srs "$proj" --tr $res -o $out $list
shp_list=$(echo $list | awk '{ for (i=NF; i>1; i--) printf("%s ",$i); print $1; }')
echo $shp_list
make_stripindex.sh $list
#Modify Hood timestamp to be 20070915
