#! /bin/bash

#This was extracted from make_mos.sh
#At some point, clean up so make_mos.sh just calls this modular script

ncpu=$(cat /proc/cpuinfo | egrep "core id|physical id" | tr -d "\n" | sed s/physical/\\nphysical/g | grep -v ^$ | sort | uniq | wc -l)
threads=$((ncpu-1))
res=32
lowres=100
tol=0.001

#site=hma
#proj='+proj=aea +lat_1=25 +lat_2=47 +lat_0=36 +lon_0=85 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs '
site=conus
proj='+proj=aea +lat_1=36 +lat_2=49 +lat_0=43 +lon_0=-115 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs '

ts=`date +%Y%m%d`

#out=mos/${site}_${ts}_mos/mos_${res}m/${site}_${ts}_mos_${res}m
out=mos/mos_2007-2010/mos_2007-2010

if [ ! -d $(dirname $out) ] ; then 
    mkdir -p $(dirname $out)
    lfs setstripe -c $threads $(dirname $out)
fi

#list=$(ls *track/*00/dem*/*-DEM_${res}m.tif)
list="$@"
#list=$(echo $list | tr ' ' '\n' | sort -n -t'/' -k 3)
#list=$(echo $list | tr ' ' '\n' | sort -n -r -t'/' -k 7)

echo $list | tr ' ' '\n' > ${out}_input_DEM_list.txt
echo $(echo $list | wc -w) input DEMs

parallel -j $threads "if [ ! -e {.}_${lowres}m.tif ] ; then gdalwarp -overwrite -r cubic -t_srs \"$proj\" -tr $lowres $lowres -dstnodata -9999 {} {.}_${lowres}m.tif; fi" ::: $list
raster2shp.py -merge_fn ${out}_stripindex.shp $(echo $list | sed "s/.tif/_${lowres}m.tif/g") 
#echo "Removing intermediate files"
#eval rm $(echo $list | sed "s/DEM_${res}m.tif/DEM_${res}m_${lowres}m.{shp,shx,dbf,prj,tif}/")
ogr2ogr -simplify $tol ${out}_stripindex_simp.shp ${out}_stripindex.shp
rm ${out}_stripindex.{shp,dbf,prj,shx}
if [ ! -e ${out}_stripindex_simp.kml ] ; then 
    ogr2ogr -f KML ${out}_stripindex_simp.kml ${out}_stripindex_simp.shp
fi
