#! /bin/bash 

#qsub -I -q devel -lselect=1:model=bro,walltime=2:00:00

count=true
index=true
res=32
ncpu=$(cat /proc/cpuinfo | egrep "core id|physical id" | tr -d "\n" | sed s/physical/\\nphysical/g | grep -v ^$ | sort | uniq | wc -l)
threads=$((ncpu-1))
#threads=28
tilesize=100000
mos=~/src/Tools/dem_mosaic_validtiles.py

ts=`date +%Y%m%d`

out=hma_${res}m_mos_${ts}/hma_${res}m
out=hma_${res}m_mos_20170601/hma_${res}m

#hma
proj='+proj=aea +lat_1=25 +lat_2=47 +lat_0=36 +lon_0=85 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs '

echo "Identifying input DEMs"
list=$(ls */*/*00/dem*/*-DEM_${res}m.tif */*00/dem*/*-DEM_${res}m.tif)
#list=$(ls */*/*/dem*/*-DEM_8m_trans.tif | grep -v QB)
#parallel -j $threads "if [ ! -e {.}_aea.tif ] ; then gdalwarp -overwrite -r cubic -t_srs \"$proj\" -tr $res $res -dstnodata -9999 {} {.}_aea.tif; fi" ::: $list
#list=$(echo $list | sed "s/DEM_${res}m.tif/DEM_${res}m_aea.tif/g")
#list=$(ls */*/dem*/*-DEM_${res}m_aea.tif)

if [ ! -d $(dirname $out) ] ; then 
    mkdir $(dirname $out)
    lfs setstripe -c $threads $(dirname $out)
fi

echo $list | tr ' ' '\n' > ${out}_input_DEM_list.txt
echo $(echo $list | wc -w) input DEMs

if [ ! -e $out.vrt ] ; then
    $mos --threads $threads --tr $res --t_srs "$proj" --georef_tile_size=$tilesize -o $out $list
fi

if (( "$res" == "32" )) ; then
    #Should run these in parallel
    gdal_opt='-co TILED=YES -co COMPRESS=LZW -co BIGTIFF=YES'
    lowres=100
    if [ ! -e ${out}_${lowres}m.tif ] ; then 
        echo "Preparing lowres $lowres m mosaic"
        gdalwarp -tr $lowres $lowres $gdal_opt $out.vrt ${out}_${lowres}m.tif
        gdaladdo_ro.sh ${out}_${lowres}m.tif
        hs.sh ${out}_${lowres}m.tif
        gdaladdo_ro.sh ${out}_${lowres}m_az315.tif
    fi
    if $count ; then
        if [ ! -e ${out}_count.vrt ] ; then
            echo "Preparing countmap"
            $mos --stat count --threads $threads --tr $res --t_srs "$proj" --georef_tile_size=$tilesize -o $out $list
        fi
        if [ ! -e ${out}_count.vrt ] ; then
            echo "Preparing lowres $lowres m countmap"
            gdalwarp -tr $lowres $lowres $gdal_opt ${out}_count.vrt ${out}_count_${lowres}m.tif
            gdaladdo_ro.sh ${out}_count_${lowres}m.tif
        fi
    fi
    if $index ; then
        parallel -j $threads "if [ ! -e {.}_${lowres}m.tif ] ; then gdalwarp -overwrite -r cubic -t_srs \"$proj\" -tr $lowres $lowres -dstnodata -9999 {} {.}_${lowres}m.tif; fi" ::: $list
        #Need to fix raster2shp, doesn't like subdir in merge_fn
        #raster2shp.py -merge_fn ${out}_stripindex.shp $(echo $list | sed 's/DEM_32m.tif/DEM_32m_100m.tif/g') 
        raster2shp.py $(echo $list | sed "s/DEM_${res}m.tif/DEM_${res}m_100m.tif/g")
        for i in $list
        do
            echo rm $(echo $i | sed "s/DEM_${res}m.tif/DEM_${res}m_100m.{shp,shx,dbf,prj}/")
        done
        ogr2ogr -f KML merge.kml merge.shp
    fi
fi
