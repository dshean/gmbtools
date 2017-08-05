#! /bin/bash 

#qsub -I -q devel -lselect=1:model=bro,walltime=2:00:00
#cd /nobackup/deshean/hma
#~/src/conus/conus/hma_mos.sh

#Set open file limit
#Default is 2048
ulimit -n 65536

res=32
#res=8
count=true
index=true
tileindex=true

ncpu=$(cat /proc/cpuinfo | egrep "core id|physical id" | tr -d "\n" | sed s/physical/\\nphysical/g | grep -v ^$ | sort | uniq | wc -l)
threads=$((ncpu-1))
tilesize=100000
mos=~/src/Tools/dem_mosaic_validtiles.py
#Simplify tolerance in decimal degrees
tol=0.001

ts=`date +%Y%m%d`

site=hma
proj='+proj=aea +lat_1=25 +lat_2=47 +lat_0=36 +lon_0=85 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs '

#site=conus
#proj='+proj=aea +lat_1=36 +lat_2=49 +lat_0=43 +lon_0=-115 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs '

out=mos/${site}_${ts}_mos/mos_${res}m/${site}_${ts}_mos_${res}m

#Should add option to split annually
echo "Identifying input DEMs"
#list=$(ls */*/*00/dem*/*-DEM_${res}m.tif */*00/dem*/*-DEM_${res}m.tif)
list=$(ls */stereo/*00/dem*/*-DEM_${res}m.tif */*00/dem*/*-DEM_${res}m.tif)
#list=$(ls */*/*/dem*/*-DEM_8m_trans.tif | grep -v QB)
#parallel -j $threads "if [ ! -e {.}_aea.tif ] ; then gdalwarp -overwrite -r cubic -t_srs \"$proj\" -tr $res $res -dstnodata -9999 {} {.}_aea.tif; fi" ::: $list
#list=$(echo $list | sed "s/DEM_${res}m.tif/DEM_${res}m_aea.tif/g")
#list=$(ls */*/dem*/*-DEM_${res}m_aea.tif)

if [ ! -d $(dirname $out) ] ; then 
    mkdir -p $(dirname $out)
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
        gdaladdo_ro.sh ${out}_${lowres}m_hs_az315.tif
    fi
    if $count ; then
        if [ ! -e ${out}_count.vrt ] ; then
            echo "Preparing countmap"
            $mos --stat count --threads $threads --tr $res --t_srs "$proj" --georef_tile_size=$tilesize -o $out $list
        fi
        if [ ! -e ${out}_count_${lowres}m.tif ] ; then
            echo "Preparing lowres $lowres m countmap"
            gdalwarp -tr $lowres $lowres $gdal_opt ${out}_count.vrt ${out}_count_${lowres}m.tif
            gdaladdo_ro.sh ${out}_count_${lowres}m.tif
        fi
    fi
    if $index ; then
        if [ ! -e ${out}_stripindex.shp ] ; then
            parallel -j $threads "if [ ! -e {.}_${lowres}m.tif ] ; then gdalwarp -overwrite -r cubic -t_srs \"$proj\" -tr $lowres $lowres -dstnodata -9999 {} {.}_${lowres}m.tif; fi" ::: $list
            raster2shp.py -merge_fn ${out}_stripindex.shp $(echo $list | sed "s/DEM_${res}m.tif/DEM_${res}m_${lowres}m.tif/g") 
            echo "Removing intermediate files"
            eval rm $(echo $list | sed "s/DEM_${res}m.tif/DEM_${res}m_${lowres}m.{shp,shx,dbf,prj,tif}/")
            ogr2ogr -simplify $tol ${out}_stripindex_simp.shp ${out}_stripindex.shp
            rm ${out}_stripindex.{shp,dbf,prj,shx}
            if [ ! -e ${out}_stripindex_simp.kml ] ; then 
                ogr2ogr -f KML ${out}_stripindex_simp.kml ${out}_stripindex_simp.shp
            fi
        fi
    fi
    #This is a hack for now, should add this functionality to dem_mosaic_validtiles.py
    if $tileindex ; then 
        raster2shp.py -merge_fn ${out}_tileindex.shp $(gdalinfo $out.vrt | grep tif)
        ogr2ogr -simplify $tol ${out}_tileindex_simp.shp ${out}_tileindex.shp
        echo "Removing intermediate files"
        rm ${out}_tileindex.{shp,dbf,prj,shx}
        eval rm $(gdalinfo $out.vrt | grep tif | sed 's/tif/{shp,dbf,prj,shx}/')
        if [ ! -e ${out}_tileindex_simp.kml ] ; then 
            ogr2ogr -f KML ${out}_tileindex_simp.kml ${out}_tileindex_simp.shp
        fi
    fi
fi

#Set permissions
chmod -R 755 mos
chmod 644 $out/*
