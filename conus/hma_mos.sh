#! /bin/bash 

#qsub -I -q devel -lselect=1:model=bro,walltime=2:00:00

res=8
threads=28
tilesize=100000

out=hma_${res}m_mos_20170410/hma_${res}m

#hma
proj='+proj=aea +lat_1=25 +lat_2=47 +lat_0=36 +lon_0=85 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs '

list=$(ls */*/*00/dem*/*-DEM_${res}m.tif */*00/dem*/*-DEM_${res}m.tif)
echo $list | tr ' ' '\n' > ${out}_input_DEM_list.txt
#list=$(ls */*/*/dem*/*-DEM_8m_trans.tif | grep -v QB)
#parallel -j $threads "if [ ! -e {.}_aea.tif ] ; then gdalwarp -overwrite -r cubic -t_srs \"$proj\" -tr $res $res -dstnodata -9999 {} {.}_aea.tif; fi" ::: $list
#list=$(echo $list | sed "s/DEM_${res}m.tif/DEM_${res}m_aea.tif/g")
#list=$(ls */*/dem*/*-DEM_${res}m_aea.tif)

mkdir $(dirname $out); lfs setstripe -c $threads $(dirname $out)

~/src/Tools/dem_mosaic_validtiles.py --threads $threads --tr $res --t_srs "$proj" --georef_tile_size=$tilesize -o $out $list
#~/src/Tools/dem_mosaic_validtiles.py --stat count --threads $threads --tr $res --t_srs "$proj" --georef_tile_size=$tilesize -o $out $list

if (( "$res" == "32" )) ; then
    #Should run these in parallel
    gdal_opt='-co TILED=YES -co COMPRESS=LZW -co BIGTIFF=YES'
    lowres=100
    #gdal_translate $gdal_opt $out.vrt $out.tif
    #gdaladdo_ro.sh $out.tif
    #gdal_translate $gdal_opt ${out}_count.vrt ${out}_count.tif
    #gdaladdo_ro.sh ${out}_count.tif
    gdalwarp -tr $lowres $lowres $gdal_opt $out.vrt ${out}_${lowres}m.tif
    gdaladdo_ro.sh ${out}_${lowres}m.tif
    gdalwarp -tr $lowres $lowres $gdal_opt ${out}_count.vrt ${out}_count_${lowres}m.tif
    gdaladdo_ro.sh ${out}_count_${lowres}m.tif

    parallel -j $threads "if [ ! -e {.}_${lowres}m.tif ] ; then gdalwarp -overwrite -r cubic -t_srs \"$proj\" -tr $lowres $lowres -dstnodata -9999 {} {.}_${lowres}m.tif; fi" ::: $list
    #Need to fix raster2shp, doesn't like subdir in merge_fn
    #raster2shp.py -merge_fn ${out}_stripindex.shp $(echo $list | sed 's/DEM_32m.tif/DEM_32m_100m.tif/g') 
    raster2shp.py $(echo $list | sed 's/DEM_32m.tif/DEM_32m_100m.tif/g') 
    ogr2ogr -f KML merge.kml merge.shp
fi
