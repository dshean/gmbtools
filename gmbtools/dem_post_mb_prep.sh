#! /bin/bash

#After running dem_post_parallel.pbs, which generates stacks for all rgi polygons

set -e

gdal_opt='-co COMPRESS=LZW -co TILED=YES -co BIGTIFF=IF_SAFER'

topdir=/nobackupp8/deshean/hma/aster/dsm
shp=dem_align_ASTER_round2_index_2000-2018_aea.shp

cd ${shp%.*}_stack

#Generate products for consistent date
echo "Generating list of npz (finished stacks)"
if [ ! -e ${shp%.*}_npz_fn_list.txt ] ; then
    find . -name '*.npz' > ${shp%.*}_npz_fn_list.txt
fi

echo "Running stack_interp.py for all npz"
parallel --arg-file ${shp%.*}_npz_fn_list.txt 'if [ ! -e {.}_20180531.tif ] ; then ~/src/gmbtools/gmbtools/stack_interp.py {}; fi'

#Should pass these into stack_interp.py
mos_list="20000211 20000531 20090531 20180531 nmad trend"
#mos_list+="trend_shpclip intercept diff diff_shpclip"
echo "Processing the following layers:"
echo $mos_list

#difference SRTM and 20000211 grid
##parallel --workdir $topdir/${shp%.*}_stack --sshloginfile $PBS_NODEFILE "compute_diff.py {} $rpcdem" ::: [0-9]*/*20000211.tif
##find . -name '*diff.tif' | parallel --workdir $topdir/${shp%.*}_stack --sshloginfile $PBS_NODEFILE 'if [ ! -e {.}_shpclip.tif ] ; then ~/src/pygeotools/pygeotools/clip_raster_by_shp.py -extent raster {} rgi ; fi'

echo "Generating fn lists for each layer"
parallel "if [ ! -e ${shp%.*}_{}_fn_list.txt ] ; then find . -name *{}.tif > ${shp%.*}_{}_fn_list.txt ; fi" ::: $mos_list

echo "Generating vrt for each layer"
parallel "if [ ! -e ${shp%.*}_{}_mos.vrt ] ; then gdalbuildvrt -r cubic -input_file_list ${shp%.*}_{}_fn_list.txt ${shp%.*}_{}_mos.vrt; fi" ::: $mos_list

#Regenerate tiles, much more efficient than working with vrt containing 90K files
echo "Retiling vrts"
parallel "mkdir ${shp%.*}_{}_mos_retile; gdal_retile.py -r cubic $gdal_opt -ps 2048 2048 -targetDir ${shp%.*}_{}_mos_retile ${shp%.*}_{}_mos.vrt; " ::: $mos_list

echo "Generating new vrt for retile"
#Set nodata
ndv=-9999
#Shouldn't need -r cubic here, as these are regular grid
parallel "gdalbuildvrt -srcnodata $ndv -vrtnodata $ndv ${shp%.*}_{}_mos_retile.vrt ${shp%.*}_{}_mos_retile/*.tif" ::: $mos_list 

echo "Converting to tif"
parallel "gdalwarp $gdal_opt ${shp%.*}_{}_mos_retile.vrt ${shp%.*}_{}_mos_retile.tif; gdaladdo_ro.sh ${shp%.*}_{}_mos_retile.tif" ::: trend nmad 

#Now mb_parallel.py
