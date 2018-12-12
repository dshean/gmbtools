#! /bin/bash

#After running dem_post_parallel.pbs, which generates stacks for all rgi polygons

set -e

gdal_opt='-co COMPRESS=LZW -co TILED=YES -co BIGTIFF=IF_SAFER'

#ASTER
#topdir=/nobackupp8/deshean/hma/aster/dsm
#shp=dem_align_ASTER_round2_index_2000-2018_aea.shp
#dt_list="20000211 20000531 20090531 20180531"
#WV
#topdir=/nobackup/deshean/hma/dem_coreg
#shp=dem_align_WV_index_2007-2018_aea.shp
#dt_list="20080101 20180101"
#Combined
topdir=/nobackupp8/deshean/hma/combined_aster_wv
shp=dem_align_ASTER_WV_index_2000-2018_aea.shp
dt_list="20000211 20000531 20090531 20180531"

dt_list="$(echo $dt_list | tr ' ' '\n' | sort -n)"
last_dt=$(echo $dt_list | awk '{print $NF}')
echo $dt_list

mos_list+="$dt_list trend trend_filt nmad"
#mos_list+=" trend_shpclip intercept diff diff_shpclip"

echo "Processing the following layers:"
echo $mos_list

cd $topdir/${shp%.*}_stack

#Generate products for consistent date
echo "Generating list of npz (finished stacks)"
if [ ! -e ${shp%.*}_npz_fn_list.txt ] ; then
    find . -name '*.npz' > ${shp%.*}_npz_fn_list.txt
fi

echo "Running stack_interp.py for all npz"
cat ${shp%.*}_npz_fn_list.txt | parallel "if [ ! -e {.}_$last_dt.tif ] ; then ~/src/gmbtools/gmbtools/stack_interp.py {} "$dt_list" ; fi"
#Overwrite
#cat ${shp%.*}_npz_fn_list.txt | parallel --progress "~/src/gmbtools/gmbtools/stack_interp.py {} "$dt_list""
#Parallel (now that filtering is implemented)
#qsub ~/src/gmbtools/stack_interp_parallel.pbs

#difference SRTM and 20000211 grid
##parallel --workdir $topdir/${shp%.*}_stack --sshloginfile $PBS_NODEFILE "compute_diff.py {} $rpcdem" ::: [0-9]*/*20000211.tif
##find . -name '*diff.tif' | parallel --workdir $topdir/${shp%.*}_stack --sshloginfile $PBS_NODEFILE 'if [ ! -e {.}_shpclip.tif ] ; then ~/src/pygeotools/pygeotools/clip_raster_by_shp.py -extent raster {} rgi ; fi'

echo "Generating fn lists for each layer"
parallel "if [ ! -e ${shp%.*}_{}_fn_list.txt ] ; then find . -name *{}.tif > ${shp%.*}_{}_fn_list.txt ; fi" ::: $mos_list

echo "Generating vrt for each layer"
parallel "if [ ! -e ${shp%.*}_{}_mos.vrt ] ; then gdalbuildvrt -r cubic -input_file_list ${shp%.*}_{}_fn_list.txt ${shp%.*}_{}_mos.vrt; fi" ::: $mos_list

#Regenerate tiles, much more efficient than working with vrt containing 90K files
#Want to select tilesize so that we have a reasonable number of tiles, maybe ~1K
tilesize=2048
echo "Retiling vrts"
parallel "if [ ! -d ${shp%.*}_{}_mos_retile ] ; then mkdir ${shp%.*}_{}_mos_retile; fi; gdal_retile.py -r cubic $gdal_opt -ps $tilesize $tilesize -targetDir ${shp%.*}_{}_mos_retile ${shp%.*}_{}_mos.vrt; " ::: $mos_list

echo "Generating new vrt for retile"
#Set nodata - gdal_retile.py doesn't propagate
#ndv=-9999
ndv=$(gdalinfo ${shp%.*}_${last_dt}_mos.vrt | grep NoData | awk -F'=' '{print $NF}')
#Shouldn't need -r cubic here, as these are regular grid
parallel "gdalbuildvrt -srcnodata $ndv -vrtnodata $ndv ${shp%.*}_{}_mos_retile.vrt ${shp%.*}_{}_mos_retile/*.tif" ::: $mos_list 

echo "Converting to tif"
parallel "gdalwarp $gdal_opt ${shp%.*}_{}_mos_retile.vrt ${shp%.*}_{}_mos_retile.tif; gdaladdo_ro.sh ${shp%.*}_{}_mos_retile.tif" ::: trend nmad 

echo "Extracting subsampled version"
parallel "gdal_translate $gdal_opt -outsize 25% 25% ${shp%.*}_{}_mos_retile.tif ${shp%.*}_{}_mos_retile_sub4.tif" ::: trend nmad 

#Now mb_parallel.py
