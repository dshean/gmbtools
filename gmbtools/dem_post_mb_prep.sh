#! /bin/bash

#After running dem_post_parallel.pbs
#These take <1s, inefficient for sshd on large block, just check out single node and chug
#pbs_rfe or
#qsub -I -q devel -lselect=1:model=bro

#make_mos.sh

gdal_opt='-co COMPRESS=LZW -co TILED=YES -co BIGTIFF=IF_SAFER'

#shp=aster_align_index_2009-2018_aea.shp

#Generate products for consistent dates
#parallel '~/src/gmbtools/gmbtools/stack_interp.py {}' ::: */*.npz
#find . -name '*.npz' > npz_fn_list.txt
#parallel --arg-file npz_fn_list.txt 'if [ ! -e {.}_20180531.tif ] ; then ~/src/gmbtools/gmbtools/stack_interp.py {}; fi'

#difference SRTM and 20000211 grid
#parallel --workdir $topdir/${shp%.*}_stack --sshloginfile $PBS_NODEFILE "compute_dz.py {} $rpcdem" ::: [0-9]*/*20000211.tif

#find . -name '*eul.tif' | parallel --workdir $topdir/${shp%.*}_stack --sshloginfile $PBS_NODEFILE 'if [ ! -e {.}_shpclip.tif ] ; then ~/src/pygeotools/pygeotools/clip_raster_by_shp.py -extent raster {} rgi ; fi'

#Build vrt mosaics
#parallel "gdalbuildvrt ${shp%.*}_stack/${shp%.*}_mos_{}.vrt ${shp%.*}_stack/[0-9]*/*{}.tif" ::: 20000211 20000531 20090531 20180531 nmad trend trend_shpclip intercept eul eul_shpclip
#cd ${shp%.*}_stack
#parallel "find . -name *{}.tif > {}_fn_list.txt; gdalbuildvrt -input_file_list {}_fn_list.txt ${shp%.*}_mos_{}.vrt" ::: 20000211 20000531 20090531 20180531 nmad trend 

#Export tif
#Regenerate tiles, much more efficient than large vrt
#parallel "mkdir {.}; gdal_retile.py -r cubic $gdal_opt -ps 2048 2048 -targetDir {.} {}" ::: *vrt
#ndv=-9999
#parallel "gdalbuildvrt -srcnodata $ndv -vrtnodata $ndv {.}_tiled.vrt ${shp%.*}_mos_{}/*.tif" :::  20000211  20000531 20090531 20180531 nmad trend
#Set nodata
#parallel --workdir $topdir/${shp%.*}_stack --sshloginfile $PBS_NODEFILE "gdalwarp -r cubic $gdal_opt ${shp%.*}_mos_{}.vrt ${shp%.*}_mos_trend.tif; gdaladdo_ro.sh ${shp%.*}_mos_{}.tif" ::: trend trend_shpclip eul eul_shpclip

#Now mb_parallel.py
