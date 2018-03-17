#! /bin/bash

#Run on tpfe1 or `qsub -I -q devel -lselect=1:model=bro,walltime=2:00:00`
topdir=/nobackupp8/deshean/conus_combined
#topdir=/nobackupp8/deshean/hma

cd $topdir
#make_mos.sh 
#Get latest strip index
strip_shp=$(ls -tr mos/latest/*mos_32m/*mos_32m_stripindex_simp.shp | tail -1)

#Polygons to rerun
#site_shp=shp/rainier_rerun_wgs84.shp
#site_shp=sites/khumbu/khumbu_clip.shp
#site_shp=shp/hma_sites_20171004.shp
site_shp=shp/conus_sites_20180125_wgs84.shp

#This was a quick workaround to site_query.py
#ogr2ogr -clipsrc $site_shp ${site_shp%.*}_stripindex.shp $strip_shp
#ogrinfo -al ${site_shp%.*}_stripindex.shp | grep name | awk '{print $NF}' | awk -F'/' '{print $3}' > ${site_shp%.*}_stripindex_pairlist.txt

#This will create lists of overlapping files, stacks and mosaics
#site_query.py $site_shp $strip_shp

cd sites
for site in $(ls -d *)
do
    echo $site
    dem=$(ls $site/stack_all/${site}_stack_all-tile-0.tif)
    #Might need to update site_rerun_demproc.sh to handle relative paths
    #dem_rerun=$(ls $site/stack_all/${site}_stack_all-tile-0_dzfilt_0.00-*.00_gaussfill-tile-0.tif 1> /dev/null 2>&1)
    dem_rerun=$(ls $site/stack_all/${site}_stack_all-tile-0_dzfilt_0.00-*.00_gaussfill-tile-0.tif)
    #if $dem_rerun ; then 
    if [ ! -e $dem_rerun ] ; then 
        echo "Running site_rerun_demprep.sh to create filled DEM for rerun"
        #Pad crop_extent
        dem_rerun=$(site_rerun_demprep.sh $dem | tail -1)
    fi
    dem_fnlist=$(cat $site/site_${site}.csv | grep tif)

    #This creates subdir for alongtrack/crosstrack and pairname (HMA)
    #dem_pairlist=$(cat $site/site_${site}.csv | grep tif | awk -F'/' '{print $5 "/" $6}')

    #This creates subdir for pair name (CONUS)
    dem_pairlist=$(cat $site/site_${site}.csv | grep tif | awk -F'/' '{print $5}')

    if [ ! -d $site/rerun ] ; then
        mkdir $site/rerun
    fi
    pushd $site/rerun
    for pair in $dem_pairlist
    do
        if [ ! -d $pair ] ; then 
            mkdir -p $pair
        fi
        #HMA
        #r100_fn_list=$(ls $topdir/$pair/*r100.{tif,xml})
        #CONUS
        r100_fn_list=$(ls ../../../$pair/1*r100.{tif,xml})
        for i in $r100_fn_list
        do
            ln -sv ../$i $pair/ 
        done
    done
    pushd 

    #Run site_rerun_lfetransfer.sh to stage ntf/xml and/or r100.tif/xml
    #ssh lfe
    #/nobackup/deshean/src/gmbtools/gmbtools/site_rerun_lfetransfer.sh /nobackup/deshean/hma/sites2/$site/rerun

    #Currently, need to edit new_dg_stereo.sh, singlepair.pbs, dg_stereo_qsub.sh
    #Update rpcdem with filled dem_rerun product
    #Update crop_extent (~/src/demtools/get_extent.py $dem_rerun)

    #dg_stereo_qsub.sh

    #When complete

    #Manually remove bad DEMs

    #Create weighted-average mosaic to use during dz filter
    #Include NED/SRTM for sanity
    ref=/nobackup/deshean/rpcdem/ned13/ned13_tiles_glac24k_115kmbuff.vrt 
    lidar_ref=/nobackup/deshean/rpcdem/lidar/conus_lidar_1m.vrt

    dem_mosaic -o $site/rerun/stack_32m/mos_32m $site/rerun/*00/dem*/*DEM_32m.tif
    dem_mosaic --first-dem-as-reference -o $site/rerun/stack_32m/mos_32m_ned $site/rerun/stack_32m/mos_32m-tile-0.tif $site/rerun/stack_32m/mos_32m-tile-0.tif $site/rerun/stack_32m/mos_32m-tile-0.tif $ref

    #Filter clearly bogus pixels
    #ref=$site/rerun/stack_32m/*med.tif
    ref=$site/rerun/stack_32m/mos_32m_ned-tile-0.tif
    parallel --delay 0.5 "filter.py {} -filt dz -param $ref -200 200" ::: $site/rerun/*00/dem*/*DEM_32m.tif
    parallel --delay 0.5 "filter.py {} -filt dz -param $ref -200 200" ::: $site/rerun/*00/dem*/*DEM_8m.tif
    parallel --delay 0.5 "filter.py {} -filt dz -param $ref -200 200" ::: $site/rerun/*00/dem*/*DEM_2m.tif

    #Create preliminary 32m stack and gallery
    res=32
    #stack_dir=stack_${res}m
    #ext=DEM_${res}m.tif
    stack_dir=stack_${res}m_filt
    ext="DEM_${res}m_dzfilt_-200_200.tif"
    make_stack.py -o $site/rerun/$stack_dir --med $site/rerun/*00/dem*/*$ext
    ~/src/demtools/stack_extract.py $site/rerun/$stack_dir/*npz
    pushd $site/rerun/$stack_dir/*extract; dem_gallery.py *$ext ; mv dem_gallery.png ${site}_dem_gallery.png; pushd

    #ll -S *00/dem*/*DEM_32m.tif
    #mkdir bad

    #Prepare TOA
    #parallel --jobs 16 --delay 2 'toa.sh {}' ::: *00
    parallel --jobs 16 --delay 2 'toa.sh {}' ::: */rerun/*00

    #Prepare ref

    res=8
    re=''
    ext="*${re}*DEM_${res}m_dzfilt_-200_200.tif"

    list=$(ls $site/rerun/*00/dem*/$ext)
    parallel 'dem_mask.py --filter='none' --no_icemask --toa {}' ::: $list
    #Note: some of these will be empty, so no _ref.tif output
    dem_mosaic -o $site/rerun/mos_${res}m_all_toamask $(echo $list | sed 's/.tif/_ref.tif/g' | xargs ls)
    dem_mosaic --count -o $site/rerun/mos_${res}m_all_toamask $(echo $list | sed 's/.tif/_ref.tif/g' | xargs ls)

    dem_mosaic -o $site/rerun/mos_${res}m_all $list
    #Make summer-only stack to use as ref
    #Aug through Oct
    re='20[01][0-9][01][089][0-9][0-9]_'
    #June through Oct
    #re='[01][06789][0-9][0-9]'
    list=$(ls $site/rerun/*00/dem*/*${re}*DEM_${res}m.tif)
    dem_mosaic -o $site/rerun/mos_${res}m_Aug-Oct $list
    #Embed lidar?
    dem_mosaic -o $site/rerun/mos_${res}m_lidar --first $lidar_ref $site/rerun/mos_${res}m_Aug-Oct-tile-0.tif $site/rerun/mos_${res}m_all-tile-0.tif

    #coregister

    parallel --delay 2 "warptool.py -outdir {//} -te {} -tr {} -t_srs {} $lidar_ref" ::: */rerun/stack_32m/mos_32m-tile-0.tif 

    #dem_coreg_all.sh

    #make_stack.py -o stack_32m_trans $site/*00/dem*/*DEM_32m_trans.tif
    #dem_gallery.py $site/*00/dem*/*DEM_32m_trans.tif

    #Create 32m stack, identify outliers
    #anomaly_maps.py *npz
   
    #make_stack.py -o stack_2m_clean *00/dem*/*align/*DEM.tif

    #Generate new orthoimages
    #parallel --jobs 14 --delay 1 --verbose --progress 'ortho_proc.sh {}' ::: *00
    echo
done
