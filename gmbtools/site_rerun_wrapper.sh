#! /bin/bash

#On lfe - Khumbu rerun, 9/12/17

#cd ~/hma
#nbdir=/nobackupp8/deshean/hma/sites/khumbu/rerun
#pushd $nbdir
#pairlist=$(ls -d *00)
#pushd
#for pair in $pairlist; do echo $pair; p=$(find 2* -name $pair) ; shiftc -R $p $nbdir ; done
#Crosstrack
#pushd $nbdir/crosstrack
#pairlist=$(ls -d *00)
#for pair in $pairlist ; do id1=$(echo $pair | awk -F'_' '{print $3}') ; id2=$(echo $pair | awk -F'_' '{print $4}') ; echo $pair ; shiftc -R $id1.r100* $id2.r100* /nobackupp8/deshean/hma/sites/khumbu/rerun/crosstrack/$pair/ ; done

#Create list of IDs to reprocess

#This creates dir of links, avoids find, can just find pairname directly
#ids=$(cat good_list_dem.txt  | awk -F'_' '{print $3}')
#for id in $ids; do pair=$(ls -d /nobackupp8/deshean/conus_dir/*${id}* | awk -F'/' '{print $NF}'); echo $pair >> pairlist; done
#shiftc -L -r -d --include 'r100.tif' --include 'r100.xml' $(cat pairlist | sed 's#^#/nobackup/deshean/conus_dir/#')

#for i in $idlist
#do
#    i=$(find conus[1-5] -name "$id.r100.tif")
#    #i=conus2/WV01_20151014_1020010043334500_1020010044474A00
#    rsync -av --include='*/' --include='*r100*' --exclude='*' $i /nobackupp8/deshean/conus/scg_rerun/
#done

#Run on tpfe1 or `qsub -I -q devel -lselect=1:model=bro,walltime=2:00:00`
#topdir=/nobackupp8/deshean/conus_combined
topdir=/nobackupp8/deshean/hma
cd $topdir
#make_mos.sh 
#Get latest strip index
#strip_shp=$(ls -tr mos/latest/mos_32m/*mos_32m_stripindex_simp.shp | tail -1)
strip_shp=$(ls -tr mos/*_mos/mos_32m/*mos_32m_stripindex_simp.shp | tail -1)

#Polygons to rerun
#site_shp=shp/rainier_rerun_wgs84.shp
#site_shp=sites/khumbu/khumbu_clip.shp
site_shp=shp/hma_sites_20171004.shp

#This was a quick workaround to site_query.py
#ogr2ogr -clipsrc $site_shp ${site_shp%.*}_stripindex.shp $strip_shp
#ogrinfo -al ${site_shp%.*}_stripindex.shp | grep name | awk '{print $NF}' | awk -F'/' '{print $3}' > ${site_shp%.*}_stripindex_pairlist.txt

#This will create lists of overlapping files, stacks and mosaics
#site_query.py $site_shp $strip_shp

cd sites2
for site in $(ls -d *)
do
    echo $site
    dem=$(ls $site/stack_all/${site}_stack_all-tile-0.tif)
    #Might need to update site_rerun_demproc.sh to handle relative paths
    dem_rerun=$(ls $site/stack_all/${site}_stack_all-tile-0_dzfilt_0.00-*.00_gaussfill-tile-0.tif 1> /dev/null 2>&1)
    #if [ ! -e $dem_rerun ] ; then 
    if $dem_rerun ; then 
        echo "Running site_rerun_demprep.sh to create filled DEM for rerun"
        dem_rerun=$(site_rerun_demprep.sh $dem | tail -1)
    fi
    dem_fnlist=$(cat $site/site_${site}.csv | grep tif)
    #dem_pairlist=$(cat $site/site_${site}.csv | grep tif | awk -F'/' '{print $5}')
    dem_pairlist=$(cat $site/site_${site}.csv | grep tif | awk -F'/' '{print $5 "/" $6}')

    if [ ! -d $site/rerun ] ; then
        mkdir $site/rerun
    fi
    pushd $site/rerun
    for pair in $dem_pairlist
    do
        mkdir -p $pair
        #HMA
        #r100_fn_list=$(ls $topdir/$pair/*r100.{tif,xml})
        #CONUS
        #r100_fn_list=$(ls ../../../$pair/1*r100.{tif,xml})
        #for i in $r100_fn_list
        #do
        #    ln -s ../$i $pair/ 
        #done
    done
    pushd 

    #Currently, need to edit new_dg_stereo.sh
    #Update rpcdem
    #Update crop_extent (~/src/demtools/get_extent.py $dem_rerun)

    #dg_stereo_qsub.sh

    #When complete
    #coregister

    #dem_gallery.py $site/*00/dem*/*DEM_32m.tif
    #make_stack.py $site/*00/dem*/*DEM_32m.tif

    #Generate new orthoimages
    #parallel --jobs 14 --delay 1 --verbose --progress 'ortho_proc.sh {}' ::: *00
    echo
done
