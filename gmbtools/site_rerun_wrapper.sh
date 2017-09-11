#! /bin/bash

#Run on tpfe1 or `qsub -I -q devel -lselect=1:model=bro,walltime=2:00:00`
topdir=/nobackupp8/deshean/conus_combined
cd $topdir
#make_mos.sh 
#Get latest strip index
strip_shp=$(ls -tr mos/*/*32m/*stripindex_simp.shp | tail -1)

#Polygons to rerun
site_shp=shp/rainier_rerun_wgs84.shp
#conus_site_poly.sh $site_shp $strip_shp

cd sites
for site in $(ls -d *)
do
    dem=$(ls $site/stack_all/${site}_stack_all-tile-0.tif)
    #Might need to update site_rerun_demproc.sh to handle relative paths
    #dem_rerun=$(site_rerun_demproc.sh $dem | tail -1)
    dem_rerun=$site/stack_all/${site}_stack_all-tile-0_dzfilt_0.00-100.00_gaussfill-tile-0.tif
    dem_fnlist=$(cat $site/site_${site}.csv | grep tif)
    dem_pairlist=$(cat $site/site_${site}.csv | grep tif | awk -F'/' '{print $5}')

    if [ ! -d rerun ] ; then
        mkdir rerun
    fi
    cd rerun
    for pair in $dem_pairlist
    do
        mkdir $pair
        #r100_fn_list=$(ls $topdir/$pair/*r100.{tif,xml})
        r100_fn_list=$(ls ../../../$pair/1*r100.{tif,xml})
        for i in $r100_fn_list
        do
            ln -s ../$i $pair/ 
        done
    done
    #Currently, need to edit new_dg_stereo.sh
    #Update rpcdem
    #Update crop_extent (~/src/demtools/get_extent.py $dem_rerun)

done

exit

#Create list of IDs to reprocess

#This creates dir of links, avoids find, can just find pairname directly
#ids=$(cat good_list_dem.txt  | awk -F'_' '{print $3}')
#for id in $ids; do pair=$(ls -d /nobackupp8/deshean/conus_dir/*${id}* | awk -F'/' '{print $NF}'); echo $pair >> pairlist; done
#shiftc -L -r -d --include 'r100.tif' --include 'r100.xml' $(cat pairlist | sed 's#^#/nobackup/deshean/conus_dir/#')

for i in $idlist
do
    i=$(find conus[1-5] -name "$id.r100.tif")
    #i=conus2/WV01_20151014_1020010043334500_1020010044474A00
    rsync -av --include='*/' --include='*r100*' --exclude='*' $i /nobackupp8/deshean/conus/scg_rerun/
done

#Create rpc DEM, burning in summer 2008
dem_mosaic --priority-blending-length 5 /nobackup/deshean/conus/dem2/noca_mora_site_rerun/dem_rerun_coreg/site_poly_highcount_rect3_rerun/scg/scg_2008_summer_n1_20081001-20081001-tile-0.tif scg_2012-2016_8m_trans_mos-tile-0.tif -o scg_2012-2016_8m_trans_mos_burn_2008

qsub -I -q devel -lselect=3:model=bro,walltime=2:00:00
