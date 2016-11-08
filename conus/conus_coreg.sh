#! /bin/bash

#David Shean
#dshean@gmail.com

#Prep with toa.sh and mask_lulc.py
#Run on bro node
#qsub ~/bin/devel.pbs
#cd /nobackup/deshean/conus5
#set pc_align_wrapper threads to 2
#parallel -j 16 --delay 3 '~/src/demtools/conus/conus_coreg.sh {}' ::: */dem*/*DEM_2m.tif

#Clean up past runs
#rm -r */*/*align */*/ned*eul.tif */*/*trans.tif

#Need to create vrt with 1 arcsec over areas where 1/3 is not avail

#Use 1 arcsec NED (30 m)
neddir=/nobackup/deshean/rpcdem/ned1
ned=$neddir/ned1_tiles_glac24k_115kmbuff.vrt
nedres=30

#Use 1/3 arcsec NED (10 m)
neddir=/nobackup/deshean/rpcdem/ned13
ned=$neddir/ned13_tiles_glac24k_115kmbuff.vrt
nedres=10

dem=$1
demdir=$(dirname $dem)
dembase=$(basename $dem)
outdir=${dembase%.*}_grid_align
dembase=$(echo $dembase | awk -F'-' '{print $1}')

#Use lowres mask
demref=$(ls $demdir/*DEM_32m_ref.tif)

warptool.py -te $dem -tr $nedres -t_srs $dem -outdir $demdir $ned

demned=$demdir/$(basename $ned)
demned=${demned%.*}_warp.tif
#Mask the NED over valid pixels
apply_mask_new.py $demned $demref
demned_masked=${demned%.*}_masked.tif

#point-to-point
pc_align_wrapper.sh $demned_masked $dem

cd $demdir
log=$(ls $outdir/*log)
if [ -e $log ] ; then 
    apply_dem_translation.py ${dembase}-DEM_32m.tif $log
    apply_dem_translation.py ${dembase}-DEM_8m.tif $log
    ln -s $outdir/*DEM.tif ${dembase}-DEM_2m_trans.tif
    compute_dh.py $(basename $demned) ${dembase}-DEM_8m_trans.tif
fi

