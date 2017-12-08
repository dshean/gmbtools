#! /bin/bash

#ssh lfe

#On lfe - Khumbu rerun, 9/12/17
#nbdir=/nobackupp8/deshean/hma/sites/khumbu/rerun
nbdir=$1

cd ~/hma
#Alongtrack
if [ -d $nbdir/alongtrack ] ; then
    pushd $nbdir/alongtrack
    pairlist=$(ls -d *00)
    pushd
    for pair in $pairlist; do echo $pair; p=$(find 2* -name $pair) ; shiftc -R $p $nbdir/alongtrack ; done
fi
#Crosstrack
if [ -d $nbdir/crosstrack ] ; then
    pushd $nbdir/crosstrack
    pairlist=$(ls -d *00)
    pushd ~/hma/mono/r100_mono
    for pair in $pairlist ; do id1=$(echo $pair | awk -F'_' '{print $3}') ; id2=$(echo $pair | awk -F'_' '{print $4}') ; echo $pair ; shiftc -R $id1.r100* $id2.r100* $nbdir/crosstrack/$pair/ ; done
fi

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

