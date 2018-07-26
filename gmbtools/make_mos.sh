#! /bin/bash 

set -e

#pbs_rfe --duration 0+6 --model bro
#qsub -I -q devel -lselect=1:model=bro,walltime=2:00:00
#qsub -I -q long -lselect=1:model=bro,walltime=6:00:00
#cd /nobackup/deshean/hma
#make_mos.sh

#Want to make cascading mosaic
#gdalbuildvrt mos_8m_all mos_8m_all_trans mos_8m_summer_trans
#Note on ordering within vrt
#If there is some amount of spatial overlapping between files, the order of files appearing in the list of source matter: files that are listed at the end are the ones from which the content will be fetched. Note that nodata will be taken into account to potentially fetch data from less prioritary datasets

#Set open file limit
#Default is 2048
ulimit -n 65536

latest=false
summer=false
trans=false

#Output mosaic res
#res=2
res=8
#res=32

#Write out tif mosaic at full res
out_tif=false

count=false
stddev=false
last=true
first=true
index=false
tileindex=false

#Exclude Quickbird-2
noQB=true

#Hardcoded site and projection for now
site=hma
proj='+proj=aea +lat_1=25 +lat_2=47 +lat_0=36 +lon_0=85 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs '
#site=conus
#proj='+proj=aea +lat_1=36 +lat_2=49 +lat_0=43 +lon_0=-115 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs '
#site=fuego
#proj='EPSG:32615'

#Mosaic tile size in meters
#tilesize=20000
tilesize=100000
mos=~/src/Tools/dem_mosaic_validtiles.py
#Simplify tolerance in decimal degrees
tol=0.001

ncpu=$(cat /proc/cpuinfo | egrep "core id|physical id" | tr -d "\n" | sed s/physical/\\nphysical/g | grep -v ^$ | sort | uniq | wc -l)
threads=$((ncpu-1))

ts=`date +%Y%m%d`

#Should add option to split annually
echo "Identifying input DEMs"

ext="-DEM_${res}m"
mos_ext="${site}_mos_${res}m"

if $latest ; then
    #Note: CONUS is only current through June 2016, missing 
    re='201[5-9]'
    mos_ext+='_latest'
else
    #CONUS
    #Want to exclude 2008/2009 data
    #re='201[2-9]'
    #HMA
    #re='201[0-9]'
    re='20[0-9]'
fi
    
if $summer ; then
    #August through October
    re+='[01][089][0-9][0-9]'
    if ! $latest ; then
        #June through October
        re='[01][06789][0-9][0-9]'
    fi
    mos_ext+='_summer'
fi

if $trans ; then
    ext+='_trans'
    mos_ext+='_trans'
fi

#list=$(ls */*/*00/dem*/*-DEM_${res}m.tif */*00/dem*/*-DEM_${res}m.tif)
#list=$(ls */stereo/*00/dem*/*-DEM_${res}m.tif */*00/dem*/*-DEM_${res}m.tif)
#list=$(ls */*/*/dem*/*-DEM_8m_trans.tif | grep -v QB)

#ext="-DEM_${res}m"
#ext="-DEM_${res}m_dem_align"
#ext="-DEM_${res}m_dzfilt_-200_200"
ext="-DEM_${res}m_dzfilt_-200_200_*align"

mos_ext="${site}_mos_${res}m"
#mos_ext="${site}_mos_${res}m_dem_align"

echo $re
echo $ext
echo $mos_ext

#list=$(ls *track/*00/dem*/${re}*${ext}.tif)
list=$(ls *align/${re}*${ext}.tif)
#list=$(ls validpairs*/*00/dem*/${re}*${ext}.tif)
echo $list | wc -w

if $noQB ; then
    echo "Removing QB02"
    #list=$(echo $list | tr ' ' '\n' | grep -v QB02)
    list=$(echo $list | tr ' ' '\n' | grep -v '_101[0-9A-Z]*_101')
    echo $list | wc -w
fi

#HMA Nuth and Kaab, round 2
#list=$(ls *track/dem_coreg/*00/dem*/*-DEM_8m_dem_align/*_align_dem_align/*_align.tif)

out=mos/${site}_${ts}_mos/${mos_ext}/${mos_ext}

#Sort by date
#NOTE: Need to update sort key with increased path depths
list=$(echo $list | tr ' ' '\n' | sort -n -t'/' -k 2)
#list=$(echo $list | tr ' ' '\n' | sort -n -t'/' -k 3)
#list=$(echo $list | tr ' ' '\n' | sort -n -t'/' -k 4)
#list=$(echo $list | tr ' ' '\n' | sort -n -t'/' -k 5)

echo
echo $out
echo
echo $list
echo

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
    if $out_tif ; then
        gdalwarp $gdal_opt $out.vrt ${out}.tif
        gdaladdo_ro.sh ${out}.tif
        hs.sh ${out}.tif
        gdaladdo_ro.sh ${out}_hs_az315.tif
    fi
fi

if $last; then
    if [ ! -e ${out}_last.vrt ] ; then
        echo "Preparing lastmap"
        $mos --stat last --threads $threads --tr $res --t_srs "$proj" --georef_tile_size=$tilesize -o $out $list
        if $out_tif ; then
            gdalwarp $gdal_opt $out.vrt ${out}_last.tif
            gdaladdo_ro.sh ${out}_last.tif
            hs.sh ${out}_last.tif
            gdaladdo_ro.sh ${out}_last_hs_az315.tif
        fi
    fi
fi

if $first; then
    if [ ! -e ${out}_first.vrt ] ; then
        echo "Preparing firstmap"
        $mos --stat first --threads $threads --tr $res --t_srs "$proj" --georef_tile_size=$tilesize -o $out $list
        if $out_tif ; then
            gdalwarp $gdal_opt $out.vrt ${out}_first.tif
            gdaladdo_ro.sh ${out}_first.tif
            hs.sh ${out}_first.tif
            gdaladdo_ro.sh ${out}_first_hs_az315.tif
        fi
    fi
fi

if (( "$res" == "32" )) ; then
    #If res above is 32, make lowres products for faster browsing
    lowres=100
    #lowres=32

    #Should run these in parallel
    gdal_opt='-co TILED=YES -co COMPRESS=LZW -co BIGTIFF=YES'
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
    if $stddev; then
        if [ ! -e ${out}_stddev.vrt ] ; then
            echo "Preparing stddev map"
            $mos --stat stddev --threads $threads --tr $res --t_srs "$proj" --georef_tile_size=$tilesize -o $out $list
        fi
        if [ ! -e ${out}_stddev_${lowres}m.tif ] ; then
            echo "Preparing lowres $lowres m stddev map"
            gdalwarp -tr $lowres $lowres $gdal_opt ${out}_stddev.vrt ${out}_stddev_${lowres}m.tif
            gdaladdo_ro.sh ${out}_stddev_${lowres}m.tif
        fi
    fi
    if $last ; then
        if [ ! -e ${out}_last.vrt ] ; then
            echo "Preparing last timestamp mosaic"
            $mos --stat last --threads $threads --tr $res --t_srs "$proj" --georef_tile_size=$tilesize -o $out $list
        fi
        if [ ! -e ${out}_last_${lowres}m.tif ] ; then
            echo "Preparing lowres $lowres m lastmap"
            gdalwarp -tr $lowres $lowres $gdal_opt ${out}_last.vrt ${out}_last_${lowres}m.tif
            gdaladdo_ro.sh ${out}_last_${lowres}m.tif
        fi
    fi
    if $first ; then
        if [ ! -e ${out}_first.vrt ] ; then
            echo "Preparing first timestamp mosaic"
            $mos --stat first --threads $threads --tr $res --t_srs "$proj" --georef_tile_size=$tilesize -o $out $list
        fi
        if [ ! -e ${out}_first_${lowres}m.tif ] ; then
            echo "Preparing lowres $lowres m firstmap"
            gdalwarp -tr $lowres $lowres $gdal_opt ${out}_first.vrt ${out}_first_${lowres}m.tif
            gdaladdo_ro.sh ${out}_first_${lowres}m.tif
        fi
    fi
    if $index ; then
        if [ ! -e ${out}_stripindex.shp ] ; then
            echo "Generating strip index shp"
            #This exports bounding boxes, but doesn't show actual footprint of valid pixels
            #gdaltindex -t_srs "$proj" ${out}_stripindex.shp $list
            parallel -j $threads "if [ ! -e {.}_${lowres}m.tif ] ; then gdalwarp -overwrite -r cubic -t_srs \"$proj\" -tr $lowres $lowres -dstnodata -9999 {} {.}_${lowres}m.tif; fi" ::: $list
            raster2shp.py -merge_fn ${out}_stripindex.shp $(echo $list | sed "s/DEM_${res}m.tif/DEM_${res}m_${lowres}m.tif/g") 
            #This needs more careful testing, currently deletes source files!
            #echo "Removing intermediate files"
            #eval rm $(echo $list | sed "s/DEM_${res}m.tif/DEM_${res}m_${lowres}m.{shp,shx,dbf,prj,tif}/g")
            ogr2ogr -simplify $tol ${out}_stripindex_simp.shp ${out}_stripindex.shp
            rm ${out}_stripindex.{shp,dbf,prj,shx}
            if [ ! -e ${out}_stripindex_simp.kml ] ; then 
                ogr2ogr -f KML ${out}_stripindex_simp.kml ${out}_stripindex_simp.shp
            fi
        fi
    fi
    #This is a hack for now, should add this functionality to dem_mosaic_validtiles.py
    if $tileindex ; then 
        echo "Generating tile index shp"
        #gdaltindex ${out}_tileindex.shp $(gdalinfo $out.vrt | grep tif)
        raster2shp.py -merge_fn ${out}_tileindex.shp $(gdalinfo $out.vrt | grep tif)
        ogr2ogr -simplify $tol ${out}_tileindex_simp.shp ${out}_tileindex.shp
        echo "Removing intermediate files"
        rm ${out}_tileindex.{shp,dbf,prj,shx}
        eval rm $(gdalinfo $out.vrt | grep tif | sed 's/tif/{shp,dbf,prj,shx}/g')
        if [ ! -e ${out}_tileindex_simp.kml ] ; then 
            ogr2ogr -f KML ${out}_tileindex_simp.kml ${out}_tileindex_simp.shp
        fi
    fi
fi

#Set permissions
chmod -R 755 mos
chmod 644 ${out}*

#Create symlink pointing to latest mosaic dir
latest=mos/latest
if [ -e $latest ] ; then
    rm $latest
fi
ln -s ${site}_${ts}_mos $latest
