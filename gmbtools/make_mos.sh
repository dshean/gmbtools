#! /bin/bash 

set -e

#If interrupted
#for i in *lastindex.tif; do
#if [ ! -e $i-log-dem_mosaic*txt ] ; then rm ${i}*; fi ; done

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

#If res is better than 32, make lowres products for faster browsing
lowres=100
#lowres=32

#Write out tif mosaics at full res
fullres_tif=false

count=false
stddev=false
median=false
last=true
lastindex=true
first=true
firstindex=true
stripindex=false
tileindex=false

#Exclude Quickbird-2
noQB=true

#Default computes union from input DEMs
extent='union'

#Hardcoded site and projection for now
site=hma
proj='+proj=aea +lat_1=25 +lat_2=47 +lat_0=36 +lon_0=85 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs '
extent='-1553632.99074 -1030104.4196 1727255.00926 1268847.5804'
#site=conus
#proj='+proj=aea +lat_1=36 +lat_2=49 +lat_0=43 +lon_0=-115 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs '
#site=fuego
#proj='EPSG:32615'

#Mosaic tile size in meters
#tilesize=20000
tilesize=100000
mos=~/src/gmbtools/gmbtools/dem_mosaic_validtiles.py
#Simplify tolerance in decimal degrees
tol=0.001

ncpu=$(cat /proc/cpuinfo | egrep "core id|physical id" | tr -d "\n" | sed s/physical/\\nphysical/g | grep -v ^$ | sort | uniq | wc -l)
threads=$((ncpu-1))

ts=`date +%Y%m%d`
ts=20180724

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
#WV/GE dem_align
ext="-DEM_${res}m_dzfilt_-200_200_*align"

#ASTER
#ext="_align_dzfilt_-100_100"

mos_ext="${site}_mos_${res}m"
#mos_ext=${site}_${ts}_mos

echo $re
echo $ext
echo $mos_ext

#list=$(ls *00/dem*/${re}*${ext}.tif)
#list=$(ls *track/*00/dem*/${re}*${ext}.tif)
#WV/GE dem_align
list=$(ls *align/${re}*${ext}.tif)

#ASTER
#list=$(ls 2*/*cr_dem_align/*${re}*${ext}.tif)
echo $list | wc -w

if $noQB ; then
    echo "Removing QB02"
    #list=$(echo $list | tr ' ' '\n' | grep -v QB02)
    list=$(echo $list | tr ' ' '\n' | grep -v '_101[0-9A-Z]*_101')
    echo $list | wc -w
fi

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
j=10
echo "First $j entries":
echo $list | tr ' ' '\n' | head -$j
echo
echo "Last $j entries":
echo $list | tr ' ' '\n' | tail -n $j
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
gdal_opt='-co COMPRESS=LZW -co TILED=YES -co BIGTIFF=IF_SAFER'

#~/src/gmbtools/gmbtools/dem_mosaic_validtiles.py --threads 27 --tr 32 --t_srs '+proj=aea +lat_1=25 +lat_2=47 +lat_0=36 +lon_0=85 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ' --georef_tile_size 100000 -o mos/hma_20180731_mos 2*/*align/*100.tif

#Define common dem_mosaic_validtiles arguments
mos_arg="--threads $threads --tr $res --t_srs \"$proj\" --t_projwin \"$extent\" --georef_tile_size=$tilesize -o $out $list"

if [ ! -e $out.vrt ] ; then
    eval $mos $mos_arg
    if $fullres_tif ; then
        gdalwarp $gdal_opt $out.vrt ${out}.tif
        gdaladdo_ro.sh ${out}.tif
        #Should hs each tile, then vrt, then warp, as above
        #hs.sh ${out}.tif
        #gdaladdo_ro.sh ${out}_hs_az315.tif
    fi
fi

if $count ; then
    if [ ! -e ${out}_count.vrt ] ; then
        echo "Preparing countmap"
        eval $mos --stat count $mos_arg
        if $fullres_tif ; then
            gdalwarp $gdal_opt ${out}_count.vrt ${out}_count.tif
            gdaladdo_ro.sh ${out}_count.tif
        fi
    fi
fi

if $stddev ; then
    if [ ! -e ${out}_stddev.vrt ] ; then
        echo "Preparing stddevmap"
        eval $mos --stat stddev $mos_arg
        if $fullres_tif ; then
            gdalwarp $gdal_opt ${out}_stddev.vrt ${out}_stddev.tif
            gdaladdo_ro.sh ${out}_stddev.tif
        fi
    fi
fi

if $median ; then
    if [ ! -e ${out}_median.vrt ] ; then
        echo "Preparing medianmap"
        eval $mos --stat median $mos_arg
        if $fullres_tif ; then
            gdalwarp $gdal_opt ${out}_median.vrt ${out}_median.tif
            gdaladdo_ro.sh ${out}_median.tif
        fi
    fi
fi

if $last; then
    if [ ! -e ${out}_last.vrt ] ; then
        echo "Preparing lastmap"
        eval $mos --stat last $mos_arg
        if $fullres_tif ; then
            gdalwarp $gdal_opt ${out}_last.vrt ${out}_last.tif
            gdaladdo_ro.sh ${out}_last.tif
        fi
    fi
fi

if $lastindex; then
    if [ ! -e ${out}_lastindex_ts.vrt ] ; then
        echo "Preparing lastindex map"
        eval $mos --stat lastindex $mos_arg
        if $fullres_tif ; then
            gdalwarp $gdal_opt ${out}_lastindex_ts.vrt ${out}_lastindex_ts.tif
            gdaladdo_ro.sh ${out}_lastindex_ts.tif
        fi
    fi
fi

if $first; then
    if [ ! -e ${out}_first.vrt ] ; then
        echo "Preparing firstmap"
        eval $mos --stat first $mos_arg
        if $fullres_tif ; then
            gdalwarp $gdal_opt $out.vrt ${out}_first.tif
            gdaladdo_ro.sh ${out}_first.tif
        fi
    fi
fi

if $firstindex; then
    if [ ! -e ${out}_firstindex_ts.vrt ] ; then
        echo "Preparing firstindex map"
        eval $mos --stat firstindex $mos_arg
        if $fullres_tif ; then
            gdalwarp $gdal_opt ${out}_firstindex_ts.vrt ${out}_firstindex_ts.tif
            gdaladdo_ro.sh ${out}_firstindex_ts.tif
        fi
    fi
fi

#Generat lowres products
if (( "$res" == "32" )) ; then
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
        if [ ! -e ${out}_count_${lowres}m.tif ] ; then
            echo "Preparing lowres $lowres m countmap"
            gdalwarp -tr $lowres $lowres $gdal_opt ${out}_count.vrt ${out}_count_${lowres}m.tif
            gdaladdo_ro.sh ${out}_count_${lowres}m.tif
        fi
    fi
    if $stddev; then
        if [ ! -e ${out}_stddev_${lowres}m.tif ] ; then
            echo "Preparing lowres $lowres m stddev map"
            gdalwarp -tr $lowres $lowres $gdal_opt ${out}_stddev.vrt ${out}_stddev_${lowres}m.tif
            gdaladdo_ro.sh ${out}_stddev_${lowres}m.tif
        fi
    fi
    if $median ; then
        if [ ! -e ${out}_median_${lowres}m.tif ] ; then
            echo "Preparing lowres $lowres m medianmap"
            gdalwarp -tr $lowres $lowres $gdal_opt ${out}_median.vrt ${out}_median_${lowres}m.tif
            gdaladdo_ro.sh ${out}_median_${lowres}m.tif
        fi
    fi
    if $last ; then
        if [ ! -e ${out}_last_${lowres}m.tif ] ; then
            echo "Preparing lowres $lowres m lastmap"
            gdalwarp -tr $lowres $lowres $gdal_opt ${out}_last.vrt ${out}_last_${lowres}m.tif
            gdaladdo_ro.sh ${out}_last_${lowres}m.tif
        fi
    fi
    if $lastindex ; then
        if [ ! -e ${out}_lastindex_ts_${lowres}m.tif ] ; then
            echo "Preparing lowres $lowres m lastindexmap"
            gdalwarp -tr $lowres $lowres $gdal_opt ${out}_lastindex_ts.vrt ${out}_lastindex_ts_${lowres}m.tif
            gdaladdo_ro.sh ${out}_lastindex_ts_${lowres}m.tif
        fi
    fi
    if $first ; then
        if [ ! -e ${out}_first_${lowres}m.tif ] ; then
            echo "Preparing lowres $lowres m firstmap"
            gdalwarp -tr $lowres $lowres $gdal_opt ${out}_first.vrt ${out}_first_${lowres}m.tif
            gdaladdo_ro.sh ${out}_first_${lowres}m.tif
        fi
    fi
    if $firstindex ; then
        if [ ! -e ${out}_firstindex_ts_${lowres}m.tif ] ; then
            echo "Preparing lowres $lowres m firstindexmap"
            gdalwarp -tr $lowres $lowres $gdal_opt ${out}_firstindex_ts.vrt ${out}_firstindex_ts_${lowres}m.tif
            gdaladdo_ro.sh ${out}_firstindex_ts_${lowres}m.tif
        fi
    fi

    if $stripindex ; then
        if [ ! -e ${out}_stripindex.shp ] ; then
            echo "Generating strip index shp"
            #This exports bounding boxes, but doesn't show actual footprint of valid pixels
            #gdaltindex -t_srs "$proj" ${out}_stripindex.shp $list
            parallel -j $threads "if [ ! -e {.}_${lowres}m.tif ] ; then gdalwarp -overwrite -r cubic -t_srs \"$proj\" -tr $lowres $lowres -dstnodata -9999 {} {.}_${lowres}m.tif; fi" ::: $list
            raster2shp.py -merge_fn ${out}_stripindex.shp $(echo $list | sed "s/.tif/_${lowres}m.tif/g") 
            #This needs more careful testing, currently deletes source files!
            #echo "Removing intermediate files"
            #eval rm $(echo $list | sed "s/DEM_${res}m.tif/DEM_${res}m_${lowres}m.{shp,shx,dbf,prj,tif}/g")
            ogr2ogr -simplify $tol ${out}_stripindex_simp.shp ${out}_stripindex.shp
            rm ${out}_stripindex.{shp,dbf,prj,shx}
            if [ ! -e ${out}_stripindex_simp.kml ] ; then 
                ogr2ogr -f KML ${out}_stripindex_simp.kml ${out}_stripindex_simp.shp
            fi
            #Generate bar plot of annual count
            plot_annual_count.py ${out}_stripindex_simp.shp
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
