#! /bin/bash 

set -e

#If interrupted
#Only if tar.gz have not been generated
#for i in *-tile-[0-9]*-*.tif; do if ! ls ${i}-log-dem_mosaic*txt 1> /dev/null 2>&1 ; then rm -v ${i}*; else if ! grep -q 'Number of valid' $(ls -t ${i}-log-dem_mosaic*txt | head -1); then rm -v ${i}* ; fi ; fi ; done

#pbs_rfe --duration 0+6 --model bro
#qsub -I -q devel -lselect=1:model=bro,walltime=2:00:00
#qsub -I -q long -lselect=1:model=bro,walltime=6:00:00
#cd /nobackup/deshean/hma
#make_mos.sh

#Want to make cascading mosaic
#gdalbuildvrt mos_8m_all mos_8m_all_trans mos_8m_summer_trans
#Note on ordering within vrt
#If there is some amount of spatial overlapping between files, the order of files appearing in the list of source matter: files that are listed at the end are the ones from which the content will be fetched. Note that nodata will be taken into account to potentially fetch data from less prioritary datasets

gdal_opt='-co TILED=YES -co COMPRESS=LZW -co BIGTIFF=YES'

#Set open file limit
#Default is 2048
ulimit -n 65536

latest=false
summer=false
trans=false

#Output mosaic res
#res=2
#res=8
res=10
#res=32
#ASTER
#res=30

#If res is better than 32, make lowres products for faster browsing
#lowres=100
lowres=90
#lowres=32

#Write out tif mosaics at full res
fullres_tif=false
lowres_tif=true

#Specify the output types
statlist=""
statlist+=" wmean"
statlist+=" count"
statlist+=" stddev"
statlist+=" median medianindex"
statlist+=" nmad"
#statlist+=" last lastindex first firstindex"

#Generate strip index shp
stripindex=false

#Generate tile index shp
tileindex=false

#Exclude Quickbird-2
noQB=true

#Exclude GeoEye-1
noGE=false
GEonly=false

#Default computes union from input DEMs
extent='union'

#Hardcoded sitename, projection and extent (performance improvement)
site=hma
#proj='+proj=aea +lat_1=25 +lat_2=47 +lat_0=36 +lon_0=85 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs '

#Export 1x1 deg tiles
proj='EPSG:4326'
extent='64 24 109 45'

#1 arcsec = 1/3600 of deg = 30 m
#res=0.0002777777777778
#1/3 arcsec = 1/3600/3 = 10 m
res=0.0000925925925926
#1/9 arcsec = 1/3600/9 = 3 m
#res=0.000030864197531

lowres=$(python -c "print($res*3.*3.)")

#WV/GE extent
#extent='-1553632.99074 -1030104.4196 1727255.00926 1268847.5804'
#ASTER extent
#extent='-1604981.73315 -1094260.0 1847978.26685 1161996.0'
#Should take union of these, so we have uniform tile boundaries
#site=conus
#proj='+proj=aea +lat_1=36 +lat_2=49 +lat_0=43 +lon_0=-115 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs '
#extent='-684633.544035 -694442.824668 824774.455965 804533.175332'
#site=fuego
#proj='EPSG:32615'

#Mosaic tile size in meters
#tilesize=10000
#tilesize=100000
tilesize=1
mos=~/src/gmbtools/gmbtools/dem_mosaic_validtiles.py
#Simplify tolerance in decimal degrees
tol=0.001

#If res is 8, use all physical cores (avoid memory caching slowdows)
ncpu=$(cat /proc/cpuinfo | egrep "core id|physical id" | tr -d "\n" | sed s/physical/\\nphysical/g | grep -v ^$ | sort | uniq | wc -l)
#If res is 32, use all virtual cores
#ncpu=$(python -c 'import multiprocessing as mp; print(mp.cpu_count())')
threads=$((ncpu-1))
#threads=16

ts=`date +%Y%m%d`
#ts=20181003
#ts=20180830
#ts=20181127_GEonly
#ts=20181130_aster
#ts=20181205

#Should add option to split annually
echo "Identifying input DEMs"

#mos_ext="${site}_mos_${res}m"
#mos_ext="${site}_mos_1arcsec"
mos_ext="${site}_wvge_mos_13arcsec"

if $latest ; then
    #Note: CONUS is only current through June 2016, missing 
    re='201[5-9]'
    mos_ext+='_latest'
else
    #Exclude 2008/2009 data
    #re='201[0-9]'
    #Include all
    #re='20[0-9]'
    re=''
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

#Standard WV
#ext="-DEM_${res}m"
#list=$(ls *00/dem*/${re}*${ext}.tif)

#HMA alongtrack/crosstrack
#list=$(ls *track/*00/dem*/${re}*${ext}.tif)

#WV/GE dem_align
ext="-DEM_8m_dzfilt_-200_200_*align"
#ext="-DEM_${res}m_dzfilt_-200_200_*align"
#Reference products (masked to remove outliers)
#ext="-DEM_${res}m_dzfilt_-200_200_*align_filt.tif"
list=$(ls *track/2*align/${re}*${ext}.tif)

#ASTER
#ext="_align_dzfilt_-100_100"
#ext='_align'
#list=$(ls 2*/AST*_dem_align/*${re}*${ext}.tif)
#Round2
#list=$(ls 2*/AST*_dem_align/AST*_dem_align/*${re}*${ext}.tif)

echo $list | wc -w

if $noQB ; then
    echo "Removing QB02"
    #list=$(echo $list | tr ' ' '\n' | grep -v QB02)
    list=$(echo $list | tr ' ' '\n' | grep -v '_101[0-9A-Z]*_101')
    echo $list | wc -w
fi

#GE1 can have nasty artifacts, should isolate and then merge GE01 composite with WV composite
if $noGE ; then
    echo "Removing GE01"
    list=$(echo $list | tr ' ' '\n' | grep -v '_105[0-9A-Z]*_105')
    echo $list | wc -w
fi

if $GEonly ; then
    echo "Isolating GE01"
    list=$(echo $list | tr ' ' '\n' | grep '_105[0-9A-Z]*_105')
    echo $list | wc -w
fi

out=mos/${site}_${ts}_mos/${mos_ext}/${mos_ext}

#Sort by date
#NOTE: Need to update sort key with increased path depths
#*align/*DEM.tif 
list=$(echo $list | tr ' ' '\n' | sort -n -t'/' -k 2)
#*00/dem*/*DEM.tif
#list=$(echo $list | tr ' ' '\n' | sort -n -t'/' -k 3)
#*track/*00/dem*/*DEM.tif
#*00/dem*/*align/*DEM.tif
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
    stripecount=8
    lfs setstripe -c $stripecount $(dirname $out)
fi

echo $list | tr ' ' '\n' > ${out}_input_DEM_list.txt
echo $(echo $list | wc -w) input DEMs

#~/src/gmbtools/gmbtools/dem_mosaic_validtiles.py --threads 27 --tr 32 --t_srs '+proj=aea +lat_1=25 +lat_2=47 +lat_0=36 +lon_0=85 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ' --georef_tile_size 100000 -o mos/hma_20180731_mos 2*/*align/*100.tif

#Check for existing final products
echo $statlist
statlist_todo=''
for stat in $statlist
do
    if [ ! -e ${out}_${stat}.vrt ] && [ ! -e ${out}_${stat}_ts.vrt ] ; then
        statlist_todo+=" $stat"
    else
        echo "Found existing vrt: ${out}_${stat}.vrt"
    fi
done
echo $statlist_todo

#If some stat mosaics have not yet been generated
if [[ ! -z "${statlist_todo// }" ]] ; then
    #Run dem_mosaic_validtiles for all stats
    #Define common dem_mosaic_validtiles arguments
    #mos_arg="--threads $threads --tr $res --t_srs \"$proj\" --t_projwin \"$extent\" --georef_tile_size=$tilesize -o $out $list"
    mos_arg="--threads $threads --tr $res --t_srs \"$proj\" --t_projwin \"$extent\" --georef_tile_size=$tilesize -o $out ${out}_input_DEM_list.txt"
    echo
    echo eval $mos --stat $statlist_todo $mos_arg
    echo
    eval $mos --stat $statlist_todo $mos_arg
fi

#Mask output median or wmean to preserve pixels with count > 2, nmad < 3.0
#Should probably move this to dem_mosaic_validtiles
statlist_tomask="median wmean"
for stat in $statlist_tomask
do
    fn_list=$(ls ${out}*-${stat}.tif)
    fn_list_todo=""
    for fn in $fn_list
    do
        if [ ! -e ${fn%.*}_masked.tif ] ; then 
           fn_list_todo+=" $fn"
        fi
    done 
    if [ ! -z "$fn_list_todo" ] ; then 
        echo "Running dem_mosaic_mask for $stat"
        #These can be large, quickly fill up memory
        parallel --progress -j 16 --delay 0.1 'dem_mosaic_mask.py {}' ::: $fn_list_todo
    fi
    if [ ! -e ${out}_${stat}_masked.vrt ] ; then 
        echo "Building vrt for $stat"
        gdalbuildvrt -r cubic ${out}_${stat}_masked.vrt ${out}*-${stat}_masked.tif
    fi
    statlist+=" ${stat}_masked"
done

#Define function to convert vrt to output tif products
stat_vrt2tif() {
    in=${out}_${1}
    if $fullres_tif ; then
        if [ ! -e $in.tif ] ; then
            echo "Converting vrt to tif: $in.vrt"
            gdalwarp $gdal_opt $in.vrt $in.tif
            echo "Adding overviews: $in.tif"
            gdaladdo_ro.sh $in.tif
            #hs.sh $in.tif
            #gdaladdo_ro.sh ${in}_hs_az315.tif
        fi
    fi
    if $lowres_tif ; then
        if [ ! -e ${in}_${lowres_str}.tif ] ; then
            echo "Converting vrt to lowres tif: $in.vrt"
            gdalwarp -tr $lowres $lowres $gdal_opt $in.vrt ${in}_${lowres_str}.tif
            echo "Adding overviews: ${in}_${lowres_str}.tif"
            gdaladdo_ro.sh ${in}_${lowres_str}.tif
            #hs.sh ${in}_${lowres_str}.tif
            #gdaladdo_ro.sh ${in}_${lowres_str}_hs_az315.tif
        fi
    fi
}

#lowres_str=${lowres}m
lowres_str="sub9"

echo "Running vrt2tif"
echo $statlist
export out
export gdal_opt
export fullres_tif
export lowres_tif
export lowres
export lowres_str
export -f stat_vrt2tif
#Run tif generation and lowres in parallel for all stat mosaics
statlist=$(echo $statlist | sed -e 's/lastindex/lastindex_ts/' -e 's/firstindex/firstindex_ts/' -e 's/medianindex/medianindex_ts/')
parallel --progress --env _ stat_vrt2tif ::: $statlist 

#Generate lowres products
if (( "$res" == "32" )) ; then
    if $stripindex ; then
        if [ ! -e ${out}_stripindex.shp ] ; then
            echo "Generating strip index shp"
            #This exports bounding boxes, but doesn't show actual footprint of valid pixels
            #gdaltindex -t_srs "$proj" ${out}_stripindex.shp $list
            parallel -j $threads "if [ ! -e {.}_${lowres_str}.tif ] ; then gdalwarp -overwrite -r cubic -t_srs \"$proj\" -tr $lowres $lowres -dstnodata -9999 {} {.}_${lowres_str}.tif; fi" ::: $list
            raster2shp.py -merge_fn ${out}_stripindex.shp $(echo $list | sed "s/.tif/_${lowres_str}.tif/g") 
            #This needs more careful testing, currently deletes source files!
            #echo "Removing intermediate files"
            #eval rm $(echo $list | sed "s/DEM_${res}m.tif/DEM_${res}m_${lowres_str}.{shp,shx,dbf,prj,tif}/g")
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
