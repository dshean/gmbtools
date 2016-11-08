#! /bin/bash

#Post-process CONUS glacier DEMs

topdir=/nobackup/deshean/conus
cd $topdir

threads=32
ext=DEM_32m
#ext=DEM_32m
list=$(ls ../conus[0-9]/{WV,GE}*/*/*${ext}.tif)
#list=$(ls ../conus[0-9]/{WV,GE}*/*/*${ext}.tif | grep 201[0-9][01][7890])

#outres=100
outres=100
lowres=1000
epsg=32611

neddir=/nobackup/deshean/rpcdem/ned1
ned=$neddir/ned1_tiles_glac24k_115kmbuff.vrt
nedres=30

if true ; then
    #Post-processing individual DEMs with original UTM proj
    parallel "gdalwarp -overwrite -r cubic -t_srs EPSG:$epsg -tr $outres $outres -dstnodata -9999 {} {.}_${epsg}_${outres}m.tif" ::: $list
    list_proj=$(echo $list | sed "s/$ext.tif/${ext}_${epsg}_${outres}m.tif/g")
fi

if true ; then 
    gdalbuildvrt -vrtnodata -9999 conus_${epsg}_${outres}m.vrt $list_proj
    if [ ! -e conus_${epsg}_${outres}m-tile-0-count.tif ] ; then 
        dem_mosaic --threads $threads --count $list_proj -o conus_${epsg}_${outres}m
    fi
    if [ ! -e conus_${epsg}_${lowres}m-tile-0-count.tif ] ; then 
        #dem_mosaic --threads $threads --tr $lowres --count $list_proj -o conus_${epsg}_${lowres}m
        gdalwarp -r cubic -overwrite -dstnodata 0 -ot UInt16 -tr $lowres $lowres conus_${epsg}_${outres}m-tile-0-count.tif conus_${epsg}_${outres}m-tile-0-count_1km.tif
    fi
    if [ ! -e conus_${epsg}_${outres}m-tile-0-stddev.tif ] ; then 
        dem_mosaic --threads $threads --stddev $list_proj -o conus_${epsg}_${outres}m
    fi
    if [ ! -e conus_${epsg}_${lowres}m-tile-0-stddev.tif ] ; then 
        #dem_mosaic --threads $threads --tr $lowres --stddev $list_proj -o conus_${epsg}_${lowres}m
        gdalwarp -r cubic -overwrite -dstnodata 0 -tr $lowres $lowres conus_${epsg}_${outres}m-tile-0-stddev.tif conus_${epsg}_${outres}m-tile-0-stddev_1km.tif
    fi
    if [ ! -e conus_${epsg}_${outres}m-tile-0.tif ] ; then 
        dem_mosaic --threads $threads $list_proj -o conus_${epsg}_${outres}m
    fi
    if [ ! -e conus_${epsg}_${lowres}m-tile-0.tif ] ; then 
        #dem_mosaic --threads $threads --tr $lowres $list_proj -o conus_${epsg}_${lowres}m
        gdalwarp -r cubic -overwrite -dstnodata -9999 -tr $lowres $lowres conus_${epsg}_${outres}m-tile-0.tif conus_${epsg}_${outres}m-tile-0_1km.tif
    fi
    hs.sh conus_${epsg}_${outres}m-tile-0.tif
    hs.sh conus_${epsg}_${lowres}m-tile-0.tif

    if false ; then 
        compute_dh.py $ned conus_${epsg}_${outres}m-tile-0.tif
        mv ${ned%.*}_conus_${epsg}_${outres}m-tile-0_dz_eul.tif .
    fi
fi

ts=$(date +%Y%m%d_%H%M)
n=$(echo $list | wc -w)
outshp=shp/conus_${epsg}_${outres}m_n${n}_${ts}.shp
#This creates issues b/c individual DEMs span multiple UTM zones
#raster2shp.py $list
#Reproject to same srs
raster2shp.py $list_proj
shp_rename.sh merge.shp $outshp 

exit

#ogr2ogr -t_srs EPSG:proj ${outshp%.*}_proj.shp $outshp
#shpfltr_by_pt.sh ${outshp%.*}_proj.shp
#shpfltr_by_pt.sh $outshp

#Create 8-m tiled mosaic
#ext=DEM_8m_trans
#res=8
res=32
ext=DEM_${res}m
list=$(ls ../conus[0-9]/{WV,GE}*/*/*${ext}.tif)
dem_mosaic --threads $threads --tile-size 10000 --t_srs EPSG:32611 --tr $res -o mos_${res}m/mos_${res}m_all/mos_${res}m_all $list
gdalbuildvrt -vrtnodata -9999 -r cubic mos_${res}m/mos_${res}m_all/mos_${res}m_all.vrt mos_${res}m/mos_${res}m_all/mos_${res}m_all-tile-0.tif 
extent=$(get_extent mos_${res}m/mos_${res}m_all/mos_${res}m_all.vrt) 
#This won't work due to heterogenous projection
#gdalbuildvrt -vrtnodata 0 -r cubic all_${res}m.vrt $list
#extent=$(get_extent.py all_${res}m.vrt)

list_summer=$(ls ../conus[0-9]/{WV,GE}*/*/*${ext}.tif | grep 201[0-9][01][67890])
list_winter=$(ls ../conus[0-9]/{WV,GE}*/*/*${ext}.tif | grep -v 201[0-9][01][67890])

#These should all use the same extent, identical tiles
dem_mosaic --threads $threads --tile-size 10000 --t_projwin $extent --t_srs EPSG:32611 --tr 8 -o mos_${res}m/mos_${res}m_summer/mos_${res}m_summer $list_summer
dem_mosaic --threads $threads --tile-size 10000 --t_projwin $extent --t_srs EPSG:32611 --tr 8 -o mos_${res}m/mos_${res}m_winter/mos_${res}m_winter $list_winter

#Many tiles are empty, remove
mkdir valid
mv $(ll -Sr | grep -v 1.4M | grep -v K | awk '{print $NF}') valid/
cd valid
gdalbuildvrt -vrtnodata 0 -r cubic mos_${res}m_valid.vrt *tif
