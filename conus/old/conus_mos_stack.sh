#! /bin/bash

#Post-process hma glacier DEMs
#topdir=/nobackupp8/deshean/hma/hma1_2016dec22/stereo
#mos_prefix=hma1_stereo_trans
#topdir=/nobackupp8/deshean/hma/hma1_2016dec22/validpairs
topdir=/nobackupp8/deshean/hma/hma1_2016dec22
mos_prefix=hma1_all
lowres_dem=/nobackup/deshean/rpcdem/hma/srtm1/hma_srtm_gl1.vrt
epsg=32644

#CONUS
topdir=/nobackupp8/deshean/conus/dem2
mos_prefix=conus
lowres_dem=/nobackup/deshean/rpcdem/ned13/ned13_tiles_glac24k_115kmbuff.vrt
epsg=32611

cd $topdir

threads=32

#Define extension type
#ext=DEM_32m
#ext=DEM_32m_trans
ext=DEM_8m_trans

#Need this for hma stereo and validpairs subdir
#list=$(ls */{QB,IK,WV,GE}*/*/*${ext}.tif)
list=$(ls {QB,IK,WV,GE}*/*/*${ext}.tif)
#Summer
#list=$(ls {QB,IK,WV,GE}*/*/*${ext}.tif | grep 201[0-9][01][7890])

outres=8
#outres=100
make_lowres=false
lowres=1000

if true ; then
    #Post-processing individual DEMs with original UTM proj
    parallel -j $threads "if [ ! -e {.}_${epsg}_${outres}m.tif ] ; then gdalwarp -overwrite -r cubic -t_srs EPSG:$epsg -tr $outres $outres -dstnodata -9999 {} {.}_${epsg}_${outres}m.tif; fi" ::: $list
    list_proj=$(echo $list | sed "s/$ext.tif/${ext}_${epsg}_${outres}m.tif/g")
fi

#This creates 8-m mosaic
if [ "$outres" -eq "8" ] ; then 
    mkdir ${mos_prefix}_${epsg}_${outres}m
    dem_mosaic --threads $threads $list_proj -o ${mos_prefix}_${epsg}_${outres}m/${mos_prefix}_${epsg}_${outres}m --georef-tile-size 100000
    cd ${mos_prefix}_${epsg}_${outres}m
    gdalbuildvrt ${mos_prefix}_${epsg}_${outres}m_mos.vrt ${mos_prefix}_${epsg}_${outres}m-tile-*.tif
    exit
fi

if true ; then 
    if [ ! -e ${mos_prefix}_${epsg}_${outres}m-tile-0.tif ] ; then 
        echo "dem_mosaic"
        dem_mosaic --threads $threads $list_proj -o ${mos_prefix}_${epsg}_${outres}m
        echo "building overviews"
        gdaladdo_ro.sh ${mos_prefix}_${epsg}_${outres}m-tile-0.tif
        echo "shaded relief"
        hs.sh ${mos_prefix}_${epsg}_${outres}m-tile-0.tif
        echo "shaded relief overviews"
        gdaladdo_ro.sh ${mos_prefix}_${epsg}_${outres}m-tile-0_hs_az315.tif
    fi
    if $make_lowres ; then 
        if [ ! -e ${mos_prefix}_${epsg}_${outres}m-tile-0_${lowres}m.tif ] ; then 
            #dem_mosaic --threads $threads --tr $lowres $list_proj -o ${mos_prefix}_${epsg}_${lowres}m
            gdalwarp -r cubic -overwrite -dstnodata -9999 -tr $lowres $lowres ${mos_prefix}_${epsg}_${outres}m-tile-0.tif ${mos_prefix}_${epsg}_${outres}m-tile-0_${lowres}m.tif
            hs.sh ${mos_prefix}_${epsg}_${outres}m-tile-0_${lowres}m.tif 
        fi
    fi
fi

if false ; then 
    if [ ! -d shp ] ; then
        mkdir shp
    fi
    ts=$(date +%Y%m%d_%H%M)
    n=$(echo $list | wc -w)
    outshp=shp/${mos_prefix}_${epsg}_${outres}m_n${n}_${ts}.shp
    #This creates issues b/c individual DEMs span multiple UTM zones
    #raster2shp.py $list
    #Reproject to same srs
    raster2shp.py $list_proj
    shp_rename.sh merge.shp $outshp 
fi

if true ; then 
    #gdalbuildvrt -vrtnodata -9999 ${mos_prefix}_${epsg}_${outres}m.vrt $list_proj
    if [ ! -e ${mos_prefix}_${epsg}_${outres}m-tile-0-count.tif ] ; then 
        dem_mosaic --threads $threads --count $list_proj -o ${mos_prefix}_${epsg}_${outres}m
        gdaladdo_ro.sh ${mos_prefix}_${epsg}_${outres}m-tile-0-count.tif
    fi
    if $make_lowres ; then 
        if [ ! -e ${mos_prefix}_${epsg}_${outres}m-tile-0-count_${lowres}m.tif ] ; then 
            #dem_mosaic --threads $threads --tr $lowres --count $list_proj -o ${mos_prefix}_${epsg}_${lowres}m
            gdalwarp -r cubic -overwrite -dstnodata 0 -ot UInt16 -tr $lowres $lowres ${mos_prefix}_${epsg}_${outres}m-tile-0-count.tif ${mos_prefix}_${epsg}_${outres}m-tile-0-count_${lowres}m.tif
        fi
    fi
fi

if true ; then 
    if [ ! -e ${mos_prefix}_${epsg}_${outres}m-tile-0-stddev.tif ] ; then 
        dem_mosaic --threads $threads --stddev $list_proj -o ${mos_prefix}_${epsg}_${outres}m
        gdaladdo_ro.sh ${mos_prefix}_${epsg}_${outres}m-tile-0-stddev.tif
    fi
    if $make_lowres ; then 
        if [ ! -e ${mos_prefix}_${epsg}_${outres}m-tile-0-stddev_${lowres}m.tif ] ; then 
            #dem_mosaic --threads $threads --tr $lowres --stddev $list_proj -o ${mos_prefix}_${epsg}_${lowres}m
            gdalwarp -r cubic -overwrite -dstnodata 0 -tr $lowres $lowres ${mos_prefix}_${epsg}_${outres}m-tile-0-stddev.tif ${mos_prefix}_${epsg}_${outres}m-tile-0-stddev_${lowres}m.tif
        fi
    fi
fi 

if true ; then 
    compute_dh.py $lowres_dem ${mos_prefix}_${epsg}_${outres}m-tile-0.tif
    mv ${lowres_dem%.*}_${mos_prefix}_${epsg}_${outres}m-tile-0_dz_eul.tif .
    gdaladdo_ro.sh $(basename ${lowres_dem%.*}_${mos_prefix}_${epsg}_${outres}m-tile-0_dz_eul.tif)
fi

exit

#ogr2ogr -t_srs EPSG:proj ${outshp%.*}_proj.shp $outshp
#shpfltr_by_pt.sh ${outshp%.*}_proj.shp
#shpfltr_by_pt.sh $outshp

#Create 8-m tiled mosaic
#ext=DEM_8m_trans
#res=8
res=32
ext=DEM_${res}m
list=$(ls ../hma[0-9]/{QB,IK,WV,GE}*/*/*${ext}.tif)
dem_mosaic --threads $threads --tile-size 10000 --t_srs EPSG:32611 --tr $res -o mos_${res}m/mos_${res}m_all/mos_${res}m_all $list
gdalbuildvrt -vrtnodata -9999 -r cubic mos_${res}m/mos_${res}m_all/mos_${res}m_all.vrt mos_${res}m/mos_${res}m_all/mos_${res}m_all-tile-0.tif 
extent=$(get_extent mos_${res}m/mos_${res}m_all/mos_${res}m_all.vrt) 
#This won't work due to heterogenous projection
#gdalbuildvrt -vrtnodata 0 -r cubic all_${res}m.vrt $list
#extent=$(get_extent.py all_${res}m.vrt)

list_summer=$(ls ../hma[0-9]/{QB,IK,WV,GE}*/*/*${ext}.tif | grep 201[0-9][01][67890])
list_winter=$(ls ../hma[0-9]/{QB,IK,WV,GE}*/*/*${ext}.tif | grep -v 201[0-9][01][67890])

#These should all use the same extent, identical tiles
dem_mosaic --threads $threads --tile-size 10000 --t_projwin $extent --t_srs EPSG:32611 --tr 8 -o mos_${res}m/mos_${res}m_summer/mos_${res}m_summer $list_summer
dem_mosaic --threads $threads --tile-size 10000 --t_projwin $extent --t_srs EPSG:32611 --tr 8 -o mos_${res}m/mos_${res}m_winter/mos_${res}m_winter $list_winter

#Many tiles are empty, remove
mkdir valid
mv $(ll -Sr | grep -v 1.4M | grep -v K | awk '{print $NF}') valid/
cd valid
gdalbuildvrt -vrtnodata 0 -r cubic mos_${res}m_valid.vrt *tif
