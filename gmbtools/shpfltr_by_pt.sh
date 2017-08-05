#! /bin/bash

#Utility to isolate overlapping products for a particular site lat/lon

#Generate master shp
#raster2shp.py conus*/WV*/dem*/*DEM_32m.tif
#shp_rename.sh merge.shp conus_12_DEM_32m_20160825.shp
#ogr2ogr -t_srs EPSG:32611 conus_12_DEM_32m_20160825_32611.shp conus_12_DEM_32m_20160825.shp

#Prepare point coord for sites
#QGIS export csv with both WGS84 and UTM 11N (32611)
#cut -d, -f4,5,6 shean_site_pts_20160928_geom.csv > shean_site_pts_20160928_geom_out.csv
#cut -d, -f5,6 shean_site_pts_20160928_geom_32611.csv | paste -d, shean_site_pts_20160928_geom_out.csv - > shean_site_pts_20160928_geom_out_32611.csv

#Should really use input polygons, but points with buffer will work for now

#proj=32611
proj=32610

shp_fn=$1
lyr=${shp_fn##*/}
lyr=${lyr%.*}
site_coord_fn=shean_site_pts_20160928_geom_out_32611.csv
glac_shp=/nobackupp8/deshean/rpcdem/nlcd/24k_selection_32611.shp
outdir=sites

neddir=/nobackup/deshean/rpcdem/ned13
ned=$neddir/ned13_tiles_glac24k_115kmbuff.vrt
nedres=10
neddir=/nobackup/deshean/rpcdem/ned1
ned=$neddir/ned1_tiles_glac24k_115kmbuff.vrt
nedres=30

if [ ! -d $outdir ] ; then 
    mkdir $outdir
fi

#IFS=,
#while read -r site lon lat x y ; do
#    echo $site $lon $lat $x $y
#done < $site_coord_fn

#for line in $(cat $site_coord_fn) 
for line in $(grep scg $site_coord_fn) 
do
    site=$(echo $line | cut -d',' -f1)
    lon=$(echo $line | cut -d',' -f2)
    lat=$(echo $line | cut -d',' -f3)
    x=$(echo $line | cut -d',' -f4)
    y=$(echo $line | cut -d',' -f5)
    echo $site $lon $lat $x $y
    if [ ! -d $outdir/$site ] ; then 
        mkdir $outdir/$site
    fi
    #If input shp is WGS84
    #pad=0.000001
    #pad=0.07
    #If input is in projected coord
    pad=8000
    xmin=$(echo "$x - $pad" | bc -l)
    ymin=$(echo "$y - $pad" | bc -l)
    xmax=$(echo "$x + $pad" | bc -l)
    ymax=$(echo "$y + $pad" | bc -l)
    spat="$xmin $ymin $xmax $ymax"
    echo $spat
    #spat='-121.814399 48.776699 -121.814401 48.776701'
    out_fn=$outdir/$site/${lyr}_${site}_fltr.shp
    eval ogr2ogr -overwrite -spat $spat $out_fn $shp_fn 

    #Get filename list
    #ogrinfo -al $out_fn | grep name | awk '{print $NF}' | sed '1,2d' | sed 's/$/.tif/' > ${out_fn%.*}_fnlist.txt
    ogrinfo -al $out_fn | grep name | awk '{print $NF}' | sed '1,2d' > ${out_fn%.*}_fnlist.txt

    if [ ! -z ${out_fn%.*}_fnlist.txt ] ; then
        pushd $outdir/$site
        fnlist=$(basename ${out_fn%.*}_fnlist.txt)
        echo -n > ${fnlist%.*}_ln.txt
        for i in $(cat $fnlist) 
        do
            ln -sv ../../$i.tif . 
            echo $(basename $i).tif >> ${fnlist%.*}_ln.txt
        done
        pushd
    fi
done

cd $outdir
#SCG
yr=19580901
#extent='195873 5356342 206489 5367251'

#cat conus_32611_100m_n366_20160928_1531_scg_fltr_fnlist.txt | sed 's/DEM_32m_trans_32611_100m/DEM_8m_trans.tif/' | sed 's#^#../../#' > conus_32611_100m_n366_20160928_1531_scg_fltr_fnlist_8m.txt
#make_stack.py -outdir . -te '639187 5348720 651696 5362850' -tr 8 $(cat conus_32611_100m_n366_20160928_1531_scg_fltr_fnlist_8m.txt)
#make_stack.py -outdir . -te '639187 5348720 651696 5362850' -tr 8 $(cat conus_32611_100m_n366_20160928_1531_scg_fltr_fnlist_8m.txt | grep 201[345][01][7890])

for site in $(ls -d ./*)
do
    pushd $site
    site_fn_list=$(cat *ln.txt)
    extent=$(compute_union.py $site_fn_list)
    res=$(echo $site_fn_list | cut -f1)
    #warptool.py -te $extent -tr $res -outdir $outdir/$site $ned 
    ln -s $ned ${yr}_$(basename $ned)
    make_stack.py --no-datestack -te "$extent" -tr $res $site_fn_list ${yr}_$(basename $ned) 
    pushd
done

exit

#ogr2ogr -f "ESRI Shapefile" $out_fn $shp_fn -dialect sqlite -sql "SELECT * from $lyr WHERE ST_Intersection(POINT($x,$y), GEOMETRY)"
#sql_stmt="SELECT \* FROM $lyr WHERE ST_Intersects(GEOMETRY, \'POINT($x $y)\'::geometry)"
#ogr2ogr -progress -overwrite -sql "$sql_stmt" -dialect SQLITE -f 'ESRI Shapefile' $out_fn $shp_fn 
