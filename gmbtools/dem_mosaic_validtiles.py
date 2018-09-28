#! /usr/bin/env python
"""
Run dem_mosaic in parallel for valid tiles only
"""
#res=32
#res=8
#mkdir conus_${res}m_tile
#lfs setstripe conus_${res}m_tile --count 64
#dem_mosaic_validtiles.py --tr $res --t_projwin 'union' --t_srs '+proj=aea +lat_1=36 +lat_2=49 +lat_0=43 +lon_0=-115 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ' --georef_tile_size 100000 -o conus_${res}m_tile/conus_${res}m *00/*/*DEM_${res}m.tif

#res=8
#mkdir hma_${res}m_tile
#lfs setstripe hma_${res}m_tile --count 64
#dem_mosaic_validtiles.py --tr $res --t_projwin 'union' --t_srs '+proj=aea +lat_1=25 +lat_2=47 +lat_0=36 +lon_0=85 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ' --georef_tile_size 100000 -o hma_${res}m_tile/hma_${res}m */*00/*/*DEM_${res}m.tif

import os
import sys
import glob
import argparse
import math
import time
import subprocess
import tarfile
import pickle
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor

from osgeo import gdal, ogr

from pygeotools.lib import geolib, warplib, iolib

from dem_mosaic_index_ts import make_dem_mosaic_index_ts

#Hack to work around file open limit
#Set this in shell with `ulimit -n 65536` before running
#import resource
#resource.setrlimit(resource.RLIMIT_NOFILE,(resource.RLIM_INFINITY, resource.RLIM_INFINITY))

def getparser():
    stat_choices = ['first', 'firstindex', 'last', 'lastindex', 'min', 'max', 'mean', 'stddev', 'count', 'median', 'medianindex', 'nmad', 'wmean']
    parser = argparse.ArgumentParser(description='Wrapper for dem_mosaic that will only write valid tiles')
    parser.add_argument('--tr', default='min', help='Output resolution (default: %(default)s)')
    parser.add_argument('--t_projwin', default='union', help='Output extent (default: %(default)s)')
    parser.add_argument('--t_srs', default='first', help='Output projection (default: %(default)s)')
    parser.add_argument('--georef_tile_size', type=float, default=100000., help='Output tile width (meters)')
    parser.add_argument('--threads', type=int, default=iolib.cpu_count(logical=False), help='Number of simultaneous dem_mosaic processes to run')
    parser.add_argument('--stat', type=str, nargs='*', default=None, choices=stat_choices, \
            help='Specify space-delimited list of output statistics to pass to dem_mosaic (e..g, "count stddev", default: wmean)')
    parser.add_argument('-o', type=str, default=None, help='Output mosaic prefix')
    parser.add_argument('src_fn_list', type=str, nargs='+', help='Input filenames (img1.tif img2.tif ...)')
    return parser

def main():
    parser = getparser()
    args = parser.parse_args()

    stat_list = ['wmean',]
    if args.stat is not None:
        if isinstance(args.stat, str):
            stat_list = args.stat.split()
        else:
            stat_list = args.stat

    print("The following mosaics will be generated:")
    print(stat_list)

    #Tile dimensions in output projected units (meters)
    #Assume square
    tile_width = args.georef_tile_size
    tile_height = tile_width

    #This is number of simultaneous processes, each with one thread
    threads = args.threads

    o = args.o
    if o is None:
        o = 'mos_%im/mos' % tr
    odir = os.path.dirname(o)
    #If dirname is empty, use prefix for new directory
    if not odir:
        odir = o
        o = os.path.join(odir, o)
    if not os.path.exists(odir): 
        os.makedirs(odir)
    iolib.setstripe(odir, threads)

    out_pickle_fn = o+'_tile_dict.pkl'
    """
    if os.path.exists(out_pickle_fn):
        with open(out_pickle_fn, 'wb') as f:
            tile_dict = pickle.load(f)
    else:
    """

    #Input filelist
    fn_list = args.src_fn_list
    #Sort?
    #Might hit OS open file limit here
    print("Loading input datasets")
    print("Note: this could take several minutes depending on number of inputs and I/O performance")
    ds_list = []
    for n, fn in enumerate(fn_list):
        if (n % 100 == 0):
            print('%i of %i done' % (n, len(fn_list)))
        ds_list.append(gdal.Open(fn))

    #Mosaic t_srs
    print("\nParsing t_srs")
    t_srs = warplib.parse_srs(args.t_srs, ds_list)
    print(t_srs.ExportToProj4())

    #Mosaic res
    print("\nParsing tr")
    tr = warplib.parse_res(args.tr, ds_list, t_srs=t_srs) 
    print(tr)

    #Mosaic extent 
    #xmin, ymin, xmax, ymax
    print("Determining t_projwin (bounding box for inputs)")
    t_projwin = warplib.parse_extent(args.t_projwin, ds_list, t_srs=t_srs) 
    print(t_projwin)
    #This could trim off some fraction of a pixel around margins
    t_projwin = geolib.extent_round(t_projwin, tr)
    mos_xmin, mos_ymin, mos_xmax, mos_ymax = t_projwin

    #Compute extent geom for all input datsets
    print("Computing extent geom for all input datasets")
    input_geom_dict = OrderedDict()
    for ds in ds_list:
        ds_geom = geolib.ds_geom(ds, t_srs)
        ds_fn = ds.GetFileList()[0]
        #Could use filename as key here
        input_geom_dict[ds_fn] = geolib.geom_dup(ds_geom)
        ds = None

    ds_list = None

    #Mosaic tile size
    #Should have float extent and tile dim here
    ntiles_w = int(math.ceil((mos_xmax - mos_xmin)/tile_width))
    ntiles_h = int(math.ceil((mos_ymax - mos_ymin)/tile_height))
    ntiles = ntiles_w * ntiles_h
    print("%i (%i cols x %i rows) tiles required for full mosaic" % (ntiles, ntiles_w, ntiles_h))
    #Use this for zero-padding of tile number
    ntiles_digits = len(str(ntiles))

    print("Computing extent geom for all output tiles")
    tile_dict = OrderedDict()
    for i in range(ntiles_w):
        for j in range(ntiles_h):
            tilenum = j*ntiles_w + i
            tile_xmin = mos_xmin + i*tile_width
            tile_xmax = mos_xmin + (i+1)*tile_width
            tile_ymax = mos_ymax - j*tile_height
            tile_ymin = mos_ymax - (j+1)*tile_height
            #Corner coord needed for geom
            x = [tile_xmin, tile_xmax, tile_xmax, tile_xmin, tile_xmin]
            y = [tile_ymax, tile_ymax, tile_ymin, tile_ymin, tile_ymax]
            tile_geom_wkt = 'POLYGON(({0}))'.format(', '.join(['{0} {1}'.format(*a) for a in zip(x,y)]))
            tile_geom = ogr.CreateGeometryFromWkt(tile_geom_wkt)
            tile_geom.AssignSpatialReference(t_srs)
            #tile_dict[tilenum] = tile_geom
            tile_dict[tilenum] = {}
            tile_dict[tilenum]['geom'] = tile_geom
            tile_dict[tilenum]['extent'] = [tile_xmin, tile_ymin, tile_xmax, tile_ymax]

    print("Computing valid intersections between input dataset geom and tile geom")
    for tilenum in sorted(tile_dict.keys()):
        print('%i of %i' % (tilenum, len(tile_dict.keys())))
        tile_geom = tile_dict[tilenum]['geom']
        tile_dict_fn = []
        for ds_fn, ds_geom in input_geom_dict.iteritems():
            if tile_geom.Intersects(ds_geom):
                tile_dict_fn.append(ds_fn)
                #Write out shp for debugging
                #geolib.geom2shp(tile_geom, 'tile_%03i.shp' % tilenum)
        if tile_dict_fn:
            tile_dict[tilenum]['fn_list'] = tile_dict_fn
    
    #Write out dictionary with list of fn for each tile
    print("Writing out tile dictionary")
    with open(out_pickle_fn, 'wb') as f:
        pickle.dump(tile_dict, f)

    out_tile_list = []
    for tilenum in tile_dict.keys():
        if 'fn_list' in tile_dict[tilenum]:
            out_tile_list.append(tilenum)
        else:
            del tile_dict[tilenum]

    print("%i valid output tiles" % len(out_tile_list))
    out_tile_list.sort()
    out_tile_list = list(set(out_tile_list))
    ni = max([len(str(i)) for i in out_tile_list])
    out_tile_list_str = ' '.join(map(str, out_tile_list))
    print(out_tile_list_str)

    delay = 0.001
    outf = open(os.devnull, 'w') 
    #outf = open('%s-log-dem_mosaic-tile-%i.log' % (o, tile), 'w')

    #Should run the tiles with the largest file count first, as they will likely take longer
    tile_dict = OrderedDict(sorted(tile_dict.items(), key=lambda item: len(item[1]['fn_list']), reverse=True))
    out_tile_list = tile_dict.keys()

    #Reorganized to run all tiles for all stats in parallel
    #So we're not waiting for one tile to finish 
    with ThreadPoolExecutor(max_workers=threads) as executor:
        print("Running dem_mosaic in parallel with %i threads" % threads)
        for n, tile in enumerate(out_tile_list):
            #print('%i of %i tiles: %i' % (n+1, len(out_tile_list), tile))
            tile_fn_base = '%s-tile-%0*i.tif' % (o, ni, tile)
            for stat in stat_list:
                tile_fn = os.path.splitext(tile_fn_base)[0]+'-%s.tif' % stat
                dem_mosaic_args = {'fn_list':tile_dict[tile]['fn_list'], 'o':tile_fn, 'tr':tr, 't_srs':t_srs, \
                        't_projwin':tile_dict[tile]['extent'], 'threads':1, 'stat':stat}
                if not os.path.exists(tile_fn):
                    #This passes only files that intersect the tile, but issues with dem_mosaic reducing bounding box
                    #dem_mosaic_args[0] = tile_dict[tile]['fn_list']
                    cmd = geolib.get_dem_mosaic_cmd(**dem_mosaic_args)
                    #print(cmd)
                    executor.submit(subprocess.call, cmd, stdout=outf, stderr=subprocess.STDOUT)
            time.sleep(delay)

    #Now aggegate into stats
    for stat in stat_list:
        tile_fn_list = []
        for n, tile in enumerate(out_tile_list):
            tile_fn_base = '%s-tile-%0*i.tif' % (o, ni, tile)
            tile_fn = os.path.splitext(tile_fn_base)[0]+'-%s.tif' % stat
            if os.path.exists(tile_fn):
                tile_fn_list.append(tile_fn)
        print("\nMosaic type: %s" % stat)
        #Convert dem_mosaic index files to timestamp arrays
        if stat in ['lastindex', 'firstindex', 'medianindex']:
            print("Running dem_mosaic_index_ts in parallel with %i threads" % threads)
            from multiprocessing import Pool
            pool = Pool(processes=threads)
            results = pool.map(make_dem_mosaic_index_ts, tile_fn_list)
            pool.close()
            #results.wait()
            #Update filenames
            tile_fn_list = [os.path.splitext(tile_fn)[0]+'_ts.tif' for tile_fn in tile_fn_list]

        print("\nCreating vrt of valid tiles")
        #tile_fn_list = glob.glob(o+'-tile-*.tif')
        vrt_fn = o+'.vrt'
        if stat is not None:
            vrt_fn = os.path.splitext(vrt_fn)[0]+'_%s.vrt' % stat
            if stat in ['lastindex', 'firstindex', 'medianindex']:
                vrt_fn = os.path.splitext(vrt_fn)[0]+'_ts.vrt'
        cmd = ['gdalbuildvrt', vrt_fn] 
        vrt_fn_list = []
        for tile_fn in tile_fn_list:
            if os.path.exists(tile_fn):
                vrt_fn_list.append(tile_fn)
            else:
                print("Missing file: %s" % tile_fn)
        cmd.extend(vrt_fn_list)
        print(cmd)
        subprocess.call(cmd)

        #Should create tile index shp/kml from tile_geom

        #This cleans up all of the log txt files (potentially 1000s of files)
        #Want to preserve these, as they contain list of DEMs that went into each tile
        log_fn_list = glob.glob(o+'*%s.tif-log-dem_mosaic-*.txt' % stat)
        print("\nCleaning up %i dem_mosaic log files" % len(log_fn_list))
        if stat is not None:
            tar_fn = o+'_%s_dem_mosaic_log.tar.gz' % stat
        else:
            tar_fn = o+'_dem_mosaic_log.tar.gz'
        with tarfile.open(tar_fn, "w:gz") as tar:
            for log_fn in log_fn_list:
                tar.add(log_fn)
        for log_fn in log_fn_list:
            os.remove(log_fn)

    outf = None

if __name__ == "__main__":
    main()
