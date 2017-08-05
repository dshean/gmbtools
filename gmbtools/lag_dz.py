#! /usr/bin/env python

#David Shean
#dshean@gmail.com

#This script takes in two dems and a displacement field to compute elevation differences

#Input dem1 dem2 [dem1_dem2-F.tif]

#To do:
#Clean up stats computation using new functions
#Check to make sure files are clipped to same extent - if not, clip.  Move the clipping function to a central library.
#Fix vrt issues - always write absolue path to vrt
#Move writing out files to new function
#Move each elev diff approach to new function

import sys
import os
import re
import subprocess
import glob

import gdal
import numpy as np

from pygeotools.lib import iolib
from pygeotools.lib import malib
from pygeotools.lib import geolib
from pygeotools.lib import warplib
from pygeotools.lib import filtlib
#from lib import glaclib

#This is output ndv, avoid using 0 for differences
diffndv = -9999

dem1_fn = sys.argv[1]
dem2_fn = sys.argv[2]

if dem1_fn == dem2_fn:
    sys.exit('Input filenames are identical')

fn_list = [dem1_fn, dem2_fn]

#This converts external velocity grids to pixel disparities
#Lingering issues for inputs with variable timestamps
vel_input = False 

#Note: velocity is expected to be a 2-band file
#for i in *mos_vx.tif ; do
#base=$(echo $i | awk -F'_v' '{print $1}')
#gdal_merge.py -separate -a_nodata -2000000000 -o ${base}_vxy.tif ${base}_vx.tif ${base}_vy.tif 
#done

print "Warping DEMs to same res/extent/proj"
#Might want to limit intersection to dem1 and dem2, as vmap could be significantly smaller
disp_ds = None
if len(sys.argv) == 4:
    disp_fn = sys.argv[3]
    if 'track' in disp_fn or 'tsx' in disp_fn:
        vel_input = True
    fn_list.append(disp_fn)
    dem1_ds, dem2_ds, disp_ds = warplib.memwarp_multi_fn(fn_list, extent='intersection', res='max')
else:
    dem1_ds, dem2_ds = warplib.memwarp_multi_fn(fn_list, extent='intersection', res='max')

outdir = os.path.split(dem1_fn)[0]
outprefix = os.path.splitext(os.path.split(dem1_fn)[1])[0]+'_'+os.path.splitext(os.path.split(dem2_fn)[1])[0]

#Load input DEMs into masked arrays
print "Loading input DEMs into masked arrays"
dem1 = iolib.ds_getma(dem1_ds, 1)
dem2 = iolib.ds_getma(dem2_ds, 1)

print

adj = ''
if '-adj' in dem1_fn:
    adj = '-adj' 
dem1_fn_base = re.sub(adj, '', os.path.splitext(dem1_fn)[0]) 
dem2_fn_base = re.sub(adj, '', os.path.splitext(dem2_fn)[0]) 

#dem1_geoid = geolib.dem_geoid(dem1_fn)
#dem2_geoid = geolib.dem_geoid(dem1_fn)

#Want to move this timestamp array handling to timelib
#Need to handle case where one ts array is present and other is constant ts
rates = True 
#Scale output to m/yr
if rates:
    #Attempt to load timestamp arrays (for mosaics)
    t1_fn = dem1_fn_base+'_ts.tif'
    t2_fn = dem2_fn_base+'_ts.tif'
    if os.path.exists(t1_fn) and os.path.exists(t2_fn):
        constant_dt = False
        print "Preparing timestamp arrays"
        t1_ds, t2_ds = warplib.memwarp_multi_fn([t1_fn, t2_fn], extent=dem1_ds, res=dem1_ds)
        print "Loading timestamps into masked arrays"
        t1 = iolib.ds_getma(t1_ds)
        t2 = iolib.ds_getma(t2_ds)
        #Compute dt in days
        t_factor = t2 - t1
        t_factor /= 365.25
    else:
        from datetime import datetime, timedelta
        from pygeotools.lib import timelib 
        t1 = timelib.fn_getdatetime(dem1_fn)
        t2 = timelib.fn_getdatetime(dem2_fn)
        if t1 is not None and t2 is not None and t1 != t2:  
            constant_dt = True
            dt = t2 - t1
            #Might be better to do this with dateutil - not sure about leap years
            #from dateutil.relativedelta import relativedelta
            #dt = relativedelta(dt1, dt2)) 
            #dt.years
            year = timedelta(days=365.25)
            t_factor = abs(dt.total_seconds()/year.total_seconds()) 
            print "Time differences is %s, dh/%0.3f" % (dt, t_factor)
        else:
            print "Unable to extract timestamps for input images"
            rates = False

#Less restrictive basename
#dem1_fn_base = os.path.splitext(dem1_fn)[0].split('_mos_')[0]+'_mos'
#dem2_fn_base = os.path.splitext(dem2_fn)[0].split('_mos_')[0]+'_mos'
#This works for SPIRIT products as well
dem1_fn_base = os.path.join(os.path.split(dem1_fn)[0], os.path.split(dem1_fn)[-1][0:13]) 
dem2_fn_base = os.path.join(os.path.split(dem2_fn)[0], os.path.split(dem2_fn)[-1][0:13]) 

def findfile(base_fn, ext):
    fn = glob.glob(base_fn+'*'+ext)
    if not fn:
        fn = None
    else:
        fn = fn[0]
    return fn

tidecorr = True 
if tidecorr:
    #dem1_tide_fn = dem1_fn_base+'_tidecorr.tif'
    #dem2_tide_fn = dem2_fn_base+'_tidecorr.tif'
    dem1_tide_fn = findfile(dem1_fn_base, 'tidecorr.tif')
    if not dem1_tide_fn:
        dem1_tide_fn = findfile(dem1_fn_base, 'tidemodel_smooth_full_clip.tif')
    dem2_tide_fn = findfile(dem2_fn_base, 'tidecorr.tif')
    if not dem2_tide_fn:
        dem2_tide_fn = findfile(dem2_fn_base, 'tidemodel_smooth_full_clip.tif')
    if dem1_tide_fn and dem2_tide_fn:
        dem1_tide_ds, dem2_tide_ds = warplib.memwarp_multi_fn([dem1_tide_fn, dem2_tide_fn], extent=dem1_ds, res=dem1_ds)
        dem1_tide = iolib.ds_getma(dem1_tide_ds)
        dem2_tide = iolib.ds_getma(dem2_tide_ds)
        tide_diff = dem2_tide - dem1_tide
        dst_fn = os.path.join(outdir, outprefix+'_tide_diff.tif')
        iolib.writeGTiff(tide_diff, dst_fn, dem1_ds, ndv=diffndv)
        #These values are tide prediction, to remove, want to subtract from observed elevation
        #Need to fill with 0 to prevent clipping to floating ice
        dem1 -= dem1_tide.filled(0)
        dem2 -= dem2_tide.filled(0)

firnair = True 
#This is constant value for PIG
dem1_firnair = 15.0
dem2_firnair = dem1_firnair 
firnair_diff = dem2_firnair - dem1_firnair 
if firnair:
    #dem1_firnair_fn = dem1_fn_base+'_racmo_FirnAir.tif'
    #dem2_firnair_fn = dem2_fn_base+'_racmo_FirnAir.tif'
    dem1_firnair_fn = findfile(dem1_fn_base, 'racmo_FirnAir.tif')
    dem2_firnair_fn = findfile(dem2_fn_base, 'racmo_FirnAir.tif')
    if dem1_firnair_fn and dem2_firnair_fn:
        dem1_firnair_ds, dem2_firnair_ds = warplib.memwarp_multi_fn([dem1_firnair_fn, dem2_firnair_fn], extent=dem1_ds, res=dem1_ds)
        dem1_firnair = iolib.ds_getma(dem1_firnair_ds)
        dem2_firnair = iolib.ds_getma(dem2_firnair_ds)
        firnair_diff = dem2_firnair - dem1_firnair 
        dst_fn = os.path.join(outdir, outprefix+'_FirnAir_diff.tif')
        iolib.writeGTiff(firnair_diff, dst_fn, dem1_ds, ndv=diffndv)
        #These values are positive, total firn air content, want to subtract
        #dem1 -= dem1_firnair
        #dem2 -= dem2_firnair

zs = True 
zs_diff = 0
if zs:
    #dem1_zs_fn = dem1_fn_base+'_racmo_zs.tif'
    #dem2_zs_fn = dem2_fn_base+'_racmo_zs.tif'
    dem1_zs_fn = findfile(dem1_fn_base, 'racmo_zs.tif')
    dem2_zs_fn = findfile(dem2_fn_base, 'racmo_zs.tif')
    if dem1_zs_fn and dem2_zs_fn:
        dem1_zs_ds, dem2_zs_ds = warplib.memwarp_multi_fn([dem1_zs_fn, dem2_zs_fn], extent=dem1_ds, res=dem1_ds)
        dem1_zs = iolib.ds_getma(dem1_zs_ds)
        dem2_zs = iolib.ds_getma(dem2_zs_ds)
        zs_diff = dem2_zs - dem1_zs
        dst_fn = os.path.join(outdir, outprefix+'_zs_diff.tif')
        iolib.writeGTiff(zs_diff, dst_fn, dem1_ds, ndv=diffndv)
        smb_diff = zs_diff - firnair_diff 
        #dst_fn = os.path.join(outdir, outprefix+'_smb_diff.tif')
        #iolib.writeGTiff(smb_diff, dst_fn, dem1_ds, ndv=diffndv)

#Check to make sure inputs actually intersect
#Masked pixels are True
if not np.any(~dem1.mask*~dem2.mask):
    sys.exit("No valid overlap between input data")

#Compute common mask
print "Generating common mask"
common_mask = malib.common_mask([dem1, dem2])

#Compute relative elevation difference with Eulerian approach 
print "Computing elevation difference with Eulerian approach"
diff_euler = np.ma.array(dem2-dem1, mask=common_mask)

#Add output that has difference filter applied
#See dem_align

if True:
    print "Eulerian elevation difference stats:"
    diff_euler_stats = malib.print_stats(diff_euler)
    diff_euler_med = diff_euler_stats[5]

if True:
    print "Writing Eulerian elevation difference map"
    dst_fn = os.path.join(outdir, outprefix+'_dz_eul.tif')
    print dst_fn
    iolib.writeGTiff(diff_euler, dst_fn, dem1_ds, ndv=diffndv)
    if rates:
        print "Writing Eulerian rate map"
        dst_fn = os.path.join(outdir, outprefix+'_dz_eul_rate.tif')
        print dst_fn
        iolib.writeGTiff(diff_euler/t_factor, dst_fn, dem1_ds, ndv=diffndv)

if False:
    print "Writing Eulerian relative elevation difference map"
    diff_euler_rel = diff_euler - diff_euler_med
    dst_fn = os.path.join(outdir, outprefix+'_dz_eul_rel.tif')
    print dst_fn
    iolib.writeGTiff(diff_euler_rel, dst_fn, dem1_ds, ndv=diffndv)

if False:
    print "Writing out DEM2 with median elevation difference removed"
    dst_fn = os.path.splitext(dem2_fn)[0]+'_med'+diff_euler_med+'.tif'
    print dst_fn
    iolib.writeGTiff(dem2 - diff_euler_med, dst_fn, dem1_ds, ndv=diffndv)

if False:
    print "Writing Eulerian elevation difference percentage map"
    diff_euler_perc = 100.0*diff_euler/dem1
    dst_fn = os.path.join(outdir, outprefix+'_dz_eul_perc.tif')
    print dst_fn
    iolib.writeGTiff(diff_euler_perc, dst_fn, dem1_ds, ndv=diffndv)

#def compute_dh_vs_z(ref_dem, dh, nbins=20, binwidth=None):
if False:
    print "Compute dh/dt vs. elevation"
    nbins = 20
    #binwidth = None
    binwidth = 50
    ref_dem = dem1
    dh = diff_euler
    if rates:
        dh = diff_euler/t_factor
    min, max = malib.calcperc(ref_dem)
    #If binwidth specified, override
    if binwidth is not None:
        nbins = None
    if nbins:
        edges = np.linspace(min, max, num=nbins)
        binwidth = edges[1] - edges[0]
    else:
        edges = np.arange(min, max, binwidth)
    diff_euler_bincenters = []
    diff_euler_binvals = []
    diff_euler_bincount = []
    for i in range(edges.size-1):
        print "%i of %i: %0.1f to %0.1f m" % (i, edges.size - 1, edges[i], edges[i+1])
        idx = np.logical_and((ref_dem > edges[i]).data, (ref_dem <= edges[i+1]).data)
        bincenter = edges[i]+((edges[i+1]-edges[i])/2.0)
        diff_euler_bincenters.append(bincenter)
        dh_idx = dh[idx]
        diff_euler_binvals.append(np.median(dh_idx))
        diff_euler_bincount.append(dh_idx.count())
    diff_euler_binvals = np.ma.array(diff_euler_binvals)
    diff_euler_bincount = np.ma.array(diff_euler_bincount)
    res = geolib.get_res(dem1_ds, square=True)[0]
    diff_euler_binarea = (diff_euler_bincount*res**2)/1E6
    #Threshold area for bin in km^2
    binarea_thresh = 0.1
    binarea_mask = ~((diff_euler_binarea > binarea_thresh).data)
    diff_euler_binvals = np.ma.array(diff_euler_binvals, mask=binarea_mask)
    diff_euler_binarea = np.ma.array(diff_euler_binarea, mask=binarea_mask)
    edges = np.ma.array(edges[:-1], mask=binarea_mask)

    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(5, 10))
    fig.add_subplot(211)
    #Plot negative in Red, postive in Blue
    plt.bar(edges, diff_euler_binvals, width=binwidth, color='k')
    if rates:
        plt.ylabel('Median Elevation Change Rate (m/yr)')
    else:
        plt.ylabel('Median Elevation Change (m)')
    plt.xlabel('Elevation (m)')
    fig.add_subplot(212)
    plt.bar(edges, diff_euler_binarea, width=binwidth, color='k')
    plt.ylabel('Area (km^2)')
    plt.xlabel('Elevation (m)')
    dst_fn = os.path.join(outdir, outprefix+'_dz_vs_z.pdf')
    plt.tight_layout()
    fig.savefig(dst_fn)
    plt.show()

if disp_ds is not None:
    #Load disparity values into masked arrays
    #Note: need to represent these as integers, as they will be used as indices for an array
    print "Loading disparity maps into masked arrays"
    disp_x_f_in = iolib.ds_getma(disp_ds, 1)
    disp_y_f_in = iolib.ds_getma(disp_ds, 2)

    #If the input disp_ds (not velocities) has been resampled from original pixel dimensions, must alter values
    #disp_scaling = 0.57/32.0
    #disp_x_f_in *= disp_scaling 
    #disp_y_f_in *= disp_scaling 

    disp_mask = malib.common_mask([disp_x_f_in, disp_y_f_in])
    disp_x_f = np.ma.array(disp_x_f_in, mask=disp_mask)
    disp_y_f = np.ma.array(disp_y_f_in, mask=disp_mask)
    disp_m_f = np.ma.sqrt(disp_x_f**2 + disp_y_f**2)

    x_res, y_res = geolib.get_res(disp_ds)
    srs = geolib.get_ds_srs(disp_ds)
    proj_scale_factor = 1.0
    if srs.IsSame(geolib.nps_srs) or srs.IsSame(geolib.sps_srs):
        proj_scale_factor = geolib.scale_ps_ds(disp_ds)
    print "Projection scale factor: ", proj_scale_factor

    #This is an ugly hack to iterate when external velocity map is used - want to clean up
    t_factor_lag = t_factor
    if vel_input:
        niter = 4
    else:
        niter = 1
    while niter > 0:
        #Convert input velocity grids to pixel disparities
        if vel_input:
            #Note minus sign for disp_y_f
            #Ian's products have y positive upward (map coordinates)
            #ASP maps have y postive downward (image coordinates)
            #!!!!!
            #Note: applying t_factor here is in eulerial reference frame (common dem1+dem2 mask)
            #Need to know pixel offsets up front to compute time difference at each pixel and scale m/yr to px
            #If t1 and t2 are constant, no problem, but for mosaics w/ different timestamps this is problematic
            #Temporary hack is to use eulerian t_factor as guess, then iterate
            #Can also propagate day by day until t2 day matches
            #!!!!!
            disp_x_f = disp_x_f_in * t_factor_lag/(x_res*proj_scale_factor)
            disp_y_f = -disp_y_f_in * t_factor_lag/(y_res*proj_scale_factor)
            #disp_mask_dem2 = malib.common_mask([disp_mask, dem2.mask])
            #disp_x_f.mask = disp_mask_dem2
            #disp_y_f.mask = disp_mask_dem2
            disp_mask = malib.common_mask([disp_x_f, disp_y_f])
            disp_x_f = np.ma.array(disp_x_f, mask=disp_mask)
            disp_y_f = np.ma.array(disp_y_f, mask=disp_mask)
            disp_m_f = np.ma.sqrt(disp_x_f**2 + disp_y_f**2)

        #Convert to integer indices
        disp_x = np.ma.around(disp_x_f).astype('int16')
        disp_y = np.ma.around(disp_y_f).astype('int16')

        #(dx,dy) are values in disparity map
        #(x0,y0) are original locations in dem1
        #(x1,y1) are displaced locations of same feature in dem2
        #(x1,y1) = (x0+dx,y0+dy)

        #Generate index arrays - these are needed to go from relative disparities to absolute indices
        idx = np.indices(dem1.shape, dtype='int16')

        #Compute absolute indices for relative disparity offsets
        print "Computing offset indices using disparity values"

        #This stores indices for expected feature locations in dem2
        #Clip to dimensions of input DEM
        #Apply dem1 mask, since everything is in dem1 coords, must have valid elevations
        disp_x_abs_fwd = np.ma.masked_outside(idx[1]+disp_x, 0, dem1.shape[1]-1)
        disp_x_abs_fwd = np.ma.array(disp_x_abs_fwd, mask=dem1.mask)
        disp_y_abs_fwd = np.ma.masked_outside(idx[0]+disp_y, 0, dem1.shape[0]-1)
        disp_y_abs_fwd = np.ma.array(disp_y_abs_fwd, mask=dem1.mask)
        if np.any(disp_x_abs_fwd.mask != disp_y_abs_fwd.mask):
            disp_fwd_mask = malib.common_mask([disp_x_abs_fwd, disp_y_abs_fwd])
            disp_x_abs_fwd.mask[:] = disp_fwd_mask
            disp_y_abs_fwd.mask[:] = disp_fwd_mask

        #Compute time differences for actual displacements
        if constant_dt:
            t_factor_lag = t_factor
            niter = 0
        else:
            t2_fwd = t2[np.ma.filled(disp_y_abs_fwd,0),np.ma.filled(disp_x_abs_fwd,0)]
            t_factor_lag = (t2_fwd - t1)/365.25
            niter -= 1
            print "%i iterations remain" % niter

    print "Computing velocities from input disparities"
    #This contains x and y components
    #dem1_floatation = glaclib.freeboard_thickness(dem1)
    #dem1_grad = np.gradient(dem1)
    #dem1_floatation_grad = np.gradient(dem1_floatation)
    u = disp_x_f*x_res*proj_scale_factor/t_factor_lag
    v = disp_y_f*y_res*proj_scale_factor/t_factor_lag
    v_mag = np.ma.sqrt(u**2 + v**2) 
    #v_mag_f = malib.robust_spread_fltr(v_mag)
    v_stats = malib.print_stats(v_mag_f)
    dst_fn = os.path.join(outdir, outprefix+'_vmag.tif')
    iolib.writeGTiff(v_mag_f.astype(float), dst_fn, dem1_ds, ndv=diffndv)
    #dem1_grad_adv = disp_x_abs_fwd*x_res * dem1_grad[1] + disp_y_abs_fwd*y_res * dem1_grad[0]
    #dem1_grad_adv = u * dem1_grad[1] + v * dem1_grad[0]
    #dem1_grad_adv = u*t_factor * dem1_grad[1] + v*t_factor * dem1_grad[0]
    #dem1_floatation_grad_adv = u*t_factor * dem1_floatation_grad[1] + v*t_factor * dem1_floatation_grad[0]

    #fwd_idx = (~disp_x_abs_fwd.mask).nonzero()
    #dem1_fwd_obs = np.ma.empty_like(dem1)
    #dem1_fwd_obs[fwd_idx] = dem2[disp_y_abs_fwd.compressed(),disp_x_abs_fwd.compressed()]
   
    #Extract values from dem2 for new positions
    dem1_fwd_obs = dem2[np.ma.filled(disp_y_abs_fwd,0),np.ma.filled(disp_x_abs_fwd,0)]
    dem1_fwd_obs = np.ma.masked_equal(dem1_fwd_obs, dem2.fill_value)
    #This is necessary b/c we are filling index arrays with 0, so existing nodata pixels will have valid values

    #NOTE: this was changed on 9/25/15 due to error in the Jak datset, but it didn't cause problems with earlier PIG data

    disp_m = np.percentile(np.sqrt(disp_x**2 + disp_y**2), 98)
    #disp_m = 59
    disp_m = 29

    dem1_smooth = filtlib.gauss_fltr_astropy(dem1,2*disp_m) 
    dem2_smooth = filtlib.gauss_fltr_astropy(dem2,2*disp_m) 

    #dem1_fwd_obs.mask[:] = disp_x_abs_fwd.mask
    dem1_fwd_obs = np.ma.array(dem1_fwd_obs, mask=disp_x_abs_fwd.mask)
    #Expected elevation values assuming steady-state advection along dem1 surface
    #dem1_fwd_exp = dem1[np.ma.filled(disp_y_abs_fwd,0),np.ma.filled(disp_x_abs_fwd,0)]
    dem1_fwd_exp = dem1_smooth[np.ma.filled(disp_y_abs_fwd,0),np.ma.filled(disp_x_abs_fwd,0)]
    dem1_fwd_exp = np.ma.masked_equal(dem1_fwd_exp, dem1.fill_value)
    #dem1_fwd_exp.mask[:] = disp_x_abs_fwd.mask
    dem1_fwd_exp = np.ma.array(dem1_fwd_exp, mask=disp_x_abs_fwd.mask)

    #This stores indices for expected feature locations in dem1
    disp_x_abs_inv = np.ma.masked_outside(idx[1]-disp_x, 0, dem1.shape[1]-1)
    disp_y_abs_inv = np.ma.masked_outside(idx[0]-disp_y, 0, dem1.shape[0]-1)
    if np.any(disp_x_abs_inv.mask != disp_y_abs_inv.mask):
        disp_inv_mask = malib.common_mask([disp_x_abs_inv, disp_y_abs_inv])
        disp_x_abs_inv.mask[:] = disp_inv_mask
        disp_y_abs_inv.mask[:] = disp_inv_mask

    #Extract values from dem1 for inverse disparities
    dem2_inv_obs = dem1[np.ma.filled(disp_y_abs_inv,0),np.ma.filled(disp_x_abs_inv,0)]
    dem2_inv_obs = np.ma.masked_equal(dem2_inv_obs, dem1.fill_value)
    dem2_inv_obs.mask[:] = disp_x_abs_inv.mask
    #Expected elevation values assuming steady-state advection along dem2 surface
    dem2_inv_exp= dem2[np.ma.filled(disp_y_abs_inv,0),np.ma.filled(disp_x_abs_inv,0)]
    dem2_inv_exp = np.ma.masked_equal(dem2_inv_exp, dem2.fill_value)
    dem2_inv_exp.mask[:] = disp_x_abs_inv.mask

    #Maintain all grids in dem1 coordinates
    diff_lag_fwd_obs = dem1_fwd_obs - dem1
    diff_lag_fwd_exp = dem1_fwd_exp - dem1
    #lag residual dh/dt for advection
    diff_fwd_exp = dem1_fwd_exp - dem1_smooth
    #diff_fwd_resid = dem1_fwd_obs - dem1_fwd_exp 
    diff_fwd_resid = diff_lag_fwd_obs - diff_fwd_exp
    diff_euler_resid = diff_euler - diff_fwd_resid 

    #lag residual dh/dt for translation
    diff_lag_inv_obs = dem2 - dem2_inv_obs
    diff_lag_inv_exp = dem2 - dem2_inv_exp
    diff_inv_resid = dem2_inv_obs - dem2_inv_exp 

    #import ipdb; ipdb.set_trace()

    #Now 
    diff_lag_inv_obs_orig = diff_lag_inv_obs[np.ma.filled(disp_y_abs_fwd,0),np.ma.filled(disp_x_abs_fwd,0)]
    #diff_lag_inv_obs_orig.mask[:] = diff_lag_inv_obs.mask[np.ma.filled(disp_y_abs_fwd,0),np.ma.filled(disp_x_abs_fwd,0)] 
    #diff_lag_fwd_obs_orig = diff_lag_fwd_obs[np.ma.filled(disp_y_abs_inv,0),np.ma.filled(disp_x_abs_inv,0)]
    #diff_lag_fwd_obs_orig.mask[:] = diff_lag_fwd_obs.mask[np.ma.filled(disp_y_abs_inv,0),np.ma.filled(disp_x_abs_inv,0)] 

    #This compares forward and inverse differences in dem1 coordinates
    #Should be 0 if consistent
    #Should have limited extent (masks propagate)
    diff_fwd_inv_dem1 = diff_lag_fwd_obs - diff_lag_inv_obs_orig 
    #diff_fwd_inv_dem2 = diff_lag_inv_obs - diff_lag_fwd_obs_orig 
    
    if True:
        print "Lagrangian elevation difference stats:"
        diff_lag_stats = malib.print_stats(diff_lag_fwd_obs)
        diff_lag_fwd_obs_med = diff_lag_stats[5]

    if True:
        print "Writing Lagrangian elevation difference map"
        dst_fn = os.path.join(outdir, outprefix+'_dz_lag.tif')
        print dst_fn
        iolib.writeGTiff(diff_lag_fwd_obs, dst_fn, dem1_ds, ndv=diffndv)
        if rates:
            print "Writing Lagrangian rate map"
            dst_fn = os.path.join(outdir, outprefix+'_dz_lag_rate.tif')
            print dst_fn
            iolib.writeGTiff(diff_lag_fwd_obs/t_factor_lag, dst_fn, dem1_ds, ndv=diffndv)
    
    if False:
        print "Writing Lagrangian relative elevation difference map"
        diff_lag_fwd_obs_rel = diff_lag_fwd_obs - diff_lag_fwd_obs_med
        dst_fn = os.path.join(outdir, outprefix+'_dz_lag_rel.tif')
        iolib.writeGTiff(diff_lag_fwd_obs_rel, dst_fn, dem1_ds, ndv=diffndv)

    if True:
        print "Writing predicted elevation difference map for advection"
        dst_fn = os.path.join(outdir, outprefix+'_dz_lag_adv_exp.tif')
        iolib.writeGTiff(diff_lag_fwd_exp, dst_fn, dem1_ds, ndv=diffndv)
    
    if True:
        print "Writing residual elevation difference map for advection"
        dst_fn = os.path.join(outdir, outprefix+'_dz_lag_adv_resid.tif')
        iolib.writeGTiff(diff_fwd_resid, dst_fn, dem1_ds, ndv=diffndv)

    if True:
        print "Writing residual elevation difference map for advection"
        dst_fn = os.path.join(outdir, outprefix+'_dz_lag_adv_exp_smooth.tif')
        iolib.writeGTiff(diff_fwd_exp, dst_fn, dem1_ds, ndv=diffndv)

    if True:
        print "Writing residual eulerian difference using predicted eulerian" 
        dst_fn = os.path.join(outdir, outprefix+'_dz_eul_adv_resid.tif')
        iolib.writeGTiff(diff_euler_resid, dst_fn, dem1_ds, ndv=diffndv)
    
    if True:
        print "Computing velocity divergence"
        #disp_x_div = np.gradient(disp_x_f)[1]
        #disp_y_div = np.gradient(disp_y_f)[0]
        #disp_x_div = (disp_x_f[np.ma.filled(disp_y_abs_fwd,0),np.ma.filled(disp_x_abs_fwd,0)] - disp_x_f)/disp_x_f
        #disp_y_div = (disp_y_f[np.ma.filled(disp_y_abs_fwd,0),np.ma.filled(disp_x_abs_fwd,0)] - disp_y_f)/disp_y_f
        #vdiv = disp_x_div + disp_y_div
        #u_div = np.gradient(u/x_res)[1]
        #v_div = np.gradient(v/y_res)[0]
        #if Ux is increasing, dUx/dx should be positive
        #u_div = (u - u[np.ma.filled(disp_y_abs_fwd,0),np.ma.filled(disp_x_abs_fwd,0)])/(disp_x_f*x_res)
        #v_div = (v - v[np.ma.filled(disp_y_abs_fwd,0),np.ma.filled(disp_x_abs_fwd,0)])/(disp_y_f*y_res)
        u_div = (u[np.ma.filled(disp_y_abs_fwd,0),np.ma.filled(disp_x_abs_fwd,0)] - u)/(disp_x_f*x_res)
        v_div = (v[np.ma.filled(disp_y_abs_fwd,0),np.ma.filled(disp_x_abs_fwd,0)] - v)/(disp_y_f*y_res)
        vdiv = u_div + v_div
        vdiv_f = malib.robust_spread_fltr(vdiv)
        #Could probably just filter by -0.2 < vdiv < 0.2
        vdiv_stats = malib.print_stats(vdiv_f)
        dst_fn = os.path.join(outdir, outprefix+'_vdiv.tif')
        iolib.writeGTiff(vdiv_f.astype(float), dst_fn, dem1_ds, ndv=diffndv)

    #Want to clip these products to floating regions 
    #Compute with bed?

    if False:
        print "Computing basal melt rate"
        #Calculate height to use for floatation assumption (remove tide, remove firnair)
        #dem1_corr = dem1 - dem1_tide.filled(0) - dem1_firnair
        #Now remove tides up front
        dem1_corr = dem1 - dem1_firnair
        #Compute ms - height change to due to smb and firn 
        #If positive, net accumulation from t1 to t2
        #Should probalby compute this using lagrangian ref frame 
        #ms = dem2_zs - dem1_zs
        ms = zs_diff
        #Ice flow divergence H*(dux/dx + dvy/dy)
        #If positive, net outward flux from box
        #We are assuming hydrostatic equilibrium here, so OK to use height rather than thickness
        #We scale the mb_rate later
        H_vdiv = dem1_corr*vdiv_f
        #Compute mb - height change due to basal mb
        mb = diff_lag_fwd_obs + H_vdiv - ms
        if vel_input:
            mb_f = mb
        else:
            mb_f = malib.robust_spread_fltr(mb)
        #Inverting here, so positive rate means melting of bottom
        mb_rate = glaclib.freeboard_thickness(-mb_f,clip=False)/t_factor_lag
        mb_stats = malib.print_stats(mb_rate)
        dst_fn = os.path.join(outdir, outprefix+'_meltrate.tif')
        iolib.writeGTiff(mb_rate.astype(float), dst_fn, dem1_ds, ndv=diffndv)
