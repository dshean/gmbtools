#! /usr/bin/env python

import sys
import os
import glob
import numpy as np
from pygeotools.lib import malib, iolib

topdir = '/nobackup/deshean/conus/prism'

#Generate products for June-Sept, and Oct-May periods
monthly_dir = '/nobackupp8/deshean/conus/prism/normals/monthly'
monthly_ppt_fn = np.array(glob.glob(os.path.join(monthly_dir,'ppt/*bil.bil')))
monthly_ppt_fn.sort()
monthly_tmean_fn = np.array(glob.glob(os.path.join(monthly_dir,'tmean/*bil.bil')))
monthly_tmean_fn.sort()

summer_idx = np.array([5, 6, 7, 8])
winter_idx = np.array([9, 10, 11, 0, 1, 2, 3, 4])

monthly_ppt_summer_fn = os.path.join(monthly_dir, 'PRISM_ppt_30yr_normal_800mM2_06-09_summer.npz')
ppt_summer = malib.DEMStack(monthly_ppt_fn[summer_idx], stack_fn=monthly_ppt_summer_fn, sort=False, med=True, datestack=False, trend=False)
ppt_summer_cum = ppt_summer.ma_stack.sum(axis=0)
iolib.writeGTiff(ppt_summer_cum, os.path.splitext(monthly_ppt_summer_fn)[0]+'_cum.tif', ppt_summer.get_ds())

monthly_ppt_winter_fn = os.path.join(monthly_dir, 'PRISM_ppt_30yr_normal_800mM2_10-05_winter.npz')
ppt_winter = malib.DEMStack(monthly_ppt_fn[winter_idx], stack_fn=monthly_ppt_winter_fn, sort=False, med=True, datestack=False, trend=False)
ppt_winter_cum = ppt_winter.ma_stack.sum(axis=0)
iolib.writeGTiff(ppt_winter_cum, os.path.splitext(monthly_ppt_winter_fn)[0]+'_cum.tif', ppt_winter.get_ds())

monthly_tmean_summer_fn = os.path.join(monthly_dir, 'PRISM_tmean_30yr_normal_800mM2_06-09_summer.npz')
tmean_summer = malib.DEMStack(monthly_tmean_fn[summer_idx], stack_fn=monthly_tmean_summer_fn, sort=False, med=True, datestack=False, trend=False)

monthly_tmean_winter_fn = os.path.join(monthly_dir, 'PRISM_tmean_30yr_normal_800mM2_10-05_winter.npz')
tmean_winter = malib.DEMStack(monthly_tmean_fn[winter_idx], stack_fn=monthly_tmean_winter_fn, sort=False, med=True, datestack=False, trend=False)

#Use 800 m normals to scale 4 km resolution to each glacier
#Compute PDD
#Function to adjust DEM surface
