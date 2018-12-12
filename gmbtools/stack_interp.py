#! /usr/bin/env python

#Interpolate stack dhdt to specified timestamps

import os
import sys
import numpy as np
from datetime import datetime
from pygeotools.lib import timelib, iolib, filtlib

#These are from TanDEM-X masked
#Not sure why negative outliers
#zlim = (-1013,8388)
#Ellispoidal heights at sea level should be above -100 m
zlim = (-200,8388)

#Input stack npz
stack_fn=sys.argv[1]

#User-provided timestamps with format YYYYMMDD
if len(sys.argv) > 2:
    dt_list_str = sys.argv[2:]
    if len(dt_list_str) == 1:
        dt_list_str = dt_list_str[0].split(' ')
    dt_list = [timelib.fn_getdatetime_list(dt_str)[0] for dt_str in dt_list_str]
else:
    #SRTM, then systematic ASTER timestamps
    dt_list = [datetime(2000,2,11), datetime(2000,5,31), datetime(2009,5,31), datetime(2018,5,31)]

#Use tif on disk if available
trend_fn=os.path.splitext(stack_fn)[0]+'_trend.tif'
intercept_fn=os.path.splitext(stack_fn)[0]+'_intercept.tif'
#Otherwise load stack and compute trend/intercept if necessary

trend_ds = iolib.fn_getds(trend_fn)
trend = iolib.ds_getma(trend_ds)/365.25
intercept = iolib.fn_getma(intercept_fn)

#Can vectorize
#dt_list_o = timelib.dt2o(dt_list)
#z_list = trend*dt_list_o[:,None,None]+intercept

filter=True

if filter:
    #Could remove outliers in trend at this phase
    print(trend.count())
    trend_filt = filtlib.mad_fltr(trend, n=3.5)
    trend_filt = filtlib.rolling_fltr(trend_filt, size=5, circular=False)
    trend_filt = filtlib.gauss_fltr_astropy(trend_filt, size=5, fill_interior=True)
    out_fn=os.path.splitext(trend_fn)[0]+'_filt.tif'
    print("Writing out: %s" % out_fn)
    iolib.writeGTiff(trend_filt*365.25, out_fn, trend_ds)

    #Update intercept using new filtered slope values
    dt_pivot = timelib.dt2o(datetime(2009, 5, 31))
    intercept_filt = dt_pivot * (trend - trend_filt) + intercept
    out_fn=os.path.splitext(intercept_fn)[0]+'_filt.tif'
    print("Writing out: %s" % out_fn)
    iolib.writeGTiff(intercept_filt*365.25, out_fn, trend_ds)

    trend = trend_filt
    intercept = intercept_filt

for dt in dt_list:
    dt_o = timelib.dt2o(dt)
    z = trend*dt_o+intercept
    #Remove any values outsize global limits
    #Could also do local range filter here
    z = filtlib.range_fltr(z, zlim)
    print(z.count())
    out_fn=os.path.splitext(stack_fn)[0]+'_%s.tif' % dt.strftime('%Y%m%d')
    print("Writing out: %s" % out_fn)
    iolib.writeGTiff(z, out_fn, trend_ds)
