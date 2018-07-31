#! /usr/bin/env python

#Interpolate stack dhdt to specified timestamps

import os
import sys
import numpy as np
from datetime import datetime
from pygeotools.lib import timelib, iolib

#SRTM, then systematic timestamps
dt_list = [datetime(2000,2,11), datetime(2000,5,31), datetime(2009,5,31), datetime(2018,5,31)]
stack_fn=sys.argv[1]

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

for dt in dt_list:
    dt_o = timelib.dt2o(dt)
    z = trend*dt_o+intercept
    out_fn=os.path.splitext(stack_fn)[0]+'_%s.tif' % dt.strftime('%Y%m%d')
    print("Writing out: %s" % out_fn)
    iolib.writeGTiff(z, out_fn, trend_ds)
