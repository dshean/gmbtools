#! /usr/bin/env python
"""
Utility to fetch SNOTEL data

Uses ulmo library and CUAHSI db
"""

import sys
from datetime import datetime

import pytz
import numpy as np
import matplotlib.pyplot as plt

import ulmo
from ulmo.util import convert_datetime

from pygeotools.lib import malib

#Get all site locations
def get_latlon():
    sites = ulmo.cuahsi.wof.get_sites(wsdlurl)
    lon = []
    lat = []
    code = []
    z = []
    for k,v in sites.iteritems():
        lon.append(float(v['location']['longitude']))
        lat.append(float(v['location']['latitude']))
        code.append(int(v['code']))
        z.append(float(v['elevation_m']))

    csv_fn = 'snotel_latlon.csv'
    out = np.array(zip(code,lat,lon))
    np.savetxt(csv_fn, out, delimiter=',', fmt='%i,%0.5f,%0.5f')
    return out

#Get datetime series
def get_series_dt(series, strptime_fmt='%Y-%m-%dT%H:%M:%S'):
    #ts = [convert_datetime(vd['date_time_utc']).replace(tzinfo=pytz.utc) for vd in series['values']]
    ts = [datetime.strptime(vd['date_time_utc'], strptime_fmt) for vd in series['values']]
    return np.array(ts, dtype=np.datetime64)

#Get value series
def get_series_val(series):
    # Create a clean timeseries list of (dt,val) tuples
    val = [float(vd['value']) for vd in series['values']]
    val = np.ma.masked_equal(val, -9999)
    val = np.ma.masked_equal(val, 0.0)
    return val

#URL for query
wsdlurl = "http://worldwater.byu.edu/interactive/snotel/services/index.php/cuahsi_1_1.asmx?WSDL" 

#Set start/end date range for data
dt_start = datetime(2007,10,1)
dt_end = datetime.now()

sitename = 'baker'
#Baker site codes
#Want to use output from swe.py, which searches for a given DEM extent
sitecodelist = [999, 909, 1011, 910]
vlist = ['SNWD', 'WTEQ']
#Incremental precip, cumulative precip
#SNWD is cm
#PRCP (mm), PREC (mm), SNWD (cm), TAVG, TMAX, TMIN, WTEQ (mm)
#v = 'SNOTEL:SNWD'

#DEM stack, can be used to plot lines on SNOTEL time series
stack_fn = '/Volumes/SHEAN_1TB_SSD/site_poly_highcount_rect3/baker/swe/20130911_1938_1030010027BE9000_1030010026900000-DEM_8m_trans_20160604_1941_104001001D940300_104001001CB1B100-DEM_8m_trans_stack_22.npz'
stack = malib.DEMStack(stack_fn=stack_fn)
dem_dt = stack.date_list

#Accuracy of measurements, in cm
#https://www.wcc.nrcs.usda.gov/snotel/snotel_sensors.html
sigma_factor = 3
snwd_precision = sigma_factor*1.27/100.
wteq_precision = sigma_factor*0.254/100.

d = {}
#Should split this up with multiprocessing rather than running serially - one job per site
for sitecode in sitecodelist:
    print('Processing site: %i' % sitecode)

    sitekey = 'SNOTEL:%i' % sitecode
    #site = ulmo.cuahsi.wof.get_site_info(wsdlurl, sitekey)

    #Get first variable, use to set dates
    v = vlist[0]
    sitev = 'SNOTEL:%s' % v
    print(sitev)
    series = ulmo.cuahsi.wof.get_values(wsdlurl, sitekey, sitev, start=dt_start, end=dt_end)
    dt = get_series_dt(series)
    d[sitecode] = {'dt':dt}
    val = get_series_val(series)
    d[sitecode][v] = val

    for v in vlist[1:]:
        sitev = 'SNOTEL:%s' % v
        print(sitev)
        series = ulmo.cuahsi.wof.get_values(wsdlurl, sitekey, sitev, start=dt_start, end=dt_end)
        #dt = series['values']['date_time_utc']
        #vals = series['values']['value']
        val = get_series_val(series)
        #Looks like these are not always updated simultaneously, make sure the records are same length
        #Should probably just query both dt and vals simultaneously, rather than assume all variables are same length
        if val.size != dt.size:
            val = val[0:dt.size]
        d[sitecode][v] = val

    #Convert SNWD to m 
    d[sitecode]['SNWD'] /= 100.
    d[sitecode]['WTEQ'] /= 1000.

    #Mask values less than instrument precision
    d[sitecode]['SNWD'] = np.ma.masked_less(d[sitecode]['SNWD'], snwd_precision)
    d[sitecode]['WTEQ'] = np.ma.masked_less(d[sitecode]['WTEQ'], wteq_precision)

    #Calculate density in g/cc
    rho = (d[sitecode]['WTEQ']/d[sitecode]['SNWD'])
    #Mask values when snow depth is small, helps avoid bogus density values
    depth_thresh = 0.2
    rho[(d[sitecode]['SNWD'] < depth_thresh)] = np.ma.masked
    d[sitecode]['Density'] = rho

print("Plotting")
vlist.append('Density')
f, axa = plt.subplots(len(vlist), 1, sharex=True, figsize=(10,7.5))
for sitecode in d.keys():
    dt = d[sitecode]['dt']
    for n,v in enumerate(vlist):
        vmed = np.ma.median(d[sitecode][v])
        #vmean = np.ma.mean(d[sitecode][v[)
        #lbl = '%s: %0.2f' % (sitecode, vmed)
        lbl = str(sitecode)
        p = axa[n].plot(dt, d[sitecode][v], marker='o', ms=1, linestyle='', label=lbl)
        axa[n].set_ylabel(vlist[n])
    axa[n].axhline(vmed, c=p[0].get_color(), linestyle=':', linewidth=0.5)

axa[n].xaxis_date()
axa[n].set_ylim(0,1.0)
#axa[0].set_title('SNOTEL')
axa[n].legend(prop={'size':8})

axa[0].set_ylabel('Snow Depth (m)')
axa[1].set_ylabel('SWE (m w.e.)')
axa[2].set_ylabel('Density (g/cc)')

#Plot lines for DEM timestamps
if stack_fn is not None:
    for dt in dem_dt:
        axa[0].axvline(dt, color='k', alpha=0.2)

plt.tight_layout()
f.autofmt_xdate()
fig_fn = '%s_SNOTEL.pdf' % sitename
plt.savefig(fig_fn, bbox_inches='tight')
fig_fn = '%s_SNOTEL.png' % sitename
plt.savefig(fig_fn, dpi=300, bbox_inches='tight')

axa[n].set_xlim(datetime(2013,8,1), datetime(2016,6,30))
fig_fn = '%s_SNOTEL_2013-2016.png' % sitename
f.set_size_inches(4,7.5)
plt.tight_layout()
plt.savefig(fig_fn, dpi=300, bbox_inches='tight')

plt.show()

"""
Alternative approach using climata lib - results in SSL Certificate errors
from climata.snotel import RegionDailyDataIO
#basin="18010202",
data = RegionDailyDataIO(
    start_date="2014-01-01",
    end_date="2014-07-01",
    parameter="WTEQ",
)

for series in data:
    print series
    for row in series.data:
        print row
"""

