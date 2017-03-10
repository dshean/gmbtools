#! /usr/bin/env python

from datetime import datetime

import pytz
import numpy as np
import matplotlib.pyplot as plt

import ulmo
from ulmo.util import convert_datetime

"""
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

#def site_series_values_to_df(series_values, variable_name):
def get_series_dt(series, strptime_fmt='%Y-%m-%dT%H:%M:%S'):
    #ts = [convert_datetime(vd['date_time_utc']).replace(tzinfo=pytz.utc) for vd in series['values']]
    ts = [datetime.strptime(vd['date_time_utc'], strptime_fmt) for vd in series['values']]
    return np.array(ts, dtype=np.datetime64)

def get_series_val(series):
    # Create a clean timeseries list of (dt,val) tuples
    val = [float(vd['value']) for vd in series['values']]
    val = np.ma.masked_equal(val, -9999)
    val = np.ma.masked_equal(val, 0.0)
    return val

wsdlurl = "http://worldwater.byu.edu/interactive/snotel/services/index.php/cuahsi_1_1.asmx?WSDL" 

#Import DEM
#Extract extent, timestamp
#Find all sites within buffered distance

dt_start = datetime(2007,10,1)
dt_end = datetime.now()

#Baker codes
sitename = 'baker'
sitecodelist = [999, 909, 1011, 910]
vlist = ['SNWD', 'WTEQ']
#Incremental precip, cumulative precip
#SNWD is cm
#PRCP (mm), PREC (mm), SNWD (cm), TAVG, TMAX, TMIN, WTEQ (mm)
#v = 'SNOTEL:SNWD'

#Accuracy of measurements, in cm
#https://www.wcc.nrcs.usda.gov/snotel/snotel_sensors.html
snwd_precision = 3*1.27/100.
wteq_precision = 3*0.254/100.

d = {}
#Should split this up with multiprocessing - one job per site
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
        d[sitecode][v] = val

    #Convert SNWD to m 
    d[sitecode]['SNWD'] /= 100.
    d[sitecode]['SNWD'] = np.ma.masked_less(d[sitecode]['SNWD'], snwd_precision)
    d[sitecode]['WTEQ'] /= 1000.
    d[sitecode]['WTEQ'] = np.ma.masked_less(d[sitecode]['WTEQ'], wteq_precision)

    #Calculate density in g/cc
    rho = (d[sitecode]['WTEQ']/d[sitecode]['SNWD'])
    d[sitecode]['Density'] = rho

print("Plotting")
vlist.append('Density')
f, axa = plt.subplots(len(vlist), 1, sharex=True, figsize=(10,7.5))
for sitecode in d.keys():
    dt = d[sitecode]['dt']
    for n,v in enumerate(vlist):
        vmed = np.ma.median(d[sitecode][v])
        #vmean = np.ma.mean(d[sitecode][v[)
        p = axa[n].plot(dt, d[sitecode][v], marker='o', ms=1, linestyle='', label='%s: %0.2f' % (sitecode, vmed))
        axa[n].set_ylabel(vlist[n])
    axa[n].axhline(vmed, c=p[0].get_color(), linestyle=':', linewidth=0.5)

axa[n].set_ylim(0,1.0)
#axa[0].set_title('SNOTEL')
axa[n].legend(prop={'size':8})
plt.tight_layout()
fig_fn = '%s_SNOTEL.pdf' % sitename
plt.savefig(fig_fn, bbox_inches='tight')

plt.show()

