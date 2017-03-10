#! /usr/bin/env python

import numpy as np

import ulmo

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

wsdlurl = "http://worldwater.byu.edu/interactive/snotel/services/index.php/cuahsi_1_1.asmx?WSDL" 
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
np.savetxt(csv_fn, zip(code,lat,lon), delimiter=',', fmt='%i,%0.5f,%0.5f')

#plt.figure()
#plt.scatter(lon,lat)
