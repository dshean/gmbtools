#! /usr/bin/env python

"""
Download SRTM-GL1 tiles for given lat/lon bounds
Pipe this to a file, then wget
"""

import os
import sys
import glob
import urllib2

file_list = glob.glob("*hgt.zip")

lon = (66, 106)
lat = (25, 47)

tile_list = []
for i in range(lon[0], lon[1]):
    for j in range(lat[0], lat[1]):
        tile_list.append('N%02iE%03i' % (j, i))

for tile in tile_list:
    print("https://cloud.sdsc.edu/v1/AUTH_opentopography/Raster/SRTM_GL1/SRTM_GL1_srtm/%s.hgt" % tile)
    #print("http://e4ftl01.cr.usgs.gov/SRTM/SRTMGL1.003/2000.02.11/%s.SRTMGL1.hgt.zip" % tile)

sys.exit()

response = urllib2.urlopen("http://e4ftl01.cr.usgs.gov/SRTM/SRTMGL1.003/2000.02.11/").readlines()
for line in response:
# <img src="/icons/compressed.gif" alt="[   ]"> <a href="N00E017.SRTMGL1.hgt.zip">N00E017.SRTMGL1.hgt.zip</a>     09-Oct-2014 11:06  6.0M
    if "/icons/compressed.gif" in line:
        if line.split('"')[5] not in file_list and 'hgt.zip' in line.split('"')[5]:
            os.system("wget -c http://e4ftl01.cr.usgs.gov/SRTM/SRTMGL1.003/2000.02.11/%s" % line.split('"')[5])
