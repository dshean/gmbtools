#! /usr/bin/env python

"""
Download SRTM tiles for given lat/lon bounds
Pipe this to a file, then wget
"""

#HMA
#lon = (66, 106)
#lat = (25, 47)
#ew = "E"
#CONUS
lon = (105, 124)
lat = (36, 49)
ew = "W"

tile_list = []
for i in range(lon[0], lon[1]):
    for j in range(lat[0], lat[1]):
        #tile_list.append('N%02i%s%03i' % (j, ew, i))
        tile_list.append('N%02i%s%03i' % (j, ew, i))

#This is SRTM-GL1
#topurl = "https://cloud.sdsc.edu/v1/AUTH_opentopography/Raster/SRTM_GL1/SRTM_GL1_srtm/"
#topurl = "http://e4ftl01.cr.usgs.gov/SRTM/SRTMGL1.003/2000.02.11/" 

#This is NASADEM provisional
topurl = "https://e4ftl01.cr.usgs.gov/provisional/MEaSUREs/NASADEM/NorthAmerica/"
#topurl = "https://e4ftl01.cr.usgs.gov/provisional/MEaSUREs/NASADEM/Eurasia/"

#This is non void-filled float
product = 'hgt_srtmOnly_R4'
#n00e072.hgt.zip
#product = 'img_8bit'
#n00e072_040_290_ss2_a_1_1.img.zip

topurl += product

for tile in tile_list:
    #print("%s/%s.hgt" % (topurl, tile))
    #print("%s/%s.SRTMGL1.hgt.zip" % (topurl, tile))
    #For Eurasia
    #print("%s/%s.hgt.zip" % (topurl, tile.lower()))
    #For NorthAmerica
    print("%s/%s.srtmOnly.hgt.zip" % (topurl, tile.lower()))

"""
import os
import glob
import urllib2
file_list = glob.glob("*hgt.zip")

response = urllib2.urlopen("http://e4ftl01.cr.usgs.gov/SRTM/SRTMGL1.003/2000.02.11/").readlines()
for line in response:
# <img src="/icons/compressed.gif" alt="[   ]"> <a href="N00E017.SRTMGL1.hgt.zip">N00E017.SRTMGL1.hgt.zip</a>     09-Oct-2014 11:06  6.0M
    if "/icons/compressed.gif" in line:
        if line.split('"')[5] not in file_list and 'hgt.zip' in line.split('"')[5]:
            os.system("wget -c http://e4ftl01.cr.usgs.gov/SRTM/SRTMGL1.003/2000.02.11/%s" % line.split('"')[5])
"""
