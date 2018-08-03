#! /usr/bin/env python

"""
Generate bar plot for annual DEM (or any feature) count from input shapefile
Expects features with date field with format YYYYMMDD
"""

import os
import sys

import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt

#fn = '/Users/dshean/Documents/UW/HMA/mos/hma_20180609_mos/hma_mos/hma_mos_32m_stripindex_simp.shp'
fn = sys.argv[1]
df = gpd.read_file(fn)

#Can also generate from a simple list of dates
#ogrinfo -al *stripindex_simp.shp | grep date | awk '{print $NF}' > hma_dem_dates.txt
#fn = 'hma_dem_dates.txt'
#df = pd.read_csv(fn)

df['dt'] = pd.to_datetime(df['date'], format='%Y%m%d')

#df['sensor'] = df['name'].str.split('/').str[1].str[0:4]
#sensors = df['sensor'].unique()
#df.groupby('sensor').count()

#https://stackoverflow.com/questions/27365467/can-pandas-plot-a-histogram-of-dates
#ax = df['dt'].groupby([df["dt"].dt.year, df["dt"].dt.month]).count().plot(kind="bar", color='k')
ax = df['dt'].groupby(df["dt"].dt.year).count().plot(kind="bar", color='k')
ax.set_ylabel('DEM Count')
ax.set_xlabel('')
#ax.set_title('Annual DEM count')

out_fn = os.path.splitext(fn)[0] + '_annual_count.pdf'
plt.savefig(out_fn)
plt.show()

#Plot individual sensors 
#sgb = df.groupby(['sensor', df['dt'].dt.year]).count()

#f, ax = plt.subplots(1, figsize=(10,7.5))
#s = 'WV01'
