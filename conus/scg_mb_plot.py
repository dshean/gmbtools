#! /usr/bin/env python

import matplotlib.pyplot as plt
from collections import OrderedDict

#The term U.S.Geological Survey "water year" in reports that deal with surface-water supply is defined as the 12-month period October 1, for any given year through September 30, of the following year. The water year is designated by the calendar year in which it ends and which includes 9 of the 12 months. Thus, the year ending September 30, 1999 is called the "1999" water year.

#Winter 2008/2009
3.26, 3.13
#Winter 2013/2014
3.64, 3.24
#Winter 2014/2015
1.07, 1.03
#Winter 2015/2016
1.56, 1.54

#Summer 2012 to Summer 2013
-1.81, -1.72
#Summer 2013 to Summer 2014
0.64, 0.63
#Summer 2014 to Summer 2015
-3.15, -3.12

mbd = OrderedDict() 
mbd[2009] = {'w':3.26}
mbd[2013] = {'a':-1.81}
mbd[2014] = {'w':3.64,'a':0.64}
mbd[2015] = {'w':1.07,'a':-3.15}
mbd[2016] = {'w':1.56,}

f, ax = plt.subplots()
bar_prop = {'linewidth':0.5, 'edgecolor':'k'}
wlbl = False
slbl = False
albl = False
for yr,v in mbd.iteritems():
    if 'w' in v.keys():
        if not wlbl:
            lbl = 'Winter'
            wlbl = True
        else:
            lbl = None
        ax.bar(yr, v['w'], color='b', label=lbl, **bar_prop)
    if 'a' in v.keys():
        if 'w' in v.keys():
            v['s'] = v['a'] - v['w']
            if not slbl:
                lbl = 'Summer'
                slbl = True
            else:
                lbl = None
            ax.bar(yr, v['s'], color='r', label=lbl, **bar_prop)
        if not albl:
            lbl = 'Annual'
            albl = True
        else:
            lbl = None
        ax.bar(yr, v['a'], color='gold', label=lbl, **bar_prop)
ax.axhline(0, color='k', lw=0.5)
ax.set_ylim(-5, 5)
ax.set_ylabel('Mass Balance (m w.e.)')
ax.set_title('Annual and Seasonal Geodetic Mass Balance')
ax.set_xlabel('Year')
ax.legend(prop={'size':10}, loc='lower left')
#ax.minorticks_on()
ax.tick_params(axis='y',which='minor',left='on')
plt.tight_layout()
fig_fn = 'scg_annual_seasonal_barchart.pdf'
plt.savefig(fig_fn, bbox_inches='tight'
