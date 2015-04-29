# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from os import path

def rangecolor(density):
    if density<100:
        return 'r'
    if density<200:
        return 'b'
    if density<300:
        return 'g'
    if density<400:
        return 'c'
    return 'm'

def normalize(density, rhomax=500):
    if density>rhomax:
        return 1
    return density/rhomax

plt.close('all')

dtformat_default = '%d.%m. %H:%M'
dtformat_snex = '%Y %d %B %H UTC'

e = EventsCollection('cases/pip2015.csv', dtformat_snex)
e.autoimport_data(autoshift=False, autobias=False, rule='6min', varinterval=True)

ax = plt.gca()
cmap = mpl.cm.gnuplot
norm = mpl.colors.Normalize(vmin=0,vmax=500)
cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm)

for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
    c.pluvio.shift_periods = -6
    plt.figure()
    colorstr = c.density().apply(rangecolor)
    colorval = c.density().apply(normalize)
    colorstr.name = 'colorstr'
    colorval.name = 'colorval'
    read.merge_multiseries(c.pipv.fits, colorstr, colorval).apply(lambda row: row.polfit.plot(color=cmap(row.colorval), linewidth=1, xmax=10), axis=1)
    plt.title(str(c.dt_start_end()[0].date()))
    plt.ylabel='fall velocity'
    plt.xlabel='diameter'