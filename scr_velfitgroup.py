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

plt.ioff()
plt.close('all')

dtformat_default = '%d.%m. %H:%M'
dtformat_snex = '%Y %d %B %H UTC'

e = EventsCollection('cases/pip2015.csv', dtformat_snex)
e.autoimport_data(autoshift=False, autobias=False, rule='6min', varinterval=True)

ax = plt.gca()
cmap = mpl.cm.gnuplot
norm = mpl.colors.Normalize(vmin=0,vmax=500)
#cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm)

for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
    c.pluvio.shift_periods = -6
    #plt.figure()
    colorstr = c.density().apply(rangecolor)
    colorval = c.density().apply(normalize)
    colorstr.name = 'colorstr'
    colorval.name = 'colorval'
    merged = read.merge_multiseries(c.pipv.fits, colorstr, colorval, c.density())
    #merged.apply(lambda row: row.polfit.plot(color=row.colorstr, linewidth=1, xmax=10), axis=1)
    groups = merged.groupby(colorstr)
    for name, group in groups:
        plt.figure()
        group.apply(lambda row: row.polfit.plot(linewidth=1, xmax=10, color=name, ), axis=1)
        rho = group.density.mean()
        if rho < 100:
            rhorange = 'rho<100'
        elif rho < 200:
            rhorange = '100<rho<200'
        elif rho < 300:
            rhorange = '200<rho<300'
        elif rho < 400:
            rhorange = '300<rho<400'
        else:
            rhorange = 'rho>400'
        dtstr = str(c.dt_start_end()[0].date())
        plt.ylim((0,3))
        plt.title('%s, %s (%s)' % (dtstr, rhorange, c.pluvio.name))
        plt.ylabel('fall velocity')
        plt.xlabel('diameter')
        plt.savefig('../results/pip2015/velfitgroups/' + c.pluvio.name + '_' + dtstr + '_' + rhorange + '.eps')