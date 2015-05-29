# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from os import path

def rangecolor(rho, rholimits=[150,300], colors=['r','b','g','c','m']):
    for i, rhomax in enumerate(rholimits):
        if rho < rhomax:
            return colors[i]
    return colors[i+1]

def rho_range_str(rho, rholimits=[150,300]):
    rhomax_old = ''
    for rhomax in rholimits:
        rhomin = rhomax_old
        rhomax_old = str(rhomax) + '<'
        if rho < rhomax:
            return '%srho<%s' % (rhomin, rhomax)
    return 'rho>%s' % rhomax

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
        rhorange = rho_range_str(group.density.mean())
        dtstr = str(c.dt_start_end()[0].date())
        plt.ylim((0,3))
        plt.title('%s, %s (%s)' % (dtstr, rhorange, c.pluvio.name))
        plt.ylabel('fall velocity')
        plt.xlabel('equivalent diameter')
        plt.savefig(read.ensure_dir('../results/pip2015/velfitgroups2/') + c.pluvio.name + '_' + dtstr + '_' + rhorange + '.eps')
