# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from os import path

def rangecolor(rho, rholimits=(150, 300), colors=('r', 'b', 'c', 'g', 'm')):
    for i, rhomax in enumerate(rholimits):
        if rho < rhomax:
            return colors[i]
    return colors[i+1]

def rho_range_str(rho, rholimits=(150, 300)):
    rhomax_old = ''
    for rhomax in rholimits:
        rhomin = rhomax_old
        rhomax_old = str(rhomax) + '<'
        if rho < rhomax:
            return '%srho<%s' % (rhomin, rhomax)
    return 'rho>%s' % rhomax

def segment_index(val, limits=(150, 300)):
    for i, maxval in enumerate(limits):
        if val < maxval:
            return i
    return i+1

def normalize(val, maxval=500):
    if val>maxval:
        return 1
    return val/maxval

plt.ioff()
plt.close('all')

dtformat_default = '%d.%m. %H:%M'
dtformat_snex = '%Y %d %B %H UTC'

e = EventsCollection('cases/pip2015.csv', dtformat_snex)
e.autoimport_data(autoshift=False, autobias=False, rule='6min', varinterval=True)

#ax = plt.gca()
#cmap = mpl.cm.gnuplot
#norm = mpl.colors.Normalize(vmin=0,vmax=500)
#cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm)

rholimits = (150, 300)
combined = False
subplots = True
outpath = '../results/pip2015/velfitgroups2/'
extra = ''
subpath = ''
if subplots:
    subpath = 'subplots/'
    combined = False
elif combined:
    subpath = 'combined/'

for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
    c.instr['pluvio'].shift_periods = -6
    if combined:
        fig = plt.figure()
    colorstr = c.density().apply(rangecolor, rholimits=rholimits)
    colorval = c.density().apply(normalize)
    colorstr.name = 'colorstr'
    colorval.name = 'colorval'
    merged = read.merge_multiseries(c.pipv.fits, colorstr, colorval, c.density())
    #merged.apply(lambda row: row.polfit.plot(color=row.colorstr, linewidth=1, xmax=10), axis=1)
    groups = merged.groupby(colorstr)
    if subplots:
        fig, axarr = plt.subplots(ncols=len(rholimits)+1, figsize=(22,6),
                                  sharex='all', sharey='all')
    for name, group in groups:
        rho = group.density.mean()
        alpha = 1
        i = segment_index(rho, limits=rholimits)
        rhorange = rho_range_str(rho, rholimits=rholimits)
        if not subplots:
            if combined:
                ax = plt.gca()
                alpha = 0.2
            else:
                fig, ax = plt.subplots()
                extra = '_' + rhorange
        else:
            ax = axarr[i]
        group.apply(lambda row: row.polfit.plot(ax=ax, linewidth=0.5, color=name, alpha=alpha), axis=1)
        ax.set_xlim([0, 10])
        dtstr = str(c.dt_start_end()[0].date())
        ax.set_ylim((0,3))
        ax.set_title('%s, %s (%s)' % (dtstr, rhorange, c.instr['pluvio'].name))
        if not (subplots and i>0):
            ax.set_ylabel('fall velocity')
        ax.set_xlabel('equivalent diameter')
        ax.grid(axis='y')
        if subplots:
            fig.subplots_adjust(hspace=0)
            plt.setp([a.get_yticklabels() for a in fig.axes[1:]], visible=False)
        elif not combined:
            plt.savefig(read.ensure_dir(outpath + subpath) + c.instr['pluvio'].name + '_' + dtstr + extra + '.eps')
    if subplots:
            fig.subplots_adjust(wspace=0)
            plt.setp([a.get_yticklabels() for a in fig.axes[1:]], visible=False)
    plt.savefig(read.ensure_dir(outpath + subpath) + c.instr['pluvio'].name + '_' + dtstr + extra + '.eps')