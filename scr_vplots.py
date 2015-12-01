# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import matplotlib.pyplot as plt
from os import path
import read

from scr_snowfall import pip2015events, test_events

debug = False

if debug:
    e = test_events()
else:
    e = pip2015events()

plt.close('all')
plt.ioff()

basepath = '../results/pip2015/paper/vfit'
if debug:
    basepath += '/test'
fname = '%Y%m%d_%H%M.eps'
date_format = '%d. %m. %Y'
time_format = '%H:%M'

extent = (0.375, 5, 0.5, 2.5)

for c in e.events.paper.values:
    c.density() # to initialize fits
    savepath = read.ensure_dir(path.join(basepath, c.dtstr('%Y%m%d')))
    fits = c.instr['pipv'].fits.polfit
    start_time = c.instr['pluvio'].start_time()
    axlist = []
    flist = []
    for i, vfit in fits.iteritems():
        flist.append(plt.figure(dpi=150, figsize=(4,3)))
        ax = vfit.plot(source_style='hex', unfiltered=True,
                       source_kws={'gridsize': 40, 'extent': extent})
        axlist.append(ax)
        ax.axis(extent)
        ax.legend(loc='lower right')
        ax.set_xlabel('Equivalent diameter (mm)')
        ax.set_ylabel('Fall velocity (ms$^{-1}$)')
        tstr_start = start_time[i].strftime(time_format)
        tstr_end = i.strftime(time_format)
        ax.set_title('%s - %s' % (tstr_start, tstr_end))
        plt.tight_layout()
        plt.savefig(os.path.join(savepath, i.strftime(fname)))
