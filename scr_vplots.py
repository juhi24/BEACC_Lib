# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt
from os import path
import read

dtformat_default = '%d.%m. %H:%M'
dtformat_snex = '%Y %d %B %H UTC'

e = EventsCollection('cases/pip2015test.csv', dtformat_snex)
e.autoimport_data(autoshift=False, autobias=False, rule='6min', varinterval=True)

plt.close('all')
plt.ioff()

basepath = '../results/pip2015/paper/vfit'
fname = '%Y%m%d_%H%M.eps'

for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
    c.instr['pluvio'].shift_periods = -6
    c.instr['pluvio'].n_combined_intervals = 2
    c.density() # to initialize fits
    savepath = read.ensure_dir(path.join(basepath, c.dtstr('%Y%m%d'), c.instr['pluvio'].name))
    fits = c.instr['pipv'].fits.polfit
    axlist = []
    flist = []
    for i, vfit in fits.iteritems():
        flist.append(plt.figure(dpi=120, figsize=(4,3)))
        ax = vfit.plot(source_style='hex', unfiltered=True)
        axlist.append(ax)
        plt.axis([None, 5, 0.5, 2.5])
        plt.legend(loc='lower right')
        plt.xlabel('Equivalent diameter (mm)')
        plt.ylabel('Fall velocity (ms$^{-1}$)')
        plt.tight_layout()
        plt.savefig(os.path.join(savepath, i.strftime(fname)))
