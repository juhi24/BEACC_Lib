# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import read
import numpy as np
import matplotlib.pyplot as plt
from os import path
import fit

from scr_snowfall import pip2015events, test_events

debug = True
unfiltered = False
tld = '.eps'

if debug:
    e = test_events()
else:
    e = pip2015events()

#plt.close('all')
plt.ion()

comb = e.events.paper.sum()
del(e)

kws = {}
fitargs = {'force_flip': False,
           'try_flip': False,
           'fitclass': fit.PolFit,
           'frac': 0.1}
if unfiltered:
    fitargs['filter_outliers'] = False
fig, axarr = comb.plot_vfits_in_density_ranges(separate=True,
                                               rholimits=(0, 150, 300, 800),
                                               source_style='kde',
                                               fitargs=fitargs,
                                               unfiltered=True,
                                               **kws)
plt.axis((0.25,2.8,0.5,1.8))

savepath = '../results/pip2015/vfits_density_ranges'
if debug:
    savepath += '/test'
read.ensure_dir(savepath)
fname = 'combined'
if unfiltered:
    fname += '_unfiltered'
plt.savefig(path.join(savepath, fname + tld))
