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

debug = False

if debug:
    e = test_events()
else:
    e = pip2015events()

#plt.close('all')
plt.ion()

comb = e.events.paper.sum()
del(e)

fig, axarr = comb.plot_vfits_in_density_ranges(separate=True,
                                          rholimits=(0, 150, 300, 800),
                                          source_style='kde',
                                          fitargs={'force_flip': False,
                                                   'try_flip': False,
                                                   'fitclass': fit.PolFit},
                                          unfiltered=True, parallel=False)
plt.axis((0.25,2.8,0.5,1.8))

savepath = '../results/pip2015/vfits_density_ranges'
if debug:
    savepath += '/test'
read.ensure_dir(savepath)
plt.savefig(path.join(savepath, 'combined.eps'))