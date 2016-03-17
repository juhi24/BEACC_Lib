# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
import numpy as np
import read
import matplotlib.pyplot as plt
from os import path
import fit

from scr_snowfall import pip2015events, test_events

debug = False
unfiltered = False
tld = '.eps'
rholimits = (0, 100, 200, 800)
#rholimits = (0, 150, 300, 800)
bnds = (0.25,2.8,0.5,1.8)


if debug:
    e = test_events()
else:
    e = pip2015events()

#plt.close('all')
plt.ioff()

comb = e.events.paper.sum()
del(e)

hextent = np.array(bnds)+(-0.1,0.1,-0.1,0.1)
kws = {'separate': True,
       'rholimits': rholimits,
       'source_style': 'hex',
       'source_kws': {'gridsize': 26, 'extent': hextent},
       'unfiltered': True}
fitargs = {'force_flip': False,
           'try_flip': False,
           'fitclass': fit.PolFit}
if unfiltered:
    fitargs['filter_outliers'] = False
fig, axarr = comb.plot_vfits_in_density_ranges(fitargs=fitargs, **kws)
merger = comb.d_0()
data = comb.instr['pipv'].good_data()
data_grouped = comb.group(data, merger)
data_fltr = data_grouped[data_grouped['D_0'] > 0.63]
fitargs['data'] = data_fltr
fig_fltr, axarr_fltr = comb.plot_vfits_in_density_ranges(fitargs=fitargs, **kws)
axarr[0].axis(bnds)
axarr_fltr[0].axis(bnds)

resultsdir = '../results/pip2015'
savepath = path.join(resultsdir, 'vfits_density_ranges')
paperpath = read.ensure_dir(path.join(resultsdir, 'paper'))
if debug:
    savepath += '/test'
read.ensure_dir(savepath)
fname = 'combined'
if unfiltered:
    fname += '_unfiltered'
fig.savefig(path.join(savepath, fname + tld))
fig_fltr.savefig(path.join(savepath, fname + '_d0fltr' + tld))
fig_fltr.savefig(path.join(paperpath, 'vfits_rho_ranges' + tld))
