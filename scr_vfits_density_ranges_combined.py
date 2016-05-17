# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
import snowfall as sf
import numpy as np
import read
import matplotlib.pyplot as plt
from os import path
import fit
import gc
import pandas as pd

from scr_snowfall import pip2015events, test_events, rholimits, param_table

debug = False
unfiltered = False
tld = '.eps'
bnds = (0.35,3.0,0.5,1.8)

def vel_fit(data, rhomin=None):
    if rhomin is not None:
        data = data.query('rhomin=={0}'.format(rhomin))
    d = np.concatenate([vfit.x_unfiltered for vfit in data.polfit])
    v = np.concatenate([vfit.y_unfiltered for vfit in data.polfit])
    return d, v

def vfits_debug(comb):
    fitargs = {'force_flip': False,
               'try_flip': False,
               'fitclass': fit.PolFit}
    rholimits=(0, 100, 200, 800)
    limslist = sf.limitslist(rholimits)
    return comb.vfits_density_range(limslist, **fitargs)

if debug:
    e = test_events()
else:
    e = pip2015events()

#plt.close('all')
plt.ioff()

data = param_table(e=e, debug=debug, rho_limits=rholimits)
data.index = data.index.droplevel()
vdata = read.merge_series(e.pluv_grouper(), e.vel_data())
table = pd.merge(vdata, pd.DataFrame(data.rhomin), left_on='group', right_index=True)

comb = e.events.paper.sum()
del(e)
gc.collect()

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
#merger = comb.d_0()
#data = comb.instr['pipv'].good_data()
#data_grouped = comb.group(data, merger)
#data_fltr = data_grouped[data_grouped['D_0'] > 0.63]
#fitargs['data'] = data_fltr
#fig_fltr, axarr_fltr = comb.plot_vfits_in_density_ranges(fitargs=fitargs, **kws)
axarr[0].axis(bnds)
#axarr_fltr[0].axis(bnds)

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
#fig_fltr.savefig(path.join(savepath, fname + '_d0fltr' + tld))
fig.savefig(path.join(paperpath, 'vfits_rho_ranges' + tld))
