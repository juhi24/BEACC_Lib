# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt
from os import path
import seaborn as sns

sns.set_style('ticks')

dtformat_default = '%d.%m. %H:%M'
dtformat_snex = '%Y %d %B %H UTC'

e = EventsCollection('cases/pip2015.csv', dtformat_snex)
e.autoimport_data(autoshift=False, autobias=False, rule='6min', varinterval=True)

plt.close('all')
#plt.ioff()
plt.ion()
kwargs = {'kde': True, 'rug': True, 'kde_kws': {'label': 'KDE'}}

def subplots(n_plots=1):
    return plt.subplots(n_plots, sharex=True, sharey=True,
                        tight_layout=False, dpi=100)

def remove_gaps(f):
    f.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

def plots(data, axd, axm, axn, label=None, title=None, **kwtitle):
    rng = (0,6)
    sns.distplot(data.D_0.dropna(), ax=axd, label=label, bins=17, 
                 hist_kws={'range':rng}, **kwargs)
    axd.set_xlim(rng)
    axd.yaxis.set_ticks(np.arange(0.4, 2.0, 0.4))
    rng = (-2, 8)
    sns.distplot(data.mu.dropna(), ax=axm, label=label, bins=20,
                       hist_kws={'range':rng}, **kwargs)
    axm.set_xlim(rng)
    axm.yaxis.set_ticks(np.arange(0.2, 0.7, 0.2))
    sns.distplot(data.N_w.dropna(), ax=axn, label=label,
                       bins=10**np.linspace(0,6,20), kde=False, rug=True)
    axn.set_xscale('log')
    if title is not None:
        for ax in (axd, axm, axn):
            ax.set_title(title, **kwtitle)
        for ax in (axd, axm):
            ax.legend().set_visible(False)

n_cases = e.events.pluvio400.count()
fd, axarrd = subplots(n_cases)
fm, axarrm = subplots(n_cases)
fn, axarrn = subplots(n_cases)

for i, c in enumerate(e.events.pluvio400.values):
    c.instr['pluvio'].shift_periods = -6
    c.instr['pluvio'].n_combined_intervals = 2
    data = read.merge_multiseries(c.d_0(), c.mu(), c.n_w())
    plots(data, axarrd[i], axarrm[i], axarrn[i], title=c.dtstr(), y=0.85,
          fontdict={'verticalalignment': 'top', 'fontsize': 10})

comb400 = e.events.pluvio400.sum()
limslist = limitslist((0, 150, 300, 800))
n_ranges = len(limslist)
fdd, axarrdd = subplots(n_ranges)
fmd, axarrmd = subplots(n_ranges)
fnd, axarrnd = subplots(n_ranges)
data400 = read.merge_multiseries(comb400.d_0(), comb400.mu(), comb400.n_w())

for i, lims in enumerate(limslist):
    data = comb400.data_in_density_range(data400, *lims)
    limitstr = '$%s < \\rho < %s$' % (lims[0], lims[1])
    plots(data, axarrdd[i], axarrmd[i], axarrnd[i], title=limitstr, y=0.9, fontdict={'verticalalignment': 'top'})

for ax in (axarrdd[-1], axarrmd[-1], axarrnd[-1]):
    ax.set_title('$\\rho > %s$' % lims[0], y=0.9, fontdict={'verticalalignment': 'top'})

for f in (fd, fm, fn, fdd, fmd, fnd):
    remove_gaps(f)

for axarr in (axarrd, axarrm, axarrn, axarrdd, axarrmd, axarrnd):
    for ax in axarr[1:]:
        sns.despine(ax=ax, top=True, left=False, right=False, bottom=False)
