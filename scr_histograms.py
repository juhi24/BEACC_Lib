# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt
from os import path
import seaborn as sns

from scr_snowfall import pip2015events

sns.set_style('ticks')

plt.close('all')
#plt.ioff()
plt.ion()
kwargs = {'kde': True, 'rug': True, 'kde_kws': {'label': 'KDE'}}
resultsdir = read.ensure_dir('../results/pip2015/hist')

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

e = pip2015events()

n_cases = e.events.paper.count()
fd, axarrd = subplots(n_cases)
fm, axarrm = subplots(n_cases)
fn, axarrn = subplots(n_cases)

for i, c in enumerate(e.events.paper.values):
    data = read.merge_multiseries(c.d_0(), c.mu(), c.n_w())
    plots(data, axarrd[i], axarrm[i], axarrn[i], title=c.dtstr(), y=0.85,
          fontdict={'verticalalignment': 'top', 'fontsize': 10})

c = e.events.paper.sum()
del(e)
limslist = limitslist((0, 150, 300, 800))
n_ranges = len(limslist)

fdd, axarrdd = subplots(n_ranges)
fmd, axarrmd = subplots(n_ranges)
fnd, axarrnd = subplots(n_ranges)
data = read.merge_multiseries(c.d_0(), c.mu(), c.n_w())

for i, lims in enumerate(limslist):
    dat = c.data_in_density_range(data, *lims)
    limitstr = '$%s < \\rho < %s$' % (lims[0], lims[1])
    plots(dat, axarrdd[i], axarrmd[i], axarrnd[i], title=limitstr, y=0.9, fontdict={'verticalalignment': 'top'})

for ax in (axarrdd[-1], axarrmd[-1], axarrnd[-1]):
    ax.set_title('$\\rho > %s$' % lims[0], y=0.9, fontdict={'verticalalignment': 'top'})

for f in (fd, fm, fn, fdd, fmd, fnd):
    remove_gaps(f)

fd.savefig(path.join(resultsdir, 'd0_cases.eps'))
fm.savefig(path.join(resultsdir, 'mu_cases.eps'))
fn.savefig(path.join(resultsdir, 'nw_cases.eps'))
fdd.savefig(path.join(resultsdir, 'd0_rho.eps'))
fmd.savefig(path.join(resultsdir, 'mu_rho.eps'))
fnd.savefig(path.join(resultsdir, 'nw_rho.eps'))

for axarr in (axarrd, axarrm, axarrn, axarrdd, axarrmd, axarrnd):
    for ax in axarr[1:]:
        sns.despine(ax=ax, top=True, left=False, right=False, bottom=False)
