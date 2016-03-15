# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
from fit import set_plot_style
import numpy as np
import matplotlib.pyplot as plt
from os import path
import seaborn as sns

from scr_snowfall import pip2015events, test_events

debug = False
rholims = (0, 100, 200, 800)
#rholims = (0, 150, 300, 800)

#sns.set_context('talk')
major_size = 8
set_plot_style(**{'xtick.major.size': major_size,
                  'ytick.major.size': major_size})

plt.close('all')
plt.ioff()
kwargs = {'kde': True, 'rug': True, 'kde_kws': {'label': 'KDE'}}
resultsdir = '../results/pip2015'
paperdir = read.ensure_dir(path.join(resultsdir, 'paper'))
savedir = path.join(resultsdir, 'hist')
if debug:
    savedir += '/test'
read.ensure_dir(savedir)

def split_hist(data, **kwargs):
    data = split_index(data, names=('first', 'second'))
    return plt.hist([data.loc['first'], data.loc['second']], stacked=True,
                    rwidth=1, **kwargs)

def subplots(n_plots=1):
    return plt.subplots(n_plots, sharex=True, sharey=True,
                        tight_layout=False, figsize=(4,5))

def plots(data, axd, axm, axn, label=None, title=None, **kwtitle):
    rng = (0,6)
    sns.distplot(data.D_0.dropna(), ax=axd, label=label, bins=12, 
                 hist_kws={'range':rng}, **kwargs)
    axd.set_xlim(rng)
    axd.yaxis.set_ticks(np.arange(0, 1.5, 0.5))
    axd.set_xlabel('$D_0$, mm')
    rng = (-2, 8)
    axd.axis([None, None, None, 1.5])
    sns.distplot(data.mu.dropna(), ax=axm, label=label, bins=20,
                       hist_kws={'range':rng}, **kwargs)
    axm.set_xlim(rng)
    axm.yaxis.set_ticks(np.arange(0, 0.7, 0.2))
    axm.set_xlabel('$\mu$')
    axm.axis([None, None, None, 0.5])
    sns.distplot(data.N_w.dropna(), ax=axn, label=label,
                       bins=10**np.linspace(0,6,20), kde=False, rug=True)
    axn.set_xscale('log')
    axn.set_xlabel('$N_w$, mm$^{-1}$m$^{-3}$')
    axn.yaxis.get_major_ticks()[-1].label1.set_visible(False)
    if title is not None:
        for ax in (axd, axm, axn):
            ax.set_title(title, **kwtitle)
        for ax in (axd, axm):
            ax.legend().set_visible(False)

def hist_data(case):
    data = read.merge_multiseries(case.d_0(), case.mu(), case.n_w())
    return d0fltr(data, apply=True)

if debug:
    e = test_events()
else:
    e = pip2015events()

n_cases = e.events.paper.count()
fd, axarrd = subplots(n_cases)
fm, axarrm = subplots(n_cases)
fn, axarrn = subplots(n_cases)

for i, c in enumerate(e.events.paper.values):
    print(c)
    data = hist_data(c)
    plots(data, axarrd[i], axarrm[i], axarrn[i], title=c.dtstr(), y=0.85,
          fontdict={'verticalalignment': 'top', 'fontsize': 10})

#c = e.events.paper.sum()
#del(e)
limslist = limitslist(rholims)
n_ranges = len(limslist)

fdd, axarrdd = subplots(n_ranges)
fmd, axarrmd = subplots(n_ranges)
fnd, axarrnd = subplots(n_ranges)
data = e.summary(col='paper', split_date=pd.datetime(2014,7,1))
data = d0fltr(data, apply=True)
titlekws = {'y': 0.85, 'fontdict': {'verticalalignment': 'top'}}

for i, (rhomin, rhomax) in enumerate(limslist):
    dat = data.query('%s<density<%s' % (rhomin, rhomax))
    limitstr = '$%s < \\rho \leq %s$' % (rhomin*read.RHO_SCALE, rhomax*read.RHO_SCALE)
    plots(dat, axarrdd[i], axarrmd[i], axarrnd[i], title=limitstr, **titlekws)

for ax in (axarrdd[-1], axarrmd[-1], axarrnd[-1]):
    ax.set_title('$\\rho > %s$' % (rhomin*read.RHO_SCALE), **titlekws)
    #ax.set_title('$\\rho > %s$' % rhomin, **titlekws)

for ax in (axarrdd[1], axarrmd[1]):
    ax.set_ylabel('PDF')
axarrnd[1].set_ylabel('Frequency')

for f in (fd, fm, fn, fdd, fmd, fnd):
    remove_subplot_gaps(f, axis='col')

tld = '.png'
savekws = {'dpi': 400}

fd.savefig(path.join(savedir, 'd0_cases' + tld), **savekws)
fm.savefig(path.join(savedir, 'mu_cases' + tld), **savekws)
fn.savefig(path.join(savedir, 'nw_cases' + tld), **savekws)
fdd.savefig(path.join(savedir, 'd0_rho' + tld), **savekws)
fmd.savefig(path.join(savedir, 'mu_rho' + tld), **savekws)
fnd.savefig(path.join(savedir, 'nw_rho' + tld), **savekws)
fdd.savefig(path.join(paperdir, 'hist_d0' + tld), **savekws)
fmd.savefig(path.join(paperdir, 'hist_mu' + tld), **savekws)
fnd.savefig(path.join(paperdir, 'hist_nw' + tld), **savekws)

for axarr in (axarrd, axarrm, axarrn, axarrdd, axarrmd, axarrnd):
    for ax in axarr[1:]:
        sns.despine(ax=ax, top=True, left=False, right=False, bottom=False)
