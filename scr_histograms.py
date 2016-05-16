# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
import snowfall as sf
import read
from fit import set_plot_style
import numpy as np
import matplotlib.pyplot as plt
from os import path
import seaborn as sns
import gc
from scipy.stats import gaussian_kde
import pandas as pd

from scr_snowfall import pip2015events, test_events, param_table, rholimits

debug = False
d0_col = 'D_0_gamma'

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

data = param_table(debug=debug)


def kde_peak(data, xlim=(-1,1), n_samples=500):
    kde = gaussian_kde(data)
    xs = np.linspace(*xlim, num=n_samples)
    ys = kde(xs)
    index = np.argmax(ys)
    return xs[index], ys[index]


def split_hist(data, **kws):
    data = sf.split_index(data, names=('first', 'second'))
    return plt.hist([data.loc['first'], data.loc['second']], stacked=True,
                    rwidth=1, **kws)


def subplots(n_plots=1):
    return plt.subplots(n_plots, sharex=True, sharey=True,
                        tight_layout=False, figsize=(4,5))


def plots(d0, mu, nw, axd, axm, axn, label=None, title=None, **kwtitle):
    rng = (0,6)
    sns.distplot(d0, ax=axd, label=label, bins=12, 
                 hist_kws={'range':rng}, **kwargs)
    axd.set_xlim(rng)
    axd.yaxis.set_ticks(np.arange(0, 1.5, 0.5))
    axd.set_xlabel('$D_0$, mm')
    rng = (-3, 7)
    axd.axis([None, None, None, 1.5])
    sns.distplot(mu, ax=axm, label=label, bins=20,
                       hist_kws={'range':rng}, **kwargs)
    axm.set_xlim(rng)
    axm.yaxis.set_ticks(np.arange(0, 0.7, 0.2))
    axm.set_xlabel('$\mu$')
    axm.axis([None, None, None, 0.5])
    sns.distplot(nw, ax=axn, label=label,
                       bins=10**np.linspace(0,6,20), kde=False, rug=True)
    axn.set_xscale('log')
    axn.set_xlabel('$N_w$, mm$^{-1}$m$^{-3}$')
    #axn.yaxis.get_major_ticks()[-1].label1.set_visible(False)
    axn.yaxis.set_ticks(np.arange(0, 140, 50))
    axn.axis([1e1, None, None, 140])
    if title is not None:
        for ax in (axd, axm, axn):
            ax.set_title(title, **kwtitle)
        for ax in (axd, axm):
            ax.legend().set_visible(False)


def hist_data(case):
    return read.merge_multiseries(case.d_0(), case.mu(), case.n_w())


def casewise_hist():
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
    del(e)
    gc.collect()
    for f in (fd, fm, fn):
        sf.remove_subplot_gaps(f, axis='col')
    fd.savefig(path.join(savedir, 'd0_cases' + tld), **savekws)
    fm.savefig(path.join(savedir, 'mu_cases' + tld), **savekws)
    fn.savefig(path.join(savedir, 'nw_cases' + tld), **savekws)


limslist = sf.limitslist(rholimits)
n_ranges = len(limslist)

fdd, axarrdd = subplots(n_ranges)
fmd, axarrmd = subplots(n_ranges)
fnd, axarrnd = subplots(n_ranges)
titlekws = {'y': 0.85, 'fontdict': {'verticalalignment': 'top'}}

stats = pd.DataFrame(index=rholimits[:-1])
kde_peak_d0 = [None, None, None]
kde_peak_mu = [None, None, None]
data_by_rho = sf.df_ranges(data, 'density', limslist).groupby('density_range_min')

for i, (rhomin, rhomax) in enumerate(limslist):
    dat = data.query('%s<density<%s' % (rhomin, rhomax))
    d0 = dat[d0_col].dropna()
    mu = dat.mu.dropna()
    nw = dat.N_w.dropna()
    kde_peak_d0[i] = kde_peak(d0, xlim=(0,6))[0]
    kde_peak_mu[i] = kde_peak(mu, xlim=(-4, 4))[0]
    limitstr = '$%s < \\rho \leq %s$' % (rhomin*read.RHO_SCALE, rhomax*read.RHO_SCALE)
    plots(d0, mu, nw, axarrdd[i], axarrmd[i], axarrnd[i], title=limitstr, **titlekws)

stats['d0_peak'] = kde_peak_d0
stats['mu_peak'] = kde_peak_mu
stats['d0_median'] = data_by_rho[d0_col].median()
stats['mu_median'] = data_by_rho.mu.median()
stats['nw_median'] = data_by_rho.N_w.median()
stats['d0_std'] = data_by_rho[d0_col].std()
stats['mu_std'] = data_by_rho.mu.std()
stats.to_csv(path.join(savedir, 'psd_stats.csv'))

for ax in (axarrdd[-1], axarrmd[-1], axarrnd[-1]):
    ax.set_title('$\\rho > %s$' % (rhomin*read.RHO_SCALE), **titlekws)
    #ax.set_title('$\\rho > %s$' % rhomin, **titlekws)

for ax in (axarrdd[1], axarrmd[1]):
    ax.set_ylabel('PDF')
axarrnd[1].set_ylabel('Frequency')

for f in (fdd, fmd, fnd):
    sf.remove_subplot_gaps(f, axis='col')
fnd.subplots_adjust(left=0.15) # prevent ylabel cutoff

tld = '.png'
savekws = {'dpi': 400}

fdd.savefig(path.join(savedir, 'd0_rho' + tld), **savekws)
fmd.savefig(path.join(savedir, 'mu_rho' + tld), **savekws)
fnd.savefig(path.join(savedir, 'nw_rho' + tld), **savekws)
fdd.savefig(path.join(paperdir, 'hist_d0' + tld), **savekws)
fmd.savefig(path.join(paperdir, 'hist_mu' + tld), **savekws)
fnd.savefig(path.join(paperdir, 'hist_nw' + tld), **savekws)

# Wat?
for axarr in (axarrdd, axarrmd, axarrnd):
    for ax in axarr[1:]:
        sns.despine(ax=ax, top=True, left=False, right=False, bottom=False)
