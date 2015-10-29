# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt
from os import path
import seaborn as sns
import fit

from scr_snowfall import pip2015events

plt.close('all')
plt.ion()

e = pip2015events()
comb = e.events.paper.sum()
del(e)

limslist = limitslist((0, 150, 300, 800))
n_ranges = len(limslist)
separate = False

basepath = '../results/pip2015/paper/psd'

def d0_nw_data(c):
    d0 = c.d_0()
    nw = c.n_w()
    nw_log = np.log10(nw)
    nw_log.name = 'log_nw'
    df = read.merge_multiseries(d0, nw, nw_log)
    return df.dropna()

def d0_nw_datarho(c, rholimits=[0,150,300,800]):
    data = d0_nw_data(c)
    return c.group_by_density(data, rholimits)

def d0_nw_paper(c, rholimits=[0,150,300,800]):
    rhopairs = limitslist(rholimits)
    data = d0_nw_datarho(c=c, rholimits=rholimits)
    fig, axarr = plt.subplots(nrows=1, ncols=3, dpi=100, figsize=(14,5),
                              tight_layout=True, sharex=True, sharey=True)
    for i, (rhomin, rhomax) in enumerate(rhopairs):
        ax = axarr[i]
        datarho = data[data.rhomin==rhomin]
        #datarho.plot(kind='scatter', x='D_0', y='N_w', ax=axarr[i], logy=True)
        #datarho.plot(kind='scatter', x='D_0', y='log_nw', ax=axarr[i])
        lfit = fit.LinFit(x=datarho.D_0, y=datarho.log_nw, xname='D_0')
        lfit.find_fit()
        lfit.plot(source_style='raw', ax=ax)
        #efit = fit.ExponentialFit(x=datarho.D_0, y=datarho.N_w)
        #efit.find_fit()
        #efit.plot(source_style='raw', ax=ax)
        ax.legend()
        ax.set_xlabel('$D_0$')
        ax.set_ylabel('$log(N_w)$')
    return fig, axarr

def d0_nw_plots(c):
    data = d0_nw_data()
    datarho = c.group_by_density(data, [0,150,300,800])
    f, ax = plt.subplots(dpi=120, tight_layout=True)
    nwcol = 'log_nw'
    df.plot(kind='hexbin', x='D_0', y=nwcol, logy=False, ax=ax)
    sns.jointplot(x='D_0', y=nwcol, data=data, kind='kde', stat_func=None,
                  joint_kws={'kde_kws':{'bw':0.001}})
    sns.jointplot(x='D_0', y=nwcol, data=data, kind='hex', stat_func=None)
    sns.jointplot(x='D_0', y=nwcol, data=data, kind='scatter', stat_func=None)
    sns.lmplot(x='D_0', y=nwcol, data=datarho, col='rhomin', col_order=[0,150,300])
    sns.lmplot(x='D_0', y=nwcol, data=datarho, hue='rhomin', hue_order=[0,150,300])
    return datarho

def psds_in_rho_range(c):
    d0 = c.d_0()
    nt = c.n_t()
    nw = c.n_w()
    #data = read.merge_multiseries(d0, nt, nw)
    n = c.intervalled(c.instr['dsd'].psd)
    y = n.div(nw, axis=0)
    d = y*0+y.columns.values
    xd = {'D':d, '$D D_0^{-1}$':d.div(d0.replace({0: np.nan}), axis=0)}
    for xlabel, x in xd.items():
        if not separate:
            fig, axarr = plt.subplots(1, n_ranges, dpi=100, sharex=True, sharey=True,
                                      tight_layout=True, figsize=(n_ranges*6, 6))
        for i, rhorange in enumerate(limslist):
            if separate:
                fig, ax = plt.subplots(dpi=100)
            else:
                ax = axarr[i]
            rhomin = rhorange[0]
            rhomax = rhorange[1]
            xpart = c.data_in_density_range(x, rhomin, rhomax)
            ypart = c.data_in_density_range(y, rhomin, rhomax)
            if xpart.empty or ypart.empty:
                ax.semilogy()
                continue
            ax.semilogy(xpart, ypart, marker='o', linestyle='None', color='black')
            ax.set_ylabel('$N_D N_w^{-1}$')
            ax.set_xlabel(xlabel)
            ax.set_title('$%s < \\rho < %s$' % (rhomin, rhomax))
            ax.grid(True)
        plt.axis([0, 5, 10e-5, 1000])
        ax.set_title('$\\rho > ' + str(rhomin) + '$')

#psds_in_rho_range(comb)
#datarho = d0_nw_plots(comb)
data = d0_nw_datarho(comb)
fig, axarr = d0_nw_paper(comb)