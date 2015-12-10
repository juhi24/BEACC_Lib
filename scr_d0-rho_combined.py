# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
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

savepath = '../results/pip2015'

#plt.close('all')
plt.ion()

def plot_d0_rho(data):
    plotkws = {'x': 'D_0_gamma',
               'y': 'density',
               'sizecol': 'count',
               #'groupby': 'case',
               #'c': 'case',
               'scale': 0.1,
               'colorbar': False,
               'xlim': [0.6,6],
               'ylim': [0,500]}
    ax = plot_pairs(data.loc['first'], c='blue', label='BAECC', **plotkws)
    plot_pairs(data.loc['second'], c='green', ax=ax, label='winter 2014-2015',
               **plotkws)
    plt.tight_layout()
    rho_d0_cols = ['density','D_0_gamma', 'count']
    rho_d0_data = data.loc[:, rho_d0_cols].dropna()
    rho_d0_data_baecc = rho_d0_data.loc['first']
    rho_d0_data_1415 = rho_d0_data.loc['second']
    #rho_d0 = fit.PolFit(x=rho_d0_data.D_0_gamma, y=rho_d0_data.density,
    #                    sigma=1/rho_d0_data['count'], xname='D_0')
    rho_d0_baecc = fit.PolFit(x=rho_d0_data_baecc.D_0_gamma,
                              y=rho_d0_data_baecc.density,
                              sigma=1/rho_d0_data_baecc['count'], xname='D_0')
    rho_d0_1415 = fit.PolFit(x=rho_d0_data_1415.D_0_gamma,
                             y=rho_d0_data_1415.density,
                             sigma=1/rho_d0_data_1415['count'], xname='D_0')
    #rho_d0.find_fit()
    #rho_d0.plot(ax=ax4)
    rho_d0_baecc.find_fit()
    rho_d0_1415.find_fit()
    rho_d0_baecc.plot(ax=ax)
    rho_d0_1415.plot(ax=ax)
    brandes = fit.PolFit(params=[178, -0.922])
    brandes.plot(ax=ax, label='Brandes et al.')
    ax.set_ylabel('$\\rho$, kg$\,$m$^{-3}$')
    ax.set_xlabel('$D_{0,\gamma}$, mm')
    plt.legend()
    return ax

#during_baecc = e.events.start<pd.datetime(2014,7,1)
#w14 = e.events[during_baecc]
#w1415 = e.events[-during_baecc]

data = e.summary(col='paper', split_date=pd.datetime(2014,7,1))
data = data.query('density<600 & count>1000 & b>0')
data_fltr = data[data.D_0 > 0.63]

figkws = {'dpi': 150, 'figsize': (5,6)}
fig = plt.figure(**figkws)
ax = plot_d0_rho(data)
fig_fltr = plt.figure(**figkws)
ax_fltr = plot_d0_rho(data_fltr)

if debug:
    savepath += '/test'
read.ensure_dir(savepath)
fig.savefig(path.join(savepath, 'rho_d0_combined.eps'))
fig_fltr.savefig(path.join(savepath, 'rho_d0_combined_d0fltr.eps'))