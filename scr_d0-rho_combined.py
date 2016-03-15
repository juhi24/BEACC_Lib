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
plt.ioff()

# dry and wet snows, from original
magono65 = fit.PolFit(params=[20, -2])
# "Fit of data from Magano and Nakamura (1965) for dry snowflakes", fit by Holroyd
magono65dry = fit.PolFit(params=[22, -1.5])
holroyd71 = fit.PolFit(params=[170, -1]) # from Brandes
# "Used by Schaller et al. (1982) for brightband modeling" from Fabry
schaller82 = fit.PolFit(params=[64, -0.65])
 # from original
muramoto95 = fit.PolFit(params=[48, -0.406])
# "Measurements from Switzerland by Barthazy (1997, personal communication)" from Fabry
barthazy97 = fit.PolFit(params=[18, -0.8])
# from Brandes
fabry99 = fit.PolFit(params=[150, -1])
# from Brandes
heymsfield04 = fit.PolFit(params=[104, -0.95])
brandes07 = fit.PolFit(params=[178, -0.922])

color1 = 'gray'
color2 = 'red'
fits_to_plot = {#magono65: {'color':color1, 'linestyle':'--', 'label':'Magono and Nakamura (1965)'},
                #magono65dry: {'label':'Magono and Nakamura (1965), dry snow'}
                #holroyd71: {'color':color1, 'linestyle':':', 'label':'Holroyd (1971)'},
                #schaller82: {'color':color1, 'linestyle':'-.', 'label':'Schaller et al. (1982)'},
                #muramoto95: {'color':color1, 'linestyle':'-', 'label':'Muramoto et al. (1995)'},
                #barthazy97: {'color':color, linestyle:'--', 'label':'Barthazy (1997)'},
                #fabry99: {'color':color2, 'linestyle':'-', 'label':'Fabry and Szyrmer (1999)'},
                #heymsfield04: {'color':color2, 'linestyle':'--', 'label':'Heymsfield et al. (2004)'},
                brandes07: {'color':'black', 'linestyle':'--', 'label':'Brandes et al. (2007)'}}

def plot_density_histogram(data, bins=60, **kws):
    ax = data.density.hist(bins=bins, **kws)
    read.rho_scale(ax.xaxis)
    ax.set_xlabel('bulk density')
    ax.set_ylabel('frequency')

def plot_d0_rho(data):
    plotkws = {'x': 'D_0_gamma',
               'y': 'density',
               'sizecol': 'count',
               #'groupby': 'case',
               #'c': 'case',
               'scale': 0.1,
               'colorbar': False,
               'xlim': [0.6,6],
               'ylim': [0,450]}
    ax = plot_pairs(data.loc['first'], c='blue', label='BAECC', **plotkws)
    plot_pairs(data.loc['second'], c='green', ax=ax, label='winter 2014-2015',
               **plotkws)
    plt.tight_layout()
    rho_d0_cols = ['density','D_0_gamma', 'count']
    rho_d0_data = data.loc[:, rho_d0_cols].dropna()
    rho_d0_data_baecc = rho_d0_data.loc['first']
    rho_d0_data_1415 = rho_d0_data.loc['second']
    rho_d0 = fit.PolFit(x=rho_d0_data.D_0_gamma, y=rho_d0_data.density,
                        sigma=1/rho_d0_data['count'], xname='D_0')
    rho_d0_baecc = fit.PolFit(x=rho_d0_data_baecc.D_0_gamma,
                              y=rho_d0_data_baecc.density,
                              sigma=1/rho_d0_data_baecc['count'], xname='D_0')
    rho_d0_1415 = fit.PolFit(x=rho_d0_data_1415.D_0_gamma,
                             y=rho_d0_data_1415.density,
                             sigma=1/rho_d0_data_1415['count'], xname='D_0')
    rho_d0.find_fit(loglog=True)
    rho_d0_baecc.find_fit(loglog=True)
    rho_d0_1415.find_fit(loglog=True)
    rho_d0_baecc.plot(ax=ax)
    rho_d0_1415.plot(ax=ax)
    rho_d0.plot(ax=ax, color='black', label='all cases: $%s$' % str(rho_d0))
    for key, kws in fits_to_plot.items():
        key.plot(ax=ax, **kws)
    read.rho_scale(ax.yaxis)
    ax.set_ylabel('$\\rho$, ' + read.RHO_UNITS)
    ax.set_xlabel('$D_{0,\gamma}$, mm')
    plt.legend()
    return ax

#during_baecc = e.events.start<pd.datetime(2014,7,1)
#w14 = e.events[during_baecc]
#w1415 = e.events[-during_baecc]

data = e.summary(col='paper', split_date=pd.datetime(2014,7,1))
data = data.query('density<600 & count>1000 & b>0')
data_fltr = data[data.D_0 > 0.63]

figkws = {'dpi': 150, 'figsize': (5,5)}
fig = plt.figure(**figkws)
ax = plot_d0_rho(data)
fig_fltr = plt.figure(**figkws)
ax_fltr = plot_d0_rho(data_fltr)

if debug:
    savepath += '/test'
paperpath = path.join(savepath, 'paper')
read.ensure_dir(paperpath)
fig.savefig(path.join(savepath, 'rho_d0_combined.eps'))
fig_fltr.savefig(path.join(savepath, 'rho_d0_combined_d0fltr.eps'))
fig_fltr.savefig(path.join(paperpath, 'd0_rho.eps'))