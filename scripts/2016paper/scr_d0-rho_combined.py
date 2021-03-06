# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
import snowfall as sf
import read
import numpy as np
import matplotlib.pyplot as plt
from os import path
import fit
from scipy.special import gamma

from scr_snowfall import param_table

debug = False
savepath = '../results/pip2015'
d0_col = 'D_0_gamma'
#plt.close('all')
plt.ioff()

fit_scale = [read.RHO_SCALE,1]

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

def prep_d0_rho(data):
    rho_d0 = fit.PolFit(x=data[d0_col], y=data.density,
                        sigma=1/data['count'], xname='D_0',
                        disp_scale=fit_scale)
    rho_d0.find_fit(loglog=True)
    return rho_d0

def prepare_d0_rho(data):
    rho_d0_cols = ['density',d0_col, 'count']
    rho_d0_data = data.loc[:, rho_d0_cols].dropna()
    rho_d0 = prep_d0_rho(rho_d0_data)
    rho_d0_baecc = prep_d0_rho(rho_d0_data.loc['first'])
    rho_d0_1415 = prep_d0_rho(rho_d0_data.loc['second'])
    return rho_d0, rho_d0_baecc, rho_d0_1415

def plot_d0_rho(data):
    plotkws = {'x': d0_col,
               'y': 'density',
               'sizecol': 'count',
               #'groupby': 'case',
               #'c': 'case',
               'scale': 0.1,
               'colorbar': False,
               'xlim': [0.5,6],
               'ylim': [0,450],
               'alpha': 0.5}
    ax = sf.plot_pairs(data.loc['first'], c=(.6, .6 , .92, .8), label='BAECC', **plotkws)
    sf.plot_pairs(data.loc['second'], c=(.2, .92, .2, .8), ax=ax, label='winter 2014-2015',
               **plotkws)
    plt.tight_layout()
    rho_d0, rho_d0_baecc, rho_d0_1415 = prepare_d0_rho(data)
    rho_d0_baecc.plot(ax=ax)
    rho_d0_1415.plot(ax=ax)
    rho_d0.plot(ax=ax, color='black', label='all cases: $%s$' % str(rho_d0))
    for key, kws in fits_to_plot.items():
        key.plot(ax=ax, **kws)
    read.rho_scale(ax.yaxis)
    ax.set_ylabel('$\\rho$, ' + read.RHO_UNITS)
    ax.set_xlabel('$D_0$, mm')
    ax.set_xticks((0.5, 1, 2, 3, 4, 5, 6))
    plt.legend()
    return ax, rho_d0, rho_d0_baecc, rho_d0_1415

def mass_dim(rho_d0, b_v=0.2):
    a_d0, b_d0 = rho_d0.params
    a_d0 = a_d0/1000*10**b_d0
    beta = b_d0 + 3
    alpha = np.pi/6*3.67**b_d0*a_d0*gamma(b_v+4)/gamma(b_v+b_d0+4)
    return fit.PolFit(params=(alpha, beta), xname='D')

data_fltr = param_table(debug=debug)

figkws = {'dpi': 150, 'figsize': (5,5)}
#fig = plt.figure(**figkws)
#ax = plot_d0_rho(data)
fig_fltr = plt.figure(**figkws)
ax_fltr, rho_d0, rho_d0_baecc, rho_d0_1415 = plot_d0_rho(data_fltr)

m_d = mass_dim(rho_d0)
m_d_baecc = mass_dim(rho_d0_baecc)
m_d_1415 = mass_dim(rho_d0_1415)

if debug:
    savepath += '/test'
paperpath = path.join(savepath, 'paper')
read.ensure_dir(paperpath)
#fig.savefig(path.join(savepath, 'rho_d0_combined.eps'))
fig_fltr.savefig(path.join(savepath, 'rho_d0_combined_d0fltr.eps'))
fig_fltr.savefig(path.join(paperpath, 'd0_rho.eps'))
fig_fltr.savefig(path.join(paperpath, 'd0_rho.png'))