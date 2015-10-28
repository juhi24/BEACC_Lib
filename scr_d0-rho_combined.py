# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt
from os import path

from scr_snowfall import pip2015events

e = pip2015events()

basedir = '/home/jussitii/results/pip2015'

plt.close('all')
plt.ioff()

during_baecc = e.events.start<pd.datetime(2014,7,1)
w14 = e.events[during_baecc]
w1415 = e.events[-during_baecc]

plotkws = {'x': 'D_0_gamma', 'y': 'density', 'sizecol': 'count', 'scale': 0.5,
           'query': 'density<600 & count>1000 & b>0', 'colorbar': False,
           'xlim': [0,6], 'ylim': [0,500]}

fig4 = plt.figure(dpi=400)
ax4 = w1415.plot_pairs(c='k', **plotkws)
w14.plot_pairs(ax=ax4, c='m', **plotkws)
plt.tight_layout()

brandes = read.PolFit(params=[178, -0.922])
brandes.plot(ax=ax4, label='Brandes et al.')

s = e.summary(col='paper')
rho_d0 = read.PolFit(x=s.D_0_gamma, y=s.density, sigma=1/s['count'], xname='D_0')
rho_d0.find_fit()
rho_d0.plot(ax=ax4)

plt.legend()
