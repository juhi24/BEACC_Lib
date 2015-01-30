# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt
from os import path

def init_dataset(e, shift=-6):
    for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
        c.dsd.store_good_data() # improve performance by storing filtered dsd tables in memory
        c.pluvio.shift_periods = shift
        c.reset() # reset memory cache after changing pluvio timeshift

dtformat_snex = '%Y %d %B %H UTC'
dtformat_davide = '%d.%m.%y %H:%M'

w1415 = EventsCollection('cases/2015.csv', dtformat_davide)
w1415.autoimport_data(autoshift=False, autobias=False, rule='5min',
                      varinterval=True, datafile=['../DATA/winter1415.h5'])

w14 = EventsCollection('cases/pip2015.csv', dtformat_snex)
w14.autoimport_data(autoshift=False, autobias=False, rule='5min',
                    varinterval=True, datafile=['../DATA/baecc.h5'])

init_dataset(w14, shift=-6)
init_dataset(w1415, shift=-5)

basedir = '/home/jussitii/results/pip2015'

plt.close('all')
plt.ioff()

e = combine_datasets(w14, w1415)
#del w1415
#del w14

fig4 = plt.figure(dpi=120)
ax4 = w1415.plot_pairs(c='k', x='D_0_gamma', y='density', sizecol='count', scale=0.5,
             query='density<600 & count>1000 & b>0', colorbar=False, xlim=[0,6], ylim=[0,500])
w14.plot_pairs(ax=ax4, c='m', x='D_0_gamma', y='density', sizecol='count', scale=0.5,
             query='density<600 & count>1000 & b>0', colorbar=False, xlim=[0,6], ylim=[0,500])
plt.tight_layout()

brandes = read.PolFit(params=[178, -0.922])
brandes.plot(ax=ax4, label='Brandes et al.')

s = e.summary()
rho_d0 = read.PolFit(x=s.D_0_gamma, y=s.density, sigma=1/s['count'], xname='D_0')
rho_d0.find_fit()
rho_d0.plot(ax=ax4)

plt.legend()
