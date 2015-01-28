# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt

dtformat_snex = '%Y %d %B %H UTC'
dtformat_davide = '%d.%m.%y %H:%M'

e = EventsCollection('cases/2015.csv', dtformat_davide)
e.autoimport_data(autoshift=False, autobias=False, rule='5min', 
                  varinterval=True, datafile=['../DATA/test_winter1415.h5'])

basedir = '/home/jussitii/results/pip2015'

plt.close('all')

for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
    c.dsd.store_good_data() # improve performance by storing filtered dsd tables in memory
    c.pluvio.shift_periods = -6
    c.reset() # reset memory cache after changing pluvio timeshift


plt.ion()
fig1 = plt.figure(dpi=120)
ax1 = plt.gca()

e.plot_pairs(ax=ax1, c='density', sizecol='count', vmin=0, vmax=300,
             query='density<600 & count>1000 & b>0', colorbar=True, 
             xlim=[0.5,2.5])

plt.tight_layout()

fig2 = plt.figure(dpi=120)
ax2 = plt.gca()

e.plot_pairs(ax=ax2, x='D_0', y='density', sizecol='count', vmin=0, vmax=800,
             query='density<600 & count>1000 & b>0', colorbar=False, xlim=[0,8], ylim=[0,600])

plt.tight_layout()

fig4 = plt.figure(dpi=120)
ax4 = plt.gca()

e.plot_pairs(ax=ax4, x='D_0_gamma', y='density', sizecol='count', vmin=0, vmax=800,
             query='density<600 & count>1000 & b>0', colorbar=False, xlim=[0,8], ylim=[0,600])


plt.tight_layout()


fig3 = plt.figure(dpi=120)
ax3 = plt.gca()

e.plot_pairs(ax=ax3, x='D_max', y='b',c='density', sizecol='count', vmin=0, vmax=300,
             query='density<600 & count>1000 & b>0', colorbar=True)


plt.tight_layout()

brandes = read.PolFit(params=[178, -0.922])
brandes.plot(ax=ax2)
brandes.plot(ax=ax4)
