# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt

dtformat_default = '%d.%m. %H:%M'
dtformat_snex = '%d %B %H UTC'

e = EventsCollection('cases/pip2015.csv', dtformat_snex)
e.autoimport_data(autoshift=False, autobias=False, rule='6min', varinterval=True)

basedir = '/home/jussitii/results/pip2015'

plt.ion()
#ax = plt.gca()
#plt.figure(dpi=120)

for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
    c.dsd.store_good_data() # improve performance by storing filtered dsd tables in memory
    c.pluvio.shift_periods = -6
    c.reset() # reset on changing pluvio timeshift

#fig = plt.figure()
ax=plt.gca()

markers = ['o', 's', '^', 'v', 'D', '*']

#for i,c in enumerate(e.events.pluvio400):
#    #c.pipv.find_fits(rule=c.rule)
#    c.plot_d0_bv(countmin=1, rhomax=600, colorbar=False, edgecolors='none',
#                 ax=ax, countscale=3, legend=True, marker=markers[i])