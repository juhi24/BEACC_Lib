# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt
from os import path

dtformat_default = '%d.%m. %H:%M'
dtformat_snex = '%Y %d %B %H UTC'

e = EventsCollection('cases/pip2015.csv', dtformat_snex)
e.autoimport_data(autoshift=False, autobias=False, rule='6min', varinterval=True)

plt.close('all')
plt.ion()

axarr = []

for i, c in enumerate(e.events.pluvio400.values):
    c.pluvio.shift_periods = -6
    c.pluvio.n_combined_intervals = 2
    plt.figure()
    c.pipv.use_flip = True
    c.clear_cache()
    axarr.append(c.density(rhomax=2000).plot(label='with flip'))
    c.pipv.use_flip = False
    c.clear_cache()
    c.density(rhomax=2000).plot(ax=axarr[i], label='no flip')
    plt.ylabel('bulk density')
    plt.legend()