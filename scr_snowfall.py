# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt
from os import path

dtformat_snex = '%Y %d %B %H UTC'
dtformat_davide = '%d.%m.%y %H:%M'

e = EventsCollection('cases/2015.csv', dtformat_davide)
e.autoimport_data(autoshift=False, autobias=False, rule='5min', 
                  varinterval=True, datafile=['../DATA/test_winter1415.h5'])

basedir = '/home/jussitii/results/pip2015'

plt.close('all')
plt.ioff()

for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
    c.dsd.store_good_data() # improve performance by storing filtered dsd tables in memory
    c.pluvio.shift_periods = -5
    c.reset() # reset memory cache after changing pluvio timeshift
    c.pipv.find_fits(rule=c.rule)
    