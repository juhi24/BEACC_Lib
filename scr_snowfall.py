# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt

dtformat_default = '%d.%m. %H:%M'
dtformat_snex = '%d %B %H UTC'

e = EventsCollection('cases/test.csv', dtformat_snex)
e.autoimport_data(autoshift=False, autobias=False, rule='6min', varinterval=True)

for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
    c.pluvio.shift_periods = -6

#c.plot_velfitcoefs(rhomax=600, countmin=2000)
