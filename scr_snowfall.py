# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt

e = EventsCollection('cases/rain.csv', '%d.%m. %H:%M')
e.autoimport_data(autoshift=False, autobias=False, rule='6min', varinterval=True)

for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
    c.pluvio.shift_periods = -6

#c.plot_velfitcoefs(rhomax=600, countmin=2000)
