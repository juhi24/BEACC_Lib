# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt

e = EventsCollection('cases/cases_of_interest.csv', '%d.%m. %H:%M')
e.autoimport_data(autoshift=False, rule='6min')

for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
    c.pluvio.shift_periods = -4
    #c.plot(pip=False)
    #c.pipv.plots(save=True)
