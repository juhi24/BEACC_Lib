# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt

e = EventsCollection('cases/cases_of_interest.csv', '%d.%m. %H:%M')
e.autoimport_data(autoshift=False, autobias=False, rule='6min', varinterval=True)

axarrarr=[]

for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
    c.pluvio.shift_periods = -6
    axarr = c.plot(label_suffix='varying intervals')
    c.varinterval = False
    c.minimize_lsq()
    c.plot(axarr=axarr, style='o--', label_suffix='6min')
    axarrarr.append(axarr)
    #c.minimize_lsq()
    #c.plot(pip=False)
    #c.pipv.plots(save=True)
