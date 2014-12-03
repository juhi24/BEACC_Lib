# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt

e = EventsCollection('cases/cases_of_interest.csv', '%d.%m. %H:%M')
e.autoimport_data(autoshift=False, autobias=False, rule='6min', varinterval=True)

for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
    c.pluvio.shift_periods = -6
    axarr = c.plot()
    c.varinterval = False
    c.minimize_lsq()
    c.plot(axarr=axarr, label_suffix='6min')
    #c.minimize_lsq()
    #c.plot(pip=False)
    #c.pipv.plots(save=True)

c=e.events.pluvio200[1]
#c.plot()
#c.varinterval=True
#c.minimize_lsq()
#c.plot()
#c.ab=(0.1,2.1)
#vels = c.pipv.v(0.375, rule=c.pluvio.grouper())
