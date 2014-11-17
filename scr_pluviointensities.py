# -*- coding: utf-8 -*-
"""
@author: jussitii
"""

from snowfall import *
import numpy as np
import matplotlib.pyplot as plt

e = EventsCollection('cases/cases_of_interest.csv', '%d.%m. %H:%M')
e.autoimport_data(autoshift=False, rule='6min')

for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
    c.pluvio.shift_periods = -4

for index, event in e.events.iterrows():
    a = []
    for p in ['pluvio200','pluvio400']:
        i = event[p].pluvio.intensity(rule='6min')
        i.name = p
        i.index.name = 'datetime'
        a.append(i)
    out=pd.DataFrame(a).T
    out.plot()
    out.to_csv('/home/jussitii/' + event.start.date().strftime('%y%m%d') + 'pluvio_intensity.csv')
