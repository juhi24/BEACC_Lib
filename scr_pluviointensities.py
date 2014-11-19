# -*- coding: utf-8 -*-
"""
@author: jussitii
"""

from snowfall import *
import os
import numpy as np
import matplotlib.pyplot as plt

e = EventsCollection('cases/cases_of_interest.csv', '%d.%m. %H:%M')
e.autoimport_data(autoshift=False, rule='6min')

home = os.curdir
if 'HOME' in os.environ:
    home = os.environ['HOME']

for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
    c.pluvio.shift_periods = -4

for index, event in e.events.iterrows():
    iarr = []
    darr = []
    dfarr = []
    for p in ['pluvio200', 'pluvio400']:
        intens = event[p].pluvio.intensity(rule='6min')
        dens = event[p].density()
        for var in [intens, dens]:    
            var.name = p
            var.index.name = 'datetime'
        iarr.append(intens)
        darr.append(dens)
    for varr in [iarr, darr]:
        df = pd.DataFrame(varr).T
        df.plot()
        dfarr.append(df)
    dtstr = event.end.date().strftime('%y%m%d')
    for i, label in enumerate(['pluvio_intensity', 'density']):
        dfarr[i].to_csv(os.path.join(home, dtstr + label + '.csv'))
