# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt
from os import path

def pip2015events():
    dtformat_default = '%d.%m.%y %H:%M'
    dtformat_snex = '%Y %d %B %H UTC'
    dtformat_paper = '%Y %b %d %H:%M'
    h5baecc_path = path.join(read.DATA_DIR, 'baecc.h5')
    h5w1415path = path.join(read.DATA_DIR, 'dec-jan1415.h5')

    e = EventsCollection('cases/pip2015.csv', dtformat_snex)
    e.autoimport_data(datafile=read.H5_PATH, autoshift=False, autobias=False,
                      rule='6min', varinterval=True)

    for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
        c.instr['pluvio'].shift_periods = -6
        c.instr['pluvio'].n_combined_intervals = 2

    e1415 = EventsCollection('cases/pip2015_14-15.csv', dtformat_paper)
    e1415.autoimport_data(datafile=h5w1415path, autoshift=False, autobias=False,
                      rule='6min', varinterval=True)

    for c in np.append(e1415.events.pluvio200.values, e1415.events.pluvio400.values):
        c.instr['pluvio'].shift_periods = -5
        c.instr['pluvio'].n_combined_intervals = 2

    e.events = e.events.append(e1415.events, ignore_index=True)
    del(e1415)

    e.events['paper'] = e.events.pluvio200
    e.events.paper[12] = e.events.pluvio400[12]
    return e

#plt.close('all')
#plt.ioff()
#plt.ion()
