# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt
from os import path

dtformat_default = '%d.%m.%y %H:%M'
dtformat_snex = '%Y %d %B %H UTC'
dtformat_paper = '%Y %b %d %H:%M'

h5baecc_path = path.join(read.DATA_DIR, 'baecc.h5')
h5w1415path = path.join(read.DATA_DIR, 'dec-jan1415.h5')

def test_events():
    return events(casesfile_baecc='cases/pip2015test.csv',
                  casesfile_1415='cases/pip2015_14-15test.csv')

def events(casesfile_baecc='cases/pip2015.csv',
           casesfile_1415='cases/pip2015_14-15.csv'):
    e = EventsCollection(casesfile_baecc, dtformat_snex)
    e.autoimport_data(datafile=read.H5_PATH, autoshift=False, autobias=False,
                      rule='6min', varinterval=True)

    for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
        c.instr['pluvio'].shift_periods = -6
        c.instr['pluvio'].n_combined_intervals = 2

    e1415 = EventsCollection(casesfile_1415, dtformat_paper)
    e1415.autoimport_data(datafile=h5w1415path, autoshift=False, autobias=False,
                      rule='6min', varinterval=True)

    for c in np.append(e1415.events.pluvio200.values, e1415.events.pluvio400.values):
        c.instr['pluvio'].shift_periods = -5
        c.instr['pluvio'].n_combined_intervals = 2

    e.events = e.events.append(e1415.events, ignore_index=True)
    del(e1415)

    e.events['paper'] = e.events.pluvio200

    e.split_index()
    e.events = before_after_col(e.events, date=pd.datetime(2014,7,1),
                                datecol='start')
    return e


def pip2015events():
    e = events()
    e.events.paper[12] = e.events.pluvio400[12]
    e.events.pluvio200[12] = None
    return e

#e = pip2015events()

#plt.close('all')
#plt.ioff()
#plt.ion()
