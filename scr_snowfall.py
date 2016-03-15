# -*- coding: utf-8 -*-
"""
@author: Jussi Tiira
"""
from snowfall import *
import numpy as np
import matplotlib.pyplot as plt
from os import path

pluvio_comb_intervals = 2

dtformat_default = '%d.%m.%y %H:%M'
dtformat_snex = '%Y %d %B %H UTC'
dtformat_paper = '%Y %b %d %H:%M'

h5baecc_path = path.join(read.DATA_DIR, 'baecc.h5')
h5nov14_path = path.join(read.DATA_DIR, '2014nov1-23.h5')
h5w1415path = path.join(read.DATA_DIR, 'dec-jan1415.h5')

def test_events():
    return events(casesfile_baecc='cases/pip2015test.csv',
                  casesfile_1415='cases/pip2015_14-15test.csv')

def pluvio_config(e, tshift_minutes, n_comb_intervals):
    for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
        c.instr['pluvio'].shift_periods = tshift_minutes
        c.instr['pluvio'].n_combined_intervals = n_comb_intervals

def extra_events(e, extra_cases_file, extra_h5_file, *pluvio_conf_args):
    ee = EventsCollection(extra_cases_file, dtformat_paper)
    ee.autoimport_data(datafile=extra_h5_file, autoshift=False, autobias=False,
                      rule='6min', varinterval=True)
    pluvio_config(ee, *pluvio_conf_args)
    e.events = e.events.append(ee.events, ignore_index=True)
    del(ee)

def events(casesfile_baecc='cases/pip2015.csv',
           casesfile_nov14='cases/pip2015_nov14.csv',
           casesfile_1415='cases/pip2015_14-15.csv'):
    e = EventsCollection(casesfile_baecc, dtformat_paper)
    e.autoimport_data(datafile=read.H5_PATH, autoshift=False, autobias=False,
                      rule='6min', varinterval=True)
    pluvio_config(e, -6, pluvio_comb_intervals)
    extra_events(e, casesfile_nov14, h5nov14_path, -5, pluvio_comb_intervals)
    extra_events(e, casesfile_1415, h5w1415path, -5, pluvio_comb_intervals)
    e.events['paper'] = e.events.pluvio200
    e.split_index()
    e.events = before_after_col(e.events, date=pd.datetime(2014,7,1),
                                datecol='start')
    return e


def pip2015events():
    e = events()
    pref400 = e.events.pluvio_pref==400
    e.events.paper[pref400] = e.events.pluvio400[pref400] # TODO
    #e.events.pluvio200[12] = None # TODO
    return e

#e = pip2015events()

#plt.close('all')
#plt.ioff()
#plt.ion()
