# -*- coding: utf-8 -*-
"""
Paper specific settings and helper functions.
@author: Jussi Tiira
"""
import snowfall as sf
import read
import numpy as np
import pandas as pd
from os import path
import gc

N_COMB_INTERVALS = 2

dtformat_default = '%d.%m.%y %H:%M'
dtformat_snex = '%Y %d %B %H UTC'
dtformat_paper = '%Y %b %d %H:%M'
#QSTR = 'density<600 & count>800 & b>0'
QSTR = 'D_0_gamma>0.6 & intensity>0.2 & count>800 & density==density'
rholimits = (0, 100, 200, 1000)
#rholimits = (0, 150, 300, 800)
home = path.expanduser('~')
resultspath = path.join(home, 'results', 'pip2015')
paperpath = path.join(resultspath, 'paper')
paths = {'results': read.ensure_dir(resultspath),
         'paper': read.ensure_dir(paperpath),
         'tables': read.ensure_dir(path.join(paperpath, 'tables'))}
files = {'h5nov14': path.join(read.DATA_DIR, '2014nov1-23.h5'),
         'h5w1415': path.join(read.DATA_DIR, 'dec-jan1415.h5'),
         'h5baecc': read.H5_PATH,
         'cbaecc': 'cases/pip2015.csv',
         'c14nov': 'cases/pip2015_nov14.csv',
         'c1415': 'cases/pip2015_14-15.csv',
         'cbaecc_test': 'cases/pip2015test.csv',
         'c1415_test': 'cases/pip2015_14-15test.csv',
         'params_cache': path.join(read.CACHE_DIR, 'param_table'+ read.MSGTLD)}


def test_events():
    return events(casesfile_baecc=files['cbaecc_test'],
                  casesfile_1415=files['c1415_test'])


def pluvio_config(e, tshift_minutes, n_comb_intervals):
    for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
        c.instr['pluvio'].shift_periods = tshift_minutes
        c.instr['pluvio'].n_combined_intervals = n_comb_intervals


def extra_events(e, extra_cases_file, extra_h5_file, *pluvio_conf_args):
    ee = sf.EventsCollection(extra_cases_file, dtformat_paper)
    ee.autoimport_data(datafile=extra_h5_file, autoshift=False, autobias=False,
                      rule='6min', varinterval=True)
    pluvio_config(ee, *pluvio_conf_args)
    e.events = e.events.append(ee.events, ignore_index=True)
    del(ee)


def events(casesfile_baecc=files['cbaecc'],
           casesfile_nov14=files['c14nov'],
           casesfile_1415=files['c1415']):
    e = sf.EventsCollection(casesfile_baecc, dtformat_paper)
    e.autoimport_data(datafile=files['h5baecc'], autoshift=False, autobias=False,
                      rule='6min', varinterval=True)
    pluvio_config(e, -6, N_COMB_INTERVALS)
    extra_events(e, casesfile_nov14, files['h5nov14'], -5, N_COMB_INTERVALS)
    extra_events(e, casesfile_1415, files['h5w1415'], -5, N_COMB_INTERVALS)
    e.events['paper'] = e.events.pluvio200
    e.split_index()
    e.events = sf.before_after_col(e.events, date=pd.datetime(2014,7,1),
                                datecol='start')
    pref400 = e.events.pluvio_pref==400
    e.events.paper[pref400] = e.events.pluvio400[pref400] # TODO
    return e


def pip2015events():
    e = events()
    #e.events.pluvio200[12] = None # TODO
    return e


def param_table(e=None, query_str=QSTR, debug=False, rho_limits=None,
                use_cache=True):
    cached_table = files['params_cache']
    if path.isfile(cached_table) and use_cache:
        return pd.read_msgpack(cached_table)
    if e is None:
        if debug:
            e = test_events()
        else:
            e = pip2015events()
    if rho_limits is None:
        rho_limits = rholimits
    data = e.summary(col='paper', split_date=pd.datetime(2014,7,1))
    del(e)
    gc.collect()
    data = sf.apply_rho_intervals(data, rho_limits)
    if len(query_str)<1:
        return data
    return data.query(query_str)
