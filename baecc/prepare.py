# -*- coding: utf-8 -*-
"""
Paper specific settings and helper functions.
@author: Jussi Tiira
"""
import gc
import numpy as np
import pandas as pd
import bisect
import baecc
from os import path
from baecc import RESULTS_DIR, DATA_DIR, USER_DIR, H5_PATH
from baecc import caching
from j24 import ensure_dir, ensure_join

N_COMB_INTERVALS = 2

dtformat_default = '%Y-%m-%d %H:%M'
dtformat_snex = '%Y %d %B %H UTC'
dtformat_paper = '%Y %b %d %H:%M'
cond = lambda df: (df.intensity>0.2) & (df.D_0_gamma>0.6) & \
                  (df.density==df.density) & (df['count']>800)
RHO_LIMITS = (0, 100, 200, 1000)
#rholimits = (0, 150, 300, 800)
resultspath = path.join(RESULTS_DIR, 'pip2015')
paperpath = path.join(resultspath, 'paper')

paths = {'results': ensure_dir(resultspath),
         'paper': paperpath,
         'tables': ensure_join(paperpath, 'tables')}
files = {'h5nov14': path.join(DATA_DIR, '2014nov1-23.h5'),
         'h5w1415': path.join(DATA_DIR, 'dec-jan1415.h5'),
         'h5baecc': H5_PATH,
         'params_cache': path.join(caching.CACHE_DIR, 'param_table'+ caching.MSGTLD)}


def cases_filepath(name):
    return path.join(USER_DIR, 'cases', name + '.csv')

def find_interval(x, limits=(0, 100, 200, 1000)):
    """Find rightmost value less than x and leftmost value greater than x."""
    i = bisect.bisect_right(limits, x)
    return limits[i-1:i+1]


def find_interval_df(s, limits):
    """Find intervals for Series s, output as a two-column DataFrame."""
    return s.apply(find_interval, limits=limits).apply(pd.Series)


def apply_rho_intervals(df, limits, rho_col='density'):
    """Add columns for density intervals."""
    data = df.copy()
    data[['rhomin', 'rhomax']] = find_interval_df(data[rho_col], limits)
    return data


def before_after_col(df, date=pd.datetime(2014,7,1), colname='winter',
                     datecol=None):
    if datecol is None:
        dates = df.index
    else:
        dates = df[datecol]
    isfirst = dates > date
    df[colname] = isfirst.astype(int)
    return df


def test_events():
    return events(casesfile_baecc=files['cbaecc_test'],
                  casesfile_1415=files['c1415_test'])


def pluvio_config(e, tshift_minutes, n_comb_intervals):
    for c in np.append(e.events.pluvio200.values, e.events.pluvio400.values):
        c.instr['pluvio'].shift_periods = tshift_minutes
        c.instr['pluvio'].n_combined_intervals = n_comb_intervals


def extra_events(e, extra_cases_file, extra_h5_file, *pluvio_conf_args):
    ee = baecc.events.EventsCollection(extra_cases_file, dtformat_default)
    ee.autoimport_data(datafile=extra_h5_file, autoshift=False, autobias=False,
                      rule='6min', varinterval=True)
    pluvio_config(ee, *pluvio_conf_args)
    e.events = e.events.append(ee.events, ignore_index=True)
    del(ee)


def events(casesname_baecc=None, casesname_nov14=None, casesname_1415=None):
    casesfile_baecc = cases_filepath(casesname_baecc)
    e = baecc.events.EventsCollection(casesfile_baecc, dtformat_default)
    e.autoimport_data(datafile=files['h5baecc'], autoshift=False, autobias=False,
                      rule='6min', varinterval=True)
    pluvio_config(e, -6, N_COMB_INTERVALS)
    #extra_events(e, casesfile_nov14, files['h5nov14'], -5, N_COMB_INTERVALS)
    if casesname_1415 is not None:
        casesfile_1415 = cases_filepath(casesname_1415)
        extra_events(e, casesfile_1415, files['h5w1415'], -5, N_COMB_INTERVALS)
    e.events['paper'] = e.events.pluvio200
    e.split_index()
    e.events = before_after_col(e.events, date=pd.datetime(2014,7,1),
                                datecol='start')
    if 'pluvio_pref' in e.events.columns:
        pref400 = e.events.pluvio_pref==400
        e.events.paper[pref400] = e.events.pluvio400[pref400].copy() # TODO
    return e


def param_table(e=None, cond=cond, debug=False, rho_limits=None,
                use_cache=True, split_date=pd.datetime(2014,7,1), **kws):
    cached_table = files['params_cache']
    if path.isfile(cached_table) and use_cache:
        return pd.read_msgpack(cached_table)
    if e is None:
        if debug:
            e = test_events(**kws)
        else:
            e = events(**kws)
    if rho_limits is None:
        rho_limits = RHO_LIMITS
    data = e.summary(col='paper', split_date=split_date)
    del(e)
    gc.collect()
    data = apply_rho_intervals(data, rho_limits)
    if cond is None:
        return data
    return data.where(cond)
