# coding: utf-8
import pandas as pd
import datetime
import netCDF4 as nc
from matplotlib import ticker
from baecc import RHO_SCALE


def rho_scale(axis):
    axis_scale(axis, RHO_SCALE)


def axis_scale(axis, scale):
    formatter = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*scale))
    axis.set_major_formatter(formatter)


def cdf_to_series(path, vname, tname='time'):
    cdf = nc.Dataset(path)
    h = cdf.variables[vname]
    t = cdf.variables[tname]
    dt = nc.num2date(t[:], t.units)
    hs = pd.Series(h[:], index=dt)
    hs.name = vname
    return hs


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def datenum2datetime(matlab_datenum):
    """Convert MATLAB datenum to datetime."""
    return datetime.datetime.fromordinal(int(matlab_datenum)) + datetime.timedelta(days=matlab_datenum % 1) - datetime.timedelta(days=366)


def merge_series(s1, s2, **kwargs):
    """Merge pandas Series and/or DataFrames on index"""
    return pd.merge(pd.DataFrame(s1), pd.DataFrame(s2),
                    left_index=True, right_index=True, **kwargs)


def merge_multiseries(s1, s2, *series, **kwargs):
    """Merge multiple pandas Series and/or DataFrames on index"""
    s_all = merge_series(s1, s2, **kwargs)
    if len(series) > 0:
        for s in series:
            s_all = merge_series(s_all, s, **kwargs)
    return s_all
